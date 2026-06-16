from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba

EPS = 1e-6


# Candidate dataclass
@dataclass
class Candidate:
    cid: str
    target_metabolite: str             # product BiGG id (e.g. "succ_e" / "EX_succ_e")
    pathway_rxns: list = field(default_factory=list)   # ReactionSpec list
    knockouts: list[str] = field(default_factory=list) # reaction ids -> bounds 0
    expr_mods: dict[str, float] = field(default_factory=dict)  # rxn_id -> rescale factor
    rbs_variants: list = field(default_factory=list)   # RBSDesign list
    provenance: dict = field(default_factory=dict)     # source + generator params

    def as_dict(self) -> dict:
        return {
            "cid": self.cid,
            "target_metabolite": self.target_metabolite,
            "pathway_rxns": self.pathway_rxns,
            "knockouts": self.knockouts,
            "expr_mods": self.expr_mods,
            "rbs_variants": self.rbs_variants,
            "provenance": self.provenance,
        }


# Curated central-carbon reactions (subsystem field empty in iML1515 JSON, so
# listed manually). Used as the knockout search space; kept small so |KO|<=2.
CENTRAL_CARBON_RXNS = [
    # glycolysis / gluconeogenesis
    "PFK", "FBA", "TPI", "GAPD", "PGK", "ENO", "PYK",
    # PPP
    "G6PDH2r", "GND", "EDA",
    # pyruvate / acetyl-CoA branch
    "PDH", "PFL", "LDH_D", "ALCD2x", "ACALD", "ACKr", "PTAr",
    # TCA + anaplerotic (key to succinate coupling)
    "CS", "ICDHyr", "AKGDH", "SUCOAS", "SUCDi", "FUM", "MDH",
    "PPC", "PPCK", "ME1", "ME2", "ICL", "MALS",
]


# util
def _target_rxn_id(model: cobra.Model, target: str) -> str:
    """Accept 'EX_succ_e' or 'succ_e' -> return the existing exchange reaction id."""
    if target in model.reactions:
        return target
    guess = f"EX_{target}" if not target.startswith("EX_") else target
    if guess in model.reactions:
        return guess
    raise KeyError(f"target exchange not in model: {target} / {guess}")


def _max_product_flux(model: cobra.Model, target_rxn: str) -> float:
    with model:
        model.objective = model.reactions.get_by_id(target_rxn)
        sol = model.optimize()
        return float(sol.objective_value) if sol.status == "optimal" else 0.0


# Strategi 1 - FSEOF
def fseof(
    model: cobra.Model,
    target: str = "EX_succ_e",
    n_steps: int = 8,
    max_fraction: float = 0.9,
    min_fraction: float = 0.1,
    rise_factor: float = 1.05,
    max_targets: int = 8,
) -> list[Candidate]:
    """Flux Scanning based on Enforced Objective Flux (FSEOF, Choi et al. 2010). theoretical_max = max product flux; enforce product flux = f*theoretical_max per step while keeping growth objective + pFBA; reactions whose |flux| rises monotonically (>= rise_factor, in >= ceil(0.6*steps) steps) -> overexpression targets (expr_mods rxn_id -> up_factor = end/start, clipped [1.5, 5]). Noise-tolerant; no cofactor analysis."""
    cands: list[Candidate] = []
    trx_id = _target_rxn_id(model, target)
    tmax = _max_product_flux(model, trx_id)
    if tmax <= EPS:
        return cands  # product not producible on this medium -> no FSEOF

    fractions = [
        min_fraction + (max_fraction - min_fraction) * i / (n_steps - 1)
        for i in range(n_steps)
    ]

    # flux_profiles[rxn_id] = list of |v| per step
    flux_profiles: dict[str, list[float]] = {}
    valid_steps = 0
    for f in fractions:
        with model:
            trx = model.reactions.get_by_id(trx_id)
            enforced = f * tmax
            # enforce product minimum = enforced (growth stays the objective)
            trx.lower_bound = enforced
            try:
                pf = pfba(model)
            except Exception:
                continue
            if pf.fluxes is None:
                continue
            valid_steps += 1
            for rid, v in pf.fluxes.items():
                flux_profiles.setdefault(rid, []).append(abs(float(v)))

    if valid_steps < 3:
        return cands

    need_rises = max(2, int(round(0.6 * (valid_steps - 1))))
    scored: list[tuple[str, float]] = []
    for rid, profile in flux_profiles.items():
        if len(profile) < 3:
            continue
        if rid == trx_id:
            continue  # the product itself, not an expression target
        if max(profile) < EPS:
            continue
        # count rises between steps
        rises = sum(
            1 for a, b in zip(profile, profile[1:])
            if b > a * rise_factor and b > EPS
        )
        drops = sum(
            1 for a, b in zip(profile, profile[1:])
            if a > b * rise_factor and a > EPS
        )
        # monotonic-rise: enough rises, few drops
        if rises >= need_rises and drops <= 1:
            start = next((x for x in profile if x > EPS), EPS)
            end = profile[-1]
            up_factor = max(1.5, min(5.0, end / start if start > EPS else 5.0))
            scored.append((rid, up_factor))

    # sort: largest up_factor first (strongest coupling signal)
    scored.sort(key=lambda kv: kv[1], reverse=True)
    for i, (rid, up) in enumerate(scored[:max_targets]):
        cands.append(Candidate(
            cid=f"FSEOF_{trx_id}_{i:02d}",
            target_metabolite=trx_id,
            expr_mods={rid: round(up, 3)},
            provenance={
                "strategy": "FSEOF",
                "target_rxn": trx_id,
                "theoretical_max": round(tmax, 4),
                "n_steps": valid_steps,
                "enforced_fraction_range": [min_fraction, max_fraction],
                "rise_factor": rise_factor,
                "up_target_rxn": rid,
                "note": "FSEOF (loop+pFBA); noise-tolerant monotonic criterion.",
            },
        ))
    return cands


# Strategi 2 - Knockout enumeration (OptKnock-lite)
def _coupling_screen(
    model: cobra.Model, trx_id: str, ko_set: tuple[str, ...],
    growth_fraction: float = 0.9, run_fva: bool = True,
) -> dict | None:
    """Knock out ko_set -> measure growth and production window at near-max growth (FVA); prod_min > EPS -> production coupled, prod_max rises -> capacity increases. Returns metrics dict, or None if lethal/infeasible. run_fva=False for fast growth-only screen."""
    with model:
        for rid in ko_set:
            try:
                r = model.reactions.get_by_id(rid)
            except KeyError:
                return None
            r.lower_bound = r.upper_bound = 0.0
        sol = model.optimize()
        if sol.status != "optimal" or (sol.objective_value or 0.0) < EPS:
            return None
        growth = float(sol.objective_value)
        if not run_fva:
            return {"growth": growth, "prod_min": None, "prod_max": None}
        try:
            fva = flux_variability_analysis(
                model, reaction_list=[trx_id],
                fraction_of_optimum=growth_fraction,
            )
        except Exception:
            return None
        pmin = float(fva.loc[trx_id, "minimum"])
        pmax = float(fva.loc[trx_id, "maximum"])
    return {"growth": growth, "prod_min": pmin, "prod_max": pmax}


# Subset for double-KO pairs (TCA/anaplerotic/fermentation - most relevant to
# succinate coupling). Bounds the FVA combinatorics.
DOUBLE_KO_SUBSET = [
    "SUCDi", "FUM", "MDH", "PPC", "PPCK", "ME1", "ME2",
    "ICL", "MALS", "PFL", "LDH_D", "ACKr", "PTAr", "ALCD2x", "AKGDH",
]


def knockout_candidates(
    model: cobra.Model,
    target: str = "EX_succ_e",
    rxn_pool: list[str] | None = None,
    max_ko: int = 2,
    double_pool: list[str] | None = None,
    growth_fraction: float = 0.9,
    min_growth_keep: float = 0.05,
    max_candidates: int = 12,
) -> list[Candidate]:
    """Enumerate single + double KOs, screen growth-coupling/capacity via FVA. Two stages: (1) fast growth-only screen drops lethal/near-lethal; (2) FVA on survivors -> prod_min (coupling) & prod_max (capacity). Emit if viable AND (prod_min > EPS [coupled] OR prod_max > WT + EPS [capacity rises]). FVA screen, not bilevel OptKnock MILP; true succinate coupling needs anaerobic + many KOs (>=5 in literature); on aerobic medium typically only capacity rises."""
    cands: list[Candidate] = []
    trx_id = _target_rxn_id(model, target)
    pool = rxn_pool if rxn_pool is not None else CENTRAL_CARBON_RXNS
    pool = [r for r in pool if r in model.reactions]
    dpool = double_pool if double_pool is not None else DOUBLE_KO_SUBSET
    dpool = [r for r in dpool if r in model.reactions]

    # WT baseline (full FVA once)
    wt = _coupling_screen(model, trx_id, (), growth_fraction)
    wt_pmin = wt["prod_min"] if wt else 0.0
    wt_pmax = wt["prod_max"] if wt else 0.0

    # --- stage 1: collect candidate KO-sets (growth-only screen) ---
    ko_sets: list[tuple[str, ...]] = [(r,) for r in pool]
    if max_ko >= 2:
        ko_sets += list(itertools.combinations(dpool, 2))

    survivors: list[tuple[str, ...]] = []
    for ko_set in ko_sets:
        m = _coupling_screen(model, trx_id, ko_set, growth_fraction, run_fva=False)
        if m is None or m["growth"] < min_growth_keep:
            continue
        survivors.append(ko_set)

    # --- stage 2: FVA on survivors only ---
    results: list[tuple[tuple[str, ...], dict]] = []
    for ko_set in survivors:
        m = _coupling_screen(model, trx_id, ko_set, growth_fraction, run_fva=True)
        if m is None:
            continue
        coupled = m["prod_min"] > EPS
        capacity_gain = m["prod_max"] > wt_pmax + EPS
        if coupled or capacity_gain:
            results.append((ko_set, m))

    # ranking: coupled first, then capacity
    results.sort(
        key=lambda km: (km[1]["prod_min"], km[1]["prod_max"], km[1]["growth"]),
        reverse=True,
    )

    for i, (ko_set, m) in enumerate(results[:max_candidates]):
        coupled = m["prod_min"] > EPS
        cands.append(Candidate(
            cid=f"KO_{trx_id}_{i:02d}",
            target_metabolite=trx_id,
            knockouts=list(ko_set),
            provenance={
                "strategy": "knockout-lite (OptKnock-style screen)",
                "target_rxn": trx_id,
                "ko_size": len(ko_set),
                "growth_at_ko": round(m["growth"], 5),
                "prod_min_at_near_max_growth": round(m["prod_min"], 5),
                "prod_max_at_near_max_growth": round(m["prod_max"], 5),
                "wt_prod_min": round(wt_pmin, 5),
                "wt_prod_max": round(wt_pmax, 5),
                "coupled": coupled,                # True = production obligatory
                "capacity_gain": round(m["prod_max"] - wt_pmax, 5),
                "growth_fraction": growth_fraction,
                "note": "FVA screen, not a bilevel OptKnock MILP. coupled=False "
                        "means only capacity rises (coupling needs anaerobic/more KOs).",
            },
        ))
    return cands


# Strategi 3 - Pathway generation (RetroPath2.0 / ATLAS) - STUB
# Curated heterologous pathway library (literature, ΔG-feasible), yield-increasing.
# ReactionSpec = {id, metabolites, lb, ub} (compatible with benchmark.apply_candidate);
# metabolites use iML1515 ids.
# PYC = pyruvate carboxylase (absent in E. coli; E. coli uses PPC) - fixes CO2 to
#       OAA, raising C4 yield (succ/mal/fum). Source: Corynebacterium/Lactococcus.
# PDC = pyruvate decarboxylase (Zymomonas) -> homoethanol route, raises ethanol yield.
_PYC_SPEC = {"id": "PYC_HET",
             "metabolites": {"pyr_c": -1.0, "co2_c": -1.0, "atp_c": -1.0, "h2o_c": -1.0,
                             "oaa_c": 1.0, "adp_c": 1.0, "pi_c": 1.0, "h_c": 1.0},
             "lb": 0.0, "ub": 1000.0}
_PDC_SPEC = {"id": "PDC_HET",
             "metabolites": {"pyr_c": -1.0, "h_c": -1.0, "acald_c": 1.0, "co2_c": 1.0},
             "lb": 0.0, "ub": 1000.0}

HETEROLOGOUS_LIBRARY: dict[str, list[tuple[str, list[dict]]]] = {
    "EX_succ_e":   [("PYC_anaplerotic", [_PYC_SPEC])],
    "EX_mal__L_e": [("PYC_anaplerotic", [_PYC_SPEC])],
    "EX_fum_e":    [("PYC_anaplerotic", [_PYC_SPEC])],
    "EX_etoh_e":   [("PDC_homoethanol", [_PDC_SPEC])],
}


def pathway_candidates(
    model: cobra.Model,
    target: str,
    **kwargs,
) -> list[Candidate]:
    """Yield-increasing heterologous pathways from curated library (ΔG-feasible); inserted via pathway_rxns (ReactionSpec), thermodynamically feasible so modal screen pass expected; annotation_conf = 1.0. Automated enumeration (RetroPath2.0 / ATLAS / novoStoic2.0) is future work."""
    trx_id = _target_rxn_id(model, target)
    cands: list[Candidate] = []
    for i, (name, specs) in enumerate(HETEROLOGOUS_LIBRARY.get(trx_id, [])):
        # skip if any substrate/product metabolite is absent from the model
        ok = all(mid in model.metabolites for spec in specs for mid in spec["metabolites"])
        if not ok:
            continue
        cands.append(Candidate(
            cid=f"PATH_{trx_id}_{name}",
            target_metabolite=trx_id,
            pathway_rxns=[dict(s) for s in specs],
            provenance={"strategy": "pathway-curated (heterolog)",
                        "name": name, "annotation_conf": 1.0,
                        "label": "tp_pathway",
                        "note": "ΔG-feasible curated heterologous reaction."},
        ))
    return cands


# orchestrator
def generate_candidates(
    model: cobra.Model,
    target: str = "EX_succ_e",
    use_fseof: bool = True,
    use_knockout: bool = True,
    use_pathway: bool = False,
    fseof_kwargs: dict | None = None,
    knockout_kwargs: dict | None = None,
) -> list[Candidate]:
    """Run selected strategies on model, return combined Candidate list. Model not mutated (all changes scoped via with)."""
    out: list[Candidate] = []
    if use_fseof:
        out.extend(fseof(model, target, **(fseof_kwargs or {})))
    if use_knockout:
        out.extend(knockout_candidates(model, target, **(knockout_kwargs or {})))
    if use_pathway:
        out.extend(pathway_candidates(model, target))
    return out
