from __future__ import annotations

import json
import os
from typing import Iterable

import cobra  # noqa: F401  (version asserted in __main__)

import annotate
from coupling import score_coupling
from generate import Candidate, HETEROLOGOUS_LIBRARY
from benchmark import apply_candidate
from modal_context import ModalContext

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(CODE_DIR, "data")
OUT_JSON = os.path.join(DATA_DIR, "loop_results.json")
OUT_FASTA = os.path.join(DATA_DIR, "loop_ranked.fasta")

GROWTH_FRACTION_DEFAULT = 0.9
THROTTLE_FRAC_DEFAULT = 0.50   # illustrative (per demo_fbr_lysine); WT throttle level
EPS = 1e-6


# Known feedback-inhibited committed steps (intervention-type router).
# Each maps a reaction_id -> the gene/biology note that makes it a committed,
# allosterically-feedback-inhibited throttle point (fbr-relaxation class).
# DHDPS (dapA/b2478): single gene, no isozyme escape, first committed step of the
# DAP/lysine branch. ASPK is excluded: it has 3 isozymes (lysC/thrA/metL), so a
# single cap reroutes through thrA/metL and is not lysine-specific.
FEEDBACK_COMMITTED_STEPS: dict[str, dict] = {
    "DHDPS": {
        "gene": "dapA (b2478)",
        "metabolite": "EX_lys__L_e",
        "why": "first committed, single-gene, lysine-specific step of the DAP/lysine "
               "branch; allosterically inhibited by L-lysine in E. coli; "
               "isozyme-free throttle.",
    },
}


# fbr throttle/relax helper (demo_fbr_lysine semantics).
# Self-contained copy of the math only; the biological justification + sweep
# live in demo_fbr_lysine.py.
def _unconstrained_capacity(model: cobra.Model, rxn_id: str,
                            growth_fraction: float) -> float:
    """Max achievable flux of rxn_id while holding biomass >= frac*growth_max (the ceiling the fbr mutant restores toward; same as demo_fbr_lysine._dhdps_unconstrained_capacity, generalised)."""
    with model:
        sol = model.optimize()
        gmax = float(sol.objective_value) if sol.objective_value is not None else 0.0
        if gmax <= EPS:
            return 0.0
        for br in (r for r in model.reactions if r.objective_coefficient != 0):
            br.lower_bound = growth_fraction * gmax
        model.objective = model.reactions.get_by_id(rxn_id)
        model.objective_direction = "max"
        sol2 = model.optimize()
        cap = float(sol2.objective_value) if sol2.status == "optimal" else 0.0
    return cap


def _make_throttled(base: cobra.Model, rxn_id: str, cap: float,
                    frac: float) -> cobra.Model:
    """Wild-type regulated model: cap the committed-step UB to frac*capacity
    (allosteric feedback inhibition on). Reuses demo_fbr_lysine._make_wt."""
    m = base.copy()
    m.reactions.get_by_id(rxn_id).upper_bound = frac * cap
    return m


def _make_relaxed(base: cobra.Model, rxn_id: str, cap: float) -> cobra.Model:
    """fbr candidate: relax the throttle back to the unconstrained capacity
    (feedback-resistant point mutation). Reuses demo_fbr_lysine._make_fbr."""
    m = base.copy()
    m.reactions.get_by_id(rxn_id).upper_bound = cap
    return m


# FASTA input parsing
def _parse_fasta(path: str) -> list[tuple[str, str]]:
    """Read a FASTA file -> list of (id, dna_sequence)."""
    out: list[tuple[str, str]] = []
    cur_id, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id is not None:
                    out.append((cur_id, "".join(buf)))
                cur_id = line[1:].split()[0] if len(line) > 1 else f"seq{len(out)}"
                buf = []
            elif line.strip():
                buf.append(line.strip())
    if cur_id is not None:
        out.append((cur_id, "".join(buf)))
    return out


def _collect_inputs(dna_inputs) -> list[tuple[str, str]]:
    """Normalise dna_inputs (FASTA path, or list of (id, dna) tuples and/or FASTA paths) into a flat list of (id, dna) tuples."""
    if isinstance(dna_inputs, str):
        return _parse_fasta(dna_inputs)
    flat: list[tuple[str, str]] = []
    for item in dna_inputs:
        if isinstance(item, str):
            flat.extend(_parse_fasta(item))
        else:
            sid, seq = item
            flat.append((str(sid), str(seq)))
    return flat


# cascade output -> intervention TYPE + Candidate
def _classify_intervention(reaction_id: str, model: cobra.Model,
                           target_id: str) -> str:
    """Route a recovered reaction to an intervention class by its biology."""
    if reaction_id in FEEDBACK_COMMITTED_STEPS:
        return "fbr-relaxation"
    if reaction_id in model.reactions:
        return "overexpression"
    # reaction not in model: heterologous if the target has a curated spec
    if HETEROLOGOUS_LIBRARY.get(target_id):
        return "heterologous"
    return "unmappable-novel"


def _build_candidate(cid: str, reaction_id: str, intervention_type: str,
                     target_id: str, confidence: float,
                     low_confidence: bool, prov_extra: dict) -> Candidate | None:
    """Map a recovered reaction -> generate.Candidate for the chosen class; returns None for 'unmappable-novel'."""
    prov = {
        "source": "loop:annotate-cascade",
        "recovered_reaction": reaction_id,
        "intervention_type": intervention_type,
        "cascade_confidence": confidence,
        "low_confidence": low_confidence,
        **prov_extra,
    }
    if intervention_type == "overexpression":
        # raise the recovered reaction's flux upper bound (expr_mods factor);
        # usually Δ=0 in stoichiometric FBA (not rate-limiting / isozyme-buffered).
        prov["note"] = ("overexpression = raise flux UB; pure overexpression is "
                        "commonly Δ=0 in stoichiometric FBA (step not rate-limiting "
                        "or isozyme-buffered).")
        return Candidate(cid=cid, target_metabolite=target_id,
                         expr_mods={reaction_id: 5.0}, provenance=prov)
    if intervention_type == "heterologous":
        specs = HETEROLOGOUS_LIBRARY.get(target_id, [])
        if not specs:
            return None
        name, rxn_specs = specs[0]
        prov["heterologous_name"] = name
        prov["note"] = ("recovered reaction not in iML1515; routed via curated "
                        "HETEROLOGOUS_LIBRARY pathway spec.")
        return Candidate(cid=cid, target_metabolite=target_id,
                         pathway_rxns=[dict(s) for s in rxn_specs], provenance=prov)
    if intervention_type == "fbr-relaxation":
        # fbr handled specially in scoring (throttle WT vs relax candidate).
        prov["throttle_reaction"] = reaction_id
        prov["note"] = ("feedback-inhibited committed step; WT = throttled "
                        "(regulation ON), candidate = throttle relaxed (fbr). "
                        "Reuses demo_fbr_lysine throttle/relax.")
        return Candidate(cid=cid, target_metabolite=target_id, provenance=prov)
    return None  # unmappable-novel


# scoring: product @ near-max growth for WT and candidate under one oracle
def _prod_at_near_max(model: cobra.Model, target_id: str,
                      growth_fraction: float) -> tuple[float, float, bool]:
    """(prod_at_near_max, growth_max, feasible) via coupling.score_coupling."""
    r = score_coupling(model, target_id, growth_fraction=growth_fraction)
    return (float(r["prod_at_near_max"]), float(r["growth_max"]), bool(r["feasible"]))


def _score_one_oracle(base: cobra.Model, cand: Candidate,
                      intervention_type: str, target_id: str,
                      growth_fraction: float, throttle_frac: float) -> dict:
    """Score WT vs candidate under ONE oracle. fbr-relaxation: WT=throttled, candidate=relaxed (plain oracle: no throttle -> Δ≈0). overexpression: WT=base, candidate=base+expr_mods. heterologous: WT=base, candidate=base+pathway rxn."""
    rxn_id = cand.provenance.get("recovered_reaction")

    if intervention_type == "fbr-relaxation":
        cap = _unconstrained_capacity(base, rxn_id, growth_fraction)
        wt = _make_throttled(base, rxn_id, cap, throttle_frac)
        cd = _make_relaxed(base, rxn_id, cap)
        wt_prod, wt_g, wt_feas = _prod_at_near_max(wt, target_id, growth_fraction)
        cd_prod, cd_g, cd_feas = _prod_at_near_max(cd, target_id, growth_fraction)
    else:
        # overexpression / heterologous: WT is the untouched base, candidate is
        # the mutated model. apply_candidate mutates a COPY (we use `with`).
        wt_prod, wt_g, wt_feas = _prod_at_near_max(base, target_id, growth_fraction)
        with base as m:
            apply_candidate(m, cand)
            cd_prod, cd_g, cd_feas = _prod_at_near_max(m, target_id, growth_fraction)

    delta = round(cd_prod - wt_prod, 6)
    growth_retained = bool(cd_g > EPS and wt_g > EPS
                           and abs(cd_g - wt_g) / max(wt_g, 1e-9) < 0.01)
    return {
        "WT": round(wt_prod, 6),
        "cand": round(cd_prod, 6),
        "delta": delta,
        "WT_growth": round(wt_g, 6),
        "cand_growth": round(cd_g, 6),
        "feasible": bool(wt_feas and cd_feas),
        "growth_retained": growth_retained,
    }


# the loop
def run_loop(target_id: str,
             dna_inputs,
             growth_fraction: float = GROWTH_FRACTION_DEFAULT,
             throttle_frac: float = THROTTLE_FRAC_DEFAULT,
             write_outputs: bool = True) -> list[dict]:
    """Run the full closed loop for one target over many DNA inputs.

    target_id: product exchange rxn id (e.g. 'EX_lys__L_e'). dna_inputs: FASTA
    path OR list of (id, dna) tuples and/or FASTA paths. growth_fraction: hold
    biomass >= frac*growth_max when scoring. throttle_frac: WT throttle level for
    fbr-relaxation (illustrative - see demo_fbr_lysine). write_outputs: write JSON+FASTA.
    Returns ranked list of per-input dicts (modal Δ desc, conf tiebreak).
    """
    os.makedirs(DATA_DIR, exist_ok=True)
    inputs = _collect_inputs(dna_inputs)

    # Build the two oracle base models ONCE (reuse for every input).
    mc = ModalContext.load()
    base_modal = mc.apply(mc.base.copy())          # modal regulatory layer ON
    polos = mc.as_polos()
    base_plain = polos.apply(polos.base.copy())    # plain stoichiometric FBA (no modal)
    # cascade annotation uses the cached model + genemap (DIAMOND if present).
    annot_model = mc.base                            # untouched bounds for annotate/GPR
    genemap = annotate.load_genemap()
    modal_layers = [k for k, v in mc.toggles.items() if v]

    results: list[dict] = []
    for sid, dna in inputs:
        res = annotate.annotate(dna, annot_model, genemap=genemap)
        top = res.reactions[0] if res.reactions else None
        recovered_rxn = top.reaction_id if top else None
        confidence = res.confidence
        low_conf = res.low_confidence

        row: dict = {
            "id": sid,
            "recovered_reaction": recovered_rxn,
            "all_recovered": [r.reaction_id for r in res.reactions],
            "confidence": confidence,
            "low_confidence": low_conf,
            "intervention_type": None,
            "plain_WT": None, "plain_cand": None, "plain_delta": None,
            "modal_WT": None, "modal_cand": None, "modal_delta": None,
            "growth_retained": None,
            "scored": False,
            "skip_reason": None,
            "provenance": {
                "dna_len_bp": len(dna),
                "n_orfs": len(res.orfs),
                "cascade_method": res.method,
                "modal_layers_active": modal_layers,
                "throttle_frac_ILLUSTRATIVE": throttle_frac,
                "growth_fraction": growth_fraction,
            },
            "dna": dna,
        }

        if recovered_rxn is None:
            row["skip_reason"] = "no reaction recovered by cascade (non-coding / no hit)"
            results.append(row)
            continue

        itype = _classify_intervention(recovered_rxn, base_modal, target_id)
        row["intervention_type"] = itype

        cid = f"LOOP_{sid}_{recovered_rxn}"
        cand = _build_candidate(cid, recovered_rxn, itype, target_id,
                                confidence, low_conf,
                                prov_extra=FEEDBACK_COMMITTED_STEPS.get(recovered_rxn, {}))
        if cand is None:
            row["skip_reason"] = (
                "unmappable-novel: reaction not in iML1515 and no curated "
                "HETEROLOGOUS_LIBRARY spec for target - not scored (no fabrication)."
            )
            results.append(row)
            continue

        # score under BOTH oracles
        plain = _score_one_oracle(base_plain, cand, itype, target_id,
                                  growth_fraction, throttle_frac)
        modal = _score_one_oracle(base_modal, cand, itype, target_id,
                                  growth_fraction, throttle_frac)

        row["plain_WT"] = plain["WT"]
        row["plain_cand"] = plain["cand"]
        row["plain_delta"] = plain["delta"]
        row["modal_WT"] = modal["WT"]
        row["modal_cand"] = modal["cand"]
        row["modal_delta"] = modal["delta"]
        row["growth_retained"] = modal["growth_retained"]
        row["scored"] = True
        row["candidate"] = cand.as_dict()
        results.append(row)

    # rank: modal Δ desc, tiebreak cascade confidence desc
    scored = [r for r in results if r["scored"]]
    unscored = [r for r in results if not r["scored"]]
    scored.sort(key=lambda r: (r["modal_delta"], r["confidence"]), reverse=True)
    for i, r in enumerate(scored, start=1):
        r["rank"] = i
    for r in unscored:
        r["rank"] = None

    ranked = scored + unscored

    if write_outputs:
        _write_json(ranked, target_id, growth_fraction, throttle_frac, modal_layers)
        _write_fasta(scored, target_id)

    return ranked


# outputs
def _write_json(ranked: list[dict], target_id: str, growth_fraction: float,
                throttle_frac: float, modal_layers: list[str]) -> None:
    payload = {
        "title": "Full closed-loop pipeline: DNA -> annotate -> Candidate -> "
                 "plain+modal FBA oracle -> ranking -> lab-ready FASTA",
        "target_id": target_id,
        "growth_fraction": growth_fraction,
        "throttle_frac_ILLUSTRATIVE": throttle_frac,
        "modal_layers_active": modal_layers,
        "caveats": [
            "Cascade confidence propagates to the final ranking; low_confidence flagged.",
            "This is in-silico re-ranking, not a viability / titre prediction.",
            "Pure-overexpression Δ=0 is a stoichiometric-FBA limitation, reported.",
            "fbr throttle magnitude is illustrative (needs a measured Ki); "
            "qualitative claim: plain Δ≈0, modal Δ>0 when the throttle binds.",
        ],
        "results": [{k: v for k, v in r.items() if k != "dna"} for r in ranked],
        "cobra_version": cobra.__version__,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(payload, f, indent=2)


def _write_fasta(scored: list[dict], target_id: str) -> None:
    """Ranked candidate DNA sequences, header carries full provenance."""
    lines: list[str] = []
    for r in scored:
        hdr = (f">rank{r['rank']}|{r['id']}|rxn={r['recovered_reaction']}"
               f"|conf={r['confidence']}|modalΔ={r['modal_delta']}"
               f"|growth={'yes' if r['growth_retained'] else 'no'} lab-ready")
        lines.append(hdr)
        seq = r["dna"]
        for i in range(0, len(seq), 70):
            lines.append(seq[i:i + 70])
    with open(OUT_FASTA, "w") as f:
        f.write("\n".join(lines) + ("\n" if lines else ""))


# pretty report
def print_report(ranked: list[dict], target_id: str) -> None:
    line = "-" * 104
    print("\n" + line)
    print("FULL CLOSED-LOOP PIPELINE - ranked candidates for target %s" % target_id)
    print(line)
    hdr = "%-6s | %-7s | %-9s | %5s | %5s | %-16s | %10s | %10s | %s"
    print(hdr % ("rank", "id", "rxn", "conf", "lowC", "intervention", "plain Δ",
                 "modal Δ", "growth"))
    print(line)
    row = "%-6s | %-7s | %-9s | %5s | %5s | %-16s | %10s | %10s | %s"
    for r in ranked:
        rank = str(r["rank"]) if r["rank"] is not None else "-"
        pdelta = "%+.5f" % r["plain_delta"] if r["plain_delta"] is not None else "   n/a"
        mdelta = "%+.5f" % r["modal_delta"] if r["modal_delta"] is not None else "   n/a"
        growth = ("yes" if r["growth_retained"] else "no") if r["scored"] else "-"
        print(row % (rank, r["id"], str(r["recovered_reaction"]),
                     "%.2f" % r["confidence"], "Y" if r["low_confidence"] else "n",
                     r["intervention_type"] or "-", pdelta, mdelta, growth))
        if not r["scored"]:
            print("       -> SKIPPED: %s" % r["skip_reason"])
    print(line)
    print("JSON  ->", OUT_JSON)
    print("FASTA ->", OUT_FASTA)
    print(line + "\n")


if __name__ == "__main__":
    # Minimal self-demo if run directly (full validation lives in test_loop.py).
    print("loop.py - run test_loop.py for the lysine validation demo.")
    print("cobra version:", cobra.__version__)
