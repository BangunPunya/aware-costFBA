from __future__ import annotations

import cobra
from cobra.flux_analysis import flux_variability_analysis

EPS = 1e-6


def _objective_reactions(model: cobra.Model) -> list[cobra.Reaction]:
    """Reaksi dengan objective_coefficient != 0 (umumnya reaksi biomassa)."""
    return [r for r in model.reactions if r.objective_coefficient != 0]


def score_coupling(model: cobra.Model, target_id: str,
                   growth_fraction: float = 0.9, run_fva: bool = True) -> dict:
    """Score growth-coupling for one product reaction.

    model: COBRA model (NOT mutated - uses `with model:` context). target_id: product
    reaction id (e.g. "EX_succ_e"). growth_fraction: fraction of growth_max to hold
    (default 0.9). run_fva: True -> FVA; False -> 2 manual LPs.

    Returns dict with EXACT keys: growth_max, feasible, prod_at_near_max,
    prod_min_guaranteed, coupling_strength, status.
    """
    out = {
        "growth_max": 0.0,
        "feasible": False,
        "prod_at_near_max": 0.0,
        "prod_min_guaranteed": 0.0,
        "coupling_strength": 0.0,
        "status": "",
    }

    # target hadir?
    try:
        model.reactions.get_by_id(target_id)
    except KeyError:
        out["status"] = "no_target"
        return out

    with model:
        # 1) growth_max pada objektif (biomassa) bawaan
        sol = model.optimize()
        out["status"] = sol.status
        growth_max = float(sol.objective_value) if sol.objective_value is not None else 0.0
        out["growth_max"] = growth_max
        out["feasible"] = (sol.status == "optimal") and (growth_max > EPS)
        if not out["feasible"]:
            return out

        # 2) tahan biomassa >= growth_fraction * growth_max
        biomass_rxns = _objective_reactions(model)
        if not biomass_rxns:
            out["status"] = "no_objective"
            return out

        if run_fva:
            # FVA pada fraction_of_optimum mengikat objektif (biomassa) saat ini ke
            # >= growth_fraction*growth_max, lalu min/max reaksi target -> tepat 2 LP.
            try:
                fva = flux_variability_analysis(
                    model, reaction_list=[target_id],
                    fraction_of_optimum=growth_fraction)
                prod_max = float(fva.loc[target_id, "maximum"])
                prod_min = float(fva.loc[target_id, "minimum"])
            except Exception as e:  # fallback ke 2 LP manual
                out["status"] = f"fva_failed:{type(e).__name__}"
                return _two_lp(model, target_id, growth_fraction, growth_max,
                               biomass_rxns, out)
        else:
            return _two_lp(model, target_id, growth_fraction, growth_max,
                           biomass_rxns, out)

    # produksi = sekresi positif; clip negatif (uptake) ke 0 untuk "guaranteed prod"
    prod_max = max(prod_max, 0.0)
    prod_min = max(prod_min, 0.0)
    out["prod_at_near_max"] = prod_max
    out["prod_min_guaranteed"] = prod_min
    out["coupling_strength"] = prod_min / max(prod_max, 1e-9)
    return out


def _two_lp(model: cobra.Model, target_id: str, growth_fraction: float,
            growth_max: float, biomass_rxns: list[cobra.Reaction],
            out: dict) -> dict:
    """Jalur 2-LP manual: kunci biomassa, optimasi target max lalu min."""
    # asumsi MVP: objektif tunggal-biomassa -> pasang lower_bound.
    for br in biomass_rxns:
        br.lower_bound = growth_fraction * growth_max
    trx = model.reactions.get_by_id(target_id)
    model.objective = trx

    model.objective_direction = "max"
    smax = model.optimize()
    prod_max = float(smax.objective_value) if smax.status == "optimal" else 0.0

    model.objective_direction = "min"
    smin = model.optimize()
    prod_min = float(smin.objective_value) if smin.status == "optimal" else 0.0

    prod_max = max(prod_max, 0.0)
    prod_min = max(prod_min, 0.0)
    out["prod_at_near_max"] = prod_max
    out["prod_min_guaranteed"] = prod_min
    out["coupling_strength"] = prod_min / max(prod_max, 1e-9)
    if not out["status"]:
        out["status"] = smin.status
    return out
