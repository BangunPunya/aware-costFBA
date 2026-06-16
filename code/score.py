from __future__ import annotations

from dataclasses import dataclass, field

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba

EPS = 1e-6


@dataclass
class ScoreVector:
    growth: float = 0.0
    feasible: bool = False
    pfba_total_flux: float | None = None
    target_id: str | None = None
    target_max: float | None = None      # theoretical max product (FBA)
    target_fva_min: float | None = None  # FVA at ~optimum growth
    target_fva_max: float | None = None
    fva_width: float | None = None       # robustness (narrow = robust)
    blocked: bool = False
    status: str = ""
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {k: v for k, v in self.__dict__.items()}


def score_model(model: cobra.Model, target_id: str | None = None,
                fva_fraction: float = 0.95, run_fva: bool = True) -> ScoreVector:
    """Score one modified model. target_id = optional product reaction (EX_* / rxn)."""
    sv = ScoreVector(target_id=target_id)

    sol = model.optimize()
    sv.status = sol.status
    sv.growth = float(sol.objective_value) if sol.objective_value is not None else 0.0
    sv.feasible = (sol.status == "optimal") and (sv.growth > EPS)

    if not sv.feasible:
        sv.notes.append(f"infeasible/no-growth (status={sol.status}, growth={sv.growth:.4g})")
        return sv

    # pFBA - metabolic load
    try:
        pf = pfba(model)
        sv.pfba_total_flux = float(pf.fluxes.abs().sum())
    except Exception as e:
        sv.notes.append(f"pFBA failed: {type(e).__name__}")

    if target_id:
        try:
            trx = model.reactions.get_by_id(target_id)
        except KeyError:
            sv.notes.append(f"target {target_id} not in model")
            return sv
        # Theoretical max product: optimize target as objective (growth free)
        with model:
            model.objective = trx
            psol = model.optimize()
            sv.target_max = float(psol.objective_value) if psol.status == "optimal" else None
        sv.blocked = (sv.target_max is None) or (abs(sv.target_max) < EPS)
        if sv.blocked:
            sv.notes.append("target blocked (max flux ~0)")
        # FVA at ~optimum growth (robustness)
        if run_fva and not sv.blocked:
            try:
                fva = flux_variability_analysis(
                    model, reaction_list=[trx], fraction_of_optimum=fva_fraction)
                sv.target_fva_min = float(fva.loc[target_id, "minimum"])
                sv.target_fva_max = float(fva.loc[target_id, "maximum"])
                sv.fva_width = sv.target_fva_max - sv.target_fva_min
            except Exception as e:
                sv.notes.append(f"FVA failed: {type(e).__name__}")
    return sv
