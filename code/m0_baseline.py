import json
import os
import sys

import cobra
from cobra.io import load_model, save_json_model

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_CACHE = os.path.join(OUT_DIR, "iML1515.json")
REPORT = os.path.join(OUT_DIR, "m0_baseline_report.json")
REF_GROWTH = 0.877  # /h, BiGG reference, glucose aerobic minimal


def get_model():
    if os.path.exists(MODEL_CACHE):
        print(f"[load] local cache {MODEL_CACHE}")
        return cobra.io.load_json_model(MODEL_CACHE)
    print("[load] downloading iML1515 from BiGG (cobra.io.load_model)...")
    m = load_model("iML1515")
    save_json_model(m, MODEL_CACHE)
    print(f"[load] saved to {MODEL_CACHE}")
    return m


def main():
    m = get_model()
    print(f"[model] id={m.id} genes={len(m.genes)} rxns={len(m.reactions)} "
          f"mets={len(m.metabolites)}")

    # Default objective & medium (iML1515 default = aerobic glucose minimal ~ M9)
    obj = str(m.objective.expression).split("-")[0].strip()
    print(f"[objective] {obj}")
    medium = {k: round(v, 4) for k, v in m.medium.items()}
    print(f"[default medium] {len(medium)} exchanges open:")
    for k, v in sorted(medium.items()):
        print(f"    {k} = -{v}")

    # FBA baseline
    sol = m.optimize()
    growth = sol.objective_value
    print(f"\n[FBA] status={sol.status}  growth={growth:.6f} /h")
    print(f"[FBA] BiGG reference ~{REF_GROWTH} /h  "
          f"delta={abs(growth - REF_GROWTH):.4f}")

    gate_ok = sol.status == "optimal" and abs(growth - REF_GROWTH) < 0.01
    print(f"[GATE M0] {'PASS' if gate_ok else 'FAIL'} "
          f"(optimal & |delta|<0.01)")

    # pFBA sanity (parsimonious flux distribution)
    obj_rxn = next(r.id for r in m.reactions if r.objective_coefficient)
    try:
        from cobra.flux_analysis import pfba
        pf = pfba(m)
        print(f"[pFBA] total |flux| = {pf.fluxes.abs().sum():.2f}  "
              f"growth={pf.fluxes[obj_rxn]:.6f}")
    except Exception as e:
        print(f"[pFBA] skip ({type(e).__name__}: {e})")

    report = {
        "model_id": m.id,
        "n_genes": len(m.genes),
        "n_reactions": len(m.reactions),
        "n_metabolites": len(m.metabolites),
        "objective": obj,
        "growth_per_h": growth,
        "ref_growth_per_h": REF_GROWTH,
        "delta": abs(growth - REF_GROWTH),
        "status": sol.status,
        "gate_pass": bool(gate_ok),
        "medium_default": medium,
        "solver": str(cobra.Configuration().solver.__name__),
    }
    with open(REPORT, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n[report] {REPORT}")
    return 0 if gate_ok else 1


if __name__ == "__main__":
    sys.exit(main())
