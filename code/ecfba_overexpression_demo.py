from __future__ import annotations
import json, os, sys
import cobra

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
MODEL = os.path.join(DATA, "eciML1515_batch_fixed.mat")
OUT_JSON = os.path.join(DATA, "ecfba_overexpression_results.json")

PRODUCT_EX = "EX_lys__L_e"      # L-lysine secretion (forward leg of GECKO irrev split)
GROWTH_FRAC = 0.5               # fix growth >= gf * growth_max, then maximise product
OE_FACTORS = [2.0, 5.0, 10.0]   # graded overexpression multipliers

# label, gene, UniProt, on lysine pathway?
TARGETS = [
    ("dapA/DHDPS", "dapA", "P0A6L2", True),   # 4-hydroxy-tetrahydrodipicolinate synthase
    ("lysC/ASPK",  "lysC", "P08660", True),   # aspartokinase III (one of 3 isozymes)
    ("gltA/CS",    "gltA", "P0ABH7", False),  # citrate synthase (off-pathway control)
]


def near_max_product(model, prod_ex: str, gf: float) -> float:
    """Fix growth >= gf*growth_max, then maximise product secretion."""
    with model:
        bio = model.reactions.get_by_id(BIOMASS)
        model.objective = bio
        gmax = model.slim_optimize()
        if gmax is None or gmax <= 1e-9:
            return float("nan")
        bio.lower_bound = gf * gmax
        model.objective = model.reactions.get_by_id(prod_ex)
        val = model.slim_optimize()
        return float(val) if val is not None else float("nan")


def overexpress(model, uniprot: str, cap: float):
    """Add a dedicated, free supply of prot_<uniprot> (engineered copies on top of the
    global proteome budget), capped at `cap`. Returns the added reaction."""
    met = model.metabolites.get_by_id(f"prot_{uniprot}")
    rxn = cobra.Reaction(f"OE_prot_{uniprot}")
    rxn.lower_bound = 0.0
    rxn.upper_bound = cap
    rxn.add_metabolites({met: 1.0})   # --> prot_X
    model.add_reactions([rxn])
    return rxn


def baseline_usage(model, uniprot: str, gf: float) -> float:
    """draw_prot_X flux at the WT product-max solution (the pool-funded usage)."""
    draw = f"draw_prot_{uniprot}"
    with model:
        bio = model.reactions.get_by_id(BIOMASS)
        model.objective = bio
        gmax = model.slim_optimize()
        bio.lower_bound = gf * gmax
        model.objective = model.reactions.get_by_id(PRODUCT_EX)
        sol = model.optimize()
        return float(sol.fluxes.get(draw, 0.0)) if sol.status == "optimal" else 0.0


def run() -> int:
    assert cobra.__version__ == "0.31.1", f"cobra {cobra.__version__} != 0.31.1"
    model = cobra.io.load_matlab_model(MODEL)
    global BIOMASS
    BIOMASS = str(model.objective.expression).split("*")[1].split(" ")[0]
    print(f"model loaded: {len(model.reactions)} rxns | biomass={BIOMASS}")

    gmax = model.slim_optimize()
    prod_wt = near_max_product(model, PRODUCT_EX, GROWTH_FRAC)
    print(f"WT: growth_max={gmax:.4f}/h | lysine_at_{GROWTH_FRAC}xgrowth={prod_wt:.5f}")

    results = []
    for label, gene, up, on_path in TARGETS:
        if f"prot_{up}" not in [m.id for m in model.metabolites]:
            print(f"[skip] {label}: prot_{up} absent"); continue
        base = baseline_usage(model, up, GROWTH_FRAC)
        rec = {"label": label, "gene": gene, "uniprot": up, "on_pathway": on_path,
               "baseline_draw": round(base, 8), "prod_wt": round(prod_wt, 6),
               "delta_by_factor": {}}
        # graded overexpression: extra copies = (f-1) * baseline usage
        for f in OE_FACTORS:
            cap = max((f - 1.0) * base, 0.0)
            with model:
                overexpress(model, up, cap)
                p = near_max_product(model, PRODUCT_EX, GROWTH_FRAC)
            rec["delta_by_factor"][f"{f:g}x"] = round(p - prod_wt, 6)
        # uncapped relaxation = "is this enzyme capacity-limiting at all?" shadow test
        with model:
            overexpress(model, up, 1000.0)
            p_unc = near_max_product(model, PRODUCT_EX, GROWTH_FRAC)
        rec["delta_uncapped"] = round(p_unc - prod_wt, 6)
        rec["capacity_limiting"] = rec["delta_uncapped"] > 1e-3
        results.append(rec)
        print(f"[{label}] baseline_draw={base:.3e} "
              f"delta(2x/5x/10x)={[rec['delta_by_factor'][k] for k in rec['delta_by_factor']]} "
              f"delta_uncapped={rec['delta_uncapped']:.5f} "
              f"capacity_limiting={rec['capacity_limiting']}")

    out = {"model": "eciML1515_batch (GECKO)", "cobra": cobra.__version__,
           "product": PRODUCT_EX, "growth_frac": GROWTH_FRAC,
           "growth_max": round(gmax, 5), "prod_wt": round(prod_wt, 6),
           "oe_factors": OE_FACTORS, "results": results,
           "note": "plain stoichiometric FBA gives delta=0 for all overexpression "
                   "(no enzyme bound exists to relax); ec-FBA scores capacity-limiting ones."}
    json.dump(out, open(OUT_JSON, "w"), indent=2)
    print(f"\nsaved -> {OUT_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(run())
