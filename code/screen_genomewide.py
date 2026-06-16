from __future__ import annotations
import argparse, json, os, time
import cobra
from modal_context import ModalContext
from benchmark import build_arm_base, EXTRA_RICH, PRODUCTS
from coupling import score_coupling
from generate import fseof, _target_rxn_id

CODE_DIR = os.path.dirname(os.path.abspath(__file__))


COUPLING_OUT = os.path.join(CODE_DIR, "data", "screen_coupling_genomewide.json")
GF = 0.5  # near-maximal growth fraction (matches ec-FBA demo)

KNOWN = {p[3:].rsplit("_", 1)[0] for p in PRODUCTS} | {"lys__L", "succ"}


def secretable_exchanges(model):
    """All organic-ish exchange reactions that can carry secretion (ub > 0)."""
    exs = []
    for r in model.reactions:
        if not (r.id.startswith("EX_") and r.id.endswith("_e")):
            continue
        if r.upper_bound <= 0:
            continue
        exs.append(r.id)
    return sorted(set(exs))


def run_coupling():
    assert cobra.__version__ == "0.31.1", f"cobra {cobra.__version__}"
    proto = ModalContext.load()
    base = proto.base
    rich_medium = {ex: 10.0 for ex in proto.medium}
    for ex in EXTRA_RICH:
        if ex in base.reactions:
            rich_medium[ex] = 10.0

    products = secretable_exchanges(base)
    print(f"[screen] {len(products)} secretable exchanges | gf={GF}", flush=True)

    plain = build_arm_base("polos", rich_medium)   # rich medium, no ATPM
    full = build_arm_base("full", rich_medium)     # M9 minimal + NGAM ATPM
    t0 = time.time()

    rows = []
    for i, ex in enumerate(products):
        try:
            sp = score_coupling(plain, ex, GF)
            sf = score_coupling(full, ex, GF)
        except Exception as e:
            print(f"  [err] {ex}: {type(e).__name__} {e}", flush=True)
            continue
        met = ex[3:].rsplit("_", 1)[0]
        rows.append({
            "exchange": ex, "metabolite": met,
            "plain_prod": round(sp["prod_at_near_max"], 4),
            "full_prod": round(sf["prod_at_near_max"], 4),
            "plain_coupling": round(sp["coupling_strength"], 4),
            "full_coupling": round(sf["coupling_strength"], 4),
            "known": met in KNOWN,
        })
        if (i + 1) % 50 == 0:
            print(f"  ...{i+1}/{len(products)} ({time.time()-t0:.0f}s)", flush=True)

    # ranking + flags
    coupled_full = [r for r in rows if r["full_coupling"] > 1e-3]
    coupled_full.sort(key=lambda r: r["full_coupling"], reverse=True)
    # plain-good-but-infeasible: high plain production, collapses under resource-aware
    infeasible_flag = [r for r in rows
                       if r["plain_prod"] > 1.0 and r["full_prod"] < 0.05 * max(r["plain_prod"], 1e-9)]
    infeasible_flag.sort(key=lambda r: r["plain_prod"], reverse=True)
    nonobvious = [r for r in coupled_full if not r["known"]]

    out = {
        "model": "iML1515", "cobra": cobra.__version__, "growth_fraction": GF,
        "n_products_screened": len(rows),
        "n_growth_coupled_full": len(coupled_full),
        "n_nonobvious_coupled": len(nonobvious),
        "n_plain_good_but_infeasible": len(infeasible_flag),
        "top20_coupled_full": coupled_full[:20],
        "top20_nonobvious_coupled": nonobvious[:20],
        "top20_plain_good_but_infeasible": infeasible_flag[:20],
        "all_rows": rows,
        "note": "plain = rich medium, no ATPM; full = M9 + NGAM (resource-aware). "
                "infeasible flag uses native model reactions.",
    }
    json.dump(out, open(COUPLING_OUT, "w"), indent=2)
    print(f"\n[done] {len(rows)} products in {time.time()-t0:.0f}s -> {COUPLING_OUT}")
    print(f"  growth-coupled (full): {len(coupled_full)} | non-obvious: {len(nonobvious)} "
          f"| plain-good-but-infeasible: {len(infeasible_flag)}")
    print("\n=== TOP 15 growth-coupled under resource-aware (M9) ===")
    for r in coupled_full[:15]:
        tag = "" if r["known"] else "  <-- NON-OBVIOUS"
        print(f"  {r['metabolite']:12s} full_coupling={r['full_coupling']:.3f} "
              f"full_prod={r['full_prod']:.2f} plain_prod={r['plain_prod']:.2f}{tag}")
    print("\n=== TOP 10 PLAIN-GOOD-BUT-INFEASIBLE (native reactions) ===")
    for r in infeasible_flag[:10]:
        print(f"  {r['metabolite']:12s} plain_prod={r['plain_prod']:.2f} "
              f"-> full_prod={r['full_prod']:.3f} (collapses under M9)")
    return 0

FSEOF_OUT = os.path.join(CODE_DIR, "data", "screen_fseof_genomewide.json")
FSEOF_PRODUCTS = ["EX_lys__L_e", "EX_succ_e"]
BIG = 10**6  # uncap target list ---> genome-wide ranking


def targets_for(model, product):
    """Full FSEOF ranked amplification targets: [(rxn_id, up_factor, genes)]."""
    cands = fseof(model, product, n_steps=10, max_targets=BIG)
    out = []
    for c in cands:
        for rid, up in c.expr_mods.items():
            genes = ""
            if rid in [r.id for r in model.reactions]:
                rx = model.reactions.get_by_id(rid)
                genes = ";".join(sorted(g.name or g.id for g in rx.genes)) if rx.genes else ""
            out.append({"rxn": rid, "up_factor": up, "genes": genes})
    return out


def run_fseof():
    assert cobra.__version__ == "0.31.1", f"cobra {cobra.__version__}"
    proto = ModalContext.load()
    base = proto.base
    rich_medium = {ex: 10.0 for ex in proto.medium}
    for ex in EXTRA_RICH:
        if ex in base.reactions:
            rich_medium[ex] = 10.0

    plain = build_arm_base("polos", rich_medium)
    full = build_arm_base("full", rich_medium)
    t0 = time.time()

    results = {}
    for prod in FSEOF_PRODUCTS:
        tp = targets_for(plain, prod)
        tf = targets_for(full, prod)
        set_p = {t["rxn"] for t in tp}
        set_f = {t["rxn"] for t in tf}
        promoted = [t for t in tf if t["rxn"] not in set_p]   # only resource-aware finds
        dropped = [t for t in tp if t["rxn"] not in set_f]    # plain-only, lost under M9
        results[prod] = {
            "n_targets_plain": len(tp), "n_targets_full": len(tf),
            "n_shared": len(set_p & set_f),
            "n_promoted_resource_aware": len(promoted),
            "n_dropped_vs_plain": len(dropped),
            "top15_full": tf[:15], "top15_plain": tp[:15],
            "promoted_top10": promoted[:10], "dropped_top10": dropped[:10],
        }
        print(f"[{prod}] plain={len(tp)} full={len(tf)} shared={len(set_p & set_f)} "
              f"promoted={len(promoted)} dropped={len(dropped)} ({time.time()-t0:.0f}s)",
              flush=True)
        print(f"  TOP5 full (resource-aware) targets:")
        for t in tf[:5]:
            print(f"    {t['rxn']:14s} up x{t['up_factor']:.2f}  genes={t['genes']}")

    out = {"model": "iML1515", "cobra": cobra.__version__,
           "method": "FSEOF (Choi et al., 2010), genome-wide, plain vs resource-aware",
           "products": FSEOF_PRODUCTS, "results": results,
           "note": "up_factor = flux end/start ratio (clip 1.5-5). plain=rich, no ATPM; "
                   "full=M9+NGAM. promoted = targets only the resource-aware oracle surfaces."}
    json.dump(out, open(FSEOF_OUT, "w"), indent=2)
    print(f"\n[done] {time.time()-t0:.0f}s -> {FSEOF_OUT}")
    return 0


UNIV = os.path.join(CODE_DIR, "data", "bigg_universal_model.json")
HET_OUT = os.path.join(CODE_DIR, "data", "screen_heterologous_audit.json")
HET_PRODUCTS = ["EX_lys__L_e", "EX_succ_e"]
HET_GF = 0.5
MARGIN = 0.5          # min plain boost to count as "looks productive"
COLLAPSE = 0.10       # full boost < COLLAPSE * plain boost ---> infeasible flag


def prod_at_floor(model, prod_ex, floor):
    """Max product secretion at growth >= floor (one LP)."""
    with model:
        model.reactions.get_by_id(BIO).lower_bound = floor
        model.objective = model.reactions.get_by_id(prod_ex)
        v = model.slim_optimize()
        return float(v) if v is not None else 0.0


def addable_reactions(univ, native_mets, native_rxns):
    out = []
    for r in univ["reactions"]:
        rid = r["id"]
        if rid in native_rxns:
            continue
        if rid.startswith(("EX_", "DM_", "SK_")):
            continue
        mets = r.get("metabolites", {})
        if not mets or len(mets) > 12:
            continue
        if all(m in native_mets for m in mets):
            out.append(r)
    return out


def add_and_score(arm, r, prod_ex, floor, baseline):
    with arm:
        rxn = cobra.Reaction(r["id"])
        rxn.lower_bound = float(r.get("lower_bound", -1000.0))
        rxn.upper_bound = float(r.get("upper_bound", 1000.0))
        arm.add_reactions([rxn])
        rxn.add_metabolites({arm.metabolites.get_by_id(m): c
                             for m, c in r["metabolites"].items()})
        return prod_at_floor(arm, prod_ex, floor) - baseline


def run_heterologous():
    assert cobra.__version__ == "0.31.1", f"cobra {cobra.__version__}"
    proto = ModalContext.load()
    base = proto.base
    native_mets = {m.id for m in base.metabolites}
    native_rxns = {r.id for r in base.reactions}
    rich_medium = {ex: 10.0 for ex in proto.medium}
    for ex in EXTRA_RICH:
        if ex in base.reactions:
            rich_medium[ex] = 10.0

    univ = json.load(open(UNIV))
    cands = addable_reactions(univ, native_mets, native_rxns)
    print(f"[audit] addable heterologous reactions: {len(cands)} (of {len(univ['reactions'])})",
          flush=True)

    plain = build_arm_base("polos", rich_medium)
    full = build_arm_base("full", rich_medium)
    global BIO
    BIO = str(plain.objective.expression).split("*")[1].split(" ")[0]
    gm_p = plain.slim_optimize(); gm_f = full.slim_optimize()
    floor_p, floor_f = HET_GF * gm_p, HET_GF * gm_f
    print(f"[arms] plain gmax={gm_p:.3f} full gmax={gm_f:.3f}", flush=True)

    results = {}
    for prod in HET_PRODUCTS:
        base_p = prod_at_floor(plain, prod, floor_p)
        base_f = prod_at_floor(full, prod, floor_f)
        print(f"\n[{prod}] baseline plain={base_p:.3f} full={base_f:.3f}", flush=True)
        flags, improvers = [], 0
        t0 = time.time()
        for i, r in enumerate(cands):
            try:
                dp = add_and_score(plain, r, prod, floor_p, base_p)
            except Exception:
                continue
            if dp <= MARGIN:
                continue
            improvers += 1
            try:
                df = add_and_score(full, r, prod, floor_f, base_f)
            except Exception:
                df = 0.0
            rec = {"rxn": r["id"], "name": r.get("name", ""),
                   "plain_boost": round(dp, 3), "full_boost": round(df, 3),
                   "infeasible": df < COLLAPSE * dp}
            flags.append(rec)
            if (i + 1) % 2000 == 0:
                print(f"  ...{i+1}/{len(cands)} improvers={improvers} ({time.time()-t0:.0f}s)",
                      flush=True)
        infeasible = [f for f in flags if f["infeasible"]]
        infeasible.sort(key=lambda f: f["plain_boost"], reverse=True)
        flags.sort(key=lambda f: f["plain_boost"], reverse=True)
        results[prod] = {
            "baseline_plain": round(base_p, 3), "baseline_full": round(base_f, 3),
            "n_improvers_plain": improvers,
            "n_flagged_infeasible": len(infeasible),
            "top15_improvers": flags[:15],
            "top15_infeasible_flags": infeasible[:15],
            "all_improvers": [{"plain_boost": f["plain_boost"], "full_boost": f["full_boost"],
                               "infeasible": f["infeasible"]} for f in flags],
        }
        print(f"  improvers(plain)={improvers}  flagged-infeasible={len(infeasible)} "
              f"({time.time()-t0:.0f}s)", flush=True)
        for f in infeasible[:8]:
            print(f"    FLAG {f['rxn']:16s} plain+{f['plain_boost']:.2f} -> full+{f['full_boost']:.2f}"
                  f"  {f['name'][:40]}", flush=True)

    out = {"model": "iML1515 + BiGG universal heterologous reactions",
           "cobra": cobra.__version__, "growth_fraction": HET_GF,
           "margin": MARGIN, "collapse_ratio": COLLAPSE,
           "n_addable": len(cands), "products": HET_PRODUCTS, "results": results,
           "note": "improver = adds plain-medium product boost > MARGIN; infeasible flag = "
                   "that boost collapses (< 10 percent) under M9 resource-aware. BiGG "
                   "universal reactions. In-silico hypothesis only."}
    json.dump(out, open(HET_OUT, "w"), indent=2)
    print(f"\n[done] -> {HET_OUT}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Unified genome-wide screens (coupling / fseof / heterologous).")
    sub = parser.add_subparsers(dest="screen", required=True)
    sub.add_parser("coupling", help="genome-wide growth-coupling screen")
    sub.add_parser("fseof", help="genome-wide FSEOF amplification-target ranking")
    sub.add_parser("heterologous", help="BiGG-universal heterologous reaction audit")
    args = parser.parse_args(argv)

    dispatch = {
        "coupling": run_coupling,
        "fseof": run_fseof,
        "heterologous": run_heterologous,
    }
    return dispatch[args.screen]()


if __name__ == "__main__":
    raise SystemExit(main())
