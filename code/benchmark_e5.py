from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor

from cobra.io import load_model

from benchmark import auroc
from benchmark_m3 import (GROWTH_FRACTION, N_WORKERS, SCALAR, _tfa_screen,
                          apply_candidate, make_fp_medium, make_fp_thermo)
from generate import Candidate, generate_candidates, pathway_candidates
from modal_context import ModalContext
from score import score_model  # noqa (ensure import path ok)
from coupling import score_coupling

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(CODE_DIR, "e5_coverage_report.json")
# iML1515_noECMDB = same E.coli but ECMDB off -> separates DATA-effect
# (E.coli +/-ECMDB, same organism) from ORGANISM-effect (E.coli-ECMDB vs non-E.coli).
# (YMDB yeast not used: concentrations absent from bulk dump + /concentrations dead ->
#  yeast metabolome data not FBA-ready, reinforces R3.)
ORGANISMS = ["iML1515", "iML1515_noECMDB", "iJO1366", "iMM904", "iYO844", "iJN1463"]
PRODUCTS_E5 = ["EX_ac_e", "EX_succ_e", "EX_etoh_e", "EX_lac__D_e", "EX_mal__L_e", "EX_fum_e"]

_DG = json.load(open(os.path.join(CODE_DIR, "data", "dG_cache.json")))


def org_modal(mid):
    """ModalContext + dG0 + conc for one organism (iML1515 = full ECMDB, else no conc)."""
    if mid in ("iML1515", "iML1515_noECMDB"):
        mc = ModalContext.load()       # medium MediaDB + ATPM NGAM + ECMDB + dG
        if mid == "iML1515_noECMDB":   # ablate ECMDB -> dG@1mM (same treatment as non-E.coli)
            mc.metabolite_conc = None
            mc.toggles["ecmdb_pool"] = False
        return mc
    m = load_model(mid)
    dg = {}
    for rid, v in _DG.items():
        v = v.get("dG0") if isinstance(v, dict) else v
        if rid in m.reactions and isinstance(v, (list, tuple)):
            dg[rid] = (float(v[0]), float(v[1]))
    atpm = 0.0
    if "ATPM" in [r.id for r in m.reactions]:
        atpm = abs(m.reactions.ATPM.lower_bound) or 1.0
    mc = ModalContext(base=m, medium=dict(m.medium), atpm=atpm, dG0=dg,
                       metabolite_conc=None,  # ECMDB E.coli-only -> MISSING
                      toggles={"medium": True, "atpm": True, "ecmdb_pool": False,
                               "tfa_dG": True, "ec_kcat": False})
    return mc


_OWK: dict = {}


def _oinit(mid):
    global _OWK
    mc = org_modal(mid)
    base = mc.base
    products = [p for p in PRODUCTS_E5 if p in base.reactions]
    full_dG0 = dict(mc.dG0 or {})
    # full = modal available; polos = rich medium + ATPM0 + no TFA
    full_model = mc.apply(base.copy())
    rich = {ex.id: 10.0 for ex in base.exchanges}
    polos_mc = ModalContext(base=base, medium=rich, atpm=0.0,
                            toggles={"medium": True, "atpm": True, "ecmdb_pool": False,
                                     "tfa_dG": False, "ec_kcat": False})
    polos_model = polos_mc.apply(base.copy())
    gm = mc.apply(base.copy())  # utk generation
    _OWK = {"mc": mc, "full": full_model, "polos": polos_model, "gm": gm,
            "dG0": full_dG0, "conc": mc.metabolite_conc, "products": products,
            "tfa": bool(mc.toggles.get("tfa_dG")), "ecmdb": bool(mc.toggles.get("ecmdb_pool"))}


def _gen_one(prod):
    gm = _OWK["gm"]
    try:
        rc = generate_candidates(gm, target=prod, use_fseof=True, use_knockout=True,
                                 fseof_kwargs={"max_targets": 2},
                                 knockout_kwargs={"max_candidates": 2})
    except Exception:
        rc = []
    tp = pathway_candidates(gm, prod)
    return prod, rc, tp


def _score_one(payload):
    arm, cid, cand = payload
    base = _OWK[arm]
    conc = _OWK["conc"] if (arm == "full" and _OWK["ecmdb"]) else None
    apply_tfa = (arm == "full") and _OWK["tfa"]
    with base as m:
        apply_candidate(m, cand)
        if apply_tfa:
            _tfa_screen(m, _OWK["dG0"], conc)
        sv = score_coupling(m, cand.target_metabolite, GROWTH_FRACTION)
    return arm, cid, sv


def run_org(mid):
    """Run benchmark for one organism; return concise advantage summary."""
    # parallel generation (per product) - pool with organism initializer
    cands, labels = [], {}
    mc0 = org_modal(mid)
    products = [p for p in PRODUCTS_E5 if p in mc0.base.reactions]

    nW = min(N_WORKERS, max(1, len(products)))
    with ProcessPoolExecutor(max_workers=nW, initializer=_oinit, initargs=(mid,)) as ex:
        for prod, rc, tp in ex.map(_gen_one, products):
            for c in rc:
                labels[c.cid] = "realistic"
            cands += rc
            for c in tp:
                labels[c.cid] = "tp_pathway"; cands.append(c)
            cm = make_fp_medium(prod); labels[cm.cid] = "fp_medium"; cands.append(cm)
            ct, _ = make_fp_thermo(prod); labels[ct.cid] = "fp_thermo"; cands.append(ct)
    prod_of = {c.cid: c.target_metabolite for c in cands}

    # parallel scoring: WT + cands x {polos, full}
    tasks = []
    for arm in ("polos", "full"):
        for p in products:
            tasks.append((arm, f"WT::{p}", Candidate(cid="WT", target_metabolite=p)))
        for c in cands:
            tasks.append((arm, c.cid, c))
    sc = {"polos": {}, "full": {}}
    with ProcessPoolExecutor(max_workers=N_WORKERS, initializer=_oinit, initargs=(mid,)) as ex:
        for arm, cid, sv in ex.map(_score_one, tasks, chunksize=4):
            sc[arm][cid] = sv

    def auc(arm):
        wt = {p: sc[arm][f"WT::{p}"][SCALAR] for p in products}
        sub = [c.cid for c in cands if labels[c.cid] in ("tp_pathway", "fp_medium", "fp_thermo")]
        s = [sc[arm][c][SCALAR] - wt[prod_of[c]] for c in sub]
        y = [1 if labels[c] == "tp_pathway" else 0 for c in sub]
        return auroc(s, y)

    a_polos, a_full = auc("polos"), auc("full")
    n_tp = sum(1 for v in labels.values() if v == "tp_pathway")
    return {"organism": mid, "n_candidates": len(cands), "n_products": len(products),
            "n_tp": n_tp, "auroc_polos": a_polos, "auroc_full": a_full,
            "modal_advantage": (a_full - a_polos) if (a_full == a_full and a_polos == a_polos) else None,
            "ecmdb": mid == "iML1515"}


def main():
    print(f"[E5] {len(ORGANISMS)} organisme × {len(PRODUCTS_E5)} produk, {N_WORKERS}w", flush=True)
    results = []
    for mid in ORGANISMS:
        try:
            r = run_org(mid)
        except Exception as e:
            r = {"organism": mid, "error": f"{type(e).__name__}: {str(e)[:100]}"}
        results.append(r)
        if "error" in r:
            print(f"   {mid}: ERROR {r['error']}", flush=True)
        else:
            print(f"   {mid:9s} ECMDB={'Y' if r['ecmdb'] else 'N'}  AUROC polos={r['auroc_polos']:.3f} "
                  f"full={r['auroc_full']:.3f}  advantage={r['modal_advantage']:+.3f}  "
                  f"(cands={r['n_candidates']}, tp={r['n_tp']})", flush=True)

    print("\n[E5] degradasi cakupan:")
    ec = [r for r in results if r.get("ecmdb")]
    non = [r for r in results if (not r.get("ecmdb")) and "error" not in r]
    if ec and non:
        import statistics
        print(f"   E.coli(+ECMDB) advantage = {ec[0]['modal_advantage']:+.3f}")
        adv = [r['modal_advantage'] for r in non if r['modal_advantage'] is not None]
        if adv:
            print(f"   non-E.coli(−ECMDB) advantage = {statistics.mean(adv):+.3f} (mean of {len(adv)})")

    out = {"organisms": ORGANISMS, "products": PRODUCTS_E5, "results": results,
           "interpretation": "modal_advantage = AUROC(tp vs fp) full − polos. ECMDB-conc "
                             "(kolam-metabolit terukur) hanya E.coli; di luar itu ΔG screen "
                             "jatuh ke 1mM default. Mengukur peluruhan manfaat-modal. "
                             "ΔG (chemistry) tetap ada lintas-organisme via abbreviasi BiGG shared."}
    with open(REPORT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[report] {REPORT}")


if __name__ == "__main__":
    main()
