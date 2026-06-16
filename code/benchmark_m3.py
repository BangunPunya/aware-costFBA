from __future__ import annotations

import json
import os
import statistics
import sys

from scipy.stats import spearmanr

from benchmark import (ARM_CFG, EPS, EXTRA_RICH, PRODUCTS, _tfa_screen,
                       apply_candidate, auroc, build_arm_base, make_fp_medium,
                       make_fp_thermo)
from coupling import score_coupling
from generate import Candidate, generate_candidates, pathway_candidates
from modal_context import ModalContext

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(CODE_DIR, "m3_benchmark_report.json")
GROWTH_FRACTION = float(sys.argv[1]) if len(sys.argv) > 1 else 0.9  # sensitivity test
SCALAR = "prod_at_near_max"   # growth-coupling design score


def score_cand(arm_base, cand: Candidate, dG0: dict, apply_tfa: bool,
               conc: dict | None) -> dict:
    with arm_base as m:
        apply_candidate(m, cand)
        nc = _tfa_screen(m, dG0, conc) if apply_tfa else 0
        sv = score_coupling(m, cand.target_metabolite, GROWTH_FRACTION)
        sv["tfa_constrained"] = nc
    return sv


# Worker ~80% CPU, capped at 12 (RAM ~0.5GB/worker). Override via BENCH_WORKERS.
N_WORKERS = int(os.environ.get("BENCH_WORKERS",
                               min(int(0.8 * (os.cpu_count() or 4)), 12)))


def _build_arm(a, rich_medium):
    base = build_arm_base("full" if a == "full_anaer" else a, rich_medium)
    if a == "full_anaer":
        try:
            base.reactions.EX_o2_e.lower_bound = 0.0
        except KeyError:
            pass
    return base


def _arm_cfg(a):
    return ARM_CFG["full" if a == "full_anaer" else a]


_WK: dict = {}


def _init_worker(arm, rich_medium, dG0, apply_tfa, conc):
    global _WK
    _WK = {"base": _build_arm(arm, rich_medium), "dG0": dG0,
           "tfa": apply_tfa, "conc": conc}


def _score_task(payload):
    cid, cand = payload
    return cid, score_cand(_WK["base"], cand, _WK["dG0"], _WK["tfa"], _WK["conc"])


_GWK: dict = {}


def _init_gen_worker():
    global _GWK
    proto = ModalContext.load()
    _GWK = {"gm": proto.apply(proto.base.copy())}


def _gen_task(prod):
    gm = _GWK["gm"]
    try:
        rc = generate_candidates(gm, target=prod, use_fseof=True, use_knockout=True,
                                 fseof_kwargs={"max_targets": 3},
                                 knockout_kwargs={"max_candidates": 3})
    except Exception:
        rc = []
    tp = pathway_candidates(gm, prod)
    return prod, rc, tp


def main():
    proto = ModalContext.load()
    print(f"[setup] M9={len(proto.medium)} atpm={proto.atpm} conc={len(proto.metabolite_conc or {})} "
          f"dG0={len(proto.dG0 or {})}", flush=True)
    rich_medium = {ex: 10.0 for ex in proto.medium}
    for ex in EXTRA_RICH:
        if ex in proto.base.reactions:
            rich_medium[ex] = 10.0

    from concurrent.futures import ProcessPoolExecutor
    cands, labels, dg_inject = [], {}, {}
    gen_workers = min(N_WORKERS, len(PRODUCTS))
    with ProcessPoolExecutor(max_workers=gen_workers, initializer=_init_gen_worker) as ex:
        for prod, rc, tp in ex.map(_gen_task, PRODUCTS):  # parallel generation (~10x)
            for c in rc:
                labels[c.cid] = "realistic"
            cands += rc
            for c in tp:
                labels[c.cid] = "tp_pathway"; cands.append(c)
            cm = make_fp_medium(prod); labels[cm.cid] = "fp_medium"; cands.append(cm)
            ct, dgi = make_fp_thermo(prod); labels[ct.cid] = "fp_thermo"
            cands.append(ct); dg_inject.update(dgi)
    counts = {v: sum(1 for x in labels.values() if x == v) for v in set(labels.values())}
    print(f"[candidates] {len(cands)}  {counts}  (gen parallel {gen_workers}w)", flush=True)

    # full_anaer = full modal BUT O2 off (anaerobic) -> tests true coupling +
    # PDC-ethanol becomes TP (active fermentation without O2).
    arms = ["polos", "medium_atpm", "full", "full_anaer"]
    full_dG0 = dict(proto.dG0 or {}); full_dG0.update(dg_inject)
    conc = proto.metabolite_conc
    scores = {a: {} for a in arms}
    wt = {a: {} for a in arms}
    wt_full = {a: {} for a in arms}   # full dict (for coupling_strength)
    prod_of = {c.cid: c.target_metabolite for c in cands}

    # Per-arm pool (worker builds 1 base model = RAM-efficient). Each candidate + WT =
    # independent task.
    from concurrent.futures import ProcessPoolExecutor
    for a in arms:
        cfg = _arm_cfg(a); uc = conc if cfg["ecmdb"] else None
        wt_tasks = [(f"WT::{p}", Candidate(cid="WT", target_metabolite=p)) for p in PRODUCTS]
        tasks = wt_tasks + [(c.cid, c) for c in cands]
        with ProcessPoolExecutor(max_workers=N_WORKERS, initializer=_init_worker,
                                 initargs=(a, rich_medium, full_dG0, cfg["tfa"], uc)) as ex:
            for cid, sv in ex.map(_score_task, tasks, chunksize=2):
                if cid.startswith("WT::"):
                    p = cid[4:]; wt[a][p] = sv[SCALAR]; wt_full[a][p] = sv
                else:
                    scores[a][cid] = sv
        print(f"[arm {a}] scored (parallel {N_WORKERS}w)", flush=True)

    cids = [c.cid for c in cands]

    def winners(arm):
        return [c for c in cids if scores[arm][c][SCALAR] > wt[arm][prod_of[c]] * 1.02 + EPS]

    def fp_among_winners(arm):
        w = winners(arm)
        if not w:
            return 0.0, 0, 0
        fp = sum(1 for c in w if labels[c].startswith("fp"))
        return fp / len(w), fp, len(w)

    def tp_recall(arm):  # fraction of tp_pathway that become winners
        tps = [c for c in cids if labels[c] == "tp_pathway"]
        if not tps:
            return float("nan"), 0, 0
        w = set(winners(arm))
        hit = sum(1 for c in tps if c in w)
        return hit / len(tps), hit, len(tps)

    # Claim-2: AUROC separating tp_pathway vs fp_* using coupling score
    def auroc_tp_vs_fp(arm):
        sub = [c for c in cids if labels[c] in ("tp_pathway", "fp_medium", "fp_thermo")]
        sc = [scores[arm][c][SCALAR] - wt[arm][prod_of[c]] for c in sub]  # lift over WT
        y = [1 if labels[c] == "tp_pathway" else 0 for c in sub]
        return auroc(sc, y)

    fpw_p, nfp_p, nw_p = fp_among_winners("polos")
    fpw_f, nfp_f, nw_f = fp_among_winners("full")
    tpr_p = tp_recall("polos"); tpr_f = tp_recall("full")
    auc_p = auroc_tp_vs_fp("polos"); auc_f = auroc_tp_vs_fp("full")

    # H1 ranking shift (coupling score)
    sp = [scores["polos"][c][SCALAR] for c in cids]
    sf = [scores["full"][c][SCALAR] for c in cids]
    rho, _ = spearmanr(sp, sf)

    H2_pass = fpw_f < fpw_p
    K2_pass = (auc_f > auc_p) and auc_f >= 0.7 and tpr_f[0] >= 0.5

    print("\n" + "=" * 72)
    print("M3 BENCHMARK - growth-coupling score + TP-enrichment (Claim-2)")
    print("=" * 72)
    print(f"skor = {SCALAR} (produksi pd ≥{GROWTH_FRACTION}·growth_max)")
    print(f"H1: Spearman(polos,full) coupling-score = {rho:.3f}")
    print(f"H2 (FP-suppression): FP-rate antar-winner  polos={fpw_p:.2f}({nfp_p}/{nw_p})  "
          f"full={fpw_f:.2f}({nfp_f}/{nw_f})  -> {'PASS' if H2_pass else 'FAIL'}")
    print(f"Klaim-2 (TP-enrichment):")
    print(f"   AUROC(tp vs fp): polos={auc_p:.3f}  full={auc_f:.3f}")
    print(f"   TP-recall (tp jadi winner): polos={tpr_p[0]:.2f}({tpr_p[1]}/{tpr_p[2]})  "
          f"full={tpr_f[0]:.2f}({tpr_f[1]}/{tpr_f[2]})")
    print(f"   -> Claim-2 {'PASS' if K2_pass else 'FAIL'} (AUROC_full>polos & >=0.7 & TP-recall>=0.5)")

    # detail TP per arm
    print("\nTP_pathway skor (lift atas WT) per arm:")
    for c in [c for c in cids if labels[c] == "tp_pathway"]:
        lifts = {a: round(scores[a][c][SCALAR] - wt[a][prod_of[c]], 3) for a in arms}
        print(f"   {c:28s} {lifts}")

    # ANAEROB (full modal, O2=0): true coupling + PDC becomes TP
    auc_an = auroc_tp_vs_fp("full_anaer"); tpr_an = tp_recall("full_anaer")
    fpw_an, nfp_an, nw_an = fp_among_winners("full_anaer")
    print("\n--- ANAEROB (full modal, O2=0) ---")
    print(f"   AUROC(tp vs fp)={auc_an:.3f}  TP-recall={tpr_an[0]:.2f}({tpr_an[1]}/{tpr_an[2]})  "
          f"FP-among-winner={fpw_an:.2f}({nfp_an}/{nw_an})")
    print("   coupling_strength WT (0=production optional, >0=OBLIGATORY) aerob(full)->anaerob:")
    for prod in ["EX_succ_e", "EX_etoh_e", "EX_for_e", "EX_lac__D_e", "EX_ac_e"]:
        ca = wt_full["full"][prod]; an = wt_full["full_anaer"][prod]
        print(f"     {prod:14s} growth {ca['growth_max']:.3f}->{an['growth_max']:.3f}  "
              f"coupling {ca['coupling_strength']:.3f}->{an['coupling_strength']:.3f}  "
              f"prod_min {ca['prod_min_guaranteed']:.2f}->{an['prod_min_guaranteed']:.2f}")
    print("   PDC-etanol lift: full=%.3f  anaerob=%.3f" % (
        scores["full"]["PATH_EX_etoh_e_PDC_homoethanol"][SCALAR] - wt["full"]["EX_etoh_e"],
        scores["full_anaer"]["PATH_EX_etoh_e_PDC_homoethanol"][SCALAR] - wt["full_anaer"]["EX_etoh_e"]))

    out = {
        "n_candidates": len(cands), "label_counts": counts,
        "score_scalar": SCALAR, "growth_fraction": GROWTH_FRACTION,
        "H1_spearman_coupling": rho,
        "H2": {"fp_winner_polos": fpw_p, "fp_winner_full": fpw_f, "pass": bool(H2_pass)},
        "Klaim2": {"auroc_tp_vs_fp_polos": auc_p, "auroc_tp_vs_fp_full": auc_f,
                   "auroc_tp_vs_fp_anaer": auc_an,
                   "tp_recall_polos": tpr_p[0], "tp_recall_full": tpr_f[0],
                   "tp_recall_anaer": tpr_an[0], "pass": bool(K2_pass)},
        "anaerobic": {"auroc": auc_an, "tp_recall": tpr_an[0], "fp_among_winners": fpw_an,
                      "coupling_strength_wt": {p: {"full": wt_full["full"][p]["coupling_strength"],
                                                   "anaer": wt_full["full_anaer"][p]["coupling_strength"]}
                                               for p in PRODUCTS}},
        "ablative_fp_among_winners": {a: fp_among_winners(a)[0] for a in arms},
        "wt_baseline": wt, "scores": {a: scores[a] for a in arms}, "labels": labels,
        "caveat": "Skor coupling = prod pd near-max growth (proxy titer realistik). TP = "
                  "library heterolog kurated (PYC/PDC), bukan enumerasi RetroPath/ATLAS. "
                  "Label proxy, bukan wetlab.",
    }
    with open(REPORT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n[report] {REPORT}")


if __name__ == "__main__":
    main()
