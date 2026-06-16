from __future__ import annotations

import copy
import json
import os
import random
from concurrent.futures import ProcessPoolExecutor

from benchmark import EPS, EXTRA_RICH, PRODUCTS, auroc
from benchmark_m3 import (GROWTH_FRACTION, N_WORKERS, SCALAR, _arm_cfg,
                          _gen_task, _init_gen_worker, _init_worker, _score_task,
                          make_fp_medium, make_fp_thermo)
from generate import Candidate
from modal_context import ModalContext

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(CODE_DIR, "e6_noise_report.json")
NOISE_LEVELS = [0.0, 0.1, 0.25, 0.5, 0.75, 1.0]
SEEDS = [0, 1, 2]


def corrupt(cands, p, rng, rxn_ids, met_ids):
    """Mis-annotate each candidate with prob p (wrong sequence-to-reaction map)."""
    out = []
    for c in cands:
        c2 = copy.deepcopy(c)
        if rng.random() < p:
            for spec in c2.pathway_rxns:           # mismatched reaction product
                prods = [m for m, co in spec["metabolites"].items() if co > 0]
                if prods:
                    coef = spec["metabolites"].pop(prods[0])
                    spec["metabolites"][rng.choice(met_ids)] = coef
            if c2.knockouts:                        # wrong GPR -> KO wrong reaction
                c2.knockouts = [rng.choice(rxn_ids) for _ in c2.knockouts]
            if c2.expr_mods:                        # wrong overexpression target
                c2.expr_mods = {rng.choice(rxn_ids): v for v in c2.expr_mods.values()}
        out.append(c2)
    return out


def score_all(cands, rich_medium, full_dG0, conc):
    """Score all candidates + WT on the full arm (parallel)."""
    cfg = _arm_cfg("full"); uc = conc if cfg["ecmdb"] else None
    wt_tasks = [(f"WT::{p}", Candidate(cid="WT", target_metabolite=p)) for p in PRODUCTS]
    tasks = wt_tasks + [(c.cid, c) for c in cands]
    out = {}
    with ProcessPoolExecutor(max_workers=N_WORKERS, initializer=_init_worker,
                             initargs=("full", rich_medium, full_dG0, cfg["tfa"], uc)) as ex:
        for cid, sv in ex.map(_score_task, tasks, chunksize=2):
            out[cid] = sv
    return out


def main():
    proto = ModalContext.load()
    rich_medium = {ex: 10.0 for ex in proto.medium}
    for ex in EXTRA_RICH:
        if ex in proto.base.reactions:
            rich_medium[ex] = 10.0
    rxn_ids = [r.id for r in proto.base.reactions if not r.id.startswith("EX_")]
    met_ids = [m.id for m in proto.base.metabolites if m.id.endswith("_c")]

    # parallel generation (reuse M3)
    cands, labels, dg_inject = [], {}, {}
    with ProcessPoolExecutor(max_workers=min(N_WORKERS, len(PRODUCTS)),
                             initializer=_init_gen_worker) as ex:
        for prod, rc, tp in ex.map(_gen_task, PRODUCTS):
            for c in rc:
                labels[c.cid] = "realistic"
            cands += rc
            for c in tp:
                labels[c.cid] = "tp_pathway"; cands.append(c)
            cm = make_fp_medium(prod); labels[cm.cid] = "fp_medium"; cands.append(cm)
            ct, dgi = make_fp_thermo(prod); labels[ct.cid] = "fp_thermo"
            cands.append(ct); dg_inject.update(dgi)
    full_dG0 = dict(proto.dG0 or {}); full_dG0.update(dg_inject)
    conc = proto.metabolite_conc
    prod_of = {c.cid: c.target_metabolite for c in cands}
    print(f"[e6] {len(cands)} cands, {len(NOISE_LEVELS)} noise levels × {len(SEEDS)} seeds, "
          f"{N_WORKERS}w", flush=True)

    def auroc_tp_vs_fp(scores, wt):
        sub = [c.cid for c in cands if labels[c.cid] in ("tp_pathway", "fp_medium", "fp_thermo")]
        sc = [scores[c][SCALAR] - wt[prod_of[c]] for c in sub]
        y = [1 if labels[c] == "tp_pathway" else 0 for c in sub]
        return auroc(sc, y)

    curve = []
    for p in NOISE_LEVELS:
        aurocs = []
        for seed in SEEDS:
            rng = random.Random(seed * 100 + int(p * 1000))
            cc = corrupt(cands, p, rng, rxn_ids, met_ids)
            sv = score_all(cc, rich_medium, full_dG0, conc)
            wt = {pr: sv[f"WT::{pr}"][SCALAR] for pr in PRODUCTS}
            scores = {c.cid: sv[c.cid] for c in cc}
            aurocs.append(auroc_tp_vs_fp(scores, wt))
        import statistics
        mean_a = statistics.mean(aurocs)
        sd_a = statistics.pstdev(aurocs) if len(aurocs) > 1 else 0.0
        curve.append({"noise_p": p, "auroc_mean": mean_a, "auroc_sd": sd_a, "auroc_seeds": aurocs})
        print(f"   p={p:.2f}  AUROC(tp vs fp) full = {mean_a:.3f} ± {sd_a:.3f}", flush=True)

    # threshold: highest p where AUROC still >= 0.7
    reliable = [pt["noise_p"] for pt in curve if pt["auroc_mean"] >= 0.7]
    threshold = max(reliable) if reliable else 0.0
    print(f"\n[E6] ambang-noise andal (AUROC≥0.7): p ≤ {threshold:.2f}")
    print(f"     AUROC p=0: {curve[0]['auroc_mean']:.3f}  ->  p=1: {curve[-1]['auroc_mean']:.3f}")

    out = {"noise_levels": NOISE_LEVELS, "seeds": SEEDS, "n_candidates": len(cands),
           "curve": curve, "reliable_threshold_p": threshold,
           "growth_fraction": GROWTH_FRACTION,
           "interpretation": "Ambang p = laju mis-anotasi maks sebelum modal gagal "
                             "memisahkan TP dari FP (AUROC<0.7). Membatasi klaim ke rezim "
                             "anotasi-andal. Label proxy, bukan wetlab."}
    with open(REPORT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[report] {REPORT}")


if __name__ == "__main__":
    main()
