from __future__ import annotations

import json
import os
from dataclasses import dataclass, field

import cobra
from scipy.stats import kendalltau, spearmanr

from generate import Candidate, generate_candidates
from modal_context import ModalContext

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(CODE_DIR, "m2_benchmark_report.json")
EPS = 1e-6

PRODUCTS = ["EX_ac_e", "EX_lac__D_e", "EX_etoh_e", "EX_for_e", "EX_succ_e",
            "EX_akg_e", "EX_pyr_e", "EX_mal__L_e", "EX_fum_e", "EX_glyclt_e"]
# rich-medium nutrients (outside M9 minimal) - uptake only in polos arm
EXTRA_RICH = ["EX_glu__L_e", "EX_asp__L_e", "EX_ala__L_e", "EX_gln__L_e",
              "EX_xyl__D_e", "EX_glyc_e", "EX_arg__L_e", "EX_ser__L_e",
              "EX_gly_e", "EX_acgam_e"]


def target_to_cyt(target_ex: str) -> str:
    """EX_succ_e -> succ_c."""
    return target_ex[3:].rsplit("_", 1)[0] + "_c"


# FP traps (independently-labelled).
def make_fp_medium(target_ex: str, rich_met="glu__L_c") -> Candidate:
    """Heterolog rxn rich_met -> target_c; collapses under M9 (no substrate)."""
    tc = target_to_cyt(target_ex)
    spec = {"id": f"TRAP_MED_{target_ex}", "metabolites": {rich_met: -1.0, tc: 1.0},
            "lb": 0.0, "ub": 1000.0}
    return Candidate(
        cid=f"FPMED_{target_ex}", target_metabolite=target_ex,
        pathway_rxns=[spec],
        provenance={"label": "fp_medium", "requires": rich_met,
                    "reason": "butuh nutrien non-M9 (glutamate) -> FP saat medium minimal"})


def make_fp_thermo(target_ex: str) -> tuple[Candidate, dict]:
    """Uphill CO2-fixation heterolog rxn; dG'° forced positive -> TFA rejects."""
    tc = target_to_cyt(target_ex)
    rid = f"TRAP_THERMO_{target_ex}"
    # uphill CO2-fix + NADH reduction; dG-infeasible route (FBA solves, dG rejects).
    spec = {"id": rid,
            "metabolites": {"co2_c": -4.0, "nadh_c": -7.0, "h_c": -6.0,
                            tc: 1.0, "nad_c": 7.0, "h2o_c": 4.0},
            "lb": 0.0, "ub": 1000.0}
    cand = Candidate(
        cid=f"FPTHERMO_{target_ex}", target_metabolite=target_ex,
        pathway_rxns=[spec],
        provenance={"label": "fp_thermo", "rxn": rid,
                    "reason": "fiksasi-CO2 uphill ΔG'°>>0 -> FP termodinamika"})
    dg_inject = {rid: (200.0, 2.0)}  # dG'° very positive (kJ/mol)
    return cand, dg_inject


# Apply candidate to model (Stage 3 - insert/KO/expr).
def apply_candidate(model: cobra.Model, cand: Candidate) -> None:
    for spec in cand.pathway_rxns:
        if spec["id"] in model.reactions:
            continue
        rxn = cobra.Reaction(spec["id"])
        rxn.lower_bound = spec.get("lb", 0.0)
        rxn.upper_bound = spec.get("ub", 1000.0)
        mets = {}
        for mid, coef in spec["metabolites"].items():
            try:
                mets[model.metabolites.get_by_id(mid)] = coef
            except KeyError:
                m = cobra.Metabolite(mid, compartment="c")
                model.add_metabolites([m])
                mets[m] = coef
        rxn.add_metabolites(mets)
        model.add_reactions([rxn])
    for rid in cand.knockouts:
        try:
            r = model.reactions.get_by_id(rid)
            r.lower_bound = r.upper_bound = 0.0
        except KeyError:
            pass
    for rid, factor in cand.expr_mods.items():
        try:
            r = model.reactions.get_by_id(rid)
            if r.upper_bound > 0:
                r.upper_bound = min(1000.0, r.upper_bound * factor)
            if r.lower_bound < 0:
                r.lower_bound = max(-1000.0, r.lower_bound * factor)
        except KeyError:
            pass


def _tfa_screen(model: cobra.Model, dG0: dict, conc: dict | None) -> int:
    """Run dG screen after candidate insertion (so trap heterolog rxns get screened)."""
    import dG as _dG
    try:
        flagged = _dG.flag_infeasible(model, dG0, conc=conc)
    except TypeError:
        flagged = _dG.flag_infeasible(model, dG0)  # legacy signature fallback
    n = 0
    for rid, direction in flagged.items():
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        if direction == "forward" and rxn.upper_bound > 0.0:
            rxn.upper_bound = 0.0; n += 1
        elif direction == "reverse" and rxn.lower_bound < 0.0:
            rxn.lower_bound = 0.0; n += 1
    return n


def score_candidate(arm_base: cobra.Model, cand: Candidate, dG0: dict,
                    apply_tfa: bool, conc: dict | None) -> dict:
    """Score target_max (objective=product) + growth feasibility."""
    with arm_base as m:
        apply_candidate(m, cand)
        n_constr = _tfa_screen(m, dG0, conc) if apply_tfa else 0
        g = m.optimize()
        growth = float(g.objective_value) if g.status == "optimal" else 0.0
        feasible = g.status == "optimal" and growth > EPS
        try:
            trx = m.reactions.get_by_id(cand.target_metabolite)
            m.objective = trx
            s = m.optimize()
            tmax = float(s.objective_value) if s.status == "optimal" else 0.0
        except KeyError:
            tmax = 0.0
    return {"growth": growth, "feasible": feasible, "target_max": tmax,
            "tfa_constrained": n_constr}


# Arm modal (ablative). Annotations (candidates+reactions) identical across arms.
# Arm config: (medium, atpm_on, apply_tfa, use_ecmdb_conc)
ARM_CFG = {
    "polos":       {"rich": True,  "atpm": False, "tfa": False, "ecmdb": False},
    "medium":      {"rich": False, "atpm": False, "tfa": False, "ecmdb": False},
    "medium_atpm": {"rich": False, "atpm": True,  "tfa": False, "ecmdb": False},
    "tfa_1mM":     {"rich": False, "atpm": True,  "tfa": True,  "ecmdb": False},
    "full":        {"rich": False, "atpm": True,  "tfa": True,  "ecmdb": True},
}


def build_arm_base(name: str, rich_medium: dict) -> cobra.Model:
    """Per-arm base: medium + ATPM only (candidate-independent constraints)."""
    cfg = ARM_CFG[name]
    mc = ModalContext.load()  # M9 medium + NGAM atpm from data
    if cfg["rich"]:
        mc.medium = rich_medium
    if not cfg["atpm"]:
        mc.atpm = 0.0
    mc.toggles = {"medium": True, "atpm": True, "ecmdb_pool": False,
                  "tfa_dG": False, "ec_kcat": False}
    return mc.apply(mc.base.copy())


# Metrics
def auroc(scores: list[float], labels: list[int]) -> float:
    """AUROC via Mann-Whitney (positive=realistic=1)."""
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return float("nan")
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def fp_rate_at_k(ranked_labels: list[str], k: int) -> float:
    top = ranked_labels[:k]
    return sum(1 for l in top if l.startswith("fp")) / max(1, len(top))


def main():
    proto = ModalContext.load()
    base_default = proto.base
    print(f"[setup] iML1515 medium-M9={len(proto.medium)} atpm={proto.atpm} "
          f"conc={len(proto.metabolite_conc or {})} dG0={len(proto.dG0 or {})}", flush=True)

    # rich medium (polos): M9 + extras (uptake 10, finite -> discriminative target_max)
    rich_medium = {ex: 10.0 for ex in proto.medium}
    for ex in EXTRA_RICH:
        if ex in base_default.reactions:
            rich_medium[ex] = 10.0
    # product secretion (ub>0) already default; product uptake stays 0

    # collect candidates (fixed annotation; generate once on M9 model)
    gen_model = proto.apply(proto.base.copy())  # M9 + atpm
    cands: list[Candidate] = []
    labels: dict[str, str] = {}
    dg_inject_all: dict = {}
    for prod in PRODUCTS:
        try:
            rc = generate_candidates(gen_model, target=prod, use_fseof=True,
                                     use_knockout=True,
                                     fseof_kwargs={"max_targets": 3},
                                     knockout_kwargs={"max_candidates": 3})
        except Exception as e:
            print(f"   gen {prod} ERR {type(e).__name__}", flush=True)
            rc = []
        for c in rc:
            labels[c.cid] = "realistic"
        cands += rc
        # traps
        cm = make_fp_medium(prod); labels[cm.cid] = "fp_medium"; cands.append(cm)
        ct, dgi = make_fp_thermo(prod); labels[ct.cid] = "fp_thermo"
        cands.append(ct); dg_inject_all.update(dgi)
    print(f"[candidates] total={len(cands)}  "
          f"realistic={sum(1 for v in labels.values() if v=='realistic')}  "
          f"fp_medium={sum(1 for v in labels.values() if v=='fp_medium')}  "
          f"fp_thermo={sum(1 for v in labels.values() if v=='fp_thermo')}", flush=True)

    # score each arm (constant annotation; only modal layer differs)
    arms = ["polos", "medium", "medium_atpm", "tfa_1mM", "full"]
    full_dG0 = dict(proto.dG0 or {}); full_dG0.update(dg_inject_all)  # +dG trap
    conc = proto.metabolite_conc
    scores = {a: {} for a in arms}   # arm -> cid -> dict
    for a in arms:
        base = build_arm_base(a, rich_medium)
        cfg = ARM_CFG[a]
        use_conc = conc if cfg["ecmdb"] else None
        nc = 0
        for c in cands:
            sv = score_candidate(base, c, full_dG0, cfg["tfa"], use_conc)
            scores[a][c.cid] = sv
            nc += sv["tfa_constrained"]
        print(f"[arm {a}] scored {len(cands)} candidates  tfa_constrained_total={nc}", flush=True)

    # WT baseline per product per arm (empty candidate)
    from generate import Candidate as _C
    wt = {a: {} for a in arms}
    for a in arms:
        base = build_arm_base(a, rich_medium)
        cfg = ARM_CFG[a]
        uc = conc if cfg["ecmdb"] else None
        for prod in PRODUCTS:
            wt[a][prod] = score_candidate(base, _C(cid="WT", target_metabolite=prod),
                                          full_dG0, cfg["tfa"], uc)["target_max"]
    prod_of = {c.cid: c.target_metabolite for c in cands}
    cids = [c.cid for c in cands]

    # H1 (global transparency) + per-product
    tm_polos = [scores["polos"][cid]["target_max"] for cid in cids]
    tm_full = [scores["full"][cid]["target_max"] for cid in cids]
    rho_g, _ = spearmanr(tm_polos, tm_full)
    tau_g, _ = kendalltau(tm_polos, tm_full)
    # per-product Spearman (intervention ranking for one product = real use-case)
    per_prod_rho = {}
    for prod in PRODUCTS:
        pc = [c for c in cids if prod_of[c] == prod]
        if len(pc) >= 3:
            a = [scores["polos"][c]["target_max"] for c in pc]
            b = [scores["full"][c]["target_max"] for c in pc]
            r, _ = spearmanr(a, b)
            per_prod_rho[prod] = r
    import statistics
    mean_prod_rho = statistics.mean([r for r in per_prod_rho.values() if r == r]) \
        if per_prod_rho else float("nan")

    # H2: "winner" = score > WT*(1+eps) -> candidate appears to improve.
    # Traps masquerade as winners under polos (inflation); modal must remove them.
    EPS_WIN = 0.02
    def winners(arm):
        return [c for c in cids
                if scores[arm][c]["target_max"] > wt[arm][prod_of[c]] * (1 + EPS_WIN) + EPS]
    def fp_among_winners(arm):
        w = winners(arm)
        if not w:
            return 0.0, 0, 0
        fp = sum(1 for c in w if labels[c].startswith("fp"))
        return fp / len(w), fp, len(w)

    fpw_polos, nfp_p, nw_p = fp_among_winners("polos")
    fpw_full, nfp_f, nw_f = fp_among_winners("full")
    # traps inflated in polos (winners) then REMOVED by modal (non-winners)
    trap_inflated_polos = [c for c in winners("polos") if labels[c].startswith("fp")]
    trap_removed_full = [c for c in trap_inflated_polos if c not in winners("full")]
    trap_catch_rate = len(trap_removed_full) / max(1, len(trap_inflated_polos))

    def topk(arm, k=10):
        return set(sorted(cids, key=lambda c: scores[arm][c]["target_max"], reverse=True)[:k])
    churn = 1 - len(topk("polos") & topk("full")) / 10
    H1_pass = (rho_g <= 0.7) and (churn >= 0.3)
    H2_pass = (fpw_full < fpw_polos) and trap_catch_rate >= 0.5

    print("\n" + "=" * 72)
    print("M2 BENCHMARK - H1 (ranking shift) & H2 (false-positive (down))")
    print("=" * 72)
    print(f"H1 global: Spearman={rho_g:.3f} Kendall={tau_g:.3f} top10-churn={churn:.2f}")
    print(f"H1 per-produk: mean Spearman={mean_prod_rho:.3f}  {{"
          + ', '.join(f'{p[3:-2]}:{r:.2f}' for p, r in per_prod_rho.items()) + "}}")
    print(f"   -> H1 {'PASS' if H1_pass else 'FAIL'}")
    print(f"H2 (winner = skor>WT, kandidat tampak-improvement):")
    print(f"   FP-rate antar-winner: polos={fpw_polos:.2f} ({nfp_p}/{nw_p})  "
          f"full={fpw_full:.2f} ({nfp_f}/{nw_f})")
    print(f"   trap-inflasi-di-polos={len(trap_inflated_polos)}  "
          f"dihapus-modal={len(trap_removed_full)}  catch-rate={trap_catch_rate:.2f}")
    print(f"   -> H2 {'PASS' if H2_pass else 'FAIL'} (FP-winner down & catch>=0.5)")

    # ablatif: FP-rate antar-winner per arm
    print("\nablatif FP-rate-antar-winner per arm:")
    for a in arms:
        f, nf, nw = fp_among_winners(a)
        print(f"    {a:14s} FP={f:.2f} ({nf}/{nw} winner)")

    out = {
        "n_candidates": len(cands),
        "label_counts": {v: sum(1 for x in labels.values() if x == v) for v in set(labels.values())},
        "H1": {"spearman_global": rho_g, "kendall_global": tau_g, "top10_churn": churn,
               "spearman_per_product": per_prod_rho, "mean_per_product": mean_prod_rho,
               "pass": bool(H1_pass)},
        "H2": {"fp_among_winners_polos": fpw_polos, "fp_among_winners_full": fpw_full,
               "winners_polos": nw_p, "winners_full": nw_f,
               "trap_inflated_polos": len(trap_inflated_polos),
               "trap_removed_full": len(trap_removed_full),
               "trap_catch_rate": trap_catch_rate, "pass": bool(H2_pass)},
        "ablative_fp_among_winners": {a: fp_among_winners(a)[0] for a in arms},
        "wt_baseline": wt, "scores": {a: scores[a] for a in arms}, "labels": labels,
        "caveat": "Label realisme = proxy (ΔG/medium by-construction), BUKAN wetlab. "
                  "H2 diukur per-produk relatif-WT karena "
                  "ranking global terkonfound flux-ceiling antar-produk. Hasil null = sah.",
    }
    with open(REPORT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n[report] {REPORT}")


if __name__ == "__main__":
    main()
