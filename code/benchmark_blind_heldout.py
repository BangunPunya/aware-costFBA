from __future__ import annotations
import json, os, random
import numpy as np
from scipy import stats
import cobra

CODE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(CODE, "data")
MODEL = os.path.join(D, "eciML1515_batch_fixed.mat")
GENE_REF = os.path.join(D, "iml1515_gene_ref.fasta")
PRODUCT = "EX_lys__L_e"
GF = 0.5
N_NEG = 60
SEED = 42

POS_GENES = ["dapA", "dapB", "dapD", "dapE", "dapF", "lysA", "lysC", "asd", "ppc", "aspC"]
EXCLUDE_FROM_NEG = set(POS_GENES) | {"thrA", "metL", "dapH", "argD", "lysB", "lysP"}


def gene_to_uniprot():
    out = {}
    for ln in open(GENE_REF):
        if ln.startswith(">"):
            p = ln[1:].strip().split("|")
            if len(p) >= 3:
                out[p[1]] = p[2]
    return out


def prod_at_floor(model, floor):
    with model:
        model.reactions.get_by_id(BIO).lower_bound = floor
        model.objective = model.reactions.get_by_id(PRODUCT)
        v = model.slim_optimize()
        return float(v) if v is not None else 0.0


def overexpress_delta(model, uniprot, floor, baseline):
    mid = f"prot_{uniprot}"
    if mid not in [m.id for m in model.metabolites]:
        return None
    with model:
        rxn = cobra.Reaction(f"OE_{uniprot}")
        rxn.lower_bound = 0.0; rxn.upper_bound = 1000.0
        model.add_reactions([rxn])
        rxn.add_metabolites({model.metabolites.get_by_id(mid): 1.0})
        return prod_at_floor(model, floor) - baseline


def auroc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    return float(np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg]))


def distribution_free_stats(pos, neg, n_perm=20000, seed=SEED):
    """Small-N robust stats on raw Delta scores: Mann-Whitney U (one-sided),
    label-permutation null for AUROC, and rank-biserial effect size. None of
    these assume normality or estimate a covariance from the positive class,
    so they hold where DeLong/parametric tests are underpowered at small N."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    obs = auroc(pos, neg)
    U, p_mwu = stats.mannwhitneyu(pos, neg, alternative="greater")
    allx = np.concatenate([pos, neg])
    lab = np.array([1] * len(pos) + [0] * len(neg))
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        perm = rng.permutation(lab)
        if auroc(allx[perm == 1], allx[perm == 0]) >= obs:
            ge += 1
    p_perm = (ge + 1) / (n_perm + 1)
    return {
        "mannwhitney_U": float(U),
        "mannwhitney_p_onesided": float(p_mwu),
        "permutation_p": float(p_perm),
        "permutation_n": int(n_perm),
        "permutation_note": f"{ge}/{n_perm} label-permutations reached observed AUROC",
        "rank_biserial_r": round(2 * obs - 1, 4),
        "n_pairwise": int(len(pos) * len(neg)),
    }


def run():
    assert cobra.__version__ == "0.31.1", f"cobra {cobra.__version__}"
    g2u = gene_to_uniprot()
    model = cobra.io.load_matlab_model(MODEL)
    global BIO
    BIO = str(model.objective.expression).split("*")[1].split(" ")[0]
    draw = {r.id.replace("draw_prot_", "") for r in model.reactions if r.id.startswith("draw_prot_")}

    gm = model.slim_optimize(); floor = GF * gm
    base = prod_at_floor(model, floor)
    print(f"[blind] eciML1515 gmax={gm:.4f} floor={floor:.4f} baseline_lysine={base:.4f}", flush=True)

    # genes that have an enzyme in the ec-model (scorable)
    scorable = [g for g, u in g2u.items() if u in draw]
    rng = random.Random(SEED)
    neg_pool = [g for g in scorable if g not in EXCLUDE_FROM_NEG]
    neg_genes = rng.sample(neg_pool, min(N_NEG, len(neg_pool)))

    def score(genes):
        out = {}
        for g in genes:
            d = overexpress_delta(model, g2u[g], floor, base)
            if d is not None:
                out[g] = round(d, 5)
        return out

    pos_scores = score([g for g in POS_GENES if g in g2u])
    neg_scores = score(neg_genes)
    pv = list(pos_scores.values()); nv = list(neg_scores.values())
    au = auroc(pv, nv)

    # stratified bootstrap CI
    rng2 = np.random.default_rng(SEED)
    pa, na = np.array(pv), np.array(nv)
    B = 20000
    boots = np.empty(B)
    for b in range(B):
        boots[b] = auroc(rng2.choice(pa, len(pa), replace=True),
                         rng2.choice(na, len(na), replace=True))
    ci = (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))

    dfree = distribution_free_stats(pv, nv)

    out = {
        "model": "eciML1515 (GECKO)", "product": PRODUCT, "growth_fraction": GF,
        "design": "blind held-out: external literature positives vs random negatives",
        "n_pos": len(pv), "n_neg": len(nv), "seed": SEED,
        "auroc_resource_aware": round(au, 4), "ci95_bootstrap": [round(ci[0], 4), round(ci[1], 4)],
        "distribution_free": dfree,
        "auroc_plain_fba": 0.5, "plain_note": "plain FBA gives Δ=0 for every overexpression (no signal)",
        "pos_scores": pos_scores, "neg_scores": neg_scores,
        "neg_scores_summary": {
            "n_zero": int(sum(1 for x in nv if abs(x) < 1e-6)),
            "n_positive": int(sum(1 for x in nv if x > 1e-6)),
            "max": round(max(nv), 5) if nv else None},
        "positives_literature": POS_GENES,
    }
    json.dump(out, open(os.path.join(D, "benchmark_blind_heldout.json"), "w"), indent=2)

    print(f"\n=== BLIND HELD-OUT (resource-aware ec-FBA) ===")
    print(f"  positives (literature, n={len(pv)}): " +
          ", ".join(f"{g}={v:+.2f}" for g, v in sorted(pos_scores.items(), key=lambda kv: -kv[1])))
    print(f"  negatives (random, n={len(nv)}): {out['neg_scores_summary']}")
    print(f"  AUROC = {au:.3f}  95% bootstrap CI [{ci[0]:.3f}, {ci[1]:.3f}]   (plain FBA = 0.5, no signal)")
    print(f"  Mann-Whitney U = {dfree['mannwhitney_U']:.0f}, one-sided p = {dfree['mannwhitney_p_onesided']:.2e}")
    print(f"  permutation p = {dfree['permutation_p']:.2e} ({dfree['permutation_note']}); "
          f"rank-biserial r = {dfree['rank_biserial_r']:.3f}; pairwise = {dfree['n_pairwise']}")
    print(f"  -> non-circular: negatives random, positives external. saved benchmark_blind_heldout.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
