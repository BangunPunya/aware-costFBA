import json, sys, numpy as np
from scipy import stats

rng = np.random.default_rng(42)
_REPORT = sys.argv[1] if len(sys.argv) > 1 else "code/m3_benchmark_report.json"
m3 = json.load(open(_REPORT))
print(f"[report] {_REPORT}")
labels, scores, wt = m3["labels"], m3["scores"], m3["wt_baseline"]
SCALAR = m3["score_scalar"]   # 'prod_at_near_max'

sub = [c for c in labels if labels[c] in ("tp_pathway", "fp_medium", "fp_thermo")]
y = np.array([1 if labels[c] == "tp_pathway" else 0 for c in sub])

def product_of(cid, arm):
    # product = the wt key contained in the candidate id (longest match wins, e.g. EX_lac__D_e)
    cand = [k for k in wt[arm] if k in cid]
    return max(cand, key=len) if cand else None

def lift_vec(arm):
    out = []
    for c in sub:
        p = product_of(c, arm)
        base = wt[arm].get(p, 0.0)
        out.append(scores[arm][c][SCALAR] - base)
    return np.array(out, float)

def auroc(sc, yy):
    pos = sc[yy == 1]; neg = sc[yy == 0]
    if len(pos) == 0 or len(neg) == 0: return np.nan
    return np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg])

print(f"n_pos={int(y.sum())}  n_neg={int((y==0).sum())}")
for arm in ("polos", "full", "full_anaer"):
    a = auroc(lift_vec(arm), y)
    rep = {"polos": m3['Klaim2']['auroc_tp_vs_fp_polos'],
           "full": m3['Klaim2']['auroc_tp_vs_fp_full'],
           "full_anaer": m3['Klaim2']['auroc_tp_vs_fp_anaer']}[arm]
    print(f"  {arm:11s} AUROC={a:.4f}  (report {rep})  {'MATCH' if abs(a-rep)<1e-6 else 'MISMATCH'}")

sc_plain = lift_vec("polos"); sc_full = lift_vec("full")
auc_p = auroc(sc_plain, y); auc_f = auroc(sc_full, y)

# stratified bootstrap (resample positives & negatives w/ replacement)
B = 20000
pos_idx = np.where(y == 1)[0]; neg_idx = np.where(y == 0)[0]
def boot(sc):
    out = np.empty(B)
    for b in range(B):
        pi = rng.choice(pos_idx, len(pos_idx), replace=True)
        ni = rng.choice(neg_idx, len(neg_idx), replace=True)
        out[b] = auroc(np.r_[sc[pi], sc[ni]], np.r_[np.ones(len(pi)), np.zeros(len(ni))])
    return out
def boot_diff():
    out = np.empty(B)
    for b in range(B):
        pi = rng.choice(pos_idx, len(pos_idx), replace=True)
        ni = rng.choice(neg_idx, len(neg_idx), replace=True)
        yy = np.r_[np.ones(len(pi)), np.zeros(len(ni))]
        out[b] = auroc(np.r_[sc_full[pi], sc_full[ni]], yy) - auroc(np.r_[sc_plain[pi], sc_plain[ni]], yy)
    return out
bp, bf, bd = boot(sc_plain), boot(sc_full), boot_diff()
ci = lambda x: (np.nanpercentile(x, 2.5), np.nanpercentile(x, 97.5))
print(f"\nBootstrap 95% CI (B={B}, seed=42):")
print(f"  plain AUROC = {auc_p:.3f}  CI [{ci(bp)[0]:.3f}, {ci(bp)[1]:.3f}]")
print(f"  full  AUROC = {auc_f:.3f}  CI [{ci(bf)[0]:.3f}, {ci(bf)[1]:.3f}]")
print(f"  diff (full-plain) = {auc_f-auc_p:+.3f}  CI [{ci(bd)[0]:+.3f}, {ci(bd)[1]:+.3f}]")
print(f"  P(diff>0) bootstrap = {np.mean(bd>0):.3f}")

# DeLong (Sun & Xu 2014) for correlated AUROCs.
def midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x); T = np.zeros(N); i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]: j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1; i = j
    T2 = np.empty(N); T2[J] = T; return T2
def delong(scA, scB):
    m = int(y.sum()); n = int((y == 0).sum())
    P = np.vstack([scA[y == 1], scB[y == 1]]); Q = np.vstack([scA[y == 0], scB[y == 0]])
    Z = np.hstack([P, Q]); k = 2
    tx = np.array([midrank(P[r]) for r in range(k)])
    ty = np.array([midrank(Q[r]) for r in range(k)])
    tz = np.array([midrank(Z[r]) for r in range(k)])
    aucs = (tz[:, :m].sum(1) / m - (m + 1) / 2) / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1 - (tz[:, m:] - ty) / m
    S = np.cov(v01) / m + np.cov(v10) / n
    return aucs, S
aucs, S = delong(sc_full, sc_plain)
L = np.array([1.0, -1.0]); var = L @ S @ L; diff = aucs[0] - aucs[1]
print(f"\nDeLong (correlated, full vs plain):")
if var <= 0:
    print(f"  diff={diff:+.4f}  variance<=0 (degenerate at n_pos=4); z/p undefined.")
else:
    se = np.sqrt(var); z = diff / se; p = 2 * (1 - stats.norm.cdf(abs(z)))
    print(f"  AUC full={aucs[0]:.4f} plain={aucs[1]:.4f} diff={diff:+.4f}")
    print(f"  SE={se:.4f}  z={z:.3f}  p={p:.4g}  95%CI [{diff-1.96*se:+.3f}, {diff+1.96*se:+.3f}]")
print(f"\nn_pos=4: CIs wide, low power; results are indicative, not inferential.")
