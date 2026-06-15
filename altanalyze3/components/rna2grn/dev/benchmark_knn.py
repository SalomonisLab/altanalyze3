#!/usr/bin/env python3.11
"""Optimize the kNN GRN predictor: distance metric, k, feature set, and whether
adding the control aggregates (Multiome/TEA) to the reference helps. Also a
cross-protocol transfer test. Leave-one-SAMPLE-out on direct rows.
"""
import time
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"


def log(*a): print(*a, flush=True)


d = np.load(f"{OUT}/matched/dataset.npz", allow_pickle=True)
Xall, Yall = d["X"], d["Y"]
genes = d["genes"].astype(str)
edge_tf = d["edge_tf"].astype(str); edge_gene = d["edge_gene"].astype(str)
sample = d["sample"].astype(str); cell_state = d["cell_state"].astype(str)
group = d["group"].astype(str)

gene_pos = {g: i for i, g in enumerate(genes)}
tf_tg = sorted(set(edge_tf) | set(edge_gene))
tf_tg_idx = np.array([gene_pos[g] for g in tf_tg])

direct = ~np.isin(group, ["Multiome_control", "TEA_control"])
ctrl = ~direct

Y = Yall[direct]; samp = sample[direct]; cs = cell_state[direct]
samples = pd.unique(samp)
GM = Y.mean(0); CS_LK_honest = None  # honest baseline computed per-fold below

# honest LOSO cellstate-mean baseline ss
def cs_mean_pred():
    P = np.zeros_like(Y)
    for s in samples:
        te = samp == s; tr = ~te
        gm = Y[tr].mean(0)
        lk = {c: Y[tr][cs[tr] == c].mean(0) for c in pd.unique(cs[tr])}
        P[te] = np.vstack([lk.get(c, gm) for c in cs[te]])
    return P
CSpred = cs_mean_pred()
ss_cs = ((Y - CSpred) ** 2).sum()
Yr_cs = Y - CSpred  # residual vs honest cellstate mean


def col_pearson(A, B):
    Az = A - A.mean(0, keepdims=True); Bz = B - B.mean(0, keepdims=True)
    num = (Az * Bz).sum(0); den = np.sqrt((Az ** 2).sum(0) * (Bz ** 2).sum(0))
    out = np.full(A.shape[1], np.nan); ok = den > 1e-12; out[ok] = num[ok] / den[ok]; return out


def feat_space(name):
    if name == "tftg":
        return tf_tg_idx
    if name == "tftg_hvg":
        # tf/target + top-1500 highly variable genes (by variance across direct rows)
        v = Xall[direct].var(0)
        hvg = np.argsort(-v)[:1500]
        return np.array(sorted(set(tf_tg_idx.tolist()) | set(hvg.tolist())))
    raise ValueError(name)


def knn_predict(Xtr, Ytr, Xte, k, metric):
    mu = Xtr.mean(0); sd = Xtr.std(0); sd[sd < 1e-8] = 1.0
    Xs = (Xtr - mu) / sd; Xq = (Xte - mu) / sd
    if metric == "cosine":
        nn = NearestNeighbors(n_neighbors=min(k, len(Xs)), metric="cosine").fit(Xs)
    else:
        nn = NearestNeighbors(n_neighbors=min(k, len(Xs)), metric="euclidean").fit(Xs)
    dist, ind = nn.kneighbors(Xq)
    w = 1.0 / (dist + 1e-6); w /= w.sum(1, keepdims=True)
    return np.einsum("ij,ijk->ik", w, Ytr[ind])


def evaluate(featname, k, metric, add_controls):
    fidx = feat_space(featname)
    Xf = Xall[:, fidx]
    Xd = Xf[direct]; Xc = Xf[ctrl]; Yc = Yall[ctrl]
    P = np.zeros_like(Y)
    for s in samples:
        te = samp == s; tr = ~te
        Xtr = Xd[tr]; Ytr = Y[tr]
        if add_controls:
            Xtr = np.vstack([Xtr, Xc]); Ytr = np.vstack([Ytr, Yc])
        P[te] = knn_predict(Xtr, Ytr, Xd[te], k, metric)
    ss_res = ((Y - P) ** 2).sum()
    skill = 1 - ss_res / ss_cs
    Pr = P - CSpred
    wedge = col_pearson(Yr_cs, Pr)
    edge = col_pearson(Y, P)
    return dict(feat=featname, k=k, metric=metric, ctrl=add_controls,
                skill_vs_CS=skill, within_cs_edge_r=np.nanmedian(wedge),
                edge_r=np.nanmedian(edge), r2=1 - ss_res / ((Y - GM) ** 2).sum())


log("=== kNN optimization (LOSO on direct rows) ===")
rows = []
for feat in ["tftg", "tftg_hvg"]:
    for metric in ["euclidean", "cosine"]:
        for k in [5, 10, 15, 20, 30]:
            for ctl in [False, True]:
                t0 = time.time()
                r = evaluate(feat, k, metric, ctl)
                rows.append(r)
                log(f"  feat={feat:9s} metric={metric:9s} k={k:2d} ctrl={str(ctl):5s} "
                    f"skill={r['skill_vs_CS']:+.3f} within_cs_r={r['within_cs_edge_r']:.3f} "
                    f"edge_r={r['edge_r']:.3f} r2={r['r2']:.3f} ({time.time()-t0:.0f}s)")
res = pd.DataFrame(rows).sort_values("skill_vs_CS", ascending=False)
res.to_csv(f"{OUT}/benchmark/knn_optimization.csv", index=False)
log("\nTOP 8:\n" + res.head(8).to_string(index=False))

# ---- cross-protocol transfer ----
log("\n=== CROSS-PROTOCOL TRANSFER (kNN tftg euclidean k=10) ===")
fidx = tf_tg_idx; Xf = Xall[:, fidx]
prot = np.where(group == "Multiome_control", "Multiome",
        np.where(group == "TEA_control", "TEA", "AML3p"))
def ss(a, b): return ((a - b) ** 2).sum()
for train_p, test_p in [("AML3p", "Multiome"), ("AML3p", "TEA"),
                        ("Multiome", "AML3p"), ("Multiome", "TEA")]:
    tri = prot == train_p; tei = prot == test_p
    if tei.sum() == 0 or tri.sum() == 0: continue
    P = knn_predict(Xf[tri], Yall[tri], Xf[tei], 10, "euclidean")
    Yt = Yall[tei]
    # baseline: global mean of train
    gm = Yall[tri].mean(0, keepdims=True)
    r2 = 1 - ss(Yt, P) / ss(Yt, Yt.mean(0, keepdims=True))
    er = np.nanmedian(col_pearson(Yt, P))
    log(f"  train={train_p:9s} -> test={test_p:9s} n_test={tei.sum():3d} edge_r={er:.3f} r2={r2:+.3f}")
log("\nwrote " + f"{OUT}/benchmark/knn_optimization.csv")
