#!/usr/bin/env python3.11
"""rna2grn LOSO benchmark (vectorized). Leave-one-SAMPLE-out.

Honest metrics: skill OVER the cell-state-mean baseline, per-edge across-sample
Pearson, within/global R2 — not just raw profile correlation.
"""
import time
import numpy as np
import pandas as pd
from functools import partial
from sklearn.linear_model import Ridge
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"


def log(*a):
    print(*a, flush=True)


d = np.load(f"{OUT}/matched/dataset.npz", allow_pickle=True)
X, Y = d["X"], d["Y"]
genes = d["genes"].astype(str)
edge_tf = d["edge_tf"].astype(str)
edge_gene = d["edge_gene"].astype(str)
sample = d["sample"].astype(str)
cell_state = d["cell_state"].astype(str)
group = d["group"].astype(str)

direct = ~np.isin(group, ["Multiome_control", "TEA_control"])
X, Y, sample, cell_state = X[direct], Y[direct], sample[direct], cell_state[direct]
log(f"LOSO rows: {X.shape[0]}  edges: {Y.shape[1]}  samples: {len(set(sample))}")

gene_pos = {g: i for i, g in enumerate(genes)}
feat_genes = sorted(set(edge_tf) | set(edge_gene))
feat_idx = np.array([gene_pos[g] for g in feat_genes])
Xf = X[:, feat_idx].astype(np.float64)
log(f"feature genes (TF u target): {len(feat_genes)}")

samples = pd.unique(sample)
n, E = Y.shape
GM = Y.mean(0, keepdims=True)
CS_LK = {c: Y[cell_state == c].mean(0) for c in pd.unique(cell_state)}


def vec_col_pearson(A, B):
    """per-column Pearson between A and B (n x E), vectorized."""
    Az = A - A.mean(0, keepdims=True)
    Bz = B - B.mean(0, keepdims=True)
    num = (Az * Bz).sum(0)
    den = np.sqrt((Az ** 2).sum(0) * (Bz ** 2).sum(0))
    out = np.full(A.shape[1], np.nan)
    ok = den > 1e-12
    out[ok] = num[ok] / den[ok]
    return out


def vec_row_pearson(A, B):
    Az = A - A.mean(1, keepdims=True)
    Bz = B - B.mean(1, keepdims=True)
    num = (Az * Bz).sum(1)
    den = np.sqrt((Az ** 2).sum(1) * (Bz ** 2).sum(1))
    out = np.full(A.shape[0], np.nan)
    ok = den > 1e-12
    out[ok] = num[ok] / den[ok]
    return out


def metrics(Ytrue, Ypred):
    pr = vec_row_pearson(Ytrue, Ypred)
    ed = vec_col_pearson(Ytrue, Ypred)
    ss_res = ((Ytrue - Ypred) ** 2).sum()
    ss_tot = ((Ytrue - GM) ** 2).sum()
    CSb = np.vstack([CS_LK[c] for c in cell_state])
    ss_cs = ((Ytrue - CSb) ** 2).sum()
    return dict(profile_r_med=np.nanmedian(pr), profile_r_mean=np.nanmean(pr),
                edge_r_med=np.nanmedian(ed), edge_r_mean=np.nanmean(ed),
                edge_r_pos_frac=float(np.nanmean(ed > 0)),
                r2_global=1 - ss_res / ss_tot, skill_vs_cs_mean=1 - ss_res / ss_cs)


def standardize(Xtr, Xte):
    mu = Xtr.mean(0); sd = Xtr.std(0); sd[sd < 1e-8] = 1.0
    return (Xtr - mu) / sd, (Xte - mu) / sd


def loso_predict(fit_predict, label=""):
    P = np.zeros_like(Y)
    t0 = time.time()
    for s in samples:
        te = sample == s; tr = ~te
        P[te] = fit_predict(Xf[tr], Y[tr], Xf[te], sample[tr], cell_state[tr])
    log(f"   {label} LOSO done in {time.time()-t0:.0f}s")
    return P


def fp_global_mean(Xtr, Ytr, Xte, st, ct):
    return np.repeat(Ytr.mean(0, keepdims=True), len(Xte), axis=0)


def loso_cs_mean():
    P = np.zeros_like(Y)
    for s in samples:
        te = sample == s; tr = ~te
        gm = Y[tr].mean(0)
        lk = {c: Y[tr][cell_state[tr] == c].mean(0) for c in pd.unique(cell_state[tr])}
        P[te] = np.vstack([lk.get(c, gm) for c in cell_state[te]])
    return P


def fp_ridge(Xtr, Ytr, Xte, st, ct, alpha=10.0):
    Xs, Xq = standardize(Xtr, Xte)
    m = Ridge(alpha=alpha).fit(Xs, Ytr)
    return m.predict(Xq)


def fp_svd_ridge(Xtr, Ytr, Xte, st, ct, k=30, alpha=5.0):
    Xs, Xq = standardize(Xtr, Xte)
    svd = TruncatedSVD(n_components=min(k, Xs.shape[0] - 1), random_state=0).fit(Xs)
    Ztr, Zte = svd.transform(Xs), svd.transform(Xq)
    m = Ridge(alpha=alpha).fit(Ztr, Ytr)
    return m.predict(Zte)


def fp_knn(Xtr, Ytr, Xte, st, ct, k=10):
    Xs, Xq = standardize(Xtr, Xte)
    nn = NearestNeighbors(n_neighbors=min(k, len(Xs))).fit(Xs)
    dist, ind = nn.kneighbors(Xq)
    w = 1.0 / (dist + 1e-6); w /= w.sum(1, keepdims=True)
    return np.einsum("ij,ijk->ik", w, Ytr[ind])


MODELS = {
    "global_mean": ("special_gm", None),
    "cellstate_mean": ("special_cs", None),
    "ridge_a1": ("fp", partial(fp_ridge, alpha=1.0)),
    "ridge_a10": ("fp", partial(fp_ridge, alpha=10.0)),
    "ridge_a30": ("fp", partial(fp_ridge, alpha=30.0)),
    "ridge_a100": ("fp", partial(fp_ridge, alpha=100.0)),
    "svd30_ridge": ("fp", partial(fp_svd_ridge, k=30, alpha=5.0)),
    "svd50_ridge": ("fp", partial(fp_svd_ridge, k=50, alpha=5.0)),
    "knn_5": ("fp", partial(fp_knn, k=5)),
    "knn_10": ("fp", partial(fp_knn, k=10)),
    "knn_15": ("fp", partial(fp_knn, k=15)),
}

rows = []
np.save(f"{OUT}/benchmark/_y_true.npy", Y)
for name, (kind, fn) in MODELS.items():
    if kind == "special_gm":
        P = loso_predict(fp_global_mean, name)
    elif kind == "special_cs":
        P = loso_cs_mean()
    else:
        P = loso_predict(fn, name)
    m = metrics(Y, P); m["model"] = name
    rows.append(m)
    np.save(f"{OUT}/benchmark/pred_{name}.npy", P.astype(np.float32))
    log(f"{name:16s} prof_r={m['profile_r_med']:.3f} edge_r_med={m['edge_r_med']:.3f} "
        f"edge_r+={m['edge_r_pos_frac']:.2f} R2={m['r2_global']:+.3f} skillVsCS={m['skill_vs_cs_mean']:+.3f}")

res = pd.DataFrame(rows)[["model", "profile_r_med", "profile_r_mean", "edge_r_med",
                          "edge_r_mean", "edge_r_pos_frac", "r2_global", "skill_vs_cs_mean"]]
res.to_csv(f"{OUT}/benchmark/loso_summary.csv", index=False)
log("\n" + res.to_string(index=False))
log("\nwrote " + f"{OUT}/benchmark/loso_summary.csv")
