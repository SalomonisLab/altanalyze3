#!/usr/bin/env python3.11
"""DECISIVE experiment: can a per-edge regulon-aware regression recover the
true AML-vs-control DIFFERENTIAL edge activities from RNA, where kNN / cell-state
mean cannot (they predict ~zero differential by construction)?

Per-edge model for edge (TF_i -> gene_j):
    activity[p,e] ~ b0 + b1*expr(gene_j) + b2*expr(TF_i) + b3*regulon_mean_i(p)
fit by ridge, vectorized across all 7,486 edges. This EXTRAPOLATES: if a regulon's
target expression rises (TF mutated/active but mRNA flat), predicted edge activity
rises beyond any reference profile.

Validation (user spec): held-out AML sample, per cell state, differential vs the
Multiome_WT control. The imputed differential must recover the true differential,
especially the INCREASED edges/TFs.
"""
import time
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"


def log(*a): print(*a, flush=True)


d = np.load(f"{OUT}/matched/dataset.npz", allow_pickle=True)
X, Y = d["X"].astype(np.float64), d["Y"].astype(np.float64)
genes = d["genes"].astype(str)
edge_tf = d["edge_tf"].astype(str); edge_gene = d["edge_gene"].astype(str)
sample = d["sample"].astype(str); cell_state = d["cell_state"].astype(str)
group = d["group"].astype(str)
gene_pos = {g: i for i, g in enumerate(genes)}
E = len(edge_tf)

# edge -> gene/TF column indices
tgt_idx = np.array([gene_pos[g] for g in edge_gene])
tf_idx = np.array([gene_pos[g] for g in edge_tf])
# regulon membership: TF -> list of target gene columns
tf_to_targets = {}
for e in range(E):
    tf_to_targets.setdefault(edge_tf[e], []).append(tgt_idx[e])
uniq_tf = list(tf_to_targets)
tf_of_edge = np.array([uniq_tf.index(t) for t in edge_tf])


def build_features(Xmat):
    """Return per-edge feature stack (n, E, F): target, TF, regulon_mean."""
    n = Xmat.shape[0]
    tgt = Xmat[:, tgt_idx]                                  # (n, E)
    tf = Xmat[:, tf_idx]                                    # (n, E)
    reg_by_tf = np.zeros((n, len(uniq_tf)))
    for ti, t in enumerate(uniq_tf):
        cols = tf_to_targets[t]
        reg_by_tf[:, ti] = Xmat[:, cols].mean(1)
    reg = reg_by_tf[:, tf_of_edge]                          # (n, E)
    return tgt, tf, reg


TGT, TF, REG = build_features(X)
masks = {"control": np.isin(group, ["Multiome_control", "TEA_control"])}
is_ctrl = masks["control"]
is_aml = ~is_ctrl

# control reference: Multiome_WT GRN + RNA features per cell state
ctrl_rows = np.where(group == "Multiome_control")[0]
ctrl_cs = {cell_state[r]: r for r in ctrl_rows}


def zfit(A):
    mu = A.mean(0); sd = A.std(0); sd[sd < 1e-8] = 1.0
    return mu, sd


def per_edge_ridge_fit(feat_list, Ytr, lam=1.0):
    """Vectorized per-edge ridge. feat_list: list of (n,E) arrays. Returns
    coefs (E, F+1) and the per-feature z-stats for prediction."""
    n = Ytr.shape[0]; F = len(feat_list)
    stats = []
    cols = [np.ones((n, E))]
    for A in feat_list:
        mu = A.mean(0); sd = A.std(0); sd[sd < 1e-8] = 1.0
        stats.append((mu, sd)); cols.append((A - mu) / sd)
    Dsg = np.stack(cols, axis=2)                           # (n, E, F+1)
    # normal equations per edge
    AtA = np.einsum("nef,neg->efg", Dsg, Dsg)              # (E, F+1, F+1)
    Aty = np.einsum("nef,ne->ef", Dsg, Ytr)                # (E, F+1)
    reg = lam * np.eye(F + 1)[None]; reg[:, 0, 0] = 0.0    # don't penalize intercept
    coefs = np.linalg.solve(AtA + reg, Aty)               # (E, F+1)
    return coefs, stats


def per_edge_predict(coefs, stats, feat_list):
    n = feat_list[0].shape[0]
    cols = [np.ones((n, E))]
    for A, (mu, sd) in zip(feat_list, stats):
        cols.append((A - mu) / sd)
    Dsg = np.stack(cols, axis=2)
    return np.einsum("nef,ef->ne", Dsg, coefs)


FEATSETS = {
    "target_only": lambda i: [TGT[i]],
    "target+TF": lambda i: [TGT[i], TF[i]],
    "target+regulon": lambda i: [TGT[i], REG[i]],
    "target+TF+regulon": lambda i: [TGT[i], TF[i], REG[i]],
}


def col_pearson(A, B):
    Az = A - A.mean(0, keepdims=True); Bz = B - B.mean(0, keepdims=True)
    num = (Az * Bz).sum(0); den = np.sqrt((Az ** 2).sum(0) * (Bz ** 2).sum(0))
    out = np.full(A.shape[1], np.nan); ok = den > 1e-12; out[ok] = num[ok] / den[ok]; return out


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    de = np.sqrt((a ** 2).sum() * (b ** 2).sum())
    return float((a * b).sum() / de) if de > 1e-12 else np.nan


def auroc_increased(delta_true, delta_pred, thr):
    """AUROC for predicting which edges are truly increased (> thr)."""
    pos = delta_true > thr
    if pos.sum() < 5 or (~pos).sum() < 5:
        return np.nan
    order = np.argsort(delta_pred)
    ranks = np.empty_like(order, dtype=float); ranks[order] = np.arange(len(order))
    # Mann-Whitney U -> AUROC
    rp = ranks[pos].sum(); n1 = pos.sum(); n0 = (~pos).sum()
    auc = (rp - n1 * (n1 - 1) / 2) / (n1 * n0)
    return float(auc)


aml_samples = pd.unique(sample[is_aml])
log(f"AML samples (LOSO): {len(aml_samples)}; edges={E}; controls={is_ctrl.sum()}")

# ---------- run per-edge models + baselines on the DIFFERENTIAL metric ----------
def eval_model_differential(predict_aml_ctrl):
    """predict_aml_ctrl(train_mask, aml_rows, ctrl_rows_for_those) -> (Yhat_aml, Yhat_ctrl).
    Returns dataframe of per-(sample,cellstate) differential metrics."""
    recs = []
    for s in aml_samples:
        te = (sample == s)
        tr = ~te
        rows = np.where(te)[0]
        # only cell states with a control match
        rows = [r for r in rows if cell_state[r] in ctrl_cs]
        if not rows:
            continue
        crows = [ctrl_cs[cell_state[r]] for r in rows]
        Yhat_aml, Yhat_ctrl = predict_aml_ctrl(tr, rows, crows)
        for k, (r, cr) in enumerate(zip(rows, crows)):
            dt = Y[r] - Y[cr]                  # true differential vs control
            dp = Yhat_aml[k] - Yhat_ctrl[k]    # predicted differential
            recs.append(dict(sample=s, cell_state=cell_state[r],
                             delta_pearson=pearson(dp, dt),
                             auroc_up=auroc_increased(dt, dp, 0.02)))
    return pd.DataFrame(recs)


# baseline predictors
def make_knn(tr_mask):
    feat_idx = np.array([gene_pos[g] for g in sorted(set(edge_tf) | set(edge_gene))])
    Xf = X[:, feat_idx]
    mu, sd = zfit(Xf[tr_mask]);
    Xs = (Xf - mu) / sd
    nn = NearestNeighbors(n_neighbors=15, metric="cosine").fit(Xs[tr_mask])
    Ytr = Y[tr_mask]
    def pred(rows):
        dist, ind = nn.kneighbors(Xs[rows]); w = 1/((dist)+1e-6); w/=w.sum(1,keepdims=True)
        return np.einsum("ij,ijk->ik", w, Ytr[ind])
    return pred


def knn_aml_ctrl(tr, rows, crows):
    p = make_knn(tr); return p(rows), p(crows)


def csmean_aml_ctrl(tr, rows, crows):
    lk = {c: Y[tr & (cell_state == c)].mean(0) for c in pd.unique(cell_state[tr])}
    gm = Y[tr].mean(0)
    A = np.vstack([lk.get(cell_state[r], gm) for r in rows])
    C = np.vstack([lk.get(cell_state[r], gm) for r in crows])
    return A, C


def make_peredge(featfn, lam=1.0):
    def f(tr, rows, crows):
        coefs, stats = per_edge_ridge_fit(featfn(tr), Y[tr], lam=lam)
        Yh_a = per_edge_predict(coefs, stats, featfn(np.array(rows)))
        Yh_c = per_edge_predict(coefs, stats, featfn(np.array(crows)))
        return Yh_a, Yh_c
    return f


results = {}
for name, fn in [("cellstate_mean", csmean_aml_ctrl), ("knn", knn_aml_ctrl)]:
    t0 = time.time(); df = eval_model_differential(fn); results[name] = df
    log(f"{name:22s} delta_pearson med={df.delta_pearson.median():.3f}  "
        f"AUROC_up med={df.auroc_up.median():.3f}  ({time.time()-t0:.0f}s)")
for name, featfn in FEATSETS.items():
    t0 = time.time(); df = eval_model_differential(make_peredge(featfn)); results["peredge_" + name] = df
    log(f"{'peredge_'+name:22s} delta_pearson med={df.delta_pearson.median():.3f}  "
        f"AUROC_up med={df.auroc_up.median():.3f}  ({time.time()-t0:.0f}s)")

summary = pd.DataFrame({k: dict(delta_pearson_med=v.delta_pearson.median(),
                                delta_pearson_mean=v.delta_pearson.mean(),
                                auroc_up_med=v.auroc_up.median())
                        for k, v in results.items()}).T
summary.to_csv(f"{OUT}/benchmark/differential_summary.csv")
log("\n" + summary.to_string())
log("\nwrote " + f"{OUT}/benchmark/differential_summary.csv")
