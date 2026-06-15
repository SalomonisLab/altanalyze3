#!/usr/bin/env python3.11
"""DECISIVE test of the actual goal: detect a TF whose regulon is perturbed in a
NEW sample, where the perturbation is NOT represented in the reference.

In-silico: take a real pseudobulk, raise the expression of ONE TF's regulon
(targets) while holding the TF's own mRNA flat (mimics a mutation that changes
TF activity, not abundance). Predict GRN for baseline vs perturbed. The perturbed
TF's edges should rank top by increased predicted activity.

Retrieval (kNN) cannot do this (it returns the nearest existing reference); a
per-edge regulon regression can, because it extrapolates from regulon expression.

Also: per-sample breakdown of the held-out AML-vs-control differential, isolating
the singleton-mutation samples (SRSF2/DNMT3A/CSF3R) that have no cousin in the
reference -- the honest leave-signature-out cases.
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
tgt_idx = np.array([gene_pos[g] for g in edge_gene])
tf_idx = np.array([gene_pos[g] for g in edge_tf])
uniq_tf = sorted(set(edge_tf))
tf_pos = {t: i for i, t in enumerate(uniq_tf)}
tf_of_edge = np.array([tf_pos[t] for t in edge_tf])
tf_targets = {t: np.array([tgt_idx[e] for e in range(E) if edge_tf[e] == t]) for t in uniq_tf}
tf_edges = {t: np.array([e for e in range(E) if edge_tf[e] == t]) for t in uniq_tf}

is_ctrl = np.isin(group, ["Multiome_control", "TEA_control"])


def regulon_mean(Xmat):
    R = np.zeros((Xmat.shape[0], len(uniq_tf)))
    for t in uniq_tf:
        R[:, tf_pos[t]] = Xmat[:, tf_targets[t]].mean(1)
    return R


def features(Xmat, Rmat=None):
    if Rmat is None:
        Rmat = regulon_mean(Xmat)
    return [Xmat[:, tgt_idx], Xmat[:, tf_idx], Rmat[:, tf_of_edge]]


def fit_peredge(Xtr, Ytr, lam=1.0):
    feats = features(Xtr)
    n = Xtr.shape[0]
    stats = []; cols = [np.ones((n, E))]
    for A in feats:
        mu = A.mean(0); sd = A.std(0); sd[sd < 1e-8] = 1.0
        stats.append((mu, sd)); cols.append((A - mu) / sd)
    Dsg = np.stack(cols, 2)
    AtA = np.einsum("nef,neg->efg", Dsg, Dsg)
    Aty = np.einsum("nef,ne->ef", Dsg, Ytr)
    reg = lam * np.eye(len(feats) + 1)[None].repeat(E, 0); reg[:, 0, 0] = 0
    coefs = np.linalg.solve(AtA + reg, Aty)
    return coefs, stats


def predict_peredge(coefs, stats, Xmat):
    feats = features(Xmat)
    cols = [np.ones((Xmat.shape[0], E))]
    for A, (mu, sd) in zip(feats, stats):
        cols.append((A - mu) / sd)
    return np.einsum("nef,ef->ef" if False else "nef,ef->ne", np.stack(cols, 2), coefs)


# ===== in-silico perturbation =====
def run_perturbation(boost=1.0, min_regulon=15, n_baselines=6, seed=0):
    rng = np.random.RandomState(seed)
    # fit per-edge model on ALL data (the shipped reference)
    coefs, stats = fit_peredge(X, Y)
    # kNN setup
    feat_idx = np.array([gene_pos[g] for g in sorted(set(edge_tf) | set(edge_gene))])
    Xf = X[:, feat_idx]; mu = Xf.mean(0); sd = Xf.std(0); sd[sd < 1e-8] = 1
    Xfs = (Xf - mu) / sd
    nn = NearestNeighbors(n_neighbors=15, metric="cosine").fit(Xfs)

    def knn_pred(Xq):
        Qf = (Xq[:, feat_idx] - mu) / sd
        dist, ind = nn.kneighbors(Qf); w = 1 / (dist + 1e-6); w /= w.sum(1, keepdims=True)
        return np.einsum("ij,ijk->ik", w, Y[ind])

    testable = [t for t in uniq_tf if len(tf_targets[t]) >= min_regulon]
    # baseline pseudobulks: random reference rows (use healthy controls + a few AML)
    base_rows = rng.choice(np.where(is_ctrl)[0], size=min(n_baselines, is_ctrl.sum()), replace=False)
    recs = []
    for t in testable:
        ranks_pe, ranks_knn = [], []
        for r in base_rows:
            x0 = X[r:r + 1].copy()
            xp = x0.copy()
            xp[0, tf_targets[t]] += boost          # raise regulon expression
            # (TF mRNA held flat: we do NOT touch tf_idx of this TF)
            for model, store in [("pe", ranks_pe), ("knn", ranks_knn)]:
                if model == "pe":
                    p0 = predict_peredge(coefs, stats, x0)[0]
                    p1 = predict_peredge(coefs, stats, xp)[0]
                else:
                    p0 = knn_pred(x0)[0]; p1 = knn_pred(xp)[0]
                dedge = p1 - p0
                tf_score = np.array([dedge[tf_edges[u]].mean() for u in uniq_tf])
                order = np.argsort(-tf_score)
                rank = int(np.where(order == tf_pos[t])[0][0]) + 1
                store.append(rank)
        recs.append(dict(TF=t, regulon=len(tf_targets[t]),
                         pe_rank_med=np.median(ranks_pe), knn_rank_med=np.median(ranks_knn)))
    df = pd.DataFrame(recs)
    log(f"\n=== IN-SILICO REGULON PERTURBATION (boost=+{boost} log1p, {len(testable)} TFs, "
        f"{len(base_rows)} baselines) ===")
    log(f"  median rank of the PERTURBED TF among {len(uniq_tf)} TFs (1 = perfect):")
    log(f"    per-edge regression : median={df.pe_rank_med.median():.0f}  "
        f"top1={100*(df.pe_rank_med==1).mean():.0f}%  top5={100*(df.pe_rank_med<=5).mean():.0f}%")
    log(f"    kNN retrieval       : median={df.knn_rank_med.median():.0f}  "
        f"top1={100*(df.knn_rank_med==1).mean():.0f}%  top5={100*(df.knn_rank_med<=5).mean():.0f}%")
    df.to_csv(f"{OUT}/benchmark/perturbation_tf_ranks.csv", index=False)
    return df


# ===== per-sample held-out AML-vs-control differential (isolate singletons) =====
def per_sample_differential():
    ctrl_rows = np.where(group == "Multiome_control")[0]
    ctrl_cs = {cell_state[r]: r for r in ctrl_rows}
    ann = pd.read_csv(f"{OUT}/matched/meta.csv")
    # map sample -> mutation via h5ad already done earlier; reuse simple heuristic from names
    aml_samples = pd.unique(sample[~is_ctrl])
    feat_idx = np.array([gene_pos[g] for g in sorted(set(edge_tf) | set(edge_gene))])
    Xf = X[:, feat_idx]

    def pearson(a, b):
        a = a - a.mean(); b = b - b.mean(); de = np.sqrt((a**2).sum()*(b**2).sum())
        return float((a*b).sum()/de) if de > 1e-12 else np.nan

    rows = []
    for s in aml_samples:
        te = sample == s; tr = ~te
        rws = [r for r in np.where(te)[0] if cell_state[r] in ctrl_cs]
        if not rws:
            continue
        crows = [ctrl_cs[cell_state[r]] for r in rws]
        # kNN
        mu = Xf[tr].mean(0); sd = Xf[tr].std(0); sd[sd < 1e-8] = 1
        Xs = (Xf - mu) / sd
        nn = NearestNeighbors(n_neighbors=15, metric="cosine").fit(Xs[tr])
        def kp(rws_):
            di, ind = nn.kneighbors(Xs[rws_]); w = 1/(di+1e-6); w/=w.sum(1,keepdims=True)
            return np.einsum("ij,ijk->ik", w, Y[tr][ind])
        Ya, Yc = kp(rws), kp(crows)
        # per-edge
        coefs, stats = fit_peredge(X[tr], Y[tr])
        Pa = predict_peredge(coefs, stats, X[np.array(rws)])
        Pc = predict_peredge(coefs, stats, X[np.array(crows)])
        dk = [pearson(Ya[i]-Yc[i], Y[rws[i]]-Y[crows[i]]) for i in range(len(rws))]
        dp = [pearson(Pa[i]-Pc[i], Y[rws[i]]-Y[crows[i]]) for i in range(len(rws))]
        rows.append(dict(sample=s, n_cs=len(rws), knn_delta_r=np.nanmedian(dk),
                         peredge_delta_r=np.nanmedian(dp)))
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/benchmark/differential_per_sample.csv", index=False)
    log("\n=== per-sample held-out AML-vs-control differential (singletons = no cousin) ===")
    singles = {"AML-7_CCHMC": "SRSF2(single)", "AML-12_CCHMC": "DNMT3A(single)", "AML-14_CCHMC": "CSF3R(single)"}
    for _, r in df.sort_values("peredge_delta_r").iterrows():
        tag = singles.get(r["sample"], "")
        log(f"  {r['sample']:22s} {tag:14s} n_cs={int(r['n_cs']):2d} "
            f"knn_delta_r={r['knn_delta_r']:.3f}  peredge_delta_r={r['peredge_delta_r']:.3f}")
    return df


if __name__ == "__main__":
    t0 = time.time()
    per_sample_differential()
    for b in [0.5, 1.0, 2.0]:
        run_perturbation(boost=b)
    log(f"\ntotal {time.time()-t0:.0f}s")
