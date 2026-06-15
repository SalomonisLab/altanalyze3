#!/usr/bin/env python3.11
"""Compare normalization denominators for the per-edge model:
  (A) all-gene CP10k  (current; requires full transcriptome)
  (B) feature-gene CP10k (self-contained; only the 2,415 TF/target genes)
  (C) housekeeping-panel CP10k (self-contained: normalize by a depth proxy of
      stable, low-variance feature genes, preserving regulon-wide signal)
Metrics: leave-one-donor-out R^2, and in-silico perturbed-TF median rank.
Feature-gene representations are derived exactly from the stored all-gene-CP10k X.
"""
import sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, "/Users/saljh8/Documents/GitHub/altanalyze3")
from altanalyze3.components.rna2grn.evaluate import _Edges

OUT = "/Users/saljh8/Dropbox/Collaborations/Grimes/Human-GRN/July-2026-simple/rna2grn"
d = np.load(f"{OUT}/matched/dataset.npz", allow_pickle=True)
X = d["X"].astype(np.float64); Y = d["Y"].astype(np.float64); genes = d["genes"].astype(str)
etf = d["edge_tf"].astype(str); eg = d["edge_gene"].astype(str)
sample = d["sample"].astype(str); group = d["group"].astype(str)
ed = _Edges(genes, etf, eg)
Xf_all = X[:, ed.feat_cols]                      # all-gene CP10k+log1p, feature genes
# pseudo-CP10k counts of feature genes (proportional to raw counts / S_all)
pc = np.expm1(Xf_all)                            # = c_g/S_all*1e4
# (B) feature-gene CP10k: renormalize over the feature-gene sum
Xf_feat = np.log1p(pc / pc.sum(1, keepdims=True) * 1e4)
# (C) housekeeping CP10k: denominator = sum over low-CV feature genes (stable depth proxy)
mu = pc.mean(0); sd = pc.std(0); cv = sd / (mu + 1e-9)
hk = (mu > np.quantile(mu, 0.5)) & (cv < np.quantile(cv, 0.25))   # high-expressed, low-variance
Xf_hk = np.log1p(pc / pc[:, hk].sum(1, keepdims=True) * 1e4)
print(f"feature genes={Xf_all.shape[1]}  housekeeping-panel size={int(hk.sum())}")

direct = ~np.isin(group, ["Multiome_control", "TEA_control"])
samples = pd.unique(sample[direct]); GM = Y[direct].mean(0, keepdims=True)


def loso_r2(Xf):
    P = np.zeros((direct.sum(), ed.E)); idx = np.where(direct)[0]; pos = {r: i for i, r in enumerate(idx)}
    for s in samples:
        te = direct & (sample == s); tr = direct & ~(sample == s)
        co, st = ed.fit(Xf[tr], Y[tr]); Ph = ed.predict(co, st, Xf[te])
        for j, r in enumerate(np.where(te)[0]): P[pos[r]] = Ph[j]
    return 1 - ((Y[direct] - P) ** 2).sum() / ((Y[direct] - GM.repeat(direct.sum(), 0)) ** 2).sum()


def perturb_rank(Xf, boost=1.0, nbase=6, minreg=15, seed=0):
    rng = np.random.RandomState(seed)
    co, st = ed.fit(Xf, Y)
    tfs = ed.tf_list; tedge = {t: np.where(ed.tf_of_edge == i)[0] for i, t in enumerate(tfs)}
    # regulon target columns in feature space
    fpos = {g: i for i, g in enumerate(ed.feat)}
    tcols = {t: np.array([fpos[g] for g in pd.unique(eg[etf == t]) if g in fpos]) for t in tfs}
    testable = [t for t in tfs if (ed.tf_of_edge == tfs.index(t)).sum() >= minreg and len(tcols[t]) >= minreg]
    base = rng.choice(np.where(~direct)[0], size=nbase, replace=False)
    ranks = []
    for t in testable:
        rr = []
        for b in base:
            x0 = Xf[b:b + 1].copy(); xp = x0.copy(); xp[0, tcols[t]] += boost
            d0 = ed.predict(co, st, x0)[0]; d1 = ed.predict(co, st, xp)[0]
            de = d1 - d0
            sc = np.array([de[tedge[u]].mean() for u in tfs])
            rr.append(int(np.where(np.argsort(-sc) == tfs.index(t))[0][0]) + 1)
        ranks.append(np.median(rr))
    return float(np.median(ranks))


for name, Xf in [("A_all_gene", Xf_all), ("B_feature_gene", Xf_feat), ("C_housekeeping", Xf_hk)]:
    t0 = time.time(); r2 = loso_r2(Xf); pr = perturb_rank(Xf)
    print(f"{name:16s} LOSO_R2={r2:+.3f}   perturbed-TF median rank={pr:.0f}/217   ({time.time()-t0:.0f}s)", flush=True)
