"""Training engine for ``rna2psi``: fit one per-event PSI model from bulk RNA counts.

Pipeline (all NA-aware; PSI missing values stay NaN and are dropped per event):

1. Normalize counts exactly as inference does -- cp10k+log1p on the all-gene library
   size, then z-score per gene with training mu/sd. Feature genes are restricted to
   protein-coding Ensembl gene IDs.
2. For each event, rank candidate feature genes by **missing-value-aware Pearson
   correlation** (pairwise-complete, computed in vector/matrix mode over the event's
   non-NaN samples) and keep the top ``n_candidates``.
3. Fit ElasticNet (internal CV) on those candidates, keep the <=``max_features`` genes
   with the largest |coefficient|, then CV-pick ElasticNet vs Ridge on that reduced set
   (5-fold held-out Spearman). The winner refit on all complete samples is stored.

Produces a bundle dict consumable by ``PerEventPSIBundle`` (see ``_impute.py``) plus a
per-event performance table.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import ElasticNetCV, RidgeCV
from sklearn.model_selection import KFold, cross_val_predict


# ----------------------------------------------------------------- normalization
def cp10k_log1p_features(counts_gene_by_sample: pd.DataFrame,
                         feature_genes: list) -> tuple:
    """counts_gene_by_sample: genes (index) x samples (columns), raw counts.
    Library size is taken over ALL genes in the matrix (inference convention); features
    are then restricted to ``feature_genes``. Returns (Z samples x genes, mu, sd, genes)."""
    totals = counts_gene_by_sample.to_numpy(dtype=np.float64).sum(axis=0)     # per sample
    feat = [g for g in feature_genes if g in counts_gene_by_sample.index]
    sub = counts_gene_by_sample.loc[feat].to_numpy(dtype=np.float64)          # genes x samples
    Xn = np.log1p(sub / totals[None, :] * 1e4).T                              # samples x genes
    mu = Xn.mean(axis=0)
    sd = Xn.std(axis=0); sd[sd < 1e-8] = 1.0
    Z = (Xn - mu[None, :]) / sd[None, :]
    return Z, mu, sd, feat


# ------------------------------------------------------ NA-aware Pearson (vector mode)
def masked_pearson_topk(P: np.ndarray, Z: np.ndarray, k: int, chunk: int = 512) -> tuple:
    """Missing-value-aware Pearson correlation of every event against every gene, using
    only each event's non-NaN samples (pairwise-complete). Vectorized via matrix products.

    P : events x samples PSI, NaN = missing.   Z : samples x genes (no missing).
    Returns (top_idx events x k gene indices, top_r events x k signed r), padded with -1 /
    NaN where an event has fewer usable genes."""
    N, S = P.shape
    G = Z.shape[1]
    k = min(k, G)
    M = np.isfinite(P).astype(np.float64)          # events x samples present-mask
    P0 = np.where(np.isfinite(P), P, 0.0)          # events x samples, 0 where absent
    Z2 = Z * Z
    top_idx = np.full((N, k), -1, dtype=np.int64)
    top_r = np.full((N, k), np.nan, dtype=np.float64)
    for a in range(0, N, chunk):
        b = min(a + chunk, N)
        Mc, Pc = M[a:b], P0[a:b]
        n = Mc.sum(axis=1)                                     # (c,)
        sumP = Pc.sum(axis=1)
        sumP2 = (Pc * Pc).sum(axis=1)
        SPE = Pc @ Z                                           # c x G
        SE = Mc @ Z                                            # c x G
        SE2 = Mc @ Z2                                          # c x G
        with np.errstate(invalid="ignore", divide="ignore"):
            inv_n = 1.0 / n
            cov = SPE - sumP[:, None] * SE * inv_n[:, None]
            varP = sumP2 - sumP * sumP * inv_n                 # (c,)
            varE = SE2 - SE * SE * inv_n[:, None]              # c x G
            denom = np.sqrt(varP[:, None] * varE)
            r = np.where((denom > 0) & (n[:, None] >= 3), cov / denom, 0.0)
        r = np.clip(r, -1.0, 1.0)
        absr = np.abs(r)
        kk = min(k, G)
        part = np.argpartition(-absr, kth=kk - 1, axis=1)[:, :kk]
        rows = np.arange(b - a)[:, None]
        order = np.argsort(-absr[rows, part], axis=1)
        sel = part[rows, order]
        top_idx[a:b, :kk] = sel
        top_r[a:b, :kk] = r[rows, sel]
    return top_idx, top_r


# ----------------------------------------------------------------- per-event fit
def _cv_spearman(est, X, y, cv):
    if X.shape[1] == 0:
        return np.nan, np.nan
    try:
        pred = cross_val_predict(est, X, y, cv=cv)
    except Exception:
        return np.nan, np.nan
    if np.std(pred) < 1e-12:
        rho = 0.0
    else:
        rho = spearmanr(pred, y).statistic
        if not np.isfinite(rho):
            rho = 0.0
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2)) or 1.0
    return float(rho), float(1.0 - ss_res / ss_tot)


def fit_one_event(y_full: np.ndarray, Z: np.ndarray, cand_idx: np.ndarray, *,
                  max_features: int, cv: int, l1_ratios, random_state: int = 0) -> dict:
    """Fit one event. y_full: (n_samples,) PSI with NaN; Z: (n_samples, n_genes)
    z-scored features; cand_idx: candidate gene indices (into Z columns)."""
    mask = np.isfinite(y_full)
    n = int(mask.sum())
    out = {"n_train": n, "sel_idx": np.zeros(0, int), "coef": np.zeros(0),
           "intercept": np.nan, "estimator": None, "cv_spearman": np.nan, "cv_r2": np.nan}
    cand_idx = cand_idx[cand_idx >= 0]
    if n < max(cv + 1, 12) or cand_idx.size == 0:
        return out
    y = y_full[mask].astype(np.float64)
    Xc = Z[np.ix_(mask, cand_idx)]
    kf = KFold(n_splits=min(cv, n), shuffle=True, random_state=random_state)
    # 1) sparsify with ElasticNet on the candidate pool -> keep top |coef| genes
    en = ElasticNetCV(l1_ratio=list(l1_ratios), n_alphas=50, cv=kf,
                      max_iter=5000, random_state=random_state)
    try:
        en.fit(Xc, y)
    except Exception:
        return out
    coef = en.coef_
    nz = np.where(np.abs(coef) > 0)[0]
    take = nz[np.argsort(-np.abs(coef[nz]))][:max_features] if nz.size else \
        np.arange(min(max_features, cand_idx.size))          # fallback: top-|r| candidates
    sel = cand_idx[take]
    Xs = Z[np.ix_(mask, sel)]
    # 2) CV-pick ElasticNet vs Ridge on the reduced (<=max_features) feature set
    en2 = ElasticNetCV(l1_ratio=list(l1_ratios), n_alphas=50, cv=kf,
                       max_iter=5000, random_state=random_state)
    rg = RidgeCV(alphas=(0.1, 1.0, 10.0, 100.0))
    rho_en, r2_en = _cv_spearman(en2, Xs, y, kf)
    rho_rg, r2_rg = _cv_spearman(rg, Xs, y, kf)
    if (np.nan_to_num(rho_rg, nan=-1) > np.nan_to_num(rho_en, nan=-1)):
        best, name, rho, r2 = rg, "ridge", rho_rg, r2_rg
    else:
        best, name, rho, r2 = en2, "elasticnet", rho_en, r2_en
    try:
        best.fit(Xs, y)
    except Exception:
        return out
    out.update(sel_idx=sel.astype(int), coef=np.asarray(best.coef_, float).ravel(),
               intercept=float(best.intercept_), estimator=name,
               cv_spearman=rho, cv_r2=r2)
    return out


# ------------------------------------------------- honest (nested) reliability
def honest_oof_spearman(P: np.ndarray, Z: np.ndarray, *, k: int, cv: int = 5,
                        random_state: int = 0, min_train: int = 12) -> tuple:
    """Leakage-free out-of-fold reliability. For each outer fold, candidate genes are
    reselected by masked Pearson on that fold's TRAIN samples only, a RidgeCV is fit on
    <=k of them, and the held-out samples are predicted. Returns per-event (spearman, r2)
    of the assembled out-of-fold predictions vs observed PSI."""
    N, S = P.shape
    oof = np.full((N, S), np.nan)
    for tr, te in KFold(n_splits=cv, shuffle=True, random_state=random_state).split(np.arange(S)):
        ti, _ = masked_pearson_topk(P[:, tr], Z[tr], k=k)
        for i in range(N):
            m_tr = tr[np.isfinite(P[i, tr])]
            m_te = te[np.isfinite(P[i, te])]
            sel = ti[i][ti[i] >= 0]
            if len(m_tr) < min_train or sel.size == 0 or len(m_te) == 0:
                continue
            rg = RidgeCV(alphas=(0.1, 1.0, 10.0, 100.0)).fit(Z[np.ix_(m_tr, sel)], P[i, m_tr])
            oof[i, m_te] = rg.predict(Z[np.ix_(m_te, sel)])
    rho = np.full(N, np.nan); r2 = np.full(N, np.nan)
    for i in range(N):
        m = np.isfinite(P[i]) & np.isfinite(oof[i])
        if m.sum() >= min_train and np.std(oof[i, m]) > 1e-9:
            s = spearmanr(oof[i, m], P[i, m]).statistic
            rho[i] = s if np.isfinite(s) else np.nan
            ss_res = float(np.sum((P[i, m] - oof[i, m]) ** 2))
            ss_tot = float(np.sum((P[i, m] - P[i, m].mean()) ** 2)) or 1.0
            r2[i] = 1.0 - ss_res / ss_tot
    return rho, r2


# ----------------------------------------------------------------- driver
def train_rna2psi(counts: pd.DataFrame, psi: pd.DataFrame, *, feature_genes: list,
                  symbol_map: dict | None = None, n_candidates: int = 15,
                  max_features: int = 5, cv: int = 5,
                  l1_ratios=(0.2, 0.5, 0.8, 1.0), random_state: int = 0,
                  progress_every: int = 500) -> tuple:
    """counts: genes x samples raw counts. psi: events x samples PSI (NaN missing).
    Columns of both must already be the shared, identically-ordered sample set.
    Returns (bundle dict, per-event performance DataFrame)."""
    assert list(counts.columns) == list(psi.columns), "counts/psi sample columns must match"
    Z, mu, sd, genes = cp10k_log1p_features(counts, feature_genes)      # samples x genes
    events = list(psi.index)
    P = psi.to_numpy(dtype=np.float64)                                  # events x samples
    print(f"[rna2psi] features={len(genes)} genes  events={len(events)}  samples={Z.shape[0]}")
    print("[rna2psi] computing NA-aware Pearson candidate genes ...", flush=True)
    top_idx, top_r = masked_pearson_topk(P, Z, k=n_candidates)

    sel_idx, coef, intercept, rows = [], [], [], []
    for i, uid in enumerate(events):
        res = fit_one_event(P[i], Z, top_idx[i], max_features=max_features, cv=cv,
                            l1_ratios=l1_ratios, random_state=random_state)
        sel_idx.append(res["sel_idx"]); coef.append(res["coef"])
        intercept.append(res["intercept"])
        rows.append({"UID": uid, "n_train": res["n_train"], "estimator": res["estimator"],
                     "cv_spearman_insample": res["cv_spearman"], "cv_r2_insample": res["cv_r2"],
                     "n_features": int(res["sel_idx"].size),
                     "feature_genes": ";".join(genes[j] for j in res["sel_idx"]),
                     "top_abs_r": float(np.nanmax(np.abs(top_r[i]))) if np.isfinite(top_r[i]).any() else np.nan})
        if progress_every and (i + 1) % progress_every == 0:
            print(f"  {i+1}/{len(events)} events fit", flush=True)

    print("[rna2psi] computing leakage-free (nested) held-out reliability ...", flush=True)
    honest_rho, honest_r2 = honest_oof_spearman(P, Z, k=max_features, cv=cv,
                                                random_state=random_state)
    perf = pd.DataFrame(rows)
    perf.insert(3, "cv_spearman", honest_rho)          # primary, trustworthy reliability
    perf.insert(4, "cv_r2", honest_r2)
    meta = {
        "modality": "psi", "gene_id": "ensembl_gene_id",
        "normalization": "cp10k_log1p", "gene_filter": "protein_coding",
        "estimator": "per-event elasticnet-or-ridge (CV-pick)",
        "max_features": max_features, "n_candidates": n_candidates,
        "n_train_samples": int(Z.shape[0]), "cv_folds": cv,
        "heldout_median_spearman": float(np.nanmedian(perf["cv_spearman"])),
        "n_imputable_sp_gt_0p3": int((perf["cv_spearman"] > 0.3).sum()),
        "n_events_valid": int(perf["n_features"].gt(0).sum()),
        "symbol_map": {g: symbol_map.get(g, "") for g in genes} if symbol_map else {},
        "per_event": {
            "cv_spearman": dict(zip(perf["UID"], perf["cv_spearman"])),
            "cv_r2": dict(zip(perf["UID"], perf["cv_r2"])),
        },
    }
    bundle = {
        "X_columns": genes, "Y_columns": events,
        "mu": mu.astype(np.float64), "sd": sd.astype(np.float64),
        "sel_idx": [s.astype(np.int64) for s in sel_idx],
        "coef": [c.astype(np.float64) for c in coef],
        "intercept": np.asarray(intercept, np.float64),
        "metadata": meta,
    }
    return bundle, perf
