"""Training and held-out evaluation for per-molecule RNA -> mass-spec imputation bundles.

Ports the validated AML recipe so any matched RNA-count + molecule table can be trained
with it. Source of truth for every rule below:

    Human-MS-impute/code/evaluate_imputation.py   (features, CV, per-target metrics)
    Human-MS-impute/code/build_imputation_model.py (final all-case fit)
    Human-MS-impute/code/make_altanalyze3_bundles.py (bundle dict)

With the AML settings (``splitter="kfold"``, 5 splits, seed 0) this module reproduces the
shipped AML model; ``validation/reproduce_aml_bundle.py`` checks it value by value.

Recipe:
  X  = log1p(CP10k) of protein-coding gene counts summed per symbol; the CP10k library
       size is the ALL-gene column sum; genes kept when TPM >= min_tpm in >= min_samples.
  per molecule: rank its prior genes by |Spearman| (training rows), take the top nfeat,
       fill the rest from the genome-wide |Spearman| ranking over the expressed pool
       (mean log-expression above the median, computed over all rows), fit Ridge(alpha)
       on z-scored genes (training mu/sd), skipping rows where the molecule is NaN.
"""
from __future__ import annotations

import gzip
import pickle
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold, KFold, LeaveOneGroupOut

SPLITTERS = ("kfold", "group_kfold", "leave_one_group_out")


def tpm_from_counts(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    """TPM of a genes x samples count matrix; ``lengths`` in bp, indexed like ``counts``."""
    L = lengths.reindex(counts.index).astype(float)
    if L.isna().any() or (L <= 0).any():
        raise ValueError(f"{int(L.isna().sum())} genes lack a length; {int((L <= 0).sum())} have length <= 0")
    rate = counts.div(L, axis=0)
    return rate / rate.sum(0) * 1e6


def cp10k_log1p_symbol(counts: pd.DataFrame, gene_annotation: pd.DataFrame, tpm: pd.DataFrame,
                       coding_only: bool = True, min_tpm: float = 5.0, min_samples: int = 3) -> pd.DataFrame:
    """Samples x symbols log1p(CP10k) matrix (evaluate_imputation.cp10k_log1p_symbol).

    counts, tpm: genes (Ensembl) x samples. gene_annotation: index Ensembl, columns
    ``gene_type`` and ``gene_name``.
    """
    C = counts
    tot = C.sum(0); tot[tot == 0] = 1.0
    gt = gene_annotation.reindex(C.index)
    coding = (gt["gene_type"].astype(str) == "protein_coding").values if coding_only else np.ones(len(C), bool)
    sym = gt["gene_name"].astype(str).values
    Csym = C[coding].groupby(sym[coding]).sum()
    T = tpm.reindex(C.index)
    Tsym = T[coding].groupby(sym[coding]).sum()
    keep = (((Tsym >= min_tpm).sum(axis=1) >= min_samples).reindex(Csym.index).fillna(False))
    Csym = Csym[keep.values]
    X = np.log1p(Csym / tot * 1e4)
    return X.T


def col_spearman(Ytr: np.ndarray, Xtr: np.ndarray) -> np.ndarray:
    """Targets x genes rank correlation; target NaNs take the column median (selection only)."""
    Yf = pd.DataFrame(Ytr); Yf = Yf.fillna(Yf.median()).fillna(0.0)
    def rank(A): return pd.DataFrame(A).rank().values
    Yr = rank(Yf.values); Xr = rank(Xtr)
    Yz = (Yr - Yr.mean(0)) / (Yr.std(0) + 1e-9)
    Xz = (Xr - Xr.mean(0)) / (Xr.std(0) + 1e-9)
    return (Yz.T @ Xz) / Yr.shape[0]


def expressed_pool(Xv: np.ndarray) -> np.ndarray:
    expr = Xv.mean(0)
    return np.where(expr > np.quantile(expr, 0.5))[0]


def select_features(acorr_t: np.ndarray, prior_idx: Sequence[int], pool: np.ndarray, nfeat: int) -> List[int]:
    chosen = sorted(prior_idx, key=lambda j: -acorr_t[j])[:nfeat]
    if len(chosen) < nfeat:
        cs = set(chosen)
        for j in pool[np.argsort(-acorr_t[pool])]:
            j = int(j)
            if j not in cs:
                chosen.append(j); cs.add(j)
                if len(chosen) == nfeat:
                    break
    return chosen


def _prior_index(prior: Mapping[str, Sequence[str]], targets: Sequence[str], gidx: Mapping[str, int]) -> List[List[int]]:
    return [[gidx[g] for g in prior.get(t, []) if g in gidx] for t in targets]


def _splits(n: int, splitter: str, n_splits: int, seed: int, groups: Optional[Sequence]):
    if splitter not in SPLITTERS:
        raise ValueError(f"splitter must be one of {SPLITTERS}")
    idx = np.arange(n)
    if splitter == "kfold":
        if groups is not None:
            raise ValueError("splitter='kfold' ignores groups; pass a group splitter or groups=None")
        return list(KFold(n_splits, shuffle=True, random_state=seed).split(idx))
    if groups is None or len(groups) != n:
        raise ValueError(f"{splitter} needs one group label per row")
    g = np.asarray(groups)
    if splitter == "group_kfold":
        return list(GroupKFold(n_splits).split(idx, groups=g))
    return list(LeaveOneGroupOut().split(idx, groups=g))


def cross_validate(X: pd.DataFrame, Y: pd.DataFrame, prior: Mapping[str, Sequence[str]], *, nfeat: int,
                   alpha: float, splitter: str = "kfold", n_splits: int = 5, seed: int = 0,
                   groups: Optional[Sequence] = None, min_fit: int = 10) -> Dict[str, object]:
    """Out-of-fold predictions (evaluate_imputation.main). X: rows x genes; Y: rows x molecules."""
    if list(X.index) != list(Y.index):
        raise ValueError("X and Y rows differ in identity or order")
    genes = list(X.columns); gidx = {g: i for i, g in enumerate(genes)}
    targets = list(Y.columns)
    Xv = X.values; Yv = Y.values
    pool = expressed_pool(Xv)
    pri_idx = _prior_index(prior, targets, gidx)
    folds = _splits(len(X), splitter, n_splits, seed, groups)
    P = np.full_like(Yv, np.nan, dtype=float)
    nprior_used = np.zeros(len(targets))
    for tr, te in folds:
        Xtr, Xte = Xv[tr], Xv[te]; Ytr = Yv[tr]
        mu = Xtr.mean(0); sd = Xtr.std(0); sd[sd < 1e-8] = 1
        Xtrs = (Xtr - mu) / sd; Xtes = (Xte - mu) / sd
        acorr = np.abs(np.nan_to_num(col_spearman(Ytr, Xtr)))
        for ti in range(len(targets)):
            chosen = select_features(acorr[ti], pri_idx[ti], pool, nfeat)
            nprior_used[ti] += min(len(pri_idx[ti]), nfeat)
            if not chosen:
                continue
            ytr = Ytr[:, ti]; ok = ~np.isnan(ytr)
            if ok.sum() < min_fit:
                continue
            m = Ridge(alpha=alpha).fit(Xtrs[ok][:, chosen], ytr[ok])
            P[te, ti] = m.predict(Xtes[:, chosen])
    return dict(oof=P, nprior_used=nprior_used, n_folds=len(folds), folds=folds, pool=pool)


def per_target_metrics(X: pd.DataFrame, Y: pd.DataFrame, prior: Mapping[str, Sequence[str]], cv: Mapping,
                       *, nfeat: int, alpha: float, min_eval: int = 10, min_insample: int = 20) -> pd.DataFrame:
    """Per-molecule held-out metrics plus the all-row top-10 gene record (evaluate_imputation.main)."""
    genes = list(X.columns); gidx = {g: i for i, g in enumerate(genes)}
    targets = list(Y.columns)
    Xv = X.values; Yv = Y.values; P = cv["oof"]; pool = cv["pool"]
    mu = Xv.mean(0); sd = Xv.std(0); sd[sd < 1e-8] = 1; Xall = (Xv - mu) / sd
    corr = np.nan_to_num(col_spearman(Yv, Xv)); acorr = np.abs(corr)
    pri_idx = _prior_index(prior, targets, gidx)
    gm = np.nanmean(Yv, 0)
    rows = []
    for ti, t in enumerate(targets):
        pri = pri_idx[ti]
        chosen = select_features(acorr[ti], pri, pool, nfeat)
        top10 = sorted(chosen, key=lambda j: -acorr[ti, j])[:10]
        feature_genes = ";".join(f"{genes[j]}({'prior' if j in pri else 'corr'}:{corr[ti, j]:+.2f})" for j in top10)
        y = Yv[:, ti]; okf = ~np.isnan(y)
        insample = np.nan
        if okf.sum() >= min_insample and chosen:
            m = Ridge(alpha=alpha).fit(Xall[okf][:, chosen], y[okf])
            insample = spearmanr(y[okf], m.predict(Xall[okf][:, chosen])).correlation
        p = P[:, ti]; ok = ~np.isnan(p) & ~np.isnan(y)
        if ok.sum() < min_eval or np.nanstd(y[ok]) < 1e-9:
            continue
        rows.append(dict(target=t, n=int(ok.sum()), n_prior=round(cv["nprior_used"][ti] / cv["n_folds"], 1),
                         pearson=np.corrcoef(y[ok], p[ok])[0, 1], spearman=spearmanr(y[ok], p[ok]).correlation,
                         r2=1 - np.nansum((y[ok] - p[ok]) ** 2) / np.nansum((y[ok] - gm[ti]) ** 2),
                         insample_spearman=insample, feature_genes=feature_genes))
    return pd.DataFrame(rows)


def fit_final(X: pd.DataFrame, Y: pd.DataFrame, prior: Mapping[str, Sequence[str]], *, nfeat: int,
              alpha: float, min_fit: int = 12) -> Dict[str, object]:
    """All-row fit (build_imputation_model.build); returns the model dict before packaging."""
    if list(X.index) != list(Y.index):
        raise ValueError("X and Y rows differ in identity or order")
    genes = list(X.columns); gidx = {g: i for i, g in enumerate(genes)}
    Xv = X.values; Yv = Y.values; targets = list(Y.columns)
    mu = Xv.mean(0); sd = Xv.std(0); sd[sd < 1e-8] = 1.0
    Xstd = (Xv - mu) / sd
    pool = expressed_pool(Xv)
    acorr = np.abs(np.nan_to_num(col_spearman(Yv, Xv)))
    pri_idx = _prior_index(prior, targets, gidx)
    sel_idx, coef, intercept, insample = [], [], [], []
    for ti in range(len(targets)):
        chosen = select_features(acorr[ti], pri_idx[ti], pool, nfeat)
        y = Yv[:, ti]; ok = ~np.isnan(y)
        if ok.sum() < min_fit or not chosen:
            sel_idx.append(np.array([], int)); coef.append(np.array([])); intercept.append(np.nan); insample.append(np.nan)
            continue
        Xc = Xstd[ok][:, chosen]
        m = Ridge(alpha=alpha).fit(Xc, y[ok])
        sel_idx.append(np.asarray(chosen, int)); coef.append(m.coef_.astype(np.float32)); intercept.append(float(m.intercept_))
        insample.append(spearmanr(y[ok], m.predict(Xc)).correlation)
    return dict(nfeat=nfeat, alpha=alpha, universe_genes=genes, mu=mu.astype(np.float32), sd=sd.astype(np.float32),
                targets=targets, sel_idx=sel_idx, coef=coef, intercept=np.array(intercept),
                insample_spearman=np.array(insample), n_cases=int(Xv.shape[0]))


def _fold_predictions(Xv, Yv, pri_idx, pool, tr, te, nfeats, alphas, min_fit):
    """Predictions for rows ``te`` from a fit on rows ``tr`` at every (nfeat, alpha); one selection pass.

    A smaller nfeat selection is a prefix of the largest one (select_features orders prior genes,
    then the genome-wide fill), so selection runs once at max(nfeats) and is sliced.
    """
    Xtr, Xte = Xv[tr], Xv[te]; Ytr = Yv[tr]
    mu = Xtr.mean(0); sd = Xtr.std(0); sd[sd < 1e-8] = 1
    Xtrs = (Xtr - mu) / sd; Xtes = (Xte - mu) / sd
    acorr = np.abs(np.nan_to_num(col_spearman(Ytr, Xtr)))
    kmax = max(nfeats)
    out = {(k, a): np.full((len(te), Yv.shape[1]), np.nan) for k in nfeats for a in alphas}
    for ti in range(Yv.shape[1]):
        order = select_features(acorr[ti], pri_idx[ti], pool, kmax)
        ytr = Ytr[:, ti]; ok = ~np.isnan(ytr)
        if not order or ok.sum() < min_fit:
            continue
        for k in nfeats:
            chosen = order[:k]
            for a in alphas:
                m = Ridge(alpha=a).fit(Xtrs[ok][:, chosen], ytr[ok])
                out[(k, a)][:, ti] = m.predict(Xtes[:, chosen])
    return out


def median_target_spearman(Y: np.ndarray, P: np.ndarray, min_eval: int = 10) -> float:
    vals = []
    for ti in range(Y.shape[1]):
        y = Y[:, ti]; p = P[:, ti]; ok = ~np.isnan(y) & ~np.isnan(p)
        if ok.sum() >= min_eval and np.nanstd(y[ok]) > 1e-9 and np.nanstd(p[ok]) > 1e-12:
            vals.append(spearmanr(y[ok], p[ok]).correlation)
    return float(np.median(vals)) if vals else np.nan


def nested_cross_validate(universes: Mapping[str, pd.DataFrame], Y: pd.DataFrame,
                          prior: Mapping[str, Sequence[str]], groups: Sequence, *,
                          nfeat_grid: Sequence, alpha_grid: Sequence[float], min_fit: int = 10,
                          min_eval: int = 10, n_jobs: int = 1) -> Dict[str, object]:
    """Leave-one-group-out outer CV; inner leave-one-group-out picks one (universe, nfeat, alpha).

    universes: name -> rows x genes matrix (same rows as Y). nfeat_grid entries are ints or "all"
    (every gene in that universe). The inner score is the median per-target Spearman of the inner
    out-of-fold predictions; ties go to the earlier grid entry. Returns outer OOF predictions,
    the chosen setting per outer fold and every inner score.
    """
    rows = list(Y.index)
    for name, X in universes.items():
        if list(X.index) != rows:
            raise ValueError(f"universe {name}: rows differ from Y in identity or order")
    g = np.asarray(groups)
    if len(g) != len(rows):
        raise ValueError("one group label per row is required")
    targets = list(Y.columns); Yv = Y.values
    prep = {}
    for name, X in universes.items():
        genes = list(X.columns); gidx = {x: i for i, x in enumerate(genes)}
        ks = sorted({len(genes) if k == "all" else min(int(k), len(genes)) for k in nfeat_grid})
        prep[name] = dict(Xv=X.values, pool=expressed_pool(X.values), ks=ks,
                          pri=_prior_index(prior, targets, gidx), n_genes=len(genes))
    grid = [(u, k, float(a)) for u in universes for k in prep[u]["ks"] for a in alpha_grid]
    outer = list(LeaveOneGroupOut().split(np.arange(len(rows)), groups=g))
    args = [(fi, otr, ote, prep, Yv, g, grid, alpha_grid, min_fit, min_eval) for fi, (otr, ote) in enumerate(outer)]
    if n_jobs == 1:
        results = [_outer_fold(*a) for a in args]
    else:
        from joblib import Parallel, delayed
        results = Parallel(n_jobs=n_jobs)(delayed(_outer_fold)(*a) for a in args)
    P = np.full(Yv.shape, np.nan)
    chosen, inner_scores = [], []
    for (fi, otr, ote, *_), (pred, ch, sc) in zip(args, results):
        P[ote] = pred
        chosen.append(ch); inner_scores.extend(sc)
    return dict(oof=P, chosen=pd.DataFrame(chosen), inner_scores=pd.DataFrame(inner_scores), grid=grid)


def _outer_fold(fi, otr, ote, prep, Yv, g, grid, alpha_grid, min_fit, min_eval):
    inner = list(LeaveOneGroupOut().split(otr, groups=g[otr]))
    oof = {c: np.full((len(otr), Yv.shape[1]), np.nan) for c in grid}
    for itr, ite in inner:
        for u, pp in prep.items():
            preds = _fold_predictions(pp["Xv"][otr], Yv[otr], pp["pri"], pp["pool"], itr, ite,
                                      pp["ks"], alpha_grid, min_fit)
            for (k, a), pr in preds.items():
                oof[(u, k, float(a))][ite] = pr
    scores = {c: median_target_spearman(Yv[otr], oof[c], min_eval) for c in grid}
    best = max(grid, key=lambda c: (np.nan_to_num(scores[c], nan=-np.inf), -grid.index(c)))
    sc = [dict(outer_fold=fi, held_out_group=str(g[ote][0]), universe=c[0], nfeat=c[1], alpha=c[2],
               inner_median_spearman=scores[c], chosen=c == best) for c in grid]
    u, k, a = best
    pp = prep[u]
    pred = _fold_predictions(pp["Xv"], Yv, pp["pri"], pp["pool"], otr, ote, [k], [a], min_fit)[(k, a)]
    ch = dict(outer_fold=fi, held_out_group=str(g[ote][0]), n_test_rows=len(ote), universe=u, nfeat=k,
              alpha=a, inner_median_spearman=scores[best])
    return pred, ch, sc


def to_bundle(model: Mapping, per_target: pd.DataFrame, metadata: Mapping) -> Dict[str, object]:
    """Bundle dict read by ``_impute.PerTargetImputeBundle`` (make_altanalyze3_bundles.build)."""
    pt = per_target.set_index("target")
    meta = dict(metadata)
    meta.setdefault("nfeat", int(model["nfeat"])); meta.setdefault("alpha", float(model["alpha"]))
    meta.setdefault("n_train_cases", int(model["n_cases"]))
    meta["heldout_median_spearman"] = round(float(pt["spearman"].median()), 3)
    meta["n_imputable_sp_gt_0p3"] = int((pt["spearman"] > 0.3).sum())
    meta["per_target"] = dict(
        heldout_spearman={t: float(v) for t, v in pt["spearman"].items()},
        heldout_r2={t: float(v) for t, v in pt["r2"].items()},
        feature_genes={t: str(v) for t, v in pt["feature_genes"].items()})
    return dict(X_columns=list(model["universe_genes"]), Y_columns=list(model["targets"]),
                mu=np.asarray(model["mu"], np.float32), sd=np.asarray(model["sd"], np.float32),
                sel_idx=[np.asarray(x, np.int32) for x in model["sel_idx"]],
                coef=[np.asarray(x, np.float32) for x in model["coef"]],
                intercept=np.asarray(model["intercept"], np.float32), metadata=meta)


def write_bundle(bundle: Mapping, path) -> None:
    with gzip.open(path, "wb") as fh:
        pickle.dump(dict(bundle), fh, protocol=4)
