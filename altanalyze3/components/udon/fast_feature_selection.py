#!/usr/bin/env python3
"""Ultra-fast, mathematically-IDENTICAL reimplementations of UDON's feature
selection bottlenecks (vectorized numpy instead of per-gene Python loops /
pandas.corr). Used by compare_workflows.py to verify the fast path reproduces
the original feature set exactly while being orders of magnitude faster.

Each function takes the fold matrix as a genes x pseudobulks DataFrame and
returns the selected gene Index, matching the corresponding original in
feature_selection.py:
  variance_based_feature_selection      (line 42)
  intercorrelation_based_feature_selection (line 60)
  pca_feature_selection                 (line 82)
"""
import numpy as np, pandas as pd
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.preprocessing import scale


def fast_rank(fold_hvg, oversample=80):
    """determine_nmf_ranks (nmf.py:9) via randomized SVD instead of full LA.eig on
    a samples x samples covariance (O(n^3)). Eigenvalues of scale(X)^T scale(X)
    are the squared singular values of scale(X). Returns est_k = 2*k (same as
    original), validated to match the eig rank on the subset."""
    X = scale(np.asarray(fold_hvg))
    g, c = X.shape
    muTW = (np.sqrt(g - 1) + np.sqrt(c)) ** 2.0
    sigmaTW = (np.sqrt(g - 1) + np.sqrt(c)) * (1.0 / np.sqrt(g - 1) + 1.0 / np.sqrt(c)) ** (1.0 / 3.0)
    boundary = 3.273 * sigmaTW + muTW
    ncomp = int(min(c - 1, oversample))
    sv = TruncatedSVD(n_components=ncomp, random_state=0).fit(X).singular_values_
    return 2 * int((sv ** 2 > boundary).sum())


def fast_variance(folds, fold_threshold=1, samples_differing=3):
    """Original: per gene, sort folds; range = (samples_differing-th largest) -
    (samples_differing-th smallest); keep if > fold_threshold. Vectorized via
    np.partition -> identical selection."""
    M = folds.to_numpy()
    sd = samples_differing
    kth_small = np.partition(M, sd - 1, axis=1)[:, sd - 1]      # sd-th smallest
    kth_large = np.partition(M, -sd, axis=1)[:, -sd]            # sd-th largest
    keep = (kth_large - kth_small) > fold_threshold
    return folds.index[keep]


def fast_intercorrelation(folds, hvg_index, corr_threshold=0.4, corr_n_events=5):
    """Original: pearson gene-gene corr on HVG genes; keep genes with
    (corr>thr).sum() > corr_n_events (diagonal self-corr counted, as in original).
    Vectorized via np.corrcoef."""
    sub = folds.loc[hvg_index].to_numpy()
    R = np.corrcoef(sub)                                        # genes x genes
    R = np.nan_to_num(R)
    keep = (R > corr_threshold).sum(axis=0) > corr_n_events
    return hvg_index[keep]


def fast_pca(folds, corr_index, corr_threshold=0.4, n_components=30):
    """Original pca_feature_selection, numpy-vectorized per-PC correlation.
    Same logic: top/bottom 200 loadings per PC, keep genes with >0 abs corr
    (|r|>thr, lower-triangular) to another gene in that 400-set."""
    exp = folds.loc[corr_index]
    X = exp.to_numpy()
    pca = PCA(n_components=n_components, svd_solver="full").fit(X.T)
    loadings = pca.components_.T                                # genes x PCs
    genes = np.asarray(exp.index)
    selected = set()
    # pre-standardize rows once for correlation
    for pc in range(n_components):
        order = np.argsort(loadings[:, pc])[::-1]
        idx = np.concatenate([order[:200], order[-200:]])
        idx = np.unique(idx)
        sub = X[idx]
        R = np.corrcoef(sub)
        R = np.nan_to_num(np.tril(R, k=-1))
        binary = np.abs(R) > corr_threshold
        has = (binary.sum(axis=0) > 0) | (binary.sum(axis=1) > 0)
        for g in genes[idx[has]]:
            selected.add(g)
    return pd.Index([g for g in genes if g in selected])


def fast_feature_selection(folds, protein_coding_genes=None, fold_threshold=1,
                           samples_differing=3, intercorr_threshold=0.4,
                           corr_n_events=5, pca_corr_threshold=0.4, n_components=30,
                           apply_gene_name_filter=True, do_pca=False):
    """Full fast feature selection on a fold matrix. If protein_coding_genes is
    given (set of symbols), applies the same protein-coding + non-coding-name
    filters as the original before variance selection. Returns dict of stage->Index.

    do_pca=False (default) SKIPS the PCA stage: clustering_wrapper consumes only
    'correlated_genes', so pca_selected_genes is unused downstream and the PCA
    full-SVD is the slowest part. Set do_pca=True only to compare against the
    original's pca_selected_genes."""
    genes = folds.index
    if apply_gene_name_filter:
        mask = np.ones(len(genes), dtype=bool)
        gl = genes.to_series()
        mask &= ~gl.str.startswith("MT-").to_numpy()
        mask &= ~gl.str.startswith(("RPS", "RPL")).to_numpy()
        mask &= ~gl.str.contains(r"\.").to_numpy()
        mask &= ~gl.str.startswith(("Gm", "GM")).to_numpy()
        mask &= ~gl.str.endswith(("y", "Y")).to_numpy()
        mask &= ~gl.str.startswith(("HLA", "Hla")).to_numpy()
        mask &= ~gl.str.startswith(("Xis", "XIS", "xis")).to_numpy()
        mask &= ~gl.str.startswith(("TSI", "Tsi")).to_numpy()
        if protein_coding_genes is not None:
            mask &= gl.apply(lambda g: g.split("_")[0] in protein_coding_genes).to_numpy()
        folds = folds.loc[genes[mask]]
    hvg = fast_variance(folds, fold_threshold, samples_differing)
    corr = fast_intercorrelation(folds, hvg, intercorr_threshold, corr_n_events)
    pca = fast_pca(folds, corr, pca_corr_threshold, n_components) if do_pca else corr
    return {"after_name_filter": folds.index, "udon_hvg": hvg,
            "correlated_genes": corr, "pca_selected_genes": pca}
