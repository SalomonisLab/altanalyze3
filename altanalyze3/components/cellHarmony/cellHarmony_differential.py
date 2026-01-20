#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
cellHarmony_differential.py

Accessory pipeline for cellHarmony h5ad outputs.

Functions:
  - load_and_merge_covariates(...)
  - compute_pseudobulk_per_population(...)
  - run_de_for_comparisons(...)
  - summarize_global_local_coreg(...)
  - build_fixed_order_heatmap(...)
  - run_goelite_for_clusters(...)
  - write_differentials_only_h5ad(...)
  - main() with CLI

Notes:
  - Requires: scanpy, anndata, numpy, pandas, matplotlib
  - Optional: rpy2 (not used by default). We default to Scanpy methods.
"""

import os
import sys
import textwrap
import math
import logging
import argparse
import datetime
import collections
from collections import Counter
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sps
from tqdm import tqdm

logging.getLogger("fontTools").setLevel(logging.WARNING)
logging.getLogger("fontTools.subset").setLevel(logging.WARNING)
from altanalyze3.components.visualization import NetPerspective
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings("ignore", message="Variable names are not unique. To make them unique, call `.var_names_make_unique`.")
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning, message=".*invalid value encountered in log2.*")

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'

diagnostic_report = False
# ------------------------------- Utilities -------------------------------- #

class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
        return len(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()

    @property
    def encoding(self):
        for stream in self.streams:
            if hasattr(stream, "encoding"):
                return stream.encoding
        return "utf-8"

def _assert_no_multiindex(df, name):
    if isinstance(df.index, pd.MultiIndex):
        raise ValueError("{} has a MultiIndex index which is unsupported. Please flatten indices.".format(name))
    if isinstance(df.columns, pd.MultiIndex):
        raise ValueError("{} has a MultiIndex columns which is unsupported. Please flatten columns.".format(name))

def _ensure_numeric_matrix(x):
    if sps.issparse(x):
        return x
    if isinstance(x, np.ndarray):
        return x
    # Convert to array
    try:
        arr = np.asarray(x)
    except Exception as e:
        raise ValueError("Failed to convert matrix to ndarray: {}".format(e))
    return arr

def _suppress_fonttools_logs():
    """Reduce noisy fontTools subset logging emitted during vector export."""
    for name in ("fontTools", "fontTools.subset"):
        logger = logging.getLogger(name)
        logger.setLevel(logging.ERROR)
        logger.propagate = False
        logger.disabled = True

def _is_categorical_series(s):
    return pd.api.types.is_categorical_dtype(s)

def _ordered_categories_from_obs(obs, col):
    # Use categorical order if present, else first-appearance order
    s = obs[col].astype(str)
    if _is_categorical_series(obs[col]):
        cats = list(obs[col].cat.categories)
        return cats
    seen = []
    for v in s:
        if v not in seen:
            seen.append(v)
    return seen

def _validate_lineage_order_candidates(candidate_order, population_labels, population_col):
    """
    Ensure lineage_order pulled from adata.uns matches the population labels used for DE.
    Returns the sanitized candidate_order list if valid, otherwise raises ValueError.
    """
    raw_list = [str(x) for x in list(candidate_order)]
    expected = [str(x) for x in population_labels]

    print(f"[DEBUG] Validating lineage_order for '{population_col}': candidate={raw_list}")
    print(f"[DEBUG] Expected population labels: {expected}")

    filtered: List[str] = []
    seen = set()
    for entry in raw_list:
        if entry in seen:
            continue
        seen.add(entry)
        filtered.append(entry)

    missing = [name for name in expected if name not in filtered]
    unexpected = [name for name in filtered if name not in expected]

    if missing:
        print(f"[INFO] lineage_order missing {len(missing)} populations (e.g. {missing[:5]})")
    if unexpected:
        print(f"[INFO] lineage_order contains {len(unexpected)} populations not present in DE (e.g. {unexpected[:5]}); ignoring them.")

    filtered = [name for name in filtered if name in expected]

    return filtered

def _extract_lineage_order(adata, population_col):
    """
    Attempt to extract and validate lineage_order from adata.uns against population_col.
    Returns a list if validation succeeds, otherwise None (with warning already printed).
    """
    if not hasattr(adata, "uns") or "lineage_order" not in adata.uns:
        return None

    try:
        populations = _ordered_categories_from_obs(adata.obs, population_col)
    except KeyError:
        print(f"[WARN] Cannot validate lineage_order: obs missing '{population_col}'.")
        return None

    try:
        return _validate_lineage_order_candidates(adata.uns["lineage_order"], populations, population_col)
    except ValueError as exc:
        print(f"[WARN] Ignoring lineage_order for '{population_col}': {exc}")
        return None
    except Exception as exc:
        print(f"[WARN] Failed to parse lineage_order for '{population_col}': {exc}")
        return None

def _pattern_to_label(pattern, populations):
    """
    Convert a list of -1/0/1 values to a descriptive string, e.g.:
    [0, -1, 1] -> 'Pop2__down--Pop3__up'
    """
    labels = []
    for val, pop_name in zip(pattern, populations):
        if val == 1:
            labels.append(f"{pop_name}__up")
        elif val == -1:
            labels.append(f"{pop_name}__down")
    return "--".join(labels) if labels else "none"

def _ln_to_log2(v):
    # Scanpy logfoldchanges are in natural log when applicable; convert to log2 safely.
    return np.array(v, dtype=float) / np.log(2.0)

def _bh_fdr(pvals):
    # Benjamini-Hochberg FDR for a 1D array-like
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = np.empty(n, dtype=float)
    ranked[order] = np.arange(1, n + 1, dtype=float)
    fdr = p * (n / ranked)
    # enforce monotonicity
    fdr_sorted = np.minimum.accumulate(fdr[order][::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(fdr_sorted, 1.0)
    return out

# --------------------------- Step 1: covariates ---------------------------- #

def load_and_merge_covariates(h5ad_path, covariate_tsv, library_col, sample_col, covariate_col):
    """
    Load an AnnData object and merge external covariate metadata by library/sample.
    Ensures .obs_names remain string-based after merge (fixes AnnData indexing errors).
    """
    import os, sys, pandas as pd, anndata as ad

    # --- Validate inputs ---
    if not os.path.exists(h5ad_path):
        print(f"[ERROR] Input h5ad not found: {h5ad_path}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(covariate_tsv):
        print(f"[ERROR] Covariate TSV not found: {covariate_tsv}", file=sys.stderr)
        sys.exit(1)

    # --- Load AnnData and enforce string-based obs_names ---
    adata = ad.read_h5ad(h5ad_path)
    if isinstance(adata.obs.index, pd.RangeIndex) or adata.obs.index.dtype != object:
        adata.obs.index = adata.obs.index.map(str)
    if not isinstance(adata.obs.index[0], str):
        adata.obs.index = adata.obs.index.astype(str)
    adata.obs_names = pd.Index(adata.obs.index, dtype=object)

    # --- Load and validate covariate table ---
    cov = pd.read_csv(covariate_tsv, sep="\t")
    _assert_no_multiindex(cov, "covariate table")

    required = [library_col, sample_col, covariate_col]
    for c in required:
        if c not in cov.columns:
            print(f"[ERROR] Covariate file missing required column '{c}'", file=sys.stderr)
            sys.exit(1)

    if library_col not in adata.obs.columns:
        print(f"[ERROR] '{library_col}' not found in adata.obs", file=sys.stderr)
        sys.exit(1)

    # --- Merge covariates by library and restore original cell index ---
    obs = adata.obs.copy()
    obs[library_col] = obs[library_col].astype(str)
    cov[library_col] = cov[library_col].astype(str)

    merged = obs.merge(
        cov[[library_col, sample_col, covariate_col]].drop_duplicates(),
        on=library_col,
        how="left"
    )
    # critical fix: restore cell barcodes as index after merge
    merged.index = adata.obs_names.astype(str)

    # --- Warn on missing covariate/sample data ---
    if merged[[sample_col, covariate_col]].isna().any().any():
        n_missing = merged[[sample_col, covariate_col]].isna().any(axis=1).sum()
        print(f"[WARN] {n_missing} cells have missing covariate/sample after merge.")

    adata.obs = merged
    return adata

# --------- Step 3: pseudobulk per (population × sample) (optional) -------- #
def compute_pseudobulk_per_population(adata, population_col, sample_col, covariate_col, min_cells, outdir):
    if population_col not in adata.obs.columns or sample_col not in adata.obs.columns:
        print("[ERROR] '{}' or '{}' not in adata.obs".format(population_col, sample_col), file=sys.stderr)
        sys.exit(1)

    use_counts = "counts" in adata.layers
    if use_counts:
        print("[INFO] Using adata.layers['counts'] for pseudobulk.")
    else:
        print("[WARN] No counts layer found; using adata.X (may be normalized/logged).")

    obs = adata.obs[[population_col, sample_col]].astype(str).copy()
    obs["group"] = obs[population_col] + "|" + obs[sample_col]
    vc = obs["group"].value_counts()
    valid_groups = vc[vc >= int(min_cells)].index.tolist()
    print("[INFO] Valid pseudobulk groups: {} (min_cells={})".format(len(valid_groups), min_cells))

    mat = adata.layers["counts"] if use_counts else adata.X
    mat = _ensure_numeric_matrix(mat)

    rows = []
    linear_rows = []
    colnames = []

    for gid in valid_groups:
        pop, samp = gid.split("|", 1)
        idx = obs.index[obs["group"] == gid]
        if idx.size == 0:
            continue

        sub = adata[adata.obs_names.isin(idx), :].copy()
        X = sub.layers["counts"] if use_counts else sub.X

        if sps.issparse(X):
            summed = np.asarray(X.sum(axis=0)).ravel()
        else:
            X = _ensure_numeric_matrix(X)
            summed = np.sum(X, axis=0)

        total = float(np.sum(summed))
        if total <= 0.0:
            cptt = np.zeros_like(summed, dtype=float)
        else:
            cptt = (summed / total) * 1.0e4
        log_cptt = np.log2(cptt + 1.0)

        rows.append(log_cptt)
        linear_rows.append(cptt)
        colnames.append(gid)

    if len(rows) == 0:
        print("[ERROR] No pseudobulk groups passed threshold.", file=sys.stderr)
        sys.exit(1)

    log_matrix = np.vstack(rows).T  # genes × groups
    lin_matrix = np.vstack(linear_rows).T
    df = pd.DataFrame(log_matrix, index=adata.var_names.astype(str), columns=colnames)
    method = "counts_log2CP10k" if use_counts else "avgExpr_like"

    # --- Ensure proper gene identifiers for pseudobulk rows ---
    if "gene_symbols" in adata.var.columns:
        df.index = adata.var["gene_symbols"].astype(str)
    elif "features" in adata.var.columns:
        df.index = adata.var["features"].astype(str)
    elif not isinstance(df.index[0], str):
        df.index = df.index.astype(str)
    print("[INFO] Pseudobulk gene index set ({} rows)".format(len(df.index)))

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    tsv_path = os.path.join(outdir, "pseudobulk_{}_{}_by_{}.tsv".format(method, population_col, sample_col))
    df.to_csv(tsv_path, sep="\t")
    print("[INFO] Wrote pseudobulk TSV: {}".format(tsv_path))

    # --- Build AnnData with metadata merged from single-cell obs ---
    # Interim h5ad with pseudobulk matrix
    pb_adata = ad.AnnData(X=df.values.T)  # obs = groups, var = genes
    pb_adata.var.index = df.index.astype(str)
    pb_adata.layers["counts"] = np.asarray(lin_matrix.T, dtype=float)

    # --- Decompose and define metadata ---
    populations = [c.split("|", 1)[0] for c in df.columns]
    samples = [c.split("|", 1)[1] for c in df.columns]

    # The obs index (obs_names) becomes the population label directly
    pb_adata.obs.index = pd.Index(populations, name=population_col)

    # Retain metadata for clarity
    pb_adata.obs.index = pd.Index(df.columns, name="GroupID")
    pb_adata.obs[population_col] = populations
    pb_adata.obs["Sample"] = samples
    pb_adata.obs["Population"] = populations
    pb_adata.uns["pseudobulk_method"] = "pseudobulk"
    pb_adata.uns["log1p"] = {"base": 2.0}
    # Preserve lineage hierarchy if present in input AnnData
    lineage_order = _extract_lineage_order(adata, population_col)
    if lineage_order is not None:
        pb_adata.uns["lineage_order"] = lineage_order
        print(f"[INFO] Preserved lineage_order from input (n={len(lineage_order)}).")
    elif "lineage_order" not in getattr(adata, "uns", {}):
        print("[WARN] lineage_order not found in input AnnData; skipping preservation.")


    # --- Attach covariate metadata from original adata.obs ---
    if covariate_col in adata.obs.columns and sample_col in adata.obs.columns:
        sample_cov = (
            adata.obs[[sample_col, covariate_col]]
            .drop_duplicates()
            .set_index(sample_col)[covariate_col]
        )
        pb_adata.obs[covariate_col] = pb_adata.obs["Sample"].map(sample_cov)
        n_missing = pb_adata.obs[covariate_col].isna().sum()
        if n_missing:
            print(f"[WARN] {n_missing} pseudobulk groups missing '{covariate_col}' after mapping.")
        else:
            print(f"[INFO] Attached covariate column '{covariate_col}' to pseudobulk metadata.")
    else:
        print(f"[WARN] Could not attach covariate '{covariate_col}' — column missing in adata.obs.")

    # --- FILTER: Remove invalid pseudobulks before DE ---
    print("[INFO] Filtering low-quality pseudobulks...")

    # 1) Drop pseudobulks with fewer than the required number of cells
    if "n_cells" in pb_adata.obs.columns:
        before = pb_adata.n_obs
        pb_adata = pb_adata[pb_adata.obs["n_cells"] >= min_cells].copy()
        after = pb_adata.n_obs
        #print(f"[DEBUG] Removed {before - after} pseudobulks with <{min_cells} cells")

    # 2) Drop cell types lacking ≥2 pseudobulks per condition
    valid_groups = []
    for pop in pb_adata.obs[population_col].unique():
        sub = pb_adata.obs[pb_adata.obs[population_col] == pop]
        counts = sub[covariate_col].value_counts()
        if all(counts.get(c, 0) >= 2 for c in counts.index):
            valid_groups.append(pop)
        else:
            #print(f"[WARN] Removing '{pop}' (insufficient pseudobulks per condition: {dict(counts)})")
            pass

    before = pb_adata.n_obs
    pb_adata = pb_adata[pb_adata.obs[population_col].isin(valid_groups)].copy()
    after = pb_adata.n_obs
    print(f"[DEBUG] Removed {before - after} pseudobulks lacking >=2 replicates per condition | {after} retained / {before} total")

    h5ad_path = os.path.join(outdir, "pseudobulk_{}_{}.h5ad".format(method, population_col))
    pb_adata.write(h5ad_path)
    print("[INFO] Wrote pseudobulk h5ad: {}".format(h5ad_path))

    # --- Optional QC check ---
    pseudobulk_value_check = False
    if pseudobulk_value_check:
        gene = "SLC39A8"
        if gene in df.index:
            print("\n[QC] Pseudobulk log2CP10k expression for", gene)
            info = df.loc[gene]
            records = []
            for gid, val in info.items():
                if gid.startswith("AT1|"):
                    pop, samp = gid.split("|", 1)
                    records.append((samp, val))
            if records:
                if "Sample" in pb_adata.obs.columns and covariate_col in pb_adata.obs.columns:
                    cond_map = dict(pb_adata.obs[["Sample", covariate_col]].drop_duplicates().values)
                else:
                    cond_map = {}
                print("{:<25s} {:<10s} {:>10s}".format("Sample", covariate_col, "log2CP10k"))
                for samp, val in sorted(records):
                    cond = cond_map.get(samp, "NA")
                    print("{:<25s} {:<10s} {:>10.3f}".format(samp, cond, val))
            else:
                print("[QC] No AT1 entries found for", gene)
        else:
            print("[QC] Gene", gene, "not found in pseudobulk matrix.")
        sys.exit("[INFO] QC check complete. Exiting early after pseudobulk verification.")

    return tsv_path, h5ad_path


# --------- Step 4: per-population DE and pooled/global/co-reg logic -------- #

def _rank_genes_scanpy(two_group_adata, groupby, case_label, control_label, method):
    # Run Scanpy DE and return names, pvals_adj, log2fc, raw pvals as Series for the case vs control
    two_for_rank = _prepare_scanpy_rank_input(two_group_adata)

    sc.tl.rank_genes_groups(two_for_rank,
                            groupby=groupby,
                            groups=[case_label],
                            reference=control_label,
                            method=method,
                            use_raw=False)
    rg = two_for_rank.uns["rank_genes_groups"]
    two_group_adata.uns["rank_genes_groups"] = rg
    names = pd.Index(rg["names"][case_label]).astype(str)
    pvals_adj = pd.Series(rg["pvals_adj"][case_label], index=names, name="fdr").astype(float)
    pvals_raw = pd.Series(rg["pvals"][case_label], index=names, name="pval").astype(float)
    logfc_ln = pd.Series(rg["logfoldchanges"][case_label], index=names, name="logfc_ln").astype(float)
    logfc = pd.Series(_ln_to_log2(logfc_ln.values), index=names, name="log2fc")
    return names, pvals_adj, logfc, pvals_raw


def _prepare_scanpy_rank_input(two_group_adata):
    """Ensure matrix is log-transformed before calling Scanpy rank_genes_groups."""
    if "log1p" in two_group_adata.uns:
        return two_group_adata

    X = two_group_adata.X
    max_val = None
    if sps.issparse(X):
        try:
            max_val = X.max()
            if hasattr(max_val, "A"):
                max_val = max_val.A.max()
            else:
                max_val = float(max_val)
        except Exception:
            max_val = None
    else:
        try:
            max_val = float(np.max(X))
        except Exception:
            max_val = None

    # Heuristic: raw counts typically exceed 20, whereas log-normalised values are ~<15
    if max_val is not None and max_val > 20:
        ad = two_group_adata.copy()
        ad.raw = None
        if not sps.issparse(ad.X):
            ad.X = np.asarray(ad.X, dtype=float)
        sc.pp.normalize_total(ad, target_sum=1e4, inplace=True)
        sc.pp.log1p(ad)
        return ad
    return two_group_adata


def _extract_expression_matrix(pop_data):
    """Return matrix, gene names, and logging metadata for fold-change computation."""
    if getattr(pop_data, "raw", None) is not None:
        matrix = pop_data.raw.X
        genes = np.asarray(pop_data.raw.var_names.astype(str))
        return matrix, genes, False, None

    if "counts" in getattr(pop_data, "layers", {}):
        matrix = pop_data.layers["counts"]
        genes = np.asarray(pop_data.var_names.astype(str))
        return matrix, genes, False, None

    matrix = pop_data.X
    genes = np.asarray(pop_data.var_names.astype(str))

    log1p_info = pop_data.uns.get("log1p", {})
    logged = bool(log1p_info)
    log_base = log1p_info.get("base", None)

    if not logged:
        max_val = None
        if sps.issparse(matrix):
            if matrix.nnz > 0:
                max_val = float(matrix.max())
        else:
            if matrix.size > 0:
                max_val = float(np.nanmax(matrix))
        if max_val is not None and max_val <= 20:
            logged = True

    return matrix, genes, logged, log_base


def _matrix_to_numpy(matrix):
    if sps.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix, dtype=float)


def _compute_pseudobulk_log2fc(pop_data, covariate_col, case_label, control_label):
    cond = pop_data.obs[covariate_col].astype(str).values
    case_mask = cond == str(case_label)
    ctrl_mask = cond == str(control_label)
    if case_mask.sum() == 0 or ctrl_mask.sum() == 0:
        return pd.DataFrame(columns=["log2fc", "case_mean", "control_mean"])

    matrix, genes, logged, log_base = _extract_expression_matrix(pop_data)
    matrix = _matrix_to_numpy(matrix)

    if logged:
        if log_base is None or np.isclose(log_base, np.e):
            linear = np.expm1(matrix)
        else:
            linear = np.power(log_base, matrix) - 1.0
    else:
        linear = matrix

    eps = 1.0
    mean_case = linear[case_mask].mean(axis=0)
    mean_ctrl = linear[ctrl_mask].mean(axis=0)
    log2fc = np.log2((mean_case + eps) / (mean_ctrl + eps))
    case_mean_log2 = np.log2(mean_case + eps)
    control_mean_log2 = np.log2(mean_ctrl + eps)
    df = pd.DataFrame(
        {
            "log2fc": log2fc.astype(float),
            "case_mean_log2": case_mean_log2.astype(float),
            "control_mean_log2": control_mean_log2.astype(float),
            "case_mean_linear": mean_case.astype(float),
            "control_mean_linear": mean_ctrl.astype(float),
        },
        index=genes.astype(str),
    )
    if "gene_symbols" in pop_data.var.columns:
        symbol_map = pop_data.var["gene_symbols"].astype(str).to_dict()
        df.index = [symbol_map.get(g, g) for g in df.index]
        df = df.groupby(level=0).mean()
    df.index.name = "gene"
    return df


def _moderated_t_test(adata, covariate_col, case_label, control_label, pop):
    """
    Limma-like moderated t-test with empirical Bayes shrinkage.
    Now performs independent filtering: expression filtering only
    affects FDR correction, not test-statistic computation.
    """
    from scipy import stats
    import numpy as np, pandas as pd
    from statsmodels.stats.multitest import multipletests

    cond = adata.obs[covariate_col].astype(str)
    case_mask = cond == str(case_label)
    ctrl_mask = cond == str(control_label)

    X_case = np.asarray(adata.X[case_mask, :].todense() if hasattr(adata.X, "todense") else adata.X[case_mask, :])
    X_ctrl = np.asarray(adata.X[ctrl_mask, :].todense() if hasattr(adata.X, "todense") else adata.X[ctrl_mask, :])
    genes = np.asarray(adata.var_names)

    n_case, n_ctrl = X_case.shape[0], X_ctrl.shape[0]
    if n_case < 2 or n_ctrl < 2:
        if diagnostic_report:
            print(f"[DEBUG] Skipping {pop} (insufficient pseudobulks per condition after case/control filtering)")
        return pd.DataFrame(columns=["gene","log2fc","t","pval","fdr"]), 0

    # --- Compute t/p for all genes first (no filtering yet) ---
    mean_case = X_case.mean(0)
    mean_ctrl = X_ctrl.mean(0)
    var_case = X_case.var(0, ddof=1)
    var_ctrl = X_ctrl.var(0, ddof=1)
    s2 = (var_case + var_ctrl) / 2

    # Empirical Bayes shrinkage
    s2_prior = np.median(s2)
    shrink_factor = 0.2
    s2_shrunk = (1 - shrink_factor) * s2 + shrink_factor * s2_prior
    se = np.sqrt(s2_shrunk / n_case + s2_shrunk / n_ctrl)
    se = np.maximum(se, 1e-12)

    diff = mean_case - mean_ctrl

    # --- NEW: persist full (unfiltered) log2FC for this population ---
    # Note: this is a side-effect; we collect these across populations later.
    if "full_log2fc" not in adata.uns:
        adata.uns["full_log2fc"] = {}
    # store as a Series aligned to adata.var_names
    adata.uns["full_log2fc"][str(pop)] = pd.Series(
        np.asarray(diff).ravel(), index=np.asarray(adata.var_names), dtype=float
    )

    tvals = diff / se

    pvals = 2 * stats.t.sf(np.abs(tvals), df=n_case + n_ctrl - 2)

    # --- Independent filtering for FDR correction ---
    expr_threshold = 0.1
    min_detected_samples = 2
    X_all = np.vstack([X_case, X_ctrl])
    detected_counts = np.sum(X_all > expr_threshold, axis=0)
    filter_mask = detected_counts >= min_detected_samples

    # Apply BH-FDR to filtered subset only
    filt_mask = filter_mask & np.isfinite(pvals)
    pvals_filt = pvals[filt_mask]
    genes_filt = genes[filt_mask]

    if len(pvals_filt) == 0:
        print(f"[WARN] {pop}: No genes passed expression filter.")
        return pd.DataFrame(columns=["gene","log2fc","t","pval","fdr"]), 0

    _, fdr_filt, _, _ = multipletests(pvals_filt, method="fdr_bh")

    # Map FDRs back to full array (untested = NaN)
    fdr = np.full_like(pvals, np.nan, dtype=float)
    fdr[filt_mask] = fdr_filt

    # --- Diagnostics ---
    tested_genes = int(np.sum(filt_mask))
    adata.uns[f"tested_genes_{pop}"] = tested_genes

    if diagnostic_report:
        print(f"[DEBUG] {pop} moderated t-test (EB):")
        print(f"  [INFO]  total genes = {len(genes)}")
        print(f"  [INFO]  tested (after filter) = {tested_genes}")
        print(f"  [INFO]  min raw p = {np.nanmin(pvals):.3e}, median = {np.nanmedian(pvals):.3e}")
        if np.isfinite(fdr).any():
            print(f"  [INFO]  min adj p = {np.nanmin(fdr):.3e}, median = {np.nanmedian(fdr):.3e}")
        if "PROSER2" in genes:
            gi = np.where(genes == "PROSER2")[0][0]
            print(f"[DEBUG] PROSER2 (moderated t): p={pvals[gi]:.3e}, FDR={fdr[gi]:.3e}")

    # --- Output for all genes ---
    df_out = pd.DataFrame({
        "gene": genes.astype(str),
        "log2fc": diff.astype(float),
        "t": tvals.astype(float),
        "pval": pvals.astype(float),
        "fdr": fdr.astype(float)
    })

    # Drop only rows with non-finite pvals for safety
    df_out = df_out.loc[np.isfinite(df_out["pval"])].copy()

    return df_out, tested_genes

def run_de_for_comparisons(adata,
                           population_col,
                           covariate_col,
                           case_label,
                           control_label,
                           method,
                           alpha,
                           fc_thresh,
                           min_cells_per_group,
                           use_rawp=False,
                           progress_callback=None):

    # Filter to cells with both labels present
    keep = adata.obs[covariate_col].isin([case_label, control_label])
    sub = adata[keep].copy()
    total_cells = sub.n_obs

    min_cells = int(min_cells_per_group)
    if total_cells < 200:
        min_cells = max(4, min_cells)

    populations = _ordered_categories_from_obs(sub.obs, population_col)
    results_rows = []
    long_stats = []
    per_pop_de = {}
    all_genes = set(sub.var_names.astype(str))
    corrected_fc = None

    if not sps.issparse(sub.X):
        sub.X = _ensure_numeric_matrix(sub.X)

    # Per-population DE
    all_fold_values = {}
    for pop in populations:
        pop_mask = sub.obs[population_col].astype(str) == str(pop)
        pop_data = sub[pop_mask].copy()

        # size checks
        grp_sizes = pop_data.obs[covariate_col].value_counts()
        n_case = int(grp_sizes.get(case_label, 0))
        n_ctrl = int(grp_sizes.get(control_label, 0))
        if n_case < min_cells or n_ctrl < min_cells:
            results_rows.append([pop, n_case, n_ctrl, 0, 0])
            if progress_callback is not None:
                progress_callback()
            continue


        # --- Perform DE for individual population ---
        two = pop_data.copy()

        # --- Choose test method based on pseudobulk flag ---
        if "pseudobulk" in adata.uns.get("pseudobulk_method", ""):
            # Use limma-like moderated t-test for pseudobulk
            df, tested_genes = _moderated_t_test(two, covariate_col, case_label, control_label, pop)
            df["gene"] = df["gene"].astype(str)
            names = df["gene"].astype(str)
            fdr = df["fdr"].astype(float)
            log2fc = df["log2fc"].astype(float)
            pvals_series = df["pval"].astype(float)

            # --- capture full log2FC vector for this population ---
            if "full_log2fc" in two.uns and str(pop) in two.uns["full_log2fc"]:
                all_fold_values[str(pop)] = two.uns["full_log2fc"][str(pop)]


        else:
            # Default Scanpy DE method
            try:
                names, fdr, log2fc, pvals_series = _rank_genes_scanpy(two, covariate_col, case_label, control_label, method)
                tested_genes = len(names)
            except ValueError as e:
                err = str(e)
                if (
                    "only contain one sample" in err
                    or "needs to be one of groupby" in err
                    or "does not exist" in err
                ):
                    if diagnostic_report:
                        print(f"[WARN] Skipping population '{pop}' (invalid group composition: {err.splitlines()[0]})")
                    results_rows.append([pop, n_case, n_ctrl, n_deg, tested_genes])
                    if progress_callback is not None:
                        progress_callback()
                    continue
                else:
                    raise

        names_idx = pd.Index(names).astype(str)
        fdr_series = pd.Series(np.asarray(fdr, dtype=float), index=names_idx)
        log2fc_series = pd.Series(np.asarray(log2fc, dtype=float), index=names_idx)
        pvals_series = pd.Series(np.asarray(pvals_series, dtype=float), index=names_idx)

        df = pd.DataFrame(
            {
                "gene": names_idx,
                "log2fc": log2fc_series.values,
                "fdr": fdr_series.values,
                "pval": pvals_series.values,
            }
        )

        if "gene_symbols" in pop_data.var.columns:
            id_map = pop_data.var["gene_symbols"].astype(str).to_dict()
            df["gene"] = df["gene"].map(lambda g: id_map.get(g, g))
            df["gene"] = df["gene"].astype(str)

        df["population"] = pop
        df["case_label"] = case_label
        df["control_label"] = control_label
        df["n_case"] = n_case
        df["n_control"] = n_ctrl

        log2fc_values = df["log2fc"].to_numpy(dtype=float)
        fdr_values = df["fdr"].to_numpy(dtype=float)
        pvals_values = df["pval"].to_numpy(dtype=float)

        if use_rawp:
            if diagnostic_report:
                print(f"[INFO] Using raw p-values for significance (alpha={alpha})")
            p_base = pvals_values
        else:
            p_base = fdr_values

        sig = (np.abs(log2fc_values) > np.log2(float(fc_thresh))) & (p_base < float(alpha))
        deg_names = pd.Index(df["gene"])[sig]
        n_deg = int(sig.sum())


        # --- Correct fold computation: true mean-based log2 fold (per gene across cells) ---
        corrected_stats = _compute_pseudobulk_log2fc(pop_data, covariate_col, case_label, control_label)
        if not corrected_stats.empty:
            all_fold_values[str(pop)] = corrected_stats["log2fc"]
            updated_fc = df["gene"].map(corrected_stats["log2fc"])
            df.loc[updated_fc.notna(), "log2fc"] = updated_fc[updated_fc.notna()]

            case_map = corrected_stats.get("case_mean_log2", corrected_stats.get("case_mean"))
            ctrl_map = corrected_stats.get("control_mean_log2", corrected_stats.get("control_mean"))
            df["case_mean_expr"] = df["gene"].map(case_map) if case_map is not None else np.nan
            df["control_mean_expr"] = df["gene"].map(ctrl_map) if ctrl_map is not None else np.nan
        else:
            if diagnostic_report:
                print(f"[WARN] No corrected fold-change computed for population {pop} (insufficient data).")
            df["case_mean_expr"] = np.nan
            df["control_mean_expr"] = np.nan

        if diagnostic_report:
            # --- DEBUGGING: PROSER2 in AT1 ---
            if str(pop) == "AT1" and "PROSER2" in df["gene"].values:
                row = df.loc[df["gene"] == "PROSER2"].iloc[0]
                print("\n[DEBUG] PROSER2 in AT1:")
                print(f"  log2FC: {row['log2fc']:.3f}")
                print(f"  FDR: {row['fdr']:.3e}")
                print(f"  case/control: {row['case']} vs {row['control']}")
                print(f"  n_case={row['n_case']}, n_control={row['n_control']}")

        # --- Optional debugging ---
        fold_value_check = False
        if fold_value_check:
            gene = "SLC39A8"
            debug_population = "AT1"

            print("\n[DEBUG] Fold check active")
            print(f"[DEBUG] Current population: '{pop}'")
            print(f"[DEBUG] corrected_fc index (first 5): {list(corrected_fc.index[:5])}")

            gene_matches = [g for g in corrected_fc.index if gene.lower() in g.lower()]
            pop_match = debug_population.lower() in str(pop).lower()

            if pop_match and len(gene_matches) > 0:
                g = gene_matches[0]
                print(f"\n[DEBUG] {g} in {pop}")
                print(f"[DEBUG] Control mean = {mean_ctrl.loc[g]:.3f}")
                print(f"[DEBUG] Case mean = {mean_case.loc[g]:.3f}")
                print(f"[DEBUG] log2FC = {corrected_fc.loc[g]:.3f}")
            else:
                print(f"[DEBUG] No direct match found — pop_match={pop_match}, gene_match_count={len(gene_matches)}")


        per_pop_de[str(pop)] = df.set_index("gene")
        # detailed only for DEGs

        # --- Defensive check to avoid KeyErrors ---
        missing_genes = [g for g in deg_names if g not in df["gene"].values]
        if len(missing_genes) == len(deg_names):
            if diagnostic_report:
                print(f"[WARN] No DE genes from {pop} match df['gene']; likely symbol mismatch ({len(df)} genes tested).")
                print(f"[DEBUG] df['gene'] example: {df['gene'].head().tolist()[:5]}")
            continue
        else:
            found_genes = [g for g in deg_names if g in df["gene"].values]
            long_stats.append(df[df["gene"].isin(found_genes)].copy())

        # Count only genes with valid (finite) t-values or fdr for tested genes
        tested_genes = int(adata.uns.get(f"tested_genes_{pop}", df["fdr"].notna().sum()))
        results_rows.append([pop, n_case, n_ctrl, n_deg, tested_genes])

        # --- progress update for this finished population ---
        if progress_callback is not None:
            progress_callback()



    summary = pd.DataFrame(results_rows,
                           columns=["population", "n_case", "n_control", "num_DEG", "tested_genes"])
    detailed = pd.concat(long_stats, axis=0) if len(long_stats) > 0 else pd.DataFrame(
        columns=[
            "gene",
            "population",
            "log2fc",
            "fdr",
            "pval",
            "case",
            "control",
            "n_case",
            "n_control",
            "case_mean_expr",
            "control_mean_expr",
        ]
    )

    # Pooled overall test pO (all cells case vs control, regardless of population)
    pooled = sub.copy()
    # Ensure both groups present
    gvc = pooled.obs[covariate_col].value_counts()
    if int(gvc.get(case_label, 0)) >= min_cells and int(gvc.get(control_label, 0)) >= min_cells:
        try:
            if method == "scanpy":
                # Single-cell mode
                names_o, fdr_o, log2fc_o = _rank_genes_scanpy(pooled, covariate_col, case_label, control_label, method)
                pooled_overall = pd.DataFrame({
                    "gene": names_o.values,
                    "log2fc": log2fc_o.values,
                    "fdr": fdr_o.values
                }).set_index("gene")
            else:
                # Pseudobulk mode
                df_o, _ = _moderated_t_test(pooled, covariate_col, case_label, control_label, "pooled_overall")
                pooled_overall = df_o.set_index("gene")[["log2fc", "fdr"]]
        except Exception as e:
            print(f"[WARN] Overall pooled DE failed: {e}")
            pooled_overall = pd.DataFrame(index=pd.Index([], dtype=str), columns=["log2fc", "fdr"])
    else:
        print("[WARN] Overall pooled DE skipped due to insufficient cells.")
        pooled_overall = pd.DataFrame(index=pd.Index([], dtype=str), columns=["log2fc", "fdr"])


    # Assemble per-gene vector over populations (for co-reg patterning)
    # vector entries: 1 if (p<0.10 & log2fc>0), -1 if (p<0.10 & log2fc<0), else 0
    pattern_mat = {}
    pop_order = populations
    for pop in pop_order:
        dfp = per_pop_de.get(str(pop), None)
        if dfp is None:
            continue
        # For *all* genes tested, we need p-values; reconstruct: if absent, set to 1, log2fc=0
        # Here we have only DEGs stored per_pop_de; build light dict from Scanpy per population
        # For robustness, mark only genes we kept; others as 0
        for g in dfp.index:
            if g not in pattern_mat:
                pattern_mat[g] = [0] * len(pop_order)
    # If no DEGs at all, pattern_mat remains empty
    # We will rebuild pattern_mat from Scanpy outputs if needed
    # Extract p<0.10 using available per_pop_de frames (their fdr is BH<alpha; we need raw p<0.10).
    # As raw p are not available from stored df, approximate: use fdr as a proxy upper bound.
    # To improve sensitivity for pattern coding, fall back to fdr<0.10 here.
    for j, pop in enumerate(pop_order):
        dfp = per_pop_de.get(str(pop), None)
        if dfp is None or dfp.empty:
            continue

        # --- define the significance metric consistent with --use_rawp ---
        if use_rawp and "pval" in dfp.columns:
            dfp["sig_metric"] = dfp["pval"]
        else:
            dfp["sig_metric"] = dfp["fdr"]

        sel = dfp.index[(dfp["sig_metric"].values < alpha)]

        sgn = np.sign(dfp.loc[sel, "log2fc"].values)
        for k, g in enumerate(sel):
            if g not in pattern_mat:
                pattern_mat[g] = [0] * len(pop_order)
            pattern_mat[g][j] = int(1 if sgn[k] > 0 else (-1 if sgn[k] < 0 else 0))

    # Patterns
    pattern_series = pd.Series({g: tuple(v) for g, v in pattern_mat.items()})
    # Compute counts excluding pure global (all 1 or all -1) and local (only one non-zero)
    def _is_global(t):
        if len(t) == 0:
            return False
        all_pos = all([x == 1 for x in t if True])
        all_neg = all([x == -1 for x in t if True])
        return all_pos or all_neg

    def _is_local(t):
        return sum([1 for x in t if x != 0]) == 1

    subset_for_freq = pattern_series[~pattern_series.apply(_is_global) & ~pattern_series.apply(_is_local)]
    pattern_freq = subset_for_freq.value_counts().sort_values(ascending=False)
    top_patterns = [list(pat) for pat in pattern_freq.index[:4]]

    # For each top pattern, pool those populations and compute pooled DE -> p(co-reg)
    coreg_tables = []
    for pat in top_patterns:
        idx_pops = [i for i, v in enumerate(pat) if v != 0]
        if len(idx_pops) == 0:
            continue
        pops_sel = [pop_order[i] for i in idx_pops]

        mask_case = (sub.obs[covariate_col] == case_label) & (sub.obs[population_col].isin(pops_sel))
        mask_ctrl = (sub.obs[covariate_col] == control_label) & (sub.obs[population_col].isin(pops_sel))
        if int(mask_case.sum()) < min_cells or int(mask_ctrl.sum()) < min_cells:
            continue
        pooled_pat = sub[mask_case | mask_ctrl].copy()
        coreg_name = _pattern_to_label(pat, pop_order)
        try:
            if method == "scanpy":
                names_c, fdr_c, log2fc_c, pvals_c = _rank_genes_scanpy(
                    pooled_pat, covariate_col, case_label, control_label, method
                )
                dfc = pd.DataFrame({
                    "gene": names_c.values,
                    "log2fc": log2fc_c.values,
                    "fdr": fdr_c.values,
                    "pval": pvals_c.values,
                })
            else:
                dfc, _ = _moderated_t_test(pooled_pat, covariate_col, case_label, control_label, coreg_name)
                keep_cols = [c for c in ["gene", "log2fc", "fdr", "pval"] if c in dfc.columns]
                dfc = dfc.loc[:, keep_cols]
            if diagnostic_report:
                print(f"[DEBUG] Co-reg pattern '{coreg_name}' using pops {pops_sel} "
                      f"(n_case={int(mask_case.sum())}, n_ctrl={int(mask_ctrl.sum())})")
            dfc["pattern"] = coreg_name
            coreg_tables.append(dfc)
        except Exception as e:
            print(f"[WARN] Co-reg pooled DE failed for pattern {pat} ({coreg_name}): {e}")
            if diagnostic_report:
                print(f"[DEBUG] Failed co-reg pops: {pops_sel}")


    coreg_df = pd.concat(coreg_tables, axis=0) if len(coreg_tables) > 0 else pd.DataFrame(
        columns=["gene", "log2fc", "fdr", "pattern"]
    )

    # Final gene assignment to group (global/local/co-regulated)
    # Start from genes significant in at least one population (original DE step)
    if detailed.empty:
        assigned = pd.DataFrame(columns=["gene", "group", "key_p", "key_log2fc", "population_or_pattern"])
    else:
        deg_genes = sorted(detailed["gene"].unique().tolist())
        # Compute global criterion: pO<0.05 and pi<0.05 for at least 2 populations with consistent sign
        # Use pooled_overall.fdr as pO; for pi, we approximate using fdr per population
        global_hits = []
        local_hits = []
        coreg_hits = []

        # Build per-pop lookup
        perpop_lookup = {}
        for pop in pop_order:
            dfp = per_pop_de.get(str(pop), None)
            if dfp is None or dfp.empty:
                continue
            perpop_lookup[pop] = dfp

        for g in deg_genes:
            # collect per-pop signs where significance metric < alpha
            s = []
            for pop in pop_order:
                dfp = perpop_lookup.get(pop, None)
                if dfp is None or (g not in dfp.index):
                    continue

                sig_val = (
                    float(dfp.loc[g, "pval"])
                    if use_rawp and "pval" in dfp.columns
                    else float(dfp.loc[g, "fdr"])
                )

                if sig_val < alpha:
                    s.append(np.sign(float(dfp.loc[g, "log2fc"])))

            consistent = (len(s) >= 2) and (all([x > 0 for x in s]) or all([x < 0 for x in s]))

            is_global = False
            key_log2fc = np.nan
            key_p = np.nan
            key_label = ""

            if (not pooled_overall.empty) and (g in pooled_overall.index):
                sig_val = (
                    float(pooled_overall.loc[g, "pval"])
                    if use_rawp and "pval" in pooled_overall.columns
                    else float(pooled_overall.loc[g, "fdr"])
                )
                if sig_val < alpha and consistent:
                    is_global = True
                    key_p = float(pooled_overall.loc[g, "fdr"])
                    key_log2fc = float(pooled_overall.loc[g, "log2fc"])
                    key_label = "overall"

            if is_global:
                global_hits.append([g, "global", key_p, key_log2fc, key_label])
                continue

            # local: only one population with DE (fdr<0.05 & |log2fc|>log2(fc_thresh))
            de_pops = []
            min_p = 1.0
            best_fc = 0.0
            best_pop = ""
            for pop in pop_order:
                dfp = perpop_lookup.get(pop, None)
                if dfp is None or (g not in dfp.index):
                    continue
                sig_val = (
                    float(dfp.loc[g, "pval"])
                    if use_rawp and "pval" in dfp.columns
                    else float(dfp.loc[g, "fdr"])
                )
                if sig_val < alpha and np.abs(float(dfp.loc[g, "log2fc"])) > np.log2(fc_thresh):
                    de_pops.append(pop)
                    pval = float(dfp.loc[g, "fdr"])
                    if pval < min_p:
                        min_p = pval
                        best_fc = float(dfp.loc[g, "log2fc"])
                        best_pop = pop

            if len(de_pops) == 1:
                local_hits.append([g, "local", min_p, best_fc, best_pop])
                continue

            # co-reg: among top patterns, if pooled comparison has smallest p
            if not coreg_df.empty and (g in coreg_df["gene"].values):
                rows = coreg_df[coreg_df["gene"] == g]
                if not rows.empty:
                    valid_rows = rows.loc[np.isfinite(rows["fdr"].astype(float))]
                    if valid_rows.empty:
                        continue  # skip if all FDRs are NaN

                    idx_min = valid_rows["fdr"].astype(float).idxmin()
                    r = valid_rows.loc[idx_min]

                    # r might be a Series (single row) or a DataFrame row
                    if isinstance(r, pd.Series):
                        fdr_val = float(r["fdr"])
                        log2fc_val = float(r["log2fc"])
                        pattern_str = str(r["pattern"])
                    else:
                        fdr_val = float(r["fdr"].iloc[0])
                        log2fc_val = float(r["log2fc"].iloc[0])
                        pattern_str = str(r["pattern"].iloc[0])

                    coreg_hits.append([g, "co-regulated", fdr_val, log2fc_val, pattern_str])

                    continue

            # If ambiguous, assign to local with smallest per-pop p (if any) else skip
            if len(de_pops) >= 2:
                # choose the pop with smallest fdr
                min_p = 1.0
                best_fc = 0.0
                best_pop = ""
                for pop in de_pops:
                    pval = (
                        float(perpop_lookup[pop].loc[g, "pval"])
                        if use_rawp and "pval" in perpop_lookup[pop].columns
                        else float(perpop_lookup[pop].loc[g, "fdr"])
                    )
                    if pval < min_p:
                        min_p = pval
                        best_fc = float(perpop_lookup[pop].loc[g, "log2fc"])
                        best_pop = pop
                local_hits.append([g, "local", min_p, best_fc, best_pop])

        assigned = pd.DataFrame(global_hits + local_hits + coreg_hits,
                                columns=["gene", "group", "key_p", "key_log2fc", "population_or_pattern"])

    # Aggregate outputs
    de_store = {
        "summary_per_population": summary,
        "detailed_deg": detailed,
        "pooled_overall": pooled_overall,
        "coreg_pooled": coreg_df,
        "assigned_groups": assigned,
        "population_order": pop_order,
        "case_label": str(case_label),
        "control_label": str(control_label),
        "method": str(method),
        "alpha": float(alpha),
        "fc_thresh": float(fc_thresh),
        "min_cells_per_group": int(min_cells),
        "per_population_deg": {k: v.copy() for k, v in per_pop_de.items()},
    }

    # --- NEW: build a complete fold-change matrix from per-population stored vectors ---
    try:
        if len(all_fold_values) > 0:
            fold_matrix = pd.DataFrame(all_fold_values)
        else:
            fold_matrix = pd.DataFrame(adata.uns.get("full_log2fc", {}))
        pop_order = de_store.get("population_order", list(fold_matrix.columns))
        fold_matrix = fold_matrix.reindex(columns=[str(p) for p in pop_order if str(p) in fold_matrix.columns])
        fold_matrix.index = fold_matrix.index.astype(str)
        de_store["fold_matrix"] = fold_matrix
        if diagnostic_report:
            print(f"[DEBUG] fold_matrix built: {fold_matrix.shape[0]} genes × {fold_matrix.shape[1]} populations")
    except Exception as e:
        print(f"[WARN] fold_matrix assembly failed: {e}")

    # Preserve lineage order from input AnnData so downstream heatmap can access it
    lineage_order = _extract_lineage_order(adata, population_col)
    if lineage_order is not None:
        de_store["lineage_order"] = lineage_order
        # print(f"[INFO] Preserved lineage_order in de_store (n={len(lineage_order)}).")

    de_store["use_rawp"] = bool(use_rawp)

    # Ensure progress bar completes at the end
    if progress_callback is not None:
        progress_callback()

    return de_store


# -------------------- Step 5: Heatmap (fixed order), TSV ------------------- #

def build_fixed_order_heatmap(
    de_store,
    outdir,
    population_col,
    fc_thresh,
    heatmap_png,
    heatmap_tsv,
    goelite_payload=None,
    show_go_terms=False,
    goelite_max_terms=30,
):
    """
    Build a fixed-order heatmap of log2 fold-changes for all populations,
    restricted to genes significant in any local comparison.
    Global and co-regulated genes are included only if also local-significant.
    Ordering:
      1) Global up/down (ranked by p ascending)
      2) Co-regulated up/down (ranked by p ascending; patterns with more locals first)
      3) Local up/down per population (ranked by p ascending)
    """

    assigned = de_store["assigned_groups"]
    detailed = de_store["detailed_deg"]
    use_rawp = bool(de_store.get("use_rawp", False))
    alpha = float(de_store.get("alpha", 0.05))

    if use_rawp and "pval" in detailed.columns:
        detailed["sig_metric"] = detailed["pval"]
        print(f"[INFO] Heatmap using raw p-values (alpha={alpha}) for inclusion.")
    else:
        detailed["sig_metric"] = detailed["fdr"]
        print(f"[INFO] Heatmap using FDR (alpha={alpha}) for inclusion.")


    pop_order = de_store["population_order"]
    pooled_overall = de_store.get("pooled_overall", pd.DataFrame())
    coreg_df = de_store.get("coreg_pooled", pd.DataFrame())

    lineage_order = None
    lineage_source = None
    if "lineage_order" in de_store:
        try:
            lineage_order = _validate_lineage_order_candidates(de_store["lineage_order"], pop_order, population_col)
            lineage_source = "de_store"
            print(f"[DEBUG] de_store lineage_order present: {de_store['lineage_order']}")
        except ValueError as exc:
            print(f"[WARN] Stored lineage_order mismatch for '{population_col}': {exc}")
    elif "adata" in de_store and hasattr(de_store["adata"], "uns") and "lineage_order" in de_store["adata"].uns:
        print(f"[DEBUG] adata.uns keys: {list(de_store['adata'].uns.keys())}")
        try:
            lineage_order = _validate_lineage_order_candidates(
                de_store["adata"].uns["lineage_order"], pop_order, population_col
            )
            lineage_source = "adata.uns"
        except ValueError as exc:
            print(f"[WARN] lineage_order from embedded adata mismatched for '{population_col}': {exc}")

    if assigned.empty or detailed.empty:
        print("[WARN] No DEGs to plot heatmap.")
        return None, None

    # ------------------------- 1. Identify local-significant genes ------------------------- #
    local_genes = set(assigned.loc[assigned["group"] == "local", "gene"])
    global_genes = set(assigned.loc[assigned["group"] == "global", "gene"])
    coreg_genes = set(assigned.loc[assigned["group"] == "co-regulated", "gene"])
    print(f"[DEBUG] Heatmap gene group counts -> local: {len(local_genes)}, global: {len(global_genes)}, co-reg: {len(coreg_genes)}")

    include_genes = set(assigned["gene"])
    if len(include_genes) == 0:
        print("[WARN] Assigned groups empty; falling back to detailed DEG table.")
        include_genes = set(detailed["gene"])
        if len(include_genes) == 0:
            print("[ERROR] Detailed DEG table empty; skipping heatmap.")
            return None, None
    print(f"[DEBUG] Heatmap include_genes count: {len(include_genes)}")

    # ------------------------- 2. Build log2FC matrix using full fold_matrix ------------------------- #
    if "fold_matrix" not in de_store:
        print("[WARN] fold_matrix not found; defaulting to partial DE table values.")
        fold_df = pd.DataFrame(index=sorted(include_genes))
        for pop in pop_order:
            dfp = detailed[detailed["population"] == str(pop)].set_index("gene")
            fold_df[pop] = dfp["log2fc"].reindex(fold_df.index)
    else:
        full_fold = de_store["fold_matrix"].copy()
        fold_df = full_fold.loc[full_fold.index.intersection(sorted(include_genes))]

    fold_df = fold_df.astype(float)
    print(f"[DEBUG] fold_df shape after filtering: {fold_df.shape[0]} genes × {fold_df.shape[1]} populations")
    extra_genes = [g for g in detailed["gene"].astype(str).unique() if g not in fold_df.index]
    if extra_genes:
        print(f"[DEBUG] Adding {len(extra_genes)} genes absent from fold_df via detailed DEG table.")
        pivot_col = "population" if "population" in detailed.columns else "population_or_pattern"
        extra_matrix = detailed[detailed["gene"].isin(extra_genes)].pivot_table(
            index="gene", columns=pivot_col, values="log2fc", fill_value=0.0
        )
        fold_df = pd.concat([fold_df, extra_matrix], axis=0)

    # ------------------------- 3. Ranking logic ------------------------- #
    global_hits = assigned.loc[assigned["group"] == "global"].copy()
    coreg_hits = assigned.loc[assigned["group"] == "co-regulated"].copy()
    local_hits = assigned.loc[assigned["group"] == "local"].copy()

    def split_up_down(df):
        up = df[df["key_log2fc"] > 0].sort_values("key_p", ascending=True)
        down = df[df["key_log2fc"] < 0].sort_values("key_p", ascending=True)
        return up, down

    # (1) Global: up then down
    g_up, g_down = split_up_down(global_hits)

    # (2) Co-regulated: order patterns by number of local pops represented
    if not coreg_hits.empty:
        coreg_hits["pattern_size"] = coreg_hits["population_or_pattern"].str.count("__") + 1
        coreg_hits = coreg_hits.sort_values(["pattern_size", "key_p"], ascending=[False, True])
    c_up, c_down = split_up_down(coreg_hits)

    # (3) Local: for each population in order, up then down by p-value
    # Build filtered lists now so block sizes reflect rows actually present in fold_df
    g_up_list   = [g for g in g_up["gene"].tolist()   if g in fold_df.index]
    g_down_list = [g for g in g_down["gene"].tolist() if g in fold_df.index]

    # (2) Co-regulated (already sorted above)
    c_up_list = []
    c_up_labels = []
    for _, row in c_up.iterrows():
        gene = str(row["gene"])
        if gene not in fold_df.index:
            continue
        pattern = str(row.get("population_or_pattern", "")).replace(":", "_")
        if pattern:
            label = f"coreg_{pattern}__up"
        else:
            label = "coreg__up"
        c_up_list.append(gene)
        c_up_labels.append(label)

    c_down_list = []
    c_down_labels = []
    for _, row in c_down.iterrows():
        gene = str(row["gene"])
        if gene not in fold_df.index:
            continue
        pattern = str(row.get("population_or_pattern", "")).replace(":", "_")
        if pattern:
            label = f"coreg_{pattern}__down"
        else:
            label = "coreg__down"
        c_down_list.append(gene)
        c_down_labels.append(label)

    # (3) Local per-pop up/down — now ordered by lineage_order if available
    ordered_local_genes = []
    local_block_sizes = []
    local_gene_blocks = []

    # Use lineage-based order for populations if available
    lineage_for_local = lineage_order if lineage_order is not None else pop_order
    lineage_for_local = [p for p in lineage_for_local if p in pop_order]
    if diagnostic_report:
        print(f"[DEBUG] lineage_for_local: {lineage_for_local}")

    for pop in lineage_for_local:
        sub = local_hits[local_hits["population_or_pattern"] == str(pop)]
        u, d = split_up_down(sub)
        u_list = [g for g in u["gene"].tolist() if g in fold_df.index]
        d_list = [g for g in d["gene"].tolist() if g in fold_df.index]
        ordered_local_genes.extend(u_list + d_list)
        local_block_sizes.append(len(u_list))  # up genes for this pop
        local_block_sizes.append(len(d_list))  # down genes for this pop
        pop_label = str(pop).replace(":", "_")
        local_gene_blocks.append((pop_label, u_list, d_list))

    # Combine master order: global → coreg → local
    ordered_genes = g_up_list + g_down_list + c_up_list + c_down_list + ordered_local_genes

    print(f"[DEBUG] Heatmap ordered gene block sizes -> global_up:{len(g_up_list)}, global_down:{len(g_down_list)}, "
          f"coreg_up:{len(c_up_list)}, coreg_down:{len(c_down_list)}, local_total:{len(ordered_local_genes)}")

    extra_genes = [g for g in fold_df.index if g not in ordered_genes]
    if extra_genes:
        print(f"[DEBUG] Appending {len(extra_genes)} unassigned genes to heatmap order.")
        ordered_genes.extend(extra_genes)

    if len(ordered_genes) == 0:
        print("[WARN] No ordered genes found after block assembly; falling back to top genes by significance.")
        fallback = (
            detailed.sort_values("sig_metric", ascending=True)
            .drop_duplicates(subset=["gene"])
            .head(min(200, len(detailed)))
        )
        ordered_genes = fallback["gene"].tolist()
        pivot_col = "population" if "population" in detailed.columns else "population_or_pattern"
        fold_df = detailed.pivot_table(index="gene", columns=pivot_col, values="log2fc", fill_value=0.0)
        fold_df = fold_df.loc[ordered_genes]
        fold_df = fold_df.reindex(columns=[col for col in pop_order if col in fold_df.columns])
        g_up_list, g_down_list, c_up_list, c_down_list = [], [], [], []
        ordered_local_genes = ordered_genes
        block_breaks = []
        local_block_sizes = [len(ordered_genes)]
        extra_genes = []
        local_gene_blocks = []
        print(f"[INFO] Heatmap fallback selected {len(ordered_genes)} genes.")

    # Compute horizontal block boundaries (end indices of each block)
    block_breaks = []
    cursor = 0
    # Global blocks
    if len(g_up_list) > 0:
        cursor += len(g_up_list); block_breaks.append(cursor)
    if len(g_down_list) > 0:
        cursor += len(g_down_list); block_breaks.append(cursor)
    # Co-reg blocks
    if len(c_up_list) > 0:
        cursor += len(c_up_list); block_breaks.append(cursor)
    if len(c_down_list) > 0:
        cursor += len(c_down_list); block_breaks.append(cursor)
    # Local blocks (up then down for each population in pop_order)
    for sz in local_block_sizes:
        if sz > 0:
            cursor += sz
            block_breaks.append(cursor)

    if extra_genes:
        cursor += len(extra_genes)
        block_breaks.append(cursor)

    # Assign row clusters to match cellHarmony differential groups/patterns
    row_clusters = []
    row_clusters.extend(["coreg_global__up"] * len(g_up_list))
    row_clusters.extend(["coreg_global__down"] * len(g_down_list))
    row_clusters.extend(c_up_labels)
    row_clusters.extend(c_down_labels)
    for pop, u_list, d_list in local_gene_blocks:
        row_clusters.extend([f"{pop}__up"] * len(u_list))
        row_clusters.extend([f"{pop}__down"] * len(d_list))
    if extra_genes:
        row_clusters.extend(["unassigned"] * len(extra_genes))
    if len(row_clusters) != len(ordered_genes):
        row_clusters = ["unassigned"] * len(ordered_genes)

    # Recompute boundaries to align with row cluster labels
    block_breaks = []
    for idx in range(1, len(row_clusters)):
        if row_clusters[idx] != row_clusters[idx - 1]:
            block_breaks.append(idx)

    # Reorder the matrix rows by the final order
    fold_df = fold_df.loc[ordered_genes]
    if len(row_clusters) != len(fold_df):
        if len(row_clusters) > len(fold_df):
            row_clusters = row_clusters[: len(fold_df)]
        else:
            row_clusters = row_clusters + ["unassigned"] * (len(fold_df) - len(row_clusters))
        block_breaks = []
        for idx in range(1, len(row_clusters)):
            if row_clusters[idx] != row_clusters[idx - 1]:
                block_breaks.append(idx)

    # ------------------------- 4. Apply lineage order if available ------------------------- #
    if lineage_order is not None:
        # --- Reorder strictly by lineage_order but keep only populations with DE results ---
        existing_pops = [p for p in lineage_order if p in fold_df.columns]
        missing_pops = [p for p in lineage_order if p not in fold_df.columns]
        print(f"[DEBUG] Existing pops for lineage order: {existing_pops}")
        #print(f"[DEBUG] Missing pops: {missing_pops}")
        if lineage_source:
            print(f"[INFO] Using lineage order from {lineage_source} (n={len(lineage_order)}).")
        else:
            print(f"[INFO] Using lineage order (n={len(lineage_order)}).")
        if missing_pops:
            print(f"[WARN] {len(missing_pops)} lineage populations missing from DE results: {missing_pops}")
        if existing_pops:
            fold_df = fold_df.reindex(columns=existing_pops)
            print(f"[INFO] Reordered heatmap columns by lineage_order (n={len(existing_pops)}).")
            #print(f"[DEBUG] fold_df columns after reorder: {list(fold_df.columns)}")
        else:
            print("[WARN] No matching lineage_order populations found in DE matrix; keeping default order.")
    else:
        print("[WARN] lineage_order not found in de_store or adata; keeping existing column order.")
        #print(f"[DEBUG] pop_order: {pop_order}")
        #print(f"[DEBUG] fold_df columns: {list(fold_df.columns)}")

    # ------------------------- 4. Export TSV ------------------------- #
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    heatmap_path = os.path.join(outdir, heatmap_tsv)

    # --- Ensure heatmap.tsv follows lineage_order exactly as the plotted heatmap ---
    if pop_order and all(p in fold_df.columns for p in pop_order):
        fold_df = fold_df.loc[:, [p for p in pop_order if p in fold_df.columns]]
    else:
        print("[WARN] lineage_order incomplete or missing columns; using current order.")
        print(f"[DEBUG] pop_order when writing TSV: {pop_order}")
        print(f"[DEBUG] fold_df columns when writing TSV: {list(fold_df.columns)}")

    heatmap_export = fold_df.copy()
    heatmap_export.index = [f"{c}:{g}" for c, g in zip(row_clusters, heatmap_export.index)]
    heatmap_export.columns = [f"{c}:{c}" for c in heatmap_export.columns]
    heatmap_export.to_csv(heatmap_path, sep="\t")
    print("[INFO] Wrote heatmap matrix TSV: {}".format(heatmap_path))


    # ------------------------- 5. Draw heatmap ------------------------- #
    go_terms_map = {}
    cluster_ranges = []
    if show_go_terms and goelite_payload and goelite_payload.get("runner") and goelite_payload.get("prepared"):
        runner = goelite_payload["runner"]
        prepared = goelite_payload["prepared"]
        cluster_gene_map = collections.OrderedDict()
        for gene, cluster in zip(ordered_genes, row_clusters):
            cluster_gene_map.setdefault(cluster, []).append(gene)

        current_cluster = None
        start = 0
        for idx, cluster in enumerate(row_clusters):
            if cluster != current_cluster:
                if current_cluster is not None:
                    cluster_ranges.append((current_cluster, start, idx))
                current_cluster = cluster
                start = idx
        if current_cluster is not None:
            cluster_ranges.append((current_cluster, start, len(row_clusters)))

        for cluster, genes in cluster_gene_map.items():
            if str(cluster) == "unassigned":
                continue
            if not genes:
                continue
            results = runner.run_prepared(genes, prepared, apply_prioritization=True)
            if not results:
                continue
            selected = [res for res in results if res.selected]
            if not selected:
                selected = results
            selected.sort(key=lambda res: (res.p_value, -abs(res.z_score)))
            term_hits = []
            for res in selected:
                node = runner.go_tree.get(res.term_id)
                name = node.name if node else ""
                if name:
                    term_hits.append((name, res.p_value))
            if term_hits:
                go_terms_map[str(cluster)] = term_hits

    def YellowBlackSky():
        cdict = {
            'red':   [(0.0, 0.0, 0.0), (0.5, 0.0, 0.1), (1.0, 1.0, 1.0)],
            'green': [(0.0, 0.0, 0.8), (0.5, 0.1, 0.0), (1.0, 1.0, 1.0)],
            'blue':  [(0.0, 0.0, 1.0), (0.5, 0.1, 0.0), (1.0, 0.0, 0.0)]
        }
        return LinearSegmentedColormap('YellowBlackSky', cdict)

    cmap = YellowBlackSky()
    vmax = max(1.0, np.log2(float(fc_thresh)) * 1.5)
    vmin = -vmax

    fixed_height = 4.2  # shorter layout to reduce heatmap height
    base_width = max(6.0, 0.25 * len(pop_order))
    ax_labels = None
    if show_go_terms:
        fig = plt.figure(figsize=(base_width, fixed_height))
        plot_top = 0.78
        plot_bottom = 0.10
        gs = fig.add_gridspec(
            1,
            3,
            width_ratios=[0.35, 0.50, 0.15],
            wspace=0.0,
            left=0.38,
            right=0.98,
            top=plot_top,
            bottom=plot_bottom,
        )
        ax_terms = fig.add_subplot(gs[0, 0], sharey=None)
        ax = fig.add_subplot(gs[0, 1], sharey=ax_terms)
        ax_labels = fig.add_subplot(gs[0, 2], sharey=ax)
    else:
        fig = plt.figure(figsize=(base_width, fixed_height))
        ax = plt.gca()
        ax_terms = None

    # --- FIX: remove white space caused by NaN folds in missing populations ---
    # Fill NaN with 0 (neutral fold) so heatmap renders continuously
    # Do NOT reorder by pop_order again; preserve lineage_order from above
    fold_df = fold_df.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    from matplotlib import colors
    from matplotlib.colors import ListedColormap
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    # --- Create normalization with custom scaling ---
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    column_clusters = list(fold_df.columns)

    # --- Draw the heatmap ---
    im = ax.imshow(fold_df.values, aspect="auto", cmap=cmap, norm=norm, interpolation="none")

    # --- Add vertical and horizontal white gridlines ---
    for x in range(1, fold_df.shape[1]):
        ax.axvline(x - 0.49, color="white", linewidth=0.39, alpha=0.9)

    for y in block_breaks:
        ax.axhline(y - 0.5, color="white", linewidth=0.35, alpha=0.9)

    ax.set_yticks([])
    ax.set_xticks(np.arange(fold_df.shape[1]))
    ax.set_xticklabels(list(fold_df.columns), rotation=45, ha="left", fontsize=5.5)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True, length=0, pad=6)
    if show_go_terms:
        fig.suptitle(
            "cellHarmony DE log2FC ({} vs {}) by {}".format(
                de_store["case_label"], de_store["control_label"], population_col
            ),
            x=0.5,
            y=0.98,
            fontsize=10,
        )
    else:
        ax.set_title(
            "cellHarmony DE log2FC ({} vs {}) by {}".format(
                de_store["case_label"], de_store["control_label"], population_col
            )
        )

    column_colors = [plt.get_cmap("tab20")(i % 20) for i in range(len(column_clusters))]
    column_cmap = ListedColormap(column_colors)
    column_to_id = {cluster: idx for idx, cluster in enumerate(column_clusters)}

    row_cluster_order = list(dict.fromkeys(row_clusters))
    row_colors = [plt.get_cmap("tab20b")(i % 20) for i in range(len(row_cluster_order))]
    row_cmap = ListedColormap(row_colors)
    row_to_id = {cluster: idx for idx, cluster in enumerate(row_cluster_order)}
    row_color_map = {cluster: row_colors[row_to_id[cluster]] for cluster in row_cluster_order}

    divider = make_axes_locatable(ax)
    ax_top = divider.append_axes("top", size="1.5%", pad=0.02)
    ax_left = divider.append_axes("left", size="3%", pad=0.015)

    col_ids = np.array([column_to_id[c] for c in column_clusters], dtype=int)[None, :]
    row_ids = np.array([row_to_id.get(c, 0) for c in row_clusters], dtype=int)[:, None]

    ax_top.imshow(col_ids, aspect="auto", cmap=column_cmap, interpolation="none")
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    ax_top.set_ylabel("")

    ax_left.imshow(row_ids, aspect="auto", cmap=row_cmap, interpolation="none")
    ax_left.set_xticks([])
    ax_left.set_yticks([])
    ax_left.set_xlabel("")
    if ax_terms is not None and ax_labels is not None:
        fig.canvas.draw()
        heatmap_bbox = ax.get_position()
        ax_terms.set_position([ax_terms.get_position().x0, heatmap_bbox.y0, ax_terms.get_position().width, heatmap_bbox.height])
        ax_labels.set_position([ax_labels.get_position().x0, heatmap_bbox.y0, ax_labels.get_position().width, heatmap_bbox.height])

    if ax_terms is not None:
        ax_terms.set_xlim(0, 1)
        ax_terms.set_ylim(ax.get_ylim())
        ax_terms.set_xticks([])
        ax_terms.set_yticks([])
        ax_terms.axis("off")
        ax_terms.patch.set_alpha(0.0)

        fig.canvas.draw()
        axis_height_px = ax.get_window_extent().height
        axis_width_px = ax_terms.get_window_extent().width
        row_height_px = axis_height_px / max(1, len(row_clusters))
        text_x = 0.96

        def _format_pvalue(pvalue):
            try:
                pval = float(pvalue)
            except Exception:
                return "p=NA"
            if pval < 1e-3:
                return "p={:.1e}".format(pval)
            return "p={:.3f}".format(pval)

        go_term_font_size = 4
        text_height_px = go_term_font_size * fig.dpi / 72.0 * 2.8
        min_rows_per_term = max(1, int(math.ceil(text_height_px / max(1.0, row_height_px))))

        cluster_entries = []
        for cluster, start, end in cluster_ranges:
            terms = go_terms_map.get(str(cluster), [])
            block_height = max(1, end - start)
            slot_count = max(1, int(block_height / min_rows_per_term))
            max_terms = min(slot_count, len(terms))
            if goelite_max_terms:
                max_terms = min(max_terms, goelite_max_terms)
            terms = terms[:max_terms] if terms else []
            cluster_entries.append(
                {
                    "cluster": cluster,
                    "start": start,
                    "end": end,
                    "block_height": block_height,
                    "terms": terms,
                    "positions": [],
                }
            )

        def _compute_positions(entry):
            terms = entry["terms"]
            if not terms:
                entry["positions"] = []
                return
            block_height = max(1, entry["end"] - entry["start"])
            available_rows = max(1, block_height - 1)
            total_height = (len(terms) - 1) * min_rows_per_term if len(terms) > 1 else 0
            top_offset = max(0.0, (available_rows - total_height) / 2.0)
            start_y = entry["start"] + 0.5 + top_offset
            positions = []
            for i in range(len(terms)):
                y = start_y + i * min_rows_per_term
                if y > entry["end"] - 0.5:
                    break
                positions.append(y)
            entry["positions"] = positions

        for entry in cluster_entries:
            _compute_positions(entry)

        for idx in range(1, len(cluster_entries)):
            prev = cluster_entries[idx - 1]
            curr = cluster_entries[idx]
            while prev["positions"] and curr["positions"]:
                prev_last = prev["positions"][-1]
                curr_first = curr["positions"][0]
                if (curr_first - prev_last) >= min_rows_per_term:
                    break
                if prev["block_height"] >= curr["block_height"]:
                    prev["terms"] = prev["terms"][:-1]
                    _compute_positions(prev)
                else:
                    curr["terms"] = curr["terms"][:-1]
                    _compute_positions(curr)

        for entry in cluster_entries:
            if not entry["positions"]:
                continue
            term_color = row_color_map.get(str(entry["cluster"]), "blue")
            for y, (term, pval) in zip(entry["positions"], entry["terms"]):
                term = "{} ({})".format(term, _format_pvalue(pval))
                ax_terms.text(
                    text_x,
                    y,
                    term,
                    transform=ax_terms.get_yaxis_transform(),
                    ha="right",
                    va="center",
                    fontsize=go_term_font_size,
                    color=term_color,
                    clip_on=False,
                )

    cax = inset_axes(
        ax,
        width="35%",
        height="3%",
        loc="lower center",
        bbox_to_anchor=(0.0, -0.06, 1.0, 1.0),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )
    cbar = plt.colorbar(im, cax=cax, orientation="horizontal")
    cbar.ax.set_xlim(vmin, vmax)
    cbar.set_label("log2 fold change", fontsize=6, labelpad=-2)
    cbar.set_ticks([])
    cbar.ax.text(-0.10, 0.32, f"{vmin:.1f}", ha="right", va="center", transform=cbar.ax.transAxes, fontsize=6)
    cbar.ax.text(1.10, 0.32, f"{vmax:.1f}", ha="left", va="center", transform=cbar.ax.transAxes, fontsize=6)
    cbar.outline.set_linewidth(0.5)

    if ax_labels is not None:
        ax_labels.set_xlim(0, 1)
        ax_labels.set_ylim(ax.get_ylim())
        ax_labels.set_xticks([])
        ax_labels.set_yticks([])
        ax_labels.axis("off")
        ax_labels.patch.set_alpha(0.0)

    first_gene_positions = {}
    for idx, cluster in enumerate(row_clusters):
        if cluster not in first_gene_positions:
            first_gene_positions[cluster] = idx
    fig.canvas.draw()
    if ax_labels is not None:
        label_axis = ax_labels
    else:
        label_axis = ax
    axis_height_px = label_axis.get_window_extent().height
    row_height_px = axis_height_px / max(1, len(row_clusters))
    one_point_rows = (fig.dpi / 72.0) / max(1.0, row_height_px)
    base_label_offset = 0.15
    gene_label_offset = base_label_offset + one_point_rows
    for cluster, pos in first_gene_positions.items():
        gene = fold_df.index[pos]
        if ax_labels is not None:
            label_x = 0.02
            label_y = pos + gene_label_offset
            ax_labels.text(
                label_x,
                label_y,
                gene,
                transform=ax_labels.get_yaxis_transform(),
                ha="left",
                va="center",
                fontsize=4,
                fontstyle="italic",
                clip_on=False,
            )
        else:
            ax.text(
                1.03,
                pos + gene_label_offset,
                gene,
                transform=ax.get_yaxis_transform(),
                ha="left",
                va="center",
                fontsize=6,
                fontstyle="italic",
                clip_on=True,
            )

    pdf_path = os.path.join(outdir, heatmap_png)

    # Save high-quality raster and vector outputs
    _suppress_fonttools_logs()
    bbox = None if show_go_terms else "tight"
    plt.savefig(pdf_path, dpi=600, bbox_inches=bbox, pad_inches=0.02)
    plt.savefig(pdf_path.replace(".png", ".pdf"), dpi=600, bbox_inches=bbox, pad_inches=0.02, transparent=True)
    plt.savefig(pdf_path.replace(".png", ".svg"), bbox_inches=bbox, transparent=True)

    plt.close(fig)
    print(f"[INFO] Wrote heatmap images: {pdf_path}, {pdf_path.replace('.png', '.pdf')}, {pdf_path.replace('.png', '.svg')}")


    return pdf_path, heatmap_path


def _normalize_goelite_species(species):
    if species is None:
        return None
    value = str(species).strip().lower()
    mapping = {
        "hs": "human",
        "human": "human",
        "homo_sapiens": "human",
        "homo sapiens": "human",
        "mm": "mouse",
        "mouse": "mouse",
        "mus_musculus": "mouse",
        "mus musculus": "mouse",
    }
    return mapping.get(value, value)

def _is_url(path):
    try:
        from urllib.parse import urlparse
        return urlparse(str(path)).scheme in ("http", "https", "ftp")
    except Exception:
        return False

def _download_goelite_resource(resource, download_dir, filename_hint=None):
    if resource is None or not _is_url(resource):
        return resource
    from urllib.parse import urlparse
    import urllib.request

    parsed = urlparse(str(resource))
    name = filename_hint or os.path.basename(parsed.path) or "goelite_resource"
    os.makedirs(download_dir, exist_ok=True)
    dest = os.path.join(download_dir, name)
    if os.path.isfile(dest) and os.path.getsize(dest) > 0:
        if str(dest).endswith(".gz"):
            try:
                with open(dest, "rb") as handle:
                    magic = handle.read(2)
                if magic != b"\x1f\x8b":
                    print(f"[WARN] GO-Elite cached file is not gzipped: {dest}")
                    os.remove(dest)
                else:
                    return dest
            except Exception:
                pass
        else:
            return dest

    try:
        if parsed.scheme in ("http", "https"):
            request = urllib.request.Request(
                resource,
                headers={"User-Agent": "AltAnalyze3-GoElite/1.0"}
            )
            with urllib.request.urlopen(request) as resp, open(dest, "wb") as out:
                out.write(resp.read())
        else:
            urllib.request.urlretrieve(resource, dest)
    except Exception as exc:
        print(f"[WARN] GO-Elite download failed for {resource}: {exc}")
        return None
    return dest

def _resolve_goelite_target(resource, download_dir, filename_hint=None):
    if resource is None:
        return None, None
    if _is_url(resource):
        try:
            from urllib.parse import urlparse
            parsed = urlparse(str(resource))
            name = filename_hint or os.path.basename(parsed.path) or "goelite_resource"
        except Exception:
            name = filename_hint or "goelite_resource"
        return os.path.join(download_dir, name), resource
    return resource, None

def _prompt_goelite_download(missing):
    if not missing:
        return True
    print("[INFO] GO-Elite resources missing locally:")
    for label, url, dest in missing:
        print(f"  - {label}: {url} -> {dest}")
    if not sys.stdin.isatty():
        print("[WARN] Non-interactive session; skipping GO-Elite download.")
        return False
    resp = input("Download GO-Elite resources now? [y/N]: ").strip().lower()
    return resp in ("y", "yes")

def run_goelite_for_clusters(de_store,
                             background_genes,
                             outdir,
                             comparison_tag,
                             species,
                             obo_path=None,
                             gaf_path=None,
                             cache_dir=None,
                             download_dir=None,
                             min_term_size=5,
                             max_term_size=2000):
    """
    Run GO-Elite enrichment per population cluster and write a combined TSV.
    """
    species_key = _normalize_goelite_species(species)
    if species_key is None:
        print("[INFO] GO-Elite skipped: no species provided.")
        return None

    from altanalyze3.components.goelite.resources import prepare_species_resources
    from altanalyze3.components.goelite.runner import GOEliteRunner, EnrichmentSettings

    per_pop = de_store.get("per_population_deg", {})
    if not per_pop:
        print("[WARN] GO-Elite skipped: no per-population DEG data found.")
        return None
    print(f"[INFO] GO-Elite: running enrichment for {len(per_pop)} populations...")

    if download_dir is None:
        download_dir = os.path.join(outdir, "_goelite_downloads")
    obo_path = _download_goelite_resource(obo_path, download_dir, "go-basic.obo")
    gaf_path = _download_goelite_resource(gaf_path, download_dir, "goa.gaf.gz")

    try:
        print("[INFO] GO-Elite: preparing resources...")
        parsed = prepare_species_resources(
            species_key,
            cache_dir=cache_dir,
            obo_path=obo_path,
            gaf_path=gaf_path,
        )
    except Exception as exc:
        print(f"[WARN] GO-Elite skipped: {exc}")
        return None
    settings = EnrichmentSettings(min_term_size=int(min_term_size), max_term_size=int(max_term_size))
    runner = GOEliteRunner(parsed, settings=settings)

    use_rawp = bool(de_store.get("use_rawp", False))
    alpha = float(de_store.get("alpha", 0.05))
    fc_thresh = float(de_store.get("fc_thresh", 1.2))
    log2_fc_thresh = np.log2(fc_thresh)

    bg = [str(g) for g in pd.Index(background_genes).dropna().unique().tolist()]
    if not bg:
        print("[WARN] GO-Elite skipped: empty background gene list.")
        return None
    prepared = runner.prepare_background(bg)
    payload = {"runner": runner, "prepared": prepared, "results": None}
    print(f"[INFO] GO-Elite: prepared {len(prepared.term_genes)} terms for testing.")

    rows = []
    go_pops = sorted(per_pop.keys())
    go_pbar = tqdm(
        total=len(go_pops),
        desc="GO-Elite per population",
        ncols=100,
        dynamic_ncols=True,
        position=0,
        leave=True,
        file=sys.stdout,
        disable=False,
    )
    for pop_name in go_pops:
        go_pbar.set_postfix_str(str(pop_name))
        df = per_pop[pop_name].copy()
        if df.empty or "log2fc" not in df.columns:
            go_pbar.update(1)
            continue

        sig_col = "pval" if use_rawp and "pval" in df.columns else "fdr"
        df["log2fc"] = pd.to_numeric(df["log2fc"], errors="coerce")
        df[sig_col] = pd.to_numeric(df[sig_col], errors="coerce")
        sig = (df[sig_col] < alpha) & (np.abs(df["log2fc"]) > log2_fc_thresh)
        genes = pd.Index(df.index.astype(str))[sig].unique().tolist()
        if len(genes) == 0:
            print(f"[INFO] GO-Elite: no significant genes for {pop_name}.")
            go_pbar.update(1)
            continue

        results = runner.run_prepared(genes, prepared, apply_prioritization=True)
        if not results:
            print(f"[INFO] GO-Elite: no enrichment results for {pop_name}.")
            go_pbar.update(1)
            continue

        query_map = {}
        for gene in genes:
            key = gene.upper()
            if key not in query_map:
                query_map[key] = gene
        query_upper = set(query_map.keys())

        for res in results:
            if not res.selected:
                continue
            term_genes = prepared.term_genes.get(res.term_id, set())
            overlap_upper = query_upper & term_genes
            overlap_genes = sorted(
                query_map[key] for key in overlap_upper if key in query_map
            )
            overlap_text = ",".join(overlap_genes)
            node = runner.go_tree.get(res.term_id)
            term_name = node.name if node else ""
            namespace = node.namespace if node else ""
            rows.append({
                "population": str(pop_name),
                "term_id": res.term_id,
                "term_name": term_name,
                "namespace": namespace,
                "p_value": res.p_value,
                "fdr": res.fdr,
                "z_score": res.z_score,
                "overlap": res.overlap,
                "term_genes": res.total_genes,
                "query_size": len(query_upper),
                "selected": res.selected,
                "overlap_genes": overlap_text,
            })
        go_pbar.update(1)
    go_pbar.close()

    if not rows:
        print("[INFO] GO-Elite produced no results across populations.")
        return payload

    os.makedirs(outdir, exist_ok=True)
    safe_tag = NetPerspective.safe_component(comparison_tag)
    out_path = os.path.join(outdir, f"GOElite_{safe_tag}.tsv")
    columns = [
        "population",
        "term_id",
        "term_name",
        "namespace",
        "p_value",
        "fdr",
        "z_score",
        "overlap",
        "term_genes",
        "query_size",
        "selected",
        "overlap_genes",
    ]
    out_df = pd.DataFrame(rows, columns=columns)
    out_df.sort_values(["population", "fdr", "p_value"], inplace=True)
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[INFO] GO-Elite results written to: {out_path}")
    payload["results"] = out_df
    return payload

# ------------------------- Step 6: differentials h5ad ---------------------- #

def write_differentials_only_h5ad(adata, de_store, out_path):
    detailed = de_store["detailed_deg"]
    if detailed.empty:
        print("[WARN] No DEGs; writing a copy of input with DE container only.")
        adx = adata.copy()
        _sanitize_de_store(de_store)
        if "cellHarmony_DE" in adx.uns:
            del adx.uns["cellHarmony_DE"]
        adx.uns["cellHarmony_DE"] = de_store
        if hasattr(adx, "var") and "_index" in adx.var.columns:
            adx.var = adx.var.drop(columns=["_index"])
        if adx.raw is not None and "_index" in adx.raw.var.columns:
            adx.raw.var = adx.raw.var.drop(columns=["_index"])
        adx.write(out_path)
        print("[INFO] Wrote differentials h5ad: {}".format(out_path))
        return

    deg_union = sorted(detailed["gene"].astype(str).unique().tolist())
    # restrict to union of DEGs (var filter); retain same obs
    # map deg gene names to var_names; skip if absent
    var_names = adata.var_names.astype(str)
    keep_mask = np.isin(var_names, np.array(deg_union, dtype=str))
    if keep_mask.sum() == 0:
        print("[WARN] None of the DEG genes match var_names; writing input with DE container only.")
        adx = adata.copy()
        _sanitize_de_store(de_store)
        if "cellHarmony_DE" in adx.uns:
            del adx.uns["cellHarmony_DE"]
        adx.uns["cellHarmony_DE"] = de_store
        if hasattr(adx, "var") and "_index" in adx.var.columns:
            adx.var = adx.var.drop(columns=["_index"])
        if adx.raw is not None and "_index" in adx.raw.var.columns:
            adx.raw.var = adx.raw.var.drop(columns=["_index"])
        adx.write(out_path)
        print("[INFO] Wrote differentials h5ad: {}".format(out_path))
        return

    adx = adata[:, keep_mask].copy()
    _sanitize_de_store(de_store)
    if adx.raw is not None:
        adx.raw = adx.raw[:, keep_mask]
    if "cellHarmony_DE" in adx.uns:
        del adx.uns["cellHarmony_DE"]
    adx.uns["cellHarmony_DE"] = de_store

    # Drop reserved column names before writing
    if hasattr(adx, "var") and "_index" in adx.var.columns:
        adx.var = adx.var.drop(columns=["_index"])
    if adx.raw is not None and "_index" in adx.raw.var.columns:
        adx.raw.var = adx.raw.var.drop(columns=["_index"])

    adx.write(out_path)
    print("[INFO] Wrote differentials h5ad: {}".format(out_path))


def _sanitize_de_store(de_store):
    assigned = de_store.get("assigned_groups")
    if isinstance(assigned, pd.DataFrame):
        for col in assigned.columns:
            assigned[col] = assigned[col].apply(lambda v: str(v) if not isinstance(v, str) else v)
        de_store["assigned_groups"] = assigned

    detailed = de_store.get("detailed_deg")
    if isinstance(detailed, pd.DataFrame):
        for col in detailed.columns:
            detailed[col] = detailed[col].apply(lambda v: str(v) if not isinstance(v, str) else v)
        de_store["detailed_deg"] = detailed

# --------------------------------- CLI ------------------------------------ #

def parse_comparisons_arg(comp_str):
    """
    Parse comparisons string like:
      "breast-cancer|wild-type;A|B"
    into list of tuples: [(case, control), ...]
    """
    comps = []
    if comp_str is None or len(comp_str.strip()) == 0:
        return comps
    parts = [p for p in comp_str.split(";") if len(p.strip()) > 0]
    for p in parts:
        bits = [b for b in p.split("|") if len(b.strip()) > 0]
        if len(bits) != 2:
            print("[ERROR] Invalid comparison '{}'; expected 'CASE|CONTROL'.".format(p), file=sys.stderr)
            sys.exit(1)
        comps.append((bits[0], bits[1]))
    return comps

def main():
    ap = argparse.ArgumentParser(description="cellHarmony differential analysis and (optional) pseudobulk generation.")
    ap.add_argument("--h5ad", required=True, help="Input h5ad with aligned annotations")
    ap.add_argument("--covariates", required=True, help="TSV mapping libraries to samples and covariate")
    ap.add_argument("--library_col", default="Library", help="obs column for library IDs (default: Library)")
    ap.add_argument("--sample_col", default="Sample", help="column in covariates for sample IDs (default: Sample)")
    ap.add_argument("--covariate_col", default="Condition", help="column in covariates for the covariate/condition (default: Condition)")
    ap.add_argument("--population_col", required=True, help="obs column with projected cell populations (matches reference ref_name)")
    ap.add_argument("--comparisons", required=True, help="Semicolon-separated CASE|CONTROL pairs, e.g. 'breast-cancer|wild-type;A|B'")
    ap.add_argument("--method", choices=["wilcoxon", "t-test", "t-test_overestim_var", "logreg"], default="wilcoxon",
                    help="Scanpy rank_genes_groups method (default: wilcoxon)")
    ap.add_argument("--alpha", type=float, default=0.05, help="FDR threshold (default: 0.05)")
    ap.add_argument("--fc", type=float, default=1.2, help="Fold-change threshold (absolute, default: 1.2)")
    ap.add_argument("--min_cells_per_group", type=int, default=20, help="Minimum cells per group per population (default: 20; auto-relaxed to 4 if total<200)")
    ap.add_argument("--make_pseudobulk", action="store_true", help="If set, compute pseudobulks per (population×sample)")
    ap.add_argument("--pseudobulk_min_cells", type=int, default=10, help="Minimum cells per pseudobulk group (default: 10)")
    ap.add_argument("--use_rawp", action="store_true", help="Use raw p-values instead of FDR for significance filtering")
    ap.add_argument("--goelite_species", default=None, help="Run GO-Elite with species (e.g. human or mouse)")
    ap.add_argument("--goelite_obo", default=None, help="GO .obo path or URL (optional)")
    ap.add_argument("--goelite_gaf", default=None, help="GOA .gaf path or URL (optional)")
    ap.add_argument("--goelite_cache_dir", default=None, help="Optional GO-Elite cache directory")
    ap.add_argument("--goelite_min_term_size", type=int, default=5, help="GO-Elite minimum term size (default: 5)")
    ap.add_argument("--goelite_max_term_size", type=int, default=2000, help="GO-Elite maximum term size (default: 2000)")
    ap.add_argument("--outdir", default="cellHarmony_DE_out", help="Output directory (default: cellHarmony_DE_out)")
    ap.add_argument("--skip_grn", action="store_true", help="Skip interaction network (GRN) generation")
    args = ap.parse_args()

    comps = parse_comparisons_arg(args.comparisons)
    if len(comps) == 0:
        print("[ERROR] No valid comparisons parsed from --comparisons.", file=sys.stderr)
        sys.exit(1)

    base_outdir = args.outdir
    os.makedirs(base_outdir, exist_ok=True)
    heatmap_dir = os.path.join(base_outdir, "heatmaps")
    deg_dir = os.path.join(base_outdir, "DEGs")
    pseudobulk_dir = os.path.join(base_outdir, "pseudobulk")
    cellfreq_dir = os.path.join(base_outdir, "cell-frequency")
    goelite_dir = os.path.join(base_outdir, "GeneSetEnrichment")
    log_dir = os.path.join(base_outdir, "logs")

    for path in (heatmap_dir, deg_dir, log_dir, cellfreq_dir):
        os.makedirs(path, exist_ok=True)
    if args.make_pseudobulk:
        os.makedirs(pseudobulk_dir, exist_ok=True)
    if args.goelite_species:
        os.makedirs(goelite_dir, exist_ok=True)
        from altanalyze3.components.goelite.resources import resolve_species_cache_dir
        normalized_species = _normalize_goelite_species(args.goelite_species)
        goelite_species_dir = resolve_species_cache_dir(normalized_species, args.goelite_cache_dir)
        goelite_download_dir = os.path.join(str(goelite_species_dir), "downloads")
        obo_local, obo_url = _resolve_goelite_target(args.goelite_obo, goelite_download_dir, "go-basic.obo")
        gaf_local, gaf_url = _resolve_goelite_target(args.goelite_gaf, goelite_download_dir, "goa.gaf.gz")
        missing = []
        if obo_url and (not os.path.isfile(obo_local) or os.path.getsize(obo_local) == 0):
            missing.append(("obo", obo_url, obo_local))
        if gaf_url and (not os.path.isfile(gaf_local) or os.path.getsize(gaf_local) == 0):
            missing.append(("gaf", gaf_url, gaf_local))
        if args.goelite_obo and not obo_url and not os.path.isfile(args.goelite_obo):
            print(f"[WARN] GO-Elite OBO path not found: {args.goelite_obo}")
            args.goelite_species = None
        if args.goelite_gaf and not gaf_url and not os.path.isfile(args.goelite_gaf):
            print(f"[WARN] GO-Elite GAF path not found: {args.goelite_gaf}")
            args.goelite_species = None
        if args.goelite_species and missing:
            if _prompt_goelite_download(missing):
                for label, url, dest in missing:
                    _download_goelite_resource(url, goelite_download_dir, os.path.basename(dest))
            else:
                args.goelite_species = None
        if args.goelite_species and not args.goelite_obo and not args.goelite_gaf:
            species_key = _normalize_goelite_species(args.goelite_species)
            try:
                from altanalyze3.components.goelite.resources import load_cached_resources, SPECIES_CONFIG, GO_ONTOLOGY_URL
                load_cached_resources(species_key, cache_dir=args.goelite_cache_dir)
            except Exception:
                from urllib.parse import urlparse
                default_obo_url = GO_ONTOLOGY_URL
                default_gaf_url = None
                if species_key in SPECIES_CONFIG:
                    default_gaf_url = SPECIES_CONFIG[species_key]["gaf_url"]
                defaults = []
                if default_obo_url:
                    obo_name = os.path.basename(urlparse(default_obo_url).path) or "go-basic.obo"
                    defaults.append(("obo", default_obo_url, os.path.join(goelite_download_dir, obo_name)))
                if default_gaf_url:
                    gaf_name = os.path.basename(urlparse(default_gaf_url).path) or "goa.gaf.gz"
                    defaults.append(("gaf", default_gaf_url, os.path.join(goelite_download_dir, gaf_name)))
                if defaults:
                    if _prompt_goelite_download(defaults):
                        for label, url, dest in defaults:
                            downloaded = _download_goelite_resource(url, goelite_download_dir, os.path.basename(dest))
                            if label == "obo":
                                args.goelite_obo = downloaded
                            elif label == "gaf":
                                args.goelite_gaf = downloaded
                    else:
                        args.goelite_species = None
        args.goelite_obo = obo_local if obo_url else args.goelite_obo
        args.goelite_gaf = gaf_local if gaf_url else args.goelite_gaf

    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    log_path = os.path.join(log_dir, f"cellHarmony-differential_{timestamp}.log")
    log_file = open(log_path, "w")
    stdout_orig, stderr_orig = sys.stdout, sys.stderr
    sys.stdout = Tee(stdout_orig, log_file)
    sys.stderr = Tee(stderr_orig, log_file)

    try:
        log_file.write("# cellHarmony differential run parameters\n")
        log_file.write(f"command = {' '.join(sys.argv)}\n")
        for key, value in sorted(vars(args).items()):
            log_file.write(f"{key} = {value}\n")
        log_file.write("\n")
        log_file.flush()

        interactions_df = None
        interaction_root = os.path.join(base_outdir, "interaction-plots")
        if not args.skip_grn:
            try:
                interactions_df = NetPerspective.load_interaction_data()
            except Exception as ex:
                print(f"[WARN] Unable to load interaction data for GRN export: {ex}")
                args.skip_grn = True
            else:
                os.makedirs(interaction_root, exist_ok=True)

        # Step 1: load + merge covariates
        print("[INFO] Loading h5ad and merging covariates...")
        adata = load_and_merge_covariates(args.h5ad, args.covariates, args.library_col, args.sample_col, args.covariate_col)

        # Optional Step 3: pseudobulk
        if args.make_pseudobulk:
            print("[INFO] Computing pseudobulks per ({} × {}).".format(args.population_col, args.sample_col))
            _, pb_h5ad = compute_pseudobulk_per_population(
                adata,
                population_col=args.population_col,
                sample_col=args.sample_col,
                covariate_col=args.covariate_col,
                min_cells=int(args.pseudobulk_min_cells),
                outdir=pseudobulk_dir
            )
            adata = ad.read_h5ad(pb_h5ad)

        # Step 4–6 per comparison
        for case_label, control_label in comps:
            tag = "{}_vs_{}".format(str(case_label), str(control_label)).replace(" ", "_")
            print("[INFO] Differential analysis: {} vs {}.".format(case_label, control_label))

            all_pops = sorted(adata.obs[args.population_col].unique())
            print(f"[INFO] Running differential expression across {len(all_pops)} populations...")

            pbar = tqdm(
                total=len(all_pops),
                desc="Running DE per population",
                ncols=100,
                dynamic_ncols=True,
                position=0,
                leave=True
            )

            def progress_callback():
                pbar.update(1)
                pbar.refresh()

            de_store = run_de_for_comparisons(
                adata=adata,
                population_col=args.population_col,
                covariate_col=args.covariate_col,
                case_label=case_label,
                control_label=control_label,
                method=args.method,
                alpha=float(args.alpha),
                fc_thresh=float(args.fc),
                min_cells_per_group=int(args.min_cells_per_group),
                use_rawp=args.use_rawp,
                progress_callback=progress_callback
            )

            pbar.close()
            sys.stdout.flush()

            pop_order = _extract_lineage_order(adata, args.population_col)
            if pop_order:
                print(f"[INFO] Using lineage order from h5ad (n={len(pop_order)} populations).")
            else:
                detailed = de_store.get("detailed_deg")
                if detailed is not None and "population_or_pattern" in detailed.columns:
                    pop_order = sorted(detailed["population_or_pattern"].astype(str).unique().tolist())
                else:
                    pop_order = []
                if hasattr(adata, "uns") and "lineage_order" in adata.uns:
                    print(f"[WARN] lineage_order present but mismatched with '{args.population_col}'; using alphabetical order.")
                else:
                    print("[WARN] 'lineage_order' not found in h5ad; using alphabetical order.")

            # ----------------------------- Cell frequency plots ----------------------------- #
            freq_subdir = os.path.join(cellfreq_dir, NetPerspective.safe_component(tag))
            os.makedirs(freq_subdir, exist_ok=True)

            freq_obs = adata.obs.loc[
                adata.obs[args.covariate_col].isin([case_label, control_label]),
                [args.population_col, args.covariate_col]
            ].dropna()

            if freq_obs.empty:
                print("[WARN] Skipping cell frequency plots: no cells found for comparison subset.")
            else:
                freq_obs[args.population_col] = freq_obs[args.population_col].astype(str)
                freq_obs[args.covariate_col] = freq_obs[args.covariate_col].astype(str)

                unique_pops = list(freq_obs[args.population_col].unique())
                if pop_order:
                    freq_order = [p for p in pop_order if p in unique_pops]
                    freq_order.extend([p for p in unique_pops if p not in freq_order])
                else:
                    freq_order = sorted(unique_pops)

                counts = (
                    freq_obs.groupby([args.population_col, args.covariate_col])
                    .size()
                    .unstack(fill_value=0)
                )

                case_name = str(case_label)
                control_name = str(control_label)
                for col in (control_name, case_name):
                    if col not in counts.columns:
                        counts[col] = 0

                counts = counts.reindex(freq_order).fillna(0)
                counts = counts.loc[counts.sum(axis=1) > 0]

                if counts.empty:
                    print("[WARN] Skipping cell frequency plots: all populations have zero cells for comparison.")
                else:
                    column_order = [control_name, case_name]
                    counts = counts.loc[:, [c for c in column_order if c in counts.columns]]

                    cluster_totals = counts.sum(axis=1).replace(0, np.nan)
                    percent_by_cluster = counts.div(cluster_totals, axis=0).fillna(0) * 100.0

                    group_totals = counts.sum(axis=0).replace(0, np.nan)
                    percent_by_group = counts.div(group_totals, axis=1).fillna(0) * 100.0

                    comp_component = NetPerspective.safe_component(tag)
                    pop_component = NetPerspective.safe_component(args.population_col)
                    stacked_path = os.path.join(freq_subdir, f"{comp_component}_stacked_{pop_component}.pdf")
                    grouped_path = os.path.join(freq_subdir, f"{comp_component}_by_condition_{pop_component}.pdf")

                    control_color = "skyblue"
                    case_color = "lightcoral"

                    fig1, ax1 = plt.subplots(figsize=(8, max(2.0, 0.4 * len(percent_by_cluster) + 2)))
                    percent_by_cluster.plot(
                        kind="barh",
                        stacked=True,
                        ax=ax1,
                        color=[control_color, case_color],
                        width=0.8,
                        edgecolor="black",
                        linewidth=0.2,
                    )
                    ax1.set_xlabel("Percent of cells per population")
                    ax1.set_ylabel(args.population_col)
                    ax1.set_title(f"Cell composition by {args.population_col}: {case_name} vs {control_name}", fontsize=12)
                    ax1.legend(
                        title="Condition",
                        bbox_to_anchor=(1.05, 1),
                        loc="upper left",
                        borderaxespad=0.0,
                        frameon=False,
                        fontsize=8,
                        title_fontsize=9,
                    )
                    plt.tight_layout(rect=[0, 0, 0.85, 1])
                    fig1.savefig(stacked_path, bbox_inches="tight")
                    plt.close(fig1)
                    print(f"[INFO] Wrote cell frequency stacked bar chart: {stacked_path}")

                    fig2, ax2 = plt.subplots(figsize=(8, max(2.0, 0.4 * len(percent_by_group) + 2)))
                    y_positions = np.arange(len(percent_by_group.index))
                    bar_height = 0.35

                    ax2.barh(
                        y_positions - bar_height / 2,
                        percent_by_group[control_name],
                        height=bar_height,
                        color=control_color,
                        edgecolor="black",
                        linewidth=0.2,
                        label=control_name,
                    )
                    ax2.barh(
                        y_positions + bar_height / 2,
                        percent_by_group[case_name],
                        height=bar_height,
                        color=case_color,
                        edgecolor="black",
                        linewidth=0.2,
                        label=case_name,
                    )
                    ax2.set_yticks(y_positions)
                    ax2.set_yticklabels(percent_by_group.index)
                    ax2.invert_yaxis()
                    ax2.set_xlabel("Percent of condition cells")
                    ax2.set_ylabel(args.population_col)
                    ax2.set_title(f"Cell frequency by condition: {case_name} vs {control_name}", fontsize=12)
                    ax2.legend(loc="best", frameon=False)
                    plt.tight_layout()
                    fig2.savefig(grouped_path, bbox_inches="tight")
                    plt.close(fig2)
                    print(f"[INFO] Wrote cell frequency comparison chart: {grouped_path}")

                    freq_tsv_path = os.path.join(freq_subdir, f"{comp_component}_cell_frequencies_{pop_component}.tsv")
                    counts_long = counts.reset_index().melt(
                        id_vars=args.population_col,
                        var_name="condition",
                        value_name="cell_count",
                    )
                    pct_group_long = percent_by_group.reset_index().melt(
                        id_vars=args.population_col,
                        var_name="condition",
                        value_name="percent_of_condition",
                    )
                    pct_cluster_long = percent_by_cluster.reset_index().melt(
                        id_vars=args.population_col,
                        var_name="condition",
                        value_name="percent_of_population",
                    )
                    freq_table = counts_long.merge(
                        pct_group_long,
                        on=[args.population_col, "condition"],
                        how="left",
                    ).merge(
                        pct_cluster_long,
                        on=[args.population_col, "condition"],
                        how="left",
                    )
                    freq_table.sort_values(
                        by=[args.population_col, "condition"],
                        inplace=True,
                    )
                    freq_table.to_csv(freq_tsv_path, sep="\t", index=False, float_format="%.4f")
                    print(f"[INFO] Wrote cell frequency table: {freq_tsv_path}")

            goelite_payload = None
            if args.goelite_species:
                background_genes = adata.var_names.astype(str)
                if "gene_symbols" in adata.var.columns:
                    background_genes = adata.var["gene_symbols"].astype(str)
                elif "features" in adata.var.columns:
                    background_genes = adata.var["features"].astype(str)

                goelite_subdir = os.path.join(goelite_dir, NetPerspective.safe_component(tag))
                print(f"[INFO] GO-Elite: starting for {tag} (output: {goelite_subdir})")
                goelite_payload = run_goelite_for_clusters(
                    de_store=de_store,
                    background_genes=background_genes,
                    outdir=goelite_subdir,
                    comparison_tag=tag,
                    species=args.goelite_species,
                    obo_path=args.goelite_obo,
                    gaf_path=args.goelite_gaf,
                    cache_dir=args.goelite_cache_dir,
                    download_dir=goelite_download_dir,
                    min_term_size=args.goelite_min_term_size,
                    max_term_size=args.goelite_max_term_size,
                )

            heat_png = "heatmap_{}_by_{}.pdf".format(tag, args.population_col)
            heat_tsv = "heatmap_{}_by_{}.tsv".format(tag, args.population_col)
            build_fixed_order_heatmap(
                de_store,
                heatmap_dir,
                args.population_col,
                float(args.fc),
                heatmap_png=heat_png,
                heatmap_tsv=heat_tsv,
                goelite_payload=goelite_payload,
                show_go_terms=bool(goelite_payload),
            )

            out_h5ad = os.path.join(deg_dir, "differentials_only_{}.h5ad".format(tag))
            write_differentials_only_h5ad(adata, de_store, out_h5ad)

            det = de_store["detailed_deg"]
            summ = de_store["summary_per_population"]
            assign = de_store["assigned_groups"]
            po = de_store["pooled_overall"]
            cr = de_store["coreg_pooled"]

            det_path = os.path.join(deg_dir, "DEG_detailed_{}.tsv".format(tag))
            summ_path = os.path.join(deg_dir, "DEG_summary_{}.tsv".format(tag))
            assign_path = os.path.join(deg_dir, "DEG_assigned_groups_{}.tsv".format(tag))
            po_path = os.path.join(deg_dir, "DEG_pooled_overall_{}.tsv".format(tag))
            cr_path = os.path.join(deg_dir, "DEG_coreg_pooled_{}.tsv".format(tag))

            if not det.empty:
                det.to_csv(det_path, sep="\t", index=False)
                print("[INFO] Wrote {}".format(det_path))
            if not summ.empty:
                summ.to_csv(summ_path, sep="\t", index=False)
                print("[INFO] Wrote {}".format(summ_path))
            if not assign.empty:
                assign.to_csv(assign_path, sep="\t", index=False)
                print("[INFO] Wrote {}".format(assign_path))
            if po is not None and not po.empty:
                po.to_csv(po_path, sep="\t")
                print("[INFO] Wrote {}".format(po_path))
            if cr is not None and not cr.empty:
                cr.to_csv(cr_path, sep="\t", index=False)
                print("[INFO] Wrote {}".format(cr_path))

            if not args.skip_grn and interactions_df is not None:
                detailed = de_store.get("detailed_deg")
                if detailed is None or detailed.empty:
                    if diagnostic_report:
                        print("[DEBUG] No detailed DEG table available for GRN export.")
                else:
                    per_pop_results = {
                        str(pop): grp.copy()
                        for pop, grp in detailed.groupby("population")
                    }

                    if not per_pop_results:
                        if diagnostic_report:
                            print("[DEBUG] No per-population DEG entries found for GRN export.")
                    else:
                        comparison_dir = os.path.join(interaction_root, NetPerspective.safe_component(tag))
                        os.makedirs(comparison_dir, exist_ok=True)
                        use_rawp = bool(de_store.get("use_rawp", False))

                        for pop_name, pop_df in per_pop_results.items():
                            if pop_df is None or pop_df.empty:
                                if diagnostic_report:
                                    print(f"[DEBUG] GRN skip {pop_name}: empty per-pop DEG table.")
                                continue

                            stats_df = pop_df.copy()
                            stats_df = stats_df.reset_index(drop=True)
                            stats_df["gene"] = stats_df["gene"].astype(str)
                            stats_df["log2fc"] = pd.to_numeric(stats_df.get("log2fc"), errors="coerce")

                            significance_column = None
                            if use_rawp and "pval" in stats_df.columns and stats_df["pval"].notna().any():
                                significance_column = "pval"
                            elif "fdr" in stats_df.columns:
                                significance_column = "fdr"

                            keep_columns = {"gene", "log2fc", "fdr", "pval"}
                            if significance_column:
                                keep_columns.add(significance_column)
                            keep_columns = [col for col in keep_columns if col in stats_df.columns]

                            selected = (
                                stats_df.loc[:, keep_columns]
                                .dropna(subset=["gene", "log2fc"])
                                .drop_duplicates(subset=["gene"])
                            )

                            if diagnostic_report:
                                print(
                                    f"[DEBUG] GRN candidate genes for {pop_name}: total={len(stats_df)}, "
                                    f"after_cleanup={selected['gene'].nunique()}"
                                )

                            if selected.empty or selected["gene"].nunique() < 2:
                                if diagnostic_report:
                                    print(
                                        f"[DEBUG] GRN skip {pop_name}: insufficient unique genes (n={selected['gene'].nunique()})."
                                    )
                                continue

                            pop_component = NetPerspective.safe_component(pop_name)
                            output_prefix = os.path.join(comparison_dir, pop_component)
                            try:
                                outputs = NetPerspective.generate_network_for_genes(
                                    selected,
                                    interactions_df,
                                    output_prefix,
                                    gene_column="gene",
                                    fold_change_column="log2fc",
                                    pval_column=significance_column if significance_column in selected.columns else None,
                                    max_genes=None,
                                )
                                print(
                                    f"[INFO] Wrote interaction network for {pop_name} ({case_label} vs {control_label}): {outputs[0]}"
                                )
                            except NetPerspective.NetworkGenerationError as ex:
                                if diagnostic_report:
                                    print(f"[DEBUG] Interaction network skipped for {pop_name}: {ex}")
                                else:
                                    print(f"[WARN] No interaction edges found for {pop_name}; skipping network plot.")
                            except ImportError as ex:
                                print(f"[WARN] Interaction plots disabled (missing dependency): {ex}")
                            except Exception as ex:
                                print(f"[WARN] Interaction network failed for {pop_name}: {ex}")
                    args.skip_grn = True
                    break

        print("[INFO] Completed.")
    finally:
        sys.stdout = stdout_orig
        sys.stderr = stderr_orig
        log_file.close()
        stdout_orig.write(f"[INFO] cellHarmony differential log written to {log_path}\n")
        stdout_orig.flush()


if __name__ == "__main__":
    main()
