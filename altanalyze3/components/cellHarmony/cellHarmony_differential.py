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
  - write_differentials_only_h5ad(...)
  - main() with CLI

Notes:
  - Requires: scanpy, anndata, numpy, pandas, matplotlib
  - Optional: rpy2 (not used by default). We default to Scanpy methods.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sps
from tqdm import tqdm
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
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'

diagnostic_report = False
# ------------------------------- Utilities -------------------------------- #

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
        colnames.append(gid)

    if len(rows) == 0:
        print("[ERROR] No pseudobulk groups passed threshold.", file=sys.stderr)
        sys.exit(1)

    M = np.vstack(rows).T  # genes × groups
    df = pd.DataFrame(M, index=adata.var_names.astype(str), columns=colnames)
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
    # Preserve lineage hierarchy if present in input AnnData
    if "lineage_order" in adata.uns:
        pb_adata.uns["lineage_order"] = adata.uns["lineage_order"]
        print(f"[INFO] Preserved lineage_order from input (n={len(adata.uns['lineage_order'])}).")
    else:
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
    # Run Scanpy DE and return names, pvals_adj, log2fc as Series for the case vs control
    sc.tl.rank_genes_groups(two_group_adata,
                            groupby=groupby,
                            groups=[case_label],
                            reference=control_label,
                            method=method)
    rg = two_group_adata.uns["rank_genes_groups"]
    names = pd.Index(rg["names"][case_label]).astype(str)
    pvals_adj = pd.Series(rg["pvals_adj"][case_label], index=names, name="fdr").astype(float)
    logfc_ln = pd.Series(rg["logfoldchanges"][case_label], index=names, name="logfc_ln").astype(float)
    logfc = pd.Series(_ln_to_log2(logfc_ln.values), index=names, name="log2fc")
    return names, pvals_adj, logfc


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
            names  = df["gene"]
            fdr    = df["fdr"]
            log2fc = df["log2fc"]
            pvals  = df["pval"] if "pval" in df.columns else pd.Series(np.nan, index=names)

            # --- capture full log2FC vector for this population ---
            if "full_log2fc" in two.uns and str(pop) in two.uns["full_log2fc"]:
                all_fold_values[str(pop)] = two.uns["full_log2fc"][str(pop)]


        else:
            # Default Scanpy DE method
            try:
                names, fdr, log2fc = _rank_genes_scanpy(two, covariate_col, case_label, control_label, method)
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

        # --- Harmonize identifiers between DE results and pseudobulk matrix ---
        if corrected_fc is not None:
            df_names = pd.Index(names).astype(str).str.upper().str.replace(r'\.\d+$', '', regex=True)
            cf_names = pd.Index(corrected_fc.index).astype(str).str.upper().str.replace(r'\.\d+$', '', regex=True)
            overlap = len(set(df_names) & set(cf_names))

            if diagnostic_report:
                print(f"[DEBUG] Harmonized ID overlap: {overlap} / {len(df_names)} genes shared with corrected_fc index")
                unmatched = list(set(df_names) - set(cf_names))
                print(f"[DEBUG] Example unmatched DE identifiers (first 10): {unmatched[:10]}")
                matched = list(set(df_names) & set(cf_names))
                print(f"[DEBUG] Example matched identifiers (first 10): {matched[:10]}")
                print(f"[DEBUG] corrected_fc index example (first 10): {list(corrected_fc.index[:10])}")
        else:
            if diagnostic_report:
                print("[WARN] corrected_fc not yet defined during harmonization.")

        df = pd.DataFrame({
            "gene": names.values,
            "log2fc": log2fc.values,
            "fdr": fdr.values
        })

        # Harmonize gene identifiers if var includes gene symbols
        if "gene_symbols" in pop_data.var.columns:
            id_map = pop_data.var["gene_symbols"].to_dict()
            df["gene"] = df["gene"].map(lambda g: id_map.get(g, g))

        df["population"] = pop
        df["case_label"] = case_label
        df["control_label"] = control_label
        df["n_case"] = n_case
        df["n_control"] = n_ctrl

        # thresholding
        if use_rawp:
            # use raw p-values from the DE frame; ensure alignment to 'names'
            p_base = pvals.values
            sig = (np.abs(log2fc.values) > np.log2(float(fc_thresh))) & (p_base < float(alpha))
            if diagnostic_report:
                print(f"[INFO] Using raw p-values for significance (alpha={alpha})")
        else:
            p_base = fdr.values
            sig = (np.abs(log2fc.values) > np.log2(float(fc_thresh))) & (p_base < float(alpha))

        deg_names = names[sig]
        n_deg = int(sig.sum())

        # collect full stats table for this population
        df = pd.DataFrame({
            "gene": names.values,
            "log2fc": log2fc.values,
            "fdr": fdr.values,
            "pval": pvals if "pvals" in locals() else np.full_like(fdr.values, np.nan),
            "population": pop,
            "case_label": case_label,
            "control_label": control_label,
            "n_case": n_case,
            "n_control": n_ctrl,
        })


        # --- Correct fold computation: true mean-based log2 fold (per gene across cells) ---
        expr = pop_data.to_df()  # rows = cells, columns = genes
        cond = pop_data.obs[covariate_col].astype(str)
        cond = cond.reindex(expr.index)

        # compute mean expression per condition
        mean_case = np.log2(np.mean(2 ** expr[cond == case_label], axis=0) + 1e-9)
        mean_ctrl = np.log2(np.mean(2 ** expr[cond == control_label], axis=0) + 1e-9)
        corrected_fc = mean_case - mean_ctrl
        corrected_fc.index.name = "gene"

        # overwrite Scanpy’s log2fc column with corrected values
        df["log2fc"] = df["gene"].map(corrected_fc)

        if diagnostic_report:
            # --- DEBUGGING: PROSER2 in AT1 ---
            if str(pop) == "AT1" and "PROSER2" in df["gene"].values:
                row = df.loc[df["gene"] == "PROSER2"].iloc[0]
                print("\n[DEBUG] PROSER2 in AT1:")
                print(f"  log2FC: {row['log2fc']:.3f}")
                print(f"  FDR: {row['fdr']:.3e}")
                print(f"  case/control: {row['case']} vs {row['control']}")
                print(f"  n_case={row['n_case']}, n_control={row['n_control']}")
                print(f"  mean_case={mean_case.get('PROSER2', np.nan):.3f}, mean_ctrl={mean_ctrl.get('PROSER2', np.nan):.3f}")

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
        columns=["gene", "population", "log2fc", "fdr", "case", "control", "n_case", "n_control"]
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
        try:
            names_c, fdr_c, log2fc_c = _rank_genes_scanpy(pooled_pat, covariate_col, case_label, control_label, method)
            dfc = pd.DataFrame({"gene": names_c.values, "log2fc": log2fc_c.values, "fdr": fdr_c.values})

            # ----- Construct human-readable pattern label -----
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

            coreg_name = _pattern_to_label(pat, pop_order)

            # ----- Perform pooled DE for this co-regulation pattern -----
            if method == "scanpy":
                names_c, fdr_c, log2fc_c = _rank_genes_scanpy(pooled_pat, covariate_col, case_label, control_label, method)
                dfc = pd.DataFrame({"gene": names_c.values, "log2fc": log2fc_c.values, "fdr": fdr_c.values})
            else:
                dfc, _ = _moderated_t_test(pooled_pat, covariate_col, case_label, control_label, coreg_name)

            dfc["pattern"] = coreg_name
            coreg_tables.append(dfc)



            coreg_tables.append(dfc)
        except Exception as e:
            print("[WARN] Co-reg pooled DE failed for pattern {}: {}".format(pat, e))


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
    if "lineage_order" in adata.uns:
        de_store["lineage_order"] = list(adata.uns["lineage_order"])
        #print(f"[INFO] Preserved lineage_order in de_store (n={len(de_store['lineage_order'])}).")

    de_store["use_rawp"] = bool(use_rawp)

    # Ensure progress bar completes at the end
    if progress_callback is not None:
        progress_callback()

    return de_store


# -------------------- Step 5: Heatmap (fixed order), TSV ------------------- #

def build_fixed_order_heatmap(de_store, outdir, population_col, fc_thresh, heatmap_png, heatmap_tsv):
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

    if assigned.empty or detailed.empty:
        print("[WARN] No DEGs to plot heatmap.")
        return None, None

    # ------------------------- 1. Identify local-significant genes ------------------------- #
    local_genes = set(assigned.loc[assigned["group"] == "local", "gene"])
    if len(local_genes) == 0:
        print("[WARN] No locally significant genes; skipping heatmap.")
        return None, None

    # Include only local genes and any global/coreg genes that overlap with local
    include_genes = set(assigned.loc[assigned["gene"].isin(local_genes), "gene"])
    assigned = assigned[assigned["gene"].isin(include_genes)].copy()

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
    c_up_list   = [g for g in c_up["gene"].tolist()   if g in fold_df.index]
    c_down_list = [g for g in c_down["gene"].tolist() if g in fold_df.index]

    # (3) Local per-pop up/down — now ordered by lineage_order if available
    ordered_local_genes = []
    local_block_sizes = []

    # Use lineage-based order for populations if available
    lineage_order = de_store.get("lineage_order", pop_order)
    lineage_order = [p for p in lineage_order if p in pop_order]

    for pop in lineage_order:
        sub = local_hits[local_hits["population_or_pattern"] == str(pop)]
        u, d = split_up_down(sub)
        u_list = [g for g in u["gene"].tolist() if g in fold_df.index]
        d_list = [g for g in d["gene"].tolist() if g in fold_df.index]
        ordered_local_genes.extend(u_list + d_list)
        local_block_sizes.append(len(u_list))  # up genes for this pop
        local_block_sizes.append(len(d_list))  # down genes for this pop

    # Combine master order: global → coreg → local
    ordered_genes = g_up_list + g_down_list + c_up_list + c_down_list + ordered_local_genes

    if len(ordered_genes) == 0:
        print("[WARN] No ordered genes found; skipping heatmap.")
        return None, None

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

    # Reorder the matrix rows by the final order
    fold_df = fold_df.loc[ordered_genes]

    # ------------------------- 4. Apply lineage order if available ------------------------- #
    lineage_order = None
    # Prefer de_store (populated during DE run) over undefined adata
    if "lineage_order" in de_store:
        lineage_order = list(de_store["lineage_order"])
        print(f"[INFO] Using lineage order from de_store (n={len(lineage_order)}).")
    elif "adata" in de_store and hasattr(de_store["adata"], "uns") and "lineage_order" in de_store["adata"].uns:
        lineage_order = list(de_store["adata"].uns["lineage_order"])
        print(f"[INFO] Using lineage order from adata.uns (n={len(lineage_order)}).")

    if lineage_order is not None:
        # --- Reorder strictly by lineage_order but keep only populations with DE results ---
        existing_pops = [p for p in lineage_order if p in fold_df.columns]
        missing_pops = [p for p in lineage_order if p not in fold_df.columns]
        if missing_pops:
            print(f"[WARN] {len(missing_pops)} lineage populations missing from DE results: {missing_pops}")
        if existing_pops:
            fold_df = fold_df.reindex(columns=existing_pops)
            print(f"[INFO] Reordered heatmap columns by lineage_order (n={len(existing_pops)}).")
        else:
            print("[WARN] No matching lineage_order populations found in DE matrix; keeping default order.")
    else:
        print("[WARN] lineage_order not found in de_store or adata; keeping existing column order.")

    # ------------------------- 4. Export TSV ------------------------- #
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    heatmap_path = os.path.join(outdir, heatmap_tsv)

    # --- Ensure heatmap.tsv follows lineage_order exactly as the plotted heatmap ---
    if pop_order and all(p in fold_df.columns for p in pop_order):
        fold_df = fold_df.loc[:, [p for p in pop_order if p in fold_df.columns]]
    else:
        print("[WARN] lineage_order incomplete or missing columns; using current order.")

    fold_df.to_csv(heatmap_path, sep="\t")
    print("[INFO] Wrote heatmap matrix TSV: {}".format(heatmap_path))


    # ------------------------- 5. Draw heatmap ------------------------- #
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

    fixed_height = 7  # adjust between 6–12 depending on display
    fig = plt.figure(figsize=(max(6.0, 0.25 * len(pop_order)), fixed_height))

    ax = plt.gca()

    # --- FIX: remove white space caused by NaN folds in missing populations ---
    # Fill NaN with 0 (neutral fold) so heatmap renders continuously
    # Do NOT reorder by pop_order again; preserve lineage_order from above
    fold_df = fold_df.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    from matplotlib import colors

    # --- Create normalization with custom scaling ---
    norm = colors.Normalize(vmin=vmin, vmax=vmax)

    # --- Draw heatmap with explicit normalization ---
    im = ax.imshow(fold_df.values, aspect="auto", cmap=cmap,
                norm=norm, interpolation="none")

    # --- Draw the heatmap ---
    im = ax.imshow(fold_df.values, aspect="auto", cmap=cmap, norm=norm, interpolation="none")

    # --- Add vertical and horizontal white gridlines ---
    for x in range(1, fold_df.shape[1]):
        ax.axvline(x - 0.5, color="white", linewidth=0.2, alpha=0.9)

    for y in block_breaks:
        ax.axhline(y - 0.5, color="white", linewidth=0.2, alpha=0.9)

    ax.set_yticks([])
    ax.set_xticks(np.arange(fold_df.shape[1]))
    ax.set_xticklabels(list(fold_df.columns), rotation=90)
    ax.set_title("cellHarmony DE log2FC ({} vs {}) by {}".format(
        de_store["case_label"], de_store["control_label"], population_col))

    # --- Create smaller colorbar (1/8 height) that respects normalization ---
    cbar = plt.colorbar(im, ax=ax, fraction=0.08, pad=0.02)
    cbar.set_label(f"log2 fold change (yellow=up, blue=down)\nScale: {vmin:.1f} to {vmax:.1f}")
    cbar.set_ticks([vmin, 0, vmax])
    cbar.ax.set_yticklabels([f"{vmin:.1f}", "0", f"{vmax:.1f}"])

    pdf_path = os.path.join(outdir, heatmap_png)

    # Save high-quality raster and vector outputs
    plt.savefig(pdf_path, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.savefig(pdf_path.replace(".png", ".pdf"), dpi=600, bbox_inches="tight", pad_inches=0.02, transparent=True)
    plt.savefig(pdf_path.replace(".png", ".svg"), bbox_inches="tight", transparent=True)

    plt.close(fig)
    print(f"[INFO] Wrote heatmap images: {pdf_path}, {pdf_path.replace('.png', '.pdf')}, {pdf_path.replace('.png', '.svg')}")


    return pdf_path, heatmap_path

# ------------------------- Step 6: differentials h5ad ---------------------- #

def write_differentials_only_h5ad(adata, de_store, out_path):
    detailed = de_store["detailed_deg"]
    if detailed.empty:
        print("[WARN] No DEGs; writing a copy of input with DE container only.")
        adx = adata.copy()
        if "cellHarmony_DE" in adx.uns:
            del adx.uns["cellHarmony_DE"]
        adx.uns["cellHarmony_DE"] = de_store
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
        if "cellHarmony_DE" in adx.uns:
            del adx.uns["cellHarmony_DE"]
        adx.uns["cellHarmony_DE"] = de_store
        adx.write(out_path)
        print("[INFO] Wrote differentials h5ad: {}".format(out_path))
        return

    adx = adata[:, keep_mask].copy()
    if "cellHarmony_DE" in adx.uns:
        del adx.uns["cellHarmony_DE"]
    adx.uns["cellHarmony_DE"] = de_store
    adx.write(out_path)
    print("[INFO] Wrote differentials h5ad: {}".format(out_path))

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
    ap.add_argument("--outdir", default="cellHarmony_DE_out", help="Output directory (default: cellHarmony_DE_out)")
    ap.add_argument("--skip_grn", action="store_true", help="Skip interaction network (GRN) generation")
    args = ap.parse_args()

    comps = parse_comparisons_arg(args.comparisons)
    if len(comps) == 0:
        print("[ERROR] No valid comparisons parsed from --comparisons.", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    interactions_df = None
    interaction_root = os.path.join(args.outdir, "interaction-plots")
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
            covariate_col=args.covariate_col,   # <-- add this new argument
            min_cells=int(args.pseudobulk_min_cells),
            outdir=args.outdir
        )
        # Replace the working dataset with pseudobulk for all downstream DE
        adata = ad.read_h5ad(pb_h5ad)


    # Step 4–6 per comparison
    for case_label, control_label in comps:
        tag = "{}_vs_{}".format(str(case_label), str(control_label)).replace(" ", "_")
        print("[INFO] Differential analysis: {} vs {}.".format(case_label, control_label))

        # Identify total populations to process for progress tracking
        all_pops = sorted(adata.obs[args.population_col].unique())

        print(f"[INFO] Running differential expression across {len(all_pops)} populations...")

        # Initialize a persistent, dynamically updating progress bar
        pbar = tqdm(
            total=len(all_pops),
            desc="Running DE per population",
            ncols=100,
            dynamic_ncols=True,
            position=0,
            leave=True
        )

        def progress_callback():
            """Increment progress bar with live refresh."""
            pbar.update(1)
            pbar.refresh()


        # Pass the full callback to the DE function
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

        # Close the progress bar after completion
        pbar.close()
        sys.stdout.flush()

        # Determine population order for heatmap using lineage hierarchy if available
        if "lineage_order" in adata.uns:
            pop_order = [str(x) for x in list(adata.uns["lineage_order"])]
            print(f"[INFO] Using lineage order from h5ad (n={len(pop_order)} populations).")
        else:
            # Fallback to alphabetical order if not found
            detailed = de_store.get("detailed_deg")
            if detailed is not None and "population_or_pattern" in detailed.columns:
                pop_order = sorted(detailed["population_or_pattern"].astype(str).unique().tolist())
            else:
                pop_order = []
            print("[WARN] 'lineage_order' not found in h5ad; using alphabetical order.")

        # Step 5: heatmap (fixed order) and TSV
        heat_png = "heatmap_{}_by_{}.pdf".format(tag, args.population_col)
        heat_tsv = "heatmap_{}_by_{}.tsv".format(tag, args.population_col)
        build_fixed_order_heatmap(de_store, args.outdir, args.population_col, float(args.fc),
                                  heatmap_png=heat_png, heatmap_tsv=heat_tsv)

        # Step 6: differentials-only h5ad
        out_h5ad = os.path.join(args.outdir, "differentials_only_{}.h5ad".format(tag))
        write_differentials_only_h5ad(adata, de_store, out_h5ad)

        # Also export summary/detailed as TSV
        det = de_store["detailed_deg"]
        summ = de_store["summary_per_population"]
        assign = de_store["assigned_groups"]
        po = de_store["pooled_overall"]
        cr = de_store["coreg_pooled"]

        det_path = os.path.join(args.outdir, "DEG_detailed_{}.tsv".format(tag))
        summ_path = os.path.join(args.outdir, "DEG_summary_{}.tsv".format(tag))
        assign_path = os.path.join(args.outdir, "DEG_assigned_groups_{}.tsv".format(tag))
        po_path = os.path.join(args.outdir, "DEG_pooled_overall_{}.tsv".format(tag))
        cr_path = os.path.join(args.outdir, "DEG_coreg_pooled_{}.tsv".format(tag))

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
            per_pop_results = de_store.get("per_population_deg", {})
            if not per_pop_results:
                if diagnostic_report:
                    print("[DEBUG] No per-population DE results supplied for GRN export.")
            else:
                comparison_dir = os.path.join(interaction_root, NetPerspective.safe_component(tag))
                os.makedirs(comparison_dir, exist_ok=True)
                fc_threshold = max(0.0, float(np.log2(float(args.fc)))) if float(args.fc) > 0 else 0.0
                use_rawp = bool(de_store.get("use_rawp", False))

                for pop_name, pop_df in per_pop_results.items():
                    if pop_df is None or pop_df.empty:
                        continue

                    stats_df = pop_df.copy()
                    stats_df = stats_df.reset_index().rename(columns={stats_df.index.name or "index": "gene"})
                    stats_df["gene"] = stats_df["gene"].astype(str)
                    stats_df["log2fc"] = pd.to_numeric(stats_df.get("log2fc"), errors="coerce")

                    significance_column = None
                    if use_rawp and "pval" in stats_df.columns and stats_df["pval"].notna().any():
                        significance_column = "pval"
                    elif "fdr" in stats_df.columns:
                        significance_column = "fdr"

                    if significance_column is None:
                        if diagnostic_report:
                            print(f"[DEBUG] Skipping {pop_name}: missing significance column for GRN export.")
                        continue

                    stats_df[significance_column] = pd.to_numeric(
                        stats_df[significance_column], errors="coerce"
                    )

                    significance_mask = stats_df[significance_column] < float(args.alpha)
                    if fc_threshold > 0:
                        significance_mask &= stats_df["log2fc"].abs() >= fc_threshold

                    if not significance_mask.any():
                        continue

                    keep_columns = {"gene", "log2fc", "fdr", "pval", significance_column}
                    keep_columns = [col for col in keep_columns if col in stats_df.columns]
                    selected = (
                        stats_df.loc[significance_mask, keep_columns]
                        .dropna(subset=["gene"])
                        .drop_duplicates(subset=["gene"])
                    )

                    if selected.empty or selected["gene"].nunique() < 2:
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
                            max_genes=75,
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
                        args.skip_grn = True
                        break

    print("[INFO] Completed.")


if __name__ == "__main__":
    main()
