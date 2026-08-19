#!/usr/bin/env python3
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import sys
import time
import warnings

import matplotlib
import numpy as np
import pandas as pd
from anndata import AnnData
from pandas.errors import PerformanceWarning
import scanpy as sc
from scipy import sparse
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
from altanalyze3.components.visualization import NetPerspective

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['pdf.compression'] = 0
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'

warnings.filterwarnings(
    'ignore',
    category=PerformanceWarning,
    message='DataFrame is highly fragmented.*',
)

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


def _progress(iterable, **kwargs):
    if tqdm is None:
        return iterable
    return tqdm(iterable, **kwargs)

class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()


def _log_step_timing(label, started_at, timings=None, key=None):
    elapsed = max(0.0, float(time.perf_counter() - started_at))
    print(f"[timing] {label}={elapsed:.2f}s")
    if timings is not None:
        timings[key or label] = round(elapsed, 2)
    return elapsed


def _write_heatmap_cache(cache_path, heatmap_df, row_clusters, column_clusters, ordered_cells):
    matrix = np.asarray(heatmap_df.to_numpy(), dtype=np.float32)
    row_ids = np.asarray([f"{c}:{g}" for c, g in zip(row_clusters, heatmap_df.index.tolist())], dtype=str)
    col_ids = np.asarray([f"{c}:{b}" for c, b in zip(column_clusters, ordered_cells)], dtype=str)
    col_barcodes = np.asarray([str(barcode) for barcode in ordered_cells], dtype=str)
    np.savez_compressed(
        cache_path,
        matrix=matrix,
        row_ids=row_ids,
        col_ids=col_ids,
        col_barcodes=col_barcodes,
    )



def YellowBlackSky():
    cdict = {
        "red": [(0.0, 0.0, 0.0), (0.5, 0.0, 0.1), (1.0, 1.0, 1.0)],
        "green": [(0.0, 0.0, 0.8), (0.5, 0.1, 0.0), (1.0, 1.0, 1.0)],
        "blue": [(0.0, 0.0, 1.0), (0.5, 0.1, 0.0), (1.0, 0.0, 0.0)],
    }
    return LinearSegmentedColormap("YellowBlackSky", cdict)


def _coerce_lineage_order(raw_order):
    if raw_order is None:
        return None
    try:
        if isinstance(raw_order, dict):
            raw_order = list(raw_order.values())
        if hasattr(raw_order, "tolist"):
            raw_order = raw_order.tolist()
        if isinstance(raw_order, str):
            return None
        return [str(x) for x in raw_order]
    except Exception:
        return None


def _resolve_cluster_order(cluster_series, lineage_order):
    cluster_series = cluster_series.dropna()
    unique_clusters = [str(c) for c in pd.unique(cluster_series)]
    unique_set = set(unique_clusters)

    if lineage_order:
        ordered = [c for c in lineage_order if c in unique_set]
        missing = [c for c in unique_clusters if c not in ordered]
        if ordered:
            return ordered + missing

    if pd.api.types.is_categorical_dtype(cluster_series):
        cats = [str(c) for c in cluster_series.cat.categories if str(c) in unique_set]
        if cats:
            return cats

    return sorted(unique_clusters)


def downsample_cells_per_group(adata, groupby, cells_per_cluster=50, seed=0, group_order=None):
    if cells_per_cluster is None or cells_per_cluster <= 0:
        return adata.copy()

    rng = np.random.default_rng(seed)
    cluster_series = adata.obs[groupby].astype(str)
    if group_order is None:
        group_order = cluster_series.value_counts().index.tolist()

    idx = []
    for group in _progress(group_order, desc="Downsampling clusters"):
        cells = adata.obs_names[cluster_series == str(group)]
        if len(cells) == 0:
            continue
        selected = rng.choice(cells, min(len(cells), cells_per_cluster), replace=False)
        idx.extend(selected.tolist())

    return adata[idx].copy()


def _select_unique_markers(pvals_df, cluster_order, top_n, effect_df=None, pval_threshold=0.05):
    if pvals_df.empty:
        return pd.DataFrame(columns=["gene", "cluster", "pval"])

    pvals_df = pvals_df.reindex(columns=cluster_order)
    pvals_df = pvals_df.replace([np.inf, -np.inf], np.nan)
    mask = ~pvals_df.index.to_series().str.contains('rik', case=False, na=False)
    pvals_df = pvals_df.loc[mask]

    if effect_df is not None:
        effect_df = effect_df.reindex(index=pvals_df.index, columns=cluster_order)
        effect_df = effect_df.replace([np.inf, -np.inf], np.nan)
        pvals_df = pvals_df.where(effect_df > 0)

    if pval_threshold is not None:
        pvals_df = pvals_df.where(pvals_df <= pval_threshold)

    pvals_df = pvals_df.loc[pvals_df.notna().any(axis=1)]
    if pvals_df.empty:
        return pd.DataFrame(columns=["gene", "cluster", "pval"])

    pvals_df = pvals_df.fillna(1.0)
    min_cluster = pvals_df.idxmin(axis=1)
    min_pval = pvals_df.min(axis=1)

    markers = pd.DataFrame(
        {
            "gene": pvals_df.index.astype(str),
            "cluster": min_cluster.values,
            "pval": min_pval.values,
        }
    )
    markers = markers.dropna(subset=["cluster"])
    markers = markers[markers["cluster"].isin(cluster_order)]
    markers["cluster"] = pd.Categorical(markers["cluster"], categories=cluster_order, ordered=True)

    if effect_df is not None:
        effect_df = effect_df.reindex(index=markers["gene"], columns=cluster_order)
        markers["effect"] = [
            effect_df.loc[g, c] if g in effect_df.index and c in effect_df.columns else np.nan
            for g, c in zip(markers["gene"], markers["cluster"])
        ]
        markers["effect"] = markers["effect"].fillna(-np.inf)
        markers = markers.sort_values(
            ["cluster", "pval", "effect", "gene"],
            ascending=[True, True, False, True],
            kind="mergesort",
        )
    else:
        markers = markers.sort_values(
            ["cluster", "pval", "gene"],
            ascending=[True, True, True],
            kind="mergesort",
        )

    selected = markers.groupby("cluster", sort=False).head(top_n)
    return selected.reset_index(drop=True)


def _select_top_markers_per_cluster(pvals_df, cluster_order, top_n, effect_df=None):
    if pvals_df is None or pvals_df.empty:
        return pd.DataFrame(columns=["gene", "cluster", "pval", "effect"])

    pvals_df = pvals_df.reindex(columns=cluster_order)
    rows = []
    for cluster in cluster_order:
        if cluster not in pvals_df.columns:
            continue
        cluster_pvals = pd.to_numeric(pvals_df[cluster], errors="coerce").fillna(1.0)
        cluster_df = pd.DataFrame(
            {
                "gene": pvals_df.index.astype(str),
                "cluster": str(cluster),
                "pval": cluster_pvals.values,
            }
        )
        if effect_df is not None and cluster in effect_df.columns:
            cluster_effect = pd.to_numeric(effect_df[cluster], errors="coerce")
            cluster_df["effect"] = cluster_effect.reindex(pvals_df.index).values
        else:
            cluster_df["effect"] = np.nan
        cluster_df = cluster_df.sort_values(
            ["pval", "effect", "gene"],
            ascending=[True, False, True],
            kind="mergesort",
        )
        rows.append(cluster_df.head(int(top_n)))
    if not rows:
        return pd.DataFrame(columns=["gene", "cluster", "pval", "effect"])
    return pd.concat(rows, ignore_index=True)



def _zscore_rows(df):
    means = df.mean(axis=1)
    stds = df.std(axis=1, ddof=0).replace(0, np.nan)
    scaled = df.sub(means, axis=0).div(stds, axis=0)
    return scaled.fillna(0.0)


def _zscore_columns(df):
    means = df.mean(axis=0)
    stds = df.std(axis=0, ddof=0).replace(0, np.nan)
    scaled = df.sub(means, axis=1).div(stds, axis=1)
    return scaled.fillna(0.0)


def _bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranks = np.arange(1, n + 1)
    q = p[order] * n / ranks
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[order] = q
    return out


def _get_expression_matrix(adata, use_raw, layer):
    if use_raw:
        if adata.raw is None:
            raise ValueError("Requested use_raw=True but adata.raw is not set.")
        return adata.raw.X, adata.raw.var_names
    if layer:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers.")
        return adata.layers[layer], adata.var_names
    return adata.X, adata.var_names


def _prepare_marker_stats_aggregates(adata, cluster_key, genes, use_raw, layer):
    ordered_genes = pd.Index(pd.Series(genes, dtype=str)).drop_duplicates().tolist()
    if not ordered_genes:
        return None

    matrix, var_names = _get_expression_matrix(adata, use_raw, layer)
    var_names = pd.Index(var_names).astype(str)
    gene_indexer = var_names.get_indexer(ordered_genes)
    valid_gene_mask = gene_indexer >= 0
    if not np.any(valid_gene_mask):
        return None

    selected_genes = pd.Index(ordered_genes)[valid_gene_mask].astype(str).tolist()
    matrix = matrix[:, gene_indexer[valid_gene_mask]]

    cluster_obs = adata.obs[cluster_key]
    cluster_series = cluster_obs.astype(str)
    unique_clusters = set(pd.unique(cluster_series))
    if pd.api.types.is_categorical_dtype(cluster_obs):
        cluster_order = [str(c) for c in cluster_obs.cat.categories if str(c) in unique_clusters]
    else:
        cluster_order = [str(c) for c in pd.unique(cluster_series)]
    if not cluster_order:
        return None

    cluster_codes = pd.Categorical(cluster_series, categories=cluster_order, ordered=True).codes
    valid_cells = cluster_codes >= 0
    if not np.any(valid_cells):
        return None

    row_idx = np.nonzero(valid_cells)[0]
    col_idx = cluster_codes[valid_cells]
    indicator = sparse.csr_matrix(
        (np.ones(len(row_idx), dtype=np.float32), (row_idx, col_idx)),
        shape=(len(cluster_series), len(cluster_order)),
    )

    sums_matrix = indicator.T.dot(matrix)
    if sparse.issparse(sums_matrix):
        sums_values = sums_matrix.toarray()
    else:
        sums_values = np.asarray(sums_matrix)
    counts_values = np.bincount(col_idx, minlength=len(cluster_order)).astype(float)
    total_sum_values = np.asarray(matrix[valid_cells, :].sum(axis=0)).ravel().astype(float)
    total_count = float(np.sum(valid_cells))

    return {
        "sums": pd.DataFrame(sums_values, index=cluster_order, columns=selected_genes),
        "counts": pd.Series(counts_values, index=cluster_order),
        "total_sum": pd.Series(total_sum_values, index=selected_genes),
        "total_count": total_count,
    }


def _compute_marker_stats(adata, cluster_key, markers_df, fdr_df, use_raw, layer, aggregates=None):
    if markers_df.empty:
        return pd.DataFrame(columns=["Gene", "Fold", "Query Exp", "Ref Exp", "FDR p-value", "cluster"])

    if aggregates is None:
        aggregates = _prepare_marker_stats_aggregates(
            adata,
            cluster_key,
            markers_df["gene"].astype(str).tolist(),
            use_raw,
            layer,
        )
    if not aggregates:
        return pd.DataFrame(columns=["Gene", "Fold", "Query Exp", "Ref Exp", "FDR p-value", "cluster"])

    sums = aggregates["sums"]
    counts = aggregates["counts"]
    total_sum = aggregates["total_sum"]
    total_count = float(aggregates["total_count"])

    sums_values = sums.to_numpy(dtype=float)
    row_clusters = pd.Categorical(markers_df["cluster"].astype(str), categories=sums.index).codes
    col_genes = pd.Categorical(markers_df["gene"].astype(str), categories=sums.columns).codes
    valid = (row_clusters >= 0) & (col_genes >= 0)
    if not np.any(valid):
        return pd.DataFrame(columns=["Gene", "Fold", "Query Exp", "Ref Exp", "FDR p-value", "cluster"])

    query_sum = np.full(len(markers_df), np.nan, dtype=float)
    query_sum[valid] = sums_values[row_clusters[valid], col_genes[valid]]

    marker_clusters = markers_df["cluster"].astype(str)
    marker_genes = markers_df["gene"].astype(str)
    query_count = counts.reindex(marker_clusters).to_numpy(dtype=float)
    total_gene_sum = total_sum.reindex(marker_genes).to_numpy(dtype=float)

    query_count = np.where(np.isfinite(query_count), query_count, 0.0)
    total_gene_sum = np.where(np.isfinite(total_gene_sum), total_gene_sum, np.nan)

    ref_sum = total_gene_sum - query_sum
    ref_count = np.maximum(total_count - query_count, 1.0)
    query_mean = query_sum / np.maximum(query_count, 1.0)
    ref_mean = ref_sum / ref_count

    eps = 1e-9
    fold = np.log2((query_mean + eps) / (ref_mean + eps))

    fdr = np.full(len(markers_df), np.nan, dtype=float)
    if fdr_df is not None and not fdr_df.empty:
        fdr_dedup = fdr_df.loc[~fdr_df.index.duplicated(keep="first")]
        fdr_dedup = fdr_dedup.loc[:, ~fdr_dedup.columns.duplicated(keep="first")]
        fdr_row = pd.Categorical(marker_genes, categories=fdr_dedup.index).codes
        fdr_col = pd.Categorical(marker_clusters, categories=fdr_dedup.columns).codes
        fdr_valid = (fdr_row >= 0) & (fdr_col >= 0)
        if np.any(fdr_valid):
            fdr_values = fdr_dedup.to_numpy(dtype=float)
            fdr[fdr_valid] = fdr_values[fdr_row[fdr_valid], fdr_col[fdr_valid]]

    out = pd.DataFrame(
        {
            "Gene": marker_genes,
            "Fold": fold,
            "Query Exp": query_mean,
            "Ref Exp": ref_mean,
            "FDR p-value": fdr,
            "cluster": marker_clusters,
        }
    )
    out = out.loc[np.isfinite(out["Query Exp"].to_numpy(dtype=float))].reset_index(drop=True)
    return out


def _order_cells_by_cluster(adata, cluster_key, cluster_order):
    ordered_cells = []
    cluster_counts = []
    cluster_series = adata.obs[cluster_key].astype(str)

    for cluster in cluster_order:
        cells = adata.obs_names[cluster_series == str(cluster)]
        if len(cells) == 0:
            continue
        ordered_cells.extend(cells.tolist())
        cluster_counts.append((str(cluster), len(cells)))

    return ordered_cells, cluster_counts


def _build_heatmap(adata, cluster_key, genes, cluster_order, use_raw, layer, include_column_zscore=True):
    if not genes:
        return pd.DataFrame(), pd.DataFrame(), [], []

    matrix, gene_reference = _get_expression_matrix(adata, use_raw, layer)
    gene_reference = pd.Index(gene_reference).astype(str)
    genes_present = [str(g) for g in genes if str(g) in gene_reference]
    if not genes_present:
        return pd.DataFrame(), pd.DataFrame(), [], []

    ordered_cells, cluster_counts = _order_cells_by_cluster(adata, cluster_key, cluster_order)
    if not ordered_cells:
        return pd.DataFrame(), pd.DataFrame(), [], []

    gene_reference = pd.Index(gene_reference).astype(str)
    cell_indexer = adata.obs_names.get_indexer(ordered_cells)
    gene_indexer = gene_reference.get_indexer(genes_present)
    if np.any(cell_indexer < 0) or np.any(gene_indexer < 0):
        raise ValueError("Failed to align selected cells or genes for heatmap construction.")

    submatrix = matrix[cell_indexer, :][:, gene_indexer]
    if sparse.issparse(submatrix):
        values = submatrix.toarray()
    else:
        values = np.asarray(submatrix)

    heatmap_raw_df = pd.DataFrame(values.T, index=genes_present, columns=ordered_cells)
    heatmap_col_df = _zscore_columns(heatmap_raw_df) if include_column_zscore else pd.DataFrame()
    heatmap_row_df = _zscore_rows(heatmap_raw_df)
    return heatmap_row_df, heatmap_col_df, cluster_counts, ordered_cells


def _build_marker_centroids(adata, cluster_key, genes, cluster_order, use_raw, layer):
    if not genes:
        return pd.DataFrame()

    expr_df = sc.get.obs_df(adata, keys=genes, use_raw=use_raw, layer=layer)
    cluster_series = adata.obs[cluster_key].astype(str)
    expr_df[cluster_key] = cluster_series.values

    ordered_clusters = [cluster for cluster in cluster_order if cluster in set(cluster_series)]
    if not ordered_clusters:
        ordered_clusters = [str(cluster) for cluster in pd.unique(cluster_series)]

    centroid_columns = {}
    for cluster in ordered_clusters:
        subset = expr_df.loc[expr_df[cluster_key] == str(cluster), genes]
        if subset.empty:
            continue
        summed = subset.sum(axis=0).astype(float)
        total_counts = float(summed.sum())
        if total_counts <= 0:
            normalized = pd.Series(0.0, index=genes)
        else:
            normalized = np.log2((summed / total_counts) * 10000.0 + 1.0)
        centroid_columns[str(cluster)] = normalized.reindex(genes)

    return pd.DataFrame(centroid_columns, index=genes)


def _build_marker_centroids_from_aggregates(genes, cluster_order, aggregates):
    if not genes or not aggregates:
        return pd.DataFrame()
    sums = aggregates.get("sums")
    if sums is None or sums.empty:
        return pd.DataFrame()

    genes = [str(g) for g in genes if str(g) in sums.columns]
    if not genes:
        return pd.DataFrame()

    ordered_clusters = [str(cluster) for cluster in cluster_order if str(cluster) in sums.index]
    if not ordered_clusters:
        ordered_clusters = [str(cluster) for cluster in sums.index]
    if not ordered_clusters:
        return pd.DataFrame()

    cluster_gene_sums = sums.reindex(index=ordered_clusters, columns=genes, fill_value=0.0).to_numpy(dtype=float)
    # The centroid is a pseudobulk CPM: log2(cluster_sum / cluster_total * 10000 + 1). The
    # formula is defined only for non-negative counts. Summing mean-centred or z-scored values
    # gives negative cluster totals, and log2 of the resulting negative ratio is NaN, which
    # would write a silently unusable cellHarmony reference.
    if np.any(cluster_gene_sums < 0):
        n_neg = int((cluster_gene_sums < 0).sum())
        raise ValueError(
            "Centroid matrix cannot be built from mean-centred or z-scored input: "
            "%d of %d cluster x gene sums are negative. The centroid step needs raw counts "
            "(it applies its own log2 CPM). Point --layer at a counts matrix."
            % (n_neg, cluster_gene_sums.size)
        )
    totals = cluster_gene_sums.sum(axis=1, keepdims=True)
    normalized = np.zeros_like(cluster_gene_sums, dtype=float)
    valid = totals[:, 0] > 0
    if np.any(valid):
        normalized[valid] = np.log2((cluster_gene_sums[valid] / totals[valid]) * 10000.0 + 1.0)
    return pd.DataFrame(normalized.T, index=genes, columns=ordered_clusters)


def _append_suffix_before_extension(path, suffix):
    root, ext = os.path.splitext(path)
    if ext:
        return f"{root}{suffix}{ext}"
    return f"{path}{suffix}"


def _strip_marker_cluster_prefix(value):
    text = str(value or "")
    if ":" not in text:
        return text
    _, suffix = text.split(":", 1)
    return suffix or text


def _cluster_color_map(cluster_order):
    palette = plt.get_cmap("tab20")
    colors = [palette(i % palette.N) for i in range(len(cluster_order))]
    return {cluster: colors[i] for i, cluster in enumerate(cluster_order)}, ListedColormap(colors)


def _axis_label_fontsize(label_count, axis):
    label_count = max(1, int(label_count))
    if axis == "x":
        if label_count <= 6:
            return 10
        if label_count <= 12:
            return 9
        if label_count <= 20:
            return 8
        if label_count <= 32:
            return 7
        if label_count <= 48:
            return 6
        if label_count <= 72:
            return 5
        if label_count <= 96:
            return 4
        return 3

    if label_count <= 20:
        return 9
    if label_count <= 40:
        return 8
    if label_count <= 80:
        return 7
    if label_count <= 140:
        return 6
    if label_count <= 220:
        return 5
    if label_count <= 320:
        return 4
    return 3


def _format_heatmap_pvalue(pvalue):
    try:
        pval = float(pvalue)
    except Exception:
        return "p=NA"
    if pval < 1e-3:
        return "p={:.1e}".format(pval)
    return "p={:.3f}".format(pval)


def _category_color_map(categories, cmap_name="tab20"):
    cmap = plt.get_cmap(cmap_name, max(1, len(categories)))
    return {str(cat): tuple(cmap(i)[:3]) for i, cat in enumerate(categories)}


def _covariate_palette_name(index):
    palettes = [
        "tab20",
        "Set3",
        "Paired",
        "Dark2",
        "Accent",
        "Pastel1",
        "Pastel2",
        "Set2",
        "tab20b",
        "tab20c",
    ]
    return palettes[int(index) % len(palettes)]


def _covariate_categories(values):
    series = pd.Series(values).astype(str)
    if pd.api.types.is_categorical_dtype(values):
        present = set(series)
        return [str(c) for c in values.cat.categories if str(c) in present]
    return [str(c) for c in pd.unique(series)]


def _build_covariate_track(covariate_df):
    if covariate_df is None or covariate_df.empty:
        return None, {}
    rows = []
    legends = {}
    for idx, col in enumerate(covariate_df.columns):
        values = covariate_df[col].astype(str)
        categories = _covariate_categories(covariate_df[col])
        color_map = _category_color_map(categories, cmap_name=_covariate_palette_name(idx))
        rows.append([color_map.get(str(v), (0.8, 0.8, 0.8)) for v in values])
        legends[str(col)] = {
            "categories": categories,
            "colors": color_map,
        }
    return np.asarray(rows, dtype=float), legends


def _filter_covariate_columns_for_heatmap(obs, ordered_cells, covariate_columns, max_categories=12):
    retained = []
    skipped = []
    for covariate in covariate_columns:
        categories = _covariate_categories(obs.loc[ordered_cells, covariate])
        if len(categories) > int(max_categories):
            skipped.append((covariate, len(categories)))
        else:
            retained.append(covariate)
    return retained, skipped


def _draw_go_term_labels(fig, ax_terms, ax, row_clusters, go_terms, row_color_map, max_terms):
    if ax_terms is None or not go_terms:
        return
    ax_terms.set_xlim(0, 1)
    ax_terms.set_ylim(ax.get_ylim())
    ax_terms.set_xticks([])
    ax_terms.set_yticks([])
    ax_terms.axis("off")
    ax_terms.patch.set_alpha(0.0)

    cluster_ranges = []
    current = None
    start = 0
    for idx, cluster in enumerate(row_clusters):
        cluster = str(cluster)
        if cluster != current:
            if current is not None:
                cluster_ranges.append((current, start, idx))
            current = cluster
            start = idx
    if current is not None:
        cluster_ranges.append((current, start, len(row_clusters)))

    fig.canvas.draw()
    axis_height_px = ax.get_window_extent().height
    row_height_px = axis_height_px / max(1, len(row_clusters))
    term_font_size = 6
    term_spacing_px = max(10.0, term_font_size * fig.dpi / 72.0 * 3.0)
    gene_scale = max(1.0, len(row_clusters) / 600.0)
    min_rows_per_term = max(1.0, term_spacing_px / max(1.0, row_height_px)) * gene_scale

    entries = []
    for cluster, start, end in cluster_ranges:
        terms = go_terms.get(str(cluster), []) or []
        block_height = max(1, end - start)
        slot_count = max(1, int(block_height / max(1.0, min_rows_per_term)))
        limit = min(len(terms), slot_count)
        if max_terms:
            limit = min(limit, int(max_terms))
        entries.append(
            {
                "cluster": str(cluster),
                "start": start,
                "end": end,
                "block_height": block_height,
                "terms": terms[:limit],
                "positions": [],
            }
        )

    def compute_positions(entry):
        terms = entry["terms"]
        if not terms:
            entry["positions"] = []
            return
        block_height = max(1, entry["end"] - entry["start"])
        available_rows = max(1.0, block_height - 1)
        total_height = (len(terms) - 1) * min_rows_per_term if len(terms) > 1 else 0.0
        total_height = min(total_height, available_rows)
        top_offset = max(0.0, (available_rows - total_height) / 2.0)
        start_y = entry["start"] + 0.5 + top_offset
        positions = []
        for i in range(len(terms)):
            y = start_y + i * min_rows_per_term
            if y > entry["end"] - 0.5:
                break
            positions.append(y)
        entry["positions"] = positions

    for entry in entries:
        compute_positions(entry)

    # Resolve collisions between adjacent cluster blocks, but NEVER prune a block to
    # zero terms. The previous loop dropped the last term from whichever block was
    # taller, with `prev` winning every tie. When all blocks are the same height --
    # the normal case, since each cluster contributes the same top_n markers -- the
    # tie always fell on `prev`, so each block in turn was stripped to empty and only
    # the final cluster kept a label. On the 42-cluster Adams run that rendered 1 of
    # 42 GO-Elite labels. Keeping a floor of one term per block preserves the intended
    # one-label-per-cluster minimum; a lone label is centred in its own block by
    # compute_positions, so it stays inside that cluster's rows.
    def _n_terms(entry):
        return len(entry["terms"])

    for idx in range(1, len(entries)):
        prev = entries[idx - 1]
        curr = entries[idx]
        while prev["positions"] and curr["positions"]:
            if (curr["positions"][0] - prev["positions"][-1]) >= min_rows_per_term:
                break
            # Prune from the block that still has the most terms; stop once both are
            # down to their last label.
            if _n_terms(prev) <= 1 and _n_terms(curr) <= 1:
                break
            if _n_terms(prev) >= _n_terms(curr):
                prev["terms"] = prev["terms"][:-1]
                compute_positions(prev)
            else:
                curr["terms"] = curr["terms"][:-1]
                compute_positions(curr)

    for entry in entries:
        color = row_color_map.get(str(entry["cluster"]), "blue")
        for y, item in zip(entry["positions"], entry["terms"]):
            if isinstance(item, (tuple, list)) and len(item) >= 2:
                label = "{} {}".format(str(item[0]), _format_heatmap_pvalue(item[1]))
            else:
                label = str(item)
            ax_terms.text(
                0.995,
                y,
                label,
                transform=ax_terms.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=term_font_size,
                color=color,
                clip_on=False,
            )


def _draw_covariate_legends(fig, covariate_legends, rect=(0.255, 0.078, 0.690, 0.115)):
    if not covariate_legends:
        return
    import matplotlib.patches as mpatches
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    ax_leg = fig.add_axes(rect)
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(0, 1)
    ax_leg.axis("off")
    fig.canvas.draw()
    axis_width_pt = max(1.0, ax_leg.get_window_extent().width * 72.0 / fig.dpi)
    legend_font = FontProperties(family=plt.rcParams.get("font.family", "sans-serif"), size=8)

    def text_width_axes(label):
        if not str(label):
            return 0.0
        return TextPath((0, 0), str(label), prop=legend_font).get_extents().width / axis_width_pt

    max_x = 1.0
    title_x = 0.245
    item_start_x = 0.275
    sw = 0.022
    sh = 0.115
    row_h = 0.190
    for row_idx, (covariate, info) in enumerate(covariate_legends.items()):
        y_center = 0.82 - row_idx * row_h
        if y_center < 0.02:
            break
        ax_leg.text(
            title_x,
            y_center,
            str(covariate),
            fontsize=8,
            ha="right",
            va="center",
            weight="bold",
            transform=ax_leg.transAxes,
        )
        x = item_start_x
        categories = list(info["categories"])
        remaining_width = max(0.05, max_x - item_start_x)
        raw_widths = []
        for cat in categories:
            label = str(cat)
            raw_widths.append(sw + 0.010 + text_width_axes(label) + 0.030)
        width_scale = min(1.0, remaining_width / max(remaining_width, sum(raw_widths)))
        min_gap = 0.010
        for cat, raw_item_w in zip(categories, raw_widths):
            label = str(cat)
            item_w = raw_item_w * width_scale
            if x + sw + min_gap > max_x:
                break
            patch = mpatches.Rectangle(
                (x, y_center - sh / 2.0),
                sw,
                sh,
                transform=ax_leg.transAxes,
                color=info["colors"][label],
                linewidth=0.25,
                clip_on=False,
            )
            ax_leg.add_patch(patch)
            ax_leg.text(x + sw + min_gap, y_center, label, fontsize=8, ha="left", va="center", transform=ax_leg.transAxes)
            x += item_w


def _plot_heatmap(
    heatmap_df,
    output_path,
    cluster_counts,
    cluster_order,
    column_clusters,
    row_clusters,
    covariate_df=None,
    go_terms=None,
    go_terms_max=30,
):
    if heatmap_df.empty:
        raise ValueError("Heatmap dataframe is empty.")

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    show_go_terms = bool(go_terms)
    covariate_track, covariate_legends = _build_covariate_track(covariate_df)

    fixed_width = 10.8 if show_go_terms else 8.3
    fixed_height = 9.0 if covariate_track is not None else 8.5
    if show_go_terms:
        fig = plt.figure(figsize=(fixed_width, fixed_height))
        gs = fig.add_gridspec(
            1,
            3,
            width_ratios=[0.35, 0.50, 0.15],
            wspace=0.0,
            left=0.06,
            right=0.98,
            top=0.82,
            bottom=0.27 if covariate_track is not None else 0.10,
        )
        ax_terms = fig.add_subplot(gs[0, 0])
        ax = fig.add_subplot(gs[0, 1], sharey=ax_terms)
        ax_labels = fig.add_subplot(gs[0, 2], sharey=ax)
    else:
        fig, ax = plt.subplots(figsize=(fixed_width, fixed_height))
        fig.subplots_adjust(left=0.08, right=0.80, top=0.82, bottom=0.27 if covariate_track is not None else 0.08)
        ax_terms = None
        ax_labels = None

    contrast = float(os.environ.get("MARKER_HEATMAP_CONTRAST", "1.0"))
    if contrast <= 0:
        contrast = 1.0
    # Row-median-centred log2 values almost never reach +/-3, so a +/-3 scale wastes most of
    # the colour range and washes real differences out. +/-1.5 is the default; override with
    # MARKER_HEATMAP_RANGE (a single number, the half-range) or MARKER_HEATMAP_CONTRAST.
    half_range = float(os.environ.get("MARKER_HEATMAP_RANGE", "1.5"))
    if half_range <= 0:
        half_range = 1.5
    vmin, vmax = -half_range / contrast, half_range / contrast
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    im = ax.imshow(
        heatmap_df.values,
        aspect="auto",
        cmap=YellowBlackSky(),
        norm=norm,
        interpolation="none",
    )

    cluster_colors, cluster_cmap = _cluster_color_map(cluster_order)
    cluster_to_id = {cluster: idx for idx, cluster in enumerate(cluster_order)}

    centers = []
    labels = []
    boundaries = []
    pos = 0
    total = sum(count for _, count in cluster_counts)
    for cluster, count in cluster_counts:
        if count <= 0:
            continue
        start = pos
        end = pos + count - 1
        labels.append(cluster)
        centers.append((start + end) / 2.0)
        pos = end + 1
        if pos < total:
            boundaries.append(end + 0.5)

    x_fontsize = _axis_label_fontsize(len(labels), axis="x")
    if centers:
        ax.set_xticks(centers)
        ax.set_xticklabels(labels, rotation=45, ha="left", fontsize=x_fontsize)
    else:
        ax.set_xticks([])

    ax.xaxis.tick_top()
    ax.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True, length=0, pad=7.5)

    gene_fontsize = _axis_label_fontsize(heatmap_df.shape[0], axis="y")
    y_positions = np.arange(heatmap_df.shape[0])
    y_labels = heatmap_df.index
    if ax_labels is None:
        ax.set_yticks(y_positions)
        ax.set_yticklabels(y_labels, fontsize=gene_fontsize)
        ax.yaxis.tick_right()
        ax.tick_params(axis="y", left=False, right=True, labelleft=False, labelright=True, pad=1)
    else:
        ax.set_yticks([])
        ax.tick_params(axis="y", left=False, right=False, labelleft=False, labelright=False)

    for boundary in boundaries:
        ax.axvline(boundary, color="white", linewidth=0.5)

    divider = make_axes_locatable(ax)
    ax_top = divider.append_axes("top", size="1.125%", pad=0.02)
    ax_left = divider.append_axes("left", size="1.5%", pad=0.02)
    col_ids = np.array([cluster_to_id[c] for c in column_clusters], dtype=int)[None, :]
    row_ids = np.array([cluster_to_id[c] for c in row_clusters], dtype=int)[:, None]

    ax_top.imshow(col_ids, aspect="auto", cmap=cluster_cmap, interpolation="none")
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    ax_top.set_ylabel("")

    ax_left.imshow(row_ids, aspect="auto", cmap=cluster_cmap, interpolation="none")
    ax_left.set_xticks([])
    ax_left.set_yticks([])
    ax_left.set_xlabel("")

    if covariate_track is not None:
        cov_rows = int(covariate_track.shape[0])
        ax_cov = divider.append_axes("bottom", size=f"{max(4.0, cov_rows * 2.0)}%", pad=0.040)
        ax_cov.imshow(covariate_track, aspect="auto", interpolation="nearest")
        ax_cov.set_xticks([])
        ax_cov.set_yticks([])
        for idx, name in enumerate(covariate_df.columns):
            ax_cov.text(
                -0.006,
                idx,
                str(name),
                transform=ax_cov.get_yaxis_transform(),
                ha="right",
                va="center",
                fontsize=7,
                weight="bold",
                clip_on=False,
            )
        ax_cov.tick_params(axis="y", length=0)
        for spine in ax_cov.spines.values():
            spine.set_visible(False)
        ax_cov.set_facecolor("none")
        for boundary in boundaries:
            ax_cov.axvline(boundary, color="white", linewidth=0.35)
        for row_boundary in np.arange(0.5, cov_rows - 0.5 + 1e-9, 1.0):
            ax_cov.axhline(row_boundary, color="white", linewidth=0.35)
        _draw_covariate_legends(fig, covariate_legends)
        cbar_anchor_ax = ax_cov
    else:
        cbar_anchor_ax = ax

    if ax_terms is not None and ax_labels is not None:
        fig.canvas.draw()
        heatmap_bbox = ax.get_position()
        ax_terms.set_position([ax_terms.get_position().x0, heatmap_bbox.y0, ax_terms.get_position().width, heatmap_bbox.height])
        ax_labels.set_position([ax_labels.get_position().x0, heatmap_bbox.y0, ax_labels.get_position().width, heatmap_bbox.height])
        ax_labels.set_xlim(0, 1)
        ax_labels.set_ylim(ax.get_ylim())
        ax_labels.set_xticks([])
        ax_labels.set_yticks([])
        ax_labels.axis("off")
        ax_labels.patch.set_alpha(0.0)

        fig.canvas.draw()
        row_height_px = ax_labels.get_window_extent().height / max(1, heatmap_df.shape[0])
        row_height_pt = row_height_px * 72.0 / fig.dpi
        fitted_gene_fontsize = max(2.5, min(gene_fontsize, row_height_pt * 0.80))
        for pos, gene in enumerate(heatmap_df.index.astype(str)):
            ax_labels.text(
                0.02,
                pos,
                gene,
                transform=ax_labels.get_yaxis_transform(),
                ha="left",
                va="center",
                fontsize=fitted_gene_fontsize,
                fontstyle="italic",
                clip_on=False,
            )
        _draw_go_term_labels(fig, ax_terms, ax, row_clusters, go_terms, cluster_colors, go_terms_max)

    fig.canvas.draw()
    heatmap_bbox = ax.get_position()
    anchor_bbox = cbar_anchor_ax.get_position()
    cbar_width = heatmap_bbox.width * 0.35
    cbar_height = 0.012
    cbar_x = heatmap_bbox.x0 + (heatmap_bbox.width - cbar_width) / 2.0
    cbar_y = max(0.02, anchor_bbox.y0 - (0.026 if covariate_track is not None else 0.055))
    cax = fig.add_axes([cbar_x, cbar_y, cbar_width, cbar_height])
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.ax.set_xlim(vmin, vmax)
    cbar.set_label("Norm Exp (z-score)", fontsize=8, labelpad=2)
    cbar.set_ticks([])
    # :g not :.0f -- the default half-range is 1.5, which :.0f would print as "-2"/"2"
    cbar.ax.text(-0.08, 0.5, f"{vmin:g}", ha="right", va="center", transform=cbar.ax.transAxes)
    cbar.ax.text(1.08, 0.5, f"{vmax:g}", ha="left", va="center", transform=cbar.ax.transAxes)

    _save_figure(fig, output_path)
    plt.close(fig)


def _save_figure(fig, output_path):
    """Write the heatmap to `output_path` AND to a sibling `.svg`, always.

    Illustrator opens the SVG with the text still text, because `svg.fonttype='none'` keeps glyphs
    as characters instead of converting them to outlines, which matches what `pdf.fonttype=42`
    does for the PDF. Returns the paths written.

    An SVG failure never costs the caller the PDF: the primary write happens first, and a failed
    SVG is reported rather than raised.
    """
    written = [output_path]
    fig.savefig(output_path, facecolor="white", transparent=False)
    base, ext = os.path.splitext(output_path)
    svg_path = base + ".svg"
    if ext.lower() != ".svg":
        prior = plt.rcParams.get("svg.fonttype")
        try:
            plt.rcParams["svg.fonttype"] = "none"      # keep text editable
            fig.savefig(svg_path, facecolor="white", transparent=False)
            written.append(svg_path)
            print(f"Saved marker heatmap SVG to: {svg_path}")
        except Exception as exc:                        # never lose the primary figure over the SVG
            print(f"[WARN] SVG export failed for {svg_path}: {exc}")
        finally:
            plt.rcParams["svg.fonttype"] = prior
    return written


def _load_marker_finder():
    try:
        from cellHarmony import markerFinder as ch_marker_finder

        return ch_marker_finder
    except ImportError:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        components_dir = os.path.dirname(script_dir)
        if components_dir not in sys.path:
            sys.path.insert(0, components_dir)
        try:
            from cellHarmony import markerFinder as ch_marker_finder

            return ch_marker_finder
        except ImportError as exc:
            raise ImportError(
                "Unable to import cellHarmony.markerFinder. Ensure the altanalyze3 components directory is on PYTHONPATH."
            ) from exc


def _markerfinder_stats(adata, cluster_key, use_raw, layer, scale_data=False,
                        scale_factor=None, validate_scaling=True):
    """MarkerFinder correlations for a cells x genes matrix.

    ``validate_scaling`` defaults to True. MarkerFinder correlates each gene against a 0/1
    cluster indicator, which is sensitive to per-cell sequencing depth, so raw counts and
    log1p-without-normalization layers give markers ranked partly by cell depth.
    """
    marker_finder = _load_marker_finder()
    matrix, gene_names = _get_expression_matrix(adata, use_raw, layer)
    clusters = adata.obs[cluster_key].astype(str).tolist()
    r_df, p_df = marker_finder.marker_finder(
        matrix,
        clusters,
        gene_names=gene_names,
        scale_data=scale_data,
        scale_factor=scale_factor,
        validate_scaling=validate_scaling,
    )
    return p_df, r_df


def _generate_single_marker_network(cluster, selected, interactions_df, output_prefix):
    try:
        pdf_path, png_path, tsv_path = NetPerspective.generate_network_for_genes(
            selected,
            interactions_df,
            output_prefix,
            gene_column="gene",
            fold_change_column="log2fc",
            pval_column="fdr" if "fdr" in selected.columns else None,
            max_genes=None,
        )
    except NetPerspective.NetworkGenerationError as exc:
        print(f"[WARN] No marker interaction edges found for {cluster}; skipping network export. {exc}")
        return None
    except ImportError as exc:
        print(f"[WARN] Marker network export disabled: {exc}")
        return None
    return {
        "id": os.path.basename(output_prefix),
        "population": str(cluster),
        "pdf": pdf_path,
        "png": png_path,
        "tsv": tsv_path,
    }


def _export_marker_networks(marker_stats, out_dir, network_top_n=1000, network_jobs=1, species=None):
    if marker_stats is None or marker_stats.empty:
        return []
    try:
        interactions_df = NetPerspective.load_interaction_data(species=species)
    except Exception as exc:
        print(f"[WARN] Marker network export skipped: {exc}")
        return []

    out_dir = str(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    networks = []
    used_ids = set()
    working = marker_stats.copy()
    working["cluster"] = working["cluster"].astype(str)
    working["Gene"] = working["Gene"].astype(str)
    working["Fold"] = pd.to_numeric(working.get("Fold"), errors="coerce")
    working["FDR p-value"] = pd.to_numeric(working.get("FDR p-value"), errors="coerce")

    jobs = []
    for index, (cluster, cluster_df) in enumerate(working.groupby("cluster"), start=1):
        selected = (
            cluster_df.loc[:, [col for col in ("Gene", "Fold", "FDR p-value") if col in cluster_df.columns]]
            .dropna(subset=["Gene", "Fold"])
            .sort_values(["FDR p-value", "Fold", "Gene"], ascending=[True, False, True])
            .drop_duplicates(subset=["Gene"])
            .head(int(network_top_n))
            .rename(columns={"Gene": "gene", "Fold": "log2fc", "FDR p-value": "fdr"})
        )
        if selected.empty or selected["gene"].nunique() < 2:
            continue

        network_id = NetPerspective.safe_component(str(cluster), fallback=f"network_{index}")
        if network_id in used_ids:
            network_id = f"{network_id}_{index}"
        used_ids.add(network_id)
        output_prefix = os.path.join(out_dir, network_id)
        jobs.append((index, cluster, selected, output_prefix))

    if not jobs:
        return networks

    max_workers = max(1, int(network_jobs or 1))
    if max_workers <= 1:
        for _, cluster, selected, output_prefix in jobs:
            network = _generate_single_marker_network(cluster, selected, interactions_df, output_prefix)
            if network is not None:
                networks.append(network)
        return networks

    indexed_outputs = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_map = {
            executor.submit(_generate_single_marker_network, cluster, selected, interactions_df, output_prefix): idx
            for idx, cluster, selected, output_prefix in jobs
        }
        for future in as_completed(future_map):
            idx = future_map[future]
            try:
                network = future.result()
            except Exception as exc:
                print(f"[WARN] Marker network export task failed: {exc}")
                continue
            if network is not None:
                indexed_outputs.append((idx, network))
    indexed_outputs.sort(key=lambda x: x[0])
    networks.extend([entry for _, entry in indexed_outputs])
    return networks


def generate_marker_heatmap_from_adata(
    adata,
    *,
    cluster_key,
    out,
    top_n=5,
    markers_tsv=None,
    heatmap_tsv=None,
    heatmap_column_tsv=None,
    heatmap_cache=None,
    marker_method="scanpy",
    method="wilcoxon",
    use_raw=False,
    layer=None,
    scale_data=False,
    scale_factor=None,
    validate_scaling=True,
    cells_per_cluster=50,
    seed=0,
    export_networks=False,
    network_top_n=1000,
    network_jobs=1,
    species=None,
    write_heatmap_tsv=True,
    write_expression_tsv=True,
    write_heatmap_cache=True,
    pval_threshold=0.001,
    covariate_columns=None,
    go_terms=None,
    go_terms_max=30,
):
    total_started = time.perf_counter()
    timings = {}
    out = str(out)
    out_dir = os.path.dirname(out) or "."
    base_name = os.path.splitext(os.path.basename(out))[0]
    markers_tsv = markers_tsv or os.path.join(out_dir, f"{base_name}_markers.tsv")
    heatmap_tsv = (
        heatmap_tsv or os.path.join(out_dir, f"{base_name}_fold_matrix.tsv")
        if write_heatmap_tsv
        else None
    )
    heatmap_column_tsv = (
        heatmap_column_tsv or os.path.join(out_dir, f"{base_name}_exp_matrix.tsv")
        if write_expression_tsv
        else None
    )
    heatmap_cache = (
        heatmap_cache or os.path.join(out_dir, f"{base_name}_fold_matrix.npz")
        if write_heatmap_cache
        else None
    )
    centroids_base = heatmap_tsv or os.path.join(out_dir, f"{base_name}_fold_matrix.tsv")
    centroids_tsv = _append_suffix_before_extension(centroids_base, ".centroids")
    redundant_markers_tsv = os.path.join(out_dir, f"{base_name}_redundant_markers.tsv")

    if cluster_key not in adata.obs:
        available_obs = ", ".join(sorted(map(str, adata.obs.columns.tolist())))
        raise KeyError(
            f"Cluster key '{cluster_key}' not found in adata.obs. "
            f"Available obs columns: {available_obs}"
        )

    if use_raw and adata.raw is None:
        raise ValueError("Requested use_raw=True but adata.raw is not set.")

    if layer and layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found in adata.layers.")

    resolve_cluster_started = time.perf_counter()
    raw_clusters = adata.obs[cluster_key]
    lineage_order = _coerce_lineage_order(adata.uns.get("lineage_order", None))
    cluster_order = _resolve_cluster_order(raw_clusters, lineage_order)
    clusters = raw_clusters.astype(str)
    print(f"[INFO] Resolved {len(cluster_order)} clusters.")
    _log_step_timing(
        "marker_heatmap.resolve_cluster_order",
        resolve_cluster_started,
        timings=timings,
        key="resolve_cluster_order",
    )

    adata.obs[cluster_key] = pd.Categorical(clusters, categories=cluster_order, ordered=True)

    effect_df = None
    fdr_df = None
    marker_core_started = time.perf_counter()
    if marker_method == "markerfinder":
        print(f"[INFO] Running markerFinder using {'raw' if use_raw else (layer or 'X')} matrix.")
        pvals, effect_df = _markerfinder_stats(
            adata, cluster_key, use_raw, layer,
            scale_data=scale_data, scale_factor=scale_factor,
            validate_scaling=validate_scaling,
        )
        fdr_df = pvals.apply(_bh_fdr, axis=0)
    else:
        n_genes = adata.raw.n_vars if use_raw and adata.raw is not None else adata.n_vars
        print(f"[INFO] Running rank_genes_groups for {len(cluster_order)} clusters with {n_genes} genes.")
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method=method,
            use_raw=use_raw,
            layer=layer,
            n_genes=n_genes,
        )
        rg_df = sc.get.rank_genes_groups_df(adata, group=None)
        if rg_df.empty:
            raise ValueError("rank_genes_groups produced no results.")
        pval_col = "pvals_adj"
        if pval_col not in rg_df.columns or not rg_df[pval_col].notna().any():
            pval_col = "pvals"
        pvals = rg_df.pivot_table(index="names", columns="group", values=pval_col, aggfunc="first")
        if pval_col == "pvals_adj":
            fdr_df = pvals.copy()
        else:
            fdr_df = pvals.apply(_bh_fdr, axis=0)
        if "logfoldchanges" not in rg_df.columns:
            raise ValueError("logfoldchanges missing; cannot enforce upregulated marker selection.")
        effect_df = rg_df.pivot_table(index="names", columns="group", values="logfoldchanges", aggfunc="first")
    _log_step_timing(
        "marker_heatmap.compute_markers",
        marker_core_started,
        timings=timings,
        key="compute_markers",
    )

    unique_select_started = time.perf_counter()
    selected = _select_unique_markers(
        fdr_df,
        cluster_order,
        top_n,
        effect_df=effect_df,
        pval_threshold=pval_threshold,
    )
    print(f"[INFO] Selected {selected.shape[0]} markers after FDR/effect filtering.")
    _log_step_timing(
        "marker_heatmap.select_unique_markers",
        unique_select_started,
        timings=timings,
        key="select_unique_markers",
    )
    if selected.empty:
        raise ValueError("No markers were selected. Check inputs and parameters.")

    redundant_select_started = time.perf_counter()
    redundant_selected = _select_top_markers_per_cluster(
        fdr_df,
        cluster_order,
        250,
        effect_df=effect_df,
    )
    print(f"[INFO] Selected {redundant_selected.shape[0]} redundant markers for network export.")
    _log_step_timing(
        "marker_heatmap.select_redundant_markers",
        redundant_select_started,
        timings=timings,
        key="select_redundant_markers",
    )

    if cells_per_cluster and cells_per_cluster > 0:
        print(f"[INFO] Downsampling to {cells_per_cluster} cells per cluster (seed={seed}).")
    else:
        print("[INFO] Using all cells (no downsampling).")
    downsample_started = time.perf_counter()
    adata_plot = downsample_cells_per_group(
        adata,
        cluster_key,
        cells_per_cluster=cells_per_cluster,
        seed=seed,
        group_order=cluster_order,
    )
    _log_step_timing(
        "marker_heatmap.downsample_cells",
        downsample_started,
        timings=timings,
        key="downsample_cells",
    )
    print(f"[INFO] Heatmap will use {adata_plot.n_obs} cells.")
    adata_plot.obs[cluster_key] = pd.Categorical(
        adata_plot.obs[cluster_key].astype(str),
        categories=cluster_order,
        ordered=True,
    )

    selected_genes = selected["gene"].tolist()
    print(f"[INFO] Building heatmap matrix for {len(selected_genes)} markers.")
    build_heatmap_started = time.perf_counter()
    heatmap_df, heatmap_col_df, cluster_counts, ordered_cells = _build_heatmap(
        adata_plot,
        cluster_key,
        selected_genes,
        cluster_order,
        use_raw,
        layer,
        include_column_zscore=write_expression_tsv,
    )
    _log_step_timing(
        "marker_heatmap.build_heatmap_matrix",
        build_heatmap_started,
        timings=timings,
        key="build_heatmap_matrix",
    )
    if heatmap_df.empty:
        raise ValueError("Heatmap data is empty after filtering genes.")

    column_clusters = adata_plot.obs[cluster_key].astype(str).loc[ordered_cells].tolist()
    row_clusters = selected.set_index("gene").loc[heatmap_df.index, "cluster"].astype(str).tolist()
    covariate_df = None
    if covariate_columns:
        covariate_columns = [str(c).strip() for c in covariate_columns if str(c).strip()]
        covariate_columns = [c for c in dict.fromkeys(covariate_columns) if c in adata_plot.obs]
        if covariate_columns:
            retained_covariates, skipped_covariates = _filter_covariate_columns_for_heatmap(
                adata_plot.obs,
                ordered_cells,
                covariate_columns,
                max_categories=12,
            )
            if skipped_covariates:
                skipped_text = ", ".join(f"{name} ({count})" for name, count in skipped_covariates)
                print(f"[INFO] Skipping heatmap covariates with >12 displayed categories: {skipped_text}")
            if retained_covariates:
                covariate_df = adata_plot.obs.loc[ordered_cells, retained_covariates].copy()
                print(f"[INFO] Rendering heatmap covariate bars: {', '.join(retained_covariates)}")
            else:
                print("[INFO] No requested heatmap covariates passed the <=12 category display limit.")
        else:
            print("[INFO] No requested heatmap covariates were present in adata.obs.")
    print("[INFO] Rendering heatmap.")
    render_heatmap_started = time.perf_counter()
    _plot_heatmap(
        heatmap_df,
        out,
        cluster_counts,
        cluster_order,
        column_clusters,
        row_clusters,
        covariate_df=covariate_df,
        go_terms=go_terms,
        go_terms_max=go_terms_max,
    )
    _log_step_timing(
        "marker_heatmap.render_heatmap_pdf",
        render_heatmap_started,
        timings=timings,
        key="render_heatmap_pdf",
    )
    print(f"Saved marker heatmap to: {out}")

    print("[INFO] Computing marker statistics for TSV.")
    aggregate_stats_started = time.perf_counter()
    stat_genes = pd.Index(
        pd.concat(
            [
                selected["gene"].astype(str),
                redundant_selected["gene"].astype(str),
            ],
            ignore_index=True,
        )
    ).drop_duplicates().tolist()
    stats_aggregates = _prepare_marker_stats_aggregates(
        adata,
        cluster_key,
        stat_genes,
        use_raw,
        layer,
    )
    _log_step_timing(
        "marker_heatmap.prepare_stats_aggregates",
        aggregate_stats_started,
        timings=timings,
        key="prepare_stats_aggregates",
    )

    marker_stats_started = time.perf_counter()
    marker_stats = _compute_marker_stats(
        adata,
        cluster_key,
        selected,
        fdr_df,
        use_raw,
        layer,
        aggregates=stats_aggregates,
    )
    _log_step_timing(
        "marker_heatmap.compute_marker_stats_unique",
        marker_stats_started,
        timings=timings,
        key="compute_marker_stats_unique",
    )
    write_marker_stats_started = time.perf_counter()
    marker_stats.to_csv(markers_tsv, sep="\t", index=False, float_format="%.4g")
    _log_step_timing(
        "marker_heatmap.write_marker_stats_tsv",
        write_marker_stats_started,
        timings=timings,
        key="write_marker_stats_tsv",
    )
    print(f"Saved marker stats TSV to: {markers_tsv}")

    redundant_stats_started = time.perf_counter()
    redundant_marker_stats = _compute_marker_stats(
        adata,
        cluster_key,
        redundant_selected,
        fdr_df,
        use_raw,
        layer,
        aggregates=stats_aggregates,
    )
    _log_step_timing(
        "marker_heatmap.compute_marker_stats_redundant",
        redundant_stats_started,
        timings=timings,
        key="compute_marker_stats_redundant",
    )
    write_redundant_stats_started = time.perf_counter()
    redundant_marker_stats.to_csv(redundant_markers_tsv, sep="\t", index=False, float_format="%.4g")
    _log_step_timing(
        "marker_heatmap.write_redundant_marker_stats_tsv",
        write_redundant_stats_started,
        timings=timings,
        key="write_redundant_marker_stats_tsv",
    )
    print(f"Saved redundant marker stats TSV to: {redundant_markers_tsv}")

    write_heatmap_matrices_started = time.perf_counter()
    if write_heatmap_cache and heatmap_cache:
        write_cache_started = time.perf_counter()
        _write_heatmap_cache(
            heatmap_cache,
            heatmap_df,
            row_clusters,
            column_clusters,
            ordered_cells,
        )
        _log_step_timing(
            "marker_heatmap.write_heatmap_cache",
            write_cache_started,
            timings=timings,
            key="write_heatmap_cache",
        )
        print(f"Saved heatmap cache to: {heatmap_cache}")

    if write_heatmap_tsv and heatmap_tsv:
        heatmap_tsv_df = heatmap_df.copy()
        heatmap_tsv_df.index = [f"{c}:{g}" for c, g in zip(row_clusters, heatmap_tsv_df.index)]
        heatmap_tsv_df.columns = [f"{c}:{b}" for c, b in zip(column_clusters, ordered_cells)]
        heatmap_tsv_df.to_csv(heatmap_tsv, sep="\t", float_format="%.4g")
        print(f"Saved heatmap matrix TSV to: {heatmap_tsv}")
    else:
        print("[INFO] Skipping fold heatmap TSV export.")

    if write_expression_tsv and heatmap_column_tsv:
        heatmap_col_tsv_df = heatmap_col_df.copy()
        heatmap_col_tsv_df.index = [f"{c}:{g}" for c, g in zip(row_clusters, heatmap_col_tsv_df.index)]
        heatmap_col_tsv_df.columns = [f"{c}:{b}" for c, b in zip(column_clusters, ordered_cells)]
        heatmap_col_tsv_df.to_csv(heatmap_column_tsv, sep="\t", float_format="%.4g")
    else:
        print("[INFO] Skipping expression matrix TSV export.")
    _log_step_timing(
        "marker_heatmap.write_heatmap_tsvs",
        write_heatmap_matrices_started,
        timings=timings,
        key="write_heatmap_tsvs",
    )
    if write_expression_tsv and heatmap_column_tsv:
        print(f"Saved expression matrix TSV to: {heatmap_column_tsv}")

    centroid_started = time.perf_counter()
    centroid_df = _build_marker_centroids_from_aggregates(
        heatmap_df.index.tolist(),
        cluster_order,
        stats_aggregates,
    )
    centroid_df.index = [_strip_marker_cluster_prefix(gene) for gene in centroid_df.index]
    centroid_df.to_csv(centroids_tsv, sep="\t", float_format="%.4g")
    _log_step_timing(
        "marker_heatmap.build_and_write_centroids_tsv",
        centroid_started,
        timings=timings,
        key="build_and_write_centroids_tsv",
    )
    print(f"Saved centroid matrix TSV to: {centroids_tsv}")

    networks = []
    if export_networks:
        network_dir = os.path.join(out_dir, f"{base_name}_networks")
        print(f"[INFO] Exporting NetPerspective networks for up to {int(network_top_n)} markers per cluster.")
        networks_started = time.perf_counter()
        networks = _export_marker_networks(
            redundant_marker_stats,
            network_dir,
            network_top_n=network_top_n,
            network_jobs=network_jobs,
            species=species,
        )
        _log_step_timing(
            "marker_heatmap.export_networks",
            networks_started,
            timings=timings,
            key="export_networks",
        )

    _log_step_timing(
        "marker_heatmap.total",
        total_started,
        timings=timings,
        key="total",
    )

    return {
        "pdf": out,
        "markers_tsv": markers_tsv,
        "redundant_markers_tsv": redundant_markers_tsv,
        "heatmap_tsv": heatmap_tsv,
        "heatmap_column_tsv": heatmap_column_tsv,
        "heatmap_cache": heatmap_cache,
        "centroids_tsv": centroids_tsv,
        "networks": networks,
        "timings": timings,
    }


def _load_text_expression_inputs(matrix_tsv, groups_tsv, cluster_key):
    print(f"[INFO] Reading expression matrix TSV from: {matrix_tsv}")
    expr_df = pd.read_csv(matrix_tsv, sep="\t", index_col=0)
    if expr_df.empty:
        raise ValueError(f"Expression matrix '{matrix_tsv}' is empty.")

    expr_df.index = expr_df.index.astype(str)
    expr_df.columns = expr_df.columns.astype(str)
    expr_df = expr_df.loc[~expr_df.index.duplicated(keep="first")]
    expr_df = expr_df.loc[:, ~expr_df.columns.duplicated(keep="first")]
    expr_df = expr_df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    print(
        f"[INFO] Loaded expression matrix with {expr_df.shape[0]} genes and "
        f"{expr_df.shape[1]} cells."
    )

    print(f"[INFO] Reading groups file from: {groups_tsv}")
    groups_df = pd.read_csv(groups_tsv, sep="\t", header=None, dtype=str)
    if groups_df.empty:
        raise ValueError(f"Groups file '{groups_tsv}' is empty.")
    if groups_df.shape[1] < 2:
        raise ValueError(
            f"Groups file '{groups_tsv}' must contain at least 2 columns: cell barcode and group label."
        )

    label_col = 2 if groups_df.shape[1] >= 3 else 1
    groups_df = groups_df.iloc[:, [0, label_col]].copy()
    groups_df.columns = ["cell_id", cluster_key]
    groups_df["cell_id"] = groups_df["cell_id"].astype(str)
    groups_df[cluster_key] = groups_df[cluster_key].astype(str)
    groups_df = groups_df.loc[~groups_df["cell_id"].duplicated(keep="first")]

    group_cells = set(groups_df["cell_id"])
    common_cells = [cell for cell in expr_df.columns if cell in group_cells]
    if not common_cells:
        raise ValueError(
            "No overlapping cell identifiers were found between the expression matrix columns "
            "and the groups file."
        )

    missing_group_count = expr_df.shape[1] - len(common_cells)
    if missing_group_count:
        print(
            f"[INFO] Excluding {missing_group_count} matrix cells not present in the groups file."
        )
    extra_group_count = groups_df.shape[0] - len(common_cells)
    if extra_group_count:
        print(
            f"[INFO] Ignoring {extra_group_count} group rows not present in the expression matrix."
        )

    expr_df = expr_df.loc[:, common_cells]
    groups_df = groups_df.set_index("cell_id").loc[common_cells]

    obs = pd.DataFrame(index=common_cells)
    obs[cluster_key] = groups_df[cluster_key].values
    var = pd.DataFrame(index=expr_df.index.astype(str))
    adata = AnnData(X=expr_df.T.to_numpy(dtype=np.float32, copy=True), obs=obs, var=var)
    adata.obs_names = pd.Index(list(map(str, common_cells)))
    adata.var_names = pd.Index(expr_df.index.astype(str))
    print(
        f"[INFO] Built AnnData from text inputs with {adata.n_obs} cells and {adata.n_vars} genes."
    )
    return adata



def main():
    parser = argparse.ArgumentParser(
        description="Create a marker heatmap from an h5ad file or text expression inputs using unique markers per cluster."
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--h5ad", help="Input h5ad file.")
    input_group.add_argument(
        "--matrix-tsv",
        help="Log2-normalized expression matrix TSV with genes as rows and cells as columns.",
    )
    parser.add_argument(
        "--groups-tsv",
        default=None,
        help="Groups file for --matrix-tsv input. Column 1 is cell barcode, column 3 (or 2 if absent) is the cluster label.",
    )
    parser.add_argument(
        "--groups-txt",
        dest="groups_tsv",
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--cluster-key",
        default=None,
        help="obs column with cluster labels for h5ad input. Ignored for --matrix-tsv, which uses the groups file.",
    )
    parser.add_argument("--top-n", type=int, default=5, help="Markers per cluster (default 5).")
    parser.add_argument("--out", default="marker_heatmap.pdf", help="Output figure path.")
    parser.add_argument("--markers-tsv", default=None, help="Output TSV for marker statistics.")
    parser.add_argument("--heatmap-tsv", default=None, help="Output TSV for the row-scaled heatmap matrix.")
    parser.add_argument(
        "--heatmap-column-tsv",
        default=None,
        help="Output TSV for the column-scaled expression matrix.",
    )
    parser.add_argument(
        "--marker-method",
        choices=["scanpy", "markerfinder"],
        default="scanpy",
        help="Marker selection via scanpy rank_genes_groups or cellHarmony markerFinder.",
    )
    parser.add_argument("--method", default="wilcoxon", help="DE method for rank_genes_groups.")
    parser.add_argument("--use-raw", action="store_true", help="Use adata.raw for DE and heatmap.")
    parser.add_argument("--layer", default=None, help="Layer to use for DE and heatmap.")
    parser.add_argument(
        "--scale-data",
        action="store_true",
        help=(
            "Depth-normalize raw counts before MarkerFinder: "
            "log2(counts / total_counts_per_cell * scale_factor + 1). Required when .X or the "
            "chosen --layer holds raw counts."
        ),
    )
    parser.add_argument(
        "--scale-factor",
        type=float,
        default=None,
        help=(
            "Scale factor for --scale-data. Default: the median cell total, matching scanpy "
            "normalize_total(). Keep the default for targeted panels (Xenium, CosMx, MERFISH)."
        ),
    )
    parser.add_argument(
        "--no-scaling-check",
        action="store_true",
        help="Skip the MarkerFinder input-scaling check (not recommended).",
    )
    parser.add_argument(
        "--cells-per-cluster",
        type=int,
        default=50,
        help="Number of cells to display per cluster (default 50). Use 0 for all cells.",
    )
    parser.add_argument(
        "--export-networks",
        action="store_true",
        help="Export NetPerspective networks for each cluster using the top marker genes.",
    )
    parser.add_argument(
        "--network-top-n",
        type=int,
        default=1000,
        help="Top marker genes per cluster to use for network export (default 1000).",
    )
    parser.add_argument(
        "--network-jobs",
        type=int,
        default=1,
        help="Parallel workers for per-cluster network generation (default 1).",
    )
    parser.add_argument(
        "--species",
        default=None,
        help="Species for GRN interaction loading (human or mouse). Defaults to combined resources when omitted.",
    )
    parser.add_argument(
        "--skip-expression-tsv",
        action="store_true",
        help="Skip writing both expression-matrix and fold-heatmap TSV outputs.",
    )
    parser.add_argument(
        "--heatmap-cache",
        default=None,
        help="Output compressed NPZ cache for the marker heatmap matrix.",
    )
    parser.add_argument(
        "--skip-heatmap-cache",
        action="store_true",
        help="Skip writing the compressed heatmap cache NPZ.",
    )
    parser.add_argument(
        "--heatmap-covariates",
        default=None,
        help="Comma-delimited h5ad obs columns to render as compact bottom covariate bars.",
    )
    parser.add_argument("--seed", type=int, default=0, help="Random seed for cell sampling.")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out) or "."
    logs_dir = os.path.join(out_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(args.out))[0]
    log_path = os.path.join(logs_dir, f"{base_name}.log")
    markers_tsv = args.markers_tsv or os.path.join(out_dir, f"{base_name}_markers.tsv")
    heatmap_tsv = (
        args.heatmap_tsv or os.path.join(out_dir, f"{base_name}_fold_matrix.tsv")
        if not args.skip_expression_tsv
        else None
    )
    heatmap_column_tsv = (
        args.heatmap_column_tsv or os.path.join(out_dir, f"{base_name}_exp_matrix.tsv")
        if not args.skip_expression_tsv
        else None
    )
    heatmap_cache = (
        args.heatmap_cache or os.path.join(out_dir, f"{base_name}_fold_matrix.npz")
        if not args.skip_heatmap_cache
        else None
    )
    centroids_base = heatmap_tsv or os.path.join(out_dir, f"{base_name}_fold_matrix.tsv")
    centroids_tsv = _append_suffix_before_extension(centroids_base, ".centroids")
    redundant_markers_tsv = os.path.join(out_dir, f"{base_name}_redundant_markers.tsv")

    log_file = open(log_path, "w")
    tee_out = Tee(sys.stdout, log_file)
    tee_err = Tee(sys.stderr, log_file)
    original_stdout, original_stderr = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = tee_out, tee_err

    try:
        print("[INFO] Parameters:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print(f"[INFO] Log file: {log_path}")
        print(f"[INFO] Marker TSV: {markers_tsv}")
        print(f"[INFO] Redundant marker TSV: {redundant_markers_tsv}")
        if heatmap_tsv:
            print(f"[INFO] Heatmap TSV (fold matrix): {heatmap_tsv}")
        if heatmap_column_tsv:
            print(f"[INFO] Heatmap TSV (expression matrix): {heatmap_column_tsv}")
        if heatmap_cache:
            print(f"[INFO] Heatmap cache (npz): {heatmap_cache}")
        print(f"[INFO] Centroid TSV: {centroids_tsv}")

        layer = args.layer
        if args.use_raw and layer:
            print("Warning: --layer ignored when --use-raw is set.")
            layer = None

        cluster_key = args.cluster_key

        if args.h5ad:
            if not cluster_key:
                raise ValueError("--cluster-key is required when --h5ad is used.")
            print(f"[INFO] Reading AnnData from: {args.h5ad}")
            adata = sc.read_h5ad(args.h5ad)
            print(f"[INFO] Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes.")
        else:
            if not args.groups_tsv:
                raise ValueError("--groups-tsv is required when --matrix-tsv is used.")
            if args.use_raw:
                print("[WARN] --use-raw ignored for text matrix input.")
            if layer:
                print("[WARN] --layer ignored for text matrix input.")
            cluster_key = "cluster"
            adata = _load_text_expression_inputs(args.matrix_tsv, args.groups_tsv, cluster_key)
            layer = None

        if cluster_key not in adata.obs:
            available_obs = ", ".join(sorted(map(str, adata.obs.columns.tolist())))
            raise KeyError(
                f"Cluster key '{cluster_key}' not found in adata.obs. "
                f"Available obs columns: {available_obs}"
            )

        if args.use_raw and adata.raw is None:
            raise ValueError("Requested --use-raw but adata.raw is not set.")

        if layer and layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers.")

        outputs = generate_marker_heatmap_from_adata(
            adata,
            cluster_key=cluster_key,
            out=args.out,
            top_n=args.top_n,
            markers_tsv=markers_tsv,
            heatmap_tsv=heatmap_tsv,
            heatmap_column_tsv=heatmap_column_tsv,
            heatmap_cache=heatmap_cache,
            marker_method=args.marker_method,
            method=args.method,
            use_raw=args.use_raw,
            layer=layer,
            scale_data=args.scale_data,
            scale_factor=args.scale_factor,
            validate_scaling=not args.no_scaling_check,
            cells_per_cluster=args.cells_per_cluster,
            seed=args.seed,
            export_networks=args.export_networks,
            network_top_n=args.network_top_n,
            network_jobs=args.network_jobs,
            species=args.species,
            write_heatmap_tsv=not args.skip_expression_tsv,
            write_expression_tsv=not args.skip_expression_tsv,
            write_heatmap_cache=not args.skip_heatmap_cache,
            covariate_columns=[c.strip() for c in str(args.heatmap_covariates or "").split(",") if c.strip()],
        )
        timings = outputs.get("timings", {}) if isinstance(outputs, dict) else {}
        if timings:
            ordered_keys = [
                "resolve_cluster_order",
                "compute_markers",
                "select_unique_markers",
                "select_redundant_markers",
                "downsample_cells",
                "build_heatmap_matrix",
                "render_heatmap_pdf",
                "prepare_stats_aggregates",
                "compute_marker_stats_unique",
                "write_marker_stats_tsv",
                "compute_marker_stats_redundant",
                "write_redundant_marker_stats_tsv",
                "write_heatmap_cache",
                "write_heatmap_tsvs",
                "build_and_write_centroids_tsv",
                "export_networks",
                "total",
            ]
            summary = " | ".join(
                f"{k}={float(timings[k]):.2f}s" for k in ordered_keys if k in timings
            )
            if summary:
                print(f"[timing] marker_heatmap.summary {summary}")

        print(f"Saved log to: {log_path}")
    finally:
        sys.stdout, sys.stderr = original_stdout, original_stderr
        log_file.close()


if __name__ == "__main__":

    main()
