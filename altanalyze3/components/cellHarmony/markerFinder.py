"""
Marker gene discovery utilities for cellHarmony workflows.

This module extends the marker finder logic that exists in the UDON component by
adding a higher-level interface that accepts AnnData objects (or .h5ad files),
supports ordering clusters via ``adata.uns["lineage_order"]``, allows selecting
positive, negative, or bidirectional markers, and produces both tabular and
heatmap (figure + TSV) outputs.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Literal, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from matplotlib import cm
from matplotlib.colors import ListedColormap, to_hex
from scipy import sparse
from scipy.stats import t

# Default diverging colormap similar to the UDON visualizations
_N = 256
_vals = np.ones((_N * 2, 4))
_vals[:_N, 0] = np.linspace(15 / 256, 0, _N)
_vals[:_N, 1] = np.linspace(255 / 256, 0, _N)
_vals[:_N, 2] = np.linspace(255 / 256, 0, _N)
_vals[_N:, 0] = np.linspace(0, 255 / 256, _N)
_vals[_N:, 1] = np.linspace(0, 243 / 256, _N)
_vals[_N:, 2] = np.linspace(0, 15 / 256, _N)
BLUE_BLACK_YELLOW = ListedColormap(_vals)

MarkerDirection = Literal["up", "down", "both"]


@dataclass
class MarkerOutputs:
    markers: pd.DataFrame
    heatmap_values: pd.DataFrame
    cluster_assignments: pd.Series


def _ensure_output_dir(output_dir: Optional[str]) -> Optional[str]:
    if output_dir is None:
        return None
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def _coerce_sparse_matrix(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(matrix)


def _get_expression_matrix(
    adata: AnnData,
    layer: Optional[str],
    use_raw: bool,
) -> Tuple[sparse.csr_matrix, pd.Index, pd.Index]:
    if use_raw:
        if adata.raw is None:
            raise ValueError("Requested use_raw=True, but adata.raw is not set.")
        matrix = adata.raw.X
        var_names = adata.raw.var_names
    elif layer:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in AnnData.layers.")
        matrix = adata.layers[layer]
        var_names = adata.var_names
    else:
        matrix = adata.X
        var_names = adata.var_names
    return _coerce_sparse_matrix(matrix), adata.obs_names.copy(), pd.Index(var_names)


def _get_cluster_series(adata: AnnData, cluster_key: str) -> pd.Series:
    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster key '{cluster_key}' not present in adata.obs.")
    cluster_series = adata.obs[cluster_key].copy()
    cluster_series = cluster_series.dropna()
    return cluster_series.astype(str)


def _resolve_cluster_order(
    cluster_series: pd.Series,
    adata: Optional[AnnData],
    lineage_order_key: Optional[str],
    custom_order: Optional[Sequence[str]],
) -> List[str]:
    unique_clusters = pd.Index(cluster_series.unique())

    if custom_order:
        order = [c for c in custom_order if c in unique_clusters]
        if order:
            return order

    if adata is not None and lineage_order_key and lineage_order_key in adata.uns:
        raw_order = adata.uns[lineage_order_key]
        if isinstance(raw_order, dict):
            raw_order = list(raw_order.values())
        if isinstance(raw_order, Iterable) and not isinstance(raw_order, str):
            order = [c for c in raw_order if c in unique_clusters]
            if order:
                return order

    if pd.api.types.is_categorical_dtype(cluster_series):
        return [c for c in cluster_series.cat.categories if c in unique_clusters]

    return sorted(unique_clusters)


def _build_cluster_indicator(groups: Sequence[str]) -> Tuple[sparse.csr_matrix, pd.Index]:
    group_series = pd.Series(groups, dtype=str)
    categories = pd.Index(pd.unique(group_series))
    codes = pd.Categorical(group_series, categories=categories).codes
    valid_mask = codes >= 0
    row_idx = np.arange(len(group_series))[valid_mask]
    col_idx = codes[valid_mask]
    data = np.ones(len(row_idx), dtype=float)
    indicator = sparse.csr_matrix(
        (data, (row_idx, col_idx)),
        shape=(len(group_series), len(categories)),
    )
    return indicator, categories


def marker_finder(
    expression_matrix: sparse.spmatrix,
    groups: Sequence[str],
    gene_names: Optional[Sequence[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if not sparse.issparse(expression_matrix):
        expression_matrix = sparse.csr_matrix(expression_matrix)
    else:
        expression_matrix = expression_matrix.tocsr()

    n_cells, n_genes = expression_matrix.shape
    if len(groups) != n_cells:
        raise ValueError("Group labels must match the number of rows in the expression matrix.")

    indicator, cluster_names = _build_cluster_indicator(groups)

    sum_x = np.asarray(expression_matrix.sum(axis=0)).ravel()
    sum_x2 = np.asarray(expression_matrix.power(2).sum(axis=0)).ravel()
    sum_y = np.asarray(indicator.sum(axis=0)).ravel()

    sum_xy = expression_matrix.T.dot(indicator)
    if sparse.issparse(sum_xy):
        sum_xy = sum_xy.toarray()
    else:
        sum_xy = np.asarray(sum_xy)

    n = float(n_cells)
    numerator = sum_xy - np.outer(sum_x, sum_y) / n
    ssx = sum_x2 - (sum_x**2) / n
    ssy = sum_y - (sum_y**2) / n
    denom = np.sqrt(np.outer(ssx, ssy))

    with np.errstate(divide="ignore", invalid="ignore"):
        r = numerator / denom

    degrees_f = n_cells - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t_stat = r * np.sqrt(degrees_f) / np.sqrt(1 - (r**2))
        p = t.sf(np.abs(t_stat), df=degrees_f) * 2

    if gene_names is None:
        gene_names = pd.Index(range(n_genes)).astype(str)
    else:
        gene_names = pd.Index(gene_names)

    r_df = pd.DataFrame(r, index=gene_names, columns=cluster_names)
    p_df = pd.DataFrame(p, index=gene_names, columns=cluster_names)
    r_df = r_df.dropna(axis=0, how="any")
    p_df = p_df.loc[r_df.index, r_df.columns]
    return r_df, p_df


def _select_unique_markers(
    r_df: pd.DataFrame,
    p_df: pd.DataFrame,
    direction: MarkerDirection,
    rho_threshold: float,
    top_n: int,
    min_markers_per_cluster: int,
) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []

    if direction in ("up", "both"):
        top_cluster = r_df.idxmax(axis=1)
        top_r = r_df.max(axis=1)
        up_df = pd.DataFrame(
            {
                "marker": r_df.index,
                "cluster": top_cluster.values,
                "pearson_r": top_r.values,
                "direction": "up",
            }
        )
        up_df["p_value"] = [p_df.loc[g, c] for g, c in zip(r_df.index, top_cluster.values)]
        up_df["score"] = up_df["pearson_r"]
        up_df = up_df[up_df["pearson_r"] >= rho_threshold]
        frames.append(up_df)

    if direction in ("down", "both"):
        bottom_cluster = r_df.idxmin(axis=1)
        bottom_r = r_df.min(axis=1)
        down_df = pd.DataFrame(
            {
                "marker": r_df.index,
                "cluster": bottom_cluster.values,
                "pearson_r": bottom_r.values,
                "direction": "down",
            }
        )
        down_df["p_value"] = [p_df.loc[g, c] for g, c in zip(r_df.index, bottom_cluster.values)]
        down_df["score"] = down_df["pearson_r"].abs()
        down_df = down_df[down_df["pearson_r"] <= -rho_threshold]
        frames.append(down_df)

    if not frames:
        return pd.DataFrame(columns=["marker", "cluster", "pearson_r", "direction", "p_value", "score", "rank"])

    markers_df = pd.concat(frames, ignore_index=True)
    markers_df = markers_df.sort_values(["cluster", "direction", "score"], ascending=[True, True, False])
    markers_df = markers_df.groupby(["cluster", "direction"]).head(top_n)
    markers_df["rank"] = (
        markers_df.groupby(["cluster", "direction"])["score"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    markers_df = markers_df.groupby(["cluster", "direction"]).filter(lambda df: len(df) >= min_markers_per_cluster)
    return markers_df.reset_index(drop=True)


def _build_heatmap_dataframe(
    expression_matrix: sparse.spmatrix,
    obs_names: pd.Index,
    var_names: pd.Index,
    markers_df: pd.DataFrame,
    clusters: pd.Series,
    cluster_order: Sequence[str],
    z_score: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if markers_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    clusters = clusters.loc[clusters.index.intersection(obs_names)].dropna()
    if clusters.empty:
        return pd.DataFrame(), pd.DataFrame()

    cat = pd.Categorical(clusters, categories=cluster_order, ordered=True)
    ordered_clusters = pd.Series(cat, index=clusters.index).dropna().sort_values()
    if ordered_clusters.empty:
        return pd.DataFrame(), pd.DataFrame()

    markers_df = markers_df.copy()
    markers_df["cluster"] = pd.Categorical(markers_df["cluster"], categories=cluster_order, ordered=True)
    markers_df = markers_df.dropna(subset=["cluster"])
    markers_df = markers_df.sort_values(["cluster", "direction", "rank"])
    markers_df = markers_df[markers_df["marker"].isin(var_names)]
    marker_list = markers_df["marker"].tolist()
    if not marker_list:
        return pd.DataFrame(), pd.DataFrame()

    cell_indexer = obs_names.get_indexer(ordered_clusters.index)
    if np.any(cell_indexer < 0):
        raise ValueError("Failed to align clusters to expression matrix rows.")

    gene_indexer = var_names.get_indexer(marker_list)
    if np.any(gene_indexer < 0):
        raise ValueError("Failed to align markers to expression matrix columns.")

    submatrix = expression_matrix[cell_indexer, :][:, gene_indexer]
    heatmap_array = submatrix.T.toarray()
    heatmap_df = pd.DataFrame(heatmap_array, index=marker_list, columns=ordered_clusters.index)

    if z_score:
        means = heatmap_df.mean(axis=1)
        stds = heatmap_df.std(axis=1, ddof=0).replace(0, 1e-9)
        heatmap_df = heatmap_df.sub(means, axis=0).div(stds, axis=0)

    ordered_markers = markers_df.reset_index(drop=True)
    ordered_markers.index = heatmap_df.index
    return heatmap_df, ordered_markers


def _assign_group_colors(groups: Sequence[str]) -> Dict[str, str]:
    unique_groups = pd.Index(groups).unique()
    mapping = {}
    for i, name in enumerate(unique_groups):
        mapping[name] = to_hex(cm.rainbow(i / max(len(unique_groups), 1)))
    return mapping


def _plot_marker_heatmap(
    heatmap_df: pd.DataFrame,
    markers_df: pd.DataFrame,
    clusters: pd.Series,
    output_path: str,
) -> None:
    if heatmap_df.empty:
        return

    clusters = clusters.loc[heatmap_df.columns]
    cluster_colors = _assign_group_colors(clusters.values)
    cluster_to_id = pd.Series(range(len(cluster_colors)), index=list(cluster_colors.keys()))

    cell_labels_df = pd.DataFrame(
        {"cluster": [cluster_to_id[c] for c in clusters.values]},
        index=clusters.index,
    ).T
    marker_clusters = markers_df["cluster"].astype(str).values
    row_label_df = pd.DataFrame(
        {"cluster": [cluster_to_id[c] for c in marker_clusters]},
        index=markers_df.index,
    )

    plt.close("all")
    fig = plt.figure(figsize=(10, 6), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=(1, 20), height_ratios=(1, 20), wspace=0.02, hspace=0.02)
    ax_heatmap = fig.add_subplot(grid[1, 1])
    ax_col = fig.add_subplot(grid[0, 1])
    ax_row = fig.add_subplot(grid[1, 0])

    sns.heatmap(
        cell_labels_df,
        cmap=sns.color_palette(list(cluster_colors.values())),
        cbar=False,
        ax=ax_col,
        xticklabels=False,
        yticklabels=False,
    )
    ax_col.set_ylabel("")
    ax_col.set_xlabel("")

    sns.heatmap(
        heatmap_df,
        cmap=BLUE_BLACK_YELLOW,
        vmin=-3,
        vmax=3,
        ax=ax_heatmap,
        xticklabels=False,
        yticklabels=False,
        cbar=True,
        cbar_kws={"shrink": 0.5, "label": "Z-score"},
    )

    sns.heatmap(
        row_label_df[["cluster"]],
        cmap=sns.color_palette(list(cluster_colors.values())),
        cbar=False,
        ax=ax_row,
        xticklabels=False,
        yticklabels=False,
    )
    ax_row.set_ylabel("")
    ax_row.set_xlabel("")

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def find_markers_from_adata(
    adata: AnnData,
    cluster_key: str,
    *,
    output_dir: Optional[str] = None,
    layer: Optional[str] = None,
    use_raw: bool = False,
    n_markers: int = 60,
    direction: MarkerDirection = "up",
    rho_threshold: float = 0.2,
    min_markers_per_cluster: int = 3,
    lineage_order_key: str = "lineage_order",
    cluster_order: Optional[Sequence[str]] = None,
    write_outputs: bool = True,
    heatmap_filename: str = "marker_heatmap.pdf",
    marker_table_filename: str = "marker_genes.tsv",
    heatmap_table_filename: str = "marker_heatmap.tsv",
) -> MarkerOutputs:
    expression_matrix, obs_names, var_names = _get_expression_matrix(adata, layer=layer, use_raw=use_raw)
    clusters = _get_cluster_series(adata, cluster_key)

    cell_indexer = obs_names.get_indexer(clusters.index)
    if cell_indexer.size == 0 or np.any(cell_indexer < 0):
        raise ValueError("No overlapping cells between expression matrix and cluster assignments.")

    expression_matrix = expression_matrix[cell_indexer, :]
    obs_names = obs_names[cell_indexer]
    clusters = clusters.loc[obs_names]

    r_df, p_df = marker_finder(expression_matrix, clusters.tolist(), gene_names=var_names)
    markers_df = _select_unique_markers(
        r_df=r_df,
        p_df=p_df,
        direction=direction,
        rho_threshold=rho_threshold,
        top_n=n_markers,
        min_markers_per_cluster=min_markers_per_cluster,
    )

    resolved_order = _resolve_cluster_order(clusters, adata, lineage_order_key, cluster_order)
    heatmap_df, ordered_markers_df = _build_heatmap_dataframe(
        expression_matrix,
        obs_names,
        var_names,
        markers_df,
        clusters,
        resolved_order,
    )

    if write_outputs:
        if output_dir is None:
            raise ValueError("output_dir must be provided when write_outputs=True.")
        output_dir = _ensure_output_dir(output_dir)
        markers_df.to_csv(os.path.join(output_dir, marker_table_filename), sep="\t", index=False)
        heatmap_df.to_csv(os.path.join(output_dir, heatmap_table_filename), sep="\t")
        if not heatmap_df.empty and not ordered_markers_df.empty:
            _plot_marker_heatmap(
                heatmap_df=heatmap_df,
                markers_df=ordered_markers_df,
                clusters=clusters,
                output_path=os.path.join(output_dir, heatmap_filename),
            )

    return MarkerOutputs(markers=markers_df, heatmap_values=heatmap_df, cluster_assignments=clusters)


def run_marker_finder_on_file(
    h5ad_path: str,
    cluster_key: str,
    output_dir: str,
    **kwargs,
) -> MarkerOutputs:
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path)
    return find_markers_from_adata(
        adata,
        cluster_key=cluster_key,
        output_dir=output_dir,
        **kwargs,
    )


def build_arg_parser():
    import argparse

    parser = argparse.ArgumentParser(description="Run marker discovery on an AnnData (.h5ad) file.")
    parser.add_argument("--h5ad", required=True, help="Path to the input .h5ad file.")
    parser.add_argument("--cluster-key", required=True, help="Column in adata.obs containing cluster labels.")
    parser.add_argument("--output-dir", required=True, help="Directory to write marker and heatmap outputs.")
    parser.add_argument("--layer", help="Name of adata layer to use instead of .X.")
    parser.add_argument(
        "--use-raw",
        action="store_true",
        help="Use adata.raw.X as the expression matrix (ignored if --layer is provided).",
    )
    parser.add_argument("--top-n", type=int, default=60, help="Number of markers per cluster (per direction).")
    parser.add_argument(
        "--direction",
        choices=["up", "down", "both"],
        default="up",
        help="Select positive (up), negative (down), or both sets of markers.",
    )
    parser.add_argument(
        "--rho-threshold",
        type=float,
        default=0.2,
        help="Absolute Pearson correlation threshold for marker selection.",
    )
    parser.add_argument(
        "--min-markers",
        type=int,
        default=3,
        help="Minimum number of retained markers per cluster/direction.",
    )
    parser.add_argument(
        "--lineage-order-key",
        default="lineage_order",
        help="Name of adata.uns entry that stores the preferred cluster order.",
    )
    parser.add_argument("--heatmap-filename", default="marker_heatmap.pdf", help="Output figure filename.")
    parser.add_argument("--markers-tsv", default="marker_genes.tsv", help="Output TSV for markers.")
    parser.add_argument("--heatmap-tsv", default="marker_heatmap.tsv", help="Output TSV for the heatmap matrix.")
    return parser


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    run_marker_finder_on_file(
        h5ad_path=args.h5ad,
        cluster_key=args.cluster_key,
        output_dir=args.output_dir,
        layer=args.layer,
        use_raw=args.use_raw,
        n_markers=args.top_n,
        direction=args.direction,
        rho_threshold=args.rho_threshold,
        min_markers_per_cluster=args.min_markers,
        lineage_order_key=args.lineage_order_key,
        heatmap_filename=args.heatmap_filename,
        marker_table_filename=args.markers_tsv,
        heatmap_table_filename=args.heatmap_tsv,
    )


if __name__ == "__main__":
    main()
