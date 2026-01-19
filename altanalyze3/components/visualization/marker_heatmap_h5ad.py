#!/usr/bin/env python3
import argparse
import os
import sys
import warnings

import matplotlib
import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning
import scanpy as sc
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
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


def _compute_marker_stats(adata, cluster_key, markers_df, fdr_df, use_raw, layer):
    if markers_df.empty:
        return pd.DataFrame(columns=["Gene", "Fold", "Query Exp", "Ref Exp", "FDR p-value", "cluster"])

    genes = markers_df["gene"].tolist()
    expr_df = sc.get.obs_df(adata, keys=genes, use_raw=use_raw, layer=layer)
    cluster_series = adata.obs[cluster_key].astype(str)
    expr_df[cluster_key] = cluster_series.values

    sums = expr_df.groupby(cluster_key)[genes].sum()
    counts = expr_df.groupby(cluster_key).size()
    total_sum = expr_df[genes].sum(axis=0)
    total_count = float(expr_df.shape[0])

    rows = []
    eps = 1e-9
    for row in _progress(
        markers_df.itertuples(index=False),
        total=len(markers_df),
        desc="Computing marker statistics",
    ):
        gene = row.gene
        cluster = str(row.cluster)
        if cluster not in sums.index:
            continue
        query_sum = float(sums.loc[cluster, gene])
        query_count = float(counts.loc[cluster])
        query_mean = query_sum / max(query_count, 1.0)
        ref_sum = float(total_sum[gene] - query_sum)
        ref_count = max(total_count - query_count, 1.0)
        ref_mean = ref_sum / ref_count
        fold = float(np.log2((query_mean + eps) / (ref_mean + eps)))
        fdr = np.nan
        if fdr_df is not None and gene in fdr_df.index and cluster in fdr_df.columns:
            fdr = float(fdr_df.loc[gene, cluster])
        rows.append(
            {
                "Gene": gene,
                "Fold": fold,
                "Query Exp": query_mean,
                "Ref Exp": ref_mean,
                "FDR p-value": fdr,
                "cluster": cluster,
            }
        )

    return pd.DataFrame(rows)


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


def _build_heatmap(adata, cluster_key, genes, cluster_order, use_raw, layer):
    if not genes:
        return pd.DataFrame(), pd.DataFrame(), [], []

    if use_raw:
        if adata.raw is None:
            raise ValueError("Requested use_raw=True but adata.raw is not set.")
        gene_reference = adata.raw.var_names
    else:
        gene_reference = adata.var_names

    genes_present = [g for g in genes if g in gene_reference]
    if not genes_present:
        return pd.DataFrame(), pd.DataFrame(), [], []

    expr_df = sc.get.obs_df(adata, keys=genes_present, use_raw=use_raw, layer=layer)
    ordered_cells, cluster_counts = _order_cells_by_cluster(adata, cluster_key, cluster_order)
    if not ordered_cells:
        return pd.DataFrame(), pd.DataFrame(), [], []

    expr_df = expr_df.loc[ordered_cells]
    heatmap_raw_df = expr_df.T
    heatmap_col_df = _zscore_columns(heatmap_raw_df)
    heatmap_row_df = _zscore_rows(heatmap_raw_df)
    return heatmap_row_df, heatmap_col_df, cluster_counts, ordered_cells


def _cluster_color_map(cluster_order):
    palette = plt.get_cmap("tab20")
    colors = [palette(i % palette.N) for i in range(len(cluster_order))]
    return {cluster: colors[i] for i, cluster in enumerate(cluster_order)}, ListedColormap(colors)


def _plot_heatmap(heatmap_df, output_path, cluster_counts, cluster_order, column_clusters, row_clusters):
    if heatmap_df.empty:
        raise ValueError("Heatmap dataframe is empty.")

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    cluster_count = max(1, len(cluster_counts))
    fixed_height = 7
    fig_w = max(6.0, 0.25 * cluster_count)
    fig, ax = plt.subplots(figsize=(fig_w, fixed_height))

    vmin, vmax = -3.0, 3.0
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

    if centers:
        ax.set_xticks(centers)
        ax.set_xticklabels(labels, rotation=45, ha="left", fontsize=7)
    else:
        ax.set_xticks([])

    ax.xaxis.tick_top()
    ax.tick_params(axis="x", bottom=False, top=True, labelbottom=False, labeltop=True, length=0, pad=8)

    n_genes = heatmap_df.shape[0]
    if n_genes <= 25:
        gene_fontsize = 9
    elif n_genes <= 60:
        gene_fontsize = 8
    elif n_genes <= 120:
        gene_fontsize = 7
    elif n_genes <= 200:
        gene_fontsize = 6
    else:
        gene_fontsize = 5

    ax.set_yticks(np.arange(heatmap_df.shape[0]))
    ax.set_yticklabels(heatmap_df.index, fontsize=gene_fontsize)
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", left=False, right=True, labelleft=False, labelright=True, pad=2)

    for boundary in boundaries:
        ax.axvline(boundary, color="white", linewidth=0.5)

    divider = make_axes_locatable(ax)
    ax_top = divider.append_axes("top", size="1.5%", pad=0.02)
    ax_left = divider.append_axes("left", size="3%", pad=0.02)
    cax = inset_axes(
        ax,
        width="35%",
        height="3%",
        loc="lower center",
        bbox_to_anchor=(0.0, -0.06, 1.0, 1.0),
        bbox_transform=ax.transAxes,
        borderpad=0,
    )

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

    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.ax.set_xlim(vmin, vmax)
    cbar.set_label("Norm Exp (z-score)")
    cbar.set_ticks([])
    cbar.ax.text(-0.08, 0.5, f"{vmin:.0f}", ha="right", va="center", transform=cbar.ax.transAxes)
    cbar.ax.text(1.08, 0.5, f"{vmax:.0f}", ha="left", va="center", transform=cbar.ax.transAxes)

    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


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


def _markerfinder_stats(adata, cluster_key, use_raw, layer):
    marker_finder = _load_marker_finder()
    matrix, gene_names = _get_expression_matrix(adata, use_raw, layer)
    clusters = adata.obs[cluster_key].astype(str).tolist()
    r_df, p_df = marker_finder.marker_finder(matrix, clusters, gene_names=gene_names)
    return p_df, r_df



def main():
    parser = argparse.ArgumentParser(
        description="Create a marker heatmap from an h5ad file using unique markers per cluster."
    )
    parser.add_argument("--h5ad", required=True, help="Input h5ad file.")
    parser.add_argument("--cluster-key", required=True, help="obs column with cluster labels.")
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
        "--cells-per-cluster",
        type=int,
        default=50,
        help="Number of cells to display per cluster (default 50). Use 0 for all cells.",
    )
    parser.add_argument("--seed", type=int, default=0, help="Random seed for cell sampling.")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out) or "."
    logs_dir = os.path.join(out_dir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(args.out))[0]
    log_path = os.path.join(logs_dir, f"{base_name}.log")
    markers_tsv = args.markers_tsv or os.path.join(out_dir, f"{base_name}_markers.tsv")
    heatmap_tsv = args.heatmap_tsv or os.path.join(out_dir, f"{base_name}_fold_matrix.tsv")
    heatmap_column_tsv = args.heatmap_column_tsv or os.path.join(out_dir, f"{base_name}_exp_matrix.tsv")

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
        print(f"[INFO] Heatmap TSV (fold matrix): {heatmap_tsv}")
        print(f"[INFO] Heatmap TSV (expression matrix): {heatmap_column_tsv}")

        layer = args.layer
        if args.use_raw and layer:
            print("Warning: --layer ignored when --use-raw is set.")
            layer = None

        print(f"[INFO] Reading AnnData from: {args.h5ad}")
        adata = sc.read_h5ad(args.h5ad)
        print(f"[INFO] Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes.")

        if args.cluster_key not in adata.obs:
            raise KeyError(f"Cluster key '{args.cluster_key}' not found in adata.obs.")

        if args.use_raw and adata.raw is None:
            raise ValueError("Requested --use-raw but adata.raw is not set.")

        if layer and layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers.")

        clusters = adata.obs[args.cluster_key].astype(str)
        lineage_order = _coerce_lineage_order(adata.uns.get("lineage_order", None))
        cluster_order = _resolve_cluster_order(clusters, lineage_order)
        print(f"[INFO] Resolved {len(cluster_order)} clusters.")

        adata.obs[args.cluster_key] = pd.Categorical(clusters, categories=cluster_order, ordered=True)

        effect_df = None
        fdr_df = None
        if args.marker_method == "markerfinder":
            print(
                f"[INFO] Running markerFinder using {'raw' if args.use_raw else (layer or 'X')} matrix."
            )
            pvals, effect_df = _markerfinder_stats(adata, args.cluster_key, args.use_raw, layer)
            fdr_df = pvals.apply(_bh_fdr, axis=0)
        else:
            n_genes = adata.raw.n_vars if args.use_raw and adata.raw is not None else adata.n_vars
            print(f"[INFO] Running rank_genes_groups for {len(cluster_order)} clusters with {n_genes} genes.")

            sc.tl.rank_genes_groups(
                adata,
                groupby=args.cluster_key,
                method=args.method,
                use_raw=args.use_raw,
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

        selected = _select_unique_markers(pvals, cluster_order, args.top_n, effect_df=effect_df, pval_threshold=0.05)
        print(f"[INFO] Selected {selected.shape[0]} markers after filtering.")

        if selected.empty:
            raise ValueError("No markers were selected. Check inputs and parameters.")

        if args.cells_per_cluster and args.cells_per_cluster > 0:
            print(
                f"[INFO] Downsampling to {args.cells_per_cluster} cells per cluster (seed={args.seed})."
            )
        else:
            print("[INFO] Using all cells (no downsampling).")
        adata_plot = downsample_cells_per_group(
            adata,
            args.cluster_key,
            cells_per_cluster=args.cells_per_cluster,
            seed=args.seed,
            group_order=cluster_order,
        )
        print(f"[INFO] Heatmap will use {adata_plot.n_obs} cells.")
        adata_plot.obs[args.cluster_key] = pd.Categorical(
            adata_plot.obs[args.cluster_key].astype(str),
            categories=cluster_order,
            ordered=True,
        )

        selected_genes = selected["gene"].tolist()
        print(f"[INFO] Building heatmap matrix for {len(selected_genes)} markers.")
        heatmap_df, heatmap_col_df, cluster_counts, ordered_cells = _build_heatmap(
            adata_plot,
            args.cluster_key,
            selected_genes,
            cluster_order,
            args.use_raw,
            layer,
        )

        if heatmap_df.empty:
            raise ValueError("Heatmap data is empty after filtering genes.")

        column_clusters = adata_plot.obs[args.cluster_key].astype(str).loc[ordered_cells].tolist()
        row_clusters = selected.set_index("gene").loc[heatmap_df.index, "cluster"].astype(str).tolist()
        print("[INFO] Rendering heatmap.")
        _plot_heatmap(heatmap_df, args.out, cluster_counts, cluster_order, column_clusters, row_clusters)
        print(f"Saved marker heatmap to: {args.out}")

        print("[INFO] Computing marker statistics for TSV.")
        marker_stats = _compute_marker_stats(
            adata,
            args.cluster_key,
            selected,
            fdr_df,
            args.use_raw,
            layer,
        )
        marker_stats.to_csv(markers_tsv, sep="\t", index=False)
        print(f"Saved marker stats TSV to: {markers_tsv}")

        heatmap_tsv_df = heatmap_df.copy()
        heatmap_tsv_df.index = [f"{c}:{g}" for c, g in zip(row_clusters, heatmap_tsv_df.index)]
        heatmap_tsv_df.columns = [f"{c}:{b}" for c, b in zip(column_clusters, ordered_cells)]
        heatmap_tsv_df.to_csv(heatmap_tsv, sep="\t")
        print(f"Saved heatmap matrix TSV to: {heatmap_tsv}")

        heatmap_col_tsv_df = heatmap_col_df.copy()
        heatmap_col_tsv_df.index = [f"{c}:{g}" for c, g in zip(row_clusters, heatmap_col_tsv_df.index)]
        heatmap_col_tsv_df.columns = [f"{c}:{b}" for c, b in zip(column_clusters, ordered_cells)]
        heatmap_col_tsv_df.to_csv(heatmap_column_tsv, sep="\t")
        print(f"Saved expression matrix TSV to: {heatmap_column_tsv}")

        print(f"Saved log to: {log_path}")
    finally:
        sys.stdout, sys.stderr = original_stdout, original_stderr
        log_file.close()


if __name__ == "__main__":

    main()
