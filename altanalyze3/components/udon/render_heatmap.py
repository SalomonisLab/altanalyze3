#!/usr/bin/env python3
"""Render a raster PNG (+ compact raster PDF) from a UDON marker_heatmap.txt
(AltAnalyze FinalMarkerHeatmap format: row 'column_clusters-flat' = per-column
cluster, col 'row_clusters-flat' = per-row cluster; values already row-median
centered). The default plot_markers_df makes a per-cell VECTOR PDF that balloons
to >200 MB at this scale; imshow rasterization keeps it small.

Usage: python3 render_heatmap.py <marker_heatmap.txt> <out_prefix> [max_cols]
"""
import sys, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
CYAN_YELLOW = LinearSegmentedColormap.from_list("cyan_yellow", ["#00FFFF", "#FFFF00"])


def load_marker_heatmap(path):
    df = pd.read_csv(path, sep="\t", index_col=0, low_memory=False)
    col_clusters = df.loc["column_clusters-flat"]            # per-column cluster (first data row)
    row_clusters = df["row_clusters-flat"]                   # per-row cluster (first data col)
    mat = df.drop(index="column_clusters-flat").drop(columns="row_clusters-flat")
    mat = mat.astype(float)
    rc = row_clusters.drop("column_clusters-flat")
    cc = col_clusters.drop("row_clusters-flat")
    return mat, rc.astype(str), cc.astype(str)


def cluster_boundaries(labels):
    labels = list(labels)
    return [i for i in range(1, len(labels)) if labels[i] != labels[i - 1]]


def render(path, prefix, max_cols=4000):
    mat, rc, cc = load_marker_heatmap(path)
    print(f"matrix: {mat.shape[0]} markers x {mat.shape[1]} pseudobulks; "
          f"{rc.nunique()} row-clusters, {cc.nunique()} col-clusters")
    M = mat.to_numpy()
    # downsample columns for a legible figure if very wide (keeps cluster order)
    if M.shape[1] > max_cols:
        idx = np.linspace(0, M.shape[1] - 1, max_cols).astype(int)
        M = M[:, idx]; cc = cc.iloc[idx]
        print(f"  downsampled columns to {max_cols} for display (full matrix in marker_heatmap.txt)")
    vmax = float(np.nanpercentile(np.abs(M), 98)) or 1.0
    h = max(6, M.shape[0] * 0.045); w = 16
    fig, ax = plt.subplots(figsize=(w, h))
    im = ax.imshow(M, aspect="auto", cmap=CYAN_YELLOW, vmin=-vmax, vmax=vmax, interpolation="nearest")
    for b in cluster_boundaries(cc):
        ax.axvline(b - 0.5, color="#000000", lw=0.3)
    for b in cluster_boundaries(rc):
        ax.axhline(b - 0.5, color="#000000", lw=0.3)
    ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(mat.index, fontsize=2)
    ax.set_xticks([]); ax.set_xlabel(f"{M.shape[1]} pseudobulks (cluster-ordered)")
    ax.set_ylabel(f"{M.shape[0]} marker genes")
    ax.set_title("UDON marker heatmap (row-median centered; cyan=low, yellow=high)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.012, pad=0.01); cbar.ax.tick_params(labelsize=6)
    fig.tight_layout()
    png = prefix + ".png"; pdf = prefix + ".pdf"
    fig.savefig(png, dpi=200)
    fig.savefig(pdf, dpi=200)            # imshow stays raster -> small PDF
    plt.close(fig)
    print("wrote", png, "and", pdf)


if __name__ == "__main__":
    render(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 4000)
