#!/usr/bin/env python3

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# -----------------------------------------------------------
# GLOBAL MPL SETTINGS
# -----------------------------------------------------------
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['figure.facecolor'] = 'white'


# ===========================================================
# ARGUMENTS
# ===========================================================
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--umap-key", required=True)
    ap.add_argument("--cluster-key", required=True)
    ap.add_argument("--library-key", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--library-merge-tsv", default=None)
    ap.add_argument("--cluster-color-tsv", default="cluster_colors.tsv")
    return ap.parse_args()

"""
python umap_libraries.py \
  --h5ad cellHarmony_lite_with_author_clusters.h5ad \
  --umap-key X_umap_leiden \
  --cluster-key author-clusters \
  --library-key Library \
  --library-merge-tsv library_merge.tsv \
  --outdir umap_by_library
"""
# ===========================================================
# DATA LOADING / PREP
# ===========================================================
def load_library_merge(path):
    if path is None:
        return {}
    df = pd.read_csv(path, sep="\t", header=None, names=["old", "new"])
    return dict(zip(df.old, df.new))


def compute_umap_limits(coords):
    x = coords[:, 0]
    y = coords[:, 1]
    pad_x = 0.08 * (x.max() - x.min())
    pad_y = 0.08 * (y.max() - y.min())
    return (
        (x.min() - pad_x, x.max() + pad_x),
        (y.min() - pad_y, y.max() + pad_y),
    )


# ===========================================================
# COLOR HANDLING
# ===========================================================
def make_color_map(categories, cmap_name="tab20"):
    cats = list(categories)
    cmap = plt.cm.get_cmap(cmap_name, len(cats))
    return {c: cmap(i) for i, c in enumerate(cats)}


def export_color_map(color_map, out_tsv, label):
    df = pd.DataFrame(
        [(k, mcolors.to_hex(v)) for k, v in color_map.items()],
        columns=[label, "color"],
    )
    df.to_csv(out_tsv, sep="\t", index=False)


# ===========================================================
# PLOTTING
# ===========================================================
def plot_umap(
    ax,
    x,
    y,
    labels,
    color_map,
    title,
    xlim,
    ylim,
    legend_title,
):
    for lab, col in color_map.items():
        mask = labels == lab
        if mask.sum() == 0:
            continue
        ax.scatter(
            x[mask],
            y[mask],
            s=4,
            c=[col],
            linewidths=0,
            alpha=0.8,
            label=lab,
        )

    ax.set_title(title, fontsize=10)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    ax.set_xlabel("UMAP-X")
    ax.set_ylabel("UMAP-Y")

    ax.legend(
        title=legend_title,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=2,
        fontsize=7,
        title_fontsize=8,
    )


# ===========================================================
# MAIN
# ===========================================================
def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)

    coords = adata.obsm[args.umap_key]
    x_all = coords[:, 0]
    y_all = coords[:, 1]
    xlim, ylim = compute_umap_limits(coords)

    clusters = adata.obs[args.cluster_key].astype(str)

    lib_merge = load_library_merge(args.library_merge_tsv)
    libraries = (
        adata.obs[args.library_key]
        .astype(str)
        .replace(lib_merge)
    )

    # ---- CLUSTER COLORS ----
    cluster_colors = make_color_map(clusters.unique())
    export_color_map(cluster_colors, args.cluster_color_tsv, "cluster")

    # ---- LIBRARY COLORS ----
    library_colors = make_color_map(libraries.unique(), cmap_name="tab20")

    # =======================================================
    # COMBINED (BY CLUSTER)
    # =======================================================
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_umap(
        ax,
        x_all,
        y_all,
        clusters.values,
        cluster_colors,
        "Combined",
        xlim,
        ylim,
        "Cluster",
    )
    fig.savefig(outdir / "UMAP_combined.pdf", bbox_inches="tight")
    plt.close(fig)

    # =======================================================
    # PER-LIBRARY (BY CLUSTER)
    # =======================================================
    for lib in sorted(libraries.unique()):
        mask = libraries == lib
        fig, ax = plt.subplots(figsize=(6, 6))
        plot_umap(
            ax,
            x_all[mask],
            y_all[mask],
            clusters[mask].values,
            cluster_colors,
            lib,
            xlim,
            ylim,
            "Cluster",
        )
        fig.savefig(outdir / f"UMAP_{lib}.pdf", bbox_inches="tight")
        plt.close(fig)

    # =======================================================
    # ALWAYS PRODUCED: UMAP BY LIBRARY
    # =======================================================
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_umap(
        ax,
        x_all,
        y_all,
        libraries.values,
        library_colors,
        "UMAP by Library",
        xlim,
        ylim,
        "Library",
    )
    fig.savefig(outdir / "UMAP_by_library.pdf", bbox_inches="tight")
    plt.close(fig)

    print("UMAPs written to:", outdir)


if __name__ == "__main__":
    main()
