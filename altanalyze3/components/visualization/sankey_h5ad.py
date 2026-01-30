#!/usr/bin/env python3

import argparse
import anndata as ad
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px


# ===========================================================
# ARGUMENTS
# ===========================================================
def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--cluster-key", required=True)
    ap.add_argument("--library-key", required=True)
    ap.add_argument("--library-order-tsv", required=True)
    ap.add_argument("--cluster-color-tsv", default=None)

    ap.add_argument("--flow-mode", choices=["polygon", "bezier"], default="bezier")
    ap.add_argument("--draw-block-borders", action="store_true")
    ap.add_argument("--draw-flow-borders", action="store_true")
    ap.add_argument("--border-color", default="DARK GREY")
    ap.add_argument("--border-width", type=float, default=0.5)

    ap.add_argument("--output-prefix", required=True)

    return ap.parse_args()

"""
python alluvial_plot.py \
  --h5ad cellHarmony_lite_with_author_clusters.h5ad \
  --cluster-key author-clusters \
  --library-key Library \
  --library-order-tsv library_order.tsv \
  --cluster-color-tsv cluster_colors.tsv \
  --flow-mode bezier \
  --draw-block-borders \
  --draw-flow-borders \
  --output-prefix alluvial_with_modes
"""

# ===========================================================
# UTILITIES
# ===========================================================
def load_library_order(path):
    return (
        pd.read_csv(path, sep="\t", header=None)[0]
        .astype(str)
        .tolist()
    )


def load_cluster_colors(path, clusters):
    if path is None:
        palette = px.colors.qualitative.Vivid
        return {cl: palette[i % len(palette)] for i, cl in enumerate(clusters)}

    df = pd.read_csv(path, sep="\t")
    if not {"cluster", "color"}.issubset(df.columns):
        raise ValueError("Color TSV must contain columns: cluster, color")

    colors = dict(zip(df["cluster"].astype(str), df["color"]))
    missing = [c for c in clusters if c not in colors]
    if missing:
        raise ValueError("Missing colors for clusters: " + ", ".join(missing))
    return colors


def extract_lineage_order(adata, clusters):
    raw = adata.uns.get("lineage_order", None)
    if raw is None:
        return clusters

    try:
        lineage = [str(x) for x in raw]
    except Exception:
        return clusters

    ordered = [c for c in lineage if c in clusters]
    missing = [c for c in clusters if c not in ordered]

    if missing:
        print("[INFO] Missing from lineage_order:", ", ".join(missing))

    return ordered + missing


def bezier_curve(p0, p1, p2, p3, n=40):
    t = np.linspace(0, 1, n)
    return (
        (1 - t)**3 * p0 +
        3*(1 - t)**2 * t * p1 +
        3*(1 - t) * t**2 * p2 +
        t**3 * p3
    )


# ===========================================================
# CORE LOGIC (REMOTE CALLABLE)
# ===========================================================
def run_alluvial(
    h5ad,
    cluster_key,
    library_key,
    library_order,
    cluster_colors,
    flow_mode,
    draw_block_borders,
    draw_flow_borders,
    border_color,
    border_width,
    output_prefix,
):
    adata = ad.read_h5ad(h5ad)

    df = adata.obs[[library_key, cluster_key]].copy()
    df[library_key] = pd.Categorical(
        df[library_key],
        categories=library_order,
        ordered=True,
    )

    clusters = sorted(df[cluster_key].astype(str).unique())
    clusters = extract_lineage_order(adata, clusters)

    count_matrix = (
        df.groupby([library_key, cluster_key])
          .size()
          .unstack(fill_value=0)
          .reindex(index=library_order, columns=clusters, fill_value=0)
    )

    fraction_matrix = count_matrix.div(count_matrix.sum(axis=1), axis=0)

    fig = go.Figure()
    x_positions = np.arange(len(library_order))
    cluster_pos = {}

    # -------------------------------------------------------
    # STACKED BARS
    # -------------------------------------------------------
    for i, lib in enumerate(library_order):
        bottom = 0.0
        for cl in clusters:
            height = fraction_matrix.loc[lib, cl]
            top = bottom + height
            cluster_pos[(lib, cl)] = (bottom, top)

            fig.add_shape(
                type="rect",
                x0=i - 0.35,
                x1=i + 0.35,
                y0=bottom,
                y1=top,
                fillcolor=cluster_colors[cl],
                line=dict(
                    width=border_width if draw_block_borders else 0,
                    color=border_color,
                ),
            )
            bottom = top

    # -------------------------------------------------------
    # LABELS
    # -------------------------------------------------------
    for cl in clusters:
        lib_max = fraction_matrix[cl].idxmax()
        b, t = cluster_pos[(lib_max, cl)]
        fig.add_annotation(
            x=library_order.index(lib_max),
            y=(b + t) / 2,
            text=cl,
            showarrow=False,
            font=dict(size=15, color="black"),
        )

    # -------------------------------------------------------
    # FLOWS
    # -------------------------------------------------------
    alpha = 0.55
    bez_n = 50
    push = 0.25

    for k in range(len(library_order) - 1):
        L = library_order[k]
        R = library_order[k + 1]

        xL, xR = k + 0.35, k + 1 - 0.35
        cL, cR = xL + push, xR - push

        for cl in clusters:
            Lb, Lt = cluster_pos[(L, cl)]
            Rb, Rt = cluster_pos[(R, cl)]

            if Lt == Lb and Rt == Rb:
                continue

            rgba = cluster_colors[cl].replace(
                "rgb", "rgba"
            ).replace(")", f",{alpha})")

            if flow_mode == "polygon":
                X = [xL, xR, xR, xL]
                Y = [Lt, Rt, Rb, Lb]
            else:
                top_x = bezier_curve(xL, cL, cR, xR, bez_n)
                top_y = bezier_curve(Lt, Lt, Rt, Rt, bez_n)
                bot_x = bezier_curve(xR, cR, cL, xL, bez_n)
                bot_y = bezier_curve(Rb, Rb, Lb, Lb, bez_n)
                X = np.concatenate([top_x, bot_x])
                Y = np.concatenate([top_y, bot_y])

            fig.add_trace(go.Scatter(
                x=X,
                y=Y,
                fill="toself",
                mode="lines",
                line=dict(
                    width=border_width if draw_flow_borders else 0,
                    color=border_color,
                ),
                fillcolor=rgba,
                hoverinfo="skip",
                showlegend=False,
            ))

    # -------------------------------------------------------
    # LEGEND
    # -------------------------------------------------------
    for cl in clusters:
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode="markers",
            marker=dict(size=12, color=cluster_colors[cl]),
            name=cl,
        ))

    # -------------------------------------------------------
    # LAYOUT
    # -------------------------------------------------------
    fig.update_layout(
        title="Alluvial View of Cell Population Composition Across Libraries",
        width=1500,
        height=900,
        plot_bgcolor="white",
        xaxis=dict(
            tickmode="array",
            tickvals=x_positions,
            ticktext=library_order,
            range=[-1.6, len(library_order) - 0.25],
        ),
        yaxis=dict(range=[0, 1], showticklabels=False),
    )

    html = output_prefix + ".html"
    svg = output_prefix + ".svg"

    fig.write_html(html)
    fig.write_image(svg, format="svg", width=1500, height=900)

    print("Saved:", html)
    print("Saved:", svg)


# ===========================================================
# ENTRY POINT
# ===========================================================
def main():
    args = parse_args()

    library_order = load_library_order(args.library_order_tsv)

    adata_tmp = ad.read_h5ad(args.h5ad, backed="r")
    clusters_tmp = adata_tmp.obs[args.cluster_key].astype(str).unique().tolist()

    cluster_colors = load_cluster_colors(
        args.cluster_color_tsv,
        clusters_tmp,
    )

    run_alluvial(
        h5ad=args.h5ad,
        cluster_key=args.cluster_key,
        library_key=args.library_key,
        library_order=library_order,
        cluster_colors=cluster_colors,
        flow_mode=args.flow_mode,
        draw_block_borders=args.draw_block_borders,
        draw_flow_borders=args.draw_flow_borders,
        border_color=args.border_color,
        border_width=args.border_width,
        output_prefix=args.output_prefix,
    )


if __name__ == "__main__":
    main()
