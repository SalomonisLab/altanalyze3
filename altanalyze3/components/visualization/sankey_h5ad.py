#!/usr/bin/env python3
import anndata as ad
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import matplotlib as mpl
import matplotlib.pyplot as plt

# ===========================================================
# USER INPUTS
# ===========================================================
INPUT_H5AD = "cellHarmony_lite_with_author_clusters.h5ad"
OBS_CLUSTER_FIELD = "author-clusters"
OBS_LIBRARY_FIELD = "Library"

# Choose flow mode: "polygon" or "bezier"
flow_mode = "bezier"

# Optional: draw black outlines around bars and flows
DRAW_BLOCK_BORDERS = True
DRAW_FLOW_BORDERS  = True
BORDER_COLOR = "DARK GREY"
BORDER_WIDTH = 0.5

LIBRARY_ORDER = [
    "LWAT-CD45pos",
    "EWAT-CD45pos",
    "HINI-LWAT-CD45pos",
    "HINI-EWAT-CD45pos",
]

OUTPUT = "alluvial_with_modes.html"

# ===========================================================
# LOAD DATA
# ===========================================================
adata = ad.read_h5ad(INPUT_H5AD)
df = adata.obs[[OBS_LIBRARY_FIELD, OBS_CLUSTER_FIELD]].copy()

df[OBS_LIBRARY_FIELD] = pd.Categorical(
    df[OBS_LIBRARY_FIELD],
    categories=LIBRARY_ORDER,
    ordered=True
)

# ===========================================================
# BUILD COUNT + FRACTION MATRICES
# ===========================================================
count_matrix = (
    df.groupby([OBS_LIBRARY_FIELD, OBS_CLUSTER_FIELD])
      .size()
      .unstack(fill_value=0)
      .reindex(index=LIBRARY_ORDER, columns=sorted(df[OBS_CLUSTER_FIELD].unique()), fill_value=0)
)

fraction_matrix = count_matrix.div(count_matrix.sum(axis=1), axis=0)

# ===========================================================
# CLUSTER ORDERING (LINEAGE ORDER IF AVAILABLE)
# ===========================================================
raw_lineage = adata.uns.get("lineage_order", None)

clusters = list(sorted(df[OBS_CLUSTER_FIELD].unique()))
lineage = None

if raw_lineage is not None:
    try:
        if isinstance(raw_lineage, (list, tuple)):
            lineage = [str(x) for x in raw_lineage]
        elif hasattr(raw_lineage, "shape") and len(raw_lineage.shape) == 1:
            lineage = [str(x) for x in raw_lineage]
        elif hasattr(raw_lineage, "tolist"):
            lineage = [str(x) for x in raw_lineage.tolist()]
        elif hasattr(raw_lineage, "__iter__") and not isinstance(raw_lineage, str):
            lineage = [str(x) for x in raw_lineage]
    except:
        lineage = None

if lineage:
    lineage_filtered = [c for c in lineage if c in clusters]
    missing = [c for c in clusters if c not in lineage_filtered]
    clusters = lineage_filtered + missing

    print("\n[INFO] Using lineage_order from adata.uns:")
    for c in clusters:
        print("   ", c)
else:
    print("[WARNING] Using alphabetical ordering for clusters.")

# ===========================================================
# COLORS
# ===========================================================
palette = px.colors.qualitative.Vivid
#palette = px.colors.qualitative.Dark24
colors = {cl: palette[i % len(palette)] for i, cl in enumerate(clusters)}

# ===========================================================
# FIGURE SETUP
# ===========================================================
fig = go.Figure()
x_positions = np.arange(len(LIBRARY_ORDER))
cluster_positions = {}

# ===========================================================
# DRAW STACKED BARS (WITH OPTIONAL BLACK BORDERS)
# ===========================================================
for lib_idx, lib in enumerate(LIBRARY_ORDER):
    bottom = 0.0
    for cl in clusters:
        height = fraction_matrix.loc[lib, cl]
        top = bottom + height

        cluster_positions[(lib, cl)] = (bottom, top)

        fig.add_shape(
            type="rect",
            x0=lib_idx - 0.35,
            x1=lib_idx + 0.35,
            y0=bottom,
            y1=top,
            fillcolor=colors[cl],
            line=dict(
                width=BORDER_WIDTH if DRAW_BLOCK_BORDERS else 0,
                color=BORDER_COLOR
            )
        )

        bottom = top

# ===========================================================
# LABELS PLACED ON TOP OF THE FLOWS (ALWAYS VISIBLE)
# ===========================================================
for cl in clusters:
    lib_max = fraction_matrix[cl].idxmax()
    bottom, top = cluster_positions[(lib_max, cl)]
    mid_y = (bottom + top) / 2

    lib_idx = LIBRARY_ORDER.index(lib_max)
    mid_x = lib_idx

    fig.add_annotation(
        x=mid_x,
        y=mid_y,
        text=cl,
        showarrow=False,
        font=dict(size=15, color="black"),
        xanchor="center",
        yanchor="middle"
    )

# ===========================================================
# FLOW RENDERING (POLYGON + BEZIER) WITH OPTIONAL BORDERS
# ===========================================================
def bezier_curve(p0, p1, p2, p3, n=40):
    t = np.linspace(0, 1, n)
    return (
        (1 - t)**3 * p0 +
        3*(1 - t)**2 * t * p1 +
        3*(1 - t) * t**2 * p2 +
        t**3 * p3
    )

flow_alpha = 0.55
bezier_n = 50
bezier_push = 0.25

for k in range(len(LIBRARY_ORDER) - 1):
    L = LIBRARY_ORDER[k]
    R = LIBRARY_ORDER[k + 1]

    xL = k + 0.35
    xR = k + 1 - 0.35

    ctrlL = xL + bezier_push
    ctrlR = xR - bezier_push

    for cl in clusters:
        L_bottom, L_top = cluster_positions[(L, cl)]
        R_bottom, R_top = cluster_positions[(R, cl)]

        if L_top == L_bottom and R_top == R_bottom:
            continue

        rgba = colors[cl].replace("rgb", "rgba").replace(")", f",{flow_alpha})")

        # ---- POLYGON MODE ----
        if flow_mode == "polygon":
            X = [xL, xR, xR, xL]
            Y = [L_top, R_top, R_bottom, L_bottom]

            fig.add_trace(go.Scatter(
                x=X,
                y=Y,
                fill="toself",
                mode="lines",
                line=dict(
                    width=BORDER_WIDTH if DRAW_FLOW_BORDERS else 0,
                    color=BORDER_COLOR
                ),
                fillcolor=rgba,
                hoverinfo="skip",
                showlegend=False
            ))
            continue

        # ---- BEZIER MODE ----
        if flow_mode == "bezier":
            top_x = bezier_curve(xL, ctrlL, ctrlR, xR, n=bezier_n)
            top_y = bezier_curve(L_top, L_top, R_top, R_top, n=bezier_n)

            bot_x = bezier_curve(xR, ctrlR, ctrlL, xL, n=bezier_n)
            bot_y = bezier_curve(R_bottom, R_bottom, L_bottom, L_bottom, n=bezier_n)

            X = np.concatenate([top_x, bot_x])
            Y = np.concatenate([top_y, bot_y])

            fig.add_trace(go.Scatter(
                x=X,
                y=Y,
                fill="toself",
                mode="lines",
                line=dict(
                    width=BORDER_WIDTH if DRAW_FLOW_BORDERS else 0,
                    color=BORDER_COLOR
                ),
                fillcolor=rgba,
                hoverinfo="skip",
                showlegend=False
            ))
            continue

        raise ValueError(f"Unknown flow_mode={flow_mode}")

# ===========================================================
# LEGEND
# ===========================================================
for cl in clusters:
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode="markers",
        marker=dict(size=12, color=colors[cl]),
        name=cl
    ))

# ===========================================================
# LAYOUT
# ===========================================================
fig.update_layout(
    title="Alluvial View of Cell Population Composition Across Libraries",
    width=1500,
    height=900,
    plot_bgcolor="white",
    xaxis=dict(
        tickmode="array",
        tickvals=x_positions,
        ticktext=LIBRARY_ORDER,
        range=[-1.6, len(LIBRARY_ORDER) - 0.25]
    ),
    yaxis=dict(range=[0, 1], showticklabels=False),
)

fig.write_html(OUTPUT)
import plotly.io as pio

#OUTPUT_PDF = OUTPUT.replace(".html", ".pdf")
#pio.write_image(fig, OUTPUT_PDF, format="pdf", width=1500, height=900, scale=2)
#print("Saved PDF:", OUTPUT_PDF)

OUTPUT_SVG = OUTPUT.replace(".html", ".svg")
fig.write_image(OUTPUT_SVG, format="svg", width=1500, height=900)
print("Saved SVG:", OUTPUT_SVG)

