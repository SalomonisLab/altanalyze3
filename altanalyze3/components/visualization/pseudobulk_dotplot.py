#!/usr/bin/env python3
"""Replicate (pseudobulk) dot plots: one panel per (gene, cell state), one column per group, one point per replicate.

Input is a pseudobulk h5ad written by altanalyze3.components.aggregate.pseudobulk_h5ad (X = counts / total counts
per pseudobulk; obs carries the cluster and sample columns). Each panel draws log2(1 + CP10k) of one gene in one cell
state for every sample of every group, a 2-point mean line per group, and n per group. Text in a panel's
`mark_<group>` column (from the caller's own test, never recomputed here) is printed above that group's column.
Points are vector EllipseCollection circles, every line is a 2-point path, and clipping is off.

python -m altanalyze3.components.visualization.pseudobulk_dotplot \
  --pseudobulk-h5ad pb.h5ad --cluster-col Hs-MarrowAtlas-L3M --sample-col replicate_unit \
  --groups-tsv groups.tsv --group-order patient,MISTRG,control --panels-tsv panels.tsv --out dotplots.pdf
"""
from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
from matplotlib.colors import to_rgba

DEFAULT_COLORS = ["#C0392B", "#2471A3", "#5D6D7E", "#B9770E", "#117A65"]


def apply_rcparams():
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
    plt.rcParams["axes.linewidth"] = 0.5
    plt.rcParams["xtick.major.width"] = 0.5
    plt.rcParams["ytick.major.width"] = 0.5


def load_values(pb_path, cluster_col, sample_col, panels, groups):
    """Long table of log2(1 + CP10k) per (panel, sample); a sample without a pseudobulk in the state is kept as NaN."""
    pb = ad.read_h5ad(pb_path)
    for col in (cluster_col, sample_col):
        if col not in pb.obs:
            raise KeyError(f"{col} not in pseudobulk obs")
    missing_genes = sorted(set(panels.gene) - set(pb.var_names))
    if missing_genes:
        raise KeyError(f"genes absent from the pseudobulk matrix: {missing_genes}")
    key = pb.obs[cluster_col].astype(str) + "|" + pb.obs[sample_col].astype(str)
    if not key.is_unique:
        raise ValueError("more than one pseudobulk per (cell state, sample)")
    row_of = pd.Series(np.arange(pb.n_obs), index=key.to_numpy())
    rows = []
    for p in panels.itertuples():
        j = pb.var_names.get_loc(p.gene)
        for s in groups.itertuples():
            i = row_of.get(f"{p.cell_state}|{s.sample}")
            if i is None:
                rows.append(dict(panel=p.Index, gene=p.gene, cell_state=p.cell_state, group=s.group, sample=s.sample,
                                 n_cells=0, value=np.nan))
                continue
            x = pb.X[i, j]
            frac = float(x.toarray().ravel()[0]) if hasattr(x, "toarray") else float(x)
            rows.append(dict(panel=p.Index, gene=p.gene, cell_state=p.cell_state, group=s.group, sample=s.sample,
                             n_cells=int(pb.obs["n_cells"].iloc[i]) if "n_cells" in pb.obs else np.nan,
                             value=float(np.log2(1.0 + 1e4 * frac))))
    return pd.DataFrame(rows)


def draw(values, panels, group_order, colors, out_pdf, columns=4, title=None, caption=None,
         ylabel="log2(1 + CP10k), pseudobulk"):
    apply_rcparams()
    n = len(panels)
    ncol = min(columns, n)
    nrow = int(np.ceil(n / ncol))
    cap_lines = textwrap.wrap(caption, 150) if caption else []
    fig_h = 2.1 * nrow + 0.5 + 0.13 * len(cap_lines)
    fig = plt.figure(figsize=(1.9 * ncol + 0.6, fig_h))
    top = 1 - 0.45 / fig_h
    bottom = (0.35 + 0.13 * len(cap_lines)) / fig_h
    gs = fig.add_gridspec(nrow, ncol, left=0.09, right=0.98, top=top, bottom=bottom, hspace=0.75, wspace=0.55)
    xs = {g: i for i, g in enumerate(group_order)}
    for k, p in enumerate(panels.itertuples()):
        ax = fig.add_subplot(gs[k // ncol, k % ncol])
        v = values[values.panel == p.Index]
        finite = v.value[np.isfinite(v.value)]
        lo, hi = (float(finite.min()), float(finite.max())) if len(finite) else (0.0, 1.0)
        pad = max(0.15 * (hi - lo), 0.2)
        ymin, ymax = max(0.0, lo - pad), hi + 2.2 * pad
        offsets, facecolors = [], []
        for g in group_order:
            gv = v[(v.group == g) & np.isfinite(v.value)].sort_values("sample")
            m = len(gv)
            if m:
                dx = np.linspace(-0.18, 0.18, m) if m > 1 else np.zeros(1)
                offsets.extend(zip(xs[g] + dx, gv.value))
                facecolors.extend([to_rgba(colors[g], 0.85)] * m)
                mean = float(gv.value.mean())
                ax.plot([xs[g] - 0.28, xs[g] + 0.28], [mean, mean], color="#000000", linewidth=0.9, clip_on=False)
            mark = getattr(p, f"mark_{g}", None) if f"mark_{g}" in panels.columns else None
            if isinstance(mark, str) and mark:
                ax.text(xs[g], ymax - 0.02 * (ymax - ymin), mark, ha="center", va="top", fontsize=6, clip_on=False)
        if offsets:
            ec = EllipseCollection(widths=np.full(len(offsets), 4.5), heights=np.full(len(offsets), 4.5),
                                   angles=np.zeros(len(offsets)), units="points", offsets=np.asarray(offsets),
                                   offset_transform=ax.transData, facecolors=facecolors, edgecolors="#000000",
                                   linewidths=0.3)
            ax.add_collection(ec)
        counts = v[np.isfinite(v.value)].groupby("group").size()
        ax.set_xticks(range(len(group_order)))
        ax.set_xticklabels([f"{g}\nn={int(counts.get(g, 0))}" for g in group_order], fontsize=5.5)
        ax.set_xlim(-0.6, len(group_order) - 0.4)
        ax.set_ylim(ymin, ymax)
        ax.tick_params(axis="y", labelsize=6, length=2)
        ax.tick_params(axis="x", length=0)
        sub = getattr(p, "subtitle", "") if "subtitle" in panels.columns else ""
        ax.set_title(f"{p.gene} | {p.cell_state}", fontsize=7, pad=9 if sub else 3)
        if sub:
            ax.text(0.5, 1.01, textwrap.shorten(sub, 48), transform=ax.transAxes, ha="center", va="bottom",
                    fontsize=4.8, color="#4D4D4D", clip_on=False)
        if k % ncol == 0:
            ax.set_ylabel(ylabel, fontsize=6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.patch.set_visible(False)
        for artist in ax.get_children():
            artist.set_clip_on(False)
    if title:
        fig.text(0.5, 1 - 0.18 / fig_h, title, ha="center", va="center", fontsize=8)
    for i, line in enumerate(cap_lines):
        fig.text(0.02, (0.12 + 0.13 * (len(cap_lines) - 1 - i)) / fig_h, line, ha="left", va="bottom", fontsize=5.5)
    fig.patch.set_visible(False)
    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)
    return out_pdf


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pseudobulk-h5ad", required=True)
    ap.add_argument("--cluster-col", required=True)
    ap.add_argument("--sample-col", required=True)
    ap.add_argument("--groups-tsv", required=True, help="columns sample, group; only listed samples are drawn")
    ap.add_argument("--group-order", required=True, help="comma-separated groups, one column each, left to right")
    ap.add_argument("--group-colors", default=None, help="comma-separated #RRGGBB, one per group")
    ap.add_argument("--panels-tsv", required=True,
                    help="columns gene, cell_state; optional subtitle and mark_<group> text columns")
    ap.add_argument("--columns", type=int, default=4)
    ap.add_argument("--title", default=None)
    ap.add_argument("--caption-file", default=None, help="figure legend text printed under the panels")
    ap.add_argument("--out", required=True, help="output PDF; the plotted values go to <out>.values.tsv")
    a = ap.parse_args(argv)
    groups = pd.read_csv(a.groups_tsv, sep="\t", dtype=str)
    order = [g.strip() for g in a.group_order.split(",") if g.strip()]
    unknown = sorted(set(groups.group) - set(order))
    if unknown:
        raise ValueError(f"groups not in --group-order: {unknown}")
    hexes = [c.strip() for c in a.group_colors.split(",")] if a.group_colors else DEFAULT_COLORS[:len(order)]
    if len(hexes) != len(order) or not all(c.startswith("#") and len(c) == 7 for c in hexes):
        raise ValueError("--group-colors needs one #RRGGBB per group")
    panels = pd.read_csv(a.panels_tsv, sep="\t", dtype=str).fillna("")
    values = load_values(a.pseudobulk_h5ad, a.cluster_col, a.sample_col, panels, groups)
    caption = Path(a.caption_file).read_text().strip() if a.caption_file else None
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    values.to_csv(str(out) + ".values.tsv", sep="\t", index=False)
    draw(values, panels, order, dict(zip(order, hexes)), out, columns=a.columns, title=a.title, caption=caption)
    print(f"[pseudobulk_dotplot] {len(panels)} panels, {int(values.value.notna().sum())} of {len(values)} "
          f"(panel, sample) pseudobulks present; wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
