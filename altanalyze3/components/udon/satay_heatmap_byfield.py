#!/usr/bin/env python3
"""SATAY-UDON metadata-field-grouped heatmaps (alternative layout).

One stacked block per metadata field (covariate): the field name is the title, rows are
the cell types enriched for that field, columns are the UDON clusters, color = enrichment
p-value (red = p<0.1, black = n.s.), matching the SATAY-UDON figure style. Large/tall is OK.

Test per (field, celltype, cluster): Fisher one-sided, background = all pseudobulks of that
cell type. Pure Python (pandas/scipy/matplotlib). No R. Editable PDF fonts (Arial)."""
import os, sys, numpy as np, pandas as pd, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.stats import fisher_exact
from satay_heatmap import build_features, add_age, load_run, load_final, auto_covars, unit_fisher, PB, OUT

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['pdf.fonttype'] = 42; plt.rcParams['ps.fonttype'] = 42

MIN_PB = 8        # min pseudobulks for a (cluster x celltype) test
P_THR = 0.1       # red if p < this (n.s. = black), matching the SATAY figure


def field_matrices(cl, feat, defined, covars):
    clusters = sorted(cl["cluster"].unique(), key=lambda x: (len(x), x))
    cts_all = [ct for ct, n in cl["celltype"].value_counts().items() if n >= MIN_PB]
    out = {}
    for (key, disp, grp) in covars:
        if key not in feat: continue
        smap, dfn = feat[key], defined[key]
        M = np.ones((len(cts_all), len(clusters)))
        for i, ct in enumerate(cts_all):
            T = cl[(cl["celltype"] == ct) & (cl["Sample"].isin(dfn))]
            if len(T) < MIN_PB: continue
            fpos = T["Sample"].map(smap).fillna(0).astype(int).values
            for j, c in enumerate(clusters):
                p, _o, a, _cn, _fn = unit_fisher(T["Sample"].values, (T["cluster"] == c).values, fpos)
                M[i, j] = p                          # donor-collapsed; p=1 if < MIN_UNITS distinct donors
        keep = [i for i in range(len(cts_all)) if (M[i] < P_THR).any()]
        if keep:
            out[disp] = (np.array(cts_all)[keep].tolist(), M[keep], clusters)
    return out


def plot_byfield(mats, title, outbase):
    fields = list(mats)
    if not fields:
        print(f"  {os.path.basename(outbase)}: no enriched fields, skipped"); return
    clusters = mats[fields[0]][2]; ncol = len(clusters)
    cmap = LinearSegmentedColormap.from_list("satay_red", ["#FFB3B3", "#FF4D4D", "#E60000", "#990000"])
    cmap.set_under("#000000")                                            # p>=0.1 -> black (n.s.)
    # stack blocks top-to-bottom in ONE axes with tight, fixed gaps (no GridSpec whitespace)
    TITLE_ROWS, GAP_ROWS = 1.0, 0.55
    placements = []; y = 0.0
    for f in fields:
        cts, M, _ = mats[f]
        y += TITLE_ROWS; placements.append((f, cts, M, y)); y += len(cts) + GAP_ROWS
    total = y
    fig_h = 0.20 * total + 1.6
    fig, ax = plt.subplots(figsize=(0.34 * ncol + 4.8, fig_h))
    for (f, cts, M, y0) in placements:
        nlp = -np.log10(np.clip(M, 1e-12, 1)); nlp[M >= P_THR] = 0.0
        X = np.arange(ncol + 1); Y = y0 + np.arange(len(cts) + 1)
        ax.pcolormesh(X, Y, nlp, cmap=cmap, vmin=1.0, vmax=4.0,           # VECTOR cells -> renders in Acrobat
                      edgecolors="white", linewidth=0.7)
        for i, ct in enumerate(cts):
            ax.text(-0.25, y0 + i + 0.5, ct, ha="right", va="center", fontsize=10)
        ax.text(ncol / 2.0, y0 - 0.20, f, ha="center", va="bottom", fontsize=13, fontweight="bold")
    # cluster labels + legend in DATA coords, right under the last block (no figure-margin gap)
    for j, c in enumerate(clusters):
        ax.text(j + 0.5, total + 0.4, str(c), ha="center", va="top", rotation=90, fontsize=10)
    ax.set_xlim(-0.5, ncol + 0.2); ax.set_ylim(total + 3.4, -0.7)        # invert y -> top-down
    ax.set_xticks([]); ax.set_yticks([]); ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(title, fontsize=15, pad=12)
    leg = [Patch(fc="#990000", label="<0.0001"), Patch(fc="#E60000", label="0.001"),
           Patch(fc="#FF4D4D", label="0.01"), Patch(fc="#FFB3B3", label="0.1"),
           Patch(fc="#000000", label="n.s.")]
    ax.legend(handles=leg, title="P-value", loc="upper center", bbox_to_anchor=(ncol / 2.0, total + 2.4),
              bbox_transform=ax.transData, ncol=5, fontsize=11, title_fontsize=11, frameon=False)
    dpi_png = max(80, min(200, int(63000 / max(fig_h, 1))))              # keep PNG under the 2^16 pixel cap
    fig.savefig(outbase + ".png", dpi=dpi_png, bbox_inches="tight")
    fig.savefig(outbase + ".pdf", bbox_inches="tight"); plt.close(fig)   # vector PDF (no image stream)
    print(f"  wrote {os.path.basename(outbase)}  ({len(fields)} fields, {int(total)} rows x {ncol} clusters, {fig_h:.0f}in@{dpi_png}dpi)")


def generate_all(feat=None, defined=None):
    if feat is None:
        feat, defined = build_features()
    if ("age>=60",) not in feat:
        add_age(feat, defined)
    jobs = [("matched_sex_platform", load_run(os.path.join(PB, "UDON/matched_sex_platform/udon_core/udon_clusters.txt"))),
            ("udon_core", load_run(os.path.join(PB, "UDON/udon_core/udon_clusters.txt"))),
            ("study_aware_final", load_final(os.path.join(PB, "UDON/study_aware/final_program_assignments.tsv")))]
    for d in sorted(glob.glob(os.path.join(PB, "UDON/study_aware/per_study/*/udon_clusters.txt"))):
        jobs.append((f"per_study_{os.path.basename(os.path.dirname(d))}", load_run(d)))
    COV = auto_covars(feat, defined)                                 # automatic: every feature, no curation
    for name, cl in jobs:
        try:
            mats = field_matrices(cl, feat, defined, COV)
            plot_byfield(mats, f"SATAY-UDON by metadata field: {name}", os.path.join(OUT, f"satay_byfield_{name}"))
        except Exception as e:
            print(f"  {name}: FAILED {e}")


if __name__ == "__main__":
    generate_all()
