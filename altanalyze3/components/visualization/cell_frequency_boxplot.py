#!/usr/bin/env python3
"""Cell-state frequency per donor, two groups compared, one box plot per capture.

The analysis and figure behind Figure 2d of the aged-HSC manuscript
(``cell_frequency_<capture>_aged_vs_young.pdf``, 2026-08-11). Two commands.

``compute``
    Input: one row per cell with a sample, a cell state and a capture column (for example the
    cellHarmony barcode table ``query_to_sample_annotated.tsv``), plus a sample-to-group map
    (tab file or Excel sheet with a Sample column and a Group column). For every sample and
    capture, the percentage of that sample's cells in that capture that fall in each cell state.
    The denominator is the sample's own cells inside the capture, so every donor contributes one
    number per state.
    Test: Mann-Whitney U, two-sided, group 1 against group 2, one test per state per capture.
    Frequencies are bounded percentages, many states are rare and right-skewed and several hold
    zeros, so a t-test's normality assumption is not supportable at n of 5 to 8 per group. A
    Shapiro-Wilk p on the pooled values is written per state so the choice is auditable.
    Effect sizes: difference and ratio of medians, rank-biserial correlation. Benjamini-Hochberg
    across the states tested within each capture. States with zero variance across every sample
    are counted and reported, not dropped silently. A state needs at least ``--min-group-n``
    samples in both groups to be tested.
    The same test runs on lineage sums when a lineage map is given (``--lineage-tsv``: CellState,
    Lineage). The default map is the module's own curation of the hematopoietic order; it is
    written out so it can be corrected.

``plot``
    One figure per capture. Only the states that differ are drawn: raw p below ``--alpha``. With
    ``--min-mean-pct`` a state must also reach that mean percentage in at least one group (the
    abundance filter of the 2026-08-17 variant); dropped states go to ``dropped_states_*.tsv``.
    Group 1 on the upper bar in blue with a triangle per sample, group 2 on the lower bar in red
    with a circle per sample. The box spans the first to the third quartile, the thin line spans
    the minimum to the maximum. The y label is coloured by the higher group, and a state at
    BH-FDR below ``--alpha`` is bold with an asterisk. Rows follow the hematopoietic order.

Outputs of ``compute`` in ``--outdir``: per_sample_frequencies.tsv, cellstate_frequency_tests.tsv,
lineage_frequency_tests.tsv, lineage_assignment.tsv, compute_parameters.json.
Outputs of ``plot``: cell_frequency_<capture>_<g2>_vs_<g1>[_min<N>pct].pdf, plotted_states_<capture>.tsv,
dropped_states_<capture>.tsv (with the filter), plot_parameters.json.

Validated 2026-10-23: ``compute`` on the 2026-07-08 barcode table reproduced the 2026-08-11
per_sample_frequencies.tsv and cellstate_frequency_tests.tsv (every value), and ``plot``
reproduced the plotted state sets of the three captures.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, shapiro

try:
    from altanalyze3.components.visualization.figure_style import (
        HEMATOPOIETIC_STATE_ORDER, apply_rcparams, order_states, pdf_vector_report, strip_masks,
        vector_circles, vector_triangles)
except ImportError:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from figure_style import (HEMATOPOIETIC_STATE_ORDER, apply_rcparams, order_states,
                              pdf_vector_report, strip_masks, vector_circles, vector_triangles)

C_GROUP1 = "#5BC8F0"   # young, upper bar
C_GROUP2 = "#EE3B24"   # aged, lower bar
C_POINT = "#000000"
C_WHISKER = "#808080"

# Lineage map: an assistant curation of the hematopoietic order (2026-08-11), not a lab-supplied
# file. ``compute`` writes it to lineage_assignment.tsv; pass --lineage-tsv to replace it.
DEFAULT_LINEAGE: Dict[str, List[str]] = {
    "HSC/MPP": ["HSC-1", "HSC-2", "MPP-1", "MPP-2", "MPP-MEP", "LMPP-1", "LMPP-1-cycling",
                "LMPP-2", "Multilin-1", "Multilin-2"],
    "MEP/Erythroid": ["MEP-1", "MEP-2", "MEP-Eryth-1", "MEP-Eryth-2", "ERP-1", "ERP-2", "ERP-3",
                      "ERP-4", "ERP-5", "ERP-6", "ERP-7", "ERP-8", "Erythroblast-1",
                      "Erythroblast-2", "Erythroblast-3"],
    "Megakaryocyte": ["MKP-early", "MKP-late", "MK-Platelet"],
    "Basophil/Mast": ["BMCP-1", "BMCP-2"],
    "Lymphoid progenitor": ["CLP", "Pro-B-Early", "Pro-B-1", "Pro-B-2", "Pro-B-3",
                            "Pro-B-cycling-1", "Pro-B-cycling-2", "Pro-B-Early-cycling",
                            "Transitional-B-1", "Transitional-B-2", "pre-B"],
    "Mature B/Plasma": ["B Memory-1", "B Memory-2", "Plasma Cell"],
    "T/NK": ["CD4 TCM", "T CD4 Activated", "T CD4 Naive-1", "T CD4 Naive-2", "T CD8 Naive",
             "CD8 TEM", "MAIT", "NK-Mature"],
    "Granulocyte/GMP": ["MultiLin-GMP-1", "MultiLin-GMP-2", "MultiLin-GMP-3", "preNeu",
                        "immNeu-1", "immNeu-2"],
    "Monocyte": ["cMOP", "Mono-1", "Mono-2", "Intermediate Mono-1", "Intermediate Mono-2",
                 "Intermediate Mono-3", "Classical-Mono", "Non-Classical Mono-1",
                 "Non-Classical Mono-2", "Myeloid intermediate 1", "Myeloid intermediate 2",
                 "Myeloid intermediate 3", "Mac"],
    "Dendritic": ["MDP-1", "MDP-2", "MDP-3", "MDP-4", "MDP-5", "MDP-6", "pre-DC-1", "pre-DC-2",
                  "pre-DC-3", "cDC1", "cDC2-1", "cDC2-2", "pDC", "ASDC"],
    "Stromal": ["MSC Fibroblast-1", "MSC Fibroblast-2", "Osteoblast", "Stromal Vascular"],
}


# ------------------------------------------------------------------ statistics

def benjamini_hochberg(p) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(o[::-1]):
        r = n - rank
        prev = min(prev, p[i] * n / r)
        q[i] = prev
    return q


def test_table(freq: pd.DataFrame, feature_col: str, capture: str, g1: str, g2: str,
               min_group_n: int) -> Tuple[pd.DataFrame, int]:
    """Mann-Whitney per feature, g2 against g1, on one capture's frequency table."""
    rows = []
    n_flat = 0
    for feat, sub in freq.groupby(feature_col, observed=True):
        a = sub.loc[sub["Grp"] == g2, "pct"].to_numpy(float)   # aged
        y = sub.loc[sub["Grp"] == g1, "pct"].to_numpy(float)   # young
        if len(a) < min_group_n or len(y) < min_group_n:
            continue
        pooled = np.concatenate([a, y])
        if np.allclose(pooled, pooled[0]):
            n_flat += 1
            continue
        u, p = mannwhitneyu(a, y, alternative="two-sided")
        rb = 2.0 * u / (len(a) * len(y)) - 1.0     # +1: g2 always higher
        sw = np.nan
        if len(pooled) >= 3 and np.ptp(pooled) > 0:
            sw = float(shapiro(pooled).pvalue)
        rows.append({"capture": capture, feature_col: feat,
                     f"n_{g2}": len(a), f"n_{g1}": len(y),
                     f"median_{g2}": np.median(a), f"median_{g1}": np.median(y),
                     f"mean_{g2}": a.mean(), f"mean_{g1}": y.mean(),
                     "median_diff": np.median(a) - np.median(y),
                     "ratio_medians": (np.median(a) + 1e-9) / (np.median(y) + 1e-9),
                     "rank_biserial": rb, "U": u, "pval": p, "shapiro_p_pooled": sw})
    T = pd.DataFrame(rows)
    if T.empty:
        return T, n_flat
    T["fdr"] = benjamini_hochberg(T["pval"].values)
    T["higher_in"] = np.where(T["median_diff"] > 0, g2, g1)
    T.loc[T["median_diff"] == 0, "higher_in"] = "tie"
    T = T.sort_values("pval").reset_index(drop=True)
    return T, n_flat


def read_group_map(path: str, sheet: Optional[str], sample_col: str, group_col: str) -> Dict[str, str]:
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        w = pd.read_excel(p, sheet_name=sheet or 0, dtype=str)
    else:
        w = pd.read_csv(p, sep="\t", dtype=str)
    for c in (sample_col, group_col):
        if c not in w.columns:
            raise ValueError(f"{path} lacks column {c!r}; columns are {list(w.columns)}")
    w = w.dropna(subset=[sample_col])
    dup = w[sample_col].duplicated()
    if dup.any():
        raise ValueError(f"duplicated samples in {path}: {sorted(w.loc[dup, sample_col])}")
    return dict(zip(w[sample_col], w[group_col]))


def read_lineage(path: Optional[str]) -> Dict[str, str]:
    if path:
        t = pd.read_csv(path, sep="\t", dtype=str)
        return dict(zip(t["CellState"], t["Lineage"]))
    return {s: L for L, ss in DEFAULT_LINEAGE.items() for s in ss}


def compute(args) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    g1, g2 = [g.strip() for g in args.groups.split(",")]
    usecols = [args.state_col, args.capture_col, args.sample_col]
    d = pd.read_csv(args.annotations, sep="\t", dtype=str, usecols=usecols)
    gmap = read_group_map(args.group_map, args.group_map_sheet, args.map_sample_col, args.map_group_col)
    d["Grp"] = d[args.sample_col].map(gmap)
    h = d[d["Grp"].isin([g1, g2])].copy()
    print(f"[load] {args.annotations}: {len(d)} cells; {len(h)} in groups {g1}/{g2}; "
          f"samples {h[args.sample_col].nunique()}; states {h[args.state_col].nunique()}")
    s2l = read_lineage(args.lineage_tsv)
    h["Lineage"] = h[args.state_col].map(s2l)
    unmapped = sorted(set(h.loc[h["Lineage"].isna(), args.state_col]))
    if unmapped:
        print(f"[WARN] {len(unmapped)} cell states not in the lineage map: {unmapped}")
        h["Lineage"] = h["Lineage"].fillna("Unassigned")
    pd.DataFrame([{"CellState": s, "Lineage": L} for s, L in sorted(s2l.items())]).to_csv(
        outdir / "lineage_assignment.tsv", sep="\t", index=False)

    captures = ([c.strip() for c in args.captures.split(",")] if args.captures
                else sorted(h[args.capture_col].dropna().unique()))
    all_state, all_lin, freqs, summary = [], [], [], []
    for cap in captures:
        c = h[h[args.capture_col] == cap]
        tot = c.groupby(["Grp", args.sample_col]).size().rename("n_total")
        n1 = c.loc[c["Grp"] == g1, args.sample_col].nunique()
        n2 = c.loc[c["Grp"] == g2, args.sample_col].nunique()
        print(f"\n===== capture {cap}: {len(c)} cells, {g2} {n2} samples, {g1} {n1} samples")
        for col, store in [(args.state_col, all_state), ("Lineage", all_lin)]:
            cnt = c.groupby(["Grp", args.sample_col, col], observed=True).size().rename("n").reset_index()
            grid = (cnt.set_index(["Grp", args.sample_col, col])["n"]
                    .unstack(fill_value=0).stack().rename("n").reset_index())
            grid = grid.merge(tot.reset_index(), on=["Grp", args.sample_col])
            grid["pct"] = 100.0 * grid["n"] / grid["n_total"]
            grid["capture"] = cap
            if col == args.state_col:
                freqs.append(grid.rename(columns={args.sample_col: "Sample", args.state_col: "CellState"}))
            T, flat = test_table(grid, col, cap, g1, g2, args.min_group_n)
            nsig = int((T["fdr"] < args.alpha).sum()) if not T.empty else 0
            nraw = int((T["pval"] < args.alpha).sum()) if not T.empty else 0
            print(f"  {col:10s}: tested {len(T)}  zero-variance skipped {flat}  "
                  f"raw p<{args.alpha} {nraw}  FDR<{args.alpha} {nsig}")
            summary.append({"capture": cap, "feature": col, "tested": len(T), "zero_variance": flat,
                            "raw_sig": nraw, "fdr_sig": nsig, f"n_{g1}": int(n1), f"n_{g2}": int(n2)})
            store.append(T.rename(columns={args.state_col: "CellState"}))

    S = pd.concat(all_state, ignore_index=True)
    L = pd.concat(all_lin, ignore_index=True)
    F = pd.concat(freqs, ignore_index=True)
    S.to_csv(outdir / "cellstate_frequency_tests.tsv", sep="\t", index=False)
    L.to_csv(outdir / "lineage_frequency_tests.tsv", sep="\t", index=False)
    F.to_csv(outdir / "per_sample_frequencies.tsv", sep="\t", index=False)
    for p in ("cellstate_frequency_tests.tsv", "lineage_frequency_tests.tsv",
              "per_sample_frequencies.tsv", "lineage_assignment.tsv"):
        print(f"[write] {outdir / p}")
    nn = S[S["shapiro_p_pooled"] < 0.05]
    print(f"[diagnostic] Shapiro-Wilk rejects normality for {len(nn)} of "
          f"{int(S['shapiro_p_pooled'].notna().sum())} tested state-capture combinations")
    params = {k: v for k, v in vars(args).items() if k != "func"}
    params.update({"summary": summary, "python": sys.version})
    with open(outdir / "compute_parameters.json", "w") as fh:
        json.dump(params, fh, indent=2)


# ------------------------------------------------------------------ plot

def draw(ax, freq: pd.DataFrame, states, tests: pd.DataFrame, g1: str, g2: str, alpha: float,
         fs_ytick: float, marker_pt: float):
    order = order_states(states)
    tmap = tests.set_index("CellState")
    bh = 0.17
    labels, colors, weights = [], [], []
    xmax = 0.0
    for i, st in enumerate(order):
        yc = len(order) - 1 - i
        for grp, col, draw_pts, dy in [(g1, C_GROUP1, vector_triangles, +bh),
                                       (g2, C_GROUP2, vector_circles, -bh)]:
            v = freq[(freq["CellState"] == st) & (freq["Grp"] == grp)]["pct"].to_numpy(float)
            if v.size == 0:
                continue
            q1, q3 = np.percentile(v, [25, 75])
            lo, hi = v.min(), v.max()
            xmax = max(xmax, hi)
            if hi > lo:   # a zero-length line is an empty path in the PDF
                ax.plot([lo, hi], [yc + dy, yc + dy], color=C_WHISKER, linewidth=0.7,
                        solid_capstyle="butt", zorder=1)
            if q3 - q1 > 0:
                ax.add_patch(_rect((q1, yc + dy - bh * 0.78), q3 - q1, bh * 1.56, col))
            draw_pts(ax, v, np.full(v.size, yc + dy), marker_pt, C_POINT, zorder=3)
        r = tmap.loc[st]
        sig = bool(r["fdr"] < alpha)
        labels.append(f"{st} *" if sig else st)
        colors.append(C_GROUP2 if r["higher_in"] == g2 else C_GROUP1)
        weights.append("bold" if sig else "normal")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(labels[::-1], fontsize=fs_ytick)
    for t, c, w in zip(ax.get_yticklabels(), colors[::-1], weights[::-1]):
        t.set_color(c)
        t.set_fontweight(w)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlim(-xmax * 0.02, xmax * 1.06)
    ax.set_xlabel("percent of the sample's cells in the capture")
    ax.tick_params(axis="y", length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    return order


def _rect(xy, w, h, color):
    from matplotlib.patches import Rectangle
    return Rectangle(xy, w, h, facecolor=color, edgecolor="none", zorder=2)


def plot(args) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    src = Path(args.tables_dir)
    outdir = Path(args.outdir) if args.outdir else src
    outdir.mkdir(parents=True, exist_ok=True)
    g1, g2 = [g.strip() for g in args.groups.split(",")]
    lab1, lab2 = ([g.strip() for g in args.group_labels.split(",")] if args.group_labels
                  else (g1.capitalize(), g2.capitalize()))
    F = pd.read_csv(src / "per_sample_frequencies.tsv", sep="\t")
    S = pd.read_csv(src / "cellstate_frequency_tests.tsv", sep="\t")
    print(f"[load] {src / 'per_sample_frequencies.tsv'} rows {len(F)}")
    print(f"[load] {src / 'cellstate_frequency_tests.tsv'} rows {len(S)}")
    captures = ([c.strip() for c in args.captures.split(",")] if args.captures
                else list(dict.fromkeys(S["capture"])))
    apply_rcparams(args.font_size)
    written, plotted_summary = [], []
    suffix = f"_min{args.min_mean_pct:g}pct" if args.min_mean_pct is not None else ""
    for cap in captures:
        sig = S[(S["capture"] == cap) & (S["pval"] < args.alpha)].copy()
        n_tested = int((S["capture"] == cap).sum())
        if args.min_mean_pct is not None:
            keep = (sig[f"mean_{g2}"] >= args.min_mean_pct) | (sig[f"mean_{g1}"] >= args.min_mean_pct)
            t, dropped = sig[keep].copy(), sig[~keep].copy()
            dropped.sort_values("pval").to_csv(outdir / f"dropped_states_{cap}.tsv", sep="\t", index=False)
        else:
            t, dropped = sig, sig.iloc[0:0]
        f = F[F["capture"] == cap]
        na = int(t[f"n_{g2}"].iloc[0]) if len(t) else 0
        ny = int(t[f"n_{g1}"].iloc[0]) if len(t) else 0
        print(f"\n[{cap}] tested {n_tested} states; raw p<{args.alpha}: {len(sig)} of {n_tested}; "
              f"drawn {len(t)} of {len(sig)}; dropped {len(dropped)}; "
              f"FDR<{args.alpha} among drawn {int((t['fdr'] < args.alpha).sum())}; "
              f"{lab2} n={na}, {lab1} n={ny}")
        plotted_summary.append({"capture": cap, "tested": n_tested, "raw_sig": len(sig),
                                "drawn": len(t), "dropped": len(dropped)})
        if t.empty:
            print(f"  [{cap}] no state to draw")
            continue
        h = max(2.4, args.row_height * len(t) + 1.5)
        fig, ax = plt.subplots(figsize=(args.fig_width, h))
        order = draw(ax, f, set(t["CellState"]), t, g1, g2, args.alpha, args.font_size, args.marker_pt)
        handles = [Rectangle((0, 0), 1, 1, facecolor=C_GROUP1, edgecolor="none", label=f"{lab1} (n={ny}, triangles)"),
                   Rectangle((0, 0), 1, 1, facecolor=C_GROUP2, edgecolor="none", label=f"{lab2} (n={na}, circles)")]
        ax.legend(handles=handles, frameon=False, loc="lower right")
        filt = (f"; mean >= {args.min_mean_pct:g}% in one group" if args.min_mean_pct is not None else "")
        ax.set_title(f"{cap} capture: cell-state frequency per sample, {lab2} (n={na}) vs {lab1} (n={ny})\n"
                     f"{len(t)} of {n_tested} states at Mann-Whitney two-sided p<{args.alpha}{filt}; "
                     f"* and bold = BH-FDR<{args.alpha} within the capture; label colour = higher median group\n"
                     f"box = quartile 1 to 3, line = minimum to maximum, one mark per sample",
                     fontsize=args.font_size)
        fig.tight_layout()
        strip_masks(fig)
        p = outdir / f"cell_frequency_{cap}_{g2}_vs_{g1}{suffix}.pdf"
        fig.savefig(p, format="pdf")
        if args.svg:
            fig.savefig(p.with_suffix(".svg"), format="svg")
        plt.close(fig)
        written.append(p)
        t.assign(plot_order=[order.index(s) for s in t["CellState"]]).sort_values("plot_order").to_csv(
            outdir / f"plotted_states_{cap}.tsv", sep="\t", index=False)
    reports = {str(p): pdf_vector_report(p) for p in written}
    for p in written:
        print(f"[write] {p}")
    params = {k: v for k, v in vars(args).items() if k != "func"}
    params.update({"summary": plotted_summary, "pdf": reports, "python": sys.version})
    with open(outdir / "plot_parameters.json", "w") as fh:
        json.dump(params, fh, indent=2)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="command", required=True)
    a = sp.add_parser("compute", help="per-sample frequencies and Mann-Whitney tests")
    a.add_argument("--annotations", required=True, help="tab file, one row per cell")
    a.add_argument("--group-map", required=True, help="tab or Excel file mapping sample to group")
    a.add_argument("--group-map-sheet", default=None, help="Excel sheet name")
    a.add_argument("--map-sample-col", default="Sample")
    a.add_argument("--map-group-col", default="Group")
    a.add_argument("--sample-col", default="Sample")
    a.add_argument("--state-col", default="CellState")
    a.add_argument("--capture-col", default="SortGate")
    a.add_argument("--captures", default=None, help="comma list; default: every capture present")
    a.add_argument("--groups", default="young,aged", help="group 1 (blue, upper), group 2 (red, lower)")
    a.add_argument("--lineage-tsv", default=None, help="CellState, Lineage table replacing the default map")
    a.add_argument("--min-group-n", type=int, default=3)
    a.add_argument("--alpha", type=float, default=0.05)
    a.add_argument("--outdir", required=True)
    a.set_defaults(func=compute)

    b = sp.add_parser("plot", help="box plots from the compute tables")
    b.add_argument("--tables-dir", required=True, help="folder holding the compute outputs")
    b.add_argument("--outdir", default=None, help="default: --tables-dir")
    b.add_argument("--captures", default=None, help="comma list; default: every capture in the tests")
    b.add_argument("--groups", default="young,aged")
    b.add_argument("--group-labels", default="Young,Aged")
    b.add_argument("--alpha", type=float, default=0.05)
    b.add_argument("--min-mean-pct", type=float, default=None,
                   help="draw a state only when its mean reaches this percent in one group")
    b.add_argument("--fig-width", type=float, default=8.2)
    b.add_argument("--row-height", type=float, default=0.40)
    b.add_argument("--font-size", type=float, default=8)
    b.add_argument("--marker-pt", type=float, default=3.4, help="sample mark diameter in points")
    b.add_argument("--svg", action=argparse.BooleanOptionalAction, default=False)
    b.set_defaults(func=plot)
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
