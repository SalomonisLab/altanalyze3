#!/usr/bin/env python3
"""UMAP of one cell state coloured by donor, with the donors shaded by group (young blue, aged red).

The figure behind Figure 2c of the aged-HSC manuscript. The embedding itself comes from
AltAnalyze 2.1.4 (``AltAnalyze.py --image UMAP``: genes z-scored, ``umap.UMAP(n_neighbors=50,
min_dist=0.75, metric='correlation')``, umap-learn 0.3.10, Python 2.7), run on the MarkerFinder
genes of the k=2 subcluster split (source: https://github.com/nsalomonis/altanalyze,
``visualization_scripts/clustering.py``). This module does not compute an embedding. It has two commands.

``prepare``
    Writes the AltAnalyze2 UMAP input for one cell state: the expression matrix restricted to the
    MarkerFinder genes of the split, and the groups file that AltAnalyze2 reads beside it. It prints
    the exact AltAnalyze2 command to run.

``plot``
    Reads an AltAnalyze2 ``*-UMAP_scores.txt`` (uid, x, y; no header) and a groups table (uid,
    donor, group) and draws every cell as a vector circle coloured by donor. Donors of the first
    group take blue shades, donors of the second group take red shades, light to dark in the donor
    order given (default: sorted donor names within each group). With ``--subclusters`` it adds a
    second panel coloured by the k=2 subcluster and writes the subcluster x group cell counts.

Outputs of ``plot`` beside ``--out``:
  <stem>.pdf                       the UMAP panel(s)
  <stem>_cell_counts.tsv           cells per donor and group (and per subcluster when given)
  <stem>_colors.tsv                the hex colour of every donor
  <stem>_parameters.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, to_hex

try:
    from altanalyze3.components.visualization.figure_style import (
        apply_rcparams, pdf_vector_report, strip_masks, vector_circles)
except ImportError:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from figure_style import apply_rcparams, pdf_vector_report, strip_masks, vector_circles

# Shade ramps, light to dark. Blue for the first group (young), red for the second (aged).
BLUE_RAMP = ["#B7C9F2", "#6E8FE0", "#2F57C9", "#1030A0", "#0A1F6E"]
RED_RAMP = ["#F6C3C3", "#EE8A8A", "#E24B4B", "#C01E1E", "#8A0B0B"]
SUBCLUSTER_COLORS = ["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02"]


def shade_ramp(anchors: Sequence[str], n: int) -> List[str]:
    """n hex colours sampled light to dark along the anchor ramp."""
    cmap = LinearSegmentedColormap.from_list("ramp", anchors)
    if n == 1:
        return [to_hex(cmap(0.75))]
    return [to_hex(cmap(v)) for v in np.linspace(0.0, 1.0, n)]


def short_donor(name: str) -> str:
    """'WM72-N-1130' -> 'WM72'; names without '-' stay as they are."""
    return name.split("-")[0]


# ------------------------------------------------------------------ prepare

def prepare(args) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    groups = _read_groups(Path(args.groups), args.uid_col, args.donor_col, args.group_col)
    markers = pd.read_csv(args.marker_genes, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    if markers and markers[0].lower() == "uid":
        markers = markers[1:]
    marker_set = set(markers)
    kept, n_rows, header = [], 0, None
    for chunk in pd.read_csv(args.expression, sep="\t", index_col=0, chunksize=2000, dtype=str):
        n_rows += len(chunk)
        header = chunk.columns if header is None else header
        kept.append(chunk[chunk.index.isin(marker_set)])
    expr = pd.concat(kept)
    missing = [m for m in markers if m not in expr.index]
    if missing:
        raise ValueError(f"{len(missing)} of {len(markers)} marker genes absent from {args.expression}")
    expr = expr.loc[markers]
    exp_path = outdir / f"exp.{args.name}.txt"
    expr.to_csv(exp_path, sep="\t", index_label=expr.index.name or "UID")

    cells = pd.DataFrame({"uid": list(expr.columns)})
    g = cells.merge(groups, on="uid", how="left")
    n_missing = int(g["donor"].isna().sum())
    if n_missing:
        raise ValueError(f"{n_missing} of {len(cells)} cells have no groups row")
    label = g["donor"].map(short_donor) if args.short_donor else g["donor"]
    grp_path = outdir / f"groups.{args.name}.txt"
    pd.DataFrame({"uid": g["uid"], "a": label, "b": label}).to_csv(
        grp_path, sep="\t", index=False, header=False)
    print(f"[prepare] {args.expression}: {n_rows} genes x {expr.shape[1]} cells; "
          f"{len(markers)} MarkerFinder genes kept")
    print(f"[write] {exp_path} ({expr.shape[0]} genes x {expr.shape[1]} cells)")
    print(f"[write] {grp_path} ({len(g)} cells, {label.nunique()} donor labels)")
    print("[altanalyze2] run, from the AltAnalyze 2.1.4 checkout, with its Python 2.7:")
    print(f"  python2.7 AltAnalyze.py --input {exp_path} --image UMAP --species Hs --platform RNASeq --display False")
    print(f"  -> {outdir / 'DataPlots' / f'exp.{args.name}-UMAP_scores.txt'}")


# ------------------------------------------------------------------ plot

def _read_groups(path: Path, uid_col: str, donor_col: str, group_col: str) -> pd.DataFrame:
    g = pd.read_csv(path, sep="\t", dtype=str)
    for c in (uid_col, donor_col, group_col):
        if c not in g.columns:
            raise ValueError(f"{path} lacks column {c!r}; columns are {list(g.columns)}")
    g = g[[uid_col, donor_col, group_col]].rename(
        columns={uid_col: "uid", donor_col: "donor", group_col: "group"})
    dup = int(g["uid"].duplicated().sum())
    if dup:
        raise ValueError(f"{dup} duplicated uids in {path}")
    return g


def read_umap_scores(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, dtype={0: str})
    if df.shape[1] < 3:
        raise ValueError(f"{path}: expected uid, x, y columns, found {df.shape[1]}")
    df = df.iloc[:, :3]
    df.columns = ["uid", "x", "y"]
    if df["uid"].iloc[0].lower() in ("uid", "cell"):
        df = df.iloc[1:]
    df["x"] = df["x"].astype(float)
    df["y"] = df["y"].astype(float)
    if int(df["uid"].duplicated().sum()):
        raise ValueError(f"duplicated uids in {path}")
    return df.reset_index(drop=True)


def donor_colours(cells: pd.DataFrame, group_order: Sequence[str],
                  donor_order: Optional[Sequence[str]]) -> Dict[str, str]:
    """Hex colour per donor label: group 1 blue ramp, group 2 red ramp, light to dark."""
    ramps = [BLUE_RAMP, RED_RAMP]
    if len(group_order) != 2:
        raise ValueError("exactly two groups are supported (blue and red ramps)")
    donors_by_group: Dict[str, List[str]] = {}
    for grp in group_order:
        labels = sorted(cells.loc[cells["group"] == grp, "label"].unique())
        if donor_order:
            given = [d for d in donor_order if d in labels]
            extra = [d for d in labels if d not in given]
            labels = given + extra
        donors_by_group[grp] = labels
    out: Dict[str, str] = {}
    for grp, ramp in zip(group_order, ramps):
        labels = donors_by_group[grp]
        for lab, col in zip(labels, shade_ramp(ramp, len(labels))):
            out[lab] = col
    return out


def plot(args) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    stem = out.with_suffix("")
    coords = read_umap_scores(Path(args.coordinates))
    groups = _read_groups(Path(args.groups), args.uid_col, args.donor_col, args.group_col)
    cells = coords.merge(groups, on="uid", how="left")
    n_missing = int(cells["donor"].isna().sum())
    if n_missing:
        raise ValueError(f"{n_missing} of {len(cells)} cells in {args.coordinates} have no groups row")
    group_order = [g.strip() for g in args.group_order.split(",")]
    group_labels = ([g.strip() for g in args.group_labels.split(",")] if args.group_labels
                    else group_order)
    unknown = sorted(set(cells["group"]) - set(group_order))
    if unknown:
        raise ValueError(f"groups {unknown} in the groups file are not in --group-order {group_order}")
    cells["label"] = cells["donor"].map(short_donor) if args.short_donor else cells["donor"]
    donor_order = _read_list(args.donor_order)
    colours = donor_colours(cells, group_order, donor_order)
    cells["color"] = cells["label"].map(colours)

    counts = (cells.groupby(["group", "label"]).size().rename("n_cells").reset_index()
              .rename(columns={"label": "donor"}))
    frames = [counts]
    sub = None
    if args.subclusters:
        sub = pd.read_csv(args.subclusters, sep="\t", dtype=str)
        if "V1" in sub.columns:
            sub = sub[["uid", "V1"]].rename(columns={"V1": "subcluster"})
        else:
            sub = pd.read_csv(args.subclusters, sep="\t", header=None, dtype=str,
                              names=["uid", "V0", "subcluster"])[["uid", "subcluster"]]
        cells = cells.merge(sub, on="uid", how="left")
        n_nosub = int(cells["subcluster"].isna().sum())
        if n_nosub:
            raise ValueError(f"{n_nosub} of {len(cells)} cells have no subcluster label")
        sc = (cells.groupby(["subcluster", "group"]).size().rename("n_cells").reset_index())
        sc.insert(1, "donor", "all")
        frames.append(sc)
        scd = (cells.groupby(["subcluster", "group", "label"]).size().rename("n_cells")
               .reset_index().rename(columns={"label": "donor"}))
        frames.append(scd)
    counts_all = pd.concat(frames, ignore_index=True)
    if "subcluster" in counts_all.columns:
        counts_all["subcluster"] = counts_all["subcluster"].fillna("all")
    counts_all = counts_all[[c for c in ("subcluster", "group", "donor", "n_cells") if c in counts_all.columns]]
    counts_path = Path(f"{stem}_cell_counts.tsv")
    counts_all.to_csv(counts_path, sep="\t", index=False)
    pd.DataFrame({"donor": list(colours), "color": list(colours.values()),
                  "group": [cells.loc[cells["label"] == d, "group"].iloc[0] for d in colours]}
                 ).to_csv(f"{stem}_colors.tsv", sep="\t", index=False)

    rng = np.random.default_rng(args.seed)
    draw_order = rng.permutation(len(cells))   # no group hides another
    n_by_group = cells["group"].value_counts().to_dict()
    n_text = ", ".join(f"{lab} n={n_by_group.get(g, 0)}" for g, lab in zip(group_order, group_labels))
    title = args.title or Path(args.coordinates).stem.replace("-UMAP_scores", "").replace("exp.", "")

    apply_rcparams(7)
    n_panels = 2 if sub is not None else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(args.panel_size * n_panels + 1.4, args.panel_size))
    axes = np.atleast_1d(axes)
    pad_x = 0.06 * (cells["x"].max() - cells["x"].min())
    pad_y = 0.06 * (cells["y"].max() - cells["y"].min())
    xlim = (cells["x"].min() - pad_x, cells["x"].max() + pad_x)
    ylim = (cells["y"].min() - pad_y, cells["y"].max() + pad_y)

    ax = axes[0]
    vector_circles(ax, cells["x"].to_numpy()[draw_order], cells["y"].to_numpy()[draw_order],
                   args.point_size, cells["color"].to_numpy()[draw_order], alpha=1.0)
    ax.set_title(f"{title}: {len(cells)} cells, coloured by donor\n{n_text}", fontsize=7)
    handles = []
    for grp, lab in zip(group_order, group_labels):
        for d in [d for d in colours if cells.loc[cells['label'] == d, 'group'].iloc[0] == grp]:
            n = int((cells["label"] == d).sum())
            handles.append(Line2D([0], [0], marker="o", linestyle="none", markersize=4,
                                  markerfacecolor=colours[d], markeredgecolor="none",
                                  label=f"{d} ({lab}, n={n})"))
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              fontsize=6, handletextpad=0.4, borderaxespad=0.0)
    if args.label_groups:
        for grp, lab in zip(group_order, group_labels):
            m = cells["group"] == grp
            if m.any():
                ax.text(cells.loc[m, "x"].median(), cells.loc[m, "y"].median(), lab,
                        ha="center", va="center", fontsize=7, color="#000000", zorder=5)

    if sub is not None:
        ax2 = axes[1]
        subs = sorted(cells["subcluster"].unique())
        sub_col = {s: SUBCLUSTER_COLORS[i % len(SUBCLUSTER_COLORS)] for i, s in enumerate(subs)}
        vector_circles(ax2, cells["x"].to_numpy()[draw_order], cells["y"].to_numpy()[draw_order],
                       args.point_size, cells["subcluster"].map(sub_col).to_numpy()[draw_order])
        parts = []
        for s in subs:
            m = cells["subcluster"] == s
            parts.append(f"subcluster {s}: " + ", ".join(
                f"{lab} {int((m & (cells['group'] == g)).sum())}" for g, lab in zip(group_order, group_labels)))
        ax2.set_title(f"{title}: AltAnalyze2 k=2 subclusters\n" + "\n".join(parts), fontsize=7)
        ax2.legend(handles=[Line2D([0], [0], marker="o", linestyle="none", markersize=4,
                                   markerfacecolor=sub_col[s], markeredgecolor="none",
                                   label=f"subcluster {s} (n={int((cells['subcluster'] == s).sum())})")
                            for s in subs],
                   loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=6)

    for a in axes:
        a.set_xlim(xlim)
        a.set_ylim(ylim)
        a.set_xticks([])
        a.set_yticks([])
        a.set_aspect(args.aspect)
        a.set_xlabel("UMAP-X")
        a.set_ylabel("UMAP-Y")
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout()
    strip_masks(fig)
    fig.savefig(out, format="pdf")
    plt.close(fig)
    report = pdf_vector_report(out)
    params = {k: v for k, v in vars(args).items() if k != "func"}
    params.update({"n_cells": int(len(cells)), "n_by_group": {str(k): int(v) for k, v in n_by_group.items()},
                   "donor_colors": colours, "pdf": report, "python": sys.version})
    with open(f"{stem}_parameters.json", "w") as fh:
        json.dump(params, fh, indent=2)
    print(f"[write] {out}")
    print(f"[write] {counts_path}")
    print(f"[counts] {n_text}")
    if sub is not None:
        print(counts_all[counts_all["donor"] == "all"].to_string(index=False))


def _read_list(value: Optional[str]) -> Optional[List[str]]:
    if not value:
        return None
    p = Path(value)
    if p.exists():
        with open(p) as fh:
            return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    return [v.strip() for v in value.split(",") if v.strip()]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="command", required=True)

    a = sp.add_parser("prepare", help="write the AltAnalyze2 UMAP input for one cell state")
    a.add_argument("--expression", required=True, help="NMF-SVM/ExpressionInput/exp.NMF-MarkerFinder.txt")
    a.add_argument("--marker-genes", required=True, help="NMF-SVM/centroids/exp.MF.txt")
    a.add_argument("--groups", required=True, help="tab file with uid, donor and group columns")
    a.add_argument("--name", required=True, help="cell state name used in the output file names")
    a.add_argument("--outdir", required=True)
    a.add_argument("--uid-col", default="uid")
    a.add_argument("--donor-col", default="donor")
    a.add_argument("--group-col", default="age-group")
    a.add_argument("--short-donor", action=argparse.BooleanOptionalAction, default=True,
                   help="label donors by the token before the first '-'")
    a.set_defaults(func=prepare)

    b = sp.add_parser("plot", help="draw the donor-coloured UMAP from AltAnalyze2 coordinates")
    b.add_argument("--coordinates", required=True, help="AltAnalyze2 *-UMAP_scores.txt (uid, x, y)")
    b.add_argument("--groups", required=True, help="tab file with uid, donor and group columns")
    b.add_argument("--out", required=True, help="output PDF path")
    b.add_argument("--uid-col", default="uid")
    b.add_argument("--donor-col", default="donor")
    b.add_argument("--group-col", default="age-group")
    b.add_argument("--group-order", default="young,old", help="two group values: blue, then red")
    b.add_argument("--group-labels", default="Young,Aged", help="display names for --group-order")
    b.add_argument("--donor-order", default=None,
                   help="comma list or file: donor labels light to dark within each group")
    b.add_argument("--short-donor", action=argparse.BooleanOptionalAction, default=True)
    b.add_argument("--subclusters", default=None,
                   help="round1Metadata*.txt or FinalGroups-CellTypesFull.txt: adds a subcluster panel and counts")
    b.add_argument("--title", default=None)
    b.add_argument("--label-groups", action=argparse.BooleanOptionalAction, default=False,
                   help="write each group's display name at the median of its cells")
    b.add_argument("--aspect", choices=["equal", "auto"], default="equal",
                   help="equal keeps the embedding geometry; auto fills the square panel")
    b.add_argument("--point-size", type=float, default=2.2, help="circle diameter in points")
    b.add_argument("--panel-size", type=float, default=2.6, help="panel side in inches")
    b.add_argument("--seed", type=int, default=0, help="seed of the drawing-order shuffle")
    b.set_defaults(func=plot)
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
