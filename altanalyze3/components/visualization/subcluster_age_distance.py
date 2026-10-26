#!/usr/bin/env python3
"""Age-induced distance between the k=2 subclusters of every cell state, with a dot plot.

The analysis behind Figure 2b of the aged-HSC manuscript. AltAnalyze 2.1.4 (ICGS2 sub-clustering,
``--k 2``) split every cell state of the young + aged 10x bone-marrow atlas into two subclusters.
For each cell state this module answers one question: did the split follow donor age?

Per cell state (one directory under ``--subclusters-dir``):

1. Subcluster membership comes from ``ICGS-NMF/FinalGroups-CellTypesFull.txt`` (columns uid, V0,
   V1, no header) when it exists, else from ``NMF-SVM/SubtypeAnalyses/round1Metadata*.txt``
   (header uid, V0, V1). The subcluster of a cell is the value of its ``V1`` column.
2. Cells join the groups table (``--groups``: uid, donor, age group) on uid.
3. A subcluster is tested when it holds at least ``--min-donors`` donors and no donor contributes
   more than ``--max-donor-fraction`` of its cells.
4. Cell level: a 2x2 table of young/aged cells inside the subcluster against the rest of the cell
   state. Fisher's exact test when any cell of the table is below 5, else the chi-square test with
   Yates correction (scipy ``chi2_contingency`` default).
5. Donor level: a donor counts when it contributes more than ``--donor-min-fraction`` of the
   subcluster's cells. The 2x2 table is young/aged donors present in the subcluster against the
   young/aged donors absent from it, over all donors of the cell state. Same Fisher/chi-square rule.
6. Distance: the mean expression centroid of each subcluster over the MarkerFinder genes of the
   split (``NMF-SVM/centroids/exp.MF.txt``, log2 values from
   ``NMF-SVM/ExpressionInput/exp.NMF-MarkerFinder.txt``), Pearson r between the centroids, and
   ``distance = 1 - mean off-diagonal r``. For k=2 that is one correlation.
7. The cell state's p value is the smallest donor-level p over its tested subclusters.

The dot plot draws one point per cell state at ``distance``; colour and size encode
``-log10(min donor-level p)``.

Validated 2026-10-23 against the manuscript tables made with ``age_test3.py`` (2024-09-25) on
``.../AgedHealth-combined/ICGS2-SubClustering/k-2/SubClusters``: every count, p value, rho,
distance and -log10 p of all 87 cell states reproduced. See the validation log in the Figure 2
notebook folder.

Outputs, in ``--outdir``:
  subcluster_age_enrichment.tsv   one row per tested subcluster
  age_induced_distance.tsv        one row per cell state (the dot-plot table)
  age_induced_distance_dotplot.pdf
  parameters.json, run_log.txt
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

try:
    from altanalyze3.components.visualization.figure_style import (
        apply_rcparams, order_states, pdf_vector_report, strip_masks, vector_circles)
except ImportError:  # pragma: no cover - direct script execution
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from figure_style import (apply_rcparams, order_states, pdf_vector_report, strip_masks,
                              vector_circles)


# ------------------------------------------------------------------ inputs

def find_annotation_file(cell_type_dir: Path) -> Tuple[Path, bool]:
    """Return (path, has_header). FinalGroups-CellTypesFull.txt has no header; round1Metadata has."""
    icgs = cell_type_dir / "ICGS-NMF" / "FinalGroups-CellTypesFull.txt"
    if icgs.exists():
        return icgs, False
    sub = cell_type_dir / "NMF-SVM" / "SubtypeAnalyses"
    if sub.is_dir():
        for name in sorted(os.listdir(sub)):
            if name.startswith("round1Metadata") and name.endswith(".txt"):
                return sub / name, True
    raise FileNotFoundError(f"no FinalGroups-CellTypesFull.txt or round1Metadata*.txt in {cell_type_dir}")


def read_annotation(cell_type_dir: Path) -> pd.DataFrame:
    path, has_header = find_annotation_file(cell_type_dir)
    if has_header:
        meta = pd.read_csv(path, sep="\t")
    else:
        meta = pd.read_csv(path, sep="\t", header=None, names=["uid", "V0", "V1"])
    if "uid" not in meta.columns or "V1" not in meta.columns:
        raise ValueError(f"{path} lacks uid/V1 columns: {list(meta.columns)}")
    meta = meta[["uid", "V1"]].copy()
    meta["subcluster"] = meta["V1"].astype(str)
    return meta.drop(columns="V1")


def read_groups(path: Path, uid_col: str, donor_col: str, group_col: str) -> pd.DataFrame:
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


def read_marker_expression(cell_type_dir: Path, log) -> pd.DataFrame:
    """MarkerFinder genes (rows) x cells (columns), log2 values, read in chunks."""
    expr_path = cell_type_dir / "NMF-SVM" / "ExpressionInput" / "exp.NMF-MarkerFinder.txt"
    mf_path = cell_type_dir / "NMF-SVM" / "centroids" / "exp.MF.txt"
    markers = pd.read_csv(mf_path, sep="\t", header=None).iloc[:, 0].astype(str).tolist()
    if markers and markers[0].lower() == "uid":
        markers = markers[1:]
    marker_set = set(markers)
    kept = []
    n_rows = 0
    for chunk in pd.read_csv(expr_path, sep="\t", index_col=0, chunksize=2000):
        n_rows += len(chunk)
        kept.append(chunk[chunk.index.isin(marker_set)])
    expr = pd.concat(kept)
    missing = [m for m in markers if m not in expr.index]
    if missing:
        raise ValueError(f"{len(missing)} of {len(markers)} MarkerFinder genes absent from "
                         f"{expr_path}: {missing[:5]}")
    expr = expr.loc[markers]
    log(f"  expression {expr_path.name}: {n_rows} genes x {expr.shape[1]} cells; "
        f"{len(markers)} MarkerFinder genes kept")
    return expr


# ------------------------------------------------------------------ statistics

def two_by_two_p(table: List[List[int]]) -> Tuple[float, str]:
    """Fisher's exact test when any cell is below 5, else chi-square with Yates correction."""
    if min(table[0]) < 5 or min(table[1]) < 5:
        return float(fisher_exact(table)[1]), "fisher_exact"
    return float(chi2_contingency(table)[1]), "chi2_yates"


def valid_subclusters(data: pd.DataFrame, min_donors: int, max_donor_fraction: float) -> List[str]:
    keep = []
    for sub, d in data.groupby("subcluster", sort=False):
        counts = d["donor"].value_counts()
        if len(counts) >= min_donors and (counts / counts.sum()).max() <= max_donor_fraction:
            keep.append(sub)
    return keep


def test_subclusters(data: pd.DataFrame, subclusters: Sequence[str], group_a: str, group_b: str,
                     donor_min_fraction: float) -> pd.DataFrame:
    total_a = int((data["group"] == group_a).sum())
    total_b = int((data["group"] == group_b).sum())
    donors_a = int(data.loc[data["group"] == group_a, "donor"].nunique())
    donors_b = int(data.loc[data["group"] == group_b, "donor"].nunique())
    rows = []
    for sub in subclusters:
        d = data[data["subcluster"] == sub]
        n_a = int((d["group"] == group_a).sum())
        n_b = int((d["group"] == group_b).sum())
        cell_table = [[n_a, total_a - n_a], [n_b, total_b - n_b]]
        cell_p, cell_test = two_by_two_p(cell_table)

        counts = d["donor"].value_counts()
        valid = counts[counts > donor_min_fraction * len(d)].index
        dd = d[d["donor"].isin(valid)]
        if dd.empty:
            continue
        per_donor = dd.groupby(["donor", "group"]).size().unstack(fill_value=0)
        da = int(per_donor[group_a].gt(0).sum()) if group_a in per_donor else 0
        db = int(per_donor[group_b].gt(0).sum()) if group_b in per_donor else 0
        donor_table = [[da, donors_a - da], [db, donors_b - db]]
        donor_p, donor_test = two_by_two_p(donor_table)
        rows.append({
            "subcluster": sub, "n_cells": int(len(d)),
            f"cells_{group_a}": n_a, f"cells_{group_b}": n_b,
            "cell_level_test": cell_test, "cell_level_p": cell_p,
            f"donors_{group_a}": da, f"donors_{group_b}": db,
            f"donors_total_{group_a}": donors_a, f"donors_total_{group_b}": donors_b,
            "donor_level_test": donor_test, "donor_level_p": donor_p,
        })
    return pd.DataFrame(rows)


def centroid_correlation(expr: pd.DataFrame, data: pd.DataFrame) -> Tuple[float, pd.DataFrame]:
    """Mean centroid per subcluster over all its cells; Pearson r matrix; mean off-diagonal r."""
    centroids = {}
    for sub in data["subcluster"].unique():
        cells = data.loc[data["subcluster"] == sub, "uid"]
        missing = [c for c in cells if c not in expr.columns]
        if missing:
            raise ValueError(f"{len(missing)} of {len(cells)} cells of subcluster {sub} are absent "
                             f"from the expression matrix (first: {missing[0]})")
        centroids[sub] = expr[cells].mean(axis=1)
    cent = pd.DataFrame(centroids)
    corr = cent.corr()
    np.fill_diagonal(corr.values, np.nan)
    return float(corr.mean().mean()), corr


# ------------------------------------------------------------------ driver

def analyse_cell_type(cell_type_dir: Path, groups: pd.DataFrame, group_a: str, group_b: str,
                      min_donors: int, max_donor_fraction: float, donor_min_fraction: float,
                      log) -> Tuple[pd.DataFrame, dict]:
    meta = read_annotation(cell_type_dir)
    data = groups.merge(meta, on="uid")
    if len(data) != len(meta):
        log(f"  WARNING {len(meta) - len(data)} of {len(meta)} cells have no groups row")
    subs = valid_subclusters(data, min_donors, max_donor_fraction)
    tests = test_subclusters(data, subs, group_a, group_b, donor_min_fraction)
    expr = read_marker_expression(cell_type_dir, log)
    rho, corr = centroid_correlation(expr, data)
    summary = {
        "n_cells": int(len(data)),
        "n_subclusters": int(data["subcluster"].nunique()),
        "n_subclusters_tested": int(len(tests)),
        "average_rho": rho,
        "min_donor_p": float(tests["donor_level_p"].min()) if len(tests) else np.nan,
        "n_cells_unmatched": int(len(meta) - len(data)),
    }
    summary["distance"] = 1.0 - rho
    summary["neg_log10_p"] = (-np.log10(summary["min_donor_p"])
                              if np.isfinite(summary["min_donor_p"]) else np.nan)
    summary["status"] = "ok" if len(tests) else "no subcluster passed the donor filter"
    return tests, summary


def plot_dotplot(table: pd.DataFrame, out_pdf: Path, group_a: str, group_b: str,
                 order: Sequence[str], hline: Optional[float], cmap_name: str,
                 size_scale: float, alpha: float, figsize: Tuple[float, float]) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.patches import Rectangle

    apply_rcparams(7)
    t = table.set_index("cell_type").loc[list(order)].reset_index()
    ok = t["neg_log10_p"].notna() & t["distance"].notna()
    x = np.arange(len(t), dtype=float)
    vmax = float(np.nanmax(t["neg_log10_p"])) if ok.any() else 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)
    colors = [cmap(norm(v)) for v in t.loc[ok, "neg_log10_p"]]
    # matplotlib scatter 's' is an area in points^2; the original used s = size_scale * -log10 p
    diam = 2.0 * np.sqrt(np.maximum(size_scale * t.loc[ok, "neg_log10_p"].to_numpy(float), 1.0) / np.pi)

    fig, ax = plt.subplots(figsize=figsize)
    vector_circles(ax, x[ok.to_numpy()], t.loc[ok, "distance"].to_numpy(float), diam, colors,
                   alpha=alpha, zorder=3)
    if hline is not None:
        ax.plot([-0.5, len(t) - 0.5], [hline, hline], color="#9E9E9E", linewidth=0.5,
                linestyle=(0, (1, 2)), zorder=1)
    ax.set_xlim(-0.5, len(t) - 0.5)
    ymax = float(np.nanmax(t["distance"])) if ok.any() else 1.0
    ax.set_ylim(min(0.0, float(np.nanmin(t["distance"])) - 0.05), ymax + 0.08)
    ax.set_xticks(x)
    ax.set_xticklabels(t["cell_type"], rotation=90, fontsize=6)
    ax.set_ylabel("Age-induced distance\n(1 - Pearson r between the two subcluster centroids, MarkerFinder genes)")
    n_ok = int(ok.sum())
    ax.set_title(f"k=2 subcluster split per cell state, {group_b} vs {group_a} donors: "
                 f"{n_ok} of {len(t)} cell states tested; "
                 f"colour and size = -log10 min donor-level p (Fisher exact / chi-square)")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    # Vector colour bar: fig.colorbar embeds a raster image, so the ramp is drawn as rectangles.
    fig.tight_layout(rect=(0.0, 0.0, 0.92, 1.0))
    cax = fig.add_axes([0.935, 0.35, 0.01, 0.45])
    n_steps = 32
    edges = np.linspace(0.0, vmax, n_steps + 1)
    for lo, hi in zip(edges[:-1], edges[1:]):
        cax.add_patch(Rectangle((0.0, lo), 1.0, hi - lo, facecolor=cmap(norm((lo + hi) / 2.0)),
                                edgecolor="none"))
    cax.set_xlim(0.0, 1.0)
    cax.set_ylim(0.0, vmax)
    cax.set_xticks([])
    cax.yaxis.tick_right()
    cax.set_yticks([0.0, vmax / 2.0, vmax])
    cax.set_yticklabels([f"{v:.1f}" for v in (0.0, vmax / 2.0, vmax)])
    cax.set_title(f"{group_b} vs {group_a}\n-log10 p\n(donor level)", fontsize=6, loc="left")
    for s in ("top", "bottom", "left"):
        cax.spines[s].set_visible(False)
    strip_masks(fig)
    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)
    return pdf_vector_report(out_pdf)


def run(args) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / ("plot_log.txt" if args.plot_only else "run_log.txt")
    log_fh = open(log_path, "w")

    def log(msg: str) -> None:
        print(msg)
        log_fh.write(msg + "\n")
        log_fh.flush()

    t0 = time.time()
    if args.plot_only:
        table_path = outdir / "age_induced_distance.tsv"
        table = pd.read_csv(table_path, sep="\t")
        log(f"[plot-only] re-rendering from {table_path} ({len(table)} rows, "
            f"written {time.strftime('%Y-%m-%d %H:%M', time.localtime(table_path.stat().st_mtime))})")
        plotted = table[table["status"] == "ok"]
        order = [c for c in order_states(table["cell_type"], _read_order(args.order))
                 if c in set(plotted["cell_type"])]
        pdf_path = outdir / "age_induced_distance_dotplot.pdf"
        plot_dotplot(plotted, pdf_path, args.group_a, args.group_b, order, args.hline, args.cmap,
                     args.size_scale, args.alpha, (args.fig_width, args.fig_height))
        log(f"[write] {pdf_path}")
        log_fh.close()
        return
    sub_root = Path(args.subclusters_dir)
    groups = read_groups(Path(args.groups), args.uid_col, args.donor_col, args.group_col)
    log(f"[groups] {args.groups}: {len(groups)} cells, {groups['donor'].nunique()} donors, "
        f"groups {groups['group'].value_counts().to_dict()}")
    for g in (args.group_a, args.group_b):
        if g not in set(groups["group"]):
            raise ValueError(f"group {g!r} absent from column {args.group_col!r}")

    cell_types = sorted(d.name for d in sub_root.iterdir() if d.is_dir() and not d.name.startswith("."))
    if args.cell_types:
        wanted = [c.strip() for c in args.cell_types.split(",") if c.strip()]
        absent = [c for c in wanted if c not in cell_types]
        if absent:
            raise ValueError(f"requested cell types without a folder: {absent}")
        cell_types = wanted
    log(f"[subclusters] {sub_root}: {len(cell_types)} cell-state folders")

    all_tests, summaries, failures = [], [], []
    for i, ct in enumerate(cell_types, 1):
        log(f"[{i}/{len(cell_types)}] {ct}")
        try:
            tests, summary = analyse_cell_type(
                sub_root / ct, groups, args.group_a, args.group_b, args.min_donors,
                args.max_donor_fraction, args.donor_min_fraction, log)
        except Exception as exc:  # recorded, never swallowed
            log(f"  FAILED: {type(exc).__name__}: {exc}")
            failures.append({"cell_type": ct, "status": f"failed: {type(exc).__name__}: {exc}"})
            continue
        tests.insert(0, "cell_type", ct)
        tests["average_rho"] = summary["average_rho"]
        all_tests.append(tests)
        summaries.append({"cell_type": ct, **summary})
        log(f"  cells {summary['n_cells']}, subclusters {summary['n_subclusters']}, tested "
            f"{summary['n_subclusters_tested']}, rho {summary['average_rho']:.6f}, "
            f"min donor p {summary['min_donor_p']:.6g}, status {summary['status']}")

    enrich = pd.concat(all_tests, ignore_index=True) if all_tests else pd.DataFrame()
    table = pd.DataFrame(summaries + failures)
    order = order_states(table["cell_type"], _read_order(args.order))
    table = table.set_index("cell_type").loc[order].reset_index()
    enrich_path = outdir / "subcluster_age_enrichment.tsv"
    table_path = outdir / "age_induced_distance.tsv"
    enrich.to_csv(enrich_path, sep="\t", index=False)
    table.to_csv(table_path, sep="\t", index=False)
    log(f"[write] {enrich_path} ({len(enrich)} subcluster rows)")
    log(f"[write] {table_path} ({len(table)} cell states; "
        f"{int((table['status'] == 'ok').sum())} ok, {len(failures)} failed)")

    pdf_path = outdir / "age_induced_distance_dotplot.pdf"
    plotted = table[table["status"] == "ok"]
    report = plot_dotplot(plotted, pdf_path, args.group_a, args.group_b,
                          [c for c in order if c in set(plotted["cell_type"])], args.hline,
                          args.cmap, args.size_scale, args.alpha,
                          (args.fig_width, args.fig_height))
    log(f"[write] {pdf_path}")

    params = {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}
    params.update({"n_cell_types": len(cell_types), "n_ok": int((table["status"] == "ok").sum()),
                   "n_failed": len(failures), "pdf": report, "python": sys.version,
                   "elapsed_s": round(time.time() - t0, 1)})
    with open(outdir / "parameters.json", "w") as fh:
        json.dump(params, fh, indent=2)
    log(f"[write] {outdir / 'parameters.json'}")
    log_fh.close()


def _read_order(path: Optional[str]) -> Optional[List[str]]:
    if not path:
        return None
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--subclusters-dir", default=None,
                   help="AltAnalyze2 ICGS2 sub-clustering output: one folder per cell state")
    p.add_argument("--groups", default=None, help="tab file: uid, donor, group columns")
    p.add_argument("--outdir", required=True)
    p.add_argument("--plot-only", action="store_true",
                   help="re-draw the dot plot from age_induced_distance.tsv in --outdir")
    p.add_argument("--uid-col", default="uid")
    p.add_argument("--donor-col", default="donor")
    p.add_argument("--group-col", default="age-group")
    p.add_argument("--group-a", default="young", help="first group label in --group-col")
    p.add_argument("--group-b", default="old", help="second group label in --group-col")
    p.add_argument("--min-donors", type=int, default=4,
                   help="a subcluster is tested when it holds at least this many donors")
    p.add_argument("--max-donor-fraction", type=float, default=0.9,
                   help="and no donor contributes more than this fraction of its cells")
    p.add_argument("--donor-min-fraction", type=float, default=0.025,
                   help="donor-level test counts a donor with more than this fraction of the subcluster")
    p.add_argument("--cell-types", default=None, help="comma list to restrict the run")
    p.add_argument("--order", default=None,
                   help="file with one cell state per line for the x axis (default: hematopoietic order)")
    p.add_argument("--hline", type=float, default=None, help="dotted guide line at this distance")
    p.add_argument("--cmap", default="coolwarm")
    p.add_argument("--size-scale", type=float, default=20.0,
                   help="point area in pt^2 per unit of -log10 p (matplotlib scatter 's' scale)")
    p.add_argument("--alpha", type=float, default=0.7)
    p.add_argument("--fig-width", type=float, default=12.0)
    p.add_argument("--fig-height", type=float, default=4.0)
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    if not args.plot_only and (args.subclusters_dir is None or args.groups is None):
        build_parser().error("--subclusters-dir and --groups are required unless --plot-only")
    run(args)


if __name__ == "__main__":
    main()
