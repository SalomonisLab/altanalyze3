"""Per-cell sub-cluster labels for an AltAnalyze2 ICGS2 ``--subcluster --k 2`` run.

Why this module exists
----------------------
A k=2 sub-cluster run leaves three different answers on disk and writes none of them as a final
per-cell label file:

1. ``NMF-SVM/SubtypeAnalyses/round1Metadata<rank>.txt`` and
   ``NMF-SVM/ExpressionInput/groups.NMF-MarkerFinder.txt`` hold the NMF round-1 membership. Both are
   written inside ``NMF_Analysis.NMFAnalysis`` (altanalyze2 ``stats_scripts/NMF_Analysis.py``), before
   the SVM runs. The label is the larger of two loadings with no margin, so a cell whose loadings are
   0.51 against 0.49 is reported as a hard call.
2. ``NMF-SVM/SVMOutputs/round1SVC_decision_func.txt`` holds the supervised re-assignment that
   ``ExpandSampleClusters`` computed next. For two clusters ``LinearSVC`` gives one decision value per
   cell, and altanalyze2 never writes the predicted labels for a non-PSI platform: the line that would
   export them is commented out, and ``round1SVC_Results_max.txt`` receives the sign of that single
   value duplicated into both columns.
3. The MarkerFinder genes of the split, which is the space the sub-clusters were defined in.

``svm`` recovers answer 2 from the decision values. ``svm+centroid`` then runs nearest-centroid passes
in the MarkerFinder gene space until the labels stop changing, which resolves the cells whose SVM
margin is small. Measured on the 2024-06-10 aged-versus-young run, against the cells a UMAP separates
into two clouds: HSC-1 11 disagreements under ``svm`` and 1 under ``svm+centroid``, MPP-1 10 and 3,
Intermediate Mono-2 0 and 0.

A state whose cells form one continuous cloud has no two clouds to agree with, and there the centroid
passes move many cells rather than a few. ``labels_summary.tsv`` reports ``n_moved_centroid`` per state
so that case is visible before the labels are used.

Outputs (``--outdir``)
  <cell state>.txt       uid, V0, V1 with a header; the format ``subcluster_age_distance
                         --annotation-dir`` and ``umap_donor_groups --subclusters`` read
  labels_summary.tsv     one row per state: cells, method, cells moved at each step, convergence
  parameters.json        every parameter and input path
  run_log.txt            per-state log
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from altanalyze3.components.visualization.subcluster_age_distance import read_marker_expression
except ImportError:  # pragma: no cover - direct script execution
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from subcluster_age_distance import read_marker_expression

METHODS = ("nmf", "svm", "svm+centroid")
MAX_PASSES = 20


def read_nmf_labels(cell_type_dir: Path) -> pd.DataFrame:
    """The NMF round-1 membership: uid and the V column that holds a 1."""
    sub = cell_type_dir / "NMF-SVM" / "SubtypeAnalyses"
    hits = sorted(glob.glob(str(sub / "round1Metadata*.txt")))
    if not hits:
        raise FileNotFoundError(f"no round1Metadata*.txt in {sub}")
    meta = pd.read_csv(hits[0], sep="\t")
    v_cols = [c for c in meta.columns if c != "uid"]
    if len(v_cols) != 2:
        raise ValueError(f"{hits[0]} holds {len(v_cols)} components; this module handles k=2 only")
    flags = meta[v_cols].astype(int)
    if not bool((flags.sum(axis=1) == 1).all()):
        raise ValueError(f"{hits[0]}: every cell must carry exactly one component flag")
    return pd.DataFrame({"uid": meta["uid"], "label": flags.idxmax(axis=1).values})


def read_svm_decision(cell_type_dir: Path) -> pd.Series:
    """The one LinearSVC decision value per cell that altanalyze2 kept for a k=2 run."""
    path = cell_type_dir / "NMF-SVM" / "SVMOutputs" / "round1SVC_decision_func.txt"
    if not path.exists():
        raise FileNotFoundError(f"no SVM decision values at {path}")
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "uid"})
    values = df[df.columns[1:]].astype(float)
    if values.shape[1] > 1 and not bool((values.nunique(axis=1) == 1).all()):
        raise ValueError(f"{path}: columns differ, so this is not the k=2 single-decision case")
    return pd.Series(values.iloc[:, 0].values, index=df["uid"].astype(str), name="decision")


def svm_labels(nmf: pd.DataFrame, decision: pd.Series) -> Tuple[pd.Series, int]:
    """Sign of the decision value, oriented so each side keeps the NMF name it overlaps most.

    Returns the labels and the orientation margin, the difference in overlap between the two
    possible orientations. A margin of 0 would mean the orientation is a coin toss.
    """
    lab = nmf.set_index("uid")["label"]
    dec = decision.reindex(lab.index)
    if dec.isna().any():
        missing = int(dec.isna().sum())
        raise ValueError(f"{missing} of {len(lab)} cells have no SVM decision value")
    names = sorted(lab.unique())
    positive, negative = names[-1], names[0]
    agree = int(((dec > 0) & (lab == positive)).sum() + ((dec <= 0) & (lab == negative)).sum())
    flipped = int(((dec > 0) & (lab == negative)).sum() + ((dec <= 0) & (lab == positive)).sum())
    if flipped > agree:
        positive, negative = negative, positive
    out = pd.Series(np.where(dec > 0, positive, negative), index=lab.index, name="label")
    return out, abs(agree - flipped)


def centroid_passes(expr: pd.DataFrame, labels: pd.Series,
                    max_passes: int = MAX_PASSES) -> Tuple[pd.Series, int, bool, int]:
    """Assign each cell to the centroid it correlates with most, until the labels stop changing.

    The centroids are means over the MarkerFinder genes of the split, which is the space the
    sub-clusters were defined in.

    A cell whose MarkerFinder values are all equal has no variance, so its correlation to both
    centroids is undefined. Such a cell keeps the label it arrived with and is counted, never
    dropped and never given an arbitrary one. Measured on the 2024-06-10 run, 4 of 87 cell states
    hold at least one.

    Returns the labels, the number of passes, whether it settled, and how many cells were
    undecidable on the last pass.
    """
    cells = [c for c in labels.index if c in expr.columns]
    if len(cells) != len(labels):
        raise ValueError(f"{len(labels) - len(cells)} of {len(labels)} cells are absent from the "
                         f"expression matrix")
    X = expr[cells]
    current = labels.loc[cells]
    undecidable = 0
    for passes in range(1, max_passes + 1):
        names = sorted(current.unique())
        if len(names) < 2:
            return current, passes, True, undecidable
        centroids = pd.DataFrame({n: X[current[current == n].index].mean(axis=1) for n in names})
        r = pd.DataFrame({n: X.corrwith(centroids[n]) for n in names})
        no_call = r.isna().all(axis=1)
        undecidable = int(no_call.sum())
        nxt = current.copy()
        decidable = r.loc[~no_call]
        if len(decidable):
            nxt.loc[decidable.index] = decidable.idxmax(axis=1)
        nxt.name = "label"
        if nxt.equals(current):
            return current, passes, True, undecidable
        current = nxt
    return current, max_passes, False, undecidable


def labels_for_state(cell_type_dir: Path, method: str, log) -> Tuple[pd.DataFrame, dict]:
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, not {method!r}")
    nmf = read_nmf_labels(cell_type_dir)
    row = {"cell_type": cell_type_dir.name, "n_cells": len(nmf), "method": method,
           "n_moved_svm": 0, "n_moved_centroid": 0, "centroid_passes": 0,
           "centroid_converged": True, "n_undecidable": 0, "orientation_margin": ""}
    labels = nmf.set_index("uid")["label"]
    if method != "nmf":
        labels, margin = svm_labels(nmf, read_svm_decision(cell_type_dir))
        row["orientation_margin"] = margin
        row["n_moved_svm"] = int((labels != nmf.set_index("uid")["label"]).sum())
    if method == "svm+centroid":
        expr = read_marker_expression(cell_type_dir, log)
        before = labels
        labels, passes, converged, undecidable = centroid_passes(expr, labels)
        row["n_moved_centroid"] = int((labels != before.loc[labels.index]).sum())
        row["centroid_passes"] = passes
        row["centroid_converged"] = converged
        row["n_undecidable"] = undecidable
        if not converged:
            log(f"  WARNING centroid passes did not settle in {MAX_PASSES} passes")
        if undecidable:
            log(f"  WARNING {undecidable} of {len(labels)} cells have no variance over the "
                f"MarkerFinder genes; each kept the label it arrived with")
    names = sorted(nmf["label"].unique())
    out = pd.DataFrame({"uid": labels.index,
                        names[0]: (labels == names[0]).astype(int).values,
                        names[1]: (labels == names[1]).astype(int).values})
    row["n_" + names[0]] = int(out[names[0]].sum())
    row["n_" + names[1]] = int(out[names[1]].sum())
    return out, row


def main(argv: Optional[List[str]] = None) -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--subclusters-dir", required=True,
                   help="one folder per cell state, as ICGS2 --subcluster writes them")
    p.add_argument("--outdir", required=True)
    p.add_argument("--method", default="svm+centroid", choices=list(METHODS),
                   help="nmf: the run's round-1 labels; svm: the supervised re-assignment; "
                        "svm+centroid (default): svm then nearest-centroid passes")
    p.add_argument("--cell-types", default=None, help="comma list to restrict the run")
    args = p.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_fh = open(outdir / "run_log.txt", "w")

    def log(msg: str) -> None:
        print(msg)
        log_fh.write(msg + "\n")
        log_fh.flush()

    root = Path(args.subclusters_dir)
    states = sorted(d.name for d in root.iterdir() if d.is_dir() and not d.name.startswith("."))
    if args.cell_types:
        wanted = [c.strip() for c in args.cell_types.split(",") if c.strip()]
        absent = [c for c in wanted if c not in states]
        if absent:
            raise ValueError(f"requested cell types without a folder: {absent}")
        states = wanted
    log(f"[subclusters] {root}: {len(states)} cell-state folders; method {args.method}")

    rows, failures = [], []
    for i, state in enumerate(states, 1):
        log(f"[{i}/{len(states)}] {state}")
        try:
            table, row = labels_for_state(root / state, args.method, log)
        except Exception as exc:  # recorded, never swallowed
            log(f"  FAILED: {type(exc).__name__}: {exc}")
            failures.append({"cell_type": state, "status": f"failed: {type(exc).__name__}: {exc}"})
            continue
        table.to_csv(outdir / f"{state}.txt", sep="\t", index=False)
        rows.append(row)
        log(f"  {row['n_cells']} cells; moved by SVM {row['n_moved_svm']}, "
            f"by centroid {row['n_moved_centroid']} in {row['centroid_passes']} passes")

    summary = pd.DataFrame(rows)
    if failures:
        summary = pd.concat([summary, pd.DataFrame(failures)], ignore_index=True)
    summary.to_csv(outdir / "labels_summary.tsv", sep="\t", index=False)
    log(f"[write] {outdir / 'labels_summary.tsv'} ({len(rows)} states, {len(failures)} failed)")
    if rows:
        done = pd.DataFrame(rows)
        log(f"[totals] {int(done.n_cells.sum())} cells; SVM moved {int(done.n_moved_svm.sum())}; "
            f"centroid moved {int(done.n_moved_centroid.sum())}")
    with open(outdir / "parameters.json", "w") as fh:
        json.dump({k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
                  fh, indent=2)
    log(f"[write] {outdir / 'parameters.json'}")
    log_fh.close()
    if failures:
        raise SystemExit(f"{len(failures)} of {len(states)} cell states failed; see run_log.txt")


if __name__ == "__main__":
    main()
