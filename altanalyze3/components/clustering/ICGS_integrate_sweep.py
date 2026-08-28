"""Sweep integration parameters and score every configuration against the ground-truth standard.

A configuration is judged by whether it recovers the cell types that two or more independent
studies each resolve, not by how many states it produces. A run with many states that loses
reproducible populations is worse than a smaller run that keeps them.

Acceptance requires all of:
  - more than min_states reference states
  - zero ground-truth false negatives
  - every final state holding its required unique markers

Usage:
    python -m altanalyze3.components.clustering.ICGS_integrate_sweep \\
        --run Adams=/path --run Basil2022=/path --run PedDev=/path --run Natri2024=/path \\
        --workdir /path/sweep --min-states 80
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from typing import Dict, List

import pandas as pd

PY = sys.executable
MOD = "altanalyze3.components.clustering.ICGS_integrate"
GT = "altanalyze3.components.clustering.ICGS_integrate_groundtruth"


def run_one(runs: Dict[str, str], outdir: str, params: Dict[str, object],
            cells_per_cluster: int, final_cells: int) -> Dict[str, object]:
    cmd = [PY, "-u", "-m", MOD, "--output-dir", outdir,
           "--cells-per-cluster", str(cells_per_cluster),
           "--final-cells-per-state", str(final_cells)]
    for n, p in runs.items():
        cmd += ["--run", f"{n}={p}"]
    for k, v in params.items():
        cmd += [f"--{k.replace('_', '-')}", str(v)]
    log = os.path.join(outdir + ".log")
    os.makedirs(outdir, exist_ok=True)
    with open(log, "w") as handle:
        rc = subprocess.call(cmd, stdout=handle, stderr=subprocess.STDOUT)
    out: Dict[str, object] = dict(params)
    out["returncode"] = rc
    out["log"] = log
    states_file = os.path.join(outdir, "harmonized_states.tsv")
    markers_file = os.path.join(outdir, "final_markers.tsv")
    if rc != 0 or not os.path.exists(states_file):
        out["n_states"] = 0
        out["status"] = "failed"
        return out
    states = pd.read_csv(states_file, sep="\t")
    out["n_states"] = int(states.shape[0])
    out["n_unique_ids"] = int(states["state"].nunique())
    if os.path.exists(markers_file):
        m = pd.read_csv(markers_file, sep="\t")
        out["states_with_markers"] = int(states["state"].astype(str)
                                         .isin(set(m["top_cluster"].astype(str))).sum())
        out["n_markers"] = int(m.shape[0])
    else:
        out["states_with_markers"] = 0
        out["n_markers"] = 0
    return out


def score_groundtruth(runs: Dict[str, str], outdir: str) -> Dict[str, object]:
    cmd = [PY, "-m", GT, "--integration", outdir]
    for n, p in runs.items():
        cmd += ["--run", f"{n}={p}"]
    try:
        subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
    except subprocess.CalledProcessError as exc:
        return {"gt_total": 0, "gt_represented": 0, "gt_missing": "error", "gt_error": exc.output[-200:]}
    f = os.path.join(outdir, "GROUND_TRUTH.tsv")
    if not os.path.exists(f):
        return {"gt_total": 0, "gt_represented": 0, "gt_missing": "no file"}
    gt = pd.read_csv(f, sep="\t")
    miss = gt.loc[~gt["represented"], "cell_type"].tolist()
    return {"gt_total": int(gt.shape[0]), "gt_represented": int(gt["represented"].sum()),
            "gt_false_negatives": len(miss), "gt_missing": ";".join(miss[:8])}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", action="append", required=True, metavar="NAME=PATH")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--min-states", type=int, default=80)
    ap.add_argument("--cells-per-cluster", type=int, default=200)
    ap.add_argument("--final-cells-per-state", type=int, default=100)
    args = ap.parse_args()

    runs = dict(spec.split("=", 1) for spec in args.run)
    os.makedirs(args.workdir, exist_ok=True)

    grid: List[Dict[str, object]] = []
    # The overlap requirement is a fraction of the smaller marker set, because set sizes span
    # 1 to 449 genes. The damage floor controls how much an existing state may give up.
    # damage_floor made no difference across 1, 2 and 3, because a damaged reference state
    # drops to zero markers, so any floor triggers identically. Fixed at 3 and dropped from the
    # grid. The overlap fraction is the parameter that moves the state count.
    # survival_min_markers is the evidence a candidate must show to count as a distinct
    # population. Three is a very low bar: the Migratory DC clusters carried 19 and 26 unique
    # markers, so the interesting question is how many genuinely distinct populations survive a
    # meaningful requirement, not how few survive a trivial one.
    for frac in (0.33, 0.50):
        for min_markers in (3, 5, 10):
            for floor in (1, min_markers):
                grid.append({"nomination_overlap_fraction": frac,
                             "nomination_specificity": 0.30,
                             "damage_floor": floor,
                             "nomination_min_overlap": 3,
                             "nomination_fdr": 0.05,
                             "survival_rho": 0.4,
                             "survival_min_markers": min_markers,
                             "survival_ref_cells": 60})

    rows = []
    for i, params in enumerate(grid, start=1):
        tag = (f"fr{str(params['nomination_overlap_fraction']).replace('.', '')}_"
               f"mm{params['survival_min_markers']}_"
               f"fl{params['damage_floor']}")
        outdir = os.path.join(args.workdir, tag)
        print(f"[{i}/{len(grid)}] {tag}", flush=True)
        res = run_one(runs, outdir, params, args.cells_per_cluster, args.final_cells_per_state)
        if res.get("n_states", 0):
            res.update(score_groundtruth(runs, outdir))
        res["tag"] = tag
        res["passes"] = bool(res.get("n_states", 0) > args.min_states
                             and res.get("gt_false_negatives", 99) == 0
                             and res.get("states_with_markers", 0) == res.get("n_states", -1))
        rows.append(res)
        print(f"    states={res.get('n_states')} "
              f"gt={res.get('gt_represented')}/{res.get('gt_total')} "
              f"markers_ok={res.get('states_with_markers')} passes={res['passes']}", flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(args.workdir, "sweep_results.tsv"),
                                  sep="\t", index=False)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.workdir, "sweep_results.tsv"), sep="\t", index=False)
    ok = df.loc[df["passes"]]
    print(f"\nconfigurations passing every criterion: {ok.shape[0]} of {df.shape[0]}")
    if ok.shape[0]:
        print(ok[["tag", "n_states", "gt_represented", "gt_total"]].to_string(index=False))
    else:
        print("none passed; closest by ground truth:")
        print(df.nlargest(5, "gt_represented")[
            ["tag", "n_states", "gt_represented", "gt_total", "gt_false_negatives"]
        ].to_string(index=False))
    print(f"\nwrote {os.path.join(args.workdir, 'sweep_results.tsv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
