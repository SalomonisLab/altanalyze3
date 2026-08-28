"""Ground-truth test for a cross-study ICGS3 integration, independent of the integration itself.

The integration decides which clusters are redundant. A benchmark that reuses that decision, or
that compares annotation label strings across datasets annotated with different vocabularies,
cannot detect a wrong answer. This module builds an external standard instead.

The standard: a cell type that at least two independent HLCA-annotated studies each resolve as a
distinct cluster is a real population. The integrated reference must contain at least one state
whose marker genes identify that population. A reproducible cell type absent from the reference
is a false negative of the integration, and is reported as such.

The seed dataset is excluded from the standard when it uses a different annotation vocabulary,
because its label strings cannot be compared. PedDev annotates with `Hs-PedDev-Lung`, which
shares only 9 label strings with `HLCA`, so comparing them by name is invalid.

Usage:
    python -m altanalyze3.components.clustering.ICGS_integrate_groundtruth \\
        --integration /path/ICGS3_integration_v3 \\
        --run Adams=/path/ICGS3_Adams_corr03_ne3 \\
        --run Basil2022=/path/ICGS3_Basil2022 \\
        --run Natri2024=/path/ICGS3_Natri2024 \\
        --run PedDev=/path/ICGS3_PedDev \\
        --annotation HLCA
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, List, Set

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def cluster_annotation(run_dir: str, column: str) -> pd.DataFrame:
    """Dominant annotation label and its purity, per ICGS3 cluster."""
    cells = pd.read_csv(os.path.join(run_dir, "icgs3_cell_barcode_clusters.tsv"), sep="\t")
    if column not in cells.columns:
        return pd.DataFrame(columns=["cluster", "label", "purity", "n_cells"])
    cells["ICGS3_cluster"] = cells["ICGS3_cluster"].astype(str)
    rows = []
    for cl, grp in cells.groupby("ICGS3_cluster"):
        vc = grp[column].value_counts()
        rows.append({"cluster": cl, "label": vc.index[0],
                     "purity": float(vc.iloc[0] / grp.shape[0]), "n_cells": int(grp.shape[0])})
    return pd.DataFrame(rows)


def cluster_markers(run_dir: str, top_n: int = 50) -> Dict[str, List[str]]:
    mk = pd.read_csv(os.path.join(run_dir, "MarkerFinder", "icgs3_markers.tsv"), sep="\t")
    out: Dict[str, List[str]] = {}
    for cl, grp in mk.groupby("top_cluster"):
        out[str(cl)] = grp.sort_values("pearson_r", ascending=False)["marker"].astype(str) \
            .head(top_n).tolist()
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--integration", required=True)
    ap.add_argument("--run", action="append", required=True, metavar="NAME=PATH")
    ap.add_argument("--annotation", default="HLCA")
    ap.add_argument("--min-studies", type=int, default=2,
                    help="Studies that must independently resolve a label for it to count as a "
                         "reproducible cell type. Default 2.")
    ap.add_argument("--min-purity", type=float, default=0.5,
                    help="Dominant-label purity a cluster needs before it counts as resolving "
                         "that label. Default 0.5.")
    ap.add_argument("--overlap-fraction", type=float, default=0.25,
                    help="Shared markers required, as a fraction of the smaller of the query "
                         "and the state's marker set. A fixed count cannot serve both: the "
                         "median state carries 17 markers in total, so 10 would demand 59% of "
                         "everything it has. Default 0.25.")
    ap.add_argument("--min-overlap-floor", type=int, default=3,
                    help="Absolute floor on the shared-marker requirement. Default 3.")
    ap.add_argument("--fdr", type=float, default=0.05,
                    help="Benjamini-Hochberg adjusted enrichment FDR a state must reach for the "
                         "cell type to count as represented. Default 0.05.")
    args = ap.parse_args()

    runs = dict(spec.split("=", 1) for spec in args.run)
    ann = {n: cluster_annotation(p, args.annotation) for n, p in runs.items()}
    usable = {n: a for n, a in ann.items() if not a.empty}
    print(f"studies annotated with '{args.annotation}': {sorted(usable)}")
    print(f"studies excluded, different vocabulary: "
          f"{sorted(set(runs) - set(usable)) or 'none'}\n")

    # a reproducible cell type: resolved as a pure cluster by at least min_studies studies
    per_label_studies: Dict[str, Set[str]] = {}
    exemplar: Dict[str, List[str]] = {}
    for name, frame in usable.items():
        good = frame.loc[frame["purity"] >= args.min_purity]
        for _, row in good.iterrows():
            per_label_studies.setdefault(row["label"], set()).add(name)
            exemplar.setdefault(row["label"], []).append(f"{name}|{row['cluster']}")
    reproducible = sorted(l for l, st in per_label_studies.items()
                          if len(st) >= int(args.min_studies))
    print(f"reproducible cell types, resolved by >= {args.min_studies} studies at purity "
          f">= {args.min_purity}: {len(reproducible)}\n")

    # marker sets of the reference states
    D = args.integration
    fm = pd.read_csv(os.path.join(D, "final_markers.tsv"), sep="\t")
    state_markers: Dict[str, Set[str]] = {
        str(k): set(v) for k, v in fm.groupby("top_cluster")["marker"].apply(set).items()}
    universe = int(pd.read_csv(os.path.join(D, "harmonized_final_centroids.tsv"),
                               sep="\t", usecols=[0]).shape[0])
    markers_by_run = {n: cluster_markers(p) for n, p in runs.items()}

    rows = []
    for label in reproducible:
        query: Set[str] = set()
        for ex in exemplar[label]:
            n, cl = ex.split("|", 1)
            query |= set(markers_by_run.get(n, {}).get(cl, []))
        # score every state, then keep the best by significance rather than by raw count
        cand = []
        for state, genes in state_markers.items():
            k = len(query & genes)
            if k == 0:
                continue
            p = float(hypergeom.sf(k - 1, universe, len(query), len(genes)))
            need = max(int(args.min_overlap_floor),
                       int(np.ceil(float(args.overlap_fraction) * min(len(query), len(genes)))))
            cand.append({"state": state, "k": k, "p": p, "need": need,
                         "state_markers": len(genes)})
        if not cand:
            rows.append({"cell_type": label, "studies_resolving": len(per_label_studies[label]),
                         "exemplar_clusters": ",".join(exemplar[label]),
                         "query_markers": len(query), "best_state": "", "shared_markers": 0,
                         "state_markers": 0, "required_overlap": 0, "p_value": 1.0, "fdr": 1.0,
                         "represented": False})
            continue
        cf = pd.DataFrame(cand)
        cf["fdr"] = _bh(cf["p"].to_numpy())
        cf = cf.sort_values(["p", "k"], ascending=[True, False])
        top = cf.iloc[0]
        ok = bool(int(top["k"]) >= int(top["need"]) and float(top["fdr"]) <= float(args.fdr))
        rows.append({"cell_type": label, "studies_resolving": len(per_label_studies[label]),
                     "exemplar_clusters": ",".join(exemplar[label]),
                     "query_markers": len(query), "best_state": str(top["state"]),
                     "shared_markers": int(top["k"]), "state_markers": int(top["state_markers"]),
                     "required_overlap": int(top["need"]), "p_value": float(top["p"]),
                     "fdr": float(top["fdr"]), "represented": ok})
    res = pd.DataFrame(rows).sort_values(["represented", "shared_markers"])
    out = os.path.join(D, "GROUND_TRUTH.tsv")
    res.to_csv(out, sep="\t", index=False)

    n_ok = int(res["represented"].sum())
    print(f"reproducible cell types represented in the reference: {n_ok} of {res.shape[0]}")
    missing = res.loc[~res["represented"]]
    if missing.shape[0]:
        print(f"\nFALSE NEGATIVES, reproducible cell types no state represents:\n")
        for _, r in missing.iterrows():
            print(f"  {r['cell_type']}: resolved by {r['studies_resolving']} studies; best state "
                  f"{r['best_state'] or 'none'} shares {r['shared_markers']} of "
                  f"{r['required_overlap']} required markers, FDR {r['fdr']:.2e} "
                  f"({r['exemplar_clusters']})")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
