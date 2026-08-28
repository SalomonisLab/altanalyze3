"""Benchmarks for a cross-study ICGS3 integration result.

Each benchmark answers one question about whether the integrated reference is defensible. The
module reads only files the integration already wrote, plus the source ICGS3 runs, and writes a
report. It changes nothing.

Benchmarks
  1. Label recovery      every annotation label present in the input must reach a final state
  2. Cell recovery       what fraction of input cells a surviving state accounts for
  3. Known biology       named populations must each map to a final state
  4. Marker self-consistency  every final state must hold its required unique markers
  5. Redundancy audit    the reason each excluded cluster was excluded, by annotation label

Usage:
    python -m altanalyze3.components.clustering.ICGS_integrate_benchmark \\
        --integration /path/ICGS3_integration_v2 \\
        --run Adams=/path/ICGS3_Adams_corr03_ne3 ... \\
        --out /path/ICGS3_integration_v2/BENCHMARKS.md
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, List

import pandas as pd

KNOWN_PANELS = {
    "aberrant basaloid": ["CDH2", "ITGB6", "TP63", "MMP7", "GDF15", "EPHB2", "CDKN2A"],
    "CTHRC1 fibroblast": ["CTHRC1", "POSTN", "COMP", "FAP", "TNC", "ASPN"],
    "neuroendocrine": ["CALCA", "ASCL1", "CHGA", "GRP"],
    "mesothelium": ["ITLN1", "PRG4", "UPK3B", "CALB2"],
    "proliferating": ["MKI67", "TOP2A", "UBE2C", "BIRC5", "PCLAF"],
}


def annotation_of_clusters(run_dir: str) -> pd.Series:
    """Dominant annotation label per ICGS3 cluster, from that run's own cell table."""
    cells = pd.read_csv(os.path.join(run_dir, "icgs3_cell_barcode_clusters.tsv"), sep="\t")
    col = "HLCA" if "HLCA" in cells.columns else next(
        (c for c in cells.columns if c not in
         ("barcode", "Library", "sample", "ICGS3_cluster", "ICGS3_original_NMF_cluster",
          "ICGS3_SVM_score", "ICGS3_SVM_margin", "ICGS3_cell_state_prediction")), None)
    cells["ICGS3_cluster"] = cells["ICGS3_cluster"].astype(str)
    return cells.groupby("ICGS3_cluster")[col].agg(lambda s: s.value_counts().idxmax())


def cells_per_cluster(run_dir: str) -> pd.Series:
    cells = pd.read_csv(os.path.join(run_dir, "icgs3_cell_barcode_clusters.tsv"),
                        sep="\t", usecols=["ICGS3_cluster"])
    return cells["ICGS3_cluster"].astype(str).value_counts()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--integration", required=True)
    ap.add_argument("--run", action="append", required=True, metavar="NAME=PATH")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    D = args.integration
    runs = dict(spec.split("=", 1) for spec in args.run)
    ann = {n: annotation_of_clusters(p) for n, p in runs.items()}
    sizes = {n: cells_per_cluster(p) for n, p in runs.items()}

    states = pd.read_csv(os.path.join(D, "harmonized_states.tsv"), sep="\t")
    dec = pd.read_csv(os.path.join(D, "nomination_decisions.tsv"), sep="\t")
    markers = pd.read_csv(os.path.join(D, "final_markers.tsv"), sep="\t")

    dec["cluster"] = dec["cluster"].astype(str)
    dec["label"] = [ann[r].get(c, "?") for r, c in zip(dec["step"], dec["cluster"])]
    dec["n_cells"] = [int(sizes[r].get(c, 0)) for r, c in zip(dec["step"], dec["cluster"])]

    # every label carried by any state member, and by any cluster the state absorbed as redundant
    state_label: Dict[str, List[str]] = {}
    for _, row in states.iterrows():
        labels = []
        for m in str(row["members"]).split(","):
            if "|" in m:
                d, cl = m.split("|", 1)
                if d in ann and cl in ann[d].index:
                    labels.append(ann[d][cl])
        state_label[str(row["state"])] = labels
    represented = {l for v in state_label.values() for l in v}
    # a redundant cluster is represented by the state it was found redundant with
    for _, r in dec.loc[dec["action"] == "excluded_redundant"].iterrows():
        tgt = str(r.get("nearest_state_r025", "") or "")
        if tgt and tgt in state_label and str(r["label"]) != "?":
            represented.add(str(r["label"]))

    all_labels = {l for a in ann.values() for l in a.values}
    missing = sorted(all_labels - represented)

    total_cells = int(sum(int(s.sum()) for s in sizes.values()))
    kept_cells = int(dec.loc[dec["action"] == "added_new_state", "n_cells"].sum())
    seed_cells = total_cells - int(dec["n_cells"].sum())
    lost_cells = int(dec.loc[dec["action"] == "excluded_redundant", "n_cells"].sum())

    lines: List[str] = []
    add = lines.append
    add("# Benchmarks for the integrated reference\n")
    add(f"Integration directory: `{D}`\n")

    add("## 1. Label recovery\n")
    add(f"Annotation labels dominant in at least one input cluster: {len(all_labels)}.")
    add(f"Labels represented in the final reference: {len(represented & all_labels)}.")
    add(f"Labels absent: {len(missing)}.\n")
    if missing:
        add("Absent labels and the number of input clusters carrying them:\n")
        for lab in missing:
            n = int((dec["label"] == lab).sum())
            add(f"- {lab}: {n} excluded cluster(s)")
        add("")

    add("## 2. Cell recovery\n")
    add(f"Cells in all input clusters: {total_cells}.")
    add(f"Cells in seed clusters, admitted without testing: {seed_cells}.")
    add(f"Cells in clusters admitted as new states: {kept_cells}.")
    add(f"Cells in clusters excluded: {lost_cells} "
        f"({100.0 * lost_cells / max(total_cells, 1):.1f}% of all input cells).\n")

    add("## 3. Known biology\n")
    for name, panel in KNOWN_PANELS.items():
        hit = markers.loc[markers["marker"].isin(panel)]
        if hit.empty:
            add(f"- {name}: ABSENT, none of {len(panel)} panel genes is a marker of any state")
        else:
            top = hit["top_cluster"].value_counts()
            add(f"- {name}: {top.iloc[0]} of {len(panel)} panel genes mark state {top.index[0]}")
    add("")

    add("## 4. Marker self-consistency\n")
    with_markers = set(markers["top_cluster"].astype(str))
    n_ok = int(states["state"].astype(str).isin(with_markers).sum())
    add(f"States holding at least one marker in the final analysis: {n_ok} of {states.shape[0]}.\n")

    add("## 5. Redundancy audit\n")
    audit = dec.loc[dec["action"] == "excluded_redundant"].copy()
    audit["reason"] = audit["withdrawn_reason"].fillna("").str.split("(").str[0].str.strip()
    add("Exclusions by stage and reason:\n")
    tally = audit.groupby(["step", "reason"]).size().reset_index(name="clusters")
    for _, r in tally.iterrows():
        add(f"- {r['step']}: {r['clusters']} cluster(s), {r['reason'][:70]}")
    add("")

    report = "\n".join(lines)
    out = args.out or os.path.join(D, "BENCHMARKS.md")
    with open(out, "w") as handle:
        handle.write(report)
    print(report)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
