"""Build shared genomic CNV intervals from per-cell CNV intervals.

This is intentionally interval-first. Chromosome-arm calls are useful summaries
for karyotype comparison, but the primary CNV unit here is a genomic segment
whose support comes from overlapping per-cell/state intervals. Different cell
states can therefore support the same CNV even when expressed genes shift the
observed start/end genes.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd


LOSS_CALLS = {"loss", "homozygous_loss", "deletion", "het_loss", "hom_loss"}
GAIN_CALLS = {"gain", "amplification", "duplication", "amp"}


def _normal_call(value: str) -> str:
    text = str(value).strip().lower()
    if text in LOSS_CALLS:
        return "loss"
    if text in GAIN_CALLS:
        return "gain"
    return text


def _load_regions(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    required = {"CellBarcode", "chr", "start", "end", "call"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    if "sample" not in frame.columns:
        frame["sample"] = "sample"
    if "cell_state" not in frame.columns:
        frame["cell_state"] = "unknown"
    frame = frame.copy()
    frame["CellBarcode"] = frame["CellBarcode"].astype(str)
    frame["sample"] = frame["sample"].astype(str)
    frame["cell_state"] = frame["cell_state"].astype(str)
    frame["chr"] = frame["chr"].astype(str)
    frame["call_group"] = frame["call"].map(_normal_call)
    frame = frame.loc[frame["call_group"].isin(["loss", "gain"])].copy()
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce").astype("Int64")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce").astype("Int64")
    frame = frame.dropna(subset=["start", "end"]).copy()
    frame["start"] = frame["start"].astype(np.int64)
    frame["end"] = frame["end"].astype(np.int64)
    bad = frame["end"] < frame["start"]
    if bad.any():
        swapped = frame.loc[bad, "start"].copy()
        frame.loc[bad, "start"] = frame.loc[bad, "end"]
        frame.loc[bad, "end"] = swapped
    return frame.reset_index(drop=True)


def _load_cells(path: Optional[Path]) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=["CellBarcode", "sample", "cell_state"])
    frame = pd.read_csv(path, sep="\t")
    required = {"CellBarcode"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    if "sample" not in frame.columns:
        frame["sample"] = "sample"
    if "cell_state" not in frame.columns:
        frame["cell_state"] = "unknown"
    frame["CellBarcode"] = frame["CellBarcode"].astype(str)
    frame["sample"] = frame["sample"].astype(str)
    frame["cell_state"] = frame["cell_state"].astype(str)
    return frame


def _merge_per_cell_intervals(frame: pd.DataFrame, max_gap_bp: int) -> pd.DataFrame:
    rows: list[dict] = []
    sort_cols = ["sample", "chr", "call_group", "CellBarcode", "start", "end"]
    for keys, group in frame.sort_values(sort_cols).groupby(["sample", "chr", "call_group", "CellBarcode"], sort=False):
        sample, chrom, call_group, cell = keys
        state = str(group["cell_state"].iloc[0])
        start_gene = group["start_gene"].iloc[0] if "start_gene" in group.columns else ""
        end_gene = group["end_gene"].iloc[-1] if "end_gene" in group.columns else ""
        cur_start: Optional[int] = None
        cur_end: Optional[int] = None
        starts: list[int] = []
        ends: list[int] = []
        scores: list[float] = []
        for row in group.itertuples(index=False):
            start = int(getattr(row, "start"))
            end = int(getattr(row, "end"))
            if cur_start is None:
                cur_start = start
                cur_end = end
            elif start <= int(cur_end) + max_gap_bp:
                cur_end = max(int(cur_end), end)
            else:
                rows.append(
                    {
                        "sample": sample,
                        "chr": chrom,
                        "call_group": call_group,
                        "CellBarcode": cell,
                        "cell_state": state,
                        "start": int(cur_start),
                        "end": int(cur_end),
                        "start_gene": start_gene,
                        "end_gene": end_gene,
                        "mean_signal": float(np.nanmean(scores)) if scores else np.nan,
                    }
                )
                cur_start = start
                cur_end = end
                starts.clear()
                ends.clear()
                scores.clear()
            starts.append(start)
            ends.append(end)
            if hasattr(row, "mean_log2"):
                scores.append(float(getattr(row, "mean_log2")))
            elif hasattr(row, "mean_score"):
                scores.append(float(getattr(row, "mean_score")))
        if cur_start is not None:
            rows.append(
                {
                    "sample": sample,
                    "chr": chrom,
                    "call_group": call_group,
                    "CellBarcode": cell,
                    "cell_state": state,
                    "start": int(cur_start),
                    "end": int(cur_end),
                    "start_gene": start_gene,
                    "end_gene": end_gene,
                    "mean_signal": float(np.nanmean(scores)) if scores else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _format_counts(series: pd.Series, max_items: int = 8) -> str:
    counts = series.astype(str).value_counts()
    return ";".join(f"{key}({int(value)})" for key, value in counts.head(max_items).items())


def _consensus_for_group(
    group: pd.DataFrame,
    *,
    sample: str,
    chrom: str,
    call_group: str,
    denominator: int,
    min_cells: int,
    min_fraction: float,
    merge_gap_bp: int,
) -> list[dict]:
    if group.empty:
        return []
    points = np.unique(np.r_[group["start"].to_numpy(np.int64), group["end"].to_numpy(np.int64)])
    if points.size < 2:
        return []
    point_index = pd.Series(np.arange(points.size, dtype=np.int64), index=points)
    diff = np.zeros(points.size + 1, dtype=np.int64)
    for row in group.itertuples(index=False):
        start_idx = int(point_index.loc[int(row.start)])
        end_idx = int(point_index.loc[int(row.end)])
        if end_idx <= start_idx:
            end_idx = min(start_idx + 1, points.size - 1)
        diff[start_idx] += 1
        diff[end_idx] -= 1
    coverage = np.cumsum(diff[:-1])
    threshold = max(int(min_cells), int(np.ceil(float(denominator) * float(min_fraction))))
    selected = coverage[:-1] >= threshold
    rows: list[dict] = []
    idx = 0
    while idx < selected.size:
        if not selected[idx]:
            idx += 1
            continue
        start_idx = idx
        end_idx = idx + 1
        while end_idx < selected.size and (
            selected[end_idx] or int(points[end_idx] - points[end_idx - 1]) <= merge_gap_bp
        ):
            end_idx += 1
        start = int(points[start_idx])
        end = int(points[end_idx])
        support = group.loc[(group["start"] < end) & (group["end"] > start)].copy()
        support_cells = support["CellBarcode"].nunique()
        if support_cells < threshold:
            idx = end_idx
            continue
        support_fraction = support_cells / max(denominator, 1)
        rows.append(
            {
                "sample": sample,
                "chr": chrom,
                "call": call_group,
                "start": start,
                "end": end,
                "length_bp": int(end - start),
                "n_support_cells": int(support_cells),
                "support_fraction": support_fraction,
                "denominator_cells": int(denominator),
                "n_cell_states": int(support["cell_state"].nunique()),
                "top_cell_states": _format_counts(support["cell_state"]),
                "n_source_intervals": int(support.shape[0]),
                "mean_signal": float(np.nanmean(support["mean_signal"])) if "mean_signal" in support else np.nan,
                "support_start_q10": int(np.nanquantile(support["start"], 0.10)),
                "support_start_q50": int(np.nanquantile(support["start"], 0.50)),
                "support_start_q90": int(np.nanquantile(support["start"], 0.90)),
                "support_end_q10": int(np.nanquantile(support["end"], 0.10)),
                "support_end_q50": int(np.nanquantile(support["end"], 0.50)),
                "support_end_q90": int(np.nanquantile(support["end"], 0.90)),
            }
        )
        idx = end_idx
    return rows


def build_interval_consensus(
    regions_path: Path,
    cells_path: Optional[Path],
    output_prefix: Path,
    *,
    min_cells: int,
    min_fraction: float,
    per_sample: bool,
    merge_cell_gap_bp: int,
    merge_consensus_gap_bp: int,
) -> dict:
    regions = _load_regions(regions_path)
    cells = _load_cells(cells_path)
    merged = _merge_per_cell_intervals(regions, max_gap_bp=merge_cell_gap_bp)
    if merged.empty:
        consensus = pd.DataFrame()
    else:
        if not per_sample:
            merged["sample"] = "ALL"
            if not cells.empty:
                cells = cells.copy()
                cells["sample"] = "ALL"
        if cells.empty:
            denominators = merged.groupby("sample")["CellBarcode"].nunique().to_dict()
        else:
            denominators = cells.groupby("sample")["CellBarcode"].nunique().to_dict()
        rows: list[dict] = []
        for (sample, chrom, call_group), group in merged.groupby(["sample", "chr", "call_group"], sort=False):
            denominator = int(denominators.get(sample, group["CellBarcode"].nunique()))
            rows.extend(
                _consensus_for_group(
                    group,
                    sample=str(sample),
                    chrom=str(chrom),
                    call_group=str(call_group),
                    denominator=denominator,
                    min_cells=min_cells,
                    min_fraction=min_fraction,
                    merge_gap_bp=merge_consensus_gap_bp,
                )
            )
        consensus = pd.DataFrame(rows)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(f"{output_prefix}.per_cell_merged_intervals.tsv", sep="\t", index=False)
    if consensus.empty:
        consensus = pd.DataFrame(
            columns=[
                "sample", "chr", "call", "start", "end", "length_bp",
                "n_support_cells", "support_fraction", "denominator_cells",
                "n_cell_states", "top_cell_states", "n_source_intervals",
                "mean_signal", "support_start_q10", "support_start_q50",
                "support_start_q90", "support_end_q10", "support_end_q50",
                "support_end_q90",
            ]
        )
    consensus.sort_values(
        ["sample", "support_fraction", "n_support_cells", "chr", "start"],
        ascending=[True, False, False, True, True],
    ).to_csv(f"{output_prefix}.cnv_interval_consensus.tsv", sep="\t", index=False)
    return {
        "source_intervals": int(regions.shape[0]),
        "merged_cell_intervals": int(merged.shape[0]),
        "consensus_intervals": int(consensus.shape[0]),
        "output_prefix": str(output_prefix),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build shared genomic CNV interval consensus.")
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--cells", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--min-fraction", type=float, default=0.01)
    parser.add_argument("--all-samples", action="store_true", help="Build one consensus across all samples.")
    parser.add_argument("--merge-cell-gap-bp", type=int, default=0)
    parser.add_argument("--merge-consensus-gap-bp", type=int, default=0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    summary = build_interval_consensus(
        args.regions,
        args.cells,
        args.output_prefix,
        min_cells=args.min_cells,
        min_fraction=args.min_fraction,
        per_sample=not args.all_samples,
        merge_cell_gap_bp=args.merge_cell_gap_bp,
        merge_consensus_gap_bp=args.merge_consensus_gap_bp,
    )
    print(pd.Series(summary).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
