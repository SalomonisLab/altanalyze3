"""Evaluation metrics for the experimental fastCNV candidate.

This module evaluates outputs; it does not change CNV calls. The three metric
families map to the user acceptance requirements:

1. LOY accuracy against the existing prior LOY labels.
2. Control false-positive rate on evaluated healthy/control cells.
3. Targeted AML cytogenetic detection, where only known target lesions are
   scored as required positives and other AML CNVs are not counted as false
   positives.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class MetricResult:
    metric: str
    value: float
    threshold: str
    passed: bool
    detail: str


def _read_cells(prefix: Path) -> pd.DataFrame:
    path = Path(f"{prefix}.candidate_cnv_cells.tsv")
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    if "CellBarcode" not in frame.columns or "cnv_status" not in frame.columns:
        raise ValueError(f"{path} is missing required columns.")
    frame["CellBarcode"] = frame["CellBarcode"].astype(str)
    return frame


def _read_regions(prefix: Path) -> pd.DataFrame:
    path = Path(f"{prefix}.candidate_cnv_regions.tsv")
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    if frame.empty:
        return frame
    frame["CellBarcode"] = frame["CellBarcode"].astype(str)
    frame["chr"] = frame["chr"].astype(str)
    frame["call"] = frame["call"].astype(str)
    return frame


def _safe_div(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator else float("nan")


def _bool_pass(value: float, op: str, threshold: float) -> bool:
    if not np.isfinite(value):
        return False
    if op == "<=":
        return value <= threshold
    if op == ">=":
        return value >= threshold
    raise ValueError(op)


def evaluate_loy(
    prefix: Path,
    truth_path: Path,
    *,
    min_expected_calls: int = 13000,
    max_expected_calls: int = 14000,
    min_state_cells: int = 20,
    min_state_sensitivity: float = 0.95,
    min_state_specificity: float = 0.95,
    min_state_pass_fraction: float = 0.90,
) -> tuple[list[MetricResult], pd.DataFrame]:
    cells = _read_cells(prefix)
    truth = pd.read_csv(truth_path, sep="\t", usecols=["CellBarcode", "LOY", "cell_state"])
    truth["CellBarcode"] = truth["CellBarcode"].astype(str)
    truth["LOY"] = truth["LOY"].astype(int)
    merged = truth.merge(cells[["CellBarcode", "cnv_status"]], on="CellBarcode", how="left")
    called = merged["cnv_status"].eq("CNV")
    truth_positive = merged["LOY"].eq(1)

    tp = int((called & truth_positive).sum())
    fp = int((called & ~truth_positive).sum())
    fn = int((~called & truth_positive).sum())
    tn = int((~called & ~truth_positive).sum())
    n_calls = int(called.sum())
    sensitivity = _safe_div(tp, tp + fn)
    specificity = _safe_div(tn, tn + fp)
    precision = _safe_div(tp, tp + fp)

    per_state_rows = []
    for state, group in merged.groupby("cell_state", observed=False):
        y = group["LOY"].eq(1)
        c = group["cnv_status"].eq("CNV")
        n_pos = int(y.sum())
        n_neg = int((~y).sum())
        if n_pos < min_state_cells or n_neg < min_state_cells:
            continue
        sens = _safe_div(int((c & y).sum()), n_pos)
        spec = _safe_div(int((~c & ~y).sum()), n_neg)
        per_state_rows.append(
            {
                "cell_state": state,
                "n_truth_loy": n_pos,
                "n_truth_wt": n_neg,
                "n_called": int(c.sum()),
                "sensitivity": sens,
                "specificity": spec,
                "passes_both": bool(sens >= min_state_sensitivity and spec >= min_state_specificity),
            }
        )
    per_state = pd.DataFrame(per_state_rows).sort_values(
        ["passes_both", "sensitivity", "specificity"], ascending=[True, True, True]
    )
    state_pass_fraction = (
        float(per_state["passes_both"].mean()) if not per_state.empty else float("nan")
    )

    results = [
        MetricResult(
            "LOY_call_count",
            float(n_calls),
            f"{min_expected_calls} <= calls <= {max_expected_calls}",
            min_expected_calls <= n_calls <= max_expected_calls,
            f"called={n_calls}; truth_LOY={int(truth_positive.sum())}",
        ),
        MetricResult(
            "LOY_overall_sensitivity",
            sensitivity,
            ">= 0.95",
            _bool_pass(sensitivity, ">=", 0.95),
            f"TP={tp}; FN={fn}",
        ),
        MetricResult(
            "LOY_overall_specificity",
            specificity,
            ">= 0.95",
            _bool_pass(specificity, ">=", 0.95),
            f"TN={tn}; FP={fp}",
        ),
        MetricResult(
            "LOY_precision_vs_prior",
            precision,
            "reported",
            bool(np.isfinite(precision)),
            f"TP={tp}; FP={fp}",
        ),
        MetricResult(
            "LOY_state_fraction_sens_spec_ge95",
            state_pass_fraction,
            f">= {min_state_pass_fraction}",
            _bool_pass(state_pass_fraction, ">=", min_state_pass_fraction),
            f"passing_states={int(per_state['passes_both'].sum()) if not per_state.empty else 0}; evaluated_states={per_state.shape[0]}",
        ),
    ]
    return results, per_state


def evaluate_control_false_positive(
    prefix: Path,
    *,
    max_cnv_fraction: float = 0.05,
) -> tuple[list[MetricResult], pd.DataFrame, pd.DataFrame]:
    cells = _read_cells(prefix)
    n_cells = int(cells.shape[0])
    n_cnv = int(cells["cnv_status"].eq("CNV").sum())
    cnv_fraction = _safe_div(n_cnv, n_cells)
    sample = cells.get("sample", pd.Series(["sample"] * n_cells)).astype(str)
    by_sample = (
        cells.assign(sample=sample, is_cnv=cells["cnv_status"].eq("CNV"))
        .groupby("sample", observed=False)
        .agg(n_cells=("CellBarcode", "size"), n_cnv_cells=("is_cnv", "sum"))
        .reset_index()
    )
    by_sample["cnv_fraction"] = by_sample["n_cnv_cells"] / by_sample["n_cells"].replace(0, np.nan)
    by_sample["passes_lt5pct"] = by_sample["cnv_fraction"] <= max_cnv_fraction
    by_state_sample = (
        cells.assign(sample=sample, is_cnv=cells["cnv_status"].eq("CNV"))
        .groupby(["sample", "cell_state"], observed=False)
        .agg(n_cells=("CellBarcode", "size"), n_cnv_cells=("is_cnv", "sum"))
        .reset_index()
    )
    by_state_sample["cnv_fraction"] = by_state_sample["n_cnv_cells"] / by_state_sample["n_cells"].replace(0, np.nan)
    by_state_sample["passes_lt5pct"] = by_state_sample["cnv_fraction"] <= max_cnv_fraction
    results = [
        MetricResult(
            "control_cnv_fraction",
            cnv_fraction,
            f"<= {max_cnv_fraction}",
            _bool_pass(cnv_fraction, "<=", max_cnv_fraction),
            f"CNV={n_cnv}; evaluated_control_cells={n_cells}",
        ),
        MetricResult(
            "control_samples_lt5pct_fraction",
            float(by_sample["passes_lt5pct"].mean()) if not by_sample.empty else float("nan"),
            "all samples should pass",
            bool(by_sample["passes_lt5pct"].all()) if not by_sample.empty else False,
            f"passing_samples={int(by_sample['passes_lt5pct'].sum())}; evaluated_samples={by_sample.shape[0]}",
        ),
    ]
    return results, by_sample.sort_values("cnv_fraction", ascending=False), by_state_sample.sort_values("cnv_fraction", ascending=False)


def build_state_clone_artifact(prefix: Path) -> pd.DataFrame:
    cells = _read_cells(prefix)
    regions = _read_regions(prefix)
    if regions.empty:
        return pd.DataFrame(
            columns=[
                "cell_state", "sample", "state_clone_id", "n_cells", "fraction_of_state",
                "call", "chr", "region", "mean_log2", "max_abs_score",
                "global_clone_id", "clone_confidence",
            ]
        )
    totals = cells.groupby(["sample", "cell_state"], observed=False).size().rename("n_state_cells").reset_index()
    frame = regions.copy()
    frame["region"] = frame["chr"].astype(str) + ":" + frame["start"].astype(str) + "-" + frame["end"].astype(str)
    grouped = (
        frame.groupby(["sample", "cell_state", "call", "chr", "region"], observed=False)
        .agg(
            n_cells=("CellBarcode", "nunique"),
            mean_log2=("mean_log2", "mean"),
            max_abs_score=("max_abs_log2", "max"),
            mean_confidence_z=("mean_confidence_z", "mean"),
        )
        .reset_index()
        .merge(totals, on=["sample", "cell_state"], how="left")
    )
    grouped["fraction_of_state"] = grouped["n_cells"] / grouped["n_state_cells"].replace(0, np.nan)
    grouped = grouped.sort_values(["sample", "cell_state", "fraction_of_state", "n_cells"], ascending=[True, True, False, False])
    grouped["state_clone_id"] = (
        grouped["cell_state"].astype(str)
        + "::"
        + grouped["sample"].astype(str)
        + "::"
        + grouped.groupby(["sample", "cell_state"], observed=False).cumcount().add(1).astype(str).radd("clone")
    )
    grouped["global_clone_id"] = (
        grouped["chr"].astype(str) + ":" + grouped["call"].astype(str) + ":" + grouped["region"].astype(str)
    )
    grouped["clone_confidence"] = np.where(grouped["fraction_of_state"] >= 0.40, "high", np.where(grouped["fraction_of_state"] >= 0.10, "medium", "low"))
    return grouped[
        [
            "cell_state", "sample", "state_clone_id", "n_cells", "fraction_of_state",
            "call", "chr", "region", "mean_log2", "max_abs_score", "mean_confidence_z",
            "global_clone_id", "clone_confidence",
        ]
    ]


def _target_region_match(regions: pd.DataFrame, chrom: str, direction: str, region: str) -> pd.Series:
    if regions.empty:
        return pd.Series([], dtype=bool)
    direction = str(direction).lower()
    call_ok = regions["call"].isin(["loss", "homozygous_loss"]) if direction == "loss" else regions["call"].isin(["gain", "amplification"])
    chrom_ok = regions["chr"].astype(str).eq(str(chrom))
    if not region or str(region) in {"-", "whole", "nan"}:
        return chrom_ok & call_ok
    region_value = str(region)
    gene_text = (
        regions.get("start_gene", pd.Series([""] * regions.shape[0])).astype(str)
        + ";"
        + regions.get("end_gene", pd.Series([""] * regions.shape[0])).astype(str)
    )
    # Cytobands are not in candidate output yet, so regional del(5q) style
    # truth is scored at chromosome+direction level unless cytoband columns are
    # later added to the caller output.
    return chrom_ok & call_ok & gene_text.notna()


def evaluate_aml_targets(
    output_dir: Path,
    truth_table: Path,
    *,
    min_target_cell_fraction: float = 0.40,
) -> tuple[list[MetricResult], pd.DataFrame]:
    truth = pd.read_csv(truth_table)
    truth = truth.loc[truth["detectable"].astype(str).eq("yes")].copy()
    rows = []
    for _, truth_row in truth.iterrows():
        sample = str(truth_row["Sample"])
        prefix = output_dir / sample
        cells_path = Path(f"{prefix}.candidate_cnv_cells.tsv")
        regions_path = Path(f"{prefix}.candidate_cnv_regions.tsv")
        if not cells_path.exists() or not regions_path.exists():
            rows.append(
                {
                    "Sample": sample,
                    "lesion": truth_row["lesion"],
                    "chr": truth_row["chr"],
                    "direction": truth_row["direction"],
                    "region": truth_row["region"],
                    "n_cells": 0,
                    "n_target_cells": 0,
                    "target_cell_fraction": np.nan,
                    "detected": False,
                    "detail": "missing output",
                }
            )
            continue
        cells = _read_cells(prefix)
        regions = _read_regions(prefix)
        matches = _target_region_match(
            regions,
            str(truth_row["chr"]),
            str(truth_row["direction"]),
            str(truth_row["region"]),
        )
        target_cells = set(regions.loc[matches, "CellBarcode"].astype(str)) if len(matches) else set()
        n_cells = int(cells.shape[0])
        n_target = int(cells["CellBarcode"].isin(target_cells).sum())
        frac = _safe_div(n_target, n_cells)
        rows.append(
            {
                "Sample": sample,
                "lesion": truth_row["lesion"],
                "chr": truth_row["chr"],
                "direction": truth_row["direction"],
                "region": truth_row["region"],
                "n_cells": n_cells,
                "n_target_cells": n_target,
                "target_cell_fraction": frac,
                "detected": bool(frac >= min_target_cell_fraction),
                "detail": "other AML CNVs ignored for this metric",
            }
        )
    detail = pd.DataFrame(rows)
    detected = int(detail["detected"].sum()) if not detail.empty else 0
    total = int(detail.shape[0])
    results = [
        MetricResult(
            "AML_target_lesion_detection_fraction",
            _safe_div(detected, total),
            f"each detectable lesion >= {min_target_cell_fraction} sample cells",
            detected == total and total > 0,
            f"detected={detected}; evaluated_detectable_lesions={total}",
        )
    ]
    return results, detail


def _results_to_frame(results: Iterable[MetricResult]) -> pd.DataFrame:
    return pd.DataFrame([result.__dict__ for result in results])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Evaluate candidate fastCNV outputs.")
    parser.add_argument("--loy-prefix", type=Path)
    parser.add_argument("--loy-truth", type=Path)
    parser.add_argument("--control-prefix", type=Path)
    parser.add_argument("--aml-output-dir", type=Path)
    parser.add_argument("--aml-truth-table", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    all_results: list[MetricResult] = []
    if args.loy_prefix and args.loy_truth:
        results, per_state = evaluate_loy(args.loy_prefix, args.loy_truth)
        all_results.extend(results)
        per_state.to_csv(f"{args.output_prefix}.loy_per_state.tsv", sep="\t", index=False)
    if args.control_prefix:
        results, by_sample, by_state_sample = evaluate_control_false_positive(args.control_prefix)
        all_results.extend(results)
        by_sample.to_csv(f"{args.output_prefix}.control_by_sample.tsv", sep="\t", index=False)
        by_state_sample.to_csv(f"{args.output_prefix}.control_by_sample_state.tsv", sep="\t", index=False)
        build_state_clone_artifact(args.control_prefix).to_csv(
            f"{args.output_prefix}.control_state_clones.tsv", sep="\t", index=False
        )
    if args.aml_output_dir and args.aml_truth_table:
        results, detail = evaluate_aml_targets(args.aml_output_dir, args.aml_truth_table)
        all_results.extend(results)
        detail.to_csv(f"{args.output_prefix}.aml_targets.tsv", sep="\t", index=False)
        clone_frames = []
        for prefix_path in sorted(args.aml_output_dir.glob("*.candidate_cnv_cells.tsv")):
            prefix = Path(str(prefix_path).removesuffix(".candidate_cnv_cells.tsv"))
            frame = build_state_clone_artifact(prefix)
            if not frame.empty:
                clone_frames.append(frame)
            frame.to_csv(f"{prefix}.candidate_state_clones.tsv", sep="\t", index=False)
        if clone_frames:
            pd.concat(clone_frames, ignore_index=True).to_csv(
                f"{args.output_prefix}.aml_state_clones.tsv", sep="\t", index=False
            )
    summary = _results_to_frame(all_results)
    summary.to_csv(f"{args.output_prefix}.summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
