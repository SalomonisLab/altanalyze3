#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import random
from dataclasses import dataclass
from pathlib import Path


PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"
IN_PATH = PROCESSED / "baseline_feature_matrix.tsv"
METRICS_OUT = PROCESSED / "seed_baseline_metrics.tsv"
REPORT_OUT = PROCESSED / "seed_baseline_report.md"

TASKS = [
    "task_global_seed",
    "task_membrane_seed",
    "task_kinase_seed",
    "task_tf_seed",
    "task_high_confidence_seed",
]

FEATURE_COLUMNS = [
    "feature_protein_ratio",
    "feature_protein_ratio_distance",
    "feature_protein_delta_norm",
    "feature_transcript_delta_norm",
    "feature_reference_uniprot_links",
    "feature_alternative_uniprot_links",
    "feature_alt_has_protein",
    "feature_ref_alt_same_protein_length",
    "feature_ref_is_membrane",
    "feature_alt_is_membrane",
    "feature_membrane_mismatch",
    "feature_ref_has_signal",
    "feature_alt_has_signal",
    "feature_signal_mismatch",
    "feature_ref_is_kinase",
    "feature_alt_is_kinase",
    "feature_kinase_mismatch",
    "feature_ref_is_tf",
    "feature_alt_is_tf",
    "feature_tf_mismatch",
    "feature_gene_has_appris_principal_1",
    "feature_gene_clinvar_splice",
    "feature_gene_hpa_extracellular",
    "feature_is_high_confidence",
    "feature_is_membrane_pair",
    "feature_is_surface_pair",
    "feature_is_kinase_pair",
    "feature_is_tf_pair",
    "feature_reference_source_mane",
    "feature_reference_source_appris",
    "feature_reference_source_uniprot_longest",
    "feature_reference_source_longest_protein",
    "feature_annotation_retention",
]

STRUCTURAL_FEATURE_COLUMNS = [
    "feature_protein_ratio",
    "feature_protein_ratio_distance",
    "feature_protein_delta_norm",
    "feature_transcript_delta_norm",
    "feature_alt_has_protein",
    "feature_ref_alt_same_protein_length",
    "feature_gene_has_appris_principal_1",
    "feature_gene_clinvar_splice",
    "feature_gene_hpa_extracellular",
    "feature_is_membrane_pair",
    "feature_is_surface_pair",
    "feature_is_kinase_pair",
    "feature_is_tf_pair",
    "feature_reference_source_mane",
    "feature_reference_source_appris",
    "feature_reference_source_uniprot_longest",
    "feature_reference_source_longest_protein",
]


def _float(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0


def _sigmoid(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def _roc_auc(labels: list[int], scores: list[float]) -> float:
    pos = sum(labels)
    neg = len(labels) - pos
    if pos == 0 or neg == 0:
        return float("nan")
    ranked = sorted(zip(scores, labels), key=lambda x: x[0])
    rank_sum = 0.0
    i = 0
    rank = 1
    while i < len(ranked):
        j = i
        while j < len(ranked) and ranked[j][0] == ranked[i][0]:
            j += 1
        avg_rank = (rank + (rank + (j - i) - 1)) / 2.0
        positives_in_tie = sum(label for _, label in ranked[i:j])
        rank_sum += positives_in_tie * avg_rank
        rank += j - i
        i = j
    return (rank_sum - pos * (pos + 1) / 2.0) / (pos * neg)


def _average_precision(labels: list[int], scores: list[float]) -> float:
    pos = sum(labels)
    if pos == 0:
        return float("nan")
    ranked = sorted(zip(scores, labels), key=lambda x: x[0], reverse=True)
    tp = 0
    fp = 0
    ap = 0.0
    for _, label in ranked:
        if label == 1:
            tp += 1
            ap += tp / (tp + fp)
        else:
            fp += 1
    return ap / pos


def _brier(labels: list[int], scores: list[float]) -> float:
    return sum((score - label) ** 2 for label, score in zip(labels, scores)) / max(len(labels), 1)


def _threshold_metrics(labels: list[int], scores: list[float], threshold: float = 0.5) -> tuple[float, float]:
    tp = tn = fp = fn = 0
    for label, score in zip(labels, scores):
        pred = 1 if score >= threshold else 0
        if pred == 1 and label == 1:
            tp += 1
        elif pred == 0 and label == 0:
            tn += 1
        elif pred == 1 and label == 0:
            fp += 1
        else:
            fn += 1
    acc = (tp + tn) / max(len(labels), 1)
    tpr = tp / max(tp + fn, 1)
    tnr = tn / max(tn + fp, 1)
    bal_acc = 0.5 * (tpr + tnr)
    return acc, bal_acc


def _length_similarity(row: dict[str, str]) -> float:
    ratio = max(min(_float(row["feature_protein_ratio"]), 4.0), 0.25)
    return math.exp(-abs(math.log(ratio)))


def _structural_heuristic(row: dict[str, str]) -> float:
    return (
        0.55 * _length_similarity(row)
        + 0.25 * _float(row["feature_ref_alt_same_protein_length"])
        + 0.20 * _float(row["feature_alt_has_protein"])
    )


@dataclass
class LogisticModel:
    weights: list[float]
    bias: float
    means: list[float]
    stds: list[float]

    def score(self, features: list[float]) -> float:
        z = self.bias
        for w, x, mean, std in zip(self.weights, features, self.means, self.stds):
            z += w * ((x - mean) / std)
        return _sigmoid(z)


def _train_logistic(rows: list[dict[str, str]], feature_columns: list[str]) -> LogisticModel:
    data = [[_float(row[col]) for col in feature_columns] for row in rows]
    labels = [int(row["label_preservation_seed"]) for row in rows]
    n_features = len(feature_columns)
    means = []
    stds = []
    for j in range(n_features):
        col = [vec[j] for vec in data]
        mean = sum(col) / max(len(col), 1)
        var = sum((x - mean) ** 2 for x in col) / max(len(col), 1)
        std = math.sqrt(var) if var > 1e-12 else 1.0
        means.append(mean)
        stds.append(std)

    weights = [0.0] * n_features
    bias = 0.0
    pos = sum(labels)
    neg = len(labels) - pos
    pos_weight = (neg / pos) if pos else 1.0
    rng = random.Random(0)
    order = list(range(len(rows)))
    lr = 0.05
    reg = 1e-4
    for _epoch in range(8):
        rng.shuffle(order)
        for idx in order:
            x = [((data[idx][j] - means[j]) / stds[j]) for j in range(n_features)]
            y = labels[idx]
            p = _sigmoid(bias + sum(w * v for w, v in zip(weights, x)))
            error = p - y
            scale = pos_weight if y == 1 else 1.0
            bias -= lr * scale * error
            for j in range(n_features):
                weights[j] -= lr * (scale * error * x[j] + reg * weights[j])
    return LogisticModel(weights=weights, bias=bias, means=means, stds=stds)


def _evaluate(model_name: str, task: str, split: str, rows: list[dict[str, str]], scorer) -> dict[str, str]:
    labels = [int(row["label_preservation_seed"]) for row in rows]
    scores = [scorer(row) for row in rows]
    acc, bal_acc = _threshold_metrics(labels, scores)
    return {
        "model": model_name,
        "task": task,
        "split": split,
        "n": str(len(rows)),
        "positives": str(sum(labels)),
        "negatives": str(len(rows) - sum(labels)),
        "auroc": f"{_roc_auc(labels, scores):.6f}",
        "average_precision": f"{_average_precision(labels, scores):.6f}",
        "brier": f"{_brier(labels, scores):.6f}",
        "accuracy": f"{acc:.6f}",
        "balanced_accuracy": f"{bal_acc:.6f}",
    }


def main() -> int:
    rows = []
    with IN_PATH.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows.extend(reader)

    metric_rows: list[dict[str, str]] = []
    report_lines = ["# Seed Baseline Report", "", "Pure-Python feasibility baselines on weakly supervised preservation tasks.", ""]

    for task in TASKS:
        task_rows = [row for row in rows if row[task] == "1" and row["label_preservation_seed"] in ("0", "1")]
        train_rows = [row for row in task_rows if row["split"] == "train"]
        val_rows = [row for row in task_rows if row["split"] == "val"]
        test_rows = [row for row in task_rows if row["split"] == "test"]
        if not train_rows or not val_rows or not test_rows:
            continue

        model_all = _train_logistic(train_rows, FEATURE_COLUMNS)
        model_struct = _train_logistic(train_rows, STRUCTURAL_FEATURE_COLUMNS)
        baselines = {
            "length_similarity": _length_similarity,
            "structural_heuristic": _structural_heuristic,
            "logistic_structural": lambda row, m=model_struct: m.score([_float(row[col]) for col in STRUCTURAL_FEATURE_COLUMNS]),
            "logistic_all_features": lambda row, m=model_all: m.score([_float(row[col]) for col in FEATURE_COLUMNS]),
        }

        report_lines.append(f"## {task}")
        report_lines.append("")
        for model_name, scorer in baselines.items():
            val_metrics = _evaluate(model_name, task, "val", val_rows, scorer)
            test_metrics = _evaluate(model_name, task, "test", test_rows, scorer)
            metric_rows.extend([val_metrics, test_metrics])
            report_lines.append(
                f"- `{model_name}` val AUROC={val_metrics['auroc']} AP={val_metrics['average_precision']} | "
                f"test AUROC={test_metrics['auroc']} AP={test_metrics['average_precision']}"
            )
        report_lines.append("")

    with METRICS_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "model",
                "task",
                "split",
                "n",
                "positives",
                "negatives",
                "auroc",
                "average_precision",
                "brier",
                "accuracy",
                "balanced_accuracy",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(metric_rows)

    with REPORT_OUT.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(report_lines) + "\n")

    print(f"[ok] Wrote {METRICS_OUT} rows={len(metric_rows)}")
    print(f"[ok] Wrote {REPORT_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
