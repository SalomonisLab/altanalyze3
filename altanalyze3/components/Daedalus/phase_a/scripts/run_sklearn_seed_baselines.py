#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"
IN_PATH = PROCESSED / "baseline_feature_matrix.tsv"
METRICS_OUT = PROCESSED / "sklearn_seed_baseline_metrics.tsv"
REPORT_OUT = PROCESSED / "sklearn_seed_baseline_report.md"

TASKS = [
    "task_global_seed",
    "task_membrane_seed",
    "task_kinase_seed",
    "task_tf_seed",
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


def _safe_metric(metric_fn, y_true, y_score):
    try:
        return float(metric_fn(y_true, y_score))
    except Exception:
        return float("nan")


def _metrics_dict(model: str, task: str, split: str, y_true, y_prob, y_pred) -> dict[str, str]:
    return {
        "model": model,
        "task": task,
        "split": split,
        "n": str(len(y_true)),
        "positives": str(int(sum(y_true))),
        "negatives": str(int(len(y_true) - sum(y_true))),
        "auroc": f"{_safe_metric(roc_auc_score, y_true, y_prob):.6f}",
        "average_precision": f"{_safe_metric(average_precision_score, y_true, y_prob):.6f}",
        "brier": f"{_safe_metric(brier_score_loss, y_true, y_prob):.6f}",
        "accuracy": f"{accuracy_score(y_true, y_pred):.6f}",
        "balanced_accuracy": f"{balanced_accuracy_score(y_true, y_pred):.6f}",
    }


def _build_models() -> dict[str, Pipeline | HistGradientBoostingClassifier]:
    return {
        "logistic_structural": Pipeline(
            [
                ("impute", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scale", StandardScaler()),
                ("model", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=0)),
            ]
        ),
        "mlp_structural": Pipeline(
            [
                ("impute", SimpleImputer(strategy="constant", fill_value=0.0)),
                ("scale", StandardScaler()),
                (
                    "model",
                    MLPClassifier(
                        hidden_layer_sizes=(64, 32),
                        activation="relu",
                        alpha=1e-4,
                        batch_size=256,
                        learning_rate_init=1e-3,
                        max_iter=200,
                        early_stopping=True,
                        random_state=0,
                    ),
                ),
            ]
        ),
        "random_forest_structural": Pipeline(
            [
                ("impute", SimpleImputer(strategy="constant", fill_value=0.0)),
                (
                    "model",
                    RandomForestClassifier(
                        n_estimators=300,
                        max_depth=None,
                        min_samples_leaf=2,
                        class_weight="balanced_subsample",
                        n_jobs=-1,
                        random_state=0,
                    ),
                ),
            ]
        ),
        "hist_gbdt_structural": HistGradientBoostingClassifier(
            learning_rate=0.05,
            max_depth=6,
            max_iter=300,
            min_samples_leaf=40,
            random_state=0,
        ),
    }


def main() -> int:
    df = pd.read_csv(IN_PATH, sep="\t")
    metric_rows: list[dict[str, str]] = []
    report_lines = ["# Scikit-learn Seed Baseline Report", "", "Structural-feature baselines trained on weakly supervised seed tasks.", ""]

    for task in TASKS:
        task_df = df[(df[task] == 1) & (df["label_preservation_seed"].isin([0, 1]))].copy()
        if task_df.empty:
            continue
        train_df = task_df[task_df["split"] == "train"]
        val_df = task_df[task_df["split"] == "val"]
        test_df = task_df[task_df["split"] == "test"]
        if train_df.empty or val_df.empty or test_df.empty:
            continue

        X_train = train_df[STRUCTURAL_FEATURE_COLUMNS]
        y_train = train_df["label_preservation_seed"].astype(int)
        X_val = val_df[STRUCTURAL_FEATURE_COLUMNS]
        y_val = val_df["label_preservation_seed"].astype(int)
        X_test = test_df[STRUCTURAL_FEATURE_COLUMNS]
        y_test = test_df["label_preservation_seed"].astype(int)

        report_lines.append(f"## {task}")
        report_lines.append("")
        for model_name, model in _build_models().items():
            model.fit(X_train, y_train)
            val_prob = model.predict_proba(X_val)[:, 1]
            val_pred = (val_prob >= 0.5).astype(int)
            test_prob = model.predict_proba(X_test)[:, 1]
            test_pred = (test_prob >= 0.5).astype(int)

            val_metrics = _metrics_dict(model_name, task, "val", y_val, val_prob, val_pred)
            test_metrics = _metrics_dict(model_name, task, "test", y_test, test_prob, test_pred)
            metric_rows.extend([val_metrics, test_metrics])
            report_lines.append(
                f"- `{model_name}` val AUROC={val_metrics['auroc']} AP={val_metrics['average_precision']} | "
                f"test AUROC={test_metrics['auroc']} AP={test_metrics['average_precision']} BalAcc={test_metrics['balanced_accuracy']}"
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
