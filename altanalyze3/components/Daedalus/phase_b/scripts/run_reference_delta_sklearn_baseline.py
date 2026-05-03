#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"
CHECKPOINTS = Path(__file__).resolve().parents[1] / "checkpoints"
IN_PATH = PROCESSED / "multitask_instances.parquet"
SCHEMA_PATH = PROCESSED / "multitask_schema.json"
METRICS_TSV = PROCESSED / "reference_delta_sklearn_metrics.tsv"
METRICS_JSON = PROCESSED / "reference_delta_sklearn_metrics.json"
REPORT_PATH = CHECKPOINTS / "reference_delta_sklearn_report.md"


def _safe_metric(metric_fn, y_true, y_score) -> float:
    try:
        return float(metric_fn(y_true, y_score))
    except Exception:
        return float("nan")


def _task_metrics(frame: pd.DataFrame, probs) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    preds = (probs >= 0.5).astype(int)
    for task in sorted(frame["task"].astype(str).unique().tolist()):
        sub = frame[frame["task"].astype(str) == task].copy()
        idx = sub.index.to_numpy()
        y = sub["label_preservation_seed"].astype(int).to_numpy()
        p = probs[idx]
        yhat = preds[idx]
        rows.append(
            {
                "task": task,
                "n": int(len(sub)),
                "positives": int(y.sum()),
                "negatives": int(len(sub) - y.sum()),
                "average_precision": _safe_metric(average_precision_score, y, p),
                "auroc": _safe_metric(roc_auc_score, y, p),
                "balanced_accuracy": _safe_metric(balanced_accuracy_score, y, yhat),
            }
        )
    return rows


def main() -> int:
    with SCHEMA_PATH.open() as handle:
        schema = json.load(handle)

    df = pd.read_parquet(IN_PATH)
    numeric_columns = schema["numeric_columns"]
    categorical_columns = schema["categorical_columns"]
    feature_columns = [*numeric_columns, *categorical_columns]

    train = df[df["split"] == "train"].copy().reset_index(drop=True)
    test = df[df["split"] == "test"].copy().reset_index(drop=True)

    preprocessor = ColumnTransformer(
        [
            (
                "num",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="constant", fill_value=0.0)),
                        ("scale", StandardScaler()),
                    ]
                ),
                numeric_columns,
            ),
            (
                "cat",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="most_frequent")),
                        ("oh", OneHotEncoder(handle_unknown="ignore")),
                    ]
                ),
                categorical_columns,
            ),
        ]
    )
    model = Pipeline(
        [
            ("pre", preprocessor),
            ("clf", LogisticRegression(max_iter=500, class_weight="balanced")),
        ]
    )

    y_train = train["label_preservation_seed"].astype(int)
    y_test = test["label_preservation_seed"].astype(int)
    model.fit(train[feature_columns], y_train)
    probs = model.predict_proba(test[feature_columns])[:, 1]

    overall = {
        "task": "overall",
        "n": int(len(test)),
        "positives": int(y_test.sum()),
        "negatives": int(len(test) - y_test.sum()),
        "average_precision": _safe_metric(average_precision_score, y_test, probs),
        "auroc": _safe_metric(roc_auc_score, y_test, probs),
        "balanced_accuracy": _safe_metric(balanced_accuracy_score, y_test, (probs >= 0.5).astype(int)),
    }
    metrics = [overall, *_task_metrics(test, probs)]
    metrics_df = pd.DataFrame(metrics)

    PROCESSED.mkdir(parents=True, exist_ok=True)
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(METRICS_TSV, sep="\t", index=False)
    with METRICS_JSON.open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)

    lines = [
        "# Reference Delta sklearn Baseline",
        "",
        "Leakage-resistant logistic baseline on the Phase B multitask table.",
        "",
        "## Test metrics",
        "",
    ]
    for row in metrics:
        lines.append(
            f"- `{row['task']}` n={row['n']} AP={row['average_precision']:.6f} "
            f"AUROC={row['auroc']:.6f} BalAcc={row['balanced_accuracy']:.6f}"
        )
    with REPORT_PATH.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")

    print(metrics_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
