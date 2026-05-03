#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from scipy import sparse
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import OneHotEncoder, StandardScaler


PHASE_D = Path(__file__).resolve().parents[1]
PROCESSED = PHASE_D / "data" / "processed"
CHECKPOINTS = PHASE_D / "checkpoints"

IN_PATH = PROCESSED / "phase_d_benchmark_instances.parquet"
SCHEMA_PATH = PROCESSED / "phase_d_benchmark_schema.json"
PHASED_TORCH = CHECKPOINTS / "phase_d_multitask_metrics.json"
OUT_TSV = CHECKPOINTS / "phase_d_model_comparison_metrics.tsv"
OUT_MD = CHECKPOINTS / "phase_d_model_comparison_report.md"


def _safe_metric(fn, y_true, y_score) -> float:
    try:
        return float(fn(y_true, y_score))
    except Exception:
        return float("nan")


def _metrics(y_true: np.ndarray, probs: np.ndarray) -> dict[str, float]:
    preds = (probs >= 0.5).astype(int)
    return {
        "n": int(len(y_true)),
        "positives": int(y_true.sum()),
        "negatives": int(len(y_true) - y_true.sum()),
        "average_precision": _safe_metric(average_precision_score, y_true, probs),
        "auroc": _safe_metric(roc_auc_score, y_true, probs),
        "balanced_accuracy": _safe_metric(balanced_accuracy_score, y_true, preds),
    }


def _prepare_sparse(train_df: pd.DataFrame, val_df: pd.DataFrame, test_df: pd.DataFrame, feature_columns: list[str]) -> tuple[sparse.csr_matrix, sparse.csr_matrix, sparse.csr_matrix]:
    cat_cols = ["task_family", "reference_source", "change_class"]
    scaler = StandardScaler()
    train_num = scaler.fit_transform(train_df[feature_columns].to_numpy(dtype=np.float32))
    val_num = scaler.transform(val_df[feature_columns].to_numpy(dtype=np.float32))
    test_num = scaler.transform(test_df[feature_columns].to_numpy(dtype=np.float32))
    encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=True)
    train_cat = encoder.fit_transform(train_df[cat_cols].astype(str))
    val_cat = encoder.transform(val_df[cat_cols].astype(str))
    test_cat = encoder.transform(test_df[cat_cols].astype(str))
    x_train = sparse.hstack([sparse.csr_matrix(train_num), train_cat], format="csr")
    x_val = sparse.hstack([sparse.csr_matrix(val_num), val_cat], format="csr")
    x_test = sparse.hstack([sparse.csr_matrix(test_num), test_cat], format="csr")
    return x_train, x_val, x_test


def _prepare_dense(train_df: pd.DataFrame, val_df: pd.DataFrame, test_df: pd.DataFrame, feature_columns: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cat_cols = ["task_family", "reference_source", "change_class"]
    scaler = StandardScaler()
    train_num = scaler.fit_transform(train_df[feature_columns].to_numpy(dtype=np.float32))
    val_num = scaler.transform(val_df[feature_columns].to_numpy(dtype=np.float32))
    test_num = scaler.transform(test_df[feature_columns].to_numpy(dtype=np.float32))
    encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    train_cat = encoder.fit_transform(train_df[cat_cols].astype(str))
    val_cat = encoder.transform(val_df[cat_cols].astype(str))
    test_cat = encoder.transform(test_df[cat_cols].astype(str))
    return (
        np.concatenate([train_num, train_cat], axis=1),
        np.concatenate([val_num, val_cat], axis=1),
        np.concatenate([test_num, test_cat], axis=1),
    )


def _run_models_for_task(df: pd.DataFrame, task: str, feature_columns: list[str]) -> list[dict[str, object]]:
    sub = df[(df["task"] == task) & (df["label_state"].isin(["positive", "negative"]))].copy()
    train_df = sub[sub["split"] == "train"].copy()
    val_df = sub[sub["split"] == "val"].copy()
    test_df = sub[sub["split"] == "test"].copy()
    if (
        train_df.empty
        or val_df.empty
        or test_df.empty
        or train_df["label_binary"].nunique() < 2
        or val_df["label_binary"].nunique() < 2
        or test_df["label_binary"].nunique() < 2
    ):
        return []
    x_train_sparse, x_val_sparse, x_test_sparse = _prepare_sparse(train_df, val_df, test_df, feature_columns)
    x_train_dense, x_val_dense, x_test_dense = _prepare_dense(train_df, val_df, test_df, feature_columns)
    y_train = train_df["label_binary"].to_numpy(dtype=np.int32)
    y_val = val_df["label_binary"].to_numpy(dtype=np.int32)
    y_test = test_df["label_binary"].to_numpy(dtype=np.int32)

    rows: list[dict[str, object]] = []

    logistic = LogisticRegression(
        solver="saga",
        penalty="l2",
        C=1.0,
        max_iter=250,
        n_jobs=-1,
        class_weight="balanced",
        random_state=0,
    )
    logistic.fit(x_train_sparse, y_train)
    probs = logistic.predict_proba(x_test_sparse)[:, 1]
    rows.append({"model": "phase_d_logistic", "task": task, **_metrics(y_test, probs)})

    hgb = HistGradientBoostingClassifier(
        max_depth=6,
        learning_rate=0.05,
        max_iter=500,
        early_stopping=True,
        validation_fraction=0.15,
        random_state=0,
    )
    hgb.fit(x_train_dense, y_train)
    probs = hgb.predict_proba(x_test_dense)[:, 1]
    rows.append({"model": "phase_d_hist_gradient_boosting", "task": task, **_metrics(y_test, probs)})

    xgb_model = xgb.XGBClassifier(
        n_estimators=500,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        reg_lambda=1.0,
        min_child_weight=1.0,
        objective="binary:logistic",
        eval_metric="aucpr",
        tree_method="hist",
        random_state=0,
        n_jobs=8,
        early_stopping_rounds=30,
    )
    xgb_model.fit(x_train_sparse, y_train, eval_set=[(x_val_sparse, y_val)], verbose=False)
    probs = xgb_model.predict_proba(x_test_sparse)[:, 1]
    rows.append({"model": "phase_d_xgboost", "task": task, **_metrics(y_test, probs)})
    return rows


def main() -> int:
    df = pd.read_parquet(IN_PATH)
    with SCHEMA_PATH.open() as handle:
        schema = json.load(handle)
    feature_columns = schema["feature_columns"]

    rows: list[dict[str, object]] = []
    for task in sorted(df["task"].astype(str).unique().tolist()):
        rows.extend(_run_models_for_task(df, task, feature_columns))

    phase_d_torch = json.loads(PHASED_TORCH.read_text())["test"]
    for task, metrics in phase_d_torch.items():
        row = {"model": "phase_d_family_multitask", "task": task}
        row.update(metrics)
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(["task", "model"]).reset_index(drop=True)
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    lines = [
        "# Phase D Model Comparison",
        "",
        "These are bootstrap proxy benchmarks using labels derived from Phase D objective thresholds.",
        "",
        "Included models:",
        "- `phase_d_logistic` sanity baseline",
        "- `phase_d_hist_gradient_boosting` strong tree baseline",
        "- `phase_d_xgboost` strong tree baseline",
        "- `phase_d_family_multitask` family-aware deep multitask model",
        "",
        "## Test metrics by task",
        "",
    ]
    for task in sorted(out["task"].unique().tolist()):
        lines.append(f"### {task}")
        lines.append("")
        sub = out[out["task"] == task].copy().sort_values("average_precision", ascending=False)
        for row in sub.itertuples(index=False):
            lines.append(
                f"- `{row.model}` AP={row.average_precision:.6f} "
                f"AUROC={row.auroc:.6f} BalAcc={row.balanced_accuracy:.6f}"
            )
        lines.append("")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[ok] Wrote {OUT_TSV}")
    print(f"[ok] Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
