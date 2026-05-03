#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import OneHotEncoder, StandardScaler


ROOT = Path(__file__).resolve().parents[3]
PHASE_B = ROOT / "Daedalus" / "phase_b" / "checkpoints"
PHASE_C = Path(__file__).resolve().parents[1]
PROCESSED = PHASE_C / "data" / "processed"
CHECKPOINTS = PHASE_C / "checkpoints"

IN_PATH = PROCESSED / "biochem_multitask_instances.parquet"
SCHEMA_PATH = PROCESSED / "biochem_multitask_schema.json"
PHASEB_TORCH = PHASE_B / "reference_delta_multitask_metrics.json"
PHASEB_SKLEARN = PHASE_B / "reference_delta_sklearn_report.md"
PHASEC_TORCH = CHECKPOINTS / "biochem_multitask_metrics.json"

OUT_TSV = CHECKPOINTS / "model_comparison_metrics.tsv"
OUT_MD = CHECKPOINTS / "model_comparison_report.md"


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


def _load_phaseb_sklearn() -> dict[str, dict[str, float]]:
    pattern = re.compile(
        r"- `(?P<task>[^`]+)` n=(?P<n>\d+) AP=(?P<ap>[0-9.]+) AUROC=(?P<auroc>[0-9.]+) BalAcc=(?P<bal>[0-9.]+)"
    )
    out: dict[str, dict[str, float]] = {}
    text = PHASEB_SKLEARN.read_text()
    for match in pattern.finditer(text):
        task = match.group("task")
        if task == "overall":
            continue
        out[task] = {
            "n": int(match.group("n")),
            "average_precision": float(match.group("ap")),
            "auroc": float(match.group("auroc")),
            "balanced_accuracy": float(match.group("bal")),
        }
    return out


def _prepare_features(train_df: pd.DataFrame, val_df: pd.DataFrame, test_df: pd.DataFrame, numeric_columns: list[str]) -> tuple[sparse.csr_matrix, sparse.csr_matrix, sparse.csr_matrix]:
    cat_cols = ["family_class", "reference_source", "change_class"]
    scaler = StandardScaler()
    train_num = scaler.fit_transform(train_df[numeric_columns].to_numpy(dtype=np.float32))
    val_num = scaler.transform(val_df[numeric_columns].to_numpy(dtype=np.float32))
    test_num = scaler.transform(test_df[numeric_columns].to_numpy(dtype=np.float32))
    encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=True)
    train_cat = encoder.fit_transform(train_df[cat_cols].astype(str))
    val_cat = encoder.transform(val_df[cat_cols].astype(str))
    test_cat = encoder.transform(test_df[cat_cols].astype(str))
    x_train = sparse.hstack([sparse.csr_matrix(train_num), train_cat], format="csr")
    x_val = sparse.hstack([sparse.csr_matrix(val_num), val_cat], format="csr")
    x_test = sparse.hstack([sparse.csr_matrix(test_num), test_cat], format="csr")
    return x_train, x_val, x_test


def _run_models_for_task(df: pd.DataFrame, task: str, numeric_columns: list[str]) -> list[dict[str, object]]:
    sub = df[df["task"] == task].copy()
    train_df = sub[sub["split"] == "train"].copy()
    val_df = sub[sub["split"] == "val"].copy()
    test_df = sub[sub["split"] == "test"].copy()
    x_train, x_val, x_test = _prepare_features(train_df, val_df, test_df, numeric_columns)
    y_train = train_df["label_preservation_seed"].to_numpy(dtype=np.int32)
    y_val = val_df["label_preservation_seed"].to_numpy(dtype=np.int32)
    y_test = test_df["label_preservation_seed"].to_numpy(dtype=np.int32)

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
    logistic.fit(x_train, y_train)
    probs = logistic.predict_proba(x_test)[:, 1]
    row = {"model": "phase_c_logistic", "task": task}
    row.update(_metrics(y_test, probs))
    rows.append(row)

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
    xgb_model.fit(x_train, y_train, eval_set=[(x_val, y_val)], verbose=False)
    probs = xgb_model.predict_proba(x_test)[:, 1]
    row = {"model": "phase_c_xgboost", "task": task}
    row.update(_metrics(y_test, probs))
    rows.append(row)
    return rows


def main() -> int:
    df = pd.read_parquet(IN_PATH)
    with SCHEMA_PATH.open() as handle:
        schema = json.load(handle)
    numeric_columns = schema["numeric_columns"]

    rows: list[dict[str, object]] = []
    for task in sorted(df["task"].astype(str).unique().tolist()):
        rows.extend(_run_models_for_task(df, task, numeric_columns))

    phaseb_torch = json.loads(PHASEB_TORCH.read_text())["test"]
    for task, metrics in phaseb_torch.items():
        row = {"model": "phase_b_torch", "task": task}
        row.update(metrics)
        rows.append(row)

    phaseb_sklearn = _load_phaseb_sklearn()
    for task, metrics in phaseb_sklearn.items():
        row = {"model": "phase_b_logistic", "task": task}
        row.update(metrics)
        rows.append(row)

    phasec_torch = json.loads(PHASEC_TORCH.read_text())["test"]
    for task, metrics in phasec_torch.items():
        row = {"model": "phase_c_autoencoder", "task": task}
        row.update(metrics)
        rows.append(row)

    out = pd.DataFrame(rows).sort_values(["task", "model"]).reset_index(drop=True)
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    lines = [
        "# Phase C Model Comparison",
        "",
        "Current Phase C models are compared on the same Phase C held-out subset.",
        "Phase B rows are historical reference metrics from the broader Phase B corpus and are included for directional context only.",
        "",
        "Included models:",
        "- Phase B leakage-resistant logistic historical reference",
        "- Phase B multitask torch historical reference",
        "- Phase C rich-feature logistic baseline",
        "- Phase C rich-feature XGBoost benchmark",
        "- Phase C class-specific biochemical autoencoder",
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
