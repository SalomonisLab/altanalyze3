#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge, SGDClassifier
from sklearn.metrics import average_precision_score, balanced_accuracy_score, mean_absolute_error, r2_score, roc_auc_score


PHASE_B = Path(__file__).resolve().parents[1]
PROCESSED = PHASE_B / "data" / "processed"
IN_PATH = PROCESSED / "supervised_benchmark_instances.parquet"
SCHEMA_PATH = PROCESSED / "supervised_benchmark_schema.json"
OUT_TSV = PROCESSED / "supervised_holdout_metrics.tsv"
OUT_MD = PROCESSED / "supervised_holdout_report.md"

NON_FEATURE_COLUMNS = {
    "task",
    "benchmark_source",
    "label_binary",
    "label_continuous",
    "label_class",
    "evidence_strength",
    "cell_context",
    "notes",
    "split",
}


def _safe_metric(metric_fn, *args):
    try:
        return float(metric_fn(*args))
    except Exception:
        return float("nan")


def _build_matrix(train_df: pd.DataFrame, eval_df: pd.DataFrame, feature_columns: list[str]) -> tuple[np.ndarray, np.ndarray]:
    train_x = train_df[feature_columns].copy().reset_index(drop=True)
    eval_x = eval_df[feature_columns].copy().reset_index(drop=True)
    combined = pd.concat([train_x, eval_x], axis=0, ignore_index=True)
    encoded = pd.DataFrame(index=combined.index)
    for col in feature_columns:
        series = combined[col]
        if pd.api.types.is_numeric_dtype(series):
            encoded[col] = pd.to_numeric(series, errors="coerce").fillna(0.0).astype(np.float32)
            continue
        categories = pd.Categorical(series.fillna("__NA__"))
        encoded[col] = categories.codes.astype(np.float32)
    train_matrix = encoded.iloc[: len(train_x)].to_numpy(dtype=np.float32)
    eval_matrix = encoded.iloc[len(train_x) :].to_numpy(dtype=np.float32)
    return train_matrix, eval_matrix


def _binary_metrics(task: str, split: str, y_true: np.ndarray, prob: np.ndarray) -> dict[str, object]:
    pred = (prob >= 0.5).astype(int)
    return {
        "task": task,
        "split": split,
        "label_type": "binary",
        "n": int(len(y_true)),
        "positives": int(y_true.sum()),
        "auroc": _safe_metric(roc_auc_score, y_true, prob),
        "average_precision": _safe_metric(average_precision_score, y_true, prob),
        "balanced_accuracy": _safe_metric(balanced_accuracy_score, y_true, pred),
        "mae": float("nan"),
        "r2": float("nan"),
    }


def _continuous_metrics(task: str, split: str, y_true: np.ndarray, pred: np.ndarray) -> dict[str, object]:
    return {
        "task": task,
        "split": split,
        "label_type": "continuous",
        "n": int(len(y_true)),
        "positives": 0,
        "auroc": float("nan"),
        "average_precision": float("nan"),
        "balanced_accuracy": float("nan"),
        "mae": _safe_metric(mean_absolute_error, y_true, pred),
        "r2": _safe_metric(r2_score, y_true, pred),
    }


def main() -> int:
    df = pd.read_parquet(IN_PATH)
    with SCHEMA_PATH.open() as handle:
        schema = json.load(handle)

    if df.empty:
        pd.DataFrame(
            columns=["task", "split", "label_type", "n", "positives", "auroc", "average_precision", "balanced_accuracy", "mae", "r2"]
        ).to_csv(OUT_TSV, sep="\t", index=False)
        OUT_MD.write_text("# Supervised Holdout Baselines\n\nNo supervised benchmark rows were available.\n", encoding="utf-8")
        print(f"[ok] Wrote {OUT_TSV} rows=0")
        print(f"[ok] Wrote {OUT_MD}")
        return 0

    feature_columns = [c for c in df.columns if c not in NON_FEATURE_COLUMNS]
    rows: list[dict[str, object]] = []
    for task, task_df in sorted(df.groupby("task")):
        train_df = task_df[task_df["split"] == "train"].copy()
        val_df = task_df[task_df["split"] == "val"].copy()
        test_df = task_df[task_df["split"] == "test"].copy()
        if train_df.empty or (val_df.empty and test_df.empty):
            continue

        label_type = None
        if task_df["label_binary"].notna().any():
            label_type = "binary"
        elif task_df["label_continuous"].notna().any():
            label_type = "continuous"
        else:
            continue

        if label_type == "binary":
            train_sub = train_df[train_df["label_binary"].notna()].copy()
            if train_sub["label_binary"].nunique() < 2:
                continue
            for split_name, split_df in [("val", val_df), ("test", test_df)]:
                eval_sub = split_df[split_df["label_binary"].notna()].copy()
                if eval_sub.empty:
                    continue
                x_train, x_eval = _build_matrix(train_sub, eval_sub, feature_columns)
                y_train = train_sub["label_binary"].astype(int).to_numpy()
                y_eval = eval_sub["label_binary"].astype(int).to_numpy()
                model = SGDClassifier(
                    loss="log_loss",
                    class_weight="balanced",
                    max_iter=1000,
                    tol=1e-3,
                    random_state=0,
                )
                model.fit(x_train, y_train)
                prob = model.predict_proba(x_eval)[:, 1]
                rows.append(_binary_metrics(task, split_name, y_eval, prob))
        else:
            train_sub = train_df[train_df["label_continuous"].notna()].copy()
            if len(train_sub) < 4:
                continue
            for split_name, split_df in [("val", val_df), ("test", test_df)]:
                eval_sub = split_df[split_df["label_continuous"].notna()].copy()
                if eval_sub.empty:
                    continue
                x_train, x_eval = _build_matrix(train_sub, eval_sub, feature_columns)
                y_train = train_sub["label_continuous"].astype(float).to_numpy()
                y_eval = eval_sub["label_continuous"].astype(float).to_numpy()
                model = Ridge(alpha=1.0)
                model.fit(x_train, y_train)
                pred = model.predict(x_eval)
                rows.append(_continuous_metrics(task, split_name, y_eval, pred))

    out = pd.DataFrame(rows).sort_values(["task", "split"]) if rows else pd.DataFrame(
        columns=["task", "split", "label_type", "n", "positives", "auroc", "average_precision", "balanced_accuracy", "mae", "r2"]
    )
    out.to_csv(OUT_TSV, sep="\t", index=False)

    lines = ["# Supervised Holdout Baselines", ""]
    if out.empty:
        lines.append("No holdout-eligible supervised benchmark rows were available.")
    else:
        for task, task_rows in out.groupby("task"):
            lines.append(f"## {task}")
            lines.append("")
            for _, row in task_rows.iterrows():
                if row["label_type"] == "binary":
                    lines.append(
                        f"- `{row['split']}` n={int(row['n'])} AP={row['average_precision']:.6f} "
                        f"AUROC={row['auroc']:.6f} BalAcc={row['balanced_accuracy']:.6f}"
                    )
                else:
                    lines.append(
                        f"- `{row['split']}` n={int(row['n'])} MAE={row['mae']:.6f} R2={row['r2']:.6f}"
                    )
            lines.append("")
    OUT_MD.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")

    print(f"[ok] Wrote {OUT_TSV} rows={len(out)}")
    print(f"[ok] Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
