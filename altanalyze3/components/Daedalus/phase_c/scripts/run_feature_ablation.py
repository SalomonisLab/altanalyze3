#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import OneHotEncoder, StandardScaler


PHASE_C = Path(__file__).resolve().parents[1]
PROCESSED = PHASE_C / "data" / "processed"
CHECKPOINTS = PHASE_C / "checkpoints"
IN_PATH = PROCESSED / "biochem_multitask_instances.parquet"
SCHEMA_PATH = PROCESSED / "biochem_multitask_schema.json"
OUT_TSV = CHECKPOINTS / "feature_ablation_metrics.tsv"
OUT_MD = CHECKPOINTS / "feature_ablation_report.md"


def _group_map(numeric_columns: list[str]) -> dict[str, list[str]]:
    segment = {
        "common_prefix_len",
        "common_suffix_len",
        "ref_changed_len",
        "alt_changed_len",
        "ref_changed_frac",
        "alt_changed_frac",
        "is_n_term_change",
        "is_c_term_change",
        "is_internal_change",
        "is_complex_change",
        "is_truncation",
        "is_extension",
        "ref_change_start",
        "ref_change_end",
    }
    groups: dict[str, list[str]] = {
        "base_reference_delta": [],
        "segment_change": [],
        "full_sequence_physchem": [],
        "terminal_sequence_physchem": [],
        "motifs": [],
        "typed_overlap": [],
    }
    for col in numeric_columns:
        if col in segment:
            groups["segment_change"].append(col)
        elif col.startswith(("ref_full_", "alt_full_")):
            groups["full_sequence_physchem"].append(col)
        elif col.startswith(("ref_seq_nterm_", "ref_seq_cterm_", "alt_seq_nterm_", "alt_seq_cterm_")):
            groups["terminal_sequence_physchem"].append(col)
        elif "_motif_" in col:
            groups["motifs"].append(col)
        elif col.startswith("overlap_"):
            groups["typed_overlap"].append(col)
        else:
            groups["base_reference_delta"].append(col)
    return groups


def _metrics(y_true: np.ndarray, probs: np.ndarray) -> dict[str, float]:
    return {
        "average_precision": float(average_precision_score(y_true, probs)),
        "auroc": float(roc_auc_score(y_true, probs)),
    }


def _run_task(df: pd.DataFrame, task: str, numeric_columns: list[str], cat_cols: list[str]) -> dict[str, float]:
    sub = df[df["task"] == task].copy()
    train_df = sub[sub["split"] == "train"].copy()
    test_df = sub[sub["split"] == "test"].copy()
    scaler = StandardScaler()
    train_num = scaler.fit_transform(train_df[numeric_columns].to_numpy(dtype=np.float32))
    test_num = scaler.transform(test_df[numeric_columns].to_numpy(dtype=np.float32))
    enc = OneHotEncoder(handle_unknown="ignore", sparse_output=True)
    train_cat = enc.fit_transform(train_df[cat_cols].astype(str))
    test_cat = enc.transform(test_df[cat_cols].astype(str))
    x_train = sparse.hstack([sparse.csr_matrix(train_num), train_cat], format="csr")
    x_test = sparse.hstack([sparse.csr_matrix(test_num), test_cat], format="csr")
    y_train = train_df["label_preservation_seed"].to_numpy(dtype=np.int32)
    y_test = test_df["label_preservation_seed"].to_numpy(dtype=np.int32)
    model = LogisticRegression(
        solver="saga",
        penalty="l2",
        C=1.0,
        max_iter=250,
        n_jobs=-1,
        class_weight="balanced",
        random_state=0,
    )
    model.fit(x_train, y_train)
    probs = model.predict_proba(x_test)[:, 1]
    return _metrics(y_test, probs)


def main() -> int:
    df = pd.read_parquet(IN_PATH)
    schema = json.loads(SCHEMA_PATH.read_text())
    numeric_columns = schema["numeric_columns"]
    feature_groups = _group_map(numeric_columns)
    tasks = sorted(df["task"].astype(str).unique().tolist())

    rows: list[dict[str, object]] = []
    full_cat = ["family_class", "reference_source", "change_class"]
    base_results: dict[str, dict[str, float]] = {}
    for task in tasks:
        metrics = _run_task(df, task, numeric_columns, full_cat)
        base_results[task] = metrics
        rows.append({"task": task, "ablation": "full", **metrics})

    for group_name, group_cols in feature_groups.items():
        if not group_cols:
            continue
        keep_cols = [c for c in numeric_columns if c not in set(group_cols)]
        cat_cols = full_cat.copy()
        if group_name == "segment_change" and "change_class" in cat_cols:
            cat_cols.remove("change_class")
        for task in tasks:
            metrics = _run_task(df, task, keep_cols, cat_cols)
            rows.append({"task": task, "ablation": f"drop_{group_name}", **metrics})

    out = pd.DataFrame(rows).sort_values(["task", "ablation"]).reset_index(drop=True)
    out["delta_ap_vs_full"] = np.nan
    out["delta_auroc_vs_full"] = np.nan
    for task in tasks:
        mask = out["task"] == task
        full = out[mask & (out["ablation"] == "full")].iloc[0]
        out.loc[mask, "delta_ap_vs_full"] = out.loc[mask, "average_precision"] - full["average_precision"]
        out.loc[mask, "delta_auroc_vs_full"] = out.loc[mask, "auroc"] - full["auroc"]

    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    lines = [
        "# Phase C Feature Ablation",
        "",
        "Ablations use leakage-resistant logistic models on the Phase C rich-feature matrix.",
        "Negative deltas indicate a useful feature block.",
        "",
    ]
    for task in tasks:
        lines.append(f"## {task}")
        lines.append("")
        sub = out[out["task"] == task].copy().sort_values("delta_ap_vs_full")
        for row in sub.itertuples(index=False):
            lines.append(
                f"- `{row.ablation}` AP={row.average_precision:.6f} "
                f"(delta {row.delta_ap_vs_full:+.6f}) AUROC={row.auroc:.6f} "
                f"(delta {row.delta_auroc_vs_full:+.6f})"
            )
        lines.append("")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[ok] Wrote {OUT_TSV}")
    print(f"[ok] Wrote {OUT_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
