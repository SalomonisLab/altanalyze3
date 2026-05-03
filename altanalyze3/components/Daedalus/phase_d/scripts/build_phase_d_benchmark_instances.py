#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


PHASE_D = Path(__file__).resolve().parents[1]
PROCESSED = PHASE_D / "data" / "processed"
PHASE_A_INTERIM = Path(__file__).resolve().parents[2] / "phase_a" / "data" / "interim"
TM_NEGATIVE_LABELS = PHASE_A_INTERIM / "tm_negatives" / "tm_negative_pair_labels.tsv"

IN_PATH = PROCESSED / "phase_d_objective_scores.parquet"
REGISTRY_PATH = PROCESSED / "phase_d_task_registry.tsv"
OUT_TSV = PROCESSED / "phase_d_benchmark_instances.tsv"
OUT_PARQUET = PROCESSED / "phase_d_benchmark_instances.parquet"
OUT_SCHEMA = PROCESSED / "phase_d_benchmark_schema.json"

HIGH_POS = 0.70
HIGH_NEG = 0.30
APPLICABILITY_COLUMN_BY_FAMILY = {
    "membrane": "task_membrane_seed",
    "transcription_factor": "task_tf_seed",
}


def _label_from_score(score: float, direction: str) -> tuple[int | None, str]:
    if pd.isna(score):
        return None, "missing"
    score = float(score)
    if direction == "higher_is_better":
        if score >= HIGH_POS:
            return 1, "positive"
        if score <= HIGH_NEG:
            return 0, "negative"
        return None, "ambiguous_exclude"
    if direction == "lower_is_better":
        if score <= HIGH_NEG:
            return 1, "positive"
        if score >= HIGH_POS:
            return 0, "negative"
        return None, "ambiguous_exclude"
    raise ValueError(f"Unknown direction: {direction}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build Phase D benchmark instances from objective scores.")
    parser.add_argument(
        "--skip-tsv",
        action="store_true",
        help="Skip the large TSV export and write only parquet/schema.",
    )
    args = parser.parse_args()

    df = pd.read_parquet(IN_PATH)
    registry = pd.read_csv(REGISTRY_PATH, sep="\t")
    tm_negative_labels = None
    if TM_NEGATIVE_LABELS.exists():
        tm_negative_labels = pd.read_csv(TM_NEGATIVE_LABELS, sep="\t")
        tm_negative_labels["pair_key"] = tm_negative_labels.apply(
            lambda row: (
                str(row["species"]),
                str(row["gene_name"]),
                str(row["reference_transcript_id"]),
                str(row["alternative_transcript_id"]),
                str(row["label_task"]),
            ),
            axis=1,
        )
        tm_negative_map = {
            key: row._asdict() if hasattr(row, "_asdict") else None
            for key, row in zip(
                tm_negative_labels["pair_key"],
                tm_negative_labels.itertuples(index=False),
            )
        }
        tm_negative_keys = set(tm_negative_map)
    else:
        tm_negative_map = {}
        tm_negative_keys = set()
    objective_score_columns = set(registry["score_column"].astype(str).tolist()) | {
        "aggregate_stable_functional_score",
        "aggregate_non_stable_non_functional_score",
        "aggregate_functional_annotation",
    }
    rows = []
    id_cols = [
        "species",
        "gene_id",
        "gene_name",
        "reference_transcript_id",
        "alternative_transcript_id",
        "family_class",
        "reference_source",
        "split",
        "change_class",
        "aggregate_functional_annotation",
    ]
    feature_columns = [
        c
        for c in df.columns
        if c not in set(id_cols) | objective_score_columns
    ]
    curated_override_count = 0
    for task_row in registry.itertuples(index=False):
        score_col = task_row.score_column
        task_name = task_row.task
        direction = task_row.direction
        task_family = task_row.task_family
        sub = df[id_cols + feature_columns + [score_col]].copy()
        applicability_col = APPLICABILITY_COLUMN_BY_FAMILY.get(task_family)
        if applicability_col and applicability_col in df.columns:
            mask = pd.to_numeric(df[applicability_col], errors="coerce").fillna(0).astype(int) == 1
            sub = sub[mask].copy()
        sub["task"] = task_name
        sub["task_family"] = task_family
        sub["direction"] = direction
        sub["objective_score"] = sub[score_col]
        labels = sub["objective_score"].map(lambda s: _label_from_score(s, direction))
        sub["label_binary"] = labels.map(lambda x: x[0]).astype("Int64")
        sub["label_state"] = labels.map(lambda x: x[1])
        sub["label_source"] = "objective_threshold"
        if task_name == "cell_surface_localization" and tm_negative_keys:
            pair_keys = sub.apply(
                lambda row: (
                    str(row["species"]),
                    str(row["gene_name"]),
                    str(row["reference_transcript_id"]),
                    str(row["alternative_transcript_id"]),
                    task_name,
                ),
                axis=1,
            )
            override_mask = pair_keys.isin(tm_negative_keys)
            curated_override_count += int(override_mask.sum())
            sub.loc[override_mask, "label_binary"] = 0
            sub.loc[override_mask, "label_state"] = "negative"
            sub.loc[override_mask, "label_source"] = "curated_tm_negative"
        rows.append(sub.drop(columns=[score_col]))

    out = pd.concat(rows, axis=0, ignore_index=True)
    out.to_parquet(OUT_PARQUET, index=False)

    schema = {
        "source_scores": str(IN_PATH),
        "source_registry": str(REGISTRY_PATH),
        "label_thresholds": {
            "higher_is_better_positive": HIGH_POS,
            "higher_is_better_negative": HIGH_NEG,
            "lower_is_better_positive": HIGH_NEG,
            "lower_is_better_negative": HIGH_POS,
        },
        "categorical_columns": ["family_class", "reference_source", "change_class", "task", "task_family"],
        "id_columns": ["gene_id", "gene_name", "reference_transcript_id", "alternative_transcript_id"],
        "split_column": "split",
        "label_column": "label_binary",
        "label_state_column": "label_state",
        "score_column": "objective_score",
        "label_source_column": "label_source",
        "feature_columns": feature_columns,
        "curated_tm_negative_labels": str(TM_NEGATIVE_LABELS) if TM_NEGATIVE_LABELS.exists() else "",
    }
    OUT_SCHEMA.write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")
    print(f"[ok] Wrote {OUT_PARQUET}")
    print(f"[ok] Wrote {OUT_SCHEMA}")
    print(f"[ok] Curated TM negative overrides={curated_override_count}")
    if not args.skip_tsv:
        out.to_csv(OUT_TSV, sep="\t", index=False)
        print(f"[ok] Wrote {OUT_TSV}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
