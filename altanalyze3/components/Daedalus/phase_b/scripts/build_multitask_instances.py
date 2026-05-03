#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"
IN_PATH = PROCESSED / "reference_delta_matrix.parquet"
OUT_TSV = PROCESSED / "multitask_instances.tsv"
OUT_PARQUET = PROCESSED / "multitask_instances.parquet"
OUT_SCHEMA = PROCESSED / "multitask_schema.json"

TASK_MAP = {
    "task_global_seed": "global",
    "task_membrane_seed": "membrane",
    "task_kinase_seed": "kinase",
    "task_tf_seed": "transcription_factor",
}


def main() -> int:
    df = pd.read_parquet(IN_PATH)
    numeric_columns = [
        c
        for c in df.columns
        if c
        not in {
            "species",
            "gene_id",
            "gene_name",
            "split",
            "family_class",
            "reference_source",
            "reference_transcript_id",
            "alternative_transcript_id",
            "label_preservation_seed",
            *TASK_MAP.keys(),
        }
    ]

    rows = []
    for task_col, task_name in TASK_MAP.items():
        sub = df[(df[task_col] == 1) & (df["label_preservation_seed"].notna())].copy()
        sub["task"] = task_name
        rows.append(sub)
    out = pd.concat(rows, axis=0, ignore_index=True)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    out.to_parquet(OUT_PARQUET, index=False)

    schema = {
        "task_map": TASK_MAP,
        "categorical_columns": ["species", "family_class", "reference_source", "task"],
        "id_columns": ["gene_id", "gene_name", "reference_transcript_id", "alternative_transcript_id"],
        "split_column": "split",
        "label_column": "label_preservation_seed",
        "numeric_columns": numeric_columns,
    }
    with OUT_SCHEMA.open("w", encoding="utf-8") as handle:
        json.dump(schema, handle, indent=2)

    print(f"[ok] Wrote {OUT_TSV} rows={len(out)}")
    print(f"[ok] Wrote {OUT_PARQUET}")
    print(f"[ok] Wrote {OUT_SCHEMA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
