#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


PHASE_B = Path(__file__).resolve().parents[1]
CONFIG = PHASE_B / "config" / "task_registry.json"
PROCESSED = PHASE_B / "data" / "processed"
REFERENCE_MATRIX = PROCESSED / "reference_delta_matrix.parquet"
EXTERNAL_DEFAULT = PHASE_B / "data" / "external"
OUT_TSV = PROCESSED / "supervised_benchmark_instances.tsv"
OUT_PARQUET = PROCESSED / "supervised_benchmark_instances.parquet"
OUT_SCHEMA = PROCESSED / "supervised_benchmark_schema.json"
OUT_SUMMARY = PROCESSED / "supervised_benchmark_summary.tsv"


ID_COLUMNS = [
    "species",
    "gene_name",
    "reference_transcript_id",
    "alternative_transcript_id",
]


def _read_registry(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def _load_external_tables(external_dir: Path) -> list[tuple[Path, pd.DataFrame]]:
    tables: list[tuple[Path, pd.DataFrame]] = []
    if not external_dir.exists():
        return tables
    for path in sorted(external_dir.glob("*.tsv")):
        tables.append((path, pd.read_csv(path, sep="\t")))
    return tables


def _normalize_labels(frame: pd.DataFrame, task_spec: dict, source_name: str) -> pd.DataFrame:
    label_binary = None
    label_continuous = None
    label_class = None
    for col in task_spec.get("preferred_label_columns", []):
        if col not in frame.columns:
            continue
        if pd.api.types.is_numeric_dtype(frame[col]):
            unique = set(frame[col].dropna().astype(float).unique().tolist())
            if unique.issubset({0.0, 1.0}) and label_binary is None:
                label_binary = col
            elif label_continuous is None:
                label_continuous = col
        elif label_class is None:
            label_class = col
    out = frame.copy()
    out["task"] = task_spec["task_name"]
    out["benchmark_source"] = source_name
    out["label_binary"] = out[label_binary] if label_binary else pd.NA
    out["label_continuous"] = out[label_continuous] if label_continuous else pd.NA
    out["label_class"] = out[label_class] if label_class else pd.NA
    if "evidence_strength" not in out.columns:
        out["evidence_strength"] = pd.NA
    if "cell_context" not in out.columns:
        out["cell_context"] = pd.NA
    if "notes" not in out.columns:
        out["notes"] = pd.NA
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Build supervised benchmark instances for Phase B external tasks.")
    parser.add_argument("--external-dir", default=str(EXTERNAL_DEFAULT), help="Directory containing external benchmark TSVs.")
    args = parser.parse_args()

    registry = _read_registry(CONFIG)
    ref = pd.read_parquet(REFERENCE_MATRIX)
    ref_features = ref.drop_duplicates(subset=ID_COLUMNS)
    ref_feature_columns = [c for c in ref_features.columns if c not in ID_COLUMNS]
    ref_features = ref_features[ID_COLUMNS + ref_feature_columns]

    tables = _load_external_tables(Path(args.external_dir))
    rows: list[pd.DataFrame] = []
    for path, frame in tables:
        missing = [col for col in ID_COLUMNS if col not in frame.columns]
        if missing:
            raise ValueError(f"{path} is missing required ID columns: {missing}")
        if "task" not in frame.columns:
            raise ValueError(f"{path} must include a 'task' column")
        for task_spec in registry.get("supervised_tasks", []):
            task_name = task_spec["task_name"]
            sub = frame[frame["task"] == task_name].copy()
            if sub.empty:
                continue
            normalized = _normalize_labels(sub, task_spec, path.stem)
            merged = normalized.merge(ref_features, how="inner", on=ID_COLUMNS, suffixes=("", "_ref"))
            if not merged.empty:
                rows.append(merged)

    if rows:
        out = pd.concat(rows, axis=0, ignore_index=True)
    else:
        out = pd.DataFrame(columns=ID_COLUMNS + ["task", "benchmark_source", "label_binary", "label_continuous", "label_class", "evidence_strength", "cell_context", "notes"] + ref_feature_columns)

    out.to_csv(OUT_TSV, sep="\t", index=False)
    out.to_parquet(OUT_PARQUET, index=False)

    summary = (
        out.groupby(["task", "benchmark_source", "split"], dropna=False)
        .size()
        .reset_index(name="n_rows")
        .sort_values(["task", "benchmark_source", "split"])
    )
    summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

    schema = {
        "id_columns": ID_COLUMNS,
        "task_registry": str(CONFIG),
        "label_columns": ["label_binary", "label_continuous", "label_class"],
        "context_columns": ["family_class", "reference_source", "cell_context", "evidence_strength", "notes", "split"],
        "feature_columns": [c for c in ref_feature_columns if c != "split"],
    }
    with OUT_SCHEMA.open("w", encoding="utf-8") as handle:
        json.dump(schema, handle, indent=2)

    print(f"[ok] Wrote {OUT_TSV} rows={len(out)}")
    print(f"[ok] Wrote {OUT_PARQUET}")
    print(f"[ok] Wrote {OUT_SUMMARY}")
    print(f"[ok] Wrote {OUT_SCHEMA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
