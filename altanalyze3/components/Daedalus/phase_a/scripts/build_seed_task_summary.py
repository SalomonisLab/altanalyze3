#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"
IN_PATH = PROCESSED / "baseline_feature_matrix.tsv"
OUT_PATH = PROCESSED / "seed_task_summary.tsv"
TASKS = [
    "task_global_seed",
    "task_membrane_seed",
    "task_kinase_seed",
    "task_tf_seed",
    "task_high_confidence_seed",
]


def main() -> int:
    counts: dict[tuple[str, str], Counter] = {}
    with IN_PATH.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            label = row["label_preservation_seed"]
            for task in TASKS:
                if row[task] != "1":
                    continue
                key = (task, row["split"])
                bucket = counts.setdefault(key, Counter())
                bucket["n"] += 1
                if label == "1":
                    bucket["positive"] += 1
                elif label == "0":
                    bucket["negative"] += 1

    with OUT_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["task", "split", "n", "positive", "negative"], delimiter="\t")
        writer.writeheader()
        for (task, split), bucket in sorted(counts.items()):
            writer.writerow(
                {
                    "task": task,
                    "split": split,
                    "n": bucket["n"],
                    "positive": bucket["positive"],
                    "negative": bucket["negative"],
                }
            )
    print(f"[ok] Wrote {OUT_PATH} rows={len(counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
