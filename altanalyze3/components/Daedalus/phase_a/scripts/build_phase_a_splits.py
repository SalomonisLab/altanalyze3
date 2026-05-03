#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
from collections import Counter, defaultdict
from pathlib import Path


INTERIM = Path(__file__).resolve().parents[1] / "data" / "interim"
PAIR_IN = INTERIM / "isoform_pair_candidates.tsv"
GENE_IN = INTERIM / "benchmark_gene_sets.tsv"
PAIR_OUT = INTERIM / "isoform_pair_candidates.with_splits.tsv"
GENE_OUT = INTERIM / "benchmark_gene_sets.with_splits.tsv"
SUMMARY_OUT = INTERIM / "phase_a_split_summary.tsv"


def _assign_split(species: str, gene_id: str) -> str:
    digest = hashlib.md5(f"{species}:{gene_id}".encode("utf-8")).hexdigest()
    bucket = int(digest[:8], 16) % 100
    if bucket < 80:
        return "train"
    if bucket < 90:
        return "val"
    return "test"


def _write_with_splits(in_path: Path, out_path: Path) -> Counter:
    counts: Counter = Counter()
    with in_path.open() as in_handle, out_path.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or []) + ["split"]
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            split = _assign_split(row["species"], row["gene_id"])
            row["split"] = split
            writer.writerow(row)
            counts["total"] += 1
            counts[f"split_{split}"] += 1
            family = row.get("family_class", "other") or "other"
            counts[f"{family}_{split}"] += 1
    return counts


def _write_summary(pair_counts: Counter, gene_counts: Counter) -> None:
    rows: list[dict[str, str]] = []
    for dataset, counts in (("pairs", pair_counts), ("genes", gene_counts)):
        for key, value in sorted(counts.items()):
            if key == "total":
                metric = "total"
                subgroup = ""
            elif key.startswith("split_"):
                metric = "split_total"
                subgroup = key.replace("split_", "", 1)
            else:
                metric = "family_split"
                subgroup = key
            rows.append(
                {
                    "dataset": dataset,
                    "metric": metric,
                    "subgroup": subgroup,
                    "value": str(value),
                }
            )
    with SUMMARY_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset", "metric", "subgroup", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    pair_counts = _write_with_splits(PAIR_IN, PAIR_OUT)
    gene_counts = _write_with_splits(GENE_IN, GENE_OUT)
    _write_summary(pair_counts, gene_counts)
    print(f"[ok] Wrote {PAIR_OUT} rows={pair_counts['total']}")
    print(f"[ok] Wrote {GENE_OUT} rows={gene_counts['total']}")
    print(f"[ok] Wrote {SUMMARY_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
