#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path


INPUT_FILE = Path(__file__).resolve().parents[1] / "data" / "raw" / "variant_summary.txt.gz"
OUTPUT_FILE = Path(__file__).resolve().parents[1] / "data" / "interim" / "clinvar_splice_pathogenic.tsv"


def _contains_splice_signal(name: str) -> bool:
    lowered = name.lower()
    if any(
        token in lowered
        for token in [
            "splice",
            "+1",
            "+2",
            "+3",
            "+4",
            "+5",
            "-1",
            "-2",
            "-3",
            "-4",
            "-5",
            "acceptor",
            "donor",
            "intron",
            "ivs"
        ]
    ):
        return True
    return bool(
        re.search(
            r"c\.\d+(?:[+-]\d+)",
            lowered
        )
    )


def _is_pathogenic(significance: str) -> bool:
    lowered = significance.lower()
    return "pathogenic" in lowered


def main() -> int:
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "allele_id",
        "variation_id",
        "gene_symbol",
        "gene_id",
        "clinical_significance",
        "review_status",
        "assembly",
        "chromosome",
        "start",
        "stop",
        "name",
        "phenotype_list",
        "origin"
    ]
    seen: set[tuple[str, str]] = set()
    with gzip.open(INPUT_FILE, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        with OUTPUT_FILE.open("w", encoding="utf-8", newline="") as out_handle:
            writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in reader:
                if not _is_pathogenic(row.get("ClinicalSignificance", "")):
                    continue
                name = row.get("Name", "")
                if not _contains_splice_signal(name):
                    continue
                allele_id = row.get("AlleleID") or row.get("#AlleleID") or ""
                key = (allele_id, row.get("Assembly", ""), name)
                if key in seen:
                    continue
                seen.add(key)
                writer.writerow(
                    {
                        "allele_id": allele_id,
                        "variation_id": row.get("VariationID", ""),
                        "gene_symbol": row.get("GeneSymbol", ""),
                        "gene_id": row.get("GeneID", ""),
                        "clinical_significance": row.get("ClinicalSignificance", ""),
                        "review_status": row.get("ReviewStatus", ""),
                        "assembly": row.get("Assembly", ""),
                        "chromosome": row.get("Chromosome", ""),
                        "start": row.get("Start", ""),
                        "stop": row.get("Stop", ""),
                        "name": name,
                        "phenotype_list": row.get("PhenotypeList", ""),
                        "origin": row.get("Origin", "")
                    }
                )
    print(f"[ok] Wrote {OUTPUT_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
