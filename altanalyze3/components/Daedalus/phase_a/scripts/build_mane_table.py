#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
from pathlib import Path


RAW = Path(__file__).resolve().parents[1] / "data" / "raw"
OUT = Path(__file__).resolve().parents[1] / "data" / "interim"
MANE_IN = RAW / "MANE.GRCh38.v1.5.summary.txt.gz"
MANE_OUT = OUT / "mane_transcripts.tsv"


def _base_id(value: str) -> str:
    return value.split(".", 1)[0] if value else ""


def main() -> int:
    rows = 0
    with gzip.open(MANE_IN, "rt") as in_handle, MANE_OUT.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        writer = csv.DictWriter(
            out_handle,
            fieldnames=[
                "species",
                "gene_id",
                "gene_id_base",
                "gene_name",
                "transcript_id",
                "transcript_id_base",
                "protein_id",
                "mane_status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in reader:
            writer.writerow(
                {
                    "species": "human",
                    "gene_id": row.get("Ensembl_Gene", ""),
                    "gene_id_base": _base_id(row.get("Ensembl_Gene", "")),
                    "gene_name": row.get("symbol", ""),
                    "transcript_id": row.get("Ensembl_nuc", ""),
                    "transcript_id_base": _base_id(row.get("Ensembl_nuc", "")),
                    "protein_id": row.get("Ensembl_prot", ""),
                    "mane_status": row.get("MANE_status", ""),
                }
            )
            rows += 1
    print(f"[ok] Wrote {MANE_OUT} rows={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
