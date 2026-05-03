#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "appris"
OUT = Path(__file__).resolve().parents[1] / "data" / "interim"

PRINCIPAL_FILES = [
    ("human", "GRCh38", RAW / "appris_human_GRCh38_principal.txt"),
    ("mouse", "GRCm39", RAW / "appris_mouse_GRCm39_principal.txt")
]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / "appris_principal_isoforms.tsv"
    fields = [
        "species",
        "assembly",
        "gene_id",
        "gene_name",
        "transcript_id",
        "translation_id",
        "ccds_id",
        "principal_tag",
        "principal_score"
    ]
    rows = 0
    with out_path.open("w", encoding="utf-8", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for species, assembly, path in PRINCIPAL_FILES:
            with path.open() as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    principal_tag = row.get("APPRIS Annotation", "")
                    principal_score = ""
                    if ":" in principal_tag:
                        principal_score = principal_tag.split(":", 1)[1]
                    writer.writerow(
                        {
                            "species": species,
                            "assembly": assembly,
                            "gene_id": row.get("Gene ID", ""),
                            "gene_name": row.get("Gene name", ""),
                            "transcript_id": row.get("Transcript ID", ""),
                            "translation_id": row.get("Translation ID", ""),
                            "ccds_id": row.get("CCDS ID", ""),
                            "principal_tag": principal_tag,
                            "principal_score": principal_score
                        }
                    )
                    rows += 1
    print(f"[ok] Wrote {out_path} rows={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
