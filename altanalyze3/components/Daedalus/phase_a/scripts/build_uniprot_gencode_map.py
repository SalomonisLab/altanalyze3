#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


INTERIM = Path(__file__).resolve().parents[1] / "data" / "interim"
UNIPROT_FILE = INTERIM / "uniprot_entries.tsv"
GTF_TRANSCRIPT_FILE = INTERIM / "gencode_transcript_reference.tsv"
GTF_PROTEIN_FILE = INTERIM / "gencode_protein_reference.tsv"
OUT_FILE = INTERIM / "uniprot_gencode_map.tsv"


def _load_gencode_maps() -> tuple[dict[str, dict], dict[str, dict]]:
    transcript_map: dict[str, dict] = {}
    protein_map: dict[str, dict] = {}
    with GTF_TRANSCRIPT_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            transcript_id = row["transcript_id"]
            if transcript_id:
                transcript_map[transcript_id] = row
    with GTF_PROTEIN_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            protein_id = row["protein_id"]
            if protein_id:
                protein_map[protein_id] = row
    return transcript_map, protein_map


def main() -> int:
    transcript_map, protein_map = _load_gencode_maps()
    with UNIPROT_FILE.open() as in_handle, OUT_FILE.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = [
            "primary_accession",
            "gene_name",
            "organism",
            "ensembl_id",
            "match_type",
            "gencode_gene_id",
            "gencode_transcript_id",
            "gencode_protein_id",
            "gencode_species",
            "gencode_release"
        ]
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        rows_written = 0
        for row in reader:
            seen: set[tuple[str, str]] = set()
            for token in [x for x in row.get("ensembl_ids", "").split(";") if x]:
                if token in transcript_map:
                    meta = transcript_map[token]
                    key = (token, "transcript")
                    if key not in seen:
                        seen.add(key)
                        writer.writerow(
                            {
                                "primary_accession": row["primary_accession"],
                                "gene_name": row.get("gene_name", ""),
                                "organism": row.get("organism", ""),
                                "ensembl_id": token,
                                "match_type": "transcript",
                                "gencode_gene_id": meta.get("gene_id", ""),
                                "gencode_transcript_id": meta.get("transcript_id", ""),
                                "gencode_protein_id": meta.get("protein_id", ""),
                                "gencode_species": meta.get("species", ""),
                                "gencode_release": meta.get("release", "")
                            }
                        )
                        rows_written += 1
                if token in protein_map:
                    meta = protein_map[token]
                    key = (token, "protein")
                    if key not in seen:
                        seen.add(key)
                        writer.writerow(
                            {
                                "primary_accession": row["primary_accession"],
                                "gene_name": row.get("gene_name", ""),
                                "organism": row.get("organism", ""),
                                "ensembl_id": token,
                                "match_type": "protein",
                                "gencode_gene_id": meta.get("gene_id", ""),
                                "gencode_transcript_id": meta.get("transcript_id", ""),
                                "gencode_protein_id": meta.get("protein_id", ""),
                                "gencode_species": meta.get("species", ""),
                                "gencode_release": meta.get("release", "")
                            }
                        )
                        rows_written += 1
    print(f"[ok] Wrote {OUT_FILE} rows={rows_written}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
