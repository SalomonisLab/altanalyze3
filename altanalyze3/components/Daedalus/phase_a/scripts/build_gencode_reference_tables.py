#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
from pathlib import Path


RAW = Path(__file__).resolve().parents[1] / "data" / "raw" / "gencode"
OUT = Path(__file__).resolve().parents[1] / "data" / "interim"

RESOURCES = [
    {
        "species": "human",
        "release": "v49",
        "gtf": RAW / "gencode.v49.annotation.gtf.gz",
        "transcripts": RAW / "gencode.v49.pc_transcripts.fa.gz",
        "proteins": RAW / "gencode.v49.pc_translations.fa.gz"
    },
    {
        "species": "mouse",
        "release": "vM38",
        "gtf": RAW / "gencode.vM38.annotation.gtf.gz",
        "transcripts": RAW / "gencode.vM38.pc_transcripts.fa.gz",
        "proteins": RAW / "gencode.vM38.pc_translations.fa.gz"
    }
]


def _parse_attributes(field: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in field.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        key, _, value = item.partition(" ")
        result[key] = value.strip().strip('"')
    return result


def _parse_fasta_lengths(path: Path, is_protein: bool) -> dict[str, dict[str, str | int]]:
    result: dict[str, dict[str, str | int]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        current_header: str | None = None
        current_seq: list[str] = []
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    header_parts = current_header[1:].split("|")
                    if is_protein:
                        protein_id = header_parts[0]
                        transcript_id = header_parts[1]
                        gene_id = header_parts[2]
                        result[transcript_id] = {
                            "protein_id": protein_id,
                            "gene_id": gene_id,
                            "protein_length": len("".join(current_seq))
                        }
                    else:
                        transcript_id = header_parts[0]
                        gene_id = header_parts[1]
                        transcript_name = header_parts[4]
                        gene_name = header_parts[5]
                        result[transcript_id] = {
                            "gene_id": gene_id,
                            "transcript_name": transcript_name,
                            "gene_name": gene_name,
                            "transcript_length": len("".join(current_seq))
                        }
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            header_parts = current_header[1:].split("|")
            if is_protein:
                protein_id = header_parts[0]
                transcript_id = header_parts[1]
                gene_id = header_parts[2]
                result[transcript_id] = {
                    "protein_id": protein_id,
                    "gene_id": gene_id,
                    "protein_length": len("".join(current_seq))
                }
            else:
                transcript_id = header_parts[0]
                gene_id = header_parts[1]
                transcript_name = header_parts[4]
                gene_name = header_parts[5]
                result[transcript_id] = {
                    "gene_id": gene_id,
                    "transcript_name": transcript_name,
                    "gene_name": gene_name,
                    "transcript_length": len("".join(current_seq))
                }
    return result


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    transcript_out = OUT / "gencode_transcript_reference.tsv"
    protein_out = OUT / "gencode_protein_reference.tsv"

    transcript_fields = [
        "species",
        "release",
        "gene_id",
        "gene_name",
        "transcript_id",
        "transcript_name",
        "gene_type",
        "transcript_type",
        "chromosome",
        "strand",
        "start",
        "end",
        "transcript_length",
        "protein_id",
        "protein_length"
    ]
    protein_fields = [
        "species",
        "release",
        "gene_id",
        "transcript_id",
        "protein_id",
        "protein_length"
    ]

    with transcript_out.open("w", encoding="utf-8", newline="") as transcript_handle, \
        protein_out.open("w", encoding="utf-8", newline="") as protein_handle:
        transcript_writer = csv.DictWriter(transcript_handle, fieldnames=transcript_fields, delimiter="\t")
        protein_writer = csv.DictWriter(protein_handle, fieldnames=protein_fields, delimiter="\t")
        transcript_writer.writeheader()
        protein_writer.writeheader()

        total_transcripts = 0
        total_proteins = 0
        for resource in RESOURCES:
            transcript_lengths = _parse_fasta_lengths(resource["transcripts"], is_protein=False)
            protein_lengths = _parse_fasta_lengths(resource["proteins"], is_protein=True)
            with gzip.open(resource["gtf"], "rt", encoding="utf-8", errors="replace") as gtf_handle:
                for line in gtf_handle:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.rstrip().split("\t")
                    if len(parts) != 9 or parts[2] != "transcript":
                        continue
                    attrs = _parse_attributes(parts[8])
                    transcript_id = attrs.get("transcript_id", "")
                    gene_id = attrs.get("gene_id", "")
                    transcript_meta = transcript_lengths.get(transcript_id, {})
                    protein_meta = protein_lengths.get(transcript_id, {})
                    transcript_writer.writerow(
                        {
                            "species": resource["species"],
                            "release": resource["release"],
                            "gene_id": gene_id,
                            "gene_name": attrs.get("gene_name", transcript_meta.get("gene_name", "")),
                            "transcript_id": transcript_id,
                            "transcript_name": attrs.get("transcript_name", transcript_meta.get("transcript_name", "")),
                            "gene_type": attrs.get("gene_type", ""),
                            "transcript_type": attrs.get("transcript_type", ""),
                            "chromosome": parts[0],
                            "strand": parts[6],
                            "start": parts[3],
                            "end": parts[4],
                            "transcript_length": transcript_meta.get("transcript_length", ""),
                            "protein_id": protein_meta.get("protein_id", ""),
                            "protein_length": protein_meta.get("protein_length", "")
                        }
                    )
                    total_transcripts += 1
            for transcript_id, protein_meta in protein_lengths.items():
                protein_writer.writerow(
                    {
                        "species": resource["species"],
                        "release": resource["release"],
                        "gene_id": protein_meta.get("gene_id", ""),
                        "transcript_id": transcript_id,
                        "protein_id": protein_meta.get("protein_id", ""),
                        "protein_length": protein_meta.get("protein_length", "")
                    }
                )
                total_proteins += 1

    print(f"[ok] Wrote {transcript_out} transcripts={total_transcripts}")
    print(f"[ok] Wrote {protein_out} proteins={total_proteins}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
