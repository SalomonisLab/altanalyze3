#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import importlib.util
import json
from pathlib import Path


INTERIM = Path(__file__).resolve().parents[1] / "data" / "interim"
RAW_UNIPROT = Path(__file__).resolve().parents[1] / "data" / "raw" / "uniprot_reviewed_human_mouse.jsonl.gz"
PAIR_IN = INTERIM / "isoform_pair_candidates.with_splits.tsv"
PAIR_OUT = INTERIM / "priority_pair_subsets.tsv"
STRICT_SURFACE_OUT = INTERIM / "strict_surface_gene_set.tsv"


def _norm_species(value: str) -> str:
    value = (value or "").strip().lower()
    if value in {"homo sapiens", "human", "hs"}:
        return "human"
    if value in {"mus musculus", "mouse", "mm"}:
        return "mouse"
    return value


def _load_build_row():
    module_path = Path(__file__).resolve().parents[3] / "annotation" / "tm_negatives" / "build_isoform_catalog.py"
    spec = importlib.util.spec_from_file_location("tm_catalog_builder", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load TM catalog builder from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.build_row


def _strict_surface_gene_map() -> dict[tuple[str, str], dict[str, str]]:
    build_row = _load_build_row()
    genes: dict[tuple[str, str], dict[str, str]] = {}
    with gzip.open(RAW_UNIPROT, "rt", encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            rows = build_row(record)
            if not rows:
                continue
            qualifying = [row for row in rows if row["has_cell_membrane"] and row["has_transmembrane"]]
            if not qualifying:
                continue
            species = _norm_species(rows[0]["organism"])
            gene_name = str(rows[0]["gene_name"] or "")
            if not species or not gene_name:
                continue
            key = (species, gene_name)
            isoform_ids = sorted({str(row["isoform_id"]) for row in qualifying if row.get("isoform_id")})
            accessions = sorted({iso.split("-", 1)[0] for iso in isoform_ids if iso})
            genes[key] = {
                "species": species,
                "gene_name": gene_name,
                "qualifying_isoform_count": str(len(isoform_ids)),
                "qualifying_isoform_ids": ";".join(isoform_ids),
                "qualifying_primary_accessions": ";".join(accessions),
            }
    return genes


def main() -> int:
    strict_surface_genes = _strict_surface_gene_map()
    with STRICT_SURFACE_OUT.open("w", encoding="utf-8", newline="") as gene_handle:
        fieldnames = [
            "species",
            "gene_name",
            "qualifying_isoform_count",
            "qualifying_isoform_ids",
            "qualifying_primary_accessions",
        ]
        writer = csv.DictWriter(gene_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for _, row in sorted(strict_surface_genes.items()):
            writer.writerow(row)

    rows = 0
    with PAIR_IN.open() as in_handle, PAIR_OUT.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or []) + [
            "is_high_confidence",
            "is_strict_surface_gene",
            "is_membrane_pair",
            "is_surface_pair",
            "is_kinase_pair",
            "is_tf_pair",
            "is_clinvar_pair",
            "is_uniprot_preserved_seed",
            "is_appris_contrast_seed",
        ]
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            gene_has_appris = row.get("gene_has_appris_principal_1", "0") == "1"
            ref_has_uniprot = int(row.get("reference_uniprot_links", "0") or 0) > 0
            alt_has_uniprot = int(row.get("alternative_uniprot_links", "0") or 0) > 0
            clinvar = int(row.get("gene_clinvar_pathogenic_splice_count", "0") or 0) > 0
            strict_surface = (row.get("species", ""), row.get("gene_name", "")) in strict_surface_genes
            membrane = strict_surface
            surface = strict_surface
            kinase = row.get("family_class", "") == "kinase"
            tf = row.get("family_class", "") == "transcription_factor"
            weak = row.get("weak_label", "")
            row["is_high_confidence"] = str(int(gene_has_appris and ref_has_uniprot and alt_has_uniprot))
            row["is_strict_surface_gene"] = str(int(strict_surface))
            row["is_membrane_pair"] = str(int(membrane))
            row["is_surface_pair"] = str(int(surface))
            row["is_kinase_pair"] = str(int(kinase))
            row["is_tf_pair"] = str(int(tf))
            row["is_clinvar_pair"] = str(int(clinvar))
            row["is_uniprot_preserved_seed"] = str(int(weak == "likely_preserved"))
            row["is_appris_contrast_seed"] = str(int(weak == "reference_preferred"))
            writer.writerow(row)
            rows += 1

    print(f"[ok] Wrote {STRICT_SURFACE_OUT} rows={len(strict_surface_genes)}")
    print(f"[ok] Wrote {PAIR_OUT} rows={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
