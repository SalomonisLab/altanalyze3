#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path


INPUT_FILE = Path(__file__).resolve().parents[1] / "data" / "raw" / "uniprot_reviewed_human_mouse.jsonl.gz"
OUT_DIR = Path(__file__).resolve().parents[1] / "data" / "interim"


def _first_gene_name(genes: list[dict]) -> str:
    for gene in genes or []:
        gene_name = gene.get("geneName")
        if gene_name and gene_name.get("value"):
            return gene_name["value"]
    return ""


def _extract_ensembl_ids(xrefs: list[dict]) -> str:
    values: list[str] = []
    for xref in xrefs or []:
        if xref.get("database") != "Ensembl":
            continue
        values.append(xref.get("id", ""))
        for prop in xref.get("properties", []) or []:
            if prop.get("value"):
                values.append(prop["value"])
    return ";".join(sorted({v for v in values if v}))


def _protein_name(record: dict) -> str:
    protein_description = record.get("proteinDescription") or {}
    recommended = protein_description.get("recommendedName") or {}
    full_name = recommended.get("fullName") or {}
    if full_name.get("value"):
        return full_name["value"]
    submission_names = protein_description.get("submissionNames") or []
    for item in submission_names:
        full_name = item.get("fullName") or {}
        if full_name.get("value"):
            return full_name["value"]
    return ""


def _keywords(record: dict) -> str:
    values: list[str] = []
    for keyword in record.get("keywords", []) or []:
        if keyword.get("name"):
            values.append(keyword["name"])
    return ";".join(sorted(set(values)))


def _subcellular_locations(comments: list[dict]) -> str:
    values: list[str] = []
    for comment in comments or []:
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for location in comment.get("subcellularLocations", []) or []:
            loc = location.get("location", {})
            if loc.get("value"):
                values.append(loc["value"])
            topology = location.get("topology", {})
            if topology.get("value"):
                values.append(f"topology:{topology['value']}")
    return ";".join(sorted(set(values)))


def _interaction_partners(comments: list[dict]) -> str:
    values: list[str] = []
    for comment in comments or []:
        if comment.get("commentType") != "INTERACTION":
            continue
        for interaction in comment.get("interactions", []) or []:
            interactant = interaction.get("interactantTwo", {})
            if interactant.get("uniProtKBAccession"):
                values.append(interactant["uniProtKBAccession"])
            elif interactant.get("geneName"):
                values.append(interactant["geneName"])
    return ";".join(sorted(set(values)))


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    isoform_out = OUT_DIR / "uniprot_isoforms.tsv"
    feature_out = OUT_DIR / "uniprot_features.tsv"
    entry_out = OUT_DIR / "uniprot_entries.tsv"

    isoform_fields = [
        "primary_accession",
        "uniProtkb_id",
        "organism",
        "gene_name",
        "protein_name",
        "ensembl_ids",
        "isoform_id",
        "isoform_name",
        "isoform_synonyms",
        "isoform_sequence_status",
        "alternative_events",
        "subcellular_locations",
        "interaction_partners"
    ]
    feature_fields = [
        "primary_accession",
        "gene_name",
        "feature_type",
        "description",
        "begin",
        "end"
    ]
    entry_fields = [
        "primary_accession",
        "uniProtkb_id",
        "organism",
        "gene_name",
        "protein_name",
        "ensembl_ids",
        "subcellular_locations",
        "interaction_partners",
        "keywords",
        "keyword_count",
        "feature_count"
    ]

    with gzip.open(INPUT_FILE, "rt", encoding="utf-8") as handle, \
        entry_out.open("w", encoding="utf-8", newline="") as entry_handle, \
        isoform_out.open("w", encoding="utf-8", newline="") as isoform_handle, \
        feature_out.open("w", encoding="utf-8", newline="") as feature_handle:

        entry_writer = csv.DictWriter(entry_handle, fieldnames=entry_fields, delimiter="\t")
        isoform_writer = csv.DictWriter(isoform_handle, fieldnames=isoform_fields, delimiter="\t")
        feature_writer = csv.DictWriter(feature_handle, fieldnames=feature_fields, delimiter="\t")
        entry_writer.writeheader()
        isoform_writer.writeheader()
        feature_writer.writeheader()

        n_entries = 0
        n_isoforms = 0
        n_features = 0
        for line in handle:
            record = json.loads(line)
            primary_accession = record.get("primaryAccession", "")
            gene_name = _first_gene_name(record.get("genes", []))
            organism = (record.get("organism") or {}).get("scientificName", "")
            xrefs = record.get("uniProtKBCrossReferences", []) or []
            comments = record.get("comments", []) or []
            features = record.get("features", []) or []
            ensembl_ids = _extract_ensembl_ids(xrefs)
            protein_name = _protein_name(record)
            keywords = _keywords(record)
            subcellular_locations = _subcellular_locations(comments)
            interaction_partners = _interaction_partners(comments)

            entry_writer.writerow(
                {
                    "primary_accession": primary_accession,
                    "uniProtkb_id": record.get("uniProtkbId", ""),
                    "organism": organism,
                    "gene_name": gene_name,
                    "protein_name": protein_name,
                    "ensembl_ids": ensembl_ids,
                    "subcellular_locations": subcellular_locations,
                    "interaction_partners": interaction_partners,
                    "keywords": keywords,
                    "keyword_count": len(record.get("keywords", []) or []),
                    "feature_count": len(features)
                }
            )
            n_entries += 1

            for comment in comments:
                if comment.get("commentType") != "ALTERNATIVE PRODUCTS":
                    continue
                events = ";".join(comment.get("events", []) or [])
                for isoform in comment.get("isoforms", []) or []:
                    isoform_writer.writerow(
                        {
                            "primary_accession": primary_accession,
                            "uniProtkb_id": record.get("uniProtkbId", ""),
                            "organism": organism,
                            "gene_name": gene_name,
                            "protein_name": protein_name,
                            "ensembl_ids": ensembl_ids,
                            "isoform_id": ";".join(isoform.get("isoformIds", []) or []),
                            "isoform_name": (isoform.get("name") or {}).get("value", ""),
                            "isoform_synonyms": ";".join(
                                synonym.get("value", "")
                                for synonym in isoform.get("synonyms", []) or []
                                if synonym.get("value")
                            ),
                            "isoform_sequence_status": isoform.get("isoformSequenceStatus", ""),
                            "alternative_events": events,
                            "subcellular_locations": subcellular_locations,
                            "interaction_partners": interaction_partners
                        }
                    )
                    n_isoforms += 1

            for feature in features:
                location = feature.get("location", {}) or {}
                start = (location.get("start") or {}).get("value", "")
                end = (location.get("end") or {}).get("value", "")
                feature_writer.writerow(
                    {
                        "primary_accession": primary_accession,
                        "gene_name": gene_name,
                        "feature_type": feature.get("type", ""),
                        "description": feature.get("description", ""),
                        "begin": start,
                        "end": end
                    }
                )
                n_features += 1

    print(f"[ok] Wrote {entry_out} entries={n_entries}")
    print(f"[ok] Wrote {isoform_out} isoforms={n_isoforms}")
    print(f"[ok] Wrote {feature_out} features={n_features}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
