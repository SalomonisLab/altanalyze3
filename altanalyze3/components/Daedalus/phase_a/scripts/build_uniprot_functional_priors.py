#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path


INTERIM = Path(__file__).resolve().parents[1] / "data" / "interim"
ENTRY_FILE = INTERIM / "uniprot_entries.tsv"
FEATURE_FILE = INTERIM / "uniprot_features.tsv"
OUT_FILE = INTERIM / "uniprot_functional_priors.tsv"


def _load_feature_stats() -> dict[str, Counter]:
    stats: dict[str, Counter] = defaultdict(Counter)
    with FEATURE_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = row["primary_accession"]
            feature_type = row["feature_type"]
            description = (row.get("description") or "").lower()
            stats[accession][f"type::{feature_type}"] += 1
            if "kinase" in description:
                stats[accession]["desc::kinase"] += 1
            if "dna-binding" in description or "dna binding" in description:
                stats[accession]["desc::dna_binding"] += 1
            if "transcription" in description:
                stats[accession]["desc::transcription"] += 1
    return stats


def _contains_any(text: str, terms: list[str]) -> bool:
    lowered = text.lower()
    return any(term in lowered for term in terms)


def main() -> int:
    feature_stats = _load_feature_stats()
    with ENTRY_FILE.open() as entry_handle, OUT_FILE.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.DictReader(entry_handle, delimiter="\t")
        fieldnames = [
            "primary_accession",
            "gene_name",
            "organism",
            "protein_name",
            "is_membrane_protein",
            "has_signal_peptide",
            "tm_feature_count",
            "topological_domain_count",
            "is_kinase",
            "is_transcription_factor",
            "interaction_partner_count",
            "subcellular_locations",
            "keywords"
        ]
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            accession = row["primary_accession"]
            stats = feature_stats.get(accession, Counter())
            protein_name = row.get("protein_name", "")
            keywords = row.get("keywords", "")
            subcellular = row.get("subcellular_locations", "")
            interaction_partners = [x for x in (row.get("interaction_partners", "")).split(";") if x]

            tm_count = stats["type::Transmembrane"]
            topo_count = stats["type::Topological domain"]
            has_signal = stats["type::Signal"] > 0 or "signal" in protein_name.lower()
            is_membrane = (
                tm_count > 0
                or "membrane" in subcellular.lower()
                or "cell membrane" in keywords.lower()
                or "cell junction" in subcellular.lower()
            )
            is_kinase = (
                stats["desc::kinase"] > 0
                or _contains_any(protein_name, ["kinase"])
                or _contains_any(keywords, ["kinase", "transferase"])
            )
            is_tf = (
                stats["type::DNA binding"] > 0
                or stats["type::Zinc finger"] > 0
                or stats["desc::dna_binding"] > 0
                or _contains_any(keywords, ["transcription", "dna-binding", "homeobox", "zinc-finger"])
                or _contains_any(protein_name, ["transcription factor"])
            )

            writer.writerow(
                {
                    "primary_accession": accession,
                    "gene_name": row.get("gene_name", ""),
                    "organism": row.get("organism", ""),
                    "protein_name": protein_name,
                    "is_membrane_protein": int(is_membrane),
                    "has_signal_peptide": int(has_signal),
                    "tm_feature_count": tm_count,
                    "topological_domain_count": topo_count,
                    "is_kinase": int(is_kinase),
                    "is_transcription_factor": int(is_tf),
                    "interaction_partner_count": len(set(interaction_partners)),
                    "subcellular_locations": subcellular,
                    "keywords": keywords
                }
            )
    print(f"[ok] Wrote {OUT_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
