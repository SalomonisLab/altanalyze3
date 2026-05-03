#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path


INTERIM = Path(__file__).resolve().parents[1] / "data" / "interim"
ENTRY_FILE = INTERIM / "uniprot_entries.tsv"
FEATURE_FILE = INTERIM / "uniprot_features.tsv"
OUT_FILE = INTERIM / "uniprot_region_priors.tsv"

SMALL_AA = set("AGSCTV")


def _parse_int(value: str) -> int:
    if value in {"", None}:
        return 0
    try:
        return int(float(value))
    except ValueError:
        return 0


def _span_length(begin: str, end: str) -> int:
    start = _parse_int(begin)
    stop = _parse_int(end)
    if start <= 0 or stop <= 0 or stop < start:
        return 0
    return stop - start + 1


def _contains_any(text: str, terms: list[str]) -> bool:
    lowered = (text or "").lower()
    return any(term in lowered for term in terms)


def _load_interaction_partner_counts() -> dict[str, int]:
    counts: dict[str, int] = {}
    with ENTRY_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            partners = {x for x in (row.get("interaction_partners", "") or "").split(";") if x}
            counts[row["primary_accession"]] = len(partners)
    return counts


def main() -> int:
    interaction_counts = _load_interaction_partner_counts()
    stats: dict[str, Counter] = defaultdict(Counter)
    meta: dict[str, dict[str, str]] = {}

    with ENTRY_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            meta[row["primary_accession"]] = {
                "gene_name": row.get("gene_name", ""),
                "organism": row.get("organism", ""),
                "protein_name": row.get("protein_name", ""),
            }

    with FEATURE_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = row["primary_accession"]
            feature_type = row.get("feature_type", "")
            description = row.get("description", "") or ""
            description_l = description.lower()
            span = _span_length(row.get("begin", ""), row.get("end", ""))
            bucket = stats[accession]
            bucket[f"type::{feature_type}"] += 1
            bucket[f"span::{feature_type}"] += span

            if feature_type == "Signal":
                bucket["signal_count"] += 1
                bucket["signal_max_end"] = max(bucket["signal_max_end"], _parse_int(row.get("end", "")))
            elif feature_type == "Transmembrane":
                bucket["tm_count"] += 1
                bucket["tm_total_span"] += span
            elif feature_type == "Intramembrane":
                bucket["intramembrane_count"] += 1
                bucket["intramembrane_total_span"] += span
            elif feature_type == "Topological domain":
                if _contains_any(description, ["extracellular", "lumenal", "luminal", "external"]):
                    bucket["extracellular_topology_count"] += 1
                    bucket["extracellular_topology_aa"] += span
                if _contains_any(description, ["cytoplasmic", "cytosolic", "intracellular"]):
                    bucket["cytoplasmic_topology_count"] += 1
                    bucket["cytoplasmic_topology_aa"] += span
                if _contains_any(description, ["perinuclear", "nuclear"]):
                    bucket["nuclear_topology_count"] += 1
            elif feature_type == "Glycosylation":
                bucket["glycosylation_count"] += 1
            elif feature_type == "Disulfide bond":
                bucket["disulfide_count"] += 1
            elif feature_type == "Modified residue":
                bucket["modified_residue_count"] += 1
                if _contains_any(description, ["phospho", "phosphoserine", "phosphothreonine", "phosphotyrosine"]):
                    bucket["phospho_feature_count"] += 1
            elif feature_type in {"Binding site", "Site", "Region", "Domain", "Motif"}:
                if _contains_any(description, ["interact", "binding", "dimer", "partner", "sh2", "sh3", "pdz", "ww domain"]):
                    bucket["ppi_binding_feature_count"] += 1
                if _contains_any(description, ["phospho", "phosphoryl"]):
                    bucket["phospho_interaction_feature_count"] += 1
                if _contains_any(description, ["dna-binding", "dna binding", "homeobox", "helix-loop-helix", "bzip"]):
                    bucket["dna_region_feature_count"] += 1
            elif feature_type == "DNA binding":
                bucket["dna_binding_feature_count"] += 1
            elif feature_type == "Zinc finger":
                bucket["zinc_finger_feature_count"] += 1
            elif feature_type == "Coiled coil":
                bucket["coiled_coil_count"] += 1
            elif feature_type == "Active site":
                bucket["active_site_count"] += 1
            elif feature_type == "Domain":
                bucket["domain_count"] += 1
            elif feature_type == "Motif":
                bucket["motif_count"] += 1

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "primary_accession",
        "gene_name",
        "organism",
        "protein_name",
        "interaction_partner_count",
        "signal_count",
        "signal_max_end",
        "tm_count",
        "tm_total_span",
        "intramembrane_count",
        "intramembrane_total_span",
        "extracellular_topology_count",
        "extracellular_topology_aa",
        "cytoplasmic_topology_count",
        "cytoplasmic_topology_aa",
        "nuclear_topology_count",
        "glycosylation_count",
        "disulfide_count",
        "modified_residue_count",
        "phospho_feature_count",
        "ppi_binding_feature_count",
        "phospho_interaction_feature_count",
        "dna_binding_feature_count",
        "dna_region_feature_count",
        "zinc_finger_feature_count",
        "coiled_coil_count",
        "active_site_count",
        "domain_count",
        "motif_count",
    ]
    with OUT_FILE.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for accession in sorted(meta):
            bucket = stats.get(accession, Counter())
            info = meta[accession]
            writer.writerow(
                {
                    "primary_accession": accession,
                    "gene_name": info.get("gene_name", ""),
                    "organism": info.get("organism", ""),
                    "protein_name": info.get("protein_name", ""),
                    "interaction_partner_count": interaction_counts.get(accession, 0),
                    "signal_count": bucket.get("signal_count", 0),
                    "signal_max_end": bucket.get("signal_max_end", 0),
                    "tm_count": bucket.get("tm_count", 0),
                    "tm_total_span": bucket.get("tm_total_span", 0),
                    "intramembrane_count": bucket.get("intramembrane_count", 0),
                    "intramembrane_total_span": bucket.get("intramembrane_total_span", 0),
                    "extracellular_topology_count": bucket.get("extracellular_topology_count", 0),
                    "extracellular_topology_aa": bucket.get("extracellular_topology_aa", 0),
                    "cytoplasmic_topology_count": bucket.get("cytoplasmic_topology_count", 0),
                    "cytoplasmic_topology_aa": bucket.get("cytoplasmic_topology_aa", 0),
                    "nuclear_topology_count": bucket.get("nuclear_topology_count", 0),
                    "glycosylation_count": bucket.get("glycosylation_count", 0),
                    "disulfide_count": bucket.get("disulfide_count", 0),
                    "modified_residue_count": bucket.get("modified_residue_count", 0),
                    "phospho_feature_count": bucket.get("phospho_feature_count", 0),
                    "ppi_binding_feature_count": bucket.get("ppi_binding_feature_count", 0),
                    "phospho_interaction_feature_count": bucket.get("phospho_interaction_feature_count", 0),
                    "dna_binding_feature_count": bucket.get("dna_binding_feature_count", 0),
                    "dna_region_feature_count": bucket.get("dna_region_feature_count", 0),
                    "zinc_finger_feature_count": bucket.get("zinc_finger_feature_count", 0),
                    "coiled_coil_count": bucket.get("coiled_coil_count", 0),
                    "active_site_count": bucket.get("active_site_count", 0),
                    "domain_count": bucket.get("domain_count", 0),
                    "motif_count": bucket.get("motif_count", 0),
                }
            )
    print(f"[ok] Wrote {OUT_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
