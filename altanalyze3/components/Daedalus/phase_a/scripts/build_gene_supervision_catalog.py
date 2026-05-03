#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path


INTERIM = Path(__file__).resolve().parents[1] / "data" / "interim"


def _base_id(value: str) -> str:
    return value.split(".", 1)[0] if value else value


def _load_transcript_counts() -> dict[tuple[str, str, str], Counter]:
    counts: dict[tuple[str, str, str], Counter] = defaultdict(Counter)
    with (INTERIM / "gencode_transcript_reference.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["species"], row["gene_id"], row["gene_name"])
            counts[key]["n_transcripts"] += 1
            if row.get("protein_id"):
                counts[key]["n_protein_coding_transcripts"] += 1
    return counts


def _load_hpa() -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with (INTERIM / "hpa_localization.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[_base_id(row["ensembl_gene_id"])] = row
    return result


def _load_clinvar_by_symbol() -> Counter:
    counts: Counter = Counter()
    with (INTERIM / "clinvar_splice_pathogenic.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("gene_symbol"):
                counts[row["gene_symbol"]] += 1
    return counts


def _load_uniprot_by_gene() -> dict[tuple[str, str], Counter]:
    result: dict[tuple[str, str], Counter] = defaultdict(Counter)
    with (INTERIM / "uniprot_functional_priors.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            organism = row["organism"]
            gene_name = row["gene_name"]
            if not gene_name:
                continue
            key = (organism, gene_name)
            result[key]["n_uniprot_entries"] += 1
            result[key]["is_membrane"] = max(result[key]["is_membrane"], int(row["is_membrane_protein"]))
            result[key]["has_signal"] = max(result[key]["has_signal"], int(row["has_signal_peptide"]))
            result[key]["is_kinase"] = max(result[key]["is_kinase"], int(row["is_kinase"]))
            result[key]["is_tf"] = max(result[key]["is_tf"], int(row["is_transcription_factor"]))
            result[key]["interaction_partner_count"] = max(
                result[key]["interaction_partner_count"],
                int(row["interaction_partner_count"])
            )
    return result


def _load_uniprot_region_by_gene() -> dict[tuple[str, str], Counter]:
    result: dict[tuple[str, str], Counter] = defaultdict(Counter)
    with (INTERIM / "uniprot_region_priors.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            organism = row["organism"]
            gene_name = row["gene_name"]
            if not gene_name:
                continue
            key = (organism, gene_name)
            for field in [
                "tm_count",
                "signal_count",
                "extracellular_topology_aa",
                "cytoplasmic_topology_aa",
                "glycosylation_count",
                "disulfide_count",
                "phospho_feature_count",
                "ppi_binding_feature_count",
                "phospho_interaction_feature_count",
                "dna_binding_feature_count",
                "dna_region_feature_count",
                "zinc_finger_feature_count",
                "domain_count",
                "motif_count",
            ]:
                result[key][field] = max(result[key][field], int(row.get(field, 0) or 0))
    return result


def _load_biogrid_by_gene() -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    path = INTERIM / "biogrid_gene_interactions.tsv"
    if not path.exists():
        return result
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[(row["species"], row["gene_name"])] = row
    return result


def _load_uniprot_isoform_counts() -> dict[tuple[str, str], int]:
    counts: Counter = Counter()
    with (INTERIM / "uniprot_isoforms.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("gene_name"):
                counts[(row["organism"], row["gene_name"])] += 1
    return counts


def _load_appris_by_gene() -> dict[tuple[str, str], Counter]:
    counts: dict[tuple[str, str], Counter] = defaultdict(Counter)
    with (INTERIM / "appris_principal_isoforms.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["species"], _base_id(row["gene_id"]))
            counts[key]["n_appris_transcripts"] += 1
            principal_tag = row.get("principal_tag", "")
            if principal_tag.startswith("PRINCIPAL"):
                counts[key]["n_appris_principal_transcripts"] += 1
                if principal_tag == "PRINCIPAL:1":
                    counts[key]["has_principal_1"] = 1
    return counts


def main() -> int:
    transcript_counts = _load_transcript_counts()
    hpa = _load_hpa()
    clinvar_counts = _load_clinvar_by_symbol()
    uniprot_gene = _load_uniprot_by_gene()
    uniprot_region_gene = _load_uniprot_region_by_gene()
    uniprot_isoforms = _load_uniprot_isoform_counts()
    appris_gene = _load_appris_by_gene()
    biogrid_gene = _load_biogrid_by_gene()

    out_path = INTERIM / "gene_supervision_catalog.tsv"
    with out_path.open("w", encoding="utf-8", newline="") as out_handle:
        fieldnames = [
            "species",
            "gencode_gene_id",
            "gene_name",
            "n_transcripts",
            "n_protein_coding_transcripts",
            "n_uniprot_entries",
            "n_uniprot_isoforms",
            "is_membrane_protein",
            "has_signal_peptide",
            "is_kinase",
            "is_transcription_factor",
            "interaction_partner_count",
            "biogrid_partner_count",
            "biogrid_interaction_count",
            "biogrid_physical_partner_count",
            "biogrid_physical_interaction_count",
            "tm_feature_count",
            "signal_feature_count",
            "extracellular_topology_aa",
            "cytoplasmic_topology_aa",
            "glycosylation_count",
            "disulfide_count",
            "phospho_feature_count",
            "ppi_binding_feature_count",
            "phospho_interaction_feature_count",
            "dna_binding_feature_count",
            "dna_region_feature_count",
            "zinc_finger_feature_count",
            "domain_count",
            "motif_count",
            "n_appris_transcripts",
            "n_appris_principal_transcripts",
            "has_appris_principal_1",
            "hpa_reliability",
            "hpa_all_locations",
            "hpa_has_extracellular_annotation",
            "clinvar_pathogenic_splice_count"
        ]
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        rows = 0
        for (species, gene_id, gene_name), counts in transcript_counts.items():
            organism = "Homo sapiens" if species == "human" else "Mus musculus"
            uniprot_key = (organism, gene_name)
            base_gene_id = _base_id(gene_id)
            appris_key = (species, base_gene_id)
            uniprot_meta = uniprot_gene.get(uniprot_key, Counter())
            uniprot_region_meta = uniprot_region_gene.get(uniprot_key, Counter())
            appris_meta = appris_gene.get(appris_key, Counter())
            biogrid_meta = biogrid_gene.get((species, gene_name), {})
            hpa_row = hpa.get(base_gene_id, {})
            writer.writerow(
                {
                    "species": species,
                    "gencode_gene_id": gene_id,
                    "gene_name": gene_name,
                    "n_transcripts": counts["n_transcripts"],
                    "n_protein_coding_transcripts": counts["n_protein_coding_transcripts"],
                    "n_uniprot_entries": uniprot_meta.get("n_uniprot_entries", 0),
                    "n_uniprot_isoforms": uniprot_isoforms.get(uniprot_key, 0),
                    "is_membrane_protein": uniprot_meta.get("is_membrane", 0),
                    "has_signal_peptide": uniprot_meta.get("has_signal", 0),
                    "is_kinase": uniprot_meta.get("is_kinase", 0),
                    "is_transcription_factor": uniprot_meta.get("is_tf", 0),
                    "interaction_partner_count": uniprot_meta.get("interaction_partner_count", 0),
                    "biogrid_partner_count": biogrid_meta.get("biogrid_partner_count", 0),
                    "biogrid_interaction_count": biogrid_meta.get("biogrid_interaction_count", 0),
                    "biogrid_physical_partner_count": biogrid_meta.get("biogrid_physical_partner_count", 0),
                    "biogrid_physical_interaction_count": biogrid_meta.get("biogrid_physical_interaction_count", 0),
                    "tm_feature_count": uniprot_region_meta.get("tm_count", 0),
                    "signal_feature_count": uniprot_region_meta.get("signal_count", 0),
                    "extracellular_topology_aa": uniprot_region_meta.get("extracellular_topology_aa", 0),
                    "cytoplasmic_topology_aa": uniprot_region_meta.get("cytoplasmic_topology_aa", 0),
                    "glycosylation_count": uniprot_region_meta.get("glycosylation_count", 0),
                    "disulfide_count": uniprot_region_meta.get("disulfide_count", 0),
                    "phospho_feature_count": uniprot_region_meta.get("phospho_feature_count", 0),
                    "ppi_binding_feature_count": uniprot_region_meta.get("ppi_binding_feature_count", 0),
                    "phospho_interaction_feature_count": uniprot_region_meta.get("phospho_interaction_feature_count", 0),
                    "dna_binding_feature_count": uniprot_region_meta.get("dna_binding_feature_count", 0),
                    "dna_region_feature_count": uniprot_region_meta.get("dna_region_feature_count", 0),
                    "zinc_finger_feature_count": uniprot_region_meta.get("zinc_finger_feature_count", 0),
                    "domain_count": uniprot_region_meta.get("domain_count", 0),
                    "motif_count": uniprot_region_meta.get("motif_count", 0),
                    "n_appris_transcripts": appris_meta.get("n_appris_transcripts", 0),
                    "n_appris_principal_transcripts": appris_meta.get("n_appris_principal_transcripts", 0),
                    "has_appris_principal_1": appris_meta.get("has_principal_1", 0),
                    "hpa_reliability": hpa_row.get("reliability", ""),
                    "hpa_all_locations": hpa_row.get("all_locations", ""),
                    "hpa_has_extracellular_annotation": hpa_row.get("has_extracellular_annotation", ""),
                    "clinvar_pathogenic_splice_count": clinvar_counts.get(gene_name, 0)
                }
            )
            rows += 1

    print(f"[ok] Wrote {out_path} rows={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
