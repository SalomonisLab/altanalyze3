#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path


INTERIM = Path(__file__).resolve().parents[1] / "data" / "interim"


def _base_id(identifier: str) -> str:
    return identifier.split(".", 1)[0] if identifier else ""


def _load_appris_by_transcript() -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    with (INTERIM / "appris_principal_isoforms.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[(row["species"], _base_id(row["transcript_id"]))] = row
    return result


def _load_mane_by_transcript() -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    mane_path = INTERIM / "mane_transcripts.tsv"
    if not mane_path.exists():
        return result
    with mane_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[(row["species"], _base_id(row["transcript_id"]))] = row
    return result


def _load_uniprot_mapping() -> dict[tuple[str, str], Counter]:
    result: dict[tuple[str, str], Counter] = defaultdict(Counter)
    with (INTERIM / "uniprot_gencode_map.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            species = row["gencode_species"]
            transcript_id = row["gencode_transcript_id"]
            if not transcript_id:
                continue
            key = (species, transcript_id)
            result[key]["n_uniprot_links"] += 1
    return result


def _load_protein_priors() -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with (INTERIM / "uniprot_functional_priors.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[row["primary_accession"]] = row
    return result


def _load_region_priors() -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with (INTERIM / "uniprot_region_priors.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[row["primary_accession"]] = row
    return result


def _load_mapping_accession_by_transcript() -> dict[tuple[str, str], set[str]]:
    result: dict[tuple[str, str], set[str]] = defaultdict(set)
    with (INTERIM / "uniprot_gencode_map.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            transcript_id = row["gencode_transcript_id"]
            if transcript_id:
                result[(row["gencode_species"], transcript_id)].add(row["primary_accession"])
    return result


def _load_sequence_features() -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    with (INTERIM / "gencode_protein_sequence_features.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[(row["species"], row["transcript_id"])] = row
    return result


def main() -> int:
    appris = _load_appris_by_transcript()
    mane = _load_mane_by_transcript()
    uniprot_link_counts = _load_uniprot_mapping()
    accession_by_transcript = _load_mapping_accession_by_transcript()
    protein_priors = _load_protein_priors()
    region_priors = _load_region_priors()
    sequence_features = _load_sequence_features()

    out_path = INTERIM / "transcript_supervision_catalog.tsv"
    with (INTERIM / "gencode_transcript_reference.tsv").open() as handle, \
        out_path.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = [
            "species",
            "gene_id",
            "gene_name",
            "transcript_id",
            "transcript_name",
            "protein_id",
            "protein_length",
            "transcript_type",
            "mane_status",
            "is_mane_select",
            "appris_tag",
            "appris_score",
            "is_appris_principal",
            "n_uniprot_links",
            "is_membrane_protein",
            "has_signal_peptide",
            "is_kinase",
            "is_transcription_factor",
            "interaction_partner_count",
            "tm_feature_count",
            "signal_feature_count",
            "signal_max_end",
            "extracellular_topology_aa",
            "cytoplasmic_topology_aa",
            "glycosylation_count",
            "disulfide_count",
            "phospho_feature_count",
            "ppi_binding_feature_count",
            "dna_binding_feature_count",
            "dna_region_feature_count",
            "zinc_finger_feature_count",
            "domain_count",
            "motif_count",
            "predicted_tm_count",
            "predicted_tm_total_span",
            "predicted_tm_max_hydropathy_19",
            "predicted_signal_candidate",
            "predicted_signal_score",
            "glyco_motif_count",
            "cysteine_count",
            "cysteine_fraction"
        ]
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        rows = 0
        for row in reader:
            species = row["species"]
            transcript_id = row["transcript_id"]
            appris_row = appris.get((species, _base_id(transcript_id)), {})
            mane_row = mane.get((species, _base_id(transcript_id)), {})
            accessions = accession_by_transcript.get((species, transcript_id), set())
            membrane = 0
            signal = 0
            kinase = 0
            tf = 0
            interaction_partner_count = 0
            tm_feature_count = 0
            signal_feature_count = 0
            signal_max_end = 0
            extracellular_topology_aa = 0
            cytoplasmic_topology_aa = 0
            glycosylation_count = 0
            disulfide_count = 0
            phospho_feature_count = 0
            ppi_binding_feature_count = 0
            dna_binding_feature_count = 0
            dna_region_feature_count = 0
            zinc_finger_feature_count = 0
            domain_count = 0
            motif_count = 0
            for acc in accessions:
                prior = protein_priors.get(acc, {})
                region = region_priors.get(acc, {})
                membrane = max(membrane, int(prior.get("is_membrane_protein", 0)))
                signal = max(signal, int(prior.get("has_signal_peptide", 0)))
                kinase = max(kinase, int(prior.get("is_kinase", 0)))
                tf = max(tf, int(prior.get("is_transcription_factor", 0)))
                interaction_partner_count = max(interaction_partner_count, int(region.get("interaction_partner_count", 0)))
                tm_feature_count = max(tm_feature_count, int(region.get("tm_count", 0)))
                signal_feature_count = max(signal_feature_count, int(region.get("signal_count", 0)))
                signal_max_end = max(signal_max_end, int(region.get("signal_max_end", 0)))
                extracellular_topology_aa = max(extracellular_topology_aa, int(region.get("extracellular_topology_aa", 0)))
                cytoplasmic_topology_aa = max(cytoplasmic_topology_aa, int(region.get("cytoplasmic_topology_aa", 0)))
                glycosylation_count = max(glycosylation_count, int(region.get("glycosylation_count", 0)))
                disulfide_count = max(disulfide_count, int(region.get("disulfide_count", 0)))
                phospho_feature_count = max(phospho_feature_count, int(region.get("phospho_feature_count", 0)))
                ppi_binding_feature_count = max(ppi_binding_feature_count, int(region.get("ppi_binding_feature_count", 0)))
                dna_binding_feature_count = max(dna_binding_feature_count, int(region.get("dna_binding_feature_count", 0)))
                dna_region_feature_count = max(dna_region_feature_count, int(region.get("dna_region_feature_count", 0)))
                zinc_finger_feature_count = max(zinc_finger_feature_count, int(region.get("zinc_finger_feature_count", 0)))
                domain_count = max(domain_count, int(region.get("domain_count", 0)))
                motif_count = max(motif_count, int(region.get("motif_count", 0)))
            seq = sequence_features.get((species, transcript_id), {})
            writer.writerow(
                {
                    "species": species,
                    "gene_id": row["gene_id"],
                    "gene_name": row["gene_name"],
                    "transcript_id": transcript_id,
                    "transcript_name": row["transcript_name"],
                    "protein_id": row.get("protein_id", ""),
                    "protein_length": row.get("protein_length", ""),
                    "transcript_type": row.get("transcript_type", ""),
                    "mane_status": mane_row.get("mane_status", ""),
                    "is_mane_select": int(mane_row.get("mane_status", "") == "MANE Select"),
                    "appris_tag": appris_row.get("principal_tag", ""),
                    "appris_score": appris_row.get("principal_score", ""),
                    "is_appris_principal": int((appris_row.get("principal_tag", "") or "").startswith("PRINCIPAL")),
                    "n_uniprot_links": uniprot_link_counts.get((species, transcript_id), Counter()).get("n_uniprot_links", 0),
                    "is_membrane_protein": membrane,
                    "has_signal_peptide": signal,
                    "is_kinase": kinase,
                    "is_transcription_factor": tf,
                    "interaction_partner_count": interaction_partner_count,
                    "tm_feature_count": tm_feature_count,
                    "signal_feature_count": signal_feature_count,
                    "signal_max_end": signal_max_end,
                    "extracellular_topology_aa": extracellular_topology_aa,
                    "cytoplasmic_topology_aa": cytoplasmic_topology_aa,
                    "glycosylation_count": glycosylation_count,
                    "disulfide_count": disulfide_count,
                    "phospho_feature_count": phospho_feature_count,
                    "ppi_binding_feature_count": ppi_binding_feature_count,
                    "dna_binding_feature_count": dna_binding_feature_count,
                    "dna_region_feature_count": dna_region_feature_count,
                    "zinc_finger_feature_count": zinc_finger_feature_count,
                    "domain_count": domain_count,
                    "motif_count": motif_count,
                    "predicted_tm_count": seq.get("predicted_tm_count", 0),
                    "predicted_tm_total_span": seq.get("predicted_tm_total_span", 0),
                    "predicted_tm_max_hydropathy_19": seq.get("predicted_tm_max_hydropathy_19", 0),
                    "predicted_signal_candidate": seq.get("predicted_signal_candidate", 0),
                    "predicted_signal_score": seq.get("predicted_signal_score", 0),
                    "glyco_motif_count": seq.get("glyco_motif_count", 0),
                    "cysteine_count": seq.get("cysteine_count", 0),
                    "cysteine_fraction": seq.get("cysteine_fraction", 0),
                }
            )
            rows += 1
    print(f"[ok] Wrote {out_path} rows={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
