#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


INTERIM = Path(__file__).resolve().parents[1] / "data" / "interim"
PAIR_OUT = INTERIM / "isoform_pair_candidates.tsv"
GENE_SET_OUT = INTERIM / "benchmark_gene_sets.tsv"


def _base_id(value: str) -> str:
    return value.split(".", 1)[0] if value else ""


def _parse_int(value: str) -> int:
    if not value:
        return 0
    try:
        return int(float(value))
    except ValueError:
        return 0


def _parse_appris_rank(tag: str) -> tuple[int, int]:
    if not tag:
        return (9, 9)
    if tag.startswith("PRINCIPAL:"):
        try:
            return (0, int(tag.split(":", 1)[1]))
        except ValueError:
            return (0, 9)
    if tag.startswith("ALTERNATIVE:"):
        try:
            return (1, int(tag.split(":", 1)[1]))
        except ValueError:
            return (1, 9)
    return (9, 9)


def _reference_source(is_mane_select: bool, tag: str, protein_length: int, n_uniprot_links: int) -> str:
    if is_mane_select:
        return "MANE_SELECT"
    if tag == "PRINCIPAL:1":
        return "APPRIS_PRINCIPAL_1"
    if tag.startswith("PRINCIPAL"):
        return "APPRIS_PRINCIPAL"
    if n_uniprot_links > 0 and protein_length > 0:
        return "UNIPROT_SUPPORTED_LONGEST"
    if protein_length > 0:
        return "LONGEST_PROTEIN"
    return "LONGEST_TRANSCRIPT"


def _family_class(row: dict[str, str]) -> str:
    if row["is_kinase"] == "1":
        return "kinase"
    if row["is_transcription_factor"] == "1":
        return "transcription_factor"
    if row["is_membrane_protein"] == "1":
        return "membrane"
    return "other"


def _load_gene_meta() -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    with (INTERIM / "gene_supervision_catalog.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            result[(row["species"], _base_id(row["gencode_gene_id"]))] = row
    return result


def _load_transcripts() -> dict[tuple[str, str], list[dict[str, str]]]:
    genes: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    gencode_lengths: dict[tuple[str, str], int] = {}
    with (INTERIM / "gencode_transcript_reference.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gencode_lengths[(row["species"], row["transcript_id"])] = _parse_int(row.get("transcript_length", ""))

    with (INTERIM / "transcript_supervision_catalog.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            transcript_type = row.get("transcript_type", "")
            protein_length = _parse_int(row.get("protein_length", ""))
            transcript_id = row["transcript_id"]
            species = row["species"]
            enriched = dict(row)
            enriched["gene_id_base"] = _base_id(row["gene_id"])
            enriched["protein_length_int"] = str(protein_length)
            enriched["transcript_length_int"] = str(gencode_lengths.get((species, transcript_id), 0))
            enriched["n_uniprot_links_int"] = str(_parse_int(row.get("n_uniprot_links", "")))
            enriched["is_protein_coding_candidate"] = str(
                int(bool(row.get("protein_id")) or transcript_type == "protein_coding" or protein_length > 0)
            )
            genes[(species, enriched["gene_id_base"])].append(enriched)
    return genes


def _load_accessions_by_transcript() -> dict[tuple[str, str], set[str]]:
    result: dict[tuple[str, str], set[str]] = defaultdict(set)
    with (INTERIM / "uniprot_gencode_map.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            transcript_id = row.get("gencode_transcript_id", "")
            if transcript_id:
                result[(row["gencode_species"], transcript_id)].add(row["primary_accession"])
    return result


def main() -> int:
    gene_meta = _load_gene_meta()
    transcripts_by_gene = _load_transcripts()
    accessions_by_transcript = _load_accessions_by_transcript()

    pair_fields = [
        "species",
        "gene_id",
        "gene_name",
        "family_class",
        "reference_source",
        "reference_transcript_id",
        "alternative_transcript_id",
        "reference_appris_tag",
        "alternative_appris_tag",
        "reference_protein_id",
        "alternative_protein_id",
        "reference_protein_length",
        "alternative_protein_length",
        "protein_length_delta",
        "protein_length_ratio",
        "reference_transcript_length",
        "alternative_transcript_length",
        "transcript_length_delta",
        "reference_uniprot_links",
        "alternative_uniprot_links",
        "shared_uniprot_accessions",
        "reference_is_membrane",
        "alternative_is_membrane",
        "reference_has_signal_peptide",
        "alternative_has_signal_peptide",
        "reference_is_kinase",
        "alternative_is_kinase",
        "reference_is_transcription_factor",
        "alternative_is_transcription_factor",
        "reference_interaction_partner_count",
        "alternative_interaction_partner_count",
        "reference_tm_feature_count",
        "alternative_tm_feature_count",
        "tm_feature_count_delta",
        "reference_signal_feature_count",
        "alternative_signal_feature_count",
        "signal_feature_count_delta",
        "reference_signal_max_end",
        "alternative_signal_max_end",
        "reference_extracellular_topology_aa",
        "alternative_extracellular_topology_aa",
        "extracellular_topology_aa_delta",
        "reference_cytoplasmic_topology_aa",
        "alternative_cytoplasmic_topology_aa",
        "cytoplasmic_topology_aa_delta",
        "reference_glycosylation_count",
        "alternative_glycosylation_count",
        "glycosylation_count_delta",
        "reference_disulfide_count",
        "alternative_disulfide_count",
        "disulfide_count_delta",
        "reference_phospho_feature_count",
        "alternative_phospho_feature_count",
        "phospho_feature_count_delta",
        "reference_ppi_binding_feature_count",
        "alternative_ppi_binding_feature_count",
        "ppi_binding_feature_count_delta",
        "reference_dna_binding_feature_count",
        "alternative_dna_binding_feature_count",
        "dna_binding_feature_count_delta",
        "reference_dna_region_feature_count",
        "alternative_dna_region_feature_count",
        "dna_region_feature_count_delta",
        "reference_zinc_finger_feature_count",
        "alternative_zinc_finger_feature_count",
        "zinc_finger_feature_count_delta",
        "reference_domain_count",
        "alternative_domain_count",
        "domain_count_delta",
        "reference_motif_count",
        "alternative_motif_count",
        "motif_count_delta",
        "reference_predicted_tm_count",
        "alternative_predicted_tm_count",
        "predicted_tm_count_delta",
        "reference_predicted_tm_total_span",
        "alternative_predicted_tm_total_span",
        "predicted_tm_total_span_delta",
        "reference_predicted_signal_candidate",
        "alternative_predicted_signal_candidate",
        "reference_predicted_signal_score",
        "alternative_predicted_signal_score",
        "predicted_signal_score_delta",
        "reference_glyco_motif_count",
        "alternative_glyco_motif_count",
        "glyco_motif_count_delta",
        "reference_cysteine_count",
        "alternative_cysteine_count",
        "cysteine_count_delta",
        "gene_has_appris_principal_1",
        "gene_clinvar_pathogenic_splice_count",
        "gene_hpa_has_extracellular_annotation",
        "gene_biogrid_partner_count",
        "gene_biogrid_physical_partner_count",
        "weak_label",
        "weak_label_reason",
    ]
    gene_set_fields = [
        "species",
        "gene_id",
        "gene_name",
        "family_class",
        "n_protein_coding_transcripts",
        "has_appris_principal_1",
        "n_uniprot_entries",
        "n_uniprot_isoforms",
        "clinvar_pathogenic_splice_count",
        "hpa_has_extracellular_annotation",
        "reference_transcript_id",
        "reference_source",
        "n_pair_candidates",
    ]

    pair_rows = 0
    gene_rows = 0
    with PAIR_OUT.open("w", encoding="utf-8", newline="") as pair_handle, \
            GENE_SET_OUT.open("w", encoding="utf-8", newline="") as gene_handle:
        pair_writer = csv.DictWriter(pair_handle, fieldnames=pair_fields, delimiter="\t")
        gene_writer = csv.DictWriter(gene_handle, fieldnames=gene_set_fields, delimiter="\t")
        pair_writer.writeheader()
        gene_writer.writeheader()

        for gene_key, transcripts in sorted(transcripts_by_gene.items()):
            species, gene_id_base = gene_key
            gene_info = gene_meta.get(gene_key)
            if not gene_info:
                continue
            coding = [row for row in transcripts if row["is_protein_coding_candidate"] == "1"]
            if len(coding) < 2:
                continue

            coding.sort(
                key=lambda row: (
                    -_parse_int(row.get("is_mane_select", "0")),
                    _parse_appris_rank(row.get("appris_tag", "")),
                    -_parse_int(row.get("n_uniprot_links_int", "")),
                    -_parse_int(row.get("protein_length_int", "")),
                    -_parse_int(row.get("transcript_length_int", "")),
                    row["transcript_id"],
                )
            )
            reference = coding[0]
            family_class = _family_class(reference)
            reference_source = _reference_source(
                reference.get("is_mane_select", "0") == "1",
                reference.get("appris_tag", ""),
                _parse_int(reference.get("protein_length_int", "")),
                _parse_int(reference.get("n_uniprot_links_int", "")),
            )

            pair_count = 0
            ref_accessions = accessions_by_transcript.get((species, reference["transcript_id"]), set())
            for alt in coding[1:]:
                alt_accessions = accessions_by_transcript.get((species, alt["transcript_id"]), set())
                shared_accessions = sorted(ref_accessions & alt_accessions)
                ref_plen = _parse_int(reference["protein_length_int"])
                alt_plen = _parse_int(alt["protein_length_int"])
                ref_tlen = _parse_int(reference["transcript_length_int"])
                alt_tlen = _parse_int(alt["transcript_length_int"])
                weak_label = ""
                weak_reason = ""
                if shared_accessions:
                    weak_label = "likely_preserved"
                    weak_reason = "shared_uniprot_accession"
                elif gene_info["clinvar_pathogenic_splice_count"] != "0":
                    weak_label = "splice_sensitive_gene"
                    weak_reason = "gene_has_clinvar_pathogenic_splice"
                elif reference.get("appris_tag", "") == "PRINCIPAL:1" and not alt.get("appris_tag", "").startswith("PRINCIPAL"):
                    weak_label = "reference_preferred"
                    weak_reason = "appris_principal_vs_nonprincipal"

                pair_writer.writerow(
                    {
                        "species": species,
                        "gene_id": gene_info["gencode_gene_id"],
                        "gene_name": gene_info["gene_name"],
                        "family_class": family_class,
                        "reference_source": reference_source,
                        "reference_transcript_id": reference["transcript_id"],
                        "alternative_transcript_id": alt["transcript_id"],
                        "reference_appris_tag": reference.get("appris_tag", ""),
                        "alternative_appris_tag": alt.get("appris_tag", ""),
                        "reference_protein_id": reference.get("protein_id", ""),
                        "alternative_protein_id": alt.get("protein_id", ""),
                        "reference_protein_length": ref_plen,
                        "alternative_protein_length": alt_plen,
                        "protein_length_delta": alt_plen - ref_plen,
                        "protein_length_ratio": f"{(alt_plen / ref_plen):.6f}" if ref_plen else "",
                        "reference_transcript_length": ref_tlen,
                        "alternative_transcript_length": alt_tlen,
                        "transcript_length_delta": alt_tlen - ref_tlen,
                        "reference_uniprot_links": reference.get("n_uniprot_links_int", "0"),
                        "alternative_uniprot_links": alt.get("n_uniprot_links_int", "0"),
                        "shared_uniprot_accessions": ";".join(shared_accessions),
                        "reference_is_membrane": reference.get("is_membrane_protein", "0"),
                        "alternative_is_membrane": alt.get("is_membrane_protein", "0"),
                        "reference_has_signal_peptide": reference.get("has_signal_peptide", "0"),
                        "alternative_has_signal_peptide": alt.get("has_signal_peptide", "0"),
                        "reference_is_kinase": reference.get("is_kinase", "0"),
                        "alternative_is_kinase": alt.get("is_kinase", "0"),
                        "reference_is_transcription_factor": reference.get("is_transcription_factor", "0"),
                        "alternative_is_transcription_factor": alt.get("is_transcription_factor", "0"),
                        "reference_interaction_partner_count": reference.get("interaction_partner_count", "0"),
                        "alternative_interaction_partner_count": alt.get("interaction_partner_count", "0"),
                        "reference_tm_feature_count": reference.get("tm_feature_count", "0"),
                        "alternative_tm_feature_count": alt.get("tm_feature_count", "0"),
                        "tm_feature_count_delta": _parse_int(alt.get("tm_feature_count", "0")) - _parse_int(reference.get("tm_feature_count", "0")),
                        "reference_signal_feature_count": reference.get("signal_feature_count", "0"),
                        "alternative_signal_feature_count": alt.get("signal_feature_count", "0"),
                        "signal_feature_count_delta": _parse_int(alt.get("signal_feature_count", "0")) - _parse_int(reference.get("signal_feature_count", "0")),
                        "reference_signal_max_end": reference.get("signal_max_end", "0"),
                        "alternative_signal_max_end": alt.get("signal_max_end", "0"),
                        "reference_extracellular_topology_aa": reference.get("extracellular_topology_aa", "0"),
                        "alternative_extracellular_topology_aa": alt.get("extracellular_topology_aa", "0"),
                        "extracellular_topology_aa_delta": _parse_int(alt.get("extracellular_topology_aa", "0")) - _parse_int(reference.get("extracellular_topology_aa", "0")),
                        "reference_cytoplasmic_topology_aa": reference.get("cytoplasmic_topology_aa", "0"),
                        "alternative_cytoplasmic_topology_aa": alt.get("cytoplasmic_topology_aa", "0"),
                        "cytoplasmic_topology_aa_delta": _parse_int(alt.get("cytoplasmic_topology_aa", "0")) - _parse_int(reference.get("cytoplasmic_topology_aa", "0")),
                        "reference_glycosylation_count": reference.get("glycosylation_count", "0"),
                        "alternative_glycosylation_count": alt.get("glycosylation_count", "0"),
                        "glycosylation_count_delta": _parse_int(alt.get("glycosylation_count", "0")) - _parse_int(reference.get("glycosylation_count", "0")),
                        "reference_disulfide_count": reference.get("disulfide_count", "0"),
                        "alternative_disulfide_count": alt.get("disulfide_count", "0"),
                        "disulfide_count_delta": _parse_int(alt.get("disulfide_count", "0")) - _parse_int(reference.get("disulfide_count", "0")),
                        "reference_phospho_feature_count": reference.get("phospho_feature_count", "0"),
                        "alternative_phospho_feature_count": alt.get("phospho_feature_count", "0"),
                        "phospho_feature_count_delta": _parse_int(alt.get("phospho_feature_count", "0")) - _parse_int(reference.get("phospho_feature_count", "0")),
                        "reference_ppi_binding_feature_count": reference.get("ppi_binding_feature_count", "0"),
                        "alternative_ppi_binding_feature_count": alt.get("ppi_binding_feature_count", "0"),
                        "ppi_binding_feature_count_delta": _parse_int(alt.get("ppi_binding_feature_count", "0")) - _parse_int(reference.get("ppi_binding_feature_count", "0")),
                        "reference_dna_binding_feature_count": reference.get("dna_binding_feature_count", "0"),
                        "alternative_dna_binding_feature_count": alt.get("dna_binding_feature_count", "0"),
                        "dna_binding_feature_count_delta": _parse_int(alt.get("dna_binding_feature_count", "0")) - _parse_int(reference.get("dna_binding_feature_count", "0")),
                        "reference_dna_region_feature_count": reference.get("dna_region_feature_count", "0"),
                        "alternative_dna_region_feature_count": alt.get("dna_region_feature_count", "0"),
                        "dna_region_feature_count_delta": _parse_int(alt.get("dna_region_feature_count", "0")) - _parse_int(reference.get("dna_region_feature_count", "0")),
                        "reference_zinc_finger_feature_count": reference.get("zinc_finger_feature_count", "0"),
                        "alternative_zinc_finger_feature_count": alt.get("zinc_finger_feature_count", "0"),
                        "zinc_finger_feature_count_delta": _parse_int(alt.get("zinc_finger_feature_count", "0")) - _parse_int(reference.get("zinc_finger_feature_count", "0")),
                        "reference_domain_count": reference.get("domain_count", "0"),
                        "alternative_domain_count": alt.get("domain_count", "0"),
                        "domain_count_delta": _parse_int(alt.get("domain_count", "0")) - _parse_int(reference.get("domain_count", "0")),
                        "reference_motif_count": reference.get("motif_count", "0"),
                        "alternative_motif_count": alt.get("motif_count", "0"),
                        "motif_count_delta": _parse_int(alt.get("motif_count", "0")) - _parse_int(reference.get("motif_count", "0")),
                        "reference_predicted_tm_count": reference.get("predicted_tm_count", "0"),
                        "alternative_predicted_tm_count": alt.get("predicted_tm_count", "0"),
                        "predicted_tm_count_delta": _parse_int(alt.get("predicted_tm_count", "0")) - _parse_int(reference.get("predicted_tm_count", "0")),
                        "reference_predicted_tm_total_span": reference.get("predicted_tm_total_span", "0"),
                        "alternative_predicted_tm_total_span": alt.get("predicted_tm_total_span", "0"),
                        "predicted_tm_total_span_delta": _parse_int(alt.get("predicted_tm_total_span", "0")) - _parse_int(reference.get("predicted_tm_total_span", "0")),
                        "reference_predicted_signal_candidate": reference.get("predicted_signal_candidate", "0"),
                        "alternative_predicted_signal_candidate": alt.get("predicted_signal_candidate", "0"),
                        "reference_predicted_signal_score": reference.get("predicted_signal_score", "0"),
                        "alternative_predicted_signal_score": alt.get("predicted_signal_score", "0"),
                        "predicted_signal_score_delta": f"{(float(alt.get('predicted_signal_score', 0) or 0) - float(reference.get('predicted_signal_score', 0) or 0)):.6f}",
                        "reference_glyco_motif_count": reference.get("glyco_motif_count", "0"),
                        "alternative_glyco_motif_count": alt.get("glyco_motif_count", "0"),
                        "glyco_motif_count_delta": _parse_int(alt.get("glyco_motif_count", "0")) - _parse_int(reference.get("glyco_motif_count", "0")),
                        "reference_cysteine_count": reference.get("cysteine_count", "0"),
                        "alternative_cysteine_count": alt.get("cysteine_count", "0"),
                        "cysteine_count_delta": _parse_int(alt.get("cysteine_count", "0")) - _parse_int(reference.get("cysteine_count", "0")),
                        "gene_has_appris_principal_1": gene_info.get("has_appris_principal_1", "0"),
                        "gene_clinvar_pathogenic_splice_count": gene_info.get("clinvar_pathogenic_splice_count", "0"),
                        "gene_hpa_has_extracellular_annotation": gene_info.get("hpa_has_extracellular_annotation", ""),
                        "gene_biogrid_partner_count": gene_info.get("biogrid_partner_count", "0"),
                        "gene_biogrid_physical_partner_count": gene_info.get("biogrid_physical_partner_count", "0"),
                        "weak_label": weak_label,
                        "weak_label_reason": weak_reason,
                    }
                )
                pair_rows += 1
                pair_count += 1

            gene_writer.writerow(
                {
                    "species": species,
                    "gene_id": gene_info["gencode_gene_id"],
                    "gene_name": gene_info["gene_name"],
                    "family_class": family_class,
                    "n_protein_coding_transcripts": gene_info["n_protein_coding_transcripts"],
                    "has_appris_principal_1": gene_info["has_appris_principal_1"],
                    "n_uniprot_entries": gene_info["n_uniprot_entries"],
                    "n_uniprot_isoforms": gene_info["n_uniprot_isoforms"],
                    "clinvar_pathogenic_splice_count": gene_info["clinvar_pathogenic_splice_count"],
                    "hpa_has_extracellular_annotation": gene_info["hpa_has_extracellular_annotation"],
                    "reference_transcript_id": reference["transcript_id"],
                    "reference_source": reference_source,
                    "n_pair_candidates": pair_count,
                }
            )
            gene_rows += 1

    print(f"[ok] Wrote {PAIR_OUT} rows={pair_rows}")
    print(f"[ok] Wrote {GENE_SET_OUT} rows={gene_rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
