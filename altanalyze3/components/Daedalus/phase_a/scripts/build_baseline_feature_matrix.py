#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


INTERIM = Path(__file__).resolve().parents[1] / "data" / "interim"
PROCESSED = Path(__file__).resolve().parents[1] / "data" / "processed"
IN_PATH = INTERIM / "priority_pair_subsets.tsv"
OUT_PATH = PROCESSED / "baseline_feature_matrix.tsv"


def _int(value: str) -> int:
    if not value:
        return 0
    try:
        return int(float(value))
    except ValueError:
        return 0


def _float(value: str) -> float:
    if not value:
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def _clip(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def _annotation_retention(
    ref_membrane: int,
    alt_membrane: int,
    ref_signal: int,
    alt_signal: int,
    ref_kinase: int,
    alt_kinase: int,
    ref_tf: int,
    alt_tf: int,
    alt_has_protein: int,
) -> float:
    components: list[float] = []
    if ref_membrane:
        components.append(float(alt_membrane))
    if ref_signal:
        components.append(float(alt_signal))
    if ref_kinase:
        components.append(float(alt_kinase))
    if ref_tf:
        components.append(float(alt_tf))
    if not components:
        components.append(float(alt_has_protein))
    return sum(components) / len(components)


def main() -> int:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    rows = 0
    with IN_PATH.open() as in_handle, OUT_PATH.open("w", encoding="utf-8", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = [
            "species",
            "gene_id",
            "gene_name",
            "split",
            "family_class",
            "reference_source",
            "label_preservation_seed",
            "task_global_seed",
            "task_membrane_seed",
            "task_surface_seed",
            "task_kinase_seed",
            "task_tf_seed",
            "task_high_confidence_seed",
            "feature_protein_ratio",
            "feature_protein_ratio_distance",
            "feature_protein_delta_norm",
            "feature_transcript_delta_norm",
            "feature_reference_uniprot_links",
            "feature_alternative_uniprot_links",
            "feature_alt_has_protein",
            "feature_ref_alt_same_protein_length",
            "feature_ref_is_membrane",
            "feature_alt_is_membrane",
            "feature_membrane_mismatch",
            "feature_ref_has_signal",
            "feature_alt_has_signal",
            "feature_signal_mismatch",
            "feature_ref_is_kinase",
            "feature_alt_is_kinase",
            "feature_kinase_mismatch",
            "feature_ref_is_tf",
            "feature_alt_is_tf",
            "feature_tf_mismatch",
            "feature_gene_has_appris_principal_1",
            "feature_gene_clinvar_splice",
            "feature_gene_hpa_extracellular",
            "feature_is_high_confidence",
            "feature_is_membrane_pair",
            "feature_is_surface_pair",
            "feature_is_kinase_pair",
            "feature_is_tf_pair",
            "feature_reference_source_mane",
            "feature_reference_source_appris",
            "feature_reference_source_uniprot_longest",
            "feature_reference_source_longest_protein",
            "feature_annotation_retention",
        ]
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            label = ""
            if row["is_uniprot_preserved_seed"] == "1":
                label = "1"
            elif row["is_appris_contrast_seed"] == "1":
                label = "0"

            ref_plen = _int(row["reference_protein_length"])
            alt_plen = _int(row["alternative_protein_length"])
            ref_tlen = _int(row["reference_transcript_length"])
            alt_tlen = _int(row["alternative_transcript_length"])
            protein_ratio = _clip(_float(row["protein_length_ratio"]), 0.0, 4.0)
            ref_uniprot = _int(row["reference_uniprot_links"])
            alt_uniprot = _int(row["alternative_uniprot_links"])
            ref_membrane = _int(row["reference_is_membrane"])
            alt_membrane = _int(row["alternative_is_membrane"])
            ref_signal = _int(row["reference_has_signal_peptide"])
            alt_signal = _int(row["alternative_has_signal_peptide"])
            ref_kinase = _int(row["reference_is_kinase"])
            alt_kinase = _int(row["alternative_is_kinase"])
            ref_tf = _int(row["reference_is_transcription_factor"])
            alt_tf = _int(row["alternative_is_transcription_factor"])
            alt_has_protein = int(bool(row["alternative_protein_id"]))
            writer.writerow(
                {
                    "species": row["species"],
                    "gene_id": row["gene_id"],
                    "gene_name": row["gene_name"],
                    "split": row["split"],
                    "family_class": row["family_class"],
                    "reference_source": row["reference_source"],
                    "label_preservation_seed": label,
                    "task_global_seed": int(bool(label)),
                    "task_membrane_seed": int(bool(label) and row["is_membrane_pair"] == "1"),
                    "task_surface_seed": int(bool(label) and row["is_surface_pair"] == "1"),
                    "task_kinase_seed": int(bool(label) and row["is_kinase_pair"] == "1"),
                    "task_tf_seed": int(bool(label) and row["is_tf_pair"] == "1"),
                    "task_high_confidence_seed": int(bool(label) and row["is_high_confidence"] == "1"),
                    "feature_protein_ratio": f"{protein_ratio:.6f}",
                    "feature_protein_ratio_distance": f"{abs(1.0 - protein_ratio):.6f}",
                    "feature_protein_delta_norm": f"{abs(alt_plen - ref_plen) / max(ref_plen, 1):.6f}",
                    "feature_transcript_delta_norm": f"{abs(alt_tlen - ref_tlen) / max(ref_tlen, 1):.6f}",
                    "feature_reference_uniprot_links": ref_uniprot,
                    "feature_alternative_uniprot_links": alt_uniprot,
                    "feature_alt_has_protein": alt_has_protein,
                    "feature_ref_alt_same_protein_length": int(ref_plen == alt_plen and ref_plen > 0),
                    "feature_ref_is_membrane": ref_membrane,
                    "feature_alt_is_membrane": alt_membrane,
                    "feature_membrane_mismatch": int(ref_membrane != alt_membrane),
                    "feature_ref_has_signal": ref_signal,
                    "feature_alt_has_signal": alt_signal,
                    "feature_signal_mismatch": int(ref_signal != alt_signal),
                    "feature_ref_is_kinase": ref_kinase,
                    "feature_alt_is_kinase": alt_kinase,
                    "feature_kinase_mismatch": int(ref_kinase != alt_kinase),
                    "feature_ref_is_tf": ref_tf,
                    "feature_alt_is_tf": alt_tf,
                    "feature_tf_mismatch": int(ref_tf != alt_tf),
                    "feature_gene_has_appris_principal_1": _int(row["gene_has_appris_principal_1"]),
                    "feature_gene_clinvar_splice": int(_int(row["gene_clinvar_pathogenic_splice_count"]) > 0),
                    "feature_gene_hpa_extracellular": _int(row["gene_hpa_has_extracellular_annotation"]),
                    "feature_is_high_confidence": _int(row["is_high_confidence"]),
                    "feature_is_membrane_pair": _int(row["is_membrane_pair"]),
                    "feature_is_surface_pair": _int(row["is_surface_pair"]),
                    "feature_is_kinase_pair": _int(row["is_kinase_pair"]),
                    "feature_is_tf_pair": _int(row["is_tf_pair"]),
                    "feature_reference_source_mane": int(row["reference_source"] == "MANE_SELECT"),
                    "feature_reference_source_appris": int(row["reference_source"].startswith("APPRIS")),
                    "feature_reference_source_uniprot_longest": int(row["reference_source"] == "UNIPROT_SUPPORTED_LONGEST"),
                    "feature_reference_source_longest_protein": int(row["reference_source"] == "LONGEST_PROTEIN"),
                    "feature_annotation_retention": f"{_annotation_retention(ref_membrane, alt_membrane, ref_signal, alt_signal, ref_kinase, alt_kinase, ref_tf, alt_tf, alt_has_protein):.6f}",
                }
            )
            rows += 1

    print(f"[ok] Wrote {OUT_PATH} rows={rows}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
