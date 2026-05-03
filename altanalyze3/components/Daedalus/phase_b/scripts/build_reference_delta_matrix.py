#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PHASE_A = Path(__file__).resolve().parents[2] / "phase_a" / "data"
PHASE_B = Path(__file__).resolve().parents[1] / "data" / "processed"
IN_PATH = PHASE_A / "interim" / "priority_pair_subsets.tsv"
OUT_TSV = PHASE_B / "reference_delta_matrix.tsv"
OUT_PARQUET = PHASE_B / "reference_delta_matrix.parquet"
OUT_SCHEMA = PHASE_B / "reference_delta_schema.json"


def _flag(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().any():
        return numeric.fillna(0).astype(int)
    return series.fillna("").astype(str).str.strip().isin({"1", "true", "True"}).astype(int)


def main() -> int:
    PHASE_B.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(IN_PATH, sep="\t", low_memory=False)
    numeric = pd.DataFrame(
        {
            "reference_protein_length": df["reference_protein_length"].fillna(0).astype(float),
            "alternative_protein_length": df["alternative_protein_length"].fillna(0).astype(float),
            "reference_transcript_length": df["reference_transcript_length"].fillna(0).astype(float),
            "alternative_transcript_length": df["alternative_transcript_length"].fillna(0).astype(float),
            "reference_uniprot_links": df["reference_uniprot_links"].fillna(0).astype(float),
            "reference_is_membrane": df["reference_is_membrane"].fillna(0).astype(float),
            "reference_has_signal_peptide": df["reference_has_signal_peptide"].fillna(0).astype(float),
            "reference_is_kinase": df["reference_is_kinase"].fillna(0).astype(float),
            "reference_is_tf": df["reference_is_transcription_factor"].fillna(0).astype(float),
            "gene_has_appris_principal_1": df["gene_has_appris_principal_1"].fillna(0).astype(float),
            "gene_clinvar_pathogenic_splice_count": df["gene_clinvar_pathogenic_splice_count"].fillna(0).astype(float),
            "gene_hpa_has_extracellular_annotation": df["gene_hpa_has_extracellular_annotation"].fillna(0).astype(float),
            "gene_biogrid_partner_count": df["gene_biogrid_partner_count"].fillna(0).astype(float),
            "gene_biogrid_physical_partner_count": df["gene_biogrid_physical_partner_count"].fillna(0).astype(float),
            "protein_length_delta": df["protein_length_delta"].fillna(0).astype(float),
            "protein_length_ratio": df["protein_length_ratio"].fillna(0).astype(float),
            "transcript_length_delta": df["transcript_length_delta"].fillna(0).astype(float),
            "ref_alt_same_protein_length": (df["reference_protein_length"].fillna(0).astype(float) == df["alternative_protein_length"].fillna(0).astype(float)).astype(float),
            "alt_has_protein": (df["alternative_protein_id"].fillna("") != "").astype(float),
            "reference_tm_feature_count": df["reference_tm_feature_count"].fillna(0).astype(float),
            "reference_signal_feature_count": df["reference_signal_feature_count"].fillna(0).astype(float),
            "reference_extracellular_topology_aa": df["reference_extracellular_topology_aa"].fillna(0).astype(float),
            "reference_cytoplasmic_topology_aa": df["reference_cytoplasmic_topology_aa"].fillna(0).astype(float),
            "reference_glycosylation_count": df["reference_glycosylation_count"].fillna(0).astype(float),
            "reference_disulfide_count": df["reference_disulfide_count"].fillna(0).astype(float),
            "reference_phospho_feature_count": df["reference_phospho_feature_count"].fillna(0).astype(float),
            "reference_ppi_binding_feature_count": df["reference_ppi_binding_feature_count"].fillna(0).astype(float),
            "reference_dna_binding_feature_count": df["reference_dna_binding_feature_count"].fillna(0).astype(float),
            "reference_dna_region_feature_count": df["reference_dna_region_feature_count"].fillna(0).astype(float),
            "reference_zinc_finger_feature_count": df["reference_zinc_finger_feature_count"].fillna(0).astype(float),
            "reference_domain_count": df["reference_domain_count"].fillna(0).astype(float),
            "reference_motif_count": df["reference_motif_count"].fillna(0).astype(float),
            "alternative_predicted_tm_count": df["alternative_predicted_tm_count"].fillna(0).astype(float),
            "reference_predicted_tm_count": df["reference_predicted_tm_count"].fillna(0).astype(float),
            "predicted_tm_count_delta": df["predicted_tm_count_delta"].fillna(0).astype(float),
            "reference_predicted_tm_total_span": df["reference_predicted_tm_total_span"].fillna(0).astype(float),
            "alternative_predicted_tm_total_span": df["alternative_predicted_tm_total_span"].fillna(0).astype(float),
            "predicted_tm_total_span_delta": df["predicted_tm_total_span_delta"].fillna(0).astype(float),
            "reference_predicted_signal_candidate": df["reference_predicted_signal_candidate"].fillna(0).astype(float),
            "alternative_predicted_signal_candidate": df["alternative_predicted_signal_candidate"].fillna(0).astype(float),
            "reference_predicted_signal_score": df["reference_predicted_signal_score"].fillna(0).astype(float),
            "alternative_predicted_signal_score": df["alternative_predicted_signal_score"].fillna(0).astype(float),
            "predicted_signal_score_delta": df["predicted_signal_score_delta"].fillna(0).astype(float),
            "alternative_glyco_motif_count": df["alternative_glyco_motif_count"].fillna(0).astype(float),
            "alternative_cysteine_count": df["alternative_cysteine_count"].fillna(0).astype(float),
        }
    )
    meta = pd.DataFrame(
        {
            "species": df["species"].astype(str),
            "gene_id": df["gene_id"].astype(str),
            "gene_name": df["gene_name"].astype(str),
            "split": df["split"].astype(str),
            "family_class": df["family_class"].astype(str),
            "reference_source": df["reference_source"].astype(str),
            "reference_transcript_id": df["reference_transcript_id"].astype(str),
            "alternative_transcript_id": df["alternative_transcript_id"].astype(str),
            "task_global_seed": ((_flag(df["is_uniprot_preserved_seed"]) == 1) | (_flag(df["is_appris_contrast_seed"]) == 1)).astype(int),
            "task_membrane_seed": _flag(df["is_membrane_pair"]),
            "task_surface_seed": _flag(df["is_surface_pair"]),
            "task_kinase_seed": _flag(df["is_kinase_pair"]),
            "task_tf_seed": _flag(df["is_tf_pair"]),
            "label_preservation_seed": pd.Series(
                [
                    1 if x == 1 else 0 if y == 1 else pd.NA
                    for x, y in zip(_flag(df["is_uniprot_preserved_seed"]), _flag(df["is_appris_contrast_seed"]))
                ],
                dtype="Int64",
            ),
        }
    )
    out = pd.concat([meta, numeric], axis=1)
    out.to_csv(OUT_TSV, sep="\t", index=False)
    out.to_parquet(OUT_PARQUET, index=False)

    schema = {
        "categorical_columns": [
            "species",
            "gene_id",
            "gene_name",
            "split",
            "family_class",
            "reference_source",
            "reference_transcript_id",
            "alternative_transcript_id",
        ],
        "task_columns": [
            "task_global_seed",
            "task_membrane_seed",
            "task_surface_seed",
            "task_kinase_seed",
            "task_tf_seed",
        ],
        "label_column": "label_preservation_seed",
        "numeric_columns": list(numeric.columns),
    }
    with OUT_SCHEMA.open("w", encoding="utf-8") as handle:
        json.dump(schema, handle, indent=2)

    print(f"[ok] Wrote {OUT_TSV} rows={len(out)}")
    print(f"[ok] Wrote {OUT_PARQUET}")
    print(f"[ok] Wrote {OUT_SCHEMA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
