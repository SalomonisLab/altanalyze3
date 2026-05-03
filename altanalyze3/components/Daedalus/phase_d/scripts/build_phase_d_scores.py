#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Daedalus.phase_d.objectives import score_phase_d_objectives


PHASE_C = Path(__file__).resolve().parents[2] / "phase_c" / "data" / "processed"
PHASE_D = Path(__file__).resolve().parents[1] / "data" / "processed"

IN_PATH = PHASE_C / "biochem_delta_matrix.parquet"
OUT_SCORES_TSV = PHASE_D / "phase_d_objective_scores.tsv"
OUT_SCORES_PARQUET = PHASE_D / "phase_d_objective_scores.parquet"
OUT_INSTANCES_TSV = PHASE_D / "phase_d_task_instances.tsv"
OUT_INSTANCES_PARQUET = PHASE_D / "phase_d_task_instances.parquet"
OUT_SCHEMA = PHASE_D / "phase_d_schema.json"
OUT_TASKS = PHASE_D / "phase_d_task_registry.tsv"

TM_TASKS = [
    ("tm_secretory_pathway_compatibility_score", "secretory_pathway_compatibility", "membrane"),
    ("tm_membrane_insertion_topology_score", "membrane_insertion_topology", "membrane"),
    ("tm_extracellular_altered_segment_accessibility_score", "extracellular_altered_segment_accessibility", "membrane"),
    ("tm_folding_stability_qc_escape_score", "folding_stability_qc_escape", "membrane"),
    ("tm_cell_surface_localization_score", "cell_surface_localization", "membrane"),
    ("tm_antibody_targetability_score", "antibody_targetability", "membrane"),
    ("tm_stable_functional_score", "tm_stable_functional", "membrane"),
    ("tm_non_stable_non_functional_score", "tm_non_stable_non_functional", "membrane"),
]
TF_TASKS = [
    ("tf_nuclear_localization_score", "nuclear_localization", "transcription_factor"),
    ("tf_dbd_integrity_specificity_score", "dbd_integrity_specificity_shift", "transcription_factor"),
    ("tf_cofactor_ppi_rewiring_risk_score", "cofactor_ppi_rewiring", "transcription_factor"),
    ("tf_activation_repression_competence_score", "activation_repression_competence", "transcription_factor"),
    ("tf_dominant_negative_neomorphic_risk_score", "dominant_negative_neomorphic_behavior", "transcription_factor"),
    ("tf_overall_regulatory_functionality_score", "overall_regulatory_functionality", "transcription_factor"),
    ("tf_non_stable_non_functional_score", "tf_non_stable_non_functional", "transcription_factor"),
]
APPLICABILITY_COLUMN_BY_FAMILY = {
    "membrane": "task_membrane_seed",
    "transcription_factor": "task_tf_seed",
}


def main() -> int:
    PHASE_D.mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet(IN_PATH)
    scored = df.copy()
    scores_df = pd.DataFrame([score_phase_d_objectives(row) for row in scored.to_dict(orient="records")])
    scored = pd.concat([scored, scores_df], axis=1)
    scored.to_csv(OUT_SCORES_TSV, sep="\t", index=False)
    scored.to_parquet(OUT_SCORES_PARQUET, index=False)

    instance_rows = []
    id_cols = ["species", "gene_id", "gene_name", "reference_transcript_id", "alternative_transcript_id", "family_class", "reference_source", "split", "change_class"]
    for col, task_name, family in TM_TASKS + TF_TASKS:
        sub = scored[id_cols + [col]].copy()
        applicability_col = APPLICABILITY_COLUMN_BY_FAMILY.get(family)
        if applicability_col and applicability_col in scored.columns:
            mask = pd.to_numeric(scored[applicability_col], errors="coerce").fillna(0).astype(int) == 1
            sub = sub[mask].copy()
        sub = sub[sub[col].notna()].rename(columns={col: "objective_score"})
        sub["task"] = task_name
        sub["task_family"] = family
        instance_rows.append(sub)
    instances = pd.concat(instance_rows, axis=0, ignore_index=True)
    instances.to_csv(OUT_INSTANCES_TSV, sep="\t", index=False)
    instances.to_parquet(OUT_INSTANCES_PARQUET, index=False)

    task_registry = pd.DataFrame(
        [
            {
                "task": "secretory_pathway_compatibility",
                "task_family": "membrane",
                "score_column": "tm_secretory_pathway_compatibility_score",
                "direction": "higher_is_better",
                "description": "Compatibility with signal-driven ER entry and early secretory pathway routing.",
            },
            {
                "task": "membrane_insertion_topology",
                "task_family": "membrane",
                "score_column": "tm_membrane_insertion_topology_score",
                "direction": "higher_is_better",
                "description": "Compatibility with stable TM insertion and gross topology preservation.",
            },
            {
                "task": "extracellular_altered_segment_accessibility",
                "task_family": "membrane",
                "score_column": "tm_extracellular_altered_segment_accessibility_score",
                "direction": "higher_is_better",
                "description": "Likelihood that the altered segment is extracellularly exposed and accessible.",
            },
            {
                "task": "folding_stability_qc_escape",
                "task_family": "membrane",
                "score_column": "tm_folding_stability_qc_escape_score",
                "direction": "higher_is_better",
                "description": "Compatibility with extracellular folding, maturation, and escape from quality-control failure.",
            },
            {
                "task": "cell_surface_localization",
                "task_family": "membrane",
                "score_column": "tm_cell_surface_localization_score",
                "direction": "higher_is_better",
                "description": "Likelihood that the isoform reaches and persists at the cell surface.",
            },
            {
                "task": "antibody_targetability",
                "task_family": "membrane",
                "score_column": "tm_antibody_targetability_score",
                "direction": "higher_is_better",
                "description": "Composite accessibility and surface-persistence score for antibody targeting.",
            },
            {
                "task": "tm_stable_functional",
                "task_family": "membrane",
                "score_column": "tm_stable_functional_score",
                "direction": "higher_is_better",
                "description": "Aggregate membrane-isoform stability and functionality score.",
            },
            {
                "task": "tm_non_stable_non_functional",
                "task_family": "membrane",
                "score_column": "tm_non_stable_non_functional_score",
                "direction": "lower_is_better",
                "description": "Aggregate risk that the membrane isoform is non-stable or non-functional.",
            },
            {
                "task": "nuclear_localization",
                "task_family": "transcription_factor",
                "score_column": "tf_nuclear_localization_score",
                "direction": "higher_is_better",
                "description": "Compatibility with nuclear localization and retention.",
            },
            {
                "task": "dbd_integrity_specificity_shift",
                "task_family": "transcription_factor",
                "score_column": "tf_dbd_integrity_specificity_score",
                "direction": "higher_is_better",
                "description": "Compatibility with intact DNA-binding domain grammar and low specificity-disrupting change burden.",
            },
            {
                "task": "cofactor_ppi_rewiring",
                "task_family": "transcription_factor",
                "score_column": "tf_cofactor_ppi_rewiring_risk_score",
                "direction": "lower_is_better",
                "description": "Risk that cofactor or PPI wiring is disrupted or rewired by the isoform change.",
            },
            {
                "task": "activation_repression_competence",
                "task_family": "transcription_factor",
                "score_column": "tf_activation_repression_competence_score",
                "direction": "higher_is_better",
                "description": "Compatibility with retained activation and repression module function.",
            },
            {
                "task": "dominant_negative_neomorphic_behavior",
                "task_family": "transcription_factor",
                "score_column": "tf_dominant_negative_neomorphic_risk_score",
                "direction": "lower_is_better",
                "description": "Risk of dominant-negative or neomorphic behavior given retained localization/interface capacity but altered regulation.",
            },
            {
                "task": "overall_regulatory_functionality",
                "task_family": "transcription_factor",
                "score_column": "tf_overall_regulatory_functionality_score",
                "direction": "higher_is_better",
                "description": "Composite TF functionality score integrating localization, DBD integrity, cofactor rewiring, and activation competence.",
            },
            {
                "task": "tf_non_stable_non_functional",
                "task_family": "transcription_factor",
                "score_column": "tf_non_stable_non_functional_score",
                "direction": "lower_is_better",
                "description": "Aggregate risk that the TF isoform is non-functional or unstable as a regulator.",
            },
        ]
    )
    task_registry.to_csv(OUT_TASKS, sep="\t", index=False)

    schema = {
        "source_matrix": str(IN_PATH),
        "score_outputs": [str(OUT_SCORES_TSV), str(OUT_SCORES_PARQUET)],
        "instance_outputs": [str(OUT_INSTANCES_TSV), str(OUT_INSTANCES_PARQUET)],
        "task_registry": str(OUT_TASKS),
        "tm_tasks": [task for _, task, _ in TM_TASKS],
        "tf_tasks": [task for _, task, _ in TF_TASKS],
        "description": "Phase D objective-function task tables derived from the Phase C biochemical matrix.",
    }
    OUT_SCHEMA.write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")

    print(f"[ok] Wrote {OUT_SCORES_TSV}")
    print(f"[ok] Wrote {OUT_SCORES_PARQUET}")
    print(f"[ok] Wrote {OUT_INSTANCES_TSV}")
    print(f"[ok] Wrote {OUT_INSTANCES_PARQUET}")
    print(f"[ok] Wrote {OUT_TASKS}")
    print(f"[ok] Wrote {OUT_SCHEMA}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
