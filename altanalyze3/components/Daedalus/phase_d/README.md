# Daedalus Phase D

Phase D refactors the broad membrane/TF readout into explicit objective
functions aligned to the intended biology.

## Current objective families

### Transmembrane / surface

- `secretory_pathway_compatibility`
- `membrane_insertion_topology`
- `extracellular_altered_segment_accessibility`
- `folding_stability_qc_escape`
- `cell_surface_localization`
- `antibody_targetability`

### Transcription factor

- `nuclear_localization`
- `dbd_integrity_specificity_shift`
- `cofactor_ppi_rewiring`
- `activation_repression_competence`
- `dominant_negative_neomorphic_behavior`
- `overall_regulatory_functionality`

## What exists now

Phase D currently implements:

- explicit objective functions in:
  - `Daedalus/phase_d/objectives.py`
- a task-table builder that scores the full Phase C biochemical matrix:
  - `Daedalus/phase_d/scripts/build_phase_d_scores.py`
- a benchmark-instance builder:
  - `Daedalus/phase_d/scripts/build_phase_d_benchmark_instances.py`
- a family-aware multitask trainer:
  - `Daedalus/phase_d/scripts/train_phase_d_multitask.py`
- a model comparison harness:
  - `Daedalus/phase_d/scripts/run_phase_d_model_comparison.py`

This is an objective-function layer built on the Phase C feature substrate. It
now includes a bootstrap proxy supervised benchmark path using objective-derived
labels. It is not yet a final external-assay-trained Phase D model family.

Phase D now also includes:

- aggregate family-level annotations:
  - `stable_functional`
  - `non_stable_non_functional`
  - `indeterminate`
- an explicit training and benchmark specification:
  - `Daedalus/phase_d/TRAINING_AND_BENCHMARK_SPEC.md`
- per-task label definitions:
  - `Daedalus/phase_d/data/processed/phase_d_label_definitions.tsv`
- curated TM-negative transcript-pair labels for surface-localization benchmarking:
  - `Daedalus/phase_a/data/interim/tm_negatives/tm_negative_pair_labels.tsv`

## Build outputs

```bash
Daedalus/.venv/bin/python Daedalus/phase_d/scripts/build_phase_d_scores.py
Daedalus/.venv/bin/python Daedalus/phase_d/scripts/build_phase_d_benchmark_instances.py
Daedalus/.venv/bin/python Daedalus/phase_d/scripts/train_phase_d_multitask.py
Daedalus/.venv/bin/python Daedalus/phase_d/scripts/run_phase_d_model_comparison.py
```

Outputs:

- `Daedalus/phase_d/data/processed/phase_d_objective_scores.tsv`
- `Daedalus/phase_d/data/processed/phase_d_objective_scores.parquet`
- `Daedalus/phase_d/data/processed/phase_d_task_instances.tsv`
- `Daedalus/phase_d/data/processed/phase_d_task_instances.parquet`
- `Daedalus/phase_d/data/processed/phase_d_task_registry.tsv`
- `Daedalus/phase_d/data/processed/phase_d_benchmark_instances.tsv`
- `Daedalus/phase_d/data/processed/phase_d_benchmark_instances.parquet`
- `Daedalus/phase_d/data/processed/phase_d_benchmark_schema.json`
- `Daedalus/phase_d/data/processed/phase_d_label_definitions.tsv`
- `Daedalus/phase_d/TRAINING_AND_BENCHMARK_SPEC.md`
- `Daedalus/phase_d/checkpoints/phase_d_multitask_report.md`
- `Daedalus/phase_d/checkpoints/phase_d_multitask_metrics.json`
- `Daedalus/phase_d/checkpoints/phase_d_model_comparison_report.md`

## Current proxy benchmark status

The engineering benchmark is now complete for the Phase D proxy labels:

- `phase_d_family_multitask`
- `phase_d_logistic`
- `phase_d_hist_gradient_boosting`
- `phase_d_xgboost`

Artifacts:

- `Daedalus/phase_d/checkpoints/phase_d_multitask_report.md`
- `Daedalus/phase_d/checkpoints/phase_d_multitask_metrics.json`
- `Daedalus/phase_d/checkpoints/phase_d_model_comparison_metrics.tsv`
- `Daedalus/phase_d/checkpoints/phase_d_model_comparison_report.md`
- `Daedalus/phase_d/data/processed/phase_d_label_balance.tsv`
- `Daedalus/phase_d/checkpoints/cell_surface_localization_curated_negative_metrics.json`

Important interpretation boundary:

- many Phase D proxy tasks are extremely label-imbalanced after thresholding
- several tasks have zero negatives or zero positives in the test split
- as a result, many proxy metrics are trivially perfect and are not scientifically
  meaningful model evidence

Tasks with especially degenerate proxy test labels include:

- `secretory_pathway_compatibility`
- `nuclear_localization`
- `overall_regulatory_functionality`
- `activation_repression_competence`
- `tf_non_stable_non_functional`
- `dominant_negative_neomorphic_behavior`

One Phase D task now has explicit held-out curated negatives:

- `cell_surface_localization`

These come from UniProt isoforms with at least one TM domain that lose
`Cell membrane` localization and instead localize elsewhere, joined onto the
Daedalus transcript-pair substrate through the UniProt-GENCODE mapping and pair
tables.

Current held-out `cell_surface_localization` test composition:

- total test rows: `533`
- positives: `520`
- negatives: `13`
- curated held-out negatives: `7`

Current task-level metrics after injecting those curated negatives:

- `phase_d_logistic`
  - AP `0.999370`
  - AUROC `0.975148`
  - balanced accuracy `0.790385`
- `phase_d_hist_gradient_boosting`
  - AP `0.999830`
  - AUROC `0.993195`
  - balanced accuracy `0.692308`
- `phase_d_xgboost`
  - AP `0.999802`
  - AUROC `0.992160`
  - balanced accuracy `0.692308`
- `phase_d_family_multitask`
  - AP `0.999322`
  - AUROC `0.973669`
  - balanced accuracy `0.729808`

So the current Phase D benchmark should be used only as:

- pipeline validation
- model-family plumbing validation
- a check that the benchmark harness runs end-to-end

It should not yet be used as evidence that Phase D objectives are solved. Real
value now requires external assay labels for the new TM and TF heads.

## Interpretation boundary

The current Phase D scores are objective-function scores derived from the Phase
C biochemical feature matrix. They are designed to replace the earlier broad
readout with explicit TM and TF task definitions.

The current supervised path is a bootstrap proxy benchmark:

- labels are thresholded from the objective scores
- ambiguous rows are excluded
- holdouts remain gene-disjoint

This is useful to validate the engineering stack and compare model families, but
it is not a replacement for external assay supervision.
