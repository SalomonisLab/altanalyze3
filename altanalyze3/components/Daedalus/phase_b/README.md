# Daedalus Phase B

Phase B turns the Phase A frozen benchmark substrate into a trainable
reference-conditioned model package.

## Goals

1. Convert wide pair records into model-ready reference/alternative feature blocks
2. Expand benchmark rows into multitask training instances
3. Train a shared torch model with task-aware conditioning
4. Save reusable checkpoints and evaluation reports

## Initial tasks

- `global`
- `membrane`
- `kinase`
- `transcription_factor`
- `surface`

`surface` is included for completeness, but Phase A already indicates it is the
weakest supervision regime.

## Current state

Phase B is now operational with two baseline tracks:

- a shared torch multitask model in `scripts/train_reference_delta_multitask.py`
- a leakage-resistant sklearn logistic baseline in
  `scripts/run_reference_delta_sklearn_baseline.py`

Primary artifacts:

- `data/processed/reference_delta_matrix.parquet`
- `data/processed/multitask_instances.parquet`
- `checkpoints/reference_delta_multitask_metrics.json`
- `checkpoints/reference_delta_multitask_report.md`
- `data/processed/reference_delta_sklearn_metrics.tsv`
- `checkpoints/reference_delta_sklearn_report.md`
- `scripts/predict_query_isoform.py`
- `QUERY_WORKFLOW.md`
- `config/task_registry.json`
- `schemas/supervised_benchmark_record.schema.json`
- `scripts/build_supervised_benchmark_table.py`
- `scripts/run_supervised_holdout_baselines.py`

## Practical interpretation

- `global`, `membrane`, `kinase`, and `transcription_factor` are viable
  Phase B targets on the current weak-label corpus
- `surface` is no longer the obvious weak point after adding topology, signal,
  PPI, PTM, and cysteine/glycosylation context
- the current modeled matrix keeps reference-known annotations but excludes
  direct alternative UniProt-region annotations to reduce weak-label leakage
- the TF branch now has a defined external validation plan based on hESC
  experimental outputs rather than weak labels alone
- query-time inference now supports automatic `MANE/APPRIS` reference selection
  and explicit channel evidence for PPI, PDI/DPI, localization, DNA binding,
  kinase signaling, signal peptide retention, and TM insertion/folding support
- supervised benchmark scaffolding now exists for promoting those channel
  outputs into explicit trainable and benchmarked task families
- supervised benchmark tables now preserve the frozen Phase B holdout split and
  can be evaluated with explicit train/val/test baselines
- promoted-task proxy holdout benchmarking is now populated and reported in:
  - `data/processed/supervised_holdout_metrics.tsv`
  - `data/processed/supervised_holdout_report.md`

## Promoted-task proxy benchmark snapshot

Current proxy holdout test AUROCs:

- `signal_retention`: `0.805449`
- `tm_insertion`: `0.838737`
- `tm_fold`: `0.823234`
- `ppi_impact`: `0.636388`
- `pdi_impact`: `0.521164`
- `tf_dna_binding`: `0.521164`
- `tf_transcriptional_activity`: `0.521164`
- `tf_perturbation_response`: `0.521164`
- `localization_shift`: `0.525888`

Key caveats:

- these are proxy weak-label supervised tasks, not external assay validation
- TF proxy tasks currently share the same underlying proxy source, so the
  reported metrics are not independent
- `kinase_signaling` is currently excluded from reported holdout performance
  because the proxy labels collapsed to a single class

See `STATUS.md` for the frozen Phase B summary and `VALIDATION.md` for the
external validation plan. See `QUERY_WORKFLOW.md` for the query-time inference
path.
