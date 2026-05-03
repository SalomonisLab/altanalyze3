# Daedalus Phase B Validation

## External validation strategy

Weak-label benchmarks are only the first gate. The transcription-factor branch
will be validated against orthogonal experimental data from the hESC analyses
described in the linked paper and already precomputed internally.

## Transcription-factor validation axes

1. TF isoform PPI
2. TF isoform PDI
3. TF isoform transcriptional activation capacity
4. TF isoform cellular localization
5. TF isoform perturbation response

## Expanded supervised task set

The workflow is now being promoted from weak-label preservation-only tasks to a
broader supervised benchmark program. The intended supervised task families are:

- `ppi_impact`
- `pdi_impact`
- `localization_shift`
- `tf_dna_binding`
- `tf_transcriptional_activity`
- `tf_perturbation_response`
- `kinase_signaling`
- `signal_retention`
- `tm_insertion`
- `tm_fold`

Task definitions are tracked in:

- `Daedalus/phase_b/config/task_registry.json`

These are the tasks needed for the channel outputs to become real model claims
rather than structured evidence overlays.

## Intended mapping from model outputs to validation tasks

- `transcription_factor` preservation score:
  global reference-relative TF preservation probability
- DNA-binding features:
  evaluate against PDI and DNA-binding competent vs impaired isoforms
- PPI-binding / phospho-interaction features:
  evaluate against observed TF isoform PPIs
- localization-sensitive features:
  evaluate against nuclear vs non-nuclear localization shifts
- perturbation-linked functional preservation:
  evaluate against the magnitude and reproducibility of bootstrap-supported DEGs

## Expected validation tables

Each row should represent a `reference_isoform -> alternative_isoform`
comparison or a single alternative isoform relative to a canonical reference.

Recommended minimal fields:

- `gene_name`
- `reference_transcript_id`
- `alternative_transcript_id`
- `species`
- `cell_context`
- `ppi_preserved`
- `ppi_delta_score`
- `pdi_preserved`
- `pdi_delta_score`
- `transcriptional_activity_score`
- `localization_class`
- `localization_preserved`
- `perturbation_deg_count`
- `perturbation_effect_size`
- `bootstrap_significant`

For the broader supervised benchmark workflow, the standardized schema is:

- `Daedalus/phase_b/schemas/supervised_benchmark_record.schema.json`

## Primary external metrics

- AUROC / AUPRC for binary preservation tasks
- rank correlation for continuous activity and perturbation outputs
- calibration for preservation probabilities
- within-gene ranking:
  canonical or experimentally stronger isoforms should rank above disrupted or
  non-functional alternatives

## Benchmark assembly

External assay tables should be placed in:

- `Daedalus/phase_b/data/external/`

and normalized with:

```bash
Daedalus/.venv/bin/python Daedalus/phase_b/scripts/build_supervised_benchmark_table.py
```

That script merges assay records onto the frozen Phase B reference-delta pairs
using:

- `species`
- `gene_name`
- `reference_transcript_id`
- `alternative_transcript_id`

Outputs:

- `data/processed/supervised_benchmark_instances.tsv`
- `data/processed/supervised_benchmark_instances.parquet`
- `data/processed/supervised_benchmark_summary.tsv`
- `data/processed/supervised_benchmark_schema.json`

These benchmark instances inherit the frozen Phase B `train` / `val` / `test`
split from the underlying reference-delta matrix. That is the holdout
mechanism. External assay rows are not evaluated by random split.

To run holdout baselines:

```bash
Daedalus/.venv/bin/python Daedalus/phase_b/scripts/run_supervised_holdout_baselines.py
```

Outputs:

- `data/processed/supervised_holdout_metrics.tsv`
- `data/processed/supervised_holdout_report.md`

## Current promoted-task holdout benchmark

The supervised benchmark path has now been run end-to-end using a large
proxy-labeled benchmark assembled from the frozen Phase A corpus:

- source table:
  - `Daedalus/phase_b/data/external/proxy_supervised_tasks.tsv`
- benchmark instances:
  - `Daedalus/phase_b/data/processed/supervised_benchmark_instances.parquet`
- holdout metrics:
  - `Daedalus/phase_b/data/processed/supervised_holdout_metrics.tsv`
- holdout report:
  - `Daedalus/phase_b/data/processed/supervised_holdout_report.md`

These are still weak-label proxy tasks, not external assay validation. They are
useful for feasibility and for checking whether the promoted task heads are
learnable on the existing reference-delta matrix.

Current holdout results:

- `localization_shift`
  - `val` n=10,840 AP=0.437278 AUROC=0.834509
  - `test` n=10,837 AP=0.159922 AUROC=0.525888
- `ppi_impact`
  - `val` n=6,274 AP=0.461970 AUROC=0.818197
  - `test` n=6,333 AP=0.168298 AUROC=0.636388
- `pdi_impact`
  - `val` n=2,583 AP=0.256386 AUROC=0.708129
  - `test` n=2,685 AP=0.161066 AUROC=0.521164
- `tf_dna_binding`
  - `val` n=2,583 AP=0.256386 AUROC=0.708129
  - `test` n=2,685 AP=0.161066 AUROC=0.521164
- `tf_transcriptional_activity`
  - `val` n=2,583 AP=0.256386 AUROC=0.708129
  - `test` n=2,685 AP=0.161066 AUROC=0.521164
- `tf_perturbation_response`
  - `val` n=2,583 AP=0.256386 AUROC=0.708129
  - `test` n=2,685 AP=0.161066 AUROC=0.521164
- `signal_retention`
  - `val` n=3,390 AP=0.480827 AUROC=0.824174
  - `test` n=3,154 AP=0.459793 AUROC=0.805449
- `tm_insertion`
  - `val` n=5,241 AP=0.505154 AUROC=0.818212
  - `test` n=5,697 AP=0.433830 AUROC=0.838737
- `tm_fold`
  - `val` n=5,241 AP=0.313015 AUROC=0.819110
  - `test` n=5,697 AP=0.357609 AUROC=0.823234

Current interpretation:

- `signal_retention`, `tm_insertion`, and `tm_fold` are the strongest promoted
  membrane-oriented proxy tasks so far.
- `ppi_impact` is learnable under proxy supervision, but there is a strong
  `val` to `test` drop that needs investigation before claiming robustness.
- the TF-related proxy tasks are only modestly learnable and currently collapse
  to identical metrics because they are all derived from the same underlying
  proxy source.
- `localization_shift` is unstable across splits and should not be treated as a
  serious result until external localization data are added.
- `kinase_signaling` produced proxy benchmark rows but no usable holdout metric,
  because the current proxy labels collapsed to a single class. That proxy
  definition needs to be rebuilt before the task is scientifically meaningful.

## What this solves

This moves the TF branch beyond weak annotation-derived supervision. It provides
real external evidence for:

- preserved or altered TF protein interactions
- preserved or altered promoter/DNA engagement
- preserved or altered transcriptional activation
- preserved or altered subcellular localization
- preserved or altered downstream regulatory impact

That is the correct validation layer for a serious methods claim.

The same framework now extends that logic to:

- kinase signaling competence
- membrane signal retention
- transmembrane insertion support
- transmembrane fold support
- broader PPI-driven functional disruption
