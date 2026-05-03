# Daedalus Phase D Training And Benchmark Spec

Phase D should be trained as a family-aware multitask problem, not as a single
binary classifier.

## Core framing

Each isoform is evaluated relative to an automatically selected reference
transcript. The model should predict explicit functional sub-capacities, then
derive a family-specific aggregate annotation:

- `stable_functional`
- `non_stable_non_functional`
- `indeterminate`

This aggregate is downstream of the task heads. It should not replace the task
heads.

## Model family recommendation

Use the same benchmark ladder that worked best in prior phases where it is
applicable:

1. strong tabular baseline
   - `hist_gradient_boosting`
   - `xgboost`
2. deep comparator
   - class-specific autoencoder
   - DeepImmuno-style CNN and CNN+autoencoder
3. family-specific multitask model
   - shared encoder
   - membrane expert
   - TF expert
   - task-specific heads

Do not use logistic regression as the primary Phase D model. Keep it only as a
minimal sanity-check baseline.

## Phase D task heads

### Membrane / surface

Primary heads:

1. `secretory_pathway_compatibility`
2. `membrane_insertion_topology`
3. `folding_stability_qc_escape`
4. `cell_surface_localization`

Secondary heads:

5. `extracellular_altered_segment_accessibility`
6. `antibody_targetability`

Aggregate:

7. `tm_stable_functional`
8. `tm_non_stable_non_functional`

### Transcription factor

Primary heads:

1. `nuclear_localization`
2. `dbd_integrity_specificity_shift`
3. `activation_repression_competence`

Secondary heads:

4. `cofactor_ppi_rewiring`
5. `dominant_negative_neomorphic_behavior`

Aggregate:

6. `overall_regulatory_functionality`
7. `tf_non_stable_non_functional`

## Recommended loss

Use multitask supervised loss with optional latent regularization:

`L = sum(lambda_t * L_t) + lambda_rec * L_reconstruction + lambda_rank * L_rank`

Recommended head types:

- binary / ordinal heads for compatibility tasks
- bounded regression or ordinal classification for graded tasks
- aggregate heads can be trained either:
  - directly, if labeled
  - or derived post hoc from the primary heads

## Labeling strategy

Use three label states per task:

- `positive`
- `negative`
- `ambiguous_exclude`

Only `positive` and `negative` enter supervised training.

### Membrane task labels

| Task | Positive label | Negative label | Ambiguous / exclude |
| --- | --- | --- | --- |
| `secretory_pathway_compatibility` | Has curated or strong predicted signal-entry grammar and compatible N-terminus | Required entry signal absent or clearly disrupted with no plausible alternative entry | TM proteins not requiring classical secretory entry; unclear entry mechanism |
| `membrane_insertion_topology` | TM architecture/topology consistent with stable membrane insertion | TM anchor lost, gross TM count/order disruption, or incompatible topology | Minor TM edits without clear insertion outcome |
| `extracellular_altered_segment_accessibility` | Altered segment plausibly extracellular/luminal and exposed | Altered segment likely buried, cytoplasmic, or absent from extracellular side | Exposure unknown |
| `folding_stability_qc_escape` | Compatible with maturation, glyco/disulfide support, and no strong QC-failure evidence | Literature/assay evidence of ER retention, degradation, instability, or likely misfolding | No direct stability evidence |
| `cell_surface_localization` | Surface localization supported by curated localization / assay | Non-surface localization, intracellular retention, secretion, or degradation | Membrane but not clearly surface-localized |
| `antibody_targetability` | Surface-localized with extracellular accessible altered segment | Not surface-localized or altered segment inaccessible | Surface-localized but accessibility unresolved |
| `tm_stable_functional` | Positive on insertion/topology, folding/QC, and surface localization | Negative on any major biogenesis step with evidence of instability/non-function | Mixed task outcomes |
| `tm_non_stable_non_functional` | Strong evidence of insertion failure, QC failure, non-surface localization, or instability | Strong evidence of stable functional membrane behavior | Mixed task outcomes |

### TF task labels

| Task | Positive label | Negative label | Ambiguous / exclude |
| --- | --- | --- | --- |
| `nuclear_localization` | Nuclear localization supported by assay/annotation | Cytoplasmic or non-nuclear localization supported | Dual localization without clear dominant state |
| `dbd_integrity_specificity_shift` | DBD intact or specificity shift experimentally supported but still DNA-competent | DBD disrupted or DNA-binding lost | DBD altered but effect unresolved |
| `cofactor_ppi_rewiring` | No major cofactor rewiring or experimentally preserved interactions | Clear rewiring/loss of cofactor interaction network | Partial interaction changes |
| `activation_repression_competence` | Can still activate/repress transcription as expected | Lacks transcriptional effect or loses expected regulatory mode | Weak or context-specific activity |
| `dominant_negative_neomorphic_behavior` | Assay-supported DN/neomorphic behavior | No evidence of DN/neomorphic behavior and preserved canonical function | Suspected but not validated |
| `overall_regulatory_functionality` | Integrates nuclear localization, DBD competence, and regulatory output | Fails core localization/binding/regulatory output | Mixed behavior |
| `tf_non_stable_non_functional` | Clear evidence of failed localization, failed DNA binding, or no regulatory output | Strong evidence of preserved regulatory functionality | Mixed behavior |

## Supervision sources

### Membrane

Direct positive/negative sources:

- UniProt subcellular localization
- UniProt topological domain, signal, TM, glycosylation, disulfide features
- HPA localization
- curated literature cases of ER retention, secretion, degradation, or failed insertion
- future orthogonal surface proteomics / flow / antibody staining if available

### TF

Direct positive/negative sources:

- precomputed hESC TF isoform assays
  - TF isoform PPI
  - TF isoform PDI
  - TF isoform transcriptional activation
  - TF isoform cellular localization
  - TF isoform perturbation / DEG response
- UniProt DNA-binding and localization annotations
- curated TF-family literature where available

## Benchmark design

### Splits

Retain gene-disjoint holdouts as the primary evaluation design:

- `train`
- `val`
- `test`

Hard rule:

- no transcript from the same gene should cross splits

### Metrics

Primary:

- `average_precision`
- `AUROC`

Secondary:

- `balanced_accuracy`
- calibration error / Brier score
- thresholded precision/recall at operating points

For aggregate labels:

- macro AP across family-specific aggregate heads
- confusion matrix for:
  - `stable_functional`
  - `indeterminate`
  - `non_stable_non_functional`

### Benchmark ladder

Run in this order:

1. `logistic` sanity baseline
2. `hist_gradient_boosting`
3. `xgboost`
4. Phase C class-specific autoencoder
5. DeepImmuno-style CNN
6. DeepImmuno-style CNN+autoencoder
7. Phase D family-specific multitask model

This mirrors prior best practice:

- tree ensembles were strongest overall in Phase C
- deep models are kept as comparators and may win once labels become more specific

### Ablations

Required ablations:

- drop `base_reference_delta`
- drop `segment_change`
- drop `typed_overlap`
- drop terminal/full-sequence physicochemical blocks
- drop motif block
- drop interface-specific features for TF
- drop TM/topology-specific features for membrane tasks

Family-specific ablations:

- membrane:
  - no signal features
  - no topology features
  - no extracellular maturation features
- TF:
  - no NLS/NES features
  - no DNA-interface features
  - no activation/repression region features
  - no cofactor/PPI rewiring features

## Deliverables

Phase D should ultimately output:

- per-head probabilities or scores
- per-family aggregate stable/functional vs non-stable/non-functional labels
- holdout benchmark reports per head
- aggregate benchmark report per family
- calibrated inference outputs for query isoforms

## Current proxy benchmark outcome

The Phase D proxy benchmark has now been executed end-to-end. Outputs live in:

- `Daedalus/phase_d/checkpoints/phase_d_multitask_metrics.json`
- `Daedalus/phase_d/checkpoints/phase_d_multitask_report.md`
- `Daedalus/phase_d/checkpoints/phase_d_model_comparison_metrics.tsv`
- `Daedalus/phase_d/checkpoints/phase_d_model_comparison_report.md`

What this currently validates:

- the new Phase D task heads are wired correctly
- the benchmark tables and holdout splits build correctly
- the family-aware multitask model trains correctly
- the comparison harness runs correctly

What it does not yet validate:

- biological correctness of the Phase D tasks
- realistic discriminative performance on external assay labels

Reason:

- many proxy tasks are extremely one-sided after thresholding
- several tasks have zero negatives or zero positives in the test split
- those tasks produce trivially perfect metrics and are not scientifically
  informative

So the current proxy benchmark should be treated as:

- engineering validation
- not final scientific validation
