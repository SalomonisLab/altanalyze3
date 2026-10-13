# rna2lipid

`rna2lipid` predicts lipid abundances from RNA expression. The module supports
two model architectures and one prediction interface.

The complete-panel native-log2 release is packaged as
`rna2lipid_hs_lung_native_log2_20261005_bundle.pkl`. It retains all **202 lipid
outputs**, all **1,303 RNA inputs**, and the original lipidwise **ElasticNetCV**
training procedure. It uses the **45 authorized non-D071 training profiles** and
per-lipid median filling of the 25 missing native-log2 target entries. The release
is bound to its exact validated pickle hash by `release.py` and
`native_log2_release.json`; promotion was requested by the user on 2026-10-06.
Activation is staged pending the historical API snapshot integrity decision.

The corrected outputs are **native relative log2 abundances**. Negative values
are valid and are preserved. The linear inverse is **2**raised to the prediction,
with no subtraction of 1 and no new pseudocount. These are relative abundance
estimates, not calibrated concentrations or a complete 202-lipid MS1 calibration.
A ratio of arithmetic mean linear abundances supplies the reported fold change.
BH uses the prespecified full 202-lipid panel for lipid differentials. See
[the complete-panel comparison](artifacts/LungMAP_full202_BH_corrected_20261006/EVALUATION_SUMMARY.md).

The earlier 47-lipid ElasticNet and kernel-ridge artifacts below are historical
restricted/alternative candidates. Neither is the new release or its deployment
source. The former 50-profile production bundle remains intact for rollback as
`api.PREVIOUS_BUNDLE_PATH`.

For the completed comparison of bulk and sorted-cell normalization, combined
reference construction, differential-preservation checks, and a new candidate
predictor, see [the reference normalization workflow](REFERENCE_NORMALIZATION_WORKFLOW.md).
Its prediction command defaults to the 43 lipids meeting provisional internal
support criteria. It does not replace the production default described below.

The corrected all-reference lung-RNA/lung-lipid comparison is documented in the
[IPF candidate validation report](artifacts/IPF_candidate_validation_20261001/IPF_CANDIDATE_REPORT.md).
The old plasma-based IPF analysis is archived under `delete/IPF_plasma_20261001`.
Candidates remain separate from the production default. Reproduce with `validate_ipf_candidate.py`.

The full LungMAP evaluation candidate includes all 26,639 sample/state pseudobulks
and 230,057 metacells, seven saved IPF/control contrasts, and unchanged production
comparators. See [dataset/cell-type selection](artifacts/LungMAP_IPF_candidate_20261004/DATASET_CELLTYPE_SELECTION.md),
[both-significant overlap and fold thresholds](artifacts/LungMAP_IPF_candidate_20261004/BOTH_SIGNIFICANT_OVERLAP_REPORT.md),
and [exact RNA gene explanations](artifacts/LungMAP_IPF_candidate_20261004/RNA_GENE_EXPLANATIONS.md).
For the union of concordant lipids over any cell-type combination, separately
within each cohort, see [concordant coverage](artifacts/LungMAP_IPF_candidate_20261004/COHORT_CONCORDANT_COVERAGE.md).
This counts a lipid once if any state agrees, regardless of opposing states,
and lists every supporting cell type plus minimum covering combinations.
The [full numerical report](artifacts/LungMAP_IPF_candidate_20261004/LungMAP_IPF_CANDIDATE_REPORT.md)
retains all lipid tables, donor counts, whole-lung donor validation, and saved-result
reproduction errors. Dataset/state/threshold selection is exploratory; opposing
cell states are hypotheses rather than experimentally validated cellular effects.
Reproduction scripts are `evaluate_lungmap_ipf.py`, `evaluate_lungmap_whole_lung.py`,
`rank_lungmap_ipf_comparisons.py`, `explain_lungmap_ipf_candidate.py`,
`evaluate_lungmap_overlap_thresholds.py`, and `render_lungmap_ipf_evaluation.py`.
Reproduce the cohort unions with `summarize_lungmap_concordant_coverage.py`.
Copy-ready Adams, Natri and intersected Results/Methods drafts are in
[Adams–Natri results](artifacts/LungMAP_IPF_candidate_20261004/ADAMS_NATRI_RESULTS_AND_METHODS.md).
`write_adams_natri_results.py` also exports cross-cohort lipid intersections
and shared-cell-type support at every tested threshold.
The [best concordant cell types](artifacts/LungMAP_IPF_candidate_20261004/BEST_CONCORDANT_CELLTYPES.md)
are ranked solely by concordant counts, retaining ties and showing donor counts.
Reproduce with `summarize_lungmap_best_celltypes.py`.

| Architecture | Model | Default? |
| --- | --- | --- |
| `sparse_lipidwise` | One `ElasticNetCV` per lipid, over that lipid's top correlated genes | yes |
| `multitask` | One `MultiTaskElasticNetCV` over all lipids | no, kept for reproducibility |

`api.load_bundle` reads both. It selects the path from the bundle's own keys:
a bundle with `models` is lipid-wise, a bundle with `model` is multi-task.

## The default bundle

`api.DEFAULT_BUNDLE_PATH` selects
`rna2lipid_hs_lung_native_log2_20261005_bundle.pkl`.

- Architecture: 202 separate ElasticNetCV models with the original training settings.
- Inputs: 1,303 gene symbols. Outputs: all 202 production lipids in native relative log2 units.
- Training: 45 sorted profiles from nine lung donors, five sorted populations;
  D071 exclusion is explicitly authorized. No bulk training profiles were added.
- Target/RNA StandardScaler, Pearson ranking, top-N selection and ElasticNetCV
  grids/CV/seeds match the established procedure. No retraining occurs on promotion.
- The release has no external donor-held-out validation claim. Its disease-direction
  evaluation uses the supplied independent, unmatched lung-tissue IPF lipidomics.

`rna2lipid_hs_lung_lipidwise_bundle.pkl` is the previous 50-profile default, kept
byte-for-byte as `api.PREVIOUS_BUNDLE_PATH`. Its historical validation metrics and
constant TG(51:1) model describe that older bundle, not the promoted native-log2 one.
The independent original baseline contract remains 50/202/1,303 for provenance.

The source-confirmed [bulk unmatched IPF comparison](artifacts/bulk_IPF_CPM_harmonized_20261008/BULK_IPF_CANDIDATE_VS_PRIOR_REPORT.md)
uses **log2(CPM + 1)**, matching the training normalization established by the
supplied script and the user’s 2026-10-08 file-identity confirmation. Both fixed
202-lipid models receive identical 1,303-gene inputs and approved median fills.
At BH <0.1 in both assays with no additional fold cutoff, direction agreement
is 48/70 (68.6%) for the current model versus 23/61 (37.7%) for the prior.
On the same 41 lipids significant in the source and both models, agreement is
26/41 versus 17/41. Abby’s top-30 source-only sign metric gives 18/30 (60%) versus
12/30 (40%) after the established donor aggregation. The previous CP10k report
remains saved as a scale-sensitivity result.

Code, training inputs, source hashes and executed fold-space verification are
collected in [the normalization/training directory](workflows/lung_native_log2_20261008/README.md).
The current model uses workbook **relative log2** targets; MASSIVE raw validation
and its historical partial MS1 calibration are documented separately. The
source-confirmed CPM normalization differs from the deployed CP10k input path;
this package/evaluation does not change that deployed path.

The [fresh three-version bulk comparison](artifacts/bulk_IPF_three_version_comparison_20261008/BULK_IPF_THREE_VERSION_REPORT.md)
compares Abby’s supplied model with CPM inputs, the corrected model with its
last CP10k inputs, and that same corrected model with source-confirmed CPM
inputs. All 45 p/fold-threshold rows, every lipid and discordant list, the
same-lipid comparison, and source-only top-30 concordance are included.
The fresh rerun exactly reproduced previous full-panel predictions/statistics.
Abby’s original author-reported ~50% lacks saved outputs and is identified
separately from the supplied-model rerun under the common donor-level methods.

## The prior bundle

`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/newnormelastic_multitask_try.pkl`

One `MultiTaskElasticNetCV`, trained on 68 profiles (bulk and cell-type), same
1,303 genes and same 202 lipids. `api.LEGACY_MULTITASK_BUNDLE_PATH` points at
it. It reproduces every result published before 2026-08-30.

## Environment

`/usr/bin/python3` on this Mac runs the module: scikit-learn 1.3.0, numpy
1.24.3, pandas 2.0.3, anndata 0.9.2. Run every command from the repository root
with `PYTHONPATH` set:

```bash
cd /Users/saljh8/Documents/GitHub/altanalyze3
export PYTHONPATH=/Users/saljh8/Documents/GitHub/altanalyze3
```

## Inspect the bundle

```bash
/usr/bin/python3 -m altanalyze3.components.rna2lipid.cli model-info
```

## Train a model for a new tissue

Copy the template, fill it in, and run `train`. You change no code.

`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/configs/template_new_tissue.json`

The config declares:

1. **`data.tables`** — a list of paired RNA and lipid matrices. Each pair must
   share sample identifiers. `rna_orientation` and `lipid_orientation` say which
   axis holds the samples: `samples_by_features` or `features_by_samples`.
2. **`sample_metadata`** — how to read the donor and the group out of a sample
   identifier. `donor` drives the donor-disjoint holdout. `group` drives the
   constrained holdout. Strategies: `split_field`, `regex`, `constant` for the
   fields, and `label`, `keyword_map`, `regex` for the group.
3. **`keep_groups`** — the groups to train on. The builder drops every other
   group, and the manifest counts what each dropped group cost.
4. **`training.architecture`** — `sparse_lipidwise` or `multitask`.
5. **`training.holdout.mode`** — `constrained` (the lung protocol, not
   donor-disjoint), `donor_disjoint`, or `none`.

Then:

```bash
/usr/bin/python3 -m altanalyze3.components.rna2lipid.cli train \
  --config altanalyze3/components/rna2lipid/configs/<your_tissue>.json
```

Training writes, beside the bundle:

- `training_summary.json` — metadata, holdout metrics, resubstitution metrics
- `bundle_manifest.json` — every input path, every count, every retained fraction
- `training_predictions.tsv`, `training_per_lipid_metrics.tsv`
- `model_summary.tsv` — the chosen top-N gene count and alpha per lipid
- `nonzero_coefficients.tsv` — every retained gene and its coefficient
- `candidate_top_gene_models.tsv` — every candidate gene-set size the trainer tried

## Reproduce the lung bundle

`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/configs/lung_lipidwise.json`

**This repository does not carry the two lipid tables that config names.** The
delivered run read `Bulk_lipids_cleaned_normalized_median_527_log2.csv` and
`cell_lipids_cleaned_norm_median_286_log2.csv`. Both files live on the
collaborator's cluster under `/users/ramkd9/Lipid_Predict/`. Place both under
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/data/`
before running that config. The repository does carry the two RNA tables and
two lipid tables on a different normalization; those cannot reproduce the
shipped bundle.

## Evaluate

`evaluate` runs donor-held-out `LeaveOneGroupOut` validation against a mean
baseline and a PCA plus ridge baseline. It fits whichever architecture the
config names.

```bash
/usr/bin/python3 -m altanalyze3.components.rna2lipid.cli evaluate \
  --bundle altanalyze3/components/rna2lipid/rna2lipid_hs_lung_lipidwise_bundle.pkl \
  --config altanalyze3/components/rna2lipid/configs/lung_lipidwise.json \
  --output-dir <output directory> \
  --skip-external
```

## Predict

From a sample-by-gene CSV:

```bash
/usr/bin/python3 -m altanalyze3.components.rna2lipid.cli predict-csv \
  --input altanalyze3/components/rna2lipid/data/feature_blankreduiction.csv \
  --output <output csv> \
  --summary-json <output json>
```

Add `--transpose` when the matrix is gene-by-sample.

From a single-cell h5ad, per cell or averaged by an `.obs` column:

```bash
/usr/bin/python3 -m altanalyze3.components.rna2lipid.cli predict-h5ad \
  --input <sample.h5ad> \
  --groupby cell_type \
  --output <output csv>
```

The model is linear and the prediction path is affine, so averaging per-cell
model-scale predictions equals predicting from the group's mean expression
vector. `predict-h5ad --groupby` retains this established API behavior. scALABLE
pseudobulk differentials use arithmetic abundance means, described below.

## scALABLE

scALABLE serves the lipid modality through
`altanalyze3/components/cellHarmony/flask/pipeline.py:_build_imputed_lipid_adata`.
All five human lung references declare the bundle explicitly in
`altanalyze3/components/cellHarmony/flask/reference_config.json`:

```json
"lipids": {
  "bundle_path": "../../rna2lipid/rna2lipid_hs_lung_native_log2_20261005_bundle.pkl",
  "expression_scale": "native_relative_log2",
  "log_base": "2",
  "log_pseudocount": 0
}
```

A new tissue points `bundle_path` at its own bundle. You change no code. When a
reference names no `lipids` bundle, the builder loads `api.DEFAULT_BUNDLE_PATH`,
which is the native-log2 lung bundle.

The job log records which bundle actually ran:

```
[params] rna2lipid resolved bundle=<path> architecture=lipidwise lipids=202 matched_genes=<n>/1303
```

## Provenance

`provenance/train_sparse_lipidwise_delivered_2026-08-19.py` is the training
script as delivered. `provenance/api_delivered_2026-08-19.py` is the prediction
interface as delivered. Neither file is the supported runtime path. Both files
stay here so you can trace the shipped bundle. The older `clean_*.py` scripts
hold research provenance only.

### Cell and pseudobulk comparisons

The user can choose cells or pseudobulks. Cell comparisons use individual stored
predictions. For pseudobulks, scALABLE averages **existing linear predictions**
within each selected sample and cell state, then restores the declared log scale.
It does not run the model a second time or normalize to the feature-panel total.
For the new lung lipids this is `log2(mean(2**prediction))`. The shared aggregation
path supports verified ADT log1p encodings and linear GRN/TF scores, preserving
signed scores. Ambiguous AML target log2 encodings stop and request the target
preprocessing protocol; their inverse is not guessed. Existing cell sampling,
minimum group sizes and test settings are retained. Pseudobulk tests use sample
profiles as replicates; cell-level tests do not establish donor independence.

scALABLE converts declared natural-log CP10k RNA to log2 CP10k before new lung
lipid inference, as previously authorized. Already declared log2 RNA is unchanged.
Standalone numeric prediction calls retain their supplied-input behavior; an RNA
scale is never inferred from value ranges. The source pickle and all coefficients
are identical to the validated candidate. Output summaries record the actual
model version, release, scale and RNA conversion.

### Downstream scale and fold contract

Export, aggregation and differential analysis share
`cellHarmony/imputed_scale.py`. Imputed data must declare their source encoding;
value ranges do not select a log transform or CP10k normalization. The source
`X` is authoritative for folds, including when stale `raw` or `counts` layers
are present.

| Declared source encoding | Linear value used for abundance means |
| --- | --- |
| `native_relative_log2` | `2**X` |
| `log1p`, with declared base `b` | `b**X - 1` (natural log uses `expm1(X)`) |
| `log2p1` | `2**X - 1` |
| `log2`, with declared source pseudocount `c` | `2**X - c` |
| `linear` | `X` unchanged |

Every reported imputed log2 fold is
`log2(mean(case_linear) / mean(control_linear))`. No new pseudocount is added
to the group means. For example, means of 0.02 and 0.04 yield a twofold change;
adding 1 to these means would instead report only 1.0196-fold. Plain-log targets
use log-sum-exp for stable arithmetic means at very low baselines. Logged
pseudobulks are restored to the source encoding after averaging existing linear
predictions within sample/state. Sample profiles remain the pseudobulk replicates.

After the explicit user response on 2026-10-06 that the 0.1 filter should not
exist, **all imputed modalities use their complete supplied feature panel for
BH**, including low values and negative plain-log values. The lung lipid panel
also requires exact ordered agreement with all 202 release lipids. Nonfinite
raw p values stop the run rather than shrink its correction family. RNA keeps
its existing eligibility policy.

The raw test algorithms remain unchanged, including the existing moderated
pooled/global and co-regulated tests in cell comparisons, as the user requested.
Their inputs retain the declared source scale. Scanpy's approximate fold field
is preserved separately for audit and replaced by the exact abundance fold in
the stored results. Per-population output includes linear means and
`fold_status`. One zero mean gives an infinite directional fold; two zero means
give an undefined fold. Signed linear scores remain preserved in exports and
aggregation, but abundance-fold analysis stops for them instead of shifting or
clipping them. Current GRN predictions are nonnegative edge scores; TF activity
is their existing outgoing-edge sum, so these folds describe relative activity.

Each scALABLE differential run records `prediction_scale_audit.json` with the
source base/pseudocount, fold formula, zero policy and complete-panel BH policy.
An unresolved bare `log2` target can be displayed on its supplied scale, but
has no invented linear layer and is blocked from abundance-fold analysis.
AML target preprocessing remains pending; no inverse has been inferred for it.
Historical API snapshot approval and the unchanged integrity gate still prevent
declaring deployment activated. No scientific dataset was rerun for these code
checks.

Regression tests use known synthetic abundances across six encodings, exact
twofold changes, zero means, preserved raw tests, complete-panel BH, multi-block
disk/RAM agreement, cell and pseudobulk comparisons, and the actual local web
GRN differential routes. See `tests/test_imputed_differential_scale.py` at the
repository root and `validation/test_native_log2_deployment.py`.
