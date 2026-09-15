# rna2lipid

`rna2lipid` predicts lipid abundances from RNA expression. The module supports
two model architectures and one prediction interface.

| Architecture | Model | Default? |
| --- | --- | --- |
| `sparse_lipidwise` | One `ElasticNetCV` per lipid, over that lipid's top correlated genes | yes |
| `multitask` | One `MultiTaskElasticNetCV` over all lipids | no, kept for reproducibility |

`api.load_bundle` reads both. It selects the path from the bundle's own keys:
a bundle with `models` is lipid-wise, a bundle with `model` is multi-task.

## The default bundle

`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/rna2lipid_hs_lung_lipidwise_bundle.pkl`

- Architecture: 202 separate `ElasticNetCV` models, one per lipid.
- Inputs: 1,303 gene symbols. Outputs: 202 lipids, on a log2 scale.
- Trained on 50 sorted cell-type profiles from 10 human lung donors,
  10 profiles each for EPI, MES, MIC, END and PMX. **0 bulk profiles.**
- Nonzero genes per lipid: minimum 0, median 21, maximum 57.
- 1 of 202 lipids, `TG(51:1)`, has an all-zero model. It predicts a constant.

### Scope limit you must state when you use it

The model saw no bulk profile. On the repository's own bulk RNA table, 16.0% of
predicted values fall outside the training mean ± 3 SD of their lipid. On
cell-type profiles the same figure is 0.15%. Apply the model to cell-state or
pseudobulk profiles. Bulk input is extrapolation.

### The reported holdout number is not donor-held-out

The bundle carries `validation_global_metrics`: Pearson r 0.871, Spearman r
0.864, R² 0.754, RMSE 1.384 over 28 held-out profiles. The split rule caps each
donor at 3 profiles in the holdout set, and each donor contributes 5 profiles,
so **every donor appears in both the training and the holdout set**. Set
`holdout.mode` to `donor_disjoint` when you need a donor-held-out estimate.
See `VALIDATION.md`.

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
predictions equals predicting from the group's mean expression vector.

## scALABLE

scALABLE serves the lipid modality through
`altanalyze3/components/cellHarmony/flask/pipeline.py:_build_imputed_lipid_adata`.
All four human lung references declare the bundle explicitly in
`altanalyze3/components/cellHarmony/flask/reference_config.json`:

```json
"lipids": {
  "bundle_path": "../../rna2lipid/rna2lipid_hs_lung_lipidwise_bundle.pkl",
  "expression_scale": "log2",
  "log_base": "2"
}
```

A new tissue points `bundle_path` at its own bundle. You change no code. When a
reference names no `lipids` bundle, the builder loads `api.DEFAULT_BUNDLE_PATH`,
which is the lung bundle.

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
