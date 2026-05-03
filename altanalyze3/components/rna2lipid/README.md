# rna2lipid

`rna2lipid` predicts lipid abundances from RNA expression. The package now has
three supported surfaces:

- `predict-*` for inference from CSV/TSV or `.h5ad`
- `train` for reproducible model training from a JSON config
- `evaluate` for donor-held-out validation and external structure checks

The older `clean_*.py` scripts remain as research provenance only. They are not
the supported runtime path.

## Environment

Use `python3` and an isolated environment from the repo root:

```bash
python3 -m venv .venv-rna2lipid
.venv-rna2lipid/bin/pip install -r components/rna2lipid/requirements.txt
```

## Inspect the bundled model

The default prediction bundle is now:

`components/rna2lipid/newnormelastic_multitask_try.pkl`

The older scaled-target bundle remains available at:

`components/rna2lipid/otherelastic_multitask_try.pkl`

Inspect it with:

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli model-info
```

## Reproducible training

The canonical config is:

`components/rna2lipid/configs/default_training.json`

It defines:

- input RNA and lipid matrices
- preprocessing rules
- elastic-net hyperparameters
- default output locations
- evaluation settings

Train a reproducible bundle:

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli train
```

Optional overrides:

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli train \
  --output-bundle /tmp/rna2lipid_retrained.pkl \
  --output-dir /tmp/rna2lipid_training_report
```

Training writes:

- a pickle bundle with `model`, `scaler_x`, `X_columns`, `Y_columns`, and `metadata`
- `scaler_y` when target scaling is enabled; unscaled-target bundles may omit it or carry an unfitted placeholder for backward compatibility
- `training_summary.json`
- `bundle_manifest.json`
- `training_predictions.tsv`
- `training_per_lipid_metrics.tsv`

The current canonical training dataset resolves to `68` matched profiles, `1303`
RNA features, and `202` lipids shared across the bulk and cell-type lipid
tables.

## Evaluation

Run donor-held-out validation plus the external pediatric structure report:

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli evaluate \
  --bundle components/rna2lipid/artifacts/rna2lipid_multitask_enet_reproducible.pkl \
  --output-dir /tmp/rna2lipid_eval
```

For a fast smoke test:

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli evaluate \
  --bundle components/rna2lipid/artifacts/rna2lipid_multitask_enet_reproducible.pkl \
  --output-dir /tmp/rna2lipid_eval_smoke \
  --max-folds 2 \
  --skip-external
```

Evaluation includes:

- donor-held-out internal validation with `LeaveOneGroupOut`
- baseline comparisons against:
  - mean lipid profile
  - PCA + ridge regression
- an external pediatric pseudobulk structure-preservation report from
  `GSE161382`

The external pediatric step is not a direct lipid-truth benchmark. It checks
whether predicted lipids preserve biologically sensible separation across coarse
cell compartments and reports top differential predicted lipids.

## Output scale

`rna2lipid` now supports both historical scaled-target bundles and native-scale
bundles. The default config trains with `"target_scaling": "none"`, which keeps
predictions on the original log-transformed lipid scale so downstream
differential analyses retain realistic fold-change magnitudes.

## Predict from a sample-by-gene CSV

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli predict-csv \
  --input components/rna2lipid/data/feature_blankreduiction.csv \
  --output /tmp/rna2lipid_predictions.csv \
  --summary-json /tmp/rna2lipid_predictions.summary.json
```

If the matrix is gene-by-sample instead, add `--transpose`.

## Predict from a single-cell h5ad

Per-cell prediction:

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli predict-h5ad \
  --input sample.h5ad \
  --output /tmp/rna2lipid_per_cell.csv \
  --summary-json /tmp/rna2lipid_per_cell.summary.json
```

Grouped prediction using an `.obs` column:

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli predict-h5ad \
  --input sample.h5ad \
  --groupby cell_type \
  --output /tmp/rna2lipid_by_cell_type.csv
```

Optional flags:

- `--layer counts`
- `--gene-symbol-col gene_symbols`
- `--chunk-size 4096`

## Notes for single-cell inputs

The bundled model was trained on aggregated sample-level and cell-type-level RNA
profiles, not raw droplet counts. For future `cellHarmony-web` use, the safer
path is:

1. load one uploaded sample `.h5ad`
2. align or annotate cells
3. predict per cell only if exploratory output is needed
4. aggregate predicted lipids by aligned cell state, cell type, or another
   stable `.obs` label

Because the model is linear and the prediction path is affine, averaging
per-cell predictions is equivalent to predicting from the mean expression vector
of that group.
