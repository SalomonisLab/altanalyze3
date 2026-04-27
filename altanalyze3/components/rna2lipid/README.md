# rna2lipid

`rna2lipid` predicts lipid abundances from RNA expression using the bundled
multi-output elastic-net model in
`components/rna2lipid/otherelastic_multitask_try.pkl`.

The original `clean_*.py` scripts in this folder are still useful as research
provenance, but they are not the recommended runtime path. They are analysis
scripts with hard-coded file assumptions, plotting side effects, and no stable
command-line interface.

The supported runtime path is the new CLI:

```bash
python3 -m components.rna2lipid.cli model-info
python3 -m components.rna2lipid.cli predict-csv ...
python3 -m components.rna2lipid.cli predict-h5ad ...
```

## Environment

The system `python` on some machines in this repo still points to Python 2. Use
`python3`.

Create an isolated environment from the repo root:

```bash
python3 -m venv .venv-rna2lipid
.venv-rna2lipid/bin/pip install -r components/rna2lipid/requirements.txt
```

`scikit-learn==1.6.1` is pinned to match the version used to serialize the
bundled model.

## Inspect the bundled model

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli model-info
```

Expected high-level properties:

- `1303` input genes
- `286` output lipids
- `MultiTaskElasticNetCV` model

## Predict from a sample-by-gene CSV

Input format:

- rows are samples
- columns are gene symbols
- the first column is the sample identifier

Example:

```bash
.venv-rna2lipid/bin/python -m components.rna2lipid.cli predict-csv \
  --input components/rna2lipid/data/feature_blankreduiction.csv \
  --output /tmp/rna2lipid_predictions.csv \
  --summary-json /tmp/rna2lipid_predictions.summary.json
```

If your matrix is gene-by-sample instead, add `--transpose`.

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

- `--layer counts` to use `adata.layers["counts"]`
- `--gene-symbol-col gene_symbols` to use `adata.var["gene_symbols"]`
- `--chunk-size 4096` for larger jobs

## Notes for single-cell inputs

The bundled model was trained on aggregated sample-level and cell-type-level RNA
profiles, not raw droplet counts. It can generate per-cell predictions, but the
more defensible integration path for `cellHarmony-web` is:

1. load one uploaded sample `.h5ad`
2. align or annotate cells
3. predict per cell only if exploratory output is desired
4. aggregate predicted lipids by a stable `.obs` field such as aligned cell
   state, cell type, or sample-level subgroup

Because the model is linear and uses affine scaling, averaging per-cell
predictions is equivalent to predicting from the mean expression profile of that
group.

## Current limitations

- The bundled model artifact is present, but the legacy training script in this
  folder does not reproduce the same output shape as the bundled pickle. Treat
  retraining as a separate cleanup task.
- The legacy validation script expects `ElasticNet_model.pkl`, which is not the
  bundled artifact used by the supported CLI.
- No web integration is included here yet. This refactor is intended to make
  later `cellHarmony-web` integration straightforward.
