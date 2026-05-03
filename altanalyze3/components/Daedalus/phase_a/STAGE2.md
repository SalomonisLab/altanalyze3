# Daedalus Stage 2

## Environment

Current working interpreter:

- `Daedalus/.venv/bin/python`

Installed stack in the venv:

- `numpy 2.4.4`
- `pandas 3.0.2`
- `scipy 1.17.1`
- `scikit-learn 1.8.0`
- `pyarrow 23.0.1`
- `torch 2.11.0`

`xgboost` remains unresolved at the system level because macOS `libomp` has not been completed. Stage 2 can proceed without it.

## Stage 2 outputs

- `data/processed/sklearn_seed_baseline_metrics.tsv`
- `data/processed/sklearn_seed_baseline_report.md`
- `data/processed/torch_seed_baseline_metrics.tsv`
- `data/processed/torch_seed_baseline_report.md`
- `checkpoints/task_global_seed_demo.pt`
- `checkpoints/task_global_seed_demo.metrics.json`

## Best current structural-only test results

Scikit-learn:

- global: `hist_gbdt_structural`
  - `AP=0.819395`
  - `AUROC=0.860452`
- membrane: `hist_gbdt_structural`
  - `AP=0.812655`
  - `AUROC=0.847125`
- kinase: `hist_gbdt_structural`
  - `AP=0.869574`
  - `AUROC=0.895919`
- transcription factor: `hist_gbdt_structural`
  - `AP=0.844113`
  - `AUROC=0.872020`
- surface: `mlp_structural`
  - `AP=0.496960`
  - `AUROC=0.730348`

Torch tabular MLP:

- global
  - `AP=0.820707`
  - `AUROC=0.860386`
- membrane
  - `AP=0.818902`
  - `AUROC=0.849498`
- kinase
  - `AP=0.861197`
  - `AUROC=0.889448`
- transcription factor
  - `AP=0.838809`
  - `AUROC=0.865624`
- surface
  - `AP=0.400127`
  - `AUROC=0.526740`

## Interpretation

The structural-only benchmark remains viable and is now supported by real `scikit-learn` and `torch` baselines rather than just pure-Python prototypes.

The immediate next model-development target should remain:

- global preservation
- membrane-associated preservation
- kinase-focused preservation
- transcription-factor-focused preservation

The surface/topology task still needs better supervision before it can support a strong method claim.
