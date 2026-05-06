# training

Trains the twelve base classifiers (six classical + six residual-MLP
autoencoder variants) under gene-grouped 5-fold cross-validation, then
combines their out-of-fold predictions with a stacked logistic
meta-learner and a K-of-N consensus ensemble.

## Inputs

- `benchmark/data/uniprot_isoform_benchmark_instances.parquet`
- `benchmark/data/uniprot_isoform_benchmark_schema.json`

## Outputs (in `checkpoints/`)

| Path | Description |
|---|---|
| `kfold_classical/tm_isoform_kfold_predictions.tsv` | Per-isoform out-of-fold probabilities for the six classical models (logistic, elasticnet, adaboost, hist_gradient_boosting, xgboost, random_forest) |
| `kfold_classical/tm_isoform_kfold_fold_metrics.tsv` | Per-fold AUROC, AP, sensitivity, specificity, threshold |
| `kfold_classical/tm_isoform_kfold_summary.{tsv,json}` | Cross-fold aggregated metrics |
| `kfold_residual_mlp/...` | Same outputs for the six autoencoder variants (autoencoder, wide, deep, highdrop, strongrecon, narrow) |
| `stacking/meta_oof_predictions.tsv` | Stacked meta-learner OOF probabilities |
| `stacking/meta_threshold_sweep.tsv` | Per-threshold confusion-matrix sweep on the meta predictions |
| `stacking/meta_coefficients.tsv` | Logistic-regression coefficients of the meta-learner over the 12 base-model probabilities |
| `stacking/summary.json` | Selected operating threshold + headline metrics |
| `ensemble_consensus.tsv` | K-of-N consensus ensemble metrics for K=1..6 |
| `BENCHMARK_REPORT.md` | Final per-model TP/FN/TN/FP table at sens=85/90/95% and deployment recommendation |

## Reproduction

```bash
cd components/Daedalus

# Classical base models
python -m Daedalus.training.scripts.run_kfold_benchmark \
    --models phase_d_logistic,phase_c_elasticnet,phase_c_adaboost,phase_d_hist_gradient_boosting,phase_d_xgboost,phase_c_random_forest \
    --tune-threshold --sensitivity-floor 0.85 \
    --out-dir training/checkpoints/kfold_classical

# Autoencoder variants (slower; ~10 minutes for 6 variants × 5 folds on CPU)
python -m Daedalus.training.scripts.run_kfold_benchmark \
    --models phase_c_residual_mlp_autoencoder,residual_mlp_ae_wide,residual_mlp_ae_deep,residual_mlp_ae_highdrop,residual_mlp_ae_strongrecon,residual_mlp_ae_narrow \
    --tune-threshold --sensitivity-floor 0.85 --max-epochs 30 --patience 6 \
    --out-dir training/checkpoints/kfold_residual_mlp

# Stacking + consensus on combined OOF
python -m Daedalus.training.scripts.stacking_meta_learner
python -m Daedalus.training.scripts.ensemble_consensus \
    --out training/checkpoints/ensemble_consensus.tsv

# Confusion matrices at threshold 0.5 (optional sanity check)
python -m Daedalus.training.scripts.summarize_negative_detection \
    --kfold-dir training/checkpoints/kfold_classical
```

## Files

| File | Purpose |
|---|---|
| `scripts/run_kfold_benchmark.py` | Main k-fold runner with threshold tuning |
| `scripts/_model_comparison_lib.py` | Shared library — prepare_arrays, classical-prob helper, deep-train helper, threshold tuner. Imported by `run_kfold_benchmark.py`. |
| `scripts/stacking_meta_learner.py` | Logistic-regression meta-learner over OOF probabilities |
| `scripts/ensemble_consensus.py` | K-of-N consensus ensemble |
| `scripts/summarize_negative_detection.py` | Per-model TP/FN/TN/FP summary at fixed threshold |

## Production checkpoint

`training/checkpoints/` holds the **strict-canonical-PM** (Round 10) results
that are referenced in the grant. Earlier hyperparameter-iteration rounds
have been archived to `archive/intermediate_rounds.tar.gz`.
