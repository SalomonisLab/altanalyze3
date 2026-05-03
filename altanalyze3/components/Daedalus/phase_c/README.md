# Daedalus Phase C

Phase C is the current segment-aware biochemical modeling path beyond the Phase
B reference-delta baseline.

## Current conclusion

Use the Phase C rich biochemical matrix with `XGBoost` as the current default
classifier.

Keep the class-specific biochemical autoencoder as the deep-learning comparison
path and dimensional-reduction experiment, not the default production model.

There is now also a separate DeepImmuno-style benchmark path with:

- PCA reduction on the physicochemical blocks
- systematic ML/DL comparison
- CNN and CNN+autoencoder classifiers

## What Phase C adds

- explicit reference-vs-alternative changed-segment detection
- sequence-derived physicochemical encodings using `peptides`
- motif-aware and terminal-aware sequence summaries
- typed overlap of changed regions with UniProt feature classes
- class-specific expert branches for:
  - membrane
  - kinase
  - transcription_factor
  - other
- an autoencoder-style latent bottleneck plus task-specific supervised heads

## Build the Phase C matrix

```bash
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/build_biochem_delta_matrix.py
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/build_multitask_instances.py
```

## Train the deep model

```bash
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/train_biochem_multitask.py \
  --max-epochs 40 \
  --patience 8 \
  --batch-size 1024
```

## Benchmark models

```bash
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/run_model_comparison.py
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/run_feature_ablation.py
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/run_deepimmuno_style_benchmark.py \
  --tasks global membrane surface kinase transcription_factor
```

## Current benchmark summary

On the current Phase C subset benchmark:

- `phase_c_xgboost` is the best classifier on all five tasks
- `phase_c_autoencoder` is only competitive for the TF task
- `segment_change` features are consistently useful
- the current `typed_overlap` and several terminal/full-sequence physicochemical
  blocks are too noisy and need redesign

On the completed DeepImmuno-style five-task benchmark:

- `hist_gradient_boosting` ranks first on all five tasks
- `random_forest` ranks second on all five tasks
- strongest deep model by task:
  - `global`: `deepimmuno_cnn_autoencoder` AP `0.826830`, AUROC `0.878698`
  - `membrane`: `deepimmuno_cnn` AP `0.832216`, AUROC `0.862649`
  - `surface`: `deepimmuno_cnn` AP `0.814196`, AUROC `0.852984`
  - `kinase`: `deepimmuno_cnn_autoencoder` AP `0.796933`, AUROC `0.857060`
  - `transcription_factor`: `deepimmuno_cnn_autoencoder` AP `0.850442`,
    AUROC `0.891157`
- PCA reduction is implemented and stable across tasks:
  - `full_sequence_physchem` `74 -> 23-24`
  - `terminal_sequence_physchem` `120 -> 24`
  - `typed_overlap` `51 -> 16`
- consolidated report:
  - `Daedalus/phase_c/checkpoints/deepimmuno_style_benchmark_all_tasks.md`

See:

- `Daedalus/phase_c/STATUS.md`
- `Daedalus/phase_c/checkpoints/model_comparison_report.md`
- `Daedalus/phase_c/checkpoints/feature_ablation_report.md`
- `Daedalus/phase_c/checkpoints/deepimmuno_style_benchmark_all_tasks.md`

## Outputs

- `data/processed/biochem_delta_matrix.parquet`
- `data/processed/biochem_multitask_instances.parquet`
- `checkpoints/biochem_multitask_metrics.json`
- `checkpoints/biochem_multitask_report.md`
- `checkpoints/model_comparison_metrics.tsv`
- `checkpoints/model_comparison_report.md`
- `checkpoints/feature_ablation_metrics.tsv`
- `checkpoints/feature_ablation_report.md`
- `checkpoints/deepimmuno_style_benchmark_metrics.tsv`
- `checkpoints/deepimmuno_style_benchmark_report.md`
- `checkpoints/deepimmuno_style_benchmark_summary.json`
- `checkpoints/deepimmuno_style_benchmark_all_tasks.tsv`
- `checkpoints/deepimmuno_style_benchmark_all_tasks.md`
