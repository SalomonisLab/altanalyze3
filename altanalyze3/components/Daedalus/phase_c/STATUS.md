# Daedalus Phase C Status

## Current conclusion

The strongest currently implemented Phase C classifier is:

- rich biochemical feature matrix
- XGBoost classifier

The class-specific deep autoencoder remains implemented and benchmarked, but it
is not the default recommendation at the current stage.

## Dataset used for the current Phase C benchmark

- `40,000` reference-vs-alternative isoform pairs in the biochemical matrix
- `75,337` multitask instances after expansion into:
  - `global`
  - `membrane`
  - `surface`
  - `kinase`
  - `transcription_factor`

This is a subset benchmark for the richer Phase C biochemical pipeline, not the
full historical Phase B corpus.

## Implemented Phase C components

- segment-aware reference-vs-alternative changed-region detection
- typed biochemical overlap features from UniProt regions
- sequence-derived physicochemical descriptors from `peptides`
- terminal and full-sequence motif/physicochemical summaries
- class-specific expert autoencoder with latent bottleneck
- fast comparison baselines with logistic regression and XGBoost
- feature-block ablation benchmarking

## Benchmark commands

Build the Phase C task table:

```bash
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/build_multitask_instances.py
```

Train the deep model:

```bash
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/train_biochem_multitask.py \
  --max-epochs 40 \
  --patience 8 \
  --batch-size 1024
```

Run model comparison:

```bash
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/run_model_comparison.py
```

Run feature-block ablations:

```bash
Daedalus/.venv/bin/python Daedalus/phase_c/scripts/run_feature_ablation.py
```

## Current benchmark summary

Phase C subset test metrics:

- XGBoost:
  - `global` AP `0.878873`, AUROC `0.914486`
  - `membrane` AP `0.886331`, AUROC `0.903270`
  - `surface` AP `0.878424`, AUROC `0.895676`
  - `kinase` AP `0.881222`, AUROC `0.919501`
  - `transcription_factor` AP `0.874922`, AUROC `0.902396`

- Logistic:
  - `global` AP `0.821485`, AUROC `0.869970`
  - `membrane` AP `0.812549`, AUROC `0.850795`
  - `surface` AP `0.806671`, AUROC `0.845223`
  - `kinase` AP `0.761381`, AUROC `0.815748`
  - `transcription_factor` AP `0.827372`, AUROC `0.865278`

- Class-specific biochemical autoencoder:
  - `global` AP `0.799552`, AUROC `0.853325`
  - `membrane` AP `0.808293`, AUROC `0.836624`
  - `surface` AP `0.803444`, AUROC `0.831492`
  - `kinase` AP `0.783835`, AUROC `0.847317`
  - `transcription_factor` AP `0.845512`, AUROC `0.881434`

## Interpretation

- The richer biochemical feature layer clearly adds usable signal.
- XGBoost is currently the best-performing efficient classifier.
- The deep autoencoder is not yet the best general classifier on this feature
  matrix.
- The deep model is still worth retaining as an experimental path because it is
  competitive for the TF task and gives a dimensional bottleneck for later
  encoder work.

## DeepImmuno-style benchmark status

A second Phase C benchmark path is now implemented to mirror the key logic from
DeepImmuno:

- systematic ML/DL comparison
- PCA dimensional reduction of the high-dimensional physicochemical blocks
- explicit CNN and CNN+autoencoder classifiers

Implemented files:

- `Daedalus/phase_c/models/deepimmuno_style.py`
- `Daedalus/phase_c/scripts/run_deepimmuno_style_benchmark.py`

Tested models in this benchmark path:

- `elasticnet`
- `knn`
- `svm_rbf`
- `random_forest`
- `adaboost`
- `hist_gradient_boosting`
- `residual_mlp_autoencoder`
- `deepimmuno_cnn`
- `deepimmuno_cnn_autoencoder`

Completed tasks:

- `global`
- `membrane`
- `surface`
- `kinase`
- `transcription_factor`

Best model by task:

- `global`: `hist_gradient_boosting` AP `0.878693`, AUROC `0.916296`
- `membrane`: `hist_gradient_boosting` AP `0.883712`, AUROC `0.906119`
- `surface`: `hist_gradient_boosting` AP `0.886910`, AUROC `0.903786`
- `kinase`: `hist_gradient_boosting` AP `0.874902`, AUROC `0.909127`
- `transcription_factor`: `hist_gradient_boosting` AP `0.870329`, AUROC `0.899536`

Strongest deep model by task:

- `global`: `deepimmuno_cnn_autoencoder` AP `0.826830`, AUROC `0.878698`
- `membrane`: `deepimmuno_cnn` AP `0.832216`, AUROC `0.862649`
- `surface`: `deepimmuno_cnn` AP `0.814196`, AUROC `0.852984`
- `kinase`: `deepimmuno_cnn_autoencoder` AP `0.796933`, AUROC `0.857060`
- `transcription_factor`: `deepimmuno_cnn_autoencoder` AP `0.850442`,
  AUROC `0.891157`

PCA reduction ranges across tasks:

- `full_sequence_physchem` `74 -> 23-24`, explained variance `0.9455-0.9509`
- `terminal_sequence_physchem` `120 -> 24`, explained variance `0.8852-0.8958`
- `typed_overlap` `51 -> 16`, explained variance `0.7188-0.7747`

Current interpretation:

- PCA reduction on the physicochemical blocks is viable and implemented
- the DeepImmuno-style deep models are working and benchmarked
- on all five completed tasks, `hist_gradient_boosting` and `random_forest`
  outperform the deep CNN variants
- among the deep models, `deepimmuno_cnn_autoencoder` is strongest overall,
  but `deepimmuno_cnn` is better on `membrane` and `surface`
- the stable aggregate report is:
  - `Daedalus/phase_c/checkpoints/deepimmuno_style_benchmark_all_tasks.md`

## Ablation takeaways

Strongest useful blocks:

- `base_reference_delta`
- `segment_change`

Weak or mixed blocks in the current implementation:

- `typed_overlap`
- `terminal_sequence_physchem`
- `full_sequence_physchem`
- `motifs`

More specifically:

- dropping `base_reference_delta` causes large performance collapse for every
  task
- dropping `segment_change` consistently hurts every task
- dropping `typed_overlap` often improves membrane, surface, and TF logistic
  baselines, so the current typed-overlap implementation is too noisy
- dropping terminal/full-sequence physicochemical blocks often improves
  membrane, surface, and kinase logistic baselines, suggesting those blocks
  need pruning or more biologically constrained redesign

## Output files

- `Daedalus/phase_c/checkpoints/biochem_multitask_metrics.json`
- `Daedalus/phase_c/checkpoints/biochem_multitask_report.md`
- `Daedalus/phase_c/checkpoints/model_comparison_metrics.tsv`
- `Daedalus/phase_c/checkpoints/model_comparison_report.md`
- `Daedalus/phase_c/checkpoints/feature_ablation_metrics.tsv`
- `Daedalus/phase_c/checkpoints/feature_ablation_report.md`
