# Daedalus benchmark report

Final performance of the production stacked meta-learner and its 12 base
models on the strict-canonical-PM benchmark (3,022 (canonical,
alternative)-isoform pairs across 1,396 unique genes: 2,839
plasma-membrane-positive + 183 non-surface-negative). All metrics are
computed under gene-grouped five-fold cross-validation, with the
operating threshold tuned per-fold on validation OOF probabilities.

## Production deployment

**Stacked meta-learner @ threshold 0.79**

| Metric | Value |
|---|---|
| Sensitivity | **95.09%** |
| Specificity | **73.22%** |
| Positive predictive value (PPV) | 98.30% |
| Negative predictive value (NPV) | 47.85% |
| Accuracy | 93.82% |
| AUROC | 0.914 |
| TP / FN / TN / FP | 2,825 / 146 / 134 / 49 |

Operating-threshold sweep:

| Sens floor | Spec | Acc |
|---|---:|---:|
| 85% | 87.43% | 85.95% |
| 90% | 83.61% | 90.33% |
| **95%** | **73.22%** | 93.82% |

## Per-base-model results @ sens ≥ 95%

Sorted by specificity at the 95% sensitivity floor.

| Rank | Model | AUROC | TP | FN | TN | FP | Sens | Spec | PPV | NPV | Acc |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| — | **stacking_meta_learner** | **0.914** | 2,825 | 146 | 134 | 49 | 95.09% | **73.22%** | 98.30% | 47.85% | 93.82% |
| — | **mean_of_probs (12 models)** | 0.918 | 2,834 | 137 | 134 | 49 | 95.39% | 73.22% | 98.30% | 49.45% | 94.10% |
| 1 | residual_mlp_ae_strongrecon | 0.896 | 2,828 | 143 | 128 | 55 | 95.19% | 69.95% | 98.09% | 47.23% | 93.72% |
| 2 | phase_d_logistic | 0.880 | 2,823 | 148 | 126 | 57 | 95.02% | 68.85% | 98.02% | 45.99% | 93.50% |
| 3 | phase_c_residual_mlp_autoencoder | 0.878 | 2,827 | 144 | 126 | 57 | 95.15% | 68.85% | 98.02% | 46.67% | 93.63% |
| 4 | residual_mlp_ae_deep | 0.892 | 2,824 | 147 | 129 | 54 | 95.05% | 70.49% | 98.12% | 46.74% | 93.63% |
| 5 | residual_mlp_ae_narrow | 0.916 | 2,825 | 146 | 125 | 58 | 95.09% | 68.31% | 97.99% | 46.13% | 93.53% |
| 6 | residual_mlp_ae_wide | 0.887 | 2,827 | 144 | 123 | 60 | 95.15% | 67.21% | 97.92% | 46.07% | 93.53% |
| 7 | residual_mlp_ae_highdrop | 0.905 | 2,828 | 143 | 122 | 61 | 95.19% | 66.67% | 97.89% | 46.04% | 93.53% |
| 8 | phase_d_hist_gradient_boosting | 0.913 | 2,825 | 146 | 113 | 70 | 95.09% | 61.75% | 97.58% | 43.63% | 93.15% |
| 9 | phase_c_adaboost | 0.900 | 2,828 | 143 | 110 | 73 | 95.19% | 60.11% | 97.48% | 43.48% | 93.15% |
| — | phase_d_xgboost | 0.912 | 2,817 | 154 | 125 | 58 | 94.82% | 68.31% | 97.98% | 44.80% | 93.28% |
| — | phase_c_random_forest | 0.927 | 2,821 | 150 | 124 | 59 | 94.95% | 67.76% | 97.95% | 45.26% | 93.37% |
| — | phase_c_elasticnet | 0.792 | 2,892 | 79 | 85 | 98 | 97.34% | 46.45% | 96.72% | 51.83% | 94.39% |

The bottom three rows do not formally qualify for sens ≥ 95% AND spec ≥
60%. xgboost and random_forest each fall short of the 95% sensitivity
threshold by less than half a percentage point at their tuned operating
point; elasticnet has the lowest AUROC and falls below the 60% specificity
floor.

## Consensus ensemble (K of 6)

Six-base-model unanimity vote at threshold 0.5. Members: phase_d_logistic,
phase_c_adaboost, residual_mlp_ae_strongrecon, phase_c_elasticnet,
residual_mlp_ae_narrow, residual_mlp_ae_wide.

| Rule | TP | FN | TN | FP | Sens | Spec | Acc | Qualifies (sens≥95 spec≥60) |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|
| K≥1 | 2,802 | 37 | 73 | 110 | 98.59% | 39.89% | 95.13% | ✗ |
| K≥2 | 2,764 | 75 | 103 | 80 | 97.34% | 56.28% | 94.87% | ✗ |
| K≥3 | 2,734 | 105 | 123 | 60 | 96.30% | 67.21% | 94.55% | **✓** |
| K≥4 | 2,690 | 149 | 135 | 48 | 94.75% | 73.77% | 93.49% | ✗ (sens 94.75) |
| K≥5 | 2,641 | 198 | 143 | 40 | 93.07% | 78.14% | 92.13% | ✗ |
| K≥6 (unanimous) | 2,521 | 318 | 161 | 22 | 88.83% | 87.98% | 88.78% | ✗ |

K≥3 is the only consensus rule that qualifies; the stacked meta-learner
achieves higher specificity (73.22% vs 67.21%) at slightly lower
sensitivity (95.09% vs 96.30%).

## Why label denoising was the decisive lever

The "strict-canonical-PM" filter excludes ~419 genes whose canonical
UniProt entry annotates the canonical at both Cell Membrane and another
major compartment (Cytoplasm / Nucleus / Mitochondrion / Lysosome). These
genes contribute ambiguous reference isoforms whose alternative isoforms
are intrinsically hard to label. Removing them:

- AUROC: **0.81 → 0.91** (mean-of-probs oracle)
- spec @ sens=95%: **47% → 73%** (stacked meta-learner)
- qualifying base models (sens≥95 + spec≥60): **0 → 9 of 12**

Adding ESM2-650M embeddings and DeepLoc-style 10-class probabilities on top
of the strict-canonical filter yielded a marginal gain at the AUROC level
(+0.004 mean-of-probs) but slightly hurt the stacked spec at sens=95%
(73.22% → 70.49%) due to over-parameterization on the 183-negative panel.
Tier 1.1 (strict-canonical filter only) is therefore the production
configuration; the t33 and DeepLoc components remain in the data_ingest
stage but are not consumed by the production benchmark builder.

## Reproduction

See `../README.md`.
