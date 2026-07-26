# rna2psi — validation

## What is imputed
Per-event PSI for the 18,168 Leucegene variant/covariate splicing events (|dPSI|>0.2,
present in the 75p matrix) from bulk RNA counts, via ≤5-gene per-event ElasticNet/Ridge.

## Full-cohort result (18,168 events, 437 training samples)
Honest nested-OOF median Spearman **0.646**; **16,663** imputable (Sp>0.3), **13,497**
strong (Sp>0.5). All events got a valid ≤5-gene model (ridge 13,862 / elasticnet 4,306);
none were structurally unimputable. In-sample (biased) median 0.687.

## Feature-selection leakage was measured and removed from the reported metric
Candidate genes are chosen by correlation with the target. Scoring a model with
cross-validation **after** selecting features on all samples leaks the held-out folds into
selection and inflates performance. On a 300-event random subset:

| scoring | median Spearman | imputable Sp>0.3 | strong Sp>0.5 |
|---|---|---|---|
| selection-outside-CV (biased) | 0.677 | 291 | 253 |
| **nested OOF (honest, reported)** | **0.630** | 281 | 231 |

Median per-event inflation 0.043 (90th pct 0.105). The shipped per-event reliability
(`cv_spearman`) is the **nested out-of-fold** value: for each outer fold, masked-Pearson
feature selection is recomputed on that fold's training samples only, a RidgeCV is fit on
≤5 of them, and the held-out samples are predicted. `cv_spearman_insample` retains the
biased score for transparency.

## Why performance is high (not circular)
The targets are events pre-selected for strong differential splicing between AML
mutation/covariate subgroups; those subgroups carry distinct expression programs, so a few
protein-coding genes genuinely track PSI out-of-fold. Selection is unrestricted (any
protein-coding gene, not just splicing factors or the host gene). This is honest held-out
prediction on cells/samples never used to pick that event's features, not label reuse.

## Train/inference consistency
Features are normalized identically in training and inference: cp10k+log1p on the **all-gene**
library size, then z-score with the stored training μ/σ. Training starts from the same
**counts** matrix a query arrives as, so a new RNA-seq sample (raw counts, or cellHarmony-lite
scaled with `--normalized`) is scaled the same way. Genes absent from a query fall back to the
training mean (z=0). Predictions are clipped to the PSI range `[0,1]`; unimputable events
(too few non-NaN samples, or no usable model) are returned blank/`NaN`.

## Plumbing check (inference round-trip)
Imputing PSI back from the Leucegene counts with a 300-event bundle: 454 samples × 300
events, 100% finite, all values in [0,1], in-sample predicted-vs-observed median Pearson
0.736 (resubstitution — the honest out-of-fold number is 0.630).

## Reproduce
```bash
cd altanalyze3
python -m components.rna2psi.train_leucegene --out-dir <dir>   # full 18,168-event bundle
python -m components.rna2psi.train_leucegene --out-dir <dir> --max-events 300   # smoke
```
