# rna2psi (AML)

Impute **splicing-event PSI from bulk/pseudobulk RNA counts** for AML, on the same
plumbing as `rna2metabolite` / `rna2grn` / `rna2adt`. One independent, sparse per-event
model; feature genes are matched to the input by **Ensembl gene ID** (the native key of
the Leucegene Kallisto counts).

## Final model
- One **per-event ElasticNet-or-Ridge** model, **CV-picked per event**, with **≤5
  predictor genes/event** (no restriction to splicing factors or host genes — any
  protein-coding gene may be selected).
- **Target:** raw PSI in `[0,1]` (no logit/arcsine transform). PSI missing values are kept
  as `NaN` and dropped per event for correlation and fitting. Predictions are clipped to
  `[0,1]`; unimputable events are returned as blank/`NaN`.
- **Features:** protein-coding Ensembl genes (`Hs_Ensembl_transcript-biotypes.txt` ==
  `protein_coding`), from the Leucegene Kallisto **counts** matrix, normalized exactly as
  at inference (**cp10k+log1p** on the all-gene library size, then z-score with training
  μ/σ). Starting from counts (not TPM) makes training match query datasets that arrive as
  counts / cellHarmony-lite scaled values.
- **Feature selection:** per event, candidate genes ranked by **missing-value-aware
  Pearson correlation** (pairwise-complete, computed in vectorized/matrix "vector mode"
  over the event's non-NaN samples); ElasticNet then sparsifies to the ≤5 largest-|coef|
  genes and the estimator is CV-picked (ElasticNet vs Ridge).

## Targets (the 18,168-event catalog)
Splicing events that (a) appear in a Leucegene variant/covariate `PSI.*` file
(`.../MajorCovariates/Events-dPSI_0.1_adjp` and `.../Variants/Events-dPSI_0.1_adjp`, parent
dir only) with **|dPSI| > 0.2**, and (b) are present in the 75p PSI EventAnnotation matrix.
18,275 unique candidate UIDs → **18,168 trainable** (present in the matrix); the 107 not in
the 75p matrix have no PSI ground truth and are not modeled.

## Held-out performance (leakage-free)
Reported per-event reliability (`cv_spearman` / `cv_r2` in the performance table) is a
**nested out-of-fold** estimate: feature selection is redone inside each training fold, so
it is free of feature-selection leakage. The in-pipeline selection-outside-CV score
(`cv_spearman_insample`) is kept for transparency and runs ~0.04 (median) higher.

| subset | median honest Spearman | imputable (Sp>0.3) | strong (Sp>0.5) | estimator mix |
|---|---|---|---|---|
| 300-event validation | 0.630 | 281 | 231 | ridge 232 / EN 68 |
| **full 18,168** | **0.646** | **16,663** | **13,497** | ridge 13,862 / EN 4,306 |

All 18,168 events got a valid gene model (≤5 features; 18,144 use the full 5). 0 were
structurally unimputable. The recommended deliverable is the imputable subset (Sp>0.3);
filter events on `cv_spearman`. In-sample (selection-outside-CV) median is 0.687, i.e. the
honest number is ~0.04 lower.

## Leucegene run outputs (`~/Claude_Project/rna2psi_Leucegene/`)
- `imputed_PSI.txt` — imputed PSI, events × samples (18,168 × 454)
- `imputed_PSI.imputable_Sp0.3.txt` — recommended subset (16,663 × 454)
- `rna2psi_per_event_performance.csv` — per-event reliability (honest + in-sample, features)

## Usage
```python
from components.rna2psi import load_bundle
b = load_bundle()                                        # default bundle
res = b.predict_from_file("counts.txt", gene_axis="rows")  # genes(ENSG) x samples counts
res.predictions                                          # samples x events (PSI in [0,1])
```
```bash
python -m components.rna2psi.cli model-info
python -m components.rna2psi.cli predict --input counts.txt --output imputed_psi.txt
```
Input genes match the model by **Ensembl gene ID**; absent genes take the training mean
(z=0). Normalization (cp10k+log1p on the all-gene library size, then z-score) is applied
internally — pass `--normalized` if the input is already cp10k+log1p / cellHarmony-lite
scaled. Output is written **events × samples** (UID rows) to mirror the AltAnalyze PSI
EventAnnotation layout, blank where an event is unimputable.

## Files
- `api.py`, `cli.py`, `_impute.py` — inference engine (no training code)
- `_train.py` — training engine (masked-Pearson selection, per-event fit, honest OOF)
- `train_leucegene.py` — reproduces the Leucegene bundle
- `artifacts/rna2psi_leucegene_bundle.pkl.gz` — final bundle (18,168 events × protein-coding genes)
