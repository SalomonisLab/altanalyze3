# rna2psi (AML)

Impute **AltAnalyze2 MultiPath-PSI splicing values from bulk RNA gene counts**, and
reclassify AML driver signatures from the imputed PSI. On the same plumbing as
`rna2metabolite` / `rna2grn` / `rna2adt`; feature genes are matched to the input by
**Ensembl gene ID**.

## Default model — v2 (2026-07-09, Leucegene-focused)
- **Imputer:** one per-event **bootstrap-LASSO stability-selected Ridge** (≤5 predictor
  genes/event, selected from **all ~39,771 genes**), over **7,630 events** (union of the
  top-500-by-p-value events with `|dPSI|>0.15` from each of 23 AltAnalyze2 MultiPath-PSI
  signatures having >100 compliant events). Trained on **all 437 Leucegene AML samples**.
- **Reclassifiers:** per-mutation one-vs-rest sign-weighted marker scores (40 markers) at a
  98%-specificity operating point.
- **Full method** (mutations, thresholds, features, evaluation, performance): **`METHODS.md`**.

### Held-out assessment (BEAT-AML, cross-cohort, 98% specificity)
| mutation | N real | TP | FP | FN | sens | spec | prec |
|---|---|---|---|---|---|---|---|
| SF3B1 | 21 | 12 | 13 | 9 | 0.571 | 0.980 | 0.480 |
| ZRSR2 | 8 | 4 | 13 | 4 | 0.500 | 0.980 | 0.235 |
| SRSF2 | 41 | 16 | 12 | 25 | 0.390 | 0.981 | 0.571 |
| NPM1 | 131 | 33 | 10 | 98 | 0.252 | 0.981 | 0.767 |
| IDH2-R140Q | 60 | 8 | 12 | 52 | 0.133 | 0.980 | 0.400 |
| FLT3-ITD | 135 | 16 | 10 | 119 | 0.119 | 0.981 | 0.615 |
| TP53 | 45 | 2 | 12 | 43 | 0.044 | 0.981 | 0.143 |
| U2AF1-S34 | 17 | 0 | 13 | 17 | 0.000 | 0.980 | 0.000 |

Trained on all Leucegene; every count is on BEAT-AML (671 clinically-annotated samples).

## Usage
```python
from components.rna2psi import load_bundle
b = load_bundle()                                          # default = v2 focused bundle
res = b.predict_from_file("counts.txt", gene_axis="rows")  # genes(ENSG) x samples counts
res.predictions                                            # samples x 7,630 events (PSI in [0,1])
```
```bash
python -m components.rna2psi.cli predict --input counts.txt --output imputed_psi.txt
```
Normalization (cp10k+log1p on all-gene library size, then z-score with training μ/σ) is
applied internally; genes absent from the query take the training mean.

## Imputed PSI produced by the default model
- Leucegene: `~/Claude_Project/rna2psi_Leucegene/deploy/imputed_PSI.Leucegene.txt` (7,630 × 437)
- BEAT-AML:  `~/Claude_Project/rna2psi_Leucegene/deploy/imputed_PSI.BEATAML.txt` (7,630 × 671)
- Reclassification performance: `.../deploy/BEATAML_reclassification_performance.csv`

## Files
- `api.py`, `cli.py`, `_impute.py` — inference engine
- `_train.py`, `train_leucegene.py` — training engines
- `artifacts/rna2psi_leucegene_focused.pkl.gz` — **default imputer bundle (v2)**
- `artifacts/rna2psi_reclassifiers.pkl.gz` — per-mutation reclassifiers (v2)
- `METHODS.md` — precise methods & performance
- `legacy/v1_2026-07-07_protein_coding_marginal/` — prior v1 model (protein-coding,
  ElasticNet/Ridge, 18,168 events, honest OOF median Spearman 0.646) + its README/VALIDATION
