# rna2psi — Methods (default model v2, 2026-07-09)

Imputation of AltAnalyze2 **MultiPath-PSI** splicing values from bulk RNA gene counts, and
reclassification of AML driver signatures from the imputed PSI. Trained on Leucegene,
assessed on BEAT-AML.

## 1. Splicing quantification (input targets)
- **PSI**: AltAnalyze2 **MultiPath-PSI** event annotation, Leucegene hg38, 75%-detected
  matrix `Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation-75p.txt` (96,875 events × 437
  samples; missing values = `NA`/blank, kept as NaN).
- Event UID = AltAnalyze MultiPath-PSI `UID` (`gene:ENSG:Ex-Ey|ENSG:Ex-Ez`).

## 2. Event set — which signatures & thresholds
Per-signature differential-splicing files in
`MajorCovariates/Events-dPSI_0.1_adjp/focused/` (41 AltAnalyze2 MultiPath-PSI
mutation/fusion/subtype comparisons vs Others). Selection:
- **Compliant event**: `|dPSI| > 0.15`.
- **Signature kept** only if it has **> 100 compliant events** → **23 of 41** signatures.
- Per kept signature: **top 500 events by p-value** (adjusted `adjp`, then `rawp`), among
  the compliant events.
- **Union ∩ 75p matrix = 7,630 events** (the imputed feature space).

**23 signatures used** (AltAnalyze2 MultiPath-PSI comparisons; compliant-event count):
SNRNP40 (5152), U2AF1-S34 (3578), SRSF2-P95 (1788), RUNX1-fusions (1365), SRSF2-8AA (1251),
TP53 (1223), U2AF2 (1039), SETX (962), ZRSR2 (774), THOC5 (752), SF3B1 (717), U2AF1-Q157
(577), PML-RARA-fusions (482), NPM1-FLT3-IDH2-TET2 (459), MLL-fusions (448), CBFB-MYH11
(431), EVI1-fusions (357), HNRNPK-KD (268), IDH2-R140Q (251), NPM1 (203), FLT3-ITD (159),
FLT3-NPM1 (126), CEBPA-allelic (121).

## 3. Imputation model (per event) — features & estimator
Expression = `counts.Leucegene.txt`, **all 39,771 Ensembl genes** (no biotype filter).
- **Normalization**: cp10k + log1p (library size over all genes), then z-score with training
  μ/σ. Genes absent from a query take z = 0 (training mean).
- **Candidate genes**: top **k = 20** by missing-value-aware (pairwise-complete) Pearson
  correlation of the event's PSI against each gene, over that event's non-NaN samples.
- **Feature selection**: **bootstrap-LASSO stability selection** — 10 bootstraps of
  `Lasso(alpha=0.01)`; keep genes selected in **≥ 50%** (≥5/10); cap at **≤ 5 genes/event**.
- **Estimator**: `RidgeCV` (α ∈ {0.1, 1, 10, 100}) on the selected genes.
- **Prediction**: linear, clipped to `[0, 1]`. One independent model per event (no k-NN).
- **Training set**: **all 437 Leucegene samples** (no samples excluded).

Bundle (`artifacts/rna2psi_leucegene_focused.pkl.gz`, rna2psi format): `X_columns` (ENSG),
`Y_columns` (7,630 UIDs), `mu`, `sd`, per-event `sel_idx`/`coef`/`intercept`, `metadata`.

## 4. Reclassification model (per mutation) — one-vs-rest
For a mutation signature, from the **imputed** Leucegene PSI (3-fold out-of-fold imputed, so
the classifier is domain-matched to genuine imputation, not resubstitution):
- **Marker selection**: `oncosplice.metadataAnalysis.limmaCompute` moderated-t (mutant vs
  rest); candidate = top **200** events by `|Δmean PSI|`; keep the top **40** by training
  AUROC (clean-margin selection).
- **Score**: mean over the 40 markers of `sign(Δmean) × z(imputed PSI)` — a sign-weighted
  marker score (chosen over SVM because it transfers across cohorts).
- **Operating point**: threshold set to **98% specificity**.

Leucegene mutation labels: `Metadata-Analysis/subtype_metadata.txt` (binary subtype
membership). Reclassifier bundle: `artifacts/rna2psi_reclassifiers.pkl.gz`.

## 5. Evaluation — dataset, labels, metric
- **Assessment cohort**: **BEAT-AML** (waves 1–4), `beataml_waves1to4_counts_dbgap.txt`
  (63,677 ENSG × 707 samples; 671 with clinical annotation). No Leucegene held-out — the
  model trains on all Leucegene and is assessed cross-cohort, as specified.
- **Ground-truth labels**: BEAT-AML clinical `variantSummary` (`beataml_wv1to4_clinical.xlsx`,
  matched by `dbgap_rnaseq_sample` = `BA…R`), one-vs-rest (positive = gene present in
  variantSummary).
- **Metric**: confusion matrix (TP/FP/TN/FN), sensitivity, specificity, precision at the
  **98% specificity** operating point (39,270/39,771 model genes matched in BEAT-AML).

## 6. Performance (BEAT-AML, 98% specificity)

| mutation | Leuc train | N real | N pred | TP | FP | TN | FN | sens | spec | prec |
|---|---|---|---|---|---|---|---|---|---|---|
| SF3B1 | 8 | 21 | 25 | 12 | 13 | 637 | 9 | 0.571 | 0.980 | 0.480 |
| ZRSR2 | 5 | 8 | 17 | 4 | 13 | 650 | 4 | 0.500 | 0.980 | 0.235 |
| SRSF2 | 35 | 41 | 28 | 16 | 12 | 618 | 25 | 0.390 | 0.981 | 0.571 |
| NPM1 | 35 | 131 | 43 | 33 | 10 | 530 | 98 | 0.252 | 0.981 | 0.767 |
| IDH2-R140Q | 49 | 60 | 20 | 8 | 12 | 599 | 52 | 0.133 | 0.980 | 0.400 |
| FLT3-ITD | 36 | 135 | 26 | 16 | 10 | 526 | 119 | 0.119 | 0.981 | 0.615 |
| TP53 | 25 | 45 | 14 | 2 | 12 | 614 | 43 | 0.044 | 0.981 | 0.143 |
| U2AF1-S34 | 10 | 17 | 13 | 0 | 13 | 641 | 17 | 0.000 | 0.980 | 0.000 |

Splicing-factor (SF3B1, ZRSR2, SRSF2) and NPM1 signatures transfer Leucegene→BEAT-AML at
98% specificity with usable precision; small-training-N (U2AF1-S34) or expression-invisible
(TP53) signatures do not. The imputer (7,630 events) is reusable for any downstream PSI use.

## Provenance
- Training: Leucegene AML (hg38 AltAnalyze2 MultiPath-PSI, all 437 samples).
- Assessment: BEAT-AML waves 1–4 (671 clinically annotated samples).
- Build script: `~/Claude_Project/rna2psi_Leucegene/build_focused_model.py`.
- Prior model (v1, protein-coding + ElasticNet/Ridge, 18,168 events): `legacy/`.
