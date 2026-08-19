# rna2grn — Validation report

Imputation of GRN TF→target **edge activity** from RNA, for **identifying TFs whose
regulatory activity is perturbed** in a new sample — including the case where a TF's
activity changes (e.g. a mutation) while its own mRNA is unchanged, so the evidence
lives only in the **regulon's** expression.

All numbers below are produced by the canonical evaluation script and the `dev/`
scripts listed in **§9 Reproduce**. Output tables live in
`…/Human-GRN/July-2026-simple/rna2grn/benchmark/` and figures in `…/rna2grn/reports/`.

---

## 1. Is this ElasticNet? (method + modifications)

**Yes — it is a per-edge elastic-net regression, with the L1 mixing parameter set to
zero (`l1_ratio = 0`, i.e. an L2/ridge penalty).** It is *not* the original *single,
global, multi-task* elastic net — the change is **structural**, not the regression
family. Code: [`model.py`](model.py) (`RegulonEdgeRegressor`),
[`training.py`](training.py) (`build_bundle`), re-implemented inline for evaluation in
[`evaluate.py`](evaluate.py) (`_Edges.fit/predict`).

| | original (rna2lipid/rna2adt) | rna2grn (this module) |
|---|---|---|
| structure | **single, global, multi-task** | **one regression per edge** (7,486 independent models) |
| features | all selected genes → all targets jointly | **3 biologically-fixed covariates per edge**: target-gene expr, TF expr, TF **regulon mean** |
| penalty | elastic net (L1+L2) | **elastic net with `l1_ratio = 0` (L2/ridge)** |
| purpose | reconstruct a profile | **extrapolate**: regulon expression shift → edge-activity shift, beyond any reference |

Empirical justification (leave-one-donor-out R², 7,486 edges; `evaluate.py` model core):

- the **single global multi-task** elastic net (original structure): **R² strongly
  negative** — unstable extrapolation on held-out donors (the failure motivating the
  per-edge redesign).
- per-edge elastic net at **`l1_ratio = 0.5`: R² = 0.771** vs **`l1_ratio = 0` (ridge):
  R² = 0.773** — with only 3 covariates the L1 term is immaterial; `l1_ratio = 0` admits
  a closed-form vectorized solution (≈0.7 ms/edge vs ≈46 ms/edge for coordinate descent),
  so it is used. Setting `l1_ratio > 0` is a one-line change with no accuracy gain.

The decisive property is **structural** (per-edge + regulon features), not the penalty:
each edge's predicted activity moves with its target/regulon expression, so a perturbed
regulon raises the responsible TF's edges even when the TF's mRNA is flat. Learned
standardized coefficients: target-expr **0.017** > TF-expr **0.014** > regulon-mean
**0.006**.

---

## 2. Data and control association

Code: [`grn_io.py`](grn_io.py) (parsing/matching), [`dev/build_dataset.py`](dev/build_dataset.py)
(materializes `matched/dataset.npz`).

- **GRN** `tf_to_gene_connection_scores_…_COMBINED_AML_Multiome_RC2.txt`: 7,486 edges ×
  358 sample×cell-state columns. **RNA** `pseudobulk_counts_hashed.h5ad` (12,255
  pseudobulks keyed `Sample|CellState`).
- 354/358 GRN columns matched to an RNA pseudobulk; 4 dropped (`AML_13` — no RNA partner).
- **Control association (the key point).** There is **no multiome/TEA RNA** in the
  pseudobulk file. The two control GRN aggregates are paired, on the RNA side, with the
  **summed counts of the 12 CCHMC healthy-control donors** (`Annotation==Control &
  Dataset==CCHMC`), per cell state, then CP10k+log1p:
  `AM72, BF21, BF32, BF71, BM27, BM32, WF26, WF32, WF83, WM22, WM34, WM72`
  (10x 3′ v3, CD34/CD271/TNC-sorted BM). This **same** control RNA is connected to
  **both** control GRNs:
  - `Multiome_WT` GRN — 82 cell states (multiome-derived);
  - `RC2_TEA` GRN — 10 cell states (TEA-seq-derived): ERP-1, HSC-2, Intermediate Mono-1,
    LMPP-1, LMPP-1-cycling, LMPP-2, MEP-2, MPP-1, MPP-2, MPP-MEP.
  (Verified: the CCHMC control-RNA vector is identical across the Multiome and TEA rows
  for all 10 shared cell states.)
- Final `dataset.npz`: 341 pseudobulks (223 patient + 26 AML-CCHMC + 82 Multiome + 10
  TEA), 30 samples, 83 cell states. **100 % of 217 TFs and 2,361 target genes present in
  the RNA.**

### Global magnitude offset (must read before interpreting differentials)

Total GRN activity (Σ of 7,486 edges) per pseudobulk, **median**:

| group | assay | median Σ-activity |
|---|---|---|
| AML patient / AML-CCHMC | 3′ | 291 / 327 |
| **TEA control** | TEA-seq | **67** |
| **Multiome control** | multiome | **46** |

Cell-state-matched **AML/Multiome ratio = 4.9×**. So raw "AML−control" differentials carry
a large **assay-scale offset** (≈60 % of TFs read "up in AML" in every pseudobulk). **TEA
(67) is the closer magnitude match to AML than Multiome (46)** — confirming the lab's
preference — but covers only 10 cell states. We therefore report **both** controls and
both a **raw** and an **offset-removed (per-pseudobulk-centered)** differential.

---

## 3. Validation design

Code: [`evaluate.py`](evaluate.py) (canonical), [`dev/validate_directional.py`](dev/validate_directional.py)
(single-control variant), [`dev/perturbation_test.py`](dev/perturbation_test.py) (Test B).

- **sample = donor; pseudobulk = donor × cell-state** (the prediction unit). Hold-out is
  at the **donor** level — every pseudobulk of a held-out donor goes to test together, so
  no cell-state of a held-out donor leaks into training.
- **AML-7 excluded** from all evaluation (reserved as the raw-`.h5` real-world speed input;
  load 0.06 s, whole-sample GRN 0.32 s).
- The merged control aggregates **cannot be held out** (they are unique, merged) — they
  are **retained in training** every fold and serve as the differential baseline.
- **Grouped 5-fold** over the remaining **27 AML donors**:

| fold | held-out AML donors | retained AML donors | held-out pseudobulks |
|---|---|---|---|
| 0 | 6 | 21 | 28 |
| 1 | 6 | 21 | 67 |
| 2 | 5 | 22 | 39 |
| 3 | 5 | 22 | 35 |
| 4 | 5 | 22 | 45 |

Every donor held out once → **214 held-out pseudobulks across 40 cell states**;
**61,845** (pseudobulk, TF, control) examples scored.
- **Metric:** per held-out pseudobulk (S,c), per-TF activity = mean of that TF's imputed
  edge scores; `true_diff = TFact(S,c) − TFact(control,c)`, `imp_diff = imputed AML −
  imputed control`. Among **truly-differential TFs** (top |true_diff| strata), report the
  fraction with **correct predicted direction** (`sign(imp_diff)==sign(true_diff)`).

---

## 4. Held-out AML samples (file `benchmark/detail_per_sample.csv`)

27 donors; **Inv16 ×17, PHIP ×6, generic-AML ×3, plus AML-12/AML-14 CCHMC subtypes**
(AML-7 / SRSF2 excluded). Per-donor directional accuracy (top-25 % differential TFs,
TEA-preferred control) ranges **0.72–0.91, median 0.83**. Selected rows:

| donor | mutation | n cell-states | n diff-TF | dir-acc |
|---|---|---|---|---|
| 562_AfInv16_25F | Inv16 | 6 | 455 | 0.910 |
| 1009_AfInv16_29M | Inv16 | 6 | 434 | 0.901 |
| 1906_Aft82_P36M | PHIP | 6 | 312 | 0.897 |
| 4440_AfInv16_P24M | Inv16 | 27 | 1483 | 0.834 |
| AML-14_CCHMC | (CCHMC) | 6 | 379 | 0.757 |
| AML-12_CCHMC | (CCHMC) | 7 | 441 | 0.741 |
| 550_EInv16_38M | Inv16 | 13 | 758 | 0.719 |
| 2853_Aft36_P71F | PHIP | 10 | 260 | 0.715 |

---

## 5. Directional accuracy — dual control (file `benchmark/directional_metrics_dual.json`)

Correct-direction fraction among truly-differential TFs:

| control | metric | top 50 % | top 25 % | top 10 % |
|---|---|---|---|---|
| **TEA** (better match, 71 pb) | raw | 0.781 | 0.845 | **0.898** |
| **TEA** | **centered (TF-specific)** | **0.910** | 0.904 | 0.892 |
| Multiome (214 pb) | raw | 0.749 | 0.822 | 0.874 |
| Multiome | centered (TF-specific) | 0.907 | 0.901 | 0.848 |

Reading: (i) **TEA gives higher raw accuracy than Multiome** — its smaller assay offset is
the cleaner contrast, as the lab predicted. (ii) **Removing the per-pseudobulk offset lifts
both to ~0.90** — the raw differential was *diluted* by the global offset, and the
underlying **TF-specific** direction is recovered ~90 % of the time. (iii) Accuracy rises
with |true Δ|. All p ≈ 0 vs chance (binomial); up/down balanced (`dev/validate_directional.py`).

---

## 6. Per-cell-type variability (file `benchmark/detail_per_cellstate.csv`)

40 cell states; per-cell-state directional accuracy (top-25 % diff TFs) **range 0.72–0.90**:

- **Best:** LMPP-2 0.90, MPP-MEP 0.90 (n=23 donors), cDC1 0.89, MDP-2 0.88, Mono-1 0.87,
  Myeloid intermediate-1 0.85, cDC2-2 0.85, preNeu 0.85.
- **Lower:** rarer / lymphoid-leaning states and `MultiLin-GMP-2`, `pDC`, `Intermediate
  Mono-2` (~0.68–0.72 in the single-control run).
- High-coverage myeloid/progenitor states (MPP-MEP n=23, Intermediate Mono-1 n=22) are the
  most stable.

---

## 7. Per-TF performance (file `benchmark/detail_per_tf.csv`)

All **217 TFs scored** in every pseudobulk; 106 are frequently differential (≥20 examples).

- **Best-inferred (dir-acc = 1.0):** E2F3, USF1, TCF12, STAT1, MXD1, RXRA, RELB, ARNT,
  MAFK, ATF3, FOSB. Myeloid/AML-relevant frequently-differential TFs also inferred well:
  **RUNX1, JUN, MAFG, SPI1, JUNB, FOS**.
- **Worst-inferred (dir-acc ≤ 0.45):** USF2 0.29, KLF3 0.33, RARA 0.33, JDP2 0.34, GMEB1
  0.43, CXXC5 0.45 — mostly ubiquitous/housekeeping TFs with small differentials.

---

## 8. Top / worst inferred TFs per pseudobulk (file `benchmark/per_pseudobulk_topworst_tf.csv`)

For each held-out (donor, cell-state), the most-differential TFs ranked by whether the
imputed activity matched the measured direction (TEA control where covered, else Multiome).
214 tables (71 TEA, 143 Multiome). Worked examples:

**`1009_AfInv16_29M | Intermediate Mono-1` (control = TEA-seq):**

| | TF | measured Δ | imputed Δ | correct dir |
|---|---|---|---|---|
| top | JUND | +0.297 | +0.016 | ✓ |
| top | FOSL1 | +0.200 | +0.148 | ✓ |
| top | SREBF2 | +0.127 | +0.122 | ✓ |
| top | NFE2L2 | +0.119 | +0.065 | ✓ |
| top | STAT3 | +0.110 | +0.047 | ✓ |
| worst | CEBPD | +0.215 | −0.002 | ✗ |
| worst | SPI1 | +0.214 | −0.018 | ✗ |
| worst | CEBPB | +0.087 | −0.040 | ✗ |

**`1009_AfInv16_29M | Mono-1` (control = Multiome):** top correct = JUND, FOSL1, SREBF2,
SPI1, ELF1; worst = KLF3 (measured +0.13, imputed −0.01).

Two consistent patterns the team should note: (i) the model recovers **AP-1 / stress TFs
(JUN/JUND/FOS/FOSL1)** and many progenitor TFs reliably; (ii) it can **miss the
direction of CEBP-family and SPI1** in specific mature-monocyte pseudobulks; (iii)
imputed Δ magnitudes are **shrunk** relative to measured (ridge regularization + the
assay-offset term on the measured-control side) — **direction is recovered far better
than magnitude**, so rank/sign, not absolute Δ, is the quantity to trust.

---

## 9. Decisive out-of-distribution test (the actual goal)

Test A (above) recovers an AML signature that is partly **shared with the training AML
donors**, so a retrieval model can game it. The goal — detect a **novel** perturbation —
is tested by an **in-silico regulon perturbation** ([`dev/perturbation_test.py`](dev/perturbation_test.py)):
raise one TF's regulon in a healthy pseudobulk (TF mRNA fixed), rank the perturbed TF among
217.

| model | median rank | top-1 | top-5 |
|---|---|---|---|
| **per-edge regulon regression (shipped)** | **2** | 37 % | **83 %** |
| kNN retrieval | 74 | 1 % | 5 % |

Only the per-edge regression localizes the perturbed TF; retrieval cannot (it returns the
nearest existing reference). This is why kNN is **not** used despite scoring well on
in-distribution Test A.

---

## 9b. Lung reference (LungMAP IPF/control)

Scope: this section covers ONLY the `lung` bundle. Sections 1–9 and 10–11 describe the
`leukemia` (default) bundle. The multi-reference change leaves that bundle **bit-identical** (6 × 7,486 = 44,916 values compared, 0 differences,
across `load()`, `reference="leukemia"`, an explicit path, and `load_bundle()`).

Code: [`grn_io.py`](grn_io.py) (`mangle_cellstate`, `match_lung_columns`),
[`dev/build_lung_reference.py`](dev/build_lung_reference.py),
[`dev/validate_lung_reference.py`](dev/validate_lung_reference.py). Report:
`/Users/saljh8/Dropbox/LungMAP/GRN/rna2grn/validation/validation_report.json` (status PASS).

**Inputs** (read-only, neither modified):
`/Users/saljh8/Dropbox/LungMAP/GRN/TF_to_Gene_connection_scores_log10-NOT_ordered_clusters_ALL_GENES-threshold-1.txt`
(57,307 edges × 39 columns, mtime 2026-08-18 21:45:43) and
`/Users/saljh8/Dropbox/LungMAP/code/ILD/AltAnalyze-create-cH-reference/ExpressionInput/exp.pseudobulks-IPF-control.txt`
(33,538 genes × 85 pseudobulks, mtime 2025-11-06 18:49:26).

**Matching.** GRN columns are `<Group>.<CellState>`; pseudobulks are `<CellState>|<Group>`.
The GRN writer drops `+` then maps ` ` and `/` to `_`; `mangle_cellstate` reproduces that
exactly, so matching is string equality, **not** fuzzy suffix search, and a
non-injective mangling raises. **39/39 GRN columns matched (100 %)**; 0 unmatched;
46/85 pseudobulks have no GRN column and are not training rows. Coverage: **221/221 TFs**,
**4,984/4,984 targets**, **57,307/57,307 edges (100 %)**, 5,095 feature genes, 27 cell
states, 14/39 Control + 25/39 IPF rows.

**Normalization.** `pseudobulk_only.py:37` averages `log1p(CP10k)` over cells (layer built
by `normalize.py:79-88`), so the lung X is a **mean-of-log**, whereas the leukemia X is
**log-of-summed-counts**. The builder writes the matrix to `dataset.npz` unchanged, records the statistic as
`pseudobulk_statistic`, and `_warn_if_counts_path` raises a `RuntimeWarning` when a caller
feeds counts to this reference. Nothing is silently converted.

**Fit.** `training.build_bundle` unchanged, `ridge_lambda = 1.0`, **all 39 pseudobulks,
no holdout** (requested). Bundle 3.93 MB.

**Checks (all PASS).**

| check | result |
|---|---|
| A edge order == GRN file row order (bundle and dataset), metadata alignment, feature space, finiteness, non-negativity | 11/11 PASS |
| B1 public API on the dataset matrix == `model.predict` | max abs diff **0.000e+00** (tol 1e-9) |
| B2 public API re-reading the text file at float64 | max abs diff **1.76e-07** (tol 1e-5; float32 storage, `max\|X_text − X_npz\| = 2.2e-07`) |
| D core cross-check: `evaluate._Edges` reproduces the shipped fit | max abs diff **4.4e-16** (tol 1e-9) |
| E `RuntimeWarning` on the counts path | fires |

**Accuracy** (39 pseudobulks × 57,307 edges). Column C is the requested configuration and
is optimistic by construction. Column D is an **addition** that changed no shipped file.

| metric | C resubstitution | D leave-one-pseudobulk-out |
|---|---|---|
| R² vs global mean | 0.887 | 0.839 |
| profile r median / min | 0.940 / 0.872 | 0.907 / 0.792 |
| per-edge r median | 0.941 | 0.916 |
| mean absolute error | 0.084 | 0.097 |
| predicted / true mean | 0.281 / 0.267 | — |

**Deployment in scALABLE (cellHarmony-web).** All four human lung references expose the
`grn` modality and point at this bundle
(`altanalyze3/components/cellHarmony/flask/reference_config.json`). The pipeline builds
the query pseudobulk with the statistic the bundle declares, so the lung path averages
`query_adata.X`, which `cellHarmony_lite.py:554-556` already normalized to CP10k+log1p per
cell. Measured on 12,000 real lung cells from
`/Users/saljh8/Dropbox/LungMAP/code/ILD/output_file_with_umap.h5ad` (43 cell states, 90
samples, 1,563 pseudobulks): the builder's pseudobulk equals the training formula
(`pseudobulk_only.py:37`) to **max abs diff 1.9e-06** (float32 h5ad storage). The marrow
GRN path stays on `sum_counts` and is **bit-identical** across 8,983,200 values. Report:
`/Users/saljh8/Dropbox/LungMAP/GRN/rna2grn/validation/impute_modalities_report.json`
(status PASS), script
[`../cellHarmony/flask/dev/validate_impute_modalities.py`](../cellHarmony/flask/dev/validate_impute_modalities.py).

**Training-column provenance (partly unverified).** I reproduced the file column
`AT1|Control` two ways from
`/Users/saljh8/Dropbox/LungMAP/code/ILD/output_file_with_umap.h5ad`: a mean over
per-sample pseudobulks reaches **r = 0.99850**, a pooled cell mean reaches **r = 0.99710**.
Neither matches exactly. The residual most likely comes from the curated sample roster
(`pseudobulk_metadata_manual_annotation_1_by_Sample_Name-ILD53-removed.txt`). I have not
verified that roster. I have settled the statistic family. I have not settled the exact sample list.

**Limitations specific to the lung reference.**

1. **39 training rows** fit 4 parameters per edge. The leukemia reference has 341.
2. **No donor replication.** Each pseudobulk pools every donor in its group, so a
   Control-vs-IPF difference is a magnitude, not a significance test.
3. **46/85 pseudobulks have no measured GRN.** No training row matches those 46 profiles,
   and this report can check none of them against a truth value.
4. **Profile r carries the same high structural floor** flagged in §7 for leukemia; the
   0.94 resubstitution figure is not evidence of sample-specific accuracy.
5. **No perturbation test.** Nobody has run the §9 in-silico regulon perturbation test or
   the §8 differential-TF test on lung. I have not verified that this reference
   localizes a perturbed TF.
6. **No cross-reference test.** Do not predict lung pseudobulks with the leukemia bundle,
   or the reverse. The two edge sets and the two input statistics differ.

---

## 10. Limitations

1. **One protocol-mismatched control class.** AML (3′) vs Multiome/TEA controls → a ~5×
   assay-scale offset; use the **centered** differential and **TEA-matched** controls where
   possible. A clean disease-only differential would need a healthy **3′** GRN, absent here.
2. **In-distribution Test A is gameable** by retrieving cousin AML donors; the in-silico
   perturbation (§9) is the test that matches the use case.
3. **Magnitude is shrunk**; trust direction/rank, not absolute Δ.
4. Single-replicate pseudobulks → "differential" is a magnitude threshold, not a
   significance test (hence magnitude strata).
5. AML-12/AML-14 CCHMC subtype labels here come from a name-merged annotation lookup;
   the pseudobulk-level annotations (e.g. DNMT3A, CSF3R) are authoritative.

---

## 11. Reproduce (code pointers)

| step | code | output |
|---|---|---|
| parse GRN + match to RNA | [`grn_io.py`](grn_io.py) | — |
| build matched dataset | [`dev/build_dataset.py`](dev/build_dataset.py) | `matched/dataset.npz` |
| signal / variance diagnostics | [`dev/diagnostics.py`](dev/diagnostics.py) | console |
| model fit + bundle | [`training.py`](training.py), [`model.py`](model.py) | `rna2grn_bundle.pkl.gz` |
| **canonical evaluation (dual control + per-pseudobulk TF tables + detail)** | [`evaluate.py`](evaluate.py) | `benchmark/directional_dual_control.csv`, `directional_metrics_dual.json`, `per_pseudobulk_topworst_tf.csv`, `detail_per_sample.csv`, `detail_per_cellstate.csv`, `detail_per_tf.csv` |
| single-control directional + baselines/centering | [`dev/validate_directional.py`](dev/validate_directional.py) | `benchmark/directional_*.csv` |
| in-silico perturbation (Test B) | [`dev/perturbation_test.py`](dev/perturbation_test.py) | `benchmark/perturbation_tf_ranks.csv` |
| protocol similarity (Multiome vs TEA) | [`dev/protocol_similarity.py`](dev/protocol_similarity.py) | `reports/protocol_similarity_*.csv` |
| figures | [`dev/make_figures.py`](dev/make_figures.py) | `reports/fig1–5*.pdf` |

Run the canonical evaluation:
```bash
python -m altanalyze3.components.rna2grn.evaluate \
    --dataset .../rna2grn/matched/dataset.npz --out .../rna2grn/benchmark
```
Interpreter: `/opt/homebrew/opt/python@3.11/bin/python3.11` (numpy/scipy/sklearn only).
