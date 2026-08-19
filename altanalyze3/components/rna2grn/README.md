# rna2grn

`rna2grn` predicts **gene-regulatory-network (GRN) connection / activity scores** —
one value per TF→target-gene edge — from RNA expression, mirroring the `rna2adt`
and `rna2lipid` cross-modality imputation modules. Given a pseudobulk (a cell
state within a sample) it imputes the full 7,486-edge GRN that would otherwise
require an expensive multiome/TEA-seq GRN inference run.

## The goal: infer regulatory-edge activity, not retrieve a profile

The intended use is **perturbation detection**: give the model a new sample (e.g.
healthy-looking marrow carrying a mutation that changes one TF's *activity* but
not its mRNA level). The TF's regulon shifts at the expression level, and the
model must raise *that TF's* edge activity scores — **even beyond any profile in
the reference** — so the perturbed TF can be identified.

This rules out nearest-neighbor retrieval: a kNN model can only blend existing
reference profiles, so it can never express a perturbation that is absent from
the reference. The shipped model is therefore a **per-edge regulon regression**:

> For each edge (TF_i → gene_j), the activity score is a ridge regression on three
> standardized features of the pseudobulk — the **target gene's expression**, the
> **TF's expression**, and the **TF's regulon mean** (mean expression of all of
> TF_i's targets). Fit is vectorized across all 7,486 edges.

This is a supervised, learned form of regulon-activity inference (cf. VIPER /
AUCell / decoupleR), trained directly on the lab's GRN scores. It **extrapolates**:
when a regulon's targets rise, the predicted edge activities rise past any
reference. The learned coefficients rely most on **target expression** (mean |β|
0.017) > TF expression (0.014) > regulon mean (0.006) — exactly what lets it
respond when the TF's own mRNA is flat but its targets move.

## Decisive test: in-silico regulon perturbation

Raise one TF's regulon expression in a healthy pseudobulk (TF mRNA held flat),
predict, and ask where the perturbed TF ranks by increased predicted activity
(of 217 TFs, 1 = perfect):

| model | median rank | top-1 | top-5 |
|---|---|---|---|
| **per-edge regulon regression** | **2** | **37 %** | **83 %** |
| kNN retrieval | 74 | 1 % | 5 % |

kNN cannot localize the perturbed TF; the regression does. This is the capability
the module exists for.

## It is also a better absolute predictor (leave-one-sample-out)

| metric | per-edge regression | kNN |
|---|---|---|
| R² vs global mean | **0.76** | 0.63 |
| skill vs cell-state-mean baseline | **+0.31** | +0.23 |
| within-cell-state edge r | 0.42 | 0.42 |
| profile r | 0.97 | 0.94 |
| bundle size | **0.47 MB** | 8.5 MB |

So the regression dominates kNN on every axis and is far smaller (per-edge
coefficients only, no stored reference matrices).

## Validation on held-out AML vs control (the lab's spec)

For a held-out AML sample, the TF activities truly increased vs the Multiome_WT
control (from the real GRN) are recovered in the **imputed** GRN:

- Spearman(true, imputed) differential TF activity over 217 TFs: **median 0.51**
- precision@10 (truly top-increased TFs recovered): **median 0.50**

Well-represented signatures (Inv16) reach 0.6–0.68; the SRSF2 singleton (AML-7,
a splicing factor with subtle/indirect GRN effect and no cousin in the reference)
is the hardest at 0.26 — biologically sensible. For AML-7, the top imputed
differentially-active TFs are **SPI1 (PU.1)**, CEBPA/CEBPB/CEBPD, FOS, FLI1 — the
correct myeloid/AML regulators.

> Honest caveat: a held-out AML sample shares an "AML signature" with the *other*
> AML samples still in the reference, so this test does not isolate truly novel
> perturbations — the in-silico perturbation test does. Both are reported.

## What signal is being predicted

GRN variance is 65.5 % between cell states (RNA decodes identity) and **34.5 %
within cell state** (the hard, sample-specific rewiring). Profile correlation has
a 0.94 floor from static structure, so it is **not** the metric to optimize —
judge by skill over the cell-state-mean baseline, within-cell-state edge
correlation, and perturbation detection.

## Cross-protocol behavior (Multiome vs TEA-seq vs 3′)

Multiome and TEA-seq control GRNs for the **same cell state** correlate at Pearson
0.80 (vs 0.40 for *different* cell states) — cell identity dominates protocol. But
cross-protocol *absolute magnitudes* differ; use the model **within a protocol
family** (the AML 3′ CITE setting, matching the test data) for quantitative work.

## Data and sample matching

GRN `tf_to_gene_connection_scores_..._COMBINED_AML_Multiome_RC2.txt` (7,486 edges
× 358 sample×cell-state columns) matched to `pseudobulk_counts_hashed.h5ad`
(keyed `Sample|CellState`). 354/358 columns matched; 4 dropped (`AML_13`, no RNA
partner). Final reference: 341 pseudobulks (223 patient + 26 AML-CCHMC + 92
control aggregates), 30 samples, 83 cell states; **100 % of 217 TFs and 2,361
target genes present in the RNA**. Matching + sample-name overrides in `grn_io.py`.

## Normalization (reproducible on any new matrix)

Pseudobulk counts → **CP10k + log1p** → align to the 2,415 TF/target feature genes
→ per-edge feature standardization (stored) → regression. Reproducible from any
raw count matrix; for new data, sum cell counts per cell state → CP10k+log1p →
predict.

## Bundle

A gzip pickle (`rna2grn_bundle.pkl.gz`, ~0.5 MB): `model` (RegulonEdgeRegressor —
per-edge coefficients, standardization stats, sparse regulon-averaging matrix),
`scaler_x`/`scaler_y` = None, `X_columns` (feature genes), `Y_columns` (edge ids
`TF|Gene`), `metadata` (edge TF/Gene, TF list, training provenance).

## Usage

```python
from altanalyze3.components.rna2grn import Rna2GrnBundle
b = Rna2GrnBundle.load()

# per-cell-state GRN from an annotated h5ad
res = b.predict_from_h5ad("sample.h5ad", groupby="cell_state")
res.predictions                      # cell_state x 7486 edges

# IDENTIFY PERTURBED TFs: rank TFs more active in a query vs a control
diff = b.differential_tf_ranking(query_pred, control_pred)   # TF, delta_activity (sorted)

# per-TF activity matrix (groups x TFs)
tfa = b.tf_activity(res.predictions)

# from a raw 10x .h5 (annotate cells first for per-state GRNs)
res = b.predict_from_10x_h5("filtered_feature_bc_matrix.h5", groupby="cellharmony_label")
```

CLI (`python -m altanalyze3.components.rna2grn.cli`): `list-references`,
`model-info`, `predict-csv`, `predict-h5ad`, `predict-10x`.

## References (tissue-specific bundles)

The module ships one bundle per atlas. `leukemia` stays the default, so callers
that pass nothing keep their previous behaviour (verified bit-identical).

| reference | tissue | edges | TFs | training pseudobulks | bundle |
|---|---|---|---|---|---|
| `leukemia` (default) | human bone marrow / AML | 7,486 | 217 | 341 | `rna2grn_bundle.pkl.gz` (0.47 MB) |
| `lung` | human lung, LungMAP IPF vs control | 57,307 | 221 | 39 | `rna2grn_lung_bundle.pkl.gz` (3.93 MB) |

```python
b = Rna2GrnBundle.load(reference="lung")      # or load_bundle(reference="lung")
Rna2GrnBundle.load("/path/to/custom.pkl.gz")  # explicit path still works
```
```bash
python -m altanalyze3.components.rna2grn.cli list-references
python -m altanalyze3.components.rna2grn.cli model-info --reference lung
```
Passing `--bundle` and `--reference` together raises `ValueError` rather than
silently preferring one.

### The lung reference

Built from 39 measured GRN columns matched 39/39 to `<CellState>|<Group>` pseudobulks
(221/221 TFs and 4,984/4,984 targets present in the RNA, so 57,307/57,307 edges are
usable). **No internal validation** — all 39 pseudobulks train the model, by request.

Its pseudobulks are the **mean over cells of `log1p(CP10k)`**, not the leukemia
reference's `log1p(CP10k of summed counts)`. Same family, different statistic, so
pass `normalized=True` (or CLI `--normalized`); the counts path emits a
`RuntimeWarning` for this reference instead of silently shifting the input space.

Accuracy over 39 pseudobulks x 57,307 edges — resubstitution is the requested
configuration, leave-one-pseudobulk-out is an added honest counterpart:

| metric | resubstitution | leave-one-out |
|---|---|---|
| R² vs global mean | 0.887 | 0.839 |
| profile r (median / min) | 0.940 / 0.872 | 0.907 / 0.792 |
| per-edge r (median) | 0.941 | 0.916 |
| mean absolute error | 0.084 | 0.097 |

### Pseudobulk statistic (how inference must form its input)

A bundle declares in `metadata["pseudobulk_statistic"]` how its training pseudobulks
were formed, and `predict_from_adata` honours that declaration by default:

| statistic | meaning | bundle |
|---|---|---|
| `sum_counts` (default when a bundle declares nothing) | sum the group's counts, then CP10k + log1p | `leukemia` |
| `mean_over_cells_of_log1p_cp10k` | average the group's already CP10k+log1p rows | `lung` |

```python
b.predict_from_adata(adata, groupby="pb", layer=None)          # uses the bundle's own statistic
b.predict_from_adata(adata, groupby="pb", layer="counts",
                     pseudobulk_statistic="sum_counts")        # explicit override
```
Sending a count matrix down the mean-of-log path raises `ValueError` rather than
averaging counts as if they were log values.

Caveats specific to this reference: 39 rows fit 4 parameters per edge; each
pseudobulk pools all donors of a group, so there is no donor-level replication and a
"differential" is a magnitude, not a test; the model predicts 46/85 lung pseudobulks that have no measured
GRN, so no training row matches them. Full report:
`/Users/saljh8/Dropbox/LungMAP/GRN/rna2grn/README.md`.

## Caveats

- **Partial distillation.** GRN scores are derived from expression (and chromatin
  for multiome/TEA), so RNA→GRN is partly a fast surrogate for the GRN
  computation. The per-edge model makes this explicit and useful: it maps
  regulon-level expression to edge activity and extrapolates.
- **Profile correlation is misleading** (0.94 floor); use skill / perturbation
  detection.
- **Within-protocol** quantitative use only (see cross-protocol section).
- `KnnGrnRegressor` is retained for pure profile reconstruction but is **not**
  valid for perturbation / differential-activity detection.

## Reproduce

Dev scripts in `dev/`: `inspect_data`, `build_dataset`, `diagnostics`,
`benchmark_loso` (kNN-era), `benchmark_differential`, `perturbation_test`,
`validate_differential_tf`, `make_figures`. Rebuild the bundle:
`training.build_bundle(dataset_npz, out_path, ridge_lambda=1.0)`.

Lung reference: `dev/build_lung_reference.py` (dataset + bundle) and
`dev/validate_lung_reference.py` (checks A–E). Both take absolute-path flags and
default to the LungMAP inputs.
