# rna2metabolite (AML)

Impute **metabolite abundance from bulk/pseudobulk RNA** for AML, on the same plumbing as
`rna2grn` / `rna2adt`. AML-specific and independent of the lung `rna2lipid` models.

## Default model
- One **per-metabolite ridge** (ElasticNet `l1_ratio=0`), **1000 genes/target, α=100**.
- Trained on **CPTAC AML** (Nature Cancer 2026; 84 RNA+metabolomics cases; GDC CPTAC-3 STAR counts).
- Genes: protein-coding, GDC `tpm_unstranded ≥ 5` in ≥ 3 cases (12,416-gene universe).
- Prediction is per-molecule and linear — `ŷ_t = b_t + w_t · z[genes_t]` — **not** a per-sample
  neighbor lookup; on held-out CV it beats k-NN retrieval on the same genes (0.268 vs 0.223).
- Bundle: `artifacts/rna2metabolite_aml_NA30_20261006_bundle.pkl.gz`
  (**2,023 metabolites × 12,416 genes**). The user explicitly authorized removing
  features with **more than 30% NA values** across the original 84 training cases.
  Missingness is evaluated before filling; 25 missing values are retained and
  26 missing values are removed. Exactly 30% is retained.
- **510 targets removed**: 463 previously fitted targets plus all 47 unfitted
  targets. Retained target order, coefficients, intercepts, RNA scalers and all
  84 training-case identities are unchanged. No additional fitting or filling
  was needed: all 2,023 retained models are already fitted and finite.
- Original `artifacts/rna2metabolite_aml_bundle.pkl.gz` remains available for
  explicit rollback. `missingness_release.json` records every retained/removed
  identity and the baseline/release hashes; loading the release checks its hash,
  complete ordered panels and original ridge settings.

## Original held-out performance (5-fold CV)
| median Spearman | imputable (Sp>0.3) | strong (Sp>0.5) | classification acc / AUROC (imputable) |
|---|---|---|---|
| 0.267 | 1,084 | 313 | 0.663 / 0.718 |

This table describes the original panel. The filtered release retains historical
per-target CV values only for its retained targets and recomputes their summary
in the bundle metadata. No new CV or classification evaluation was performed.

The deliverable is the **imputable subset** (the median metabolite is not imputable). Per-molecule
held-out Spearman/R² and the genes used are in `var` of the output and in the bundle metadata.

## Usage
```python
from components.rna2metabolite import load_bundle
b = load_bundle()                                   # default bundle
res = b.predict_from_h5ad("pseudobulk_counts.h5ad") # rows = pseudobulks, var = gene symbols
res.predictions                                     # samples x metabolites
out = b.impute_anndata(adata)                       # AnnData: X=imputed, obs carried over, var=reliability
```
```bash
python -m components.rna2metabolite.cli model-info
python -m components.rna2metabolite.cli predict-h5ad --input pb.h5ad --output imp.csv --h5ad-out imp.h5ad
```
Input genes match the model by **symbol**; absent genes take the training mean (z=0). Normalization
(CP10k+log1p on the all-gene library size, then z-score with training μ/σ) is applied internally; pass
`--normalized` if the input is already CP10k+log1p.

## Files
- `api.py`, `cli.py`, `_impute.py` — inference engine (no training code)
- `artifacts/rna2metabolite_aml_NA30_20261006_bundle.pkl.gz` — default classifier
- `artifacts/rna2metabolite_aml_bundle.pkl.gz` — original rollback classifier
- `VALIDATION.md` — methods + held-out evaluation

## Verified source targets (2026-10-06)

Saved GitHub provenance points to the original source project at
`/Users/saljh8/Dropbox/Collaborations/Grimes/Human-MS-impute`.
`code/build_unique_ms_tables.py:77` imports processed workbook Tables21/22
unchanged, explicitly documenting them as log2/median-centered. The original
target loader and ridge builder introduce no additional target log transform or
target z-score. RNA CP10k/log1p/z-scoring is a separate input transformation.

All **2,533 target identities × 84 training cases** match a reconstruction of
the original workbook extraction and existing assay-selection rule. Maximum
absolute numeric disagreement is **1.78e-15**. All packaged model arrays and
identities match the original saved model after the documented packaging casts.
The normalized linear representation of a documented log2 target is `2**y`;
this establishes a relative signal scale, not an upstream absolute-intensity
inverse. Original workbook and model files were not changed.

The trace also found **47 targets with empty coefficient arrays and NaN
intercepts** in both original and packaged artifacts; **2,486 are fitted**.
Inference preserves the 47 as NaN in its prediction dataframe. Existing
`impute_anndata` selects fitted targets, while scALABLE's existing finalizer
replaces NaNs with zero. No exclusion, zero filling, retraining or policy change
was performed by this audit. A separate user question about auditing these
47 targets was answered with authorization to impute NaNs and update the model.
The subsequent observation audit found **37 targets with 1–11 measured values**
in the 84 training cases and **10 targets with no observations** in either those
cases or all 95 source cases. The original fitting rule requires 12 measured
values. The user approved the established lipid missing-value method for the 37:
fill missing entries with each target's observed median on the normalized log2
scale, then fit the original metabolite ridge procedure (alpha100, 1000 selected
genes), preserving the 2,486 existing fitted models.

The separately authorized raw-assay diagnostic found each of the ten all-missing
identities in exactly one original Table21/Table22 row. Those rows have no numeric
sample measurements, no sample-cell formulas and no alternate exact-identity
assay rows. No measurements were recovered. See
`provenance/2026-10-06/all_missing_assay_row_audit.json` and its reproducible
`audit_all_missing_assay_rows.py`. The subsequent user instruction,
“If >30% values are NA, remove that feature,” resolves the all-missing policy:
all 47 previously unfitted targets meet the authorized exclusion rule. The
original full-panel model and source tables remain unchanged; the new default
has the 2,023 explicitly retained targets. No sparse-target refitting was needed.
Do not treat substituted zeros as measured or predicted
abundances. See `unfitted_target_observation_audit.json` beside the scale audit.

See `provenance/2026-10-06/target_scale_trace.json` and its reproducible
`trace_target_scale.py`. Runtime metadata was not changed during the original
trace. The new release and scALABLE config now explicitly declare log2 with
zero offset in normalized relative space. Linear signal is `2**y`, with no
subtraction of one after centering and no lipid/metabolite panel-total scaling.
The complete missingness table is
`provenance/2026-10-06/target_missingness_policy.tsv` (all 2,533 original identities).
