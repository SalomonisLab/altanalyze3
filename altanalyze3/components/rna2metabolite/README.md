# rna2metabolite (AML)

Impute **metabolite abundance from bulk/pseudobulk RNA** for AML, on the same plumbing as
`rna2grn` / `rna2adt`. AML-specific and independent of the lung `rna2lipid` models.

## Final model (the only classifier shipped)
- One **per-metabolite ridge** (ElasticNet `l1_ratio=0`), **1000 genes/target, α=100**.
- Trained on **CPTAC AML** (Nature Cancer 2026; 84 RNA+metabolomics cases; GDC CPTAC-3 STAR counts).
- Genes: protein-coding, GDC `tpm_unstranded ≥ 5` in ≥ 3 cases (12,416-gene universe).
- Prediction is per-molecule and linear — `ŷ_t = b_t + w_t · z[genes_t]` — **not** a per-sample
  neighbor lookup; on held-out CV it beats k-NN retrieval on the same genes (0.268 vs 0.223).
- Bundle: `artifacts/rna2metabolite_aml_bundle.pkl.gz` (2,533 metabolites × 12,416 genes).

## Held-out performance (5-fold CV)
| median Spearman | imputable (Sp>0.3) | strong (Sp>0.5) | classification acc / AUROC (imputable) |
|---|---|---|---|
| 0.267 | 1,084 | 313 | 0.663 / 0.718 |

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
- `artifacts/rna2metabolite_aml_bundle.pkl.gz` — final classifier
- `VALIDATION.md` — methods + held-out evaluation
