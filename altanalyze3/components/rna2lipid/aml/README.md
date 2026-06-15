# rna2lipid · aml

AML-specific **lipid imputation from bulk/pseudobulk RNA**. This is a **separate model and
code path** from the lung/IPF `rna2lipid` multitask elastic-net model — that model and its
code (`../api.py`, `../artifacts/rna2lipid_multitask_enet_reproducible.pkl`) are untouched.

## Final model (the only classifier shipped)
- One **per-lipid ridge** (ElasticNet `l1_ratio=0`), **120 genes/target, α=100**.
- Trained on **CPTAC AML** (Nature Cancer 2026; 87 RNA+lipidomics cases; GDC CPTAC-3 STAR counts).
- Genes: protein-coding, GDC `tpm_unstranded ≥ 5` in ≥ 3 cases (12,472-gene universe; prior =
  Recon3D lipid-class subsystem genes ∪ the 1,303 lung `rna2lipid` genes).
- Per-molecule linear prediction (not per-sample neighbor lookup); beats k-NN retrieval on held-out
  CV (0.355 vs 0.308).
- Bundle: `artifacts/rna2lipid_aml_bundle.pkl.gz` (1,009 lipids × 12,472 genes).

## Held-out performance (5-fold CV)
| median Spearman | imputable (Sp>0.3) | strong (Sp>0.5) | classification acc / AUROC (imputable) |
|---|---|---|---|
| 0.355 | 625 | 245 | 0.678 / 0.723 |

Most-imputable classes: glycosphingolipids Hex2Cer/Hex3Cer, PA, PS, GM3.

## Usage
```python
from components.rna2lipid.aml import load_bundle
b = load_bundle()
res = b.predict_from_h5ad("pseudobulk_counts.h5ad")   # samples x lipids
out = b.impute_anndata(adata)                          # AnnData: X=imputed, obs carried over
```
```bash
python -m components.rna2lipid.aml.cli model-info
python -m components.rna2lipid.aml.cli predict-h5ad --input pb.h5ad --output imp.csv --h5ad-out imp.h5ad
```
Genes match by symbol; absent genes take the training mean; CP10k+log1p + z-score applied internally
(pass `--normalized` for pre-normalized input).
