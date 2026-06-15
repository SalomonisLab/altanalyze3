# rna2metabolite / rna2lipid (AML) — methods and held-out evaluation

AML-specific imputation of metabolite and lipid abundance from bulk RNA-seq (CPTAC AML,
Nature Cancer 2026). Only the **final** held-out-optimal classifiers are shipped
(metabolite: 1000 genes/target; lipid: 120 genes/target; ridge α=100).

## Data
- **RNA**: GDC `CPTAC-3` open STAR counts (matched subsets 84 metabolomics / 87 lipidomics cases).
- **Targets**: untargeted MS features collapsed to one best-detected measurement per molecule
  (HILIC+RP → 2,533 metabolites; pos+neg ESI → 1,009 lipids).

## Model
- **Per-molecule** feature set: mechanistic-prior genes (Recon3D reaction→gene for metabolites;
  Recon3D lipid-class subsystems ∪ lung `rna2lipid` genes for lipids) ranked by within-fold
  |Spearman|, then filled with top genome-wide correlated genes; no HVG.
- **Estimator**: per-molecule **ridge** (= ElasticNet `l1_ratio=0`), α=100. Sparsity is imposed by
  the gene pre-selection, so L1 is redundant; heavily-regularized ridge beats ElasticNet/lasso at
  these gene counts.
- **Normalization**: CP10k+log1p (all-gene library size); per-gene z-score with training μ/σ.
- **CV**: 5-fold over cases; all selection/standardization within the training fold (leakage-free).

## Held-out results (out-of-fold)
| | median Spearman | imputable (Sp>0.3) | Sp>0.5 | class. acc / AUROC (imputable) |
|---|---|---|---|---|
| metabolites (2,478 eval) | 0.267 | 1,084 | 313 | 0.663 / 0.718 |
| lipids (1,009) | 0.355 | 625 | 245 | 0.678 / 0.723 |

- **Not circular**: in-sample median Spearman ≈ 0.97/0.83; the out-of-fold median (0.27/0.36) is the
  honest value. The deliverable is the imputable subset (median molecule has negative R²).
- **Not retrieval**: on held-out CV the per-molecule regression beats k-NN retrieval on the same
  genes (metabolite 0.268 vs 0.223; lipid 0.355 vs 0.308) — it imputes a function of expression,
  not a per-sample neighbor average.
- **Beyond subtype**: ≤18% of MS variance is between molecular subtype; imputable molecules retain
  within-subtype Spearman ≈ 0.24–0.32 (genuine within-subtype biology).

## Limitations
1. Internal CV only (no external matched RNA+MS cohort).
2. 84–87 cases bound power; the imputable subset is what survives.
3. Cell-composition/blast-% confounds not yet regressed out; high-gene accuracy partly rides
   myeloid cell-content genes (BPI, CEACAM1, ELANE) over biosynthetic enzymes — accuracy vs mechanism.
4. Bundle trained on bulk RNA; applied to scRNA pseudobulks it can extrapolate — use rank/relative
   values and filter by `var['heldout_spearman']` and `obs['n_cells']`.

Provenance and the full analysis (build, sweeps, figures) live in the source project
`Human-MS-impute`; only the final classifiers and inference code are integrated here.
