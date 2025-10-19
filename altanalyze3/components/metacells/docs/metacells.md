# metacells/main.py

## Purpose
`metacells/main.py` constructs fast, size-controlled metacells from single-cell RNA-seq AnnData files. It prioritizes:
- Summing raw counts to create metacell expression profiles
- Flexible grouping (global, sample-restricted, or boundary-aware)
- Optional random aggregation for bootstrapped metacells
- Full provenance of member cell barcodes in the output

---

## Inputs

### Required
| Argument | Description |
|----------|-------------|
| `--input` | Input `.h5ad` file containing filtered scRNA-seq data |
| `--output` | Destination `.h5ad` file for aggregated metacells |

### Optional (Key Controls)
| Argument | Description |
|----------|-------------|
| `--mode` | `global` (default), `boundary`, or `sample` strategy |
| `--sample-column` | Name of `adata.obs` column with sample/donor IDs |
| `--cell-type-column` | Name of `adata.obs` column with cell-state labels |
| `--boundary-columns` | Explicit columns used to partition cells before clustering |
| `--restrict-sample` | Restrict analysis to specific sample IDs (repeatable) |
| `--restrict-cell-type` | Restrict to specific cell-type values (repeatable) |
| `--graph-algorithm` | `leiden`, `louvain`, `kmeans`, or `random` |
| `--target-metacell-size` | Desired cells per metacell (default `50`) |
| `--min-metacell-size` | Minimum cells per metacell after refinement |
| `--max-metacell-size` | Maximum cells per metacell after refinement |
| `--preserve-small-clusters` | Keep undersized clusters instead of merging |
| `--random-metacell-count` | Number of metacells per group when `--graph-algorithm random` |
| `--random-cells-per-metacell` | Cells sampled per random metacell (default `5`) |
| `--random-sampling-with-replacement` | Allow sampling the same cell multiple times |
| `--n-top-genes` | Number of highly-variable genes to keep (default `3000`) |
| `--n-pcs` | PCs for KNN graph construction (default `50`) |
| `--neighbor-method` | KNN backend: `auto`, `umap`, or `gauss` |
| `--neighbor-metric` | Distance metric for neighbors (default `euclidean`) |
| `--expression-layer` | Layer to aggregate instead of `X` |
| `--use-raw` | Use `adata.raw` for graph construction & aggregation |
| `--threads` | Parallel workers for Scanpy tasks |
| `--verbose` | Enable debug logging |

---

## Outputs

| File / Location | Description |
|-----------------|-------------|
| `--output` (`.h5ad`) | Metacell expression matrix with metacell metadata in `obs` and original gene metadata in `var` |
| `.uns['metacell_membership']` | DataFrame mapping each metacell ID to constituent barcodes (plus sample / cell-type annotations if available) |
| `.uns['metacell_parameters']` | Serialized CLI arguments used for the run (reproducibility) |

By default, metacell expression values are raw count sums (`sum` aggregation only).

---

## Synthetic Example Dataset
Utility script: `altanalyze3/components/tests/metacell_synthetic_data.py`

```bash
python3.11 - <<'PY'
from altanalyze3.components.tests.metacell_synthetic_data import write_synthetic_h5ad
path, _ = write_synthetic_h5ad()
print(path)
PY
```

Generates `altanalyze3/components/tests/data/metacell_synthetic.h5ad` with:
- 72 cells (3 samples × 3 cell types × 8 replicates)
- 15 genes
- Layers: `counts` (integer-scaled) and `normalized`

---

## Quick Tutorials

### 1. Global Metacells (KMeans)
```bash
python3.11 altanalyze3/components/metacells/main.py \
  --input altanalyze3/components/tests/data/metacell_synthetic.h5ad \
  --output metacells_global.h5ad \
  --target-metacell-size 12 \
  --min-metacell-size 6 \
  --max-metacell-size 18 \
  --sample-column sample_id \
  --cell-type-column cell_type \
  --graph-algorithm kmeans \
  --neighbor-method auto \
  --n-top-genes 40 \
  --n-pcs 8
```

**Expected outcome**
- ~7 metacells (`n_obs` ≈ 7) combining all 72 cells
- `.obs['metacell_size']` spans the configured min/max (e.g., `[18, 6, 7, 9, 17, 9, 6]`)
- `.uns['metacell_membership']` contains 72 rows with original barcodes and dominant annotations

---

### 2. Boundary-Aware Metacells per Sample/Cell Type
```bash
python3.11 altanalyze3/components/metacells/main.py \
  --input altanalyze3/components/tests/data/metacell_synthetic.h5ad \
  --output metacells_boundary.h5ad \
  --mode boundary \
  --sample-column sample_id \
  --cell-type-column cell_type \
  --target-metacell-size 8 \
  --min-metacell-size 4 \
  --max-metacell-size 10 \
  --graph-algorithm kmeans \
  --neighbor-method auto \
  --n-top-genes 50 \
  --n-pcs 10
```

**Expected outcome**
- One or more metacells per `(sample_id, cell_type)` subgroup (9 groups in total)
- Each metacell contains only cells from its boundary group
- Output membership table includes `sample_id` and `cell_type` columns reflecting group assignments

---

### 3. Random Aggregation for Bootstrapping
```bash
python3.11 altanalyze3/components/metacells/main.py \
  --input altanalyze3/components/tests/data/metacell_synthetic.h5ad \
  --output metacells_random.h5ad \
  --mode boundary \
  --sample-column sample_id \
  --cell-type-column cell_type \
  --graph-algorithm random \
  --random-metacell-count 5 \
  --random-cells-per-metacell 4 \
  --random-sampling-with-replacement \
  --target-metacell-size 4
```

**Expected outcome**
- For each `(sample_id, cell_type)` boundary, exactly 5 metacells are produced
- `.obs['metacell_size']` equals the requested `random-cells-per-metacell` (allowing repeats when sampling with replacement)
- Membership captures the random draws, enabling downstream resampling analyses

---

## Notes
- Sum aggregation is enforced to ensure each metacell reflects combined raw counts.
- If employing `leiden` or `louvain`, ensure Scanpy’s clustering backends (including `scikit-misc` for `kmeans`) are installed.
- The command-line interface mirrors other AltAnalyze components, allowing integration into existing pipelines or CWLs.

