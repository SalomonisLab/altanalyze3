# cellHarmony_lite.py

## Purpose
`cellHarmony_lite.py` performs reference-based cell annotation of single-cell RNA-seq data using cosine similarity. It supports:
- Merging multiple 10x Genomics `.h5` or `.mtx` files
- QC filtering
- Log2 normalization (counts per 10,000)
- Cell type/state assignment using a marker-based reference
- Export of assignments and normalized matrices

---

## Inputs

### Required
| Argument | Description |
|----------|-------------|
| `--h5dir` | Folder containing 10x `.h5` or `.mtx` files (or a single file) |
| `--refdir` | A tab-delimited file of reference expression centroids (genes x cell populations) |

### Optional
| Argument | Description |
|----------|-------------|
| `--h5ad` | Use a preprocessed `.h5ad` file instead of raw 10x files |
| `--outdir` | Output folder (default: `output/`) |
| `--cptt` | Export log2-normalized expression matrix (gene x cell) |
| `--export_h5ad` | Save combined Seurat object (`.h5ad`) |
| `--min_genes` | Min genes per cell (default=500) |
| `--min_cells` | Min cells per gene (default=0) |
| `--min_counts` | Min counts per cell (default=1000) |
| `--mit_percent` | Max mitochondrial percent cutoff (default=10%) |

---

## Output Files

| File | Description |
|------|-------------|
| `cellHarmony_lite_assignments.txt` | Cell barcode → reference cell population + alignment score |
| `CPTT_matrix.txt` | (Optional) Dense log2(CPTT) matrix |
| `combined_qc_normalized.h5ad` | (Optional) Processed h5ad with metadata |

---

## Quick Tutorial

### 1. Reference Format (`reference.txt`)
```tsv
UID	Monocyte	T-cell	B-cell
GeneA	3.2	1.1	0.2
GeneB	2.9	0.1	4.3
