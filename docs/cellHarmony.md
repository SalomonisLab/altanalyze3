# cellHarmony_lite.py

## Overview
`cellHarmony_lite.py` aligns single-cell RNA-seq datasets to a reference panel of cell-state centroids. It can ingest one or more 10x Genomics runs, re-use an existing `.h5ad`, perform light quality control, and score each query cell (or optional metacell) against the reference. Alignment supports cosine similarity (default) or a Pearson/z-difference “classic” mode, and results can be exported as text tables or AnnData objects for downstream analyses.

## Supported Input Formats
- **10x Genomics `.h5`** archives created by `cellranger count` or `cellranger arc`
- **10x Matrix Market directories** containing `matrix.mtx[.gz]`, `barcodes.tsv[.gz]`, and `features/genes.tsv[.gz]`
- **Existing `.h5ad`** objects (skips raw file assembly and QC if already processed)
- Mixed batches of the above; each run is tagged with a sample and group identifier automatically

---

## Core Command-Line Parameters
These switches cover the typical alignment workflow. Defaults are shown in parentheses.

| Parameter | Description |
|-----------|-------------|
| `--h5dir` | Directory of input runs (10x `.h5`, Matrix Market triplets, or nested folders). Use when starting from raw outputs. |
| `--h5ad` | Optional shortcut to supply a prebuilt `.h5ad` instead of `--h5dir`. When present, the file is loaded directly. |
| `--refdir` | Tab-delimited reference file (`genes × populations`) containing cellHarmony centroids used for alignment. |
| `--outdir` (`output`) | Destination folder for all reports, matrices, and AnnData exports (created if missing). |
| `--alignment_mode` (`cosine`) | Alignment algorithm. Use `cosine` for cosine similarity or `classic` for Pearson correlation with z-difference scoring. |
| `--align_cutoff` (none) | Minimum alignment score required to retain a cell/metacell in the assignment table. |
| `--export_h5ad` | Write the QC-filtered, normalized AnnData object (`combined_qc_normalized.h5ad`). |
| `--cptt` | Export a dense log-normalized expression matrix (`CPTT_matrix.txt`, genes × cells/metacells). |
| `--metacell-align` | Aggregate cells into metacells prior to alignment; helpful for very large or noisy datasets. |

> **Tip:** Provide either `--h5dir` or `--h5ad`. If both are supplied, the explicit `.h5ad` takes precedence.

---

## Accessory Parameters
Additional options are grouped by theme to keep the core workflow succinct.

### Quality Control & Normalisation
| Parameter | Description (default) |
|-----------|-----------------------|
| `--min_genes` | Minimum detected genes per cell (200). |
| `--min_cells` | Minimum cells referencing a gene before it is retained (3). |
| `--min_counts` | Minimum total counts per cell (500). |
| `--mit_percent` | Upper mitochondrial percentage threshold (10). |
| `--ambient_correct_cutoff` | If set, run SoupX ambient RNA correction with this contamination fraction (rho) **before** QC and alignment. Corrected counts and a summary TSV are written under `<outdir>/soupx/`. |
| `--generate_umap` | Run Scanpy HVG/PCA/neighbor graph steps and store UMAP coordinates plus marker rankings. |
| `--save_adata` | Persist the AnnData after optional UMAP/clustering (`combined_with_umap_and_markers.h5ad`). |
| `--unsupervised_cluster` | Compute Leiden clusters prior to exporting results. |

### Metadata & Identifier Handling
| Parameter | Description |
|-----------|-------------|
| `--append_obs` | Append values from a specified `.obs` column (e.g. donor ID) to each barcode. |
| `--gene_translation` | Two-column TSV mapping source gene IDs to symbols (e.g. Ensembl → HGNC). Applied to all input modalities—including `.h5ad`—before alignment so reference symbols match. |

### Metacell Construction
When `--metacell-align` is active, the following knobs control metacell generation. Defaults are shown in parentheses.

| Parameter | Description |
|-----------|-------------|
| `--metacell-target-size` | Desired number of cells per metacell (50). |
| `--metacell-min-size` / `--metacell-max-size` | Hard bounds on metacell size (25 / 100). |
| `--metacell-algorithm` | Clustering routine: `kmeans`, `leiden`, `louvain`, or `random` (`kmeans`). |
| `--metacell-neighbors` | Nearest neighbors for the metacell graph (30). |
| `--metacell-hvg` | Highly variable genes used during clustering (3000). |
| `--metacell-pcs` | Principal components retained for clustering (50). |
| `--metacell-random-count` / `--metacell-random-cells` | Controls for the `random` algorithm (50 metacells of 5 cells each). |
| `--metacell-random-replacement` | Sample with replacement when forming random metacells. |
| `--metacell-random-state` | Random seed for reproducibility (0). |

---

## Generated Outputs
| File | When Produced | Contents |
|------|----------------|----------|
| `cellHarmony_lite_assignments.txt` | Always | Cell (or metacell) barcode, matched reference population, alignment score. |
| `CPTT_matrix.txt` | `--cptt` | Dense log-normalized expression matrix (genes × cells/metacells). |
| `combined_qc_normalized.h5ad` | `--export_h5ad` | QC-filtered AnnData with original counts in `layers["counts"]` and normalized expression in `X`. |
| `combined_with_umap_and_markers.h5ad` | `--save_adata` or `--generate_umap` | Annotated AnnData containing UMAP coordinates, optional Leiden labels, and marker rankings. |
| `metacells.h5ad` | `--metacell-align` | AnnData describing generated metacells plus `uns["metacell_membership"]` for barcode → metacell mapping. |

All outputs are written to `--outdir` (default `output/`). Console logs summarise file loading, gene translation, QC filters, alignment statistics, and optional metacell construction.

---

## Reference File Layout Example
Supply the reference centroids via `--refdir` as shown below. Gene identifiers should match the translated symbols present in your query data.

```tsv
UID	Monocyte	T-cell	B-cell
GeneA	3.2	1.1	0.2
GeneB	2.9	0.1	4.3
```

If your dataset uses alternative identifiers (e.g. Ensembl IDs), convert them on the fly by supplying `--gene_translation` with a two-column TSV of source → symbol mappings.

---

## Typical Workflow
1. Prepare a reference centroid matrix of genes × populations.
2. Gather raw 10x outputs inside a directory (or point to an existing `.h5ad`).
3. Run:
   ```bash
   python cellHarmony_lite.py \
     --h5dir data/10x_runs \
     --refdir refs/pediatric_lung_reference.txt \
     --alignment_mode classic \
     --export_h5ad --cptt
   ```
4. Review `cellHarmony_lite_assignments.txt` for population assignments, check optional UMAP/clustering outputs, and load the exported `.h5ad` into downstream differential pipelines.

Adjust the accessory parameters to tune QC thresholds, append metadata into barcode IDs, translate gene identifiers, or condense large datasets into metacells before alignment.
