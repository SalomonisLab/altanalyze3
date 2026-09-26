# cellHarmony_lite.py

## Overview
`cellHarmony_lite.py` aligns single-cell RNA-seq datasets to a reference panel of cell-state centroids. It can ingest one or more 10x Genomics runs, re-use an existing `.h5ad`, perform light quality control, and score each query cell (or optional metacell) against the reference. Alignment supports cosine similarity (default) or a Pearson/z-difference “classic” mode, and results can be exported as text tables or AnnData objects for downstream analyses.

## Supported Input Formats
- **10x Genomics `.h5`** archives created by `cellranger count` or `cellranger arc`
- **10x Matrix Market directories** containing `matrix.mtx[.gz]`, `barcodes.tsv[.gz]`, and `features/genes.tsv[.gz]`
- **Existing `.h5ad`** objects (skips raw file assembly and QC if already processed)
- Mixed batches of the above; each run is tagged with a sample and group identifier automatically
- **altanalyze3 sparse stream** (`--sparse_stream`): a cells x genes CSR matrix, names and cell annotations read from a pipe or file, so no intermediate `.h5ad` or `.mtx` is written. `sparse_stream.py` defines the layout and checks it; `sparse_stream_from_seurat.R` streams a Seurat assay (see "Streaming a Seurat object")

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
| `--sparse_stream` | Read the input from an altanalyze3 sparse stream: a path, a named pipe, or `-` for standard input. An obs column `Library` sets the ambient correction unit. Replaces `--h5dir` and `--h5ad`. |
| `--pseudobulk_sample_col` | After alignment, sum the aligned cells into sample x cell-state pseudobulks with `aggregate.pseudobulk_h5ad` (`build_pseudobulk_from_adata`), using this obs column as the sample. Writes `pseudobulk_counts.h5ad` and no single-cell h5ad. |
| `--pseudobulk_min_cells` (10) | Minimum cells per pseudobulk. |
| `--pseudobulk_layers` | Comma-separated extra layers to sum, for example `soupx_raw` for the uncorrected counts after `--ambient_correct_cutoff`. |

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
| `logs/cellHarmony-lite_<timestamp>.log` | Always | Execution log capturing console output and command-line parameters for reproducibility. |
| `pseudobulk_counts.h5ad` | `--pseudobulk_sample_col` | Sample x cell-state pseudobulks: summed `layers["counts"]`, `X` = counts divided by each pseudobulk's total, one summed layer per `--pseudobulk_layers` name. |

### Streaming a Seurat object
R reads the RDS and writes the dgCMatrix slots unchanged; cellHarmony reads them from the pipe:

```bash
Rscript altanalyze3/components/cellHarmony/sparse_stream_from_seurat.R \
    --zip archive.zip --member seurat_object.RDS --assay RNA --slot counts \
    --obs-cols sample --library-col capture \
  | python -m altanalyze3.components.cellHarmony.cellHarmony_lite --sparse_stream - \
    --refdir Hs-MarrowAtlas-L3M.txt --outdir out --ambient_correct_cutoff auto \
    --pseudobulk_sample_col sample --pseudobulk_layers soupx_raw
```

The reader checks the magic, the trailer, the row pointers, the column index range, and each
cell's sum against the total the producer sent, so a truncated or corrupt stream stops the run.

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
