# How to use scALABLE-discover

scALABLE-discover clusters single-cell RNA data without a reference atlas. ICGS3 finds the
clusters; the Explore and Chat tabs are the ones scALABLE-web uses. The interface has three
tabs: `Run`, `Explore` and `Chat`. This version runs no differential expression and no
modality imputation.

The method behind each number is in [README.md](README.md).

## Run tab

### 1. Upload

| Rule | Value |
| --- | --- |
| File types | `.h5` (Cell Ranger) or `.h5ad` (AnnData); use the same format within a job |
| Files per job | up to 7, one per sample |
| Species | human or mouse |

1. Choose `Species`.
2. Click `Add sample` once per file. Give each row a `Sample name` and a file.
3. Click `Upload`.

H5AD import uses disk-backed arrays when the combined input exceeds 100,000 cells,
its uncompressed matrices total at least 1 GiB, or multiple H5AD files are uploaded.
The size includes `X`, layers, `raw` and stored embeddings/graphs. Multiple H5AD files
merge on disk, preserving their feature union and sample identities. Later clustering
and UMAP steps still use RAM. Completed jobs build a separate disk-backed Explore
bundle at 10,000 cells or more by default.

### 2. QC and ICGS3 clustering

| Control | Default | Meaning |
| --- | --- | --- |
| `Min genes` | 500 | drop a cell with fewer detected genes |
| `Min counts` | 1000 | drop a cell with fewer UMI counts |
| `Min cells` | 0 | drop a gene detected in fewer cells |
| `Mito %` | 15 | drop a cell at or above this mitochondrial percent |
| `Ambient RNA correction` | No | `Yes` estimates and subtracts ambient RNA per sample |
| `UMAP fitting` | Accelerated | fit representative cells using the complete MarkerFinder panel, then map every remaining cell; target 30,000 fitting cells with a minimum of 200 per cluster; select All cells for the original full feature fit |
| `Max K` | — (automatic) | optional integer of at least 2; sets ICGS3's target NMF K (`--nmf-k`); leave blank for automatic rank estimation |

QC uses scALABLE-web's shared filtering code; UMAP fitting and Max K configure ICGS3. Click
`Save QC and run`.

Max K sets the target NMF rank. ICGS3's later marker and SVM steps can yield fewer
final clusters; the field does not force the final displayed cluster count.

Accelerated UMAP includes every clustered cell in Explore and leaves clustering and marker
discovery unchanged. The fitting set contains at least 200 cells per cluster, or every cell
in smaller clusters, with additional cells sampled within larger clusters. Its budget
expands when necessary to cover every cluster. Jobs with at most 30,000 cells fit all cells automatically. The shape
of an accelerated embedding can differ from a full fit. Logged h5ad inputs without a counts
layer are passed to ICGS3 as normalized expression, with no additional log transform.
For log-normalized inputs, the shared pipeline skips count-based QC filtering; the QC
panel explicitly reports that the thresholds were not applied.

### What runs

1. QC filters the cells and normalizes the counts.
2. ICGS3 clusters the QC-retained cells. Louvain sampling keeps 10,000 cells, PageRank keeps
   5,000. UDON NMF and the first MarkerFinder pass use those 5,000. The SVM then assigns
   every QC-retained cell, and the final MarkerFinder markers, UMAP and GO-Elite BioMarkers
   labels use every clustered cell.
3. GO-Elite BioMarkers names each cluster, for example `AT2 Cells_c2` for cluster C2.
4. NetPerspective draws a marker network per cell state.
5. fastComm scores receptor-ligand communication between cell states.

While the job runs, the percentage and the line under the progress bar show the current
stage: QC, then ICGS3 steps 1 to 10, then writing the results. Step 10 reports UMAP
fitting and mapping remaining cells separately. When it finishes, the line
gives how many QC-retained cells ICGS3 placed in a cell state. ICGS3 leaves a cell
unassigned when its best SVM score is 0 or lower. An unassigned cell appears in no view.

Upload tip: scALABLE-web offers "Download unaligned QC-passed cells (h5ad)" for the cells
that passed QC but fell below its alignment cutoff. Upload that file here to cluster them.
Its counts are already ambient-corrected when the scALABLE-web job corrected them, so set
`Ambient RNA correction` to `No` for it.

## Explore tab

`Cell-state layer`, at the top of the left panel, sets how every view names cells:

| Layer | Names |
| --- | --- |
| Predicted cell states (default) | the GO-Elite BioMarkers label of each cluster, column `cell_state_predicted` |
| Clusters | C1, C2, ..., column `cluster` |

Changing the layer reloads the page. Two panels follow, each with its own plot type,
filters and `Download PDF`.

| Plot type | What it shows |
| --- | --- |
| `UMAP cell states` | every clustered cell on the UMAP, coloured by cell state |
| `Cell frequency` | each sample's fraction of cells per cell state |
| `UMAP` | one gene's value per cell on the UMAP |
| `Violin` | one gene's distribution per cell state |
| `DotPlot` | a gene set by cell state: dot size is percent expressing, colour is mean |
| `CombPlot` | a gene set by cell state and sample |
| `MarkerHeatmap` | ICGS3's MarkerFinder markers in Morpheus |
| `MarkerNetwork` | known interactions among one cell state's markers |
| `GO-Elite BioMarkers` | one cell state's enriched BioMarkers terms: z-score against FDR; hover for the overlapping genes |
| `Cell communication` | fastComm receptor-ligand signalling between cell states |
| `Pathway` | pathways holding a cell state's markers |

`Color by` on `UMAP cell states` offers the other layer, `original_NMF_cluster`, and every
other categorical cell annotation.

### Downloads

| Button | Content |
| --- | --- |
| combined_h5ad | every QC gene, normalized, for the clustered cells, with ICGS3 labels and UMAP |
| MarkerFinder ZIP | ICGS3's marker tables, centroids and marker networks |
| log | the pipeline log |

These are the default `minimal` exports. A server set to `SCALABLE_DISCOVER_EXPORTS=full`
adds cluster assignments, the ICGS3 results ZIP, the cell communication ZIP and the run
parameters.

## Chat tab

Chat works as in scALABLE-web. It answers marker, comparison-between-clusters and
where-expressed questions from this job's data. A question about a comparison between
groups answers that no comparison has run.
