# How to use scALABLE-discover

scALABLE-discover identifies transcriptionally distinct cell states and their marker genes
in human and mouse single-cell RNA-sequencing data. Upload samples in `Run`, examine
cell states and gene expression in `Explore`, and query the results in `Chat`.

The analysis uses ICGS3, a Python 3 implementation optimized for ultrafast cell-state
discovery and accessible through a web interface. It builds on the
[ICGS approach (Venkatasubramanian et al., 2020)](https://academic.oup.com/bioinformatics/article/36/12/3773/5811229),
which combined HOPACH, sparse non-negative matrix factorization, cluster fitness and SVM
classification to resolve rare and common cell states while limiting donor and batch
influences. The published atlas analyses showed that PageRank sampling preserved rare
and closely related cell types and enabled discovery of additional distinct populations.
ICGS3 replaces the legacy HOPACH workflow with sparse processing and UDON-derived feature
selection and NMF, retaining PageRank, MarkerFinder cluster fitness and SVM assignment.

[README.md](README.md) describes the methods, parameters, outputs and validation results.

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

QC applies the selected thresholds; UMAP fitting and Max K configure ICGS3. Click
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
2. ICGS3 clusters the QC-retained cells. For larger datasets, Louvain sampling selects up
   to 10,000 representative cells, followed by PageRank sampling of up to 5,000 cells.
   UDON-derived NMF and the first MarkerFinder pass use the sampled cells. The SVM then
   scores every QC-retained cell for cluster assignment. The final MarkerFinder markers,
   UMAP and GO-Elite BioMarkers labels use every clustered cell.
3. GO-Elite BioMarkers names each cluster, for example `AT2 Cells_c2` for cluster C2.
4. NetPerspective generates a marker network per cell state.
5. fastComm scores receptor-ligand communication between cell states.

While the job runs, the percentage and the line under the progress bar show the current
stage: QC, then ICGS3 steps 1 to 10, then writing the results. Step 10 reports UMAP
fitting and mapping remaining cells separately. When it finishes, the line
gives how many QC-retained cells ICGS3 placed in a cell state. ICGS3 leaves a cell
unassigned when its best SVM score is 0 or lower. An unassigned cell appears in no view.

## Explore tab

`Cell-state layer`, at the top of the left panel, sets how every view names cells:

| Layer | Names |
| --- | --- |
| Predicted cell states (default) | the GO-Elite BioMarkers label of each cluster, column `cell_state_predicted` |
| Clusters | C1, C2, ..., column `cluster` |

Changing the layer reloads the page. Each of the two plotting panels has its own plot type,
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
| `Pathway` | pathways containing a cell state's markers |

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

Chat answers questions about cell-state markers and where genes are expressed using the
job's results. Differential expression between sample groups is not computed in this
workflow; Chat reports that limitation when asked for a group comparison.
