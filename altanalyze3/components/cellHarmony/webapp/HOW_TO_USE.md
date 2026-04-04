# How to Use cellHarmony web

This guide is written for a first-time user of the `cellHarmony web` site. It
explains what files to upload, what each step does, and what the viewers are
showing at each stage of the workflow.

## What this website does

`cellHarmony web` supports two main analysis stages:

1. **Alignment to a reference atlas**
   - your uploaded single-cell data are quality controlled
   - the cells are aligned to a reference atlas cell state
   - the site then provides UMAP- and gene-expression-based viewers

2. **Differential analysis between sample groups**
   - available after alignment completes for jobs with at least 2 samples
   - compares two groups of samples within aligned cell states
   - provides an interactive heatmap, volcano plot, network, GO terms, and a
     gene-detail violin plot

## What files you upload

Upload one file per sample.

Accepted file types:

- `.h5`
- `.h5ad`

Practical expectations:

- each uploaded file should represent one sample
- each sample should contain single-cell expression data
- if you upload multiple files, they should be biologically comparable and
  intended to be processed against the same selected reference
- when you upload `.h5ad` files, compatible `.obs` metadata are retained and
  can be reused later in the interface for filtering displayed cells and for
  differential grouping

When uploading, you also provide:

- a **sample name**
- the **species**
- the **reference atlas**

The sample name is what later appears in the differential grouping interface.

## 1. Open the app

Open the web app in your browser.

Common locations:

- local run: `http://127.0.0.1:8000/`
- VM or server run: `http://<server>:8000/`
- subpath deployment: for example `https://<server>/cell-harmony/`

## 2. Upload the data

In **1. Upload**:

1. Select the **Species**.
2. Select the **Reference**.
3. Click **Add sample** for each file you want to include.
4. For each sample row:
   - enter the sample name
   - choose the `.h5` or `.h5ad` file
5. Click **Upload**.

What happens here:

- the files are copied into a new job folder on the server
- the upload progress bar shows file transfer progress
- once upload completes, a job ID is assigned
- **Panel 2** appears only after upload succeeds

At this point, no QC or alignment has been run yet. The upload step only
creates the job and stores the input files.

## 3. Run QC and alignment

In **2. QC and alignment**:

Review the QC settings:

- `Min genes`  
  default: `500`
- `Min counts`  
  default: `1000`
- `Min cells`
- `Mito %`
- `Minimum cosine similarity score`
- `Identify cell-state marker genes`
  default: `TRUE`

Then click **Save QC and run**.

What happens during this step:

1. the uploaded data are loaded
2. cells are filtered by QC thresholds
3. the filtered cells are aligned to the selected reference atlas
4. if marker analysis is enabled, cell-state marker genes are identified from the aligned dataset
5. if marker analysis is enabled, a marker heatmap and marker-driven NetPerspective networks are exported
6. approximate UMAP placement is generated relative to the reference
7. alignment outputs are written to the job folder

While it runs, Panel 2 shows:

- a **progress bar**
- a **QC cell summary bar**
- a **log window**

What those running summaries mean:

- the QC summary bar reports the number of cells:
  - before QC
  - after `min_genes`
  - after `min_counts`
  - after `mito %`
- the log window shows the live pipeline output

When alignment completes:

- **Panel 3** appears
- the lower baseline viewers become available
- if there are at least 2 samples, **Panel 4** can be used for differential analysis

### What marker analysis adds

When **Identify cell-state marker genes** is `TRUE`:

- the web pipeline identifies the top unique marker genes per aligned cell state
- a heatmap PDF and TSV outputs are exported
- NetPerspective networks are exported for each aligned cell state
- a ZIP archive of those marker outputs becomes available in **Panel 3**

This marker analysis runs after alignment and before the approximate UMAP outputs are finalized.

## 4. Review the alignment outputs

After alignment finishes, **3. Results** becomes active.

This section lets you inspect a selected gene in the aligned dataset.

Controls:

- **Gene**  
  type a gene symbol of interest
- **UMAP display**
  - `Relative to reference`
  - `Cell-state clusters`
- **Expression display**
  - `UMAP`
  - `Violin`
  - `MarkerHeatmap` when marker analysis outputs exist
  - `MarkerNetwork` when marker networks exist
- **Dot size**
  - controls the point size used in the UMAP and violin viewers
- **Filter data to display**
  - **Annotation 1** and **Annotation 2** let you choose compatible `.obs`
    metadata fields from the aligned `.h5ad`
  - **Values** lets you restrict the viewers to one selected value or `All`

### What the UMAP viewer is showing

The lower-left **Approximate UMAP** panel is based on the aligned query data
after QC and reference mapping.

Two modes are available:

- **Relative to reference**
  - shows the uploaded cells placed in the reference embedding
  - useful for seeing how the uploaded dataset maps onto the atlas structure

- **Cell-state clusters**
  - colors uploaded cells by assigned aligned cell state
  - reference-only cells remain in the background as a light grey context layer
  - can be restricted to selected subsets using **Filter data to display**

### What the Expression viewer is showing

The lower-right **Expression** panel uses the currently selected gene and the
aligned dataset after QC.

Two modes are available:

- **UMAP**
  - colors cells by normalized expression of the selected gene
  - useful for seeing where expression occurs in the aligned manifold
  - can be restricted to selected subsets using **Filter data to display**

- **Violin**
  - shows the distribution of normalized expression for that gene
  - useful for a quick expression summary across the aligned populations

Additional marker-analysis modes may also be available:

- **MarkerHeatmap**
  - opens the exported MarkerFinder heatmap in an embedded Morpheus viewer
  - uses the marker heatmap matrix produced during the alignment job
  - now respects the current **Filter data to display** barcode restrictions in the interactive viewer

- **MarkerNetwork**
  - displays the exported NetPerspective network for a selected marker-defined cell state
  - uses the **Marker cell state** dropdown that appears when marker networks are available

### How the cell-display filters work

The **Filter data to display** controls affect only what is shown in:

- the **Approximate UMAP**
- the **Expression** viewer
- the interactive **MarkerHeatmap** viewer
- the corresponding PDF downloads

They do not rerun alignment or change the saved job outputs.

Typical uses:

- display only one uploaded sample
- display only one biological subset such as a chosen cell annotation
- combine a sample-level filter with a second annotation-level filter

By default:

- **Annotation 1** starts from the sample-related field used later in
  differential grouping when available
- **Annotation 2** starts from the aligned cell-state annotation field when
  available

Only manageable categorical `.obs` fields are shown. UMAP columns and very
high-cardinality metadata are excluded from these dropdowns.

### What downloads are available in Panel 3

Depending on the job state, downloads may include:

- assignments
- combined h5ad
- marker genes ZIP
- viewer PDFs

Notes:

- **Download assignments** provides the combined assignments table with appended UMAP coordinates
- **Download marker genes ZIP** appears only when marker analysis was enabled and completed
- these files come from the completed alignment stage, not the differential stage

## 5. Run cellHarmony-differential

**4. Differential** becomes available only when:

- alignment has completed
- the job contains at least 2 samples

In this section:

1. choose the aligned cell-state field in **Cell-state aligned to**
2. choose the grouping source in **Group values from**
3. if available, choose **Comparison Type**
   - `cells`
   - `pseudobulk`
4. assign values to:
   - **Group 1 (numerator)**
   - **Group 2 (denominator)**
5. click **Run cellHarmony-differential**

Important grouping rule:

- only two-group comparisons are supported
- each group can contain one or more selected values from the chosen grouping
  field

This step is flexible for `.h5ad` uploads:

- **Cell-state aligned to** can use the reference-aligned field or another
  compatible annotation stored in `.obs`
- **Group values from** can be a different `.obs` field than the cell-state
  field
- this lets you compare values such as `condition`, `sample`, `Library`, or
  another categorical annotation while still running DE within the chosen cell
  states

What happens during differential analysis:

1. the aligned AnnData from the completed cellHarmony run is used
2. cells are grouped by the selected aligned cell-state field
3. Group 1 and Group 2 are defined from the selected **Group values from**
   field
4. Group 1 and Group 2 are compared within each selected cell state
5. differential statistics, heatmaps, GO results, and networks are generated

If **Comparison Type** is set to `pseudobulk`:

- pseudobulk profiles are computed before DE testing
- raw p-values are used in that workflow
- internal minimum-cell filtering is applied to pseudobulk groups

While this runs:

- the differential progress bar updates
- the baseline Panel 3 controls are greyed out because the interface switches
  to differential exploration mode

When it completes:

- the interactive differential explorer appears below
- a ZIP archive is available for download

## 6. Explore the differential results

The lower differential section has two panels:

- a **large left panel** for the selected visualization
- a **smaller right panel** for gene-level detail

At the top of the left panel:

- choose the visualization mode:
  - `Heatmap`
  - `Volcano`
  - `Network`
  - `GO Terms`
- choose the cell state relevant to that mode

The cell-state dropdown is restricted to only the populations that exist for the
selected visualization type.

### Heatmap view

What it is showing:

- genes selected from the detailed differential output for the chosen cell state
- fold-change patterns for those genes across aligned populations
- the color scale summarizes the differential signal for those genes

How to use it:

- switch to **Heatmap**
- choose a cell state from the dropdown
- inspect the gene rows
- click a gene row to update the right-side gene-detail panel

Use this when you want to see:

- how the selected cell state’s DE genes behave across the broader aligned atlas
- how genes from a user-defined `.h5ad` cell annotation behave across the
  populations used in the differential run

### Volcano view

What it is showing:

- one point per gene for the selected cell state
- x-axis: `log2 fold change`
- y-axis: `-log10(FDR)`

How to use it:

- switch to **Volcano**
- choose the cell state
- click a point to select a gene

Use this when you want to see:

- which genes are strongly shifted in Group 1 vs Group 2
- both effect size and significance in the same view

### Network view

What it is showing:

- an interactive network for the selected cell state
- each node is a gene
- edges represent the network relationships produced by the differential workflow

How to use it:

- switch to **Network**
- choose the cell state
- zoom, pan, and inspect the graph
- click a node to send that gene to the right-side gene-detail panel

Use this when you want to see:

- gene relationships rather than only ranked DEGs

### GO Terms view

What it is showing:

- the top GO-enrichment terms for the selected cell state
- only the top 15 terms are shown
- the hover contains the overlap genes associated with that GO term

How to use it:

- switch to **GO Terms**
- choose the cell state
- hover over a bar to inspect the overlap genes
- click the term to choose a gene from that term for the right-side panel

Use this when you want to see:

- which biological processes are enriched in the selected aligned cell state

## 7. Use the Gene Detail panel

The right-side **Gene Detail** panel is driven by the current gene selected from
the left-side visualization.

It can be updated from:

- a heatmap gene row
- a volcano point
- a network node
- a GO term gene choice

What it is showing:

- a normalized expression **violin plot**
- individual cell dots overlaid on the violins
- Group 1 vs Group 2 for the selected gene in the selected cell state

Below the violin plot, the statistics panel shows:

- `P-value`
- `FDR`
- `log2FC`

This panel answers:

- whether the selected gene is more highly expressed in Group 1 or Group 2
- how strong that change is
- how statistically significant it is in that cell state

## 8. Download the outputs

Downloads are available from multiple places.

Alignment-stage downloads:

- assignments
- combined h5ad
- UMAP outputs
- baseline viewer PDFs

Differential-stage downloads:

- differential PDF exports from the active left-side view
- the differential ZIP archive

The differential ZIP excludes `.h5ad` files from the packaged archive.

## 9. Reopen an existing job

If the server still retains the job files, the job can be revisited without
re-uploading the original data.

Typical use:

- open the web app
- load or revisit the existing job state on the same server
- continue reviewing the alignment or differential outputs

## 10. Troubleshooting

### The page looks stale or controls do not update

- do a hard refresh in the browser

### The upload completed but later panels do not appear

- wait for the server to finish updating job status
- then refresh the job state in the interface

### Differential results are empty for a view

- some visualizations only exist for specific cell states
- switch the visualization type or choose another applicable cell state

### The wrong cells or groups seem to be selected

Check that the intended fields were chosen in:

- **Filter data to display**
- **Cell-state aligned to**
- **Group values from**

If you uploaded `.h5ad` files with multiple annotation fields, the web app can
use different `.obs` columns for display filtering and for differential
analysis.

### Docker deployment changed but the app still behaves the old way

Rebuild the container:

```bash
docker compose up --build -d
```

### The app is served under a subpath such as `/cell-harmony`

Make sure `CELLHARMONY_ROOT_PATH` is set consistently for that deployment.
