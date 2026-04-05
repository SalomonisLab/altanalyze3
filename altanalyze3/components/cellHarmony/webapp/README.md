# cellHarmony web

`cellHarmony web` aligns single-cell datasets to a reference atlas, provides an
interactive Explore workspace for UMAP and gene-expression review, and supports
group-based differential analysis in a dedicated Differential workspace.

## Run locally

From the repo root:

```bash
uvicorn altanalyze3.components.cellHarmony.webapp.app:app --reload
```

Open:

```text
http://127.0.0.1:8000
```

## Interface overview

The promoted interface is organized into three top-level tabs:

- `Run`
  - upload files
  - configure QC and alignment
  - review the reference preview before upload
- `Explore`
  - inspect aligned UMAPs
  - inspect gene expression by UMAP or violin plot
  - review marker heatmaps and marker networks when marker analysis is enabled
  - download assignments, the combined h5ad, and marker ZIP outputs
- `Differential`
  - compare biological groups after alignment completes
  - inspect heatmap, volcano, network, GO terms, and gene detail views

## Accepted inputs

Upload one file per sample.

Accepted formats:

- `.h5`
- `.h5ad`

Behavior by input type:

- multiple `.h5` files support group-based differential analysis
- a single `.h5ad` can also support group-based differential analysis when the
  `.obs` metadata contain multiple biological groups
- a single `.h5` upload does not enable group differential analysis

Compatible `.obs` metadata from uploaded `.h5ad` files are preserved and reused
for:

- Explore tab filtering
- Differential cell-state selection
- Differential biological-group selection

## Alignment workflow

In `Run`:

1. choose `Species`
2. choose `Reference`
3. add one or more samples
4. upload files
5. review QC settings
6. click `Save QC and run`

QC/alignment settings include:

- `Min genes`
- `Min counts`
- `Min cells`
- `Mito %`
- `Minimum cosine similarity score`
- `Identify cell-state marker genes`

When alignment completes:

- the app switches to `Explore`
- aligned results become available immediately
- if the job supports grouped comparisons, the `Differential` tab becomes usable

## Marker analysis

When `Identify cell-state marker genes` is `TRUE`:

- markerFinder is run on the aligned dataset
- a marker heatmap PDF and TSV outputs are exported
- redundant top-ranked marker tables are exported for network generation
- NetPerspective marker networks are exported per cell state
- a ZIP archive becomes available in Explore downloads

Explore then exposes two additional expression modes:

- `MarkerHeatmap`
- `MarkerNetwork`

## Explore workspace

The Explore tab keeps aligned results available even after differential
analysis has been run.

### Approximate UMAP

The left viewer supports:

- `UMAP broad`
- `UMAP cell types`
- `Cell frequency`

Notes:

- `Cell frequency` reports normalized fractions per sample using the currently
  filtered cells
- cluster colors are kept consistent with the reference preview color mapping

### Expression

The right viewer supports:

- `UMAP`
- `Violin`
- `MarkerHeatmap` when marker outputs exist
- `MarkerNetwork` when marker networks exist

The `Marker cell state` dropdown only appears for `MarkerNetwork`.

### Explore filters

`Filter data to display` restricts the visible cells used in:

- Approximate UMAP
- Expression UMAP
- Violin
- MarkerHeatmap
- Cell frequency

These filters affect visualization only. They do not rerun alignment.

## Differential workspace

The Differential tab is used for comparisons between biological groups within
cell states.

It is enabled only when:

- two or more samples were uploaded for the job, or
- one `.h5ad` contains multiple biological groups in reusable `.obs` metadata

Controls include:

- `Cell-state aligned to`
- `Group values from`
- `Comparison Type`
  - `cells`
  - `pseudobulk` when applicable
- `Group 1 (numerator)`
- `Group 2 (denominator)`

The interactive explorer supports:

- `Heatmap`
- `Volcano`
- `Network`
- `GO Terms`

Selecting genes from those results populates the `Gene Detail` viewer on the
right.

## Downloads

Explore downloads can include:

- assignments
- combined h5ad
- marker genes ZIP

Differential downloads can include:

- differential ZIP
- PDF exports for the active differential panel
- PDF export for the selected gene detail plot

## Data location

Job data are stored under:

- `cellHarmony/webapp/jobs`
