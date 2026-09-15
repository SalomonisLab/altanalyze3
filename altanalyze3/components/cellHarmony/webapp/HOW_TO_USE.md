# How to use scALABLE

scALABLE aligns your single-cell RNA data to a reference atlas and lets you inspect the result,
compare groups, and ask questions in plain language. The interface has four tabs: `Run`,
`Explore`, `Differential` and `Chat`. Run opens first.

This guide describes the controls as they exist on 2026-09-14. The methods behind each number,
the thresholds, and the output files are in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/README.md`.

A walkthrough video of an earlier version covers upload, QC, ambient correction, plots and
differential analysis: https://vimeo.com/1179728118/c5e0c39d24

## What you can upload

| Rule | Value |
| --- | --- |
| File types | `.h5` (Cell Ranger) or `.h5ad` (AnnData) |
| Files per job | up to 7, one per sample |
| Size per upload request | up to 1 GiB |
| Species and reference | one per job |

The pipeline keeps compatible `obs` columns of an uploaded `.h5ad`. Explore filters,
Differential cell-state selection and Differential group selection reuse them.

| Upload | Group differential available |
| --- | --- |
| two or more `.h5` files | yes, samples form the groups |
| one `.h5ad` with two or more group values in `obs` | yes, an `obs` column forms the groups |
| one `.h5` file | no |

## Run tab

### 1. Upload

1. Choose `Species`, then `Reference`. The preview draws the reference UMAP with its cell-state labels.
2. Click `Add sample` once per file, and give each row a `Sample name` and a file.
3. Click `Upload`. The app copies the files into a new job folder and shows the job id.

### 2. QC and alignment

| Control | Default | Meaning |
| --- | --- | --- |
| `Min genes` | 500 | drop a cell with fewer detected genes |
| `Min counts` | 1000 | drop a cell with fewer UMI counts |
| `Min cells` | 0 | drop a gene detected in fewer cells |
| `Mito %` | 15 | drop a cell above this mitochondrial fraction |
| `Minimum cosine similarity score` | 0.4 | drop a cell whose best reference match scores below this |
| `Ambient RNA correction` | No | `Yes` estimates and subtracts ambient RNA per sample before alignment |
| `Impute modality` | None | see below; hidden when the reference declares no imputation |

Ambient RNA correction runs per uploaded sample. Inspect your data for ambient contamination
before choosing `Yes`. In many droplet datasets the estimated fraction sits near 20 percent.

`Impute modality` predicts a second data type from the aligned RNA with a model that ships with
AltAnalyze3. The menu shows only the options the chosen reference supports:

| Option | What the job gains | References |
| --- | --- | --- |
| `None` | nothing | all |
| `All available` | every modality below that the reference supports | all |
| `ADT (CITE-seq)` | surface-protein abundance per cell: 56 antibodies for human lung, 129 for human bone marrow, 103 for mouse bone marrow | human lung, human bone marrow, mouse bone marrow |
| `Lipids` | 202 lipid species per cell | human lung |
| `Metabolite (AML)` | 2,533 metabolites per sample and cell state | human bone marrow |
| `Lipid (AML)` | 1,009 lipids per sample and cell state | human bone marrow |
| `GRN (TF activity)` | one activity score per transcription factor per cell, plus TF-to-target edge scores per sample and cell state | human lung, human bone marrow |

The `Metabolite (AML)`, `Lipid (AML)` and GRN edge predictions are per pseudobulk, one value
per sample and cell state, because their models learned from bulk profiles. Every predicted
modality is a model output, not a measurement. Each model's held-out accuracy is in its own
`README.md` under `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/`
(`rna2adt`, `rna2lipid`, `rna2metabolite`, `rna2grn`).

Click `Save QC and run`.

### What runs

1. QC filters the cells and normalises the counts.
2. cellHarmony aligns each cell to the reference by cosine similarity and drops cells below the cutoff.
3. MarkerFinder finds the top 50 markers per cell state and NetPerspective draws a network per state.
4. Approximate UMAP places each cell on the reference map.
5. The selected imputation models run.
6. fastComm scores receptor-ligand communication between cell states.
7. The app writes the combined h5ad and the download files.

The Run tab shows a progress bar, the QC cell counts, and the pipeline log while the job runs.
When it finishes, the app switches to `Explore`.

## Explore tab

Explore shows two panels side by side. Each panel has its own `Select plot type`, `Modality`,
`Select gene`, `Dot size`, `Filter data to display` and `Download PDF`. Both panels read the
same job, so you can compare two views of the same cells.

### Plot types

| Plot type | What it shows | Needs |
| --- | --- | --- |
| `UMAP cell types` | every query cell on the reference map, coloured by cell state, reference cells in grey | |
| `UMAP broad` | query cells in orange over the reference in grey | |
| `Cell frequency` | for each sample, the fraction of its cells in each cell state, as stacked bars | |
| `UMAP` | one feature's value per cell on the map; cells with value 0 stay in the background | a feature |
| `Violin` | one feature's distribution per cell state, every cell as a dot, top states by mean | a feature |
| `DotPlot` | for a gene set and a grouping: dot size is the fraction of cells above zero, colour is the mean | a gene set |
| `CombPlot` | for a gene set: one bar per cell state and donor, the donor's mean value in that state | a gene set and a donor column |
| `MarkerHeatmap` | the MarkerFinder matrix in an embedded Morpheus heatmap | network access |
| `MarkerNetwork` | known interactions among one state's marker genes | a marker cell state |
| `Cell communication` | receptor-ligand signalling between cell states, seven plot types | a completed fastComm run |
| `GRN edges` | TF-to-target edges touching the chosen genes, scored per sample and cell state | the `grn` modality |

The list changes with the modality: `MarkerNetwork` and `Cell communication` appear for RNA
only, and `GRN edges` appears for GRN only.

### Modality

The `Modality` menu appears when the job imputed at least one modality. It changes what
`Select gene` offers: gene symbols for RNA, antibody names such as `Hu.CD4` for ADT, lipid or
metabolite names, or transcription factor names for GRN. Imputed modalities use a blue to
yellow to red colour ramp; RNA uses grey to red.

### UMAP options

For `UMAP cell types` two extra menus appear. `Color by` colours the cells by any `obs` column
with 2 to 60 levels. The `X` and `Y` menus choose the axes: the cellHarmony UMAP, any other
embedding stored in the h5ad, or any two float-valued `obs` columns such as a QC score. The app
hides the reference cells when the colour or the axes leave their defaults, because the
reference has no such column. A cell without a value on a chosen axis is not drawn, and the
panel reports how many it dropped.

### DotPlot and CombPlot controls

| Control | Meaning |
| --- | --- |
| `Gene set` | one or more feature names, separated by spaces, commas or semicolons; empty gives one marker gene per group |
| `Group by` | the `obs` column that forms the columns; default is the cell state |
| `Show` | the levels to keep; nothing selected keeps every level |
| `Min cells` (CombPlot) | 1, 5, 10 or 25; a donor with fewer cells in a group is left out |

In the DotPlot, dot size runs from 4 px at 0 percent to 22 px at 100 percent of cells above
zero, and colour runs from white to red with the mean. In the CombPlot, a colour strip names
the cell state of each column, one row per gene holds the bars, and each bar is the mean over
one donor's cells in that state. The title states how many donor groups the plot shows.

### Cell communication plot types

When the plot type is `Cell communication`, the `Modality` menu becomes `Plot type`, and the
`Marker cell state` menu chooses the focus state.

| Plot type | What it shows |
| --- | --- |
| `Focused incoming` | senders signalling to the focus state, up to 60 state pairs |
| `Focused outgoing` | receivers the focus state signals to |
| `Cell-state network` | every state pair, no focus |
| `Ligand-receptor dot plot` | one point per interaction of the focus state |
| `Cell-state heatmap` | sender by receiver matrix of summed scores |
| `Top interactions table` | the focus state's incoming interactions ranked by score |
| `Per-sample comparison` | per-sample totals for the top 12 partner states |

Edge width and opacity follow the summed score of the pair. The plots apply no threshold of
their own; the run's thresholds are in the README.

### GRN edges

For the `grn` modality, `GRN edges` draws every TF-to-target edge that touches a gene in the
gene set, averaged over the chosen sample and cell state. A threshold on the absolute score and a
cap of 300 edges limit the drawing.

### MarkerHeatmap and MarkerNetwork

`MarkerHeatmap` opens the MarkerFinder fold matrix in Morpheus, which loads from the Broad
Institute's server. When the job holds two MarkerFinder runs, a density menu offers
`10 metacells per cell population (random)` or `All metacells`; the two runs pick different
marker rows, so switching changes the genes shown. `MarkerNetwork` needs a `Marker cell state`.

### Filter data to display

Choose a column under `Annotation 1`, then one value, and optionally the same under
`Annotation 2`. Only cells matching both values are drawn. The filter applies to that panel
only, changes no saved output, and reruns nothing. For `Cell communication`, a value on the
sample column selects that sample's scores, and a value on the cell-state column keeps
interactions where that state sends or receives.

### Download PDF

`Download PDF` saves the current plot as a vector PDF with editable text. The MarkerHeatmap
button saves the pipeline's own heatmap PDF. A `Cell communication` network saves as SVG.

### Downloads

| File | Content |
| --- | --- |
| assignments | one row per cell: cell state, score, UMAP coordinates |
| combined h5ad | the aligned dataset with approximate UMAP |
| marker genes ZIP | MarkerFinder tables, heatmap PDF and marker networks |
| `<modality>_results.zip` | the imputed h5ad and its marker outputs, one per imputed modality |

## Differential tab

The tab works when the job has two or more `.h5` samples, or one `.h5ad` with a group column.
Otherwise it shows a message instead of the controls.

### Setup

| Control | Meaning |
| --- | --- |
| `Cell-state aligned to` | the `obs` column that defines the populations to test within |
| `Modality` | which feature matrix to test; shown when the job has more than one |
| `Group values from` | the `obs` column that holds the groups |
| `Comparison Type` | `cells` tests cells; `pseudobulk` sums each sample's cells per state and tests samples |
| `Group 1 (numerator)` | the case values |
| `Group 2 (denominator)` | the control values |

`pseudobulk` appears only for one `.h5ad` or four or more files. It needs at least 2 samples
per group within a cell state; a state with fewer is not tested and does not appear in the
views. The test and the fold and p thresholds depend on the modality; the README section
"The differential engine" lists them. The GRN modality tests TF-to-target edges, not TFs.
`Cell communication` appears as a modality once a fastComm run exists and tests interactions
per receiver state.

### Views

| View | What it shows |
| --- | --- |
| `Summary` | up- and down-regulated feature counts per cell state as diverging bars |
| `Heatmap` | the features called in the chosen state, across every state, as folds |
| `Volcano` | fold change against significance, one point per feature |
| `GO Terms` | enriched Gene Ontology terms, RNA only |
| `Network` | known interactions among the changed genes, RNA only |
| `Table` | the top interactions, cell communication only |
| `Gene Detail` | one feature's values in both groups, with its p, FDR and fold |

Click a feature in the left view to fill `Gene Detail` on the right. Every view exports to PDF
through the two `Download PDF` buttons, and the differential ZIP holds every table.

### Filter by gene

`Filter by gene` sits beside the cell-state menu. Type a feature name and press Enter, or pick
one from the list. One gene selects a set: the gene plus every gene it interacts with in that
cell state's network. Heatmap, Volcano, GO Terms and Network then show that set only. The box
follows the modality: `Filter by protein` for ADT, `Filter by TF` for GRN. A line under the
menus reports the result, for example `Filtered to FGF2 + 10 interacting genes: 60 of 2191 GO
terms shown.` With a filter on, `Download PDF` saves the filtered figure you see.

### Reading the GO Terms plot

Each point is one term. The x axis is the enrichment z-score. The y axis is the false discovery
rate, drawn so the strongest results sit lowest. GO-Elite marks a term `Representative` when it
passed `z >= 1.96`, `FDR <= 0.1`, at least 3 overlapping genes, and no parent or child term
represents it better. Other terms below the FDR cut are `Significant`. Use the representative
terms first; a term can pass the FDR cut and still fold into a broader one.

## Chat tab

Type a question and click `Ask`, or click one of the example questions the app builds from your
job. The app sends the question to an assistant service that reads it into one of 17 named
analyses, then computes the answer from your job's own data. The assistant never sees a value
and never states a number.

| Ask about | Example | Answer |
| --- | --- | --- |
| the markers of a cell state | What are the best marker genes of AT2 cells? | 25 markers and a DotPlot |
| two cell states | What distinguishes AT1 from AT2 cells? | 30 markers and a DotPlot |
| where a gene is expressed | Where is SFTPC expressed? | the top 5 states by mean and a DotPlot |
| a completed comparison in a cell state | Which genes are significant in COPD versus control in AT2 cells? | the top 25 rows and a volcano |
| pathways, regulators, cell communication | Which pathways change in AT2 cells? | the tab that holds that analysis |

Below the answer, a `Table` and `Plot` toggle appears. The plot reuses the Explore DotPlot or
the volcano, so the chat cannot draw a figure Explore cannot.

| Reply | Meaning |
| --- | --- |
| a list of states or comparisons | the question named none of them; pick one |
| `not implemented` | the analysis exists in the viewer program but not here; the reply names the missing statistic |
| `not covered` | the comparison exists but did not test that cell state; the reply lists the states it did |
| `not run` | no comparison has completed; run one in `Differential` first |
| HTTP 503 | the assistant service on port 8001 is down |

## Reset

After an upload, `Reset data` in the header returns the app to the pre-upload state without a
browser reload.

## Where the outputs are

Each job lives in one folder under
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/jobs/<job id>/`
with `uploads/`, `outputs/` and `logs/pipeline.log`. The README lists every output file. The
server deletes finished jobs older than 8 hours when a new upload arrives.
