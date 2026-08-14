# How to use scALABLE-viewer

scALABLE-viewer explores a single-cell atlas that has already been computed. You do not
upload data and you do not launch a run. The dataset opens by itself.

scALABLE, the analysis tool, is a different program. scALABLE takes your raw data,
aligns it and runs the comparisons. scALABLE-viewer serves the finished result.

---

## The three tabs

**Study** describes the dataset: title, abstract, assay, organism, technology, cell
count, and links to the raw data and the paper. Below that sit three lists: **Tools**
links to other viewers for the same study, **Downloads** lists the released files, and
**Samples** lists every sample with sorting, filtering and paging.

**Explore** draws the cells. Two panels sit side by side so you can compare two views at
once.

**Differential** shows the precomputed group comparisons.

Explore opens first. Click Study or Differential to move.

---

## Explore

Each panel has its own **Select plot type** menu:

| Plot type | What it shows |
|---|---|
| UMAP cell types | Every cell, coloured by cell state, with state names on the map |
| UMAP broad | The same map at a coarser grouping |
| UMAP | One gene's expression across the map |
| Violin | One gene's distribution within each cell state |
| DotPlot | Dot size is the fraction of cells expressing a gene; colour is the mean level |
| Cell frequency | How cell-state proportions differ between groups |
| MarkerHeatmap | The marker genes that define each cell state |
| MarkerNetwork | Known interactions among one state's marker genes |
| Cell communication | Predicted receptor-ligand signalling between cell states |

**Select gene** takes a gene symbol. Type a few letters and pick from the list.

**Dot size** sits beside the plot type and controls point size. Each panel has its own.

**Filter data to display** narrows the plot. Choose a column under Annotation 1 or
Annotation 2, then choose which values to keep. Filters apply to that panel only.

**Download PDF** saves the current plot as an editable vector PDF. Text stays text, so
you can restyle it in Illustrator or Inkscape.

---

## Differential

Pick a **contrast** to choose which comparison to view, then a **cell state**.

| View | What it shows |
|---|---|
| Heatmap | Genes changed in that cell state, across samples |
| Volcano | Fold change against significance, one point per gene |
| Gene detail | One gene's values in both groups |
| GO Terms | Enriched biological processes |
| Network | Interactions among the changed genes |

Every view exports to PDF.

### Filter by gene

**Filter by gene** sits beside the cell-state menu. Type a gene symbol and press Enter, or
pick one from the list. The list offers the genes of the view you are looking at, so every
gene it offers gives you a result.

One gene selects a gene set: the gene you typed, plus every gene it interacts with in that
cell state. All four views then show that set only.

| View | What the filter keeps |
|---|---|
| Network | The gene, the genes it connects to, and the edges between them |
| Heatmap | One row per gene in the set |
| Volcano | One point per gene in the set, each point labelled |
| GO Terms | Terms that overlap the set |
| Cell communication | Interactions that use the gene as ligand or as receptor |

A cell communication comparison has no gene-to-gene network behind it, so the filter there
matches the gene alone. Its network keeps the interactions that use the gene and the cell
states those interactions join. Its table keeps the matching rows.

scALABLE names the feature after the modality, and the box follows: `Filter by gene` for
RNA, `Filter by protein` for ADT, `Filter by TF` for a GRN. Changing the modality clears
the box, because a new modality names different features.

A line under the menus reports the result, for example
`Filtered to FGF2 + 10 interacting genes: 60 of 2191 GO terms shown.`

The box turns blue when the filter matches and red when it matches nothing. Clear the box
to see everything again.

Two limits. A gene with no interaction edge in that cell state gives a set of one gene, and
the line under the menus says so. The cell communication network draws cell states, not
genes, so the filter leaves it whole and says why.

**Download PDF** follows the filter. With a filter on, the button saves the plot you see
rather than the whole comparison.

### Reading the GO Terms plot

Each point is one biological process. The x axis is the enrichment Z-score. The y axis is
the false discovery rate, drawn so the strongest results sit lowest.

Colour carries two different meanings, so read the legend:

- **Representative** terms are the ones GO-Elite kept after it removed redundant parents
  and children. Use these first.
- **Significant** terms passed the FDR cut but were folded into a broader term, or rest on
  few genes.
- **Not significant** terms failed the FDR cut.

A term can be highly significant without being representative. GO-Elite prunes the tree so
a hundred near-identical terms do not crowd out the real signal.

---

## What the numbers mean

A cell in this atlas may be a **metacell**. A metacell sums several real cells that were
matched on their donor metadata. Metacells protect participant privacy and reduce noise.
The Study tab states the metacell count and how many donors each one mixes.

Differential testing runs on **pseudobulks**. The viewer sums each sample's cells within
one cell state, then compares those sums between groups. Each point in a comparison is
therefore a sample, not a cell, which is what keeps the statistics honest.

---

## Things that surprise people

**A cell state is missing from Differential.** A comparison needs enough samples on both
sides. A state that appears in only a few samples is skipped rather than reported at low
confidence. The comparison table lists how many samples each side had.

**MarkerHeatmap takes a moment.** It loads a large matrix into an embedded viewer. Give it
a few seconds.

**A gene shows nothing.** The gene may not be detected in this dataset. Check the spelling
against the suggestion list.

**Plot types disappear from the menu.** Reload the page. Report it if it repeats.

---

## Getting a figure out

Use **Download PDF** rather than a screenshot. The PDF holds real vector shapes and
editable text at any size. Screenshots are fixed-resolution images and cannot be edited.

---

## Where the data came from

The Study tab names the source study, the raw-data accession and the publication.
`README.md` beside this file describes how the viewer is built and served.
