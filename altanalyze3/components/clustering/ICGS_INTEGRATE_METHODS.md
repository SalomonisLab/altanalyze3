# Cross-study integration of ICGS3 cell states

Module: `/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/clustering/ICGS_integrate.py`

Below I specify the algorithm the module runs, the reason behind each choice, the tunable
thresholds, and the measured result on four human lung datasets. A reader can reproduce the run
from the commands in section 16.

---

## 1. Problem and design choice

### 1.1 What the module does

ICGS3 clusters one dataset at a time. Four independent ICGS3 runs therefore produce four separate
cluster sets that use four separate numbering schemes. The module builds one non-redundant set of
cell states from those cluster sets. A state describes a population. The module decides which
clusters describe a population the reference already holds, and which describe a population no
other dataset resolved.

### 1.2 Why the module compares clusters and not cells

Embedding-based integration methods such as Harmony, scVI and Seurat anchors place all cells in one
latent space, then shrink the distances they attribute to batch. Such a method must decide which
part of a distance comes from the assay and which part comes from biology. No labelled example
distinguishes the two, so the method decides by assumption.

The module avoids that decision. ICGS3 has already reduced each dataset to clusters with marker
genes. A marker gene set is a discrete object. Two marker gene sets either overlap or do not
overlap, and a hypergeometric test measures the overlap against an explicit null. The module never
moves a cell and never rescales an expression value between datasets.

The cost of the choice: the module cannot split a cluster that one study under-resolved. Its
resolution ceiling equals the resolution of the input runs.

### 1.3 What the module never does

The module never merges two clusters into one state. A state always holds the cells of exactly one
source cluster from exactly one dataset. When cluster B of a later dataset proves redundant with
state A, the module excludes cluster B and leaves state A untouched, with its original cells and
its original markers. The module records dataset B in the `supporting_datasets` field of state A.

Recording support rather than merging keeps every state's expression profile free of a second
study's batch signal. An averaged centroid across studies would carry both.

---

## 2. Software and dependencies

| Component | Source | Use |
|---|---|---|
| `numpy`, `pandas` | standard | arrays and tables |
| `scipy.sparse` | standard | the pooled cell-by-gene matrix, CSR format |
| `scipy.stats.hypergeom` | standard | enrichment `sf(k-1, N, n, K)` |
| `anndata` | standard | reading each run's `icgs3_result.h5ad` |
| `marker_finder_wrapper` | `altanalyze3/components/udon/markerFinder.py` | unique-marker assignment |
| `hopach` | `altanalyze3/components/clustering/hopach.py` | centroid ordering |
| `generate_marker_heatmap_from_adata` | `altanalyze3/components/visualization/marker_heatmap_h5ad.py` | heatmap PDF |

The module calls no R code and no batch-correction library.

Interpreter used for every run reported here:
`/opt/homebrew/opt/python@3.11/bin/python3.11`.

---

## 3. Required input files

The module reads only files a completed ICGS3 run already wrote. Each `--run NAME=PATH` directory
must contain the following five files.

| File, relative to the run directory | Columns the module reads | Purpose |
|---|---|---|
| `icgs3_cell_barcode_clusters.tsv` | `barcode`, `ICGS3_cluster`, `ICGS3_SVM_score` | cell-to-cluster map and ranking score |
| `MarkerFinder/icgs3_markers_all_correlations.tsv` | `marker`, `top_cluster`, `pearson_r` | unique-marker assignment, drives Gate 1 |
| `MarkerFinder/icgs3_marker_heatmap_fold_matrix.centroids.tsv` | genes by clusters | supplies the gene list of that run |
| `icgs3_result.h5ad` | `X`, `layers['counts']` | expression of the representative cells |
| `GO-Elite/icgs3_cell_state_predictions.tsv` | `cluster`, `cell_type_prediction`, `fdr` | state naming, optional |

The module raises `FileNotFoundError` when any of the first three files is absent. The module
falls back to `MarkerFinder/icgs3_markers.tsv` when the all-correlations file is absent, and warns
when an `icgs3_result.h5ad` is absent.

---

## 4. Step 1 — Select representative cells

ICGS3 finishes by training a linear support vector machine on its clusters and scoring every cell.
The score measures the distance from the decision boundary. A high score means the cell sits far
inside its own cluster and far from every other cluster centroid.

The module sorts the cells of each cluster by `ICGS3_SVM_score` in descending order and keeps the
top `--cells-per-cluster` (default 200). A cluster holding fewer than 200 cells contributes all of
its cells, and the module logs how many clusters fall short.

**Justification.** The comparison asks whether two clusters describe the same population. Cells that sit
near a decision boundary answer that question poorly, because their own dataset already treats them
as ambiguous. Taking the highest-scoring cells sharpens each cluster's profile before any
cross-dataset comparison begins.

**Measured on the four lung runs.** 39,396 cells from 214 clusters. Clusters below 200 cells:
Adams 9 of 43 (smallest 42), Basil-2022 8 of 55 (smallest 90), PedDev 20 of 63 (smallest 41),
Natri-2024 13 of 53 (smallest 19).

**Effect of changing it.** Raising the value gives each cluster more cells. More cells raise the
correlations MarkerFinder computes in Step 6, so both gates grow more permissive. Lowering the
value makes marker detection noisier and raises exclusion in Step 6. Eleven of 214 clusters already
hold fewer than 200 cells, so values above roughly 250 change little.

Output: `representative_cells.tsv`.

---

## 5. Step 2 — Build the union feature space

Each ICGS3 run selected its own informative genes and wrote only those genes to its centroid file.
The four runs selected 2,580 (Adams), 3,266 (Basil-2022), 3,748 (PedDev) and 3,103 (Natri-2024)
genes. The module takes the union of the four lists, giving 6,715 genes.

**Justification.** The intersection of the four lists would discard exactly the genes that
distinguish one dataset from another. A gene that one study found informative is not absent from
the other studies' cells. Only that study's gene-selection step omitted it. The module therefore
takes the union and reads the missing values from expression in Step 3.

Comparing two clusters on the intersection of their two parent runs' selections would also change
the gene universe for every pair, which would make the hypergeometric tests in Step 6
non-comparable across pairs. The union fixes one universe, N = 6,715, for every test.

---

## 6. Step 3 — Build the pooled expression matrix

The module opens each run's `icgs3_result.h5ad` once, selects the representative barcodes, and
writes their expression into a shared cells-by-genes CSR matrix over the 6,715 union genes. Genes
a run's h5ad does not carry stay at zero for that run's cells. The module reports the count per
run. The module copies `layers['counts']` in parallel, so the final heatmap can rescale from real
integer counts instead of a gene-subset log matrix.

**Justification for reading the h5ad rather than the centroid files.** A centroid file holds only
the genes that run selected. Reading expression instead means a gene selected by PedDev alone is
still measured in Adams cells, which is the whole point of the union space.

**Measured.** Adams contributed 8,113 cells with 6,695 of 6,715 genes present; Basil-2022 10,571
cells with 6,681; PedDev 11,091 cells with 6,704; Natri-2024 9,621 cells with 6,696. The pooled
matrix holds 39,396 cells by 6,715 genes at 14.6% non-zero density.

**Memory.** The CSR matrix stores only non-zero values. 39,396 by 6,715 at 14.6% density occupies
roughly 0.31 GB as float32 plus int32 indices. The module never densifies the full matrix. It
densifies only the row blocks that a single MarkerFinder call needs.

---

## 7. Step 4 — Recompute every cluster centroid

The module computes each cluster's centroid as the mean expression of its representative cells over
the 6,715 union genes. Every cluster in every dataset now has a profile over an identical feature
set.

**Justification.** The centroid files ICGS3 wrote span four different gene lists. A cosine or
correlation between two such centroids compares two different feature spaces and returns a number
without meaning. Recomputing from the pooled matrix removes that defect.

---

## 8. Step 5 — Order the datasets and seed the reference

The module sorts the datasets by cluster count, descending, and the first dataset initialises the
reference. Each of its clusters becomes one state. The module assigns state identifiers from a
single monotonic counter, so no identifier is ever reused.

**Measured.** Order: PedDev (63 clusters), Basil-2022 (55), Natri-2024 (53), Adams (43). PedDev
seeds 63 states.

**Justification.** The dataset that resolved the most clusters offers the finest partition of the
biology. Seeding from the finest partition makes every later cluster face the most demanding reference
available. Starting from a coarse reference would admit later clusters simply
because the reference was too coarse to contain them.

**Known consequence.** The seed dataset is never tested and never loses a cluster. PedDev keeps all
63 states. Adams enters last against 99 states and contributes 5 of 43. Section 15.3 measures what
a different order does.

### 8.1 Overriding the order

Two options override the default sort.

- `--seed-dataset NAME` names the dataset that initialises the reference. The remaining datasets
  follow in descending cluster count.
- `--dataset-order A,B,C,D` names the whole order and overrides `--seed-dataset`.

`resolve_dataset_order` validates both before the module reads any file, so a typo fails in one
second rather than after the module has built the pooled matrix. The validator raises `ValueError`
on a name absent from `--run`, on an order that omits a dataset, and on an order that repeats one.
`integrate_hierarchically` calls the same validator, because a caller can invoke it directly.

The module logs a warning when the chosen seed holds fewer clusters than a later dataset, because
a coarser seed admits later clusters that a finer seed would have absorbed.

`integration_summary.json` records `seed`, `dataset_order`, `seed_dataset_override` and
`dataset_order_override`. The module reads the resolved order off the audit table, so the record
reflects the run rather than a recomputed guess.

---

## 9. Step 6 — Gate 1, exclusion by marker gene set enrichment

Gate 1 asks one question of each incoming cluster: does the current reference already hold a state
with this cluster's marker genes?

### 9.1 Building the reference marker databases

Every state holds exactly one source cluster, so a state's marker set equals that cluster's marker
set. The module reads those markers from the source run's
`MarkerFinder/icgs3_markers_all_correlations.tsv`, restricted to the 6,715 union genes, and builds
two databases per state:

- **Database A**: genes with Pearson r > 0.25, sorted descending, capped at 200 genes.
- **Database B**: genes with Pearson r > 0.30, sorted descending, capped at 100 genes.

**Why two databases.** Database A is permissive and detects a broad resemblance. Database B is
strict and detects agreement among the strongest markers only. A cluster passing either test counts
as redundant, and the module records the overlap and FDR from both in
`nomination_decisions.tsv`.

**Measured database sizes.** Median 39 genes per state when Basil-2022 entered, 31 when Natri-2024
entered, 28 when Adams entered. The caps of 200 and 100 therefore bind rarely. The medians fall
across steps because later-added states come from clusters that already survived one gate.

### 9.2 Why the databases use ICGS3's unique-marker assignment

ICGS3's all-correlations file assigns every gene to exactly one cluster, the cluster whose indicator
vector the gene correlates with best. Genes expressed across a whole lineage therefore go to one
cluster of that lineage and to no other. ICGS3 therefore competes lineage-wide genes away inside each study
before any cross-study comparison starts.

An earlier version of this module built the databases differently. It correlated every gene against
every state indicator across the pooled matrix and kept the top 200 per state above r = 0.25. That
version gave each gene a value for every state, so lineage-wide genes entered many databases at
once. Measured consequence: the median cluster shared 30 of its top 50 genes with some reference
state, Gate 1 called 52 of 55 Basil-2022 clusters redundant, and the reference reached only 67
states. I deleted the functions implementing that version on 2026-08-21.

### 9.3 The enrichment test

For each incoming cluster the module takes its top `--nomination-query-top` (default 60) unique
markers by Pearson r, restricted to the union genes, giving a query set Q. For each state S in a
database D the module computes the overlap k = |Q ∩ D(S)| and the hypergeometric survival
probability:

```
p(S) = hypergeom.sf(k - 1, N, |Q|, |D(S)|)          N = 6715
```

`sf(k-1, ...)` gives P(X ≥ k), the probability of seeing at least k shared genes when Q is drawn at
random from the 6,715 genes. The module skips states with k = 0. The module then adjusts every p across
the states compared, by Benjamini-Hochberg, using `_bh_adjust`:

```
sort p ascending; q_i = p_i * n / i; enforce monotonicity from the largest rank down; clip at 1
```

### 9.4 The decision rule

A cluster counts as redundant when both conditions hold in database A or both hold in database B:

**Condition 1, overlap size.**

```
k >= max(--nomination-min-overlap, ceil(--nomination-overlap-fraction * min(|Q|, |D(S)|)))
        default 10                  default 0.17
```

The fixed floor of 10 states the requirement for a normal case: 10 shared genes out of a 60-gene
query. The second term handles a state whose own marker set is smaller than the query. A
state holding 20 markers cannot share 10 genes with a 60-gene query without agreeing on half of
everything it has. The second term lowers the requirement to `ceil(0.17 * 20) = 4` for that state, so
small states remain reachable. 0.17 equals 10 divided by 60, so both terms state the same
requirement at the normal query size.

**Condition 2, significance.** BH-adjusted FDR ≤ `--nomination-fdr`, default 0.05.

**Condition 3, specificity, disabled by default.**

```
margin = (k_best - k_second) / k_best  >=  --nomination-specificity     default 0.0
```

Section 14.2 gives the measurement that set this default to 0.0.

The module excludes a cluster meeting the conditions in either database, then writes the matched
state, the overlap and the FDR to `nomination_decisions.tsv`. A cluster meeting them in neither
database proceeds to Gate 2.

### 9.5 Measured separation

Across the 151 clusters tested from the three non-seed datasets:

| Group | n | Overlap with best state, r > 0.25 | | |
|---|---|---|---|---|
| | | median | IQR | range |
| Excluded by Gate 1 | 110 | 24 | 15 to 36 | 10 to 59 |
| Admitted as new states | 41 | 6 | 5 to 7 | 1 to 10 |

The two groups separate almost completely at the threshold of 10. The largest FDR among excluded
clusters is 1.40 × 10⁻⁵, well inside 0.05, so the overlap requirement and not the FDR is the binding
condition. Chance expectation for a 60-gene query against a 39-gene database over 6,715 genes is
0.35 shared genes.

---

## 10. Step 7 — Gate 2, survival under competitive marker assignment

Gate 1 asks whether the reference already describes the candidate. Gate 2 asks a different question:
can the candidate hold marker genes of its own when it competes against the entire reference at once?

### 10.1 Why a second gate is necessary

MarkerFinder assigns each gene to exactly one cluster. Adding a new cluster to a reference therefore
takes genes away from existing clusters. A candidate that survives Gate 1 might still be a slight
variant of an existing state, distinguishable by no gene at all once the two compete directly. Gate
1 cannot detect that case, because Gate 1 compares fixed marker lists rather than running the
assignment.

### 10.2 The reference panel

The module subsamples each reference state to `--survival-ref-cells` (default 60) cells, drawn with
`numpy.random.default_rng(0)` so the panel is reproducible. The module runs
`marker_finder_wrapper` once on the panel alone and records each state's baseline marker count.

**Justification for subsampling.** Only marker structure matters here, not statistical power on cell
counts. A panel of 60 cells per state keeps each MarkerFinder call to a few thousand rows. The
module runs one call per candidate, so panel size multiplies directly into runtime.

**Measured panel sizes.** 63 states and 3,756 cells for Basil-2022; 85 states and 5,076 cells for
Natri-2024; 99 states and 5,875 cells for Adams. All reference states held at least one marker at
baseline in all three steps.

### 10.3 The one-candidate-at-a-time test

The module adds **one** candidate to the panel, runs `marker_finder_wrapper` again, and applies two
rules:

**Rule 1, the candidate must earn markers.** The candidate must own at least
`--survival-min-markers` (default 1) genes at Pearson r ≥ `--survival-rho` (default 0.3).

**Rule 2, the candidate must not erase a reference state.** No reference state that met Rule 1's
requirement at baseline may fall below `--damage-floor` (default 1) markers when the candidate joins.

The module excludes a candidate failing either rule, then writes the failing rule and the numbers
into `withdrawn_reason`.

**Why one at a time, and not all candidates together.** Testing all candidates together cannot
attribute damage. Twenty new candidates spread a reference state's markers across many of them, and
a rule that then withdraws whichever took the most withdraws the wrong cluster. Measured example:
the Migratory dendritic cell clusters of Adams and Natri-2024, holding 19 and 26 unique markers of
their own, were withdrawn under the all-together test for erasing state S30. S30 is a fibroblast
state marked by CTHRC1, FAP, LUM and MXRA5, and it shares exactly one gene, CCL19, with those
clusters. Adding one candidate at a time isolates the effect of that candidate.

**Measured on the four lung runs.** Gate 2 excluded 0 of 22 Basil-2022 candidates, 0 of 14
Natri-2024 candidates and 0 of 5 Adams candidates at the default thresholds. Gate 1 does nearly all
of the exclusion at these settings, and Gate 2 acts as a check that no admitted cluster is markerless.

---

## 11. Step 8 — Order the states with HOPACH

The module runs HOPACH on the state centroids, states by 6,715 genes, with the `cosangle` metric,
`kmax=9`, `kmin=2`, `mincluster=2` and `random_state=0`.

HOPACH builds a hierarchical tree and then orders the leaves so that neighbouring leaves are
similar. AltAnalyze and ICGS heatmaps have always used that ordering. Ordering by merge order
instead would place related states arbitrarily far apart in the final figure.

**Measured.** 70 top-level groups over 104 states.

The module runs HOPACH twice, once before the survival report and once after, so the recorded order
always matches the states that survive.

Output: `hopach_state_order.tsv`, and the `hopach_position` and `hopach_cluster` columns of
`harmonized_states.tsv`.

---

## 12. Step 9 — Global survival report and final marker analysis

### 12.1 The global survival report

The module runs MarkerFinder once over every representative cell of every surviving state and counts
how many states own at least `--survival-min-markers` genes at r ≥ `--survival-rho`. Under
`--survival-mode report`, the default, the module records the states that fall short and keeps them.
Under `--survival-mode enforce`, the module merges each failing state into its nearest surviving
state by centroid cosine, subject to the one-cluster-per-dataset constraint, and repeats up to 10
rounds.

**Justification for the `report` default.** Winner-take-all assignment across 104 states starves
states whose genes a near neighbour claimed. Removing them cascades: each removal frees genes,
changes the assignment, and starves a different state in the next round. Measured result under
`enforce`: the reference collapsed below the state count of a single input dataset.

**Measured.** 102 of 104 states own at least one unique marker at r ≥ 0.3 in the full-pool report.

### 12.2 The final marker analysis

The module ranks the cells of each surviving state by `ICGS3_SVM_score`, keeps the top
`--final-cells-per-state` (default 100), and runs MarkerFinder on that reduced matrix.

**Measured.** 10,106 cells across 104 states, giving 2,728 markers across 103 states. Markers per
state: median 22, minimum 2, maximum 60. Marker Pearson r: minimum 0.300, median 0.420. One state,
S98, holds no marker at r ≥ 0.3 and stays in the table under `report` mode.

Cells by source dataset in the final analysis: PedDev 6,087, Basil-2022 2,200, Natri-2024 1,319,
Adams 500.

---

## 12A. Step 10 — Annotate every state

The module names every state by default. Two sources run, and the module adopts the better result.
`--no-annotation` turns the step off. `--species` sets the reference species and defaults to `Hs`.

### 12A.1 Source 1, reference cell type enrichment

Each source run's `icgs3_cell_barcode_clusters.tsv` carries a reference annotation per cell. The
module reads every column that is not one of the eight technical columns ICGS3 writes, so the
module finds a new reference column without a code change. `--annotation-columns` restricts the set.

For state S and reference type T the module counts

```
k = cells of S labelled T           n = cells of S carrying the column
K = all final cells labelled T      N = all final cells carrying the column
p = hypergeom.sf(k - 1, N, n, K)
```

and adjusts p across every state-by-type pair within one column, by Benjamini-Hochberg. The
universe is the final cells, so the denominator is the cell types this dataset actually contains
rather than the whole reference vocabulary.

`accuracy` is k / n, the fraction of the state's own annotated cells carrying T. The module scores precision, not
recall: a state is well named when its cells agree, and one cell type may legitimately split
across several states.

Two studies annotated with different vocabularies never enter the same test. `HLCA` covers Adams,
Basil-2022 and Natri-2024; `Hs-PedDev-Lung` covers PedDev. Pooling label strings across two
vocabularies would compare names that were never meant to match.

**Measured.** HLCA scored 41 states over 4,019 annotated cells and 50 cell types present.
`Hs-PedDev-Lung` scored 63 states over 6,087 cells and 49 cell types.

### 12A.2 Source 2, marker gene set enrichment

The module calls ICGS3's own `biomarker_enrichment` and `clean_biomarker_prediction_labels`, so
this module and the source runs call one implementation and one reference file. The default
reference is `/Users/saljh8/Documents/GitHub/altanalyze/AltDatabase/EnsMart72/goelite/Hs/gene-mapp/Ensembl-BioMarkers.txt`;
`--biomarker-file` overrides it.

**The module restricts the gene denominator.** `biomarker_enrichment` intersects every BioMarkers term
with the background before testing. The module passes the 6,715 union genes, so the universe is
the genes this integration measured, not the genome. A term keeps only its genes that this dataset
could have detected, and the test drops any term retaining fewer than 3 such genes.

`accuracy` is overlap / query size, the fraction of the state's markers inside the term.

**Name cleaning.** `clean_biomarker_prediction_labels` applies the legacy RNASeq.py rules: drop the
parenthetical citation, strip `Adult`, `Fetal`, `Embryonic`, `Embryo` and `Term`, then infer the
dominant tissue across clusters and remove it from any term that is not entirely that tissue.

`clean_biomarker_prediction_labels` records the tissue in a dictionary keyed by cluster, so the
last qualifying row per cluster casts that cluster's vote. ICGS3 passes one row per cluster, its top term. The integration
needs a cleaned name for every enriched term, because the accuracy rule can adopt a term that is
not the top hit, so the module sorts each cluster's top term last before calling the function. The
vote is then identical to ICGS3's.

**Check.** Both paths infer tissue `Lung`, and all 103 top-term names match ICGS3's own output
exactly. Cleaning turns `Adult Lung Serous (PMID...)` into `Serous` and `Adult Lung Dendritic
Cells (PMID...)` into `Dendritic Cells`.

### 12A.3 Adopting one annotation

A candidate must carry at least `--annotation-min-evidence` (default 4) matching cells or genes and
reach BH FDR at or below `--annotation-fdr` (default 0.05). Among the survivors the module takes
the highest `accuracy`, breaking ties on the smaller FDR, then the larger evidence, then the name,
so the result never depends on row order.

`state_annotations.tsv` records the adopted call and, separately, the best call from each source,
so a disagreement stays visible. `state_annotation_candidates.tsv` records every candidate.

**Measured.** 104 of 104 states named. Source 1 won 87, source 2 won 17. Median accuracy 0.97 for
source 1 and 0.75 for source 2; median evidence 97 cells and 7 genes. 90 states carry a call from
both sources.

### 12A.4 A known weakness of the accuracy rule

Accuracy divides by the state's own evidence, so a state with few markers reaches a high value
cheaply. Two states hold exactly 4 markers, all 4 fall in one term, and the resulting accuracy of
1.000 beats every reference-cell call.

| Marker count of the state | Source-2 wins |
|---|---|
| fewer than 10 markers | 8 of 17 |
| median markers, states won by source 2 | 11 |
| median markers, all states | 22 |

The biology shows what that costs. 10 of the 17 source-2 wins name a tissue absent from a
lung dataset: `Placenta Endo`, `Heart Dendritic cell`, `Stomach Smooth muscle cell`,
`Omentum Endothelial cell`, `Duodenum Enterocyte_BEST4 high`, `Esophagus Goblet cell`,
`kidney Natural killer T NKT cell`, `Pleura Unknown`, `Trachea Endothelial cell_SELE high` and
`Airway Epithelium_Plasschaert Secretory cell`. In each case a reference-cell call with a plausible
lung identity lost on accuracy alone. The rule gives S86 the name `Duodenum Enterocyte_BEST4 high`
at accuracy 1.000 from 4 genes, over `Mesothelium` at 0.44 from 44 cells.

The hypergeometric FDR already accounts for both set sizes, so it does not reward a small query.
Ranking the same surviving candidates by FDR instead of accuracy changes 16 of the 17 source-2
wins, and every replacement is a lung cell type: S10 becomes `VEC`, S102 `Classical monocytes`,
S103 `CD4 T cells`, S11 `EC venous systemic`, S21 `ASMC`, S27 `Alveolar FB`, S36 `AT0`.

The module implements the accuracy rule as specified. Changing the tie-break to FDR, or requiring a
minimum query size before accuracy counts, would remove the artifact and needs one instruction.

### 12A.5 Duplicate names

The 104 states carry 68 distinct names. Several states share one name because they are substates of
one cell type: five states carry `AM` or `Alveolar macrophages`. The `state` column remains the
unique identifier. ICGS3's own protocol appends `_c<cluster>` to guarantee uniqueness, and the
module strips that suffix because the state identifier is already a column.

---

## 13. The MarkerFinder algorithm, specified

Every gate and every marker step calls `marker_finder_wrapper` from
`altanalyze3/components/udon/markerFinder.py`. The algorithm runs as follows.

1. **Build indicator vectors.** `pandas.get_dummies(groups)` turns the cluster labels into a
   cells-by-clusters 0/1 matrix. Column c holds 1 for a cell in cluster c and 0 otherwise.
2. **Correlate.** The module computes the Pearson correlation of every gene's expression vector
   against every indicator column, across all cells.
3. **Test.** The module converts r to t with `t = r * sqrt(df / (1 - r²))`, where `df = n_cells - 2`,
   and clips r strictly inside ±1 so that r = ±1 does not produce an infinite statistic.
4. **Assign uniquely.** Each gene goes to the single cluster giving its highest r, by
   `numpy.argmax` across columns. A gene therefore marks exactly one cluster.
5. **Filter by correlation.** The module drops assignments below `rho_threshold`.
6. **Filter by cluster.** The module drops clusters holding fewer than `min_markers_per_cluster`
   surviving markers.
7. **Cap.** The module keeps the top `top_n` markers per cluster, then applies `marker_finder_rho`
   a second time.

The module calls the wrapper with `rho_threshold` and `marker_finder_rho` both set to
`--survival-rho`, `top_n` set to `--marker-top-n`, and `min_markers_per_cluster` set to
`--survival-min-markers`.

**Consequence to note.** `--survival-rho` therefore also sets the threshold of the final marker
table at `ICGS_integrate.py:1065`. Raising it to 0.4 lowers the reported marker count from 2,728 to
1,643 without changing which states exist. A reader comparing marker counts across parameter
settings must account for that coupling.

---

## 14. Parameters, defaults and effect on state removal

| Option | Default | Step | Raising it |
|---|---|---|---|
| `--seed-dataset` | most clusters | 5 | see section 15.3; a coarser seed gives fewer states |
| `--dataset-order` | descending cluster count | 5 | overrides `--seed-dataset`; changes the result |
| `--cells-per-cluster` | 200 | 1 | sharpens profiles, mildly more permissive |
| `--marker-top-n` | 60 | all | more markers reported per cluster |
| `--nomination-query-top` | 60 | 6 | larger query, more chances to overlap, more exclusion |
| `--nomination-min-overlap` | 10 | 6 | **strongest control on removal**; more exclusion |
| `--nomination-overlap-fraction` | 0.17 | 6 | raises the bar only for small reference states |
| `--nomination-fdr` | 0.05 | 6 | more exclusion; not binding at current settings |
| `--nomination-specificity` | 0.0 | 6 | **keeps cross-study duplicates apart**; fewer exclusions, more states |
| `--survival-ref-cells` | 60 | 7 | slower, more stable marker structure |
| `--survival-rho` | 0.3 | 7, 9 | more exclusion, and a smaller final marker table |
| `--survival-min-markers` | 1 | 7 | more exclusion |
| `--damage-floor` | 1 | 7 | more exclusion by the erasure rule |
| `--survival-mode` | report | 9 | `enforce` removes low-marker states and cascades |
| `--final-cells-per-state` | 100 | 9 | larger final analysis |
| `--require-multi-dataset` | off | 9 | keeps only corroborated states |
| `--no-rebuild` | off | 9 | skips Steps 9 and 10 entirely |
| `--species` | Hs | 10 | selects the BioMarkers reference; Hs or Mm |
| `--biomarker-file` | ICGS3 default | 10 | overrides the BioMarkers reference |
| `--annotation-columns` | every non-technical column | 10 | restricts source 1 |
| `--annotation-min-evidence` | 4 | 10 | fewer states named |
| `--annotation-fdr` | 0.05 | 10 | fewer states named |
| `--no-annotation` | off | 10 | skips annotation entirely |

### 14.1 Measured threshold ablation

Four configurations ran on the same four datasets on 2026-08-21. Every other parameter stayed at
its default.

| Specificity | Survival | States | Ground truth |
|---|---|---|---|
| 0.30 | ≥1 marker at r 0.3 | 114 | 31 of 31 |
| **0.00** | **≥1 marker at r 0.3** | **104** | **31 of 31** |
| 0.00 | ≥3 markers at r 0.3 | 87 | 30 of 31 |
| 0.00 | ≥1 marker at r 0.4 | 92 | 31 of 31 |

### 14.2 Why `--nomination-specificity` defaults to 0.0

Setting the margin to 0.30 kept 10 extra clusters. The module recorded each one, and a centroid
correlation against the state that absorbed it, computed over 6,715 genes and placed against all
6,441 cross-state pairs of that reference (median r 0.320, 99th percentile 0.829), gave:

| Cluster | Label | Pair r | Percentile |
|---|---|---|---|
| Basil2022 C28 | Alveolar macrophages | 0.938 | 100.0 |
| Adams C35 | AT2 | 0.910 | 99.9 |
| Basil2022 C11 | EC arterial | 0.903 | 99.8 |
| Natri2024 C21 | Ionocyte | 0.901 | 99.8 |
| Basil2022 C21 | Alveolar fibroblasts | 0.846 | 99.2 |
| Natri2024 C52 | Alveolar fibroblasts | 0.829 | 99.0 |
| Adams C4 | DC1 | 0.819 | 98.8 |
| Natri2024 C4 | EC arterial | 0.805 | 98.3 |
| Adams C26 | Pericytes | 0.625 | 92.8 |
| Adams C1 | Alveolar macrophages | 0.524 | 87.3 |

Eight of ten sit above the 98th percentile. Basil2022 C28 gave the most similar state pair in the
entire 114-state reference.

A direct marker check confirms the reading. Two clusters carry an ionocyte label: Basil-2022 C49
(119 cells, purity 0.99) and Natri-2024 C21 (128 cells, purity 0.73). At specificity 0.30 the module
kept both, and the canonical ionocyte panel split across the two states: ASCL3, ATP6V1G3, BSND,
CFTR, DMRT2, STAP1 and TMPRSS11E marked one state, while ATP6V1B1, CLCNKB, FOXI1 and HEPACAM2 marked
the other. FOXI1 is the defining ionocyte transcription factor, so neither state carried a complete
signature. At specificity 0.00 all 11 panel genes mark one state.

### 14.3 Why `--survival-min-markers` stays at 1 and `--survival-rho` at 0.3

Requiring 3 markers lost EC venous systemic, a cell type Basil-2022 C6 and Natri-2024 C6 both
resolve. Losing a population two studies independently resolve fails the ground-truth standard in
section 15.1.

Raising the correlation to 0.4 kept ground truth intact but removed 12 clusters, 11 of them under
Rule 2 rather than Rule 1. The median removed candidate held 60 unique markers of its own. In 4 of the 11 cases, a reference state holding one or two markers at baseline cast the veto. One example:
Basil-2022 C17, annotated Myofibroblasts at purity 0.99 with 283 cells and 37 unique markers, was
removed because admitting it cost state S26 its single marker. Raising `--survival-rho` thins every
reference state's marker list until single-marker states can veto strong candidates.

---

## 15. Validation

### 15.1 The external ground-truth standard

Module: `ICGS_integrate_groundtruth.py`.

The standard never reads the integration's own decisions. A cell type counts as reproducible when
at least `--min-studies` (default 2) HLCA-annotated studies each resolve it as a cluster whose
dominant label reaches `--min-purity` (default 0.5). The standard drops PedDev, because
PedDev annotates with the `Hs-PedDev-Lung` vocabulary and shares only 9 label strings with HLCA.
Comparing label strings across two vocabularies would be invalid.

For each reproducible cell type the standard pools the markers of the exemplar clusters, tests that
pooled set against every final state by hypergeometric enrichment, and requires the best state to
reach BH FDR ≤ 0.05 and an overlap of
`max(--min-overlap-floor, ceil(--overlap-fraction * min(|query|, |state markers|)))`. The adaptive
form is necessary: the median state of this reference holds 22 markers, so a fixed requirement of
10 would demand 45% of everything that state has.

**Result: 31 of 31 reproducible cell types represented. Zero false negatives.**

### 15.2 Additional benchmarks

Module: `ICGS_integrate_benchmark.py`. Results on the 104-state reference:

| Benchmark | Result |
|---|---|
| Label recovery | 74 of 74 input annotation labels represented |
| Known biology, aberrant basaloid | 6 of 7 panel genes mark S94 |
| Known biology, neuroendocrine | 4 of 4 panel genes mark S34 |
| Known biology, mesothelium | 4 of 4 panel genes mark S80 |
| Known biology, proliferating | 5 of 5 panel genes mark S74 |
| Known biology, CTHRC1 fibroblast | 2 of 6 panel genes mark S88 |
| Marker self-consistency | 103 of 104 states hold a marker |
| One cluster per dataset per state | PASS |
| Unique state identifiers | 104 of 104 |

### 15.3 Order permutation

Forcing the smallest dataset to seed the reference, with `--seed-dataset Adams`, gives:

| Entry order | Seed clusters | Final states | Ground truth |
|---|---|---|---|
| PedDev, Basil-2022, Natri-2024, Adams (default) | 63 | 104 | 31 of 31 |
| Adams, PedDev, Basil-2022, Natri-2024 | 43 | 88 | 31 of 31 |

Seeding from Adams costs 16 states, a drop of 15%. Every reproducible cell type survives both
orders, so the loss falls on states one dataset alone resolves, not on shared biology. Under the
Adams seed, Gate 1 excluded 38 of 63 PedDev clusters, 41 of 55 Basil-2022 clusters and 44 of 53
Natri-2024 clusters, against 33, 39 and 38 under the default.

The two runs together support the default rule. A seed holding 43 clusters cannot contain the
distinctions a 63-cluster partition draws, so later clusters that the finer seed would have
absorbed enter as states, while other genuine distinctions collapse into the coarser seed states.
Outputs: `/Users/saljh8/Dropbox/Transfer/ICGS3_order_seedAdams/`.

The result also sets the honest error bar on the headline number. The state count carries a
roughly 15% dependence on entry order. The set of reproducible cell types carries none.

### 15.4 Benchmarks not yet run

Two checks remain outstanding, and no claim in this document rests on them.

1. **Negative control.** Split one dataset in half, treat the halves as two datasets, and integrate.
   A correct method should recover close to the original cluster count, not double it.
2. **Leave one dataset out.** Build the reference without one dataset, then align that dataset's
   clusters to the result.

### 15.5 Residual redundancy

Nine of 5,356 state pairs in the 104-state reference exceed centroid r 0.90, with a maximum of
0.938. Some cross-study duplicates therefore survive. Two changes would address them, and neither
is implemented:

1. Require a reference state to hold at least 3 markers at baseline before Rule 2 lets it veto a
   candidate. Single-marker states currently block candidates holding 37 to 60 markers.
2. Merge state pairs above a stated centroid-correlation threshold after Step 8.

---

## 16. Reproduction

```bash
cd /Users/saljh8/Documents/GitHub/altanalyze3
PYTHONPATH=. /opt/homebrew/opt/python@3.11/bin/python3.11 -u \
  -m altanalyze3.components.clustering.ICGS_integrate \
  --run Adams=/Users/saljh8/Dropbox/Transfer/ICGS3_Adams_corr03_ne3 \
  --run Basil2022=/Users/saljh8/Dropbox/Transfer/ICGS3_Basil2022 \
  --run PedDev=/Users/saljh8/Dropbox/Transfer/ICGS3_PedDev \
  --run Natri2024=/Users/saljh8/Dropbox/Transfer/ICGS3_Natri2024 \
  --output-dir /Users/saljh8/Dropbox/Transfer/ICGS3_integration_final \
  > /Users/saljh8/Dropbox/Transfer/ICGS3_integration_final/integration.log 2>&1
```

Every threshold in section 14 is at its default, so the command needs no threshold argument.
Validation:

```bash
PYTHONPATH=. /opt/homebrew/opt/python@3.11/bin/python3.11 \
  -m altanalyze3.components.clustering.ICGS_integrate_groundtruth \
  --integration /Users/saljh8/Dropbox/Transfer/ICGS3_integration_final \
  --run Adams=/Users/saljh8/Dropbox/Transfer/ICGS3_Adams_corr03_ne3 \
  --run Basil2022=/Users/saljh8/Dropbox/Transfer/ICGS3_Basil2022 \
  --run PedDev=/Users/saljh8/Dropbox/Transfer/ICGS3_PedDev \
  --run Natri2024=/Users/saljh8/Dropbox/Transfer/ICGS3_Natri2024
```

---

## 17. Outputs

All files below sit in `/Users/saljh8/Dropbox/Transfer/ICGS3_integration_final/`.

| File | Content |
|---|---|
| `harmonized_states.tsv` | 104 states: members, source dataset, supporting datasets, cell count, GO-Elite name, HOPACH position |
| `nomination_decisions.tsv` | one row per tested cluster: action, matched state, overlap and FDR in both databases, unique markers, exclusion reason |
| `hierarchical_audit.tsv` | the same, including the 63 seed states |
| `final_markers.tsv` | 2,728 markers across 103 states from the top-100-cell analysis |
| `harmonized_markers.tsv`, `harmonized_markers_all.tsv` | markers from the full-pool analysis |
| `harmonized_final_centroids.tsv` | 6,715 genes by 104 states |
| `harmonized_centroids.tsv` | centroids in HOPACH order |
| `harmonized_cell_annotations.tsv` | 10,106 cells: state, state name, source dataset, source cluster |
| `harmonized_integrated.h5ad` | the final matrix with `layers['counts']` |
| `hopach_state_order.tsv` | state order and HOPACH group |
| `states_removed_by_survival.tsv` | states below the marker requirement, kept under `report` mode |
| `integration_summary.json` | every parameter that decided the run, plus the state counts |
| `state_annotations.tsv` | adopted annotation per state, plus the best call from each source |
| `state_annotation_candidates.tsv` | every annotation candidate with evidence, accuracy and FDR |
| `GO-Elite/icgs3_biomarker_enrichment.tsv` | source-2 enrichment, gene denominator restricted to the union genes |
| `GROUND_TRUTH.tsv` | per-cell-type representation result |
| `BENCHMARKS.md` | the benchmark report of section 15.2 |
| `MarkerFinder/harmonized_marker_heatmap.pdf` | heatmap of 10,106 cells, with a `source_dataset` covariate bar |
| `integration.log` | the full run log |

Ablation outputs sit in `/Users/saljh8/Dropbox/Transfer/ICGS3_ablation_comparison/` and in the three
`/Users/saljh8/Dropbox/Transfer/ICGS3_ablation_spec0*/` directories.

---

## 18. Result summary

The module produced 104 states from 214 input clusters. PedDev contributed 63, Basil-2022 22,
Natri-2024 14 and Adams 5. Cross-study support: 17 states corroborated by all four datasets, 11 by
three, 14 by two, and 62 seen in one dataset only. Novel additions fell from 22 to 14 to 5 across
the three non-seed datasets, which is the saturation a converging reference should show.

---

## 19. Changes from the earlier draft description

An earlier description of this work described a design the module no longer runs. The differences
matter to anyone comparing the two.

| Earlier description | Current module |
|---|---|
| Clusters assigned to states by maximum-weight bipartite matching on rank-weighted Jaccard | No matching step. Gate 1 tests each cluster independently by enrichment |
| Redundant clusters merged into the matched state | Nothing merges. The redundant cluster is excluded and the state keeps its own cells |
| Marker databases built by correlating genes against state indicators | Databases built from ICGS3's unique-marker assignment. Section 9.2 gives the measurement that forced the change |
| Query of top 50 markers, 10 shared genes required | Query of top 60, requirement `max(10, ceil(0.17 * min(|Q|, |D|)))` |
| Merging blocked when the state already held a cluster from that dataset, and the candidate retained | No merging exists, so the constraint is structural. The module asserts it and logs PASS |
| A pairwise consolidation pass merged states at Jaccard 0.30 and cosine 0.80 | Removed. `--redundancy-score` and `--redundancy-cosine` were deleted on 2026-08-21 |
| 107 states, then 106 after consolidation, 1,295 markers across 82 states at r ≥ 0.4 | 104 states, 2,728 markers across 103 states at r ≥ 0.3 |

The functions implementing the superseded design carried no caller. The module lost
`_rank_weighted_jaccard`, `state_correlations` and `build_marker_database`, and six command-line
options, on 2026-08-21. Removing them changed no output: all six result files matched byte for byte
before and after.
