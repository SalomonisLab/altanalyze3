# CellRef2 reference construction: ICGS-integrate, CellCards annotation, and UPenn validation

## Scope and provenance

This document expands the reference-construction methods summarized on the LungMAP
[Vocabulary page](https://www.lungmap.net/research/vocabulary/) and in the downloadable
[AltAnalyze3 command notebook](https://www.lungmap.net/research/vocabulary/downloads/notebook).
It focuses on the three steps that determine which transcriptional states enter CellRef2 and
which biological names they receive:

1. cross-study state consolidation with ICGS-integrate;
2. independent marker-panel annotation and rescue with LungMAP Cell Nomenclature and CellCards;
3. projection into the UPenn source collection, followed by validation and filtering of the
   final terms.

The description separates the general behavior of the AltAnalyze3 modules from the parameters
used for the CellRef2 build. This distinction matters because several optional redundancy gates
implemented in `ICGS_integrate.py` were deliberately disabled in the CellRef2 command.

## Methodological overview

The analysis treated a *cell state* as a transcriptionally coherent cluster supported by unique
marker genes. It did not concatenate all cells and batch-correct them into a common embedding.
Instead, each study was clustered or read using its supplied annotations, MarkerFinder derived a
marker signature for every input population, and ICGS-integrate admitted datasets sequentially.
The first dataset established the initial reference; every later cluster was tested against all
populations evaluated previously. Consequently, the declared dataset order was part of the
statistical design, not merely an implementation detail.

The complete state-selection sequence was:

```text
717 input clusters
  -> MarkerFinder qualification and staged ICGS-integrate
  -> 142 nonredundant candidates
  -> marker, GO-Elite, centroid, and LungMAP nomenclature review
  -> 86 curated candidates
  -> CellCards/LungMAP marker-panel rescue of 5 candidates
  -> separate addition of a Langerhans candidate
  -> 92 states projected into the UPenn collection
  -> removal of 5 states after the first projection review
  -> removal of Preterminal bronchiolar secretory after cross-dataset review
  -> 86 final CellRef2 v7 states
```

The resulting v7 reference contains 40,383 UPenn cells, 35,049 genes, and 86 cell states. The
larger CellRef2 compendium was then classified against this reference; it is not the source of the
v7 reference centroids.

## 1. Input population definition and MarkerFinder qualification

Seventeen human lung or airway scRNA-seq and snRNA-seq datasets contributed candidate states.
Author annotations were used when suitable annotations were available. Ten datasets without the
required labels were analyzed independently with ICGS3, using non-negative matrix factorization,
graph clustering, MarkerFinder, and linear-SVM reclassification. The representative command in
the command notebook used 3,000 variable genes, 30 neighbors, Leiden resolution 0.8, a fixed
random seed of 0, and up to 60 markers per cluster.

### 1.1 Definition of a unique marker

MarkerFinder correlates each gene's expression across cells with an idealized binary membership
vector for each cluster. For gene \(g\) and cluster \(c\), the score is

\[
r_{gc}=\operatorname{cor}\left(x_g, I_c\right),
\]

where \(x_g\) is the expression vector and \(I_c\) is 1 for cells in cluster \(c\) and 0 for all
other cells. Each gene is assigned to the single cluster having the greatest positive Pearson
correlation. This winner-take-all assignment prevents a lineage-wide gene from being counted as a
unique marker for every related population in the same dataset.

For CellRef2 integration, an input cluster had to possess at least three uniquely assigned markers
with Pearson \(r>0.3\). Forty-seven of the 717 input clusters did not satisfy that prerequisite.
They were passed to ICGS-integrate through `--suppress-cluster`, so they could not become reference
states. MarkerFinder was explicitly selected with `--marker-method markerfinder`; the unrelated
Scanpy marker test exposed by the heatmap command was not used. Integer counts were supplied with
`--layer counts --scale-data`, allowing MarkerFinder to perform its own depth normalization.

### 1.2 Files consumed by ICGS-integrate

Each dataset was represented as a completed ICGS3-style run. ICGS-integrate reads:

- the cell-to-cluster table and ICGS3 SVM scores;
- the all-gene MarkerFinder correlations, including each gene's winning cluster and Pearson
  correlation;
- the MarkerFinder centroid/heatmap gene list;
- the corresponding AnnData expression object and, when present, its integer-count layer; and
- optional GO-Elite cell-type predictions and source annotation columns.

These inputs preserve both evidence types needed downstream: discrete marker membership for
redundancy testing and expression profiles for competitive MarkerFinder checks and centroid
summaries.

## 2. Staged integration with ICGS-integrate

### 2.1 Representative cells and common feature universe

Within each source cluster, cells were ordered by their ICGS3 linear-SVM score and at most 200
cells were retained. High-scoring cells lie farther inside their assigned decision region and
therefore provide less ambiguous estimates of the cluster expression profile. Clusters with fewer
than 200 cells contributed all available cells.

The module formed the union of genes present in the contributing MarkerFinder centroid files. It
then re-read expression for the representative cells and rebuilt every cluster centroid in this
single union feature space. This avoids comparing centroids defined over different gene sets. It
also defines one gene universe, \(N\), for all hypergeometric overlap tests in a run.

### 2.2 Fixed reference priority

The seed and entry order were fixed explicitly. Annotated references entered before unsupervised
candidates:

1. Human Lung Cell Atlas (HLCA), the seed;
2. LungMAP Human Lung CellRef v1.1;
3. BPD/Sun;
4. TGEN interstitial lung disease;
5. COPD/Zhang;
6. UPenn pediatric/postnatal collections;
7. COVID/Deutsch;
8. unsupervised UPenn pediatric development;
9. Basil 2022;
10. Natri 2024;
11. Adams 2020;
12. combined BOS, PAM, and submucosal-gland collection;
13. ACDMPV/Guo;
14. ILD/Jaiswal;
15. BPD/Sun;
16. LAM/Olatoke; and
17. lifespan/Wang.

The seed clusters are admitted without redundancy testing. Every cluster in a later dataset is
compared with the growing reference, so a population represented in a higher-priority dataset is
preferentially retained over its lower-priority counterpart. ICGS-integrate records support from a
redundant later cluster but does not merge that cluster's cells into the earlier state. Thus each
retained state remains represented by a single source cluster during this integration stage.

### 2.3 Marker-set redundancy test

For every retained state, the implementation constructs two marker databases from its source
cluster's unique-marker table:

- database A: markers with \(r>0.25\), ranked by \(r\), capped at 200 genes;
- database B: markers with \(r>0.30\), ranked by \(r\), capped at 100 genes.

For an incoming cluster, the query is its 60 highest-correlated unique markers in the shared gene
universe. For query set \(Q\), target marker set \(T\), overlap \(k=|Q\cap T|\), and universe size
\(N\), the enrichment probability is

\[
p=P(X\geq k)=\operatorname{hypergeom.sf}(k-1;N,|T|,|Q|).
\]

Within each database, p values are Benjamini-Hochberg adjusted across targets having at least one
overlapping gene. Targets are ranked first by overlap size and then by nominal p value. In the
CellRef2 build, an incoming cluster was called redundant when either database produced a best
target with at least 12 shared genes and FDR \(\leq0.05\).

The exact build settings are important:

| Parameter | CellRef2 value | Consequence |
|---|---:|---|
| `--nomination-query-top` | 60 (module default) | tests the 60 strongest incoming markers |
| `--nomination-min-overlap` | 12 | requires at least 12 shared genes |
| `--nomination-overlap-fraction` | 0 | disables the additional relative-overlap rule |
| `--nomination-fdr` | 0.05 | controls target enrichment within each marker database |
| `--nomination-identity-min` | 0 | disables the optional top-10 identity-marker gate |
| `--redundancy-top-n` | 0 (default) | disables the optional head-to-head marker-list gate |
| `--redundancy-centroid-r` | 0 | disables centroid correlation as a redundancy condition |
| `--nomination-specificity` | 0 (default) | does not require the best hit to exceed the second hit by a margin |
| `--compare-to-excluded` | enabled | prevents a later dataset from reintroducing an already rejected population |

Centroid correlation was disabled because equivalent populations measured by single-nucleus or
10x 5-prime chemistry could have correlations of only 0.27-0.36. In this build the statistical
redundancy decision therefore rests on unique-marker overlap and its FDR, not on a latent-space or
centroid-distance cutoff.

When `--compare-to-excluded` is active, clusters excluded during an earlier integration step remain
as marker-set comparison targets. A later cluster matching one of them is also excluded. Clusters
explicitly removed by `--suppress-cluster` are qualification failures and are not used as such
targets.

### 2.4 Competitive survival test

Failure to match an earlier marker set does not by itself establish a new state. Every provisional
candidate is therefore tested separately against the current reference with MarkerFinder. Up to
60 cells per existing state are sampled with random seed 0; all available representative cells of
the candidate are added; and unique markers are reassigned competitively at \(r\geq0.3\).

A candidate must:

- have at least five cells available for the test;
- retain at least one unique marker against the current reference; and
- not erase the last qualifying marker of an existing state that has at least as much marker
  support as the candidate.

This one-candidate-at-a-time design attributes any loss of marker specificity to a particular
candidate. It also means the public summary “12 of 60 markers at FDR 0.05” describes the principal
redundancy gate, but not the entire admission procedure. The independent survival test remains
active under the published CellRef2 command.

### 2.5 Integration result and audit trail

The staged procedure reduced 717 input clusters to 142 nonredundant candidate states, a 5.0-fold
reduction. The implementation writes, among other artifacts:

- `representative_cells.tsv`, the selected cells and source clusters;
- `hierarchical_audit.tsv`, one integration action per source cluster;
- `state_retention_evidence.tsv`, the matched target, overlap, required overlap, FDR, optional-gate
  values, failed test, and decision reason;
- `harmonized_states.tsv`, retained states and cross-dataset support;
- harmonized centroids and MarkerFinder tables; and
- a harmonized AnnData object retaining source-dataset and source-cluster provenance.

These tables should be treated as the machine-readable record of automated retention. Later
biological curation is recorded separately, rather than being overwritten into the automated
decision table.

## 3. CellCards and LungMAP nomenclature annotation

### 3.1 Initial biological interpretation

The 142 automated candidates were reviewed using four complementary evidence streams:

1. the candidate's uniquely assigned marker genes;
2. GO-Elite BioMarkers enrichment over those markers;
3. its nearest retained state by centroid similarity and marker overlap; and
4. its position and expected term in the LungMAP Cell Nomenclature hierarchy.

GO-Elite used HGNC symbols as the query space and hypergeometric over-representation with ontology
graph prioritization. The notebook records `min-term-size=5`, `max-term-size=2000`, minimum z score
1.96, maximum FDR 0.1, minimum overlap 2, and delta-z 0.5. Author labels were evidence, not binding
names: marker and enrichment evidence could override a source annotation when they disagreed.

This review reduced the 142 integration candidates to 86 provisional named states. The analysis
recorded, for every decision, the closest state by centroid correlation, the closest state by
marker overlap and its hypergeometric FDR, the shared genes, the nomenclature terms on both sides,
and a written retention/removal reason.

### 3.2 Independent CellCards marker-panel test

CellCards validation was run over all 717 input clusters, not only over the retained candidates.
The validation collection comprised 331 published marker panels assembled from 307 named cell
populations across 13 LungMAP Cell Type Database datasets together with 45 CellCards cell types.
For each input cluster, the query comprised its 100 highest-correlated unique markers. Each panel
was tested for hypergeometric enrichment over the 13,601-gene universe represented by the panel
collection, with false-discovery correction across the tested terms.

Testing every input cluster served two purposes. First, it assigned external biological evidence
without reusing the integration decision. Second, it could detect a marker panel whose strongest
carrier had been removed earlier. Five removed candidates were restored on that basis, increasing
the provisional atlas from 86 to 91 states. A Langerhans candidate was then added separately,
yielding the 92-state reference used for UPenn projection.

This CellCards test is a project-level curation step named `validate_cellcards_annotation.py` in
the command notebook; it is not a public top-level AltAnalyze3 entry point in this repository.
The notebook likewise invokes the analysis-local scripts `annotate_retention_evidence.py`,
`duplicate_evidence.py`, `apply_final_removals.py`, `assign_cellcards_names.py`, and
`build_final_atlas.py`. The statistical primitives are documented, but exact reproduction of the
curation requires the panel file and decision tables archived with the CellRef2 build.

### 3.3 Controlled vocabulary assignment

After state selection, each state was mapped one-to-one to a CellCards term. A broader parent term
was used when CellCards lacked an exact state-level match. Each final record then received:

- a standardized CellRef2 long and short name;
- a broad cellular superclass and intermediate family;
- the selected CellCards term;
- a Cell Ontology label and identifier; and
- a LungMAP Cell Nomenclature identifier and aliases when available.

All 86 final states have a Cell Ontology term, and 74 have a LungMAP identifier. Using a broader
ontology term does not imply that the transcriptional states were merged: the CellRef2 short name
continues to preserve distinctions such as activation, anatomical compartment, or transitional
state that may be finer than the available ontology node.

## 4. UPenn projection, validation, and final filtering

### 4.1 Validation collection

The 92 candidate states were projected into the UPenn scRNA-seq source collection assembled by the
Morrisey laboratory. This object contains 978,892 cells from 120 libraries and 83 unique donors:

| Cohort represented in the object | Libraries | Donors within cohort | Cells |
|---|---:|---:|---:|
| Morrisey 2026 LungMAP remainder | 68 | 54 | 523,972 |
| Guo 2023 Nature Communications | 39 | 24 | 354,570 |
| Mellors 2025 JCI Insight | 13 | 6 | 100,350 |
| Combined | 120 | 83 unique | 978,892 |

One donor is represented in two cohort records, explaining why the cohort donor counts sum to 84.
The collection contains 628,702 cells from 80 10x 3-prime v3 libraries and 350,190 cells from 40
10x 3-prime v2 libraries. Fifty-nine donors are normal controls; 24 represent COPD, LAM,
CLAD-associated or cGVHD-associated bronchiolitis obliterans syndrome, IPF, PVOD, sJIA-associated
lung disease, Hermansky-Pudlak syndrome, bleomycin-associated ILD, or other ILD. Libraries
`EEM-scRNA-R26-1` and `EEM-scRNA-R26-2` were excluded from every reference build.

### 4.2 cellHarmony assignment

The UPenn cells were aligned with `run_cellharmony_lite` using only genes present in the supplied
reference (`reference_genes_only=True`), cosine alignment, no minimum alignment score, and no UMAP
generation. For a query cell vector \(q\) and reference-state vector \(s_j\), the current aligner
computes

\[
\operatorname{score}(q,s_j)=\frac{q\cdot s_j}{\lVert q\rVert_2\lVert s_j\rVert_2}
\]

over the genes shared by the query and reference, and assigns the state with the maximum score.
Because `min_alignment_score=None`, every cell surviving input QC receives its best-scoring label;
the score is retained for ranking and review rather than used as a rejection threshold.

The checked-in implementation reads the reference matrix as supplied and L2-normalizes query and
reference vectors for cosine similarity. It does **not** z-score genes inside the cosine-alignment
branch. This differs from the notebook sentence stating that cellHarmony z-scores each gene across
reference states. Reproduction should therefore archive the exact code revision and supplied
reference matrix used for the build; the executable implementation, rather than that sentence,
defines the numerical transform.

### 4.3 Two-round state filtering

Projection was evaluated in two rounds:

| Round | State count | Removed candidates | Decision basis |
|---|---:|---|---|
| v2 to v3 | 92 to 87 | Terminal and respiratory bronchiolar secretory; Goblet (Bronchioles); CCL3-positive AM; Intravascular macrophage; Myeloid dendritic CD1c+ | insufficient unique-marker support after projection into UPenn |
| v5 to v6 | 87 to 86 | Preterminal bronchiolar secretory | not consistently separable from AT2 across datasets |

These decisions tested whether a discovery-derived term remained transcriptionally distinguishable
after assigning a large, independent source collection. The first round removed states that did
not retain adequate unique-marker support in the projected UPenn cells. The second removed a
secretory candidate whose separation from AT2 was not reproducible across datasets. The notebook
does not give a single numerical marker cutoff for these six curation decisions; they should not be
described as automatic applications of the earlier three-marker or 12-gene thresholds.

Author annotations from the contributing datasets were compared with the projected final labels as
an additional validation layer. This comparison assessed concordance and recognizable biological
relationships; author labels did not automatically replace a CellRef2 assignment.

### 4.4 Construction of the final v7 reference

The final reference was rebuilt exclusively from UPenn cells. Within each of the 86 retained
states, cells were ranked by cellHarmony alignment score and up to 500 were selected. Seventy-five
states reached the 500-cell cap; the remaining 11 contributed all available cells, from 27 to 491
per state. Random seed 0 was recorded for the build. The output contains 40,383 cells by 35,049
genes.

MarkerFinder was rerun on these selected cells using integer counts, internal scaling, 60 markers
per state, and up to 100 cells per state for the heatmap. It produced 5,138 unique markers and
21,500 redundant marker entries (250 per state). The corresponding state-centroid matrix became
the deployed cellHarmony reference. Thus the final marker definitions and centroids were estimated
from UPenn cells carrying the validated labels, not copied unchanged from the heterogeneous
discovery datasets.

## 5. Reproducible commands

The essential ICGS-integrate call was:

```bash
python3.11 -m altanalyze3.components.clustering.ICGS_integrate \
  --run HLCA=/path/HLCA/ICGS3_run \
  --run CellRef=/path/CellRef/ICGS3_run \
  --run BPD=/path/BPD/ICGS3_run \
  --run TGEN=/path/TGEN/ICGS3_run \
  --run COPD=/path/COPD/ICGS3_run \
  --run PedDev=/path/PedDev/ICGS3_run \
  --run COVID=/path/COVID/ICGS3_run \
  --run L-PedDev=/path/L-PedDev/ICGS3_run \
  --run L-Basil2022=/path/L-Basil2022/ICGS3_run \
  --run L-Natri2024=/path/L-Natri2024/ICGS3_run \
  --run L-Adams=/path/L-Adams/ICGS3_run \
  --run L-BOS_PAM_SMG=/path/L-BOS_PAM_SMG/ICGS3_run \
  --run L-ACD=/path/L-ACD/ICGS3_run \
  --run L-ILD=/path/L-ILD/ICGS3_run \
  --run L-BPD=/path/L-BPD/ICGS3_run \
  --run L-LAM=/path/L-LAM/ICGS3_run \
  --run L-Lifespan=/path/L-Lifespan/ICGS3_run \
  --output-dir /path/ICGS3_atlas_17ds \
  --seed-dataset HLCA \
  --dataset-order HLCA,CellRef,BPD,TGEN,COPD,PedDev,COVID,L-PedDev,L-Basil2022,L-Natri2024,L-Adams,L-BOS_PAM_SMG,L-ACD,L-ILD,L-BPD,L-LAM,L-Lifespan \
  --cells-per-cluster 200 \
  --final-cells-per-state 100 \
  --cells-from result \
  --nomination-min-overlap 12 \
  --nomination-overlap-fraction 0 \
  --nomination-identity-min 0 \
  --nomination-fdr 0.05 \
  --redundancy-centroid-r 0 \
  --compare-to-excluded \
  --suppress-cluster "$(cat /path/disq_17ds_combined_noAT0.arg)" \
  --species Hs
```

The UPenn projection used the programmatic entry point so that no CLI default could substitute an
alignment threshold:

```python
from altanalyze3.components.cellHarmony.run_cellHarmony_lite import run_cellharmony_lite

run_cellharmony_lite(
    gene_h5ad_path="/path/UPenn_with_umap_and_markers.h5ad",
    output_dir="/path/cellharmony_UPenn",
    cellharmony_ref="/path/CellRef2_candidate_reference.txt",
    barcode_cluster_out="/path/cellharmony_UPenn/UPenn_barcode_clusters.txt",
    sample_name="UPenn",
    alignment_mode="cosine",
    min_alignment_score=None,
    reference_genes_only=True,
    generate_umap=False,
)
```

## 6. Interpretation and limitations

- ICGS-integrate selects among already resolved clusters; it cannot split a biologically mixed
  input cluster. Its resolution ceiling is therefore set by the source annotations and ICGS3
  runs.
- Reference priority affects which study supplies the retained representative of a shared state.
  The fixed order is part of the method and must be reported with the thresholds.
- Hypergeometric significance establishes non-random marker overlap, not biological identity on
  its own. The build therefore combined overlap statistics with competitive marker survival,
  CellCards panels, GO-Elite results, centroid evidence, nomenclature review, and UPenn projection.
- The CellCards rescue and final naming steps include expert curation. Their decision tables and
  exact panel resource are required to reproduce the 92-state intermediate exactly.
- UPenn is both the validation source and the exclusive source of cells in the final v7 centroid
  reference. Projection tests generalizability relative to the discovery candidates, but it is
  not a fully held-out validation of the final UPenn-derived centroids.
- Three final states originated from fewer than 40 cells in their discovery datasets and require
  additional experimental or cross-cohort validation.
- The published marker panels did not reliably separate terminal/respiratory-bronchiole secretory
  cells from the corresponding alveolar type 0 population in these data; those terms were not
  retained.

## 7. Primary implementation and records

- ICGS-integrate implementation: `components/clustering/ICGS_integrate.py`
- MarkerFinder correlation core: `components/udon/markerFinder.py`
- MarkerFinder/heatmap entry point: `components/visualization/marker_heatmap_h5ad.py`
- cellHarmony wrapper: `components/cellHarmony/run_cellHarmony_lite.py`
- cellHarmony alignment implementation: `components/cellHarmony/cellHarmony_lite.py`
- detailed general ICGS-integrate behavior and ablations:
  `components/clustering/ICGS_INTEGRATE_METHODS.md`
- public methodological summary and controlled vocabulary:
  <https://www.lungmap.net/research/vocabulary/>
- exact CellRef2 run-order command notebook:
  <https://www.lungmap.net/research/vocabulary/downloads/notebook>

