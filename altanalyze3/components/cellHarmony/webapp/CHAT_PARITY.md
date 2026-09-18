# Shared Chat protocols and query performance

Both uploaded scALABLE and the COPD viewer now execute the same implementations
in `cellHarmony/chat_protocols.py`. Implemented protocols are clinical severity
gradients, ordered-stage trends, seed-gene coexpression, cell-state composition,
annotation concordance, donor-signature heterogeneity, most-affected states and
existing GO enrichment. Markers, expression, differential, TF activity, regulatory
networks and lipid/metabolite pathways retain their shared views.

## Data and statistical scope

COPD sample analyses use the existing 4,781 real sample/state pseudobulks from
141 samples, with original cell counts and an exact study/sample clinical join.
`export_chat_pseudobulk.py` builds a memory-mapped analysis matrix and opaque
sample labels inside the viewer database. It does not rerun differentials or
modify the LungMAP database. Stage trends exclude unrecorded/not-applicable
levels and use clinical order; normal spirometry precedes GOLD I–IV.

Uploaded donor expression sums raw RNA counts by biological sample within cell
state, then normalizes to log2(1+CP10k). At least five cells per sample and eight
contributing samples are required for expression trend/correlation protocols.
Library size uses the entire transcriptome. Candidate scans cover a deterministic
spread of up to 3,000 genes, stated in the answer. These are exploratory
associations; reported correlation p-values are unadjusted. Coexpression and
stage plots show the same within-state samples used for the statistics.

Composition uses each sample's fraction of original cells, with sample-level
rank tests. Missing/conflicting categorical sample assignments are excluded.
Constant genes have undefined correlation rather than a fabricated association.
A negative-only differential signature is supported without generating NaNs.

A missing raw-count layer, clinical field, ordered category definition, second
annotation or matching GO result produces a specific missing-data response.
COPD has no second independent annotation or v8 GO enrichment in this release;
those requests report the missing inputs. This is a data limitation, not an
unimplemented Chat protocol. The v7 GO results are not substituted for v8 calls.

## Performance and routing

Recognized dataset-specific questions route locally, including common marker,
expression and differential questions. Covariate routing excludes bookkeeping
fields such as `__n_obs`; FEV1 resolves to its clinical measurement. Ambiguous or
unrecognized language still uses the configured assistant. A missing required
slot produces choices rather than an unrelated answer.

Spearman scans use vectorized average-tie ranks, avoiding a gene-by-gene
correlation matrix. Each dataset caches up to 12 sample aggregations and 64
protocol answers. Upload adapters are bounded to four cached generations and
expire on changed metadata, comparisons or artifact timestamps. Expression
caches also expire when the uploaded h5ad is replaced in place. Caches are scoped
to job and dataset; no results are shared across uploaded users' jobs.

## Reproducible validation

From the repository root:

```bash
python -m pytest tests/test_chat_parity.py tests/test_grn_web_release.py tests/test_discover_integration.py tests/test_flask_pipeline.py -q
python tests/benchmark_chat_queries.py --base http://127.0.0.1:8062 --job Hs-Lung-COPD-metacells --out /tmp/chat-benchmark.json
```

The benchmark sends nine read-only question types three times and records HTTP
latency, resolved protocol, data availability and cache hits. `--questions`
accepts a JSON list for a different dataset. Test and browser reports are saved
with the local COPD release validation artifacts. Local timings describe these
recognized simulated queries, not arbitrary assistant-model requests.

## Expression cross-modality correlations

`webapp/cross_modal.py` is shared by the uploaded analysis application and the
published viewer. Its example questions work before a differential is run:

- Which TFs have activity discordant with gene expression across cell states?
- Correlate TF activity with matching gene expression across cell states.
- Correlate ADT abundance with its cell-surface gene expression across cell states.

The default analysis correlates paired cell-state means with Spearman's rho.
Discordance means a negative rho, not a comparison of absolute values on different
measurement scales. A named state instead uses matched donor/sample means within
that state. At least three groups with five matched observations per group are
required. Constant pairs and ambiguous gene mappings are excluded and reported.
No differential statistics, p-values, or replicate claims are inferred from these
Expression plots. Imputed activity/abundance is derived from RNA and does not
provide independent experimental validation.

ADT partners use the existing curated human marrow/lung or mouse mappings, with
clone-qualified aliases and explicit multi-subunit partners. For example,
CD10 maps to MME, CD56 to NCAM1, and CD8 to CD8A and CD8B. Gene-level PTPRC values
cannot distinguish CD45 protein isoforms. Observations are joined by identifier
before averaging. Stable float64 accumulation followed by restoration of source
precision preserves ties in repeated float32 imputed values. Only matched feature
columns are read; eight bounded cached summaries make repeat questions inexpensive.

Answers open a scatter plot by default, with a pair selector and a table of
correlations. The shared interface also hides modality controls for cell-type
UMAP and MarkerNetwork, reads MarkerNetwork folds from both unique and redundant
marker tables, and scrolls wide DotPlots horizontally. Group controls retain all
cell states even when there are more than 60. Selecting a completed differential
modality restores the saved comparison with matching groups and settings.
Uploads put their job_id in the browser URL immediately; that URL reloads the
session while its files remain available.

Validation:

```bash
python -m pytest tests/test_cross_modal_expression.py tests/test_plot_semantics.py tests/test_chat_parity.py tests/test_viewer_expression_bundle.py tests/test_combplot_cells.py -q
python tests/browser/check_expression_parity.py http://127.0.0.1:8000 JOB_ID
```

### Complete, sortable correlation results

Correlation responses retain every valid matched pair. The shared Chat result
controls support sortable column headers, a sort-field menu (including absolute
rho), text search, numeric rho bounds, row counts up to 2,000 or All, and paging.
Correlation tables and scatter-plot selectors use the same filters; every table
row can open its paired scatter plot. The default 50-row page does not truncate
the underlying result set.

“Uncorrelated”, “weakly correlated”, and “least correlated” questions select a
Near zero filter and rank by smallest absolute rho. The explicit default cutoff
is |rho| <= 0.2 and can be changed. This describes weak rank association, not a
claim of statistical independence. Negative-correlation queries preselect the
Negative filter; switching to All reveals the other pairs without another query.
Constant or insufficient pairs remain undefined and are not classified as zero.

`tests/browser/check_chat_result_controls.py URL` validates paging, field sorting,
search, numeric ranges, empty matches, and plotting a pair beyond the former
50-result cutoff on either application. Validated against 217 uploaded TF pairs
and 284 published COPD TF pairs; 23 targeted Python regression tests passed.

### Best marker across modalities

“What is the best modality marker of HSC-1?” routes locally to the shared
`modality_markers.py` implementation. It reads saved unique and redundant marker
tables and ranks retained positive MarkerFinder Pearson correlations against
0/1 cell-state membership. It does not use differential comparisons or compare
fold changes between different measurement scales. Responses identify the
strongest retained marker overall and within each modality, show a bar chart,
and retain all available positive rows in the sortable/searchable table.

Unknown-labelled metabolites remain unidentified. “Best named modality marker”
and a table search for Named exclude them. Specific ADT, lipid, metabolite, or
TF-activity marker requests restrict the same lookup to that modality. Missing
marker scores are reported explicitly; neither fold values nor another
modality's markers substitute for absent Pearson scores. Published datasets use
only their retained marker statistics, so coverage can differ from uploaded jobs.

The restored marrow session has 682 retained positive HSC-1 markers across five
modalities. IFIT2/RNA ranks first (r=0.5082); the other modality leaders are
Unknown 0400/metabolite (0.2527), Cer 36:0;2O|Cer 18:0;2O/18:0/lipid (0.2448),
HLF/TF activity (0.1925), and CD73/ADT (0.1500). No new marker or differential
analysis was run. `tests/test_modality_marker_chat.py` and
`tests/browser/check_modality_marker_chat.py URL STATE` cover ranking, named
filters, scope, missing data, routing, plotting, and table access.
