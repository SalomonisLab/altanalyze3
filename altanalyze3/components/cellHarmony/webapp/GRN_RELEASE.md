# GRN release preparation

The release separates `grn` (TF-to-target edges) from `grn_tf` (summed predicted
TF activity) in the uploaded-data scALABLE app and the precomputed scALABLE viewer.
All implementation lives in altanalyze3. The COPD deployment's release catalog
remains an external dataset configuration; no COPD analyses are rerun.

## Implementation and compatibility

- `cellHarmony/modalities.py`: shared identities, aliases, capabilities, and
  read-only interpretation of legacy uploaded job artifacts.
- `cellHarmony/flask/pipeline.py`: chunked per-cell TF activity, separate
  sample x state TF and edge inputs, markers, and separate downloads.
- `cellHarmony/grn_analysis.py`: shared activity profiles, network selection,
  state-specific filtering, exact sibling-comparison matching, and local
  regulatory intent routing. RNA expression and activity folds are kept separate.
- `cellHarmony/webapp/grn_data.py`: adapter for uploaded job files and retained
  completed runs. The same groups, population field, and comparison type must
  match before sibling statistics are joined.
- `cellHarmony/flask/tasks.py`: retains completed differential runs. The UI
  selects these without recomputation. A new alignment resets old comparisons.
- `visualization/scalable_viewer/grn_network.py`: compatibility import so existing
  integrations keep their entry points.
- Shared Chat plots color numeric TF fold changes by their sign, leave absent
  changes gray, and label absent network statistics as unreported rather than zero.

`Rna2GrnBundle.tf_activity()` retains its historical mean default. New upload
analysis explicitly requests `aggregation="sum"`. Legacy standardized target-set
enrichment remains labelled as legacy; it is not silently treated as predicted
activity. Reprocess a legacy upload to create the new TF differential input.

The existing GRN differential defaults remain alpha 0.05, fold threshold 1.1,
minimum two replicates per arm, using raw p for selection. Both p and FDR remain
in the output. Integrated views default to the stored differential calls, retaining
the original analysis thresholds. Network users can explicitly add an FDR or raw-p
cutoff; the interface labels that choice. Pathways never apply a second significance
cutoff to the stored calls.

## Evaluated functionality

| Function | Uploaded scALABLE | Precomputed viewer |
| --- | --- | --- |
| TF activity and edge differential inputs | Separate sample x state matrices | Separate existing facets |
| Summary, volcano, fold heatmap, detail distribution | Shared renderers | Shared renderers |
| PDF and differential archives | Available | Existing bundle assets/renderers |
| TF Explore: UMAP, violin, dotplot, donor/group plot | Per-cell predicted activity | Existing per-cell/metacell store |
| TF marker heatmap | When markers pass selection | When supplied by bundle |
| Edge network | Sample x state edge scores | Existing edge scores |
| Regulatory Chat | Inline network with matching activity/RNA statistics | Same analysis functions |
| Activity Chat | Ranked bars and reported differential statistics | Same analysis functions |
| Differential Chat | Correct matching completed modality | Correct matching precomputed modality |
| Saved comparison selection | Retained completed runs | Precomputed comparisons |

The formerly missing donor-level and annotation Chat protocols are now shared
between uploads and the viewer; see [CHAT_PARITY.md](CHAT_PARITY.md). Inputs that
a dataset does not contain are reported explicitly. Unrecognized language still
uses the configured assistant; supported local query forms avoid that round trip.
The integrated Network uses the maximum of the two sample-arm means on
log2(1+CP10k), matching Discover, in Differentials. Explore's marker regulatory
network instead uses stored RNA marker folds, state-specific model edge scores,
and mean stored RNA expression in the selected state, without differential runs.
The model-only GRN edges view defaults to a TF and a single cell state and uses
neutral colours and TF/target shapes rather than fabricated fold changes.

## Validation

Run from the repository root:

```bash
PYTHONPATH=. python -m pytest tests/test_grn_web_release.py tests/test_flask_pipeline.py -q
```

The GRN tests upload synthetic h5ad files through the FastAPI endpoint, exercise
alignment with deterministic GRN predictions, verify independent per-cell and
sample-level inputs, execute real TF/edge differential tests, and inspect
replicate counts, fold tables, plots, PDFs, retained comparisons, assistant
responses, and legacy behavior. The existing RNA upload test is also retained.
The synthetic imputer makes the integration test independent of model downloads.

Additional local release checks use the installed lung_hybrid model on synthetic
RNA, compare chunked factor sums with direct predictions, and exercise both
deployments in Chrome. The COPD catalog is loaded in an isolated viewer process
with a temporary runtime directory for regression testing.

Validated on 2026-09-15:

- Nine GRN regression tests and the existing RNA upload/alignment test passed.
- Chrome rendered TF/edge differentials and all three regulatory Chat plot types
  for an uploaded fixture and the COPD bundle, with no JavaScript or API errors.
  Direction colors and absent statistics were also checked.
- COPD exposes 13 TF activity and 13 edge comparisons. These checks reuse the
  published tables; no COPD analyses were computed.
- A built wheel installed into an isolated target directory served its template,
  JavaScript and CSS and loaded the lung hybrid model with 63,647 edges.
  Wheel contents exclude runtime jobs and repository build/temp directories.

## Packaging

`pip install '.[web]'` supplies the web runtime extra. The wheel explicitly
includes templates, JavaScript, CSS, reference configuration, shared Python
modules, and GRN model bundles. Reference atlas files remain deployment inputs.
No release number or target branch was specified, so project version 0.1.3 is
unchanged. This prepares the working tree and validation artifacts; it does not
publish a package, tag a release, or change the deployment's pinned commit.


## Integrated Discover release, 2026-09-16

`discover_integration.py` ports LungMAP Discover's edge gates, TF expression
floor, expression/activity colors, lipid-class matching and balanced pathway
ranking. `webapp/static/integrated.js` ports the original Cytoscape layout and
WikiPathways SVG renderer, with Show, gene/edge fold, p/FDR, expression, pathway,
zoom, vector PDF and TSV controls. Explore marker queries, Differential and Chat use
these same endpoints. Missing paired analyses produce the requested notice.

Pairing requires identical comparison groups, population field and replicate
type. Missing folds remain missing in heatmaps. TF activity falls back only to
matching comparisons. Metabolites map to exact common names in the packaged
WikiPathways diagrams; unrecognized names remain unmeasured. Lipid classes use
the original longest-synonym match and mean of significant species.

The COPD deployment stages existing LungMAP v8 **pseudobulk** results with
`export_discover.py` and `stage_discover.py`. Source SQLite is opened read-only;
no COPD differential analysis runs. The group means use 4,781 true sample/state
pseudobulks, with sibling libraries collapsed before normalization. Each sample
has equal weight. Mean log2(1+CP10k) and mean CP10k are exported separately.
`integrated_pseudobulk/manifest.json` records provenance; `validation.json`
records source row-count and independent arm-mean checks. Do not substitute the
old imputed meta-sample pseudobulks or an averaged-metacell matrix.

RNA, TF, edge, lipid and ADT tables in the staged release share this evidence.
Old RNA GO/interaction assets are not attached to changed v8 calls: they were
computed from v7 and would describe different gene lists. The uploaded app
continues to generate GO/interaction results when requested. Cell-communication
results and expression stores are retained. Original release catalogs/bundles
remain available for rollback.

Additional checks:

```bash
PYTHONPATH=. python -m pytest tests/test_discover_integration.py tests/test_grn_web_release.py tests/test_flask_pipeline.py -q
```

Browser checks cover COPD and uploaded differential selection, TF/edge volcanoes,
network and pathway controls, Chat, missing-pair messages, and PDF/TSV export.
Local browser evidence is in `/tmp/scalable-integrated-browser/`. Release 0.1.3
includes shared views and pathway resources; package publication is separate.

## Integrated selection and shared platform parity

Network target ranking now checks retained differential edges and TF expression
before applying the Show limit, including marker queries. Every displayed node
is incident to a retained edge, and TFs that are also targets have a single node.
The default selection is Reported calls; explicit FDR filtering can correctly
produce an empty network when the original analysis used raw p-values.

Pathway menus include only diagrams with a regulated lipid class or metabolite.
Gene-only matches do not qualify. Both the menu and diagram use the original
differential calls without additional significance or fold filtering. Legacy
pathway API significance parameters are accepted for compatibility but ignored.

Network and Pathway views have no Interpret button or interpretation panel,
including their Explore and Chat renderings. Structured result summaries remain
available to API consumers. Both uploaded datasets and the viewer share the same
renderers, selection rules, missing-pair notices, cell-state selection and exports.

The selected outer modality owns the view: Pathway is offered for lipids and
metabolites, and Network for RNA and GRN. Embedded Measurement and Features
selectors are removed. Differential uses differential features, Explore uses
markers, and Chat uses the requested feature source. Download PDF uses the existing
panel button in Differential/Explore and the same `ghost-btn` styling in Chat.
The shared jsPDF/svg2pdf pipeline exports vector paths and editable text, including
the full pathway independent of viewport zoom, network directions, color ranges,
and comparison context. It does not rasterize the figure or offer a second SVG
download control.

The Network Show selector offers 25, 50, 100, 200, 500, 1,000 and 2,000 targets
in Differential, Explore and Chat. Uploaded saved-comparison selection applies
the returned comparison immediately. Request guards prevent late plot or detail
responses from restoring an earlier comparison, modality or cell state.

Differential heatmaps preserve full fold matrices when available and leave absent
results null. A separate black heatmap layer displays nulls with a missing-result
hover label; it does not replace them with zero. PDF exports mask missing entries
as black. The shared uploaded-data and viewer implementation covers every modality.

## Imputed marker inputs

Lipid, ADT and TF activity marker analyses use model predictions directly for
MarkerFinder correlations. RNA sequencing-depth validation remains enabled for
RNA; it is not applied to continuous predicted modalities. Negative lipid
predictions are floored at zero before analysis and export, with the count of
clipped values recorded in the prediction summary. Centroid TSVs contain
per-state arithmetic means on the
prediction scale, preserving feature names such as `PC(16:0/18:1)`. RNA centroid
normalization and its negative-input guard are unchanged.

MarkerFinder uses float64 sparse variance to avoid cancellation for nearly
constant imputed features, bounds Pearson correlations to [-1, 1], and excludes
undefined tests from FDR adjustment. An empty marker selection is nonfatal for
RNA and all predicted modalities; expression and differential analysis remain
available. Other input errors still fail explicitly. Lipid prediction preserves
model feature names, eliminating the repeated sklearn warnings without disabling
model validation. Completed uploads can be reopened with `?job_id=<job_id>`.

Marker regulatory networks select only positive RNA marker log2 folds (including when the fold threshold is 1). Marker nodes use a red-only legend; expressed TFs without a qualifying positive marker fold stay gray. Differential networks retain their signed up/down selection and color scale. Negative values are excluded, never converted to positive folds.

## Interrupted uploaded comparisons

Saved queued/processing status is reconciled with the owning worker before status
is displayed or a comparison is submitted. An orphaned run becomes retryable;
its configuration, alignment outputs and completed comparison history remain.
Active local futures and live worker processes continue to block duplicate runs.
Reloading a completed alignment resumes polling when its differential is active.
Differential working data omits unused correction layers and other imputed
modalities, and avoids copying the full matrix when all samples are selected.
