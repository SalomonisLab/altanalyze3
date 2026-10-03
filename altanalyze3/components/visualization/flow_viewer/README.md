The viewer selects populations in measured CITE-seq data and evaluates their CITE-seq → flow assignments, FlowSOM overlap, marker expression and sorting gates. It does not impute unmeasured flow antibodies or establish functional stem-cell identity.

Run from the directory containing the `altanalyze3` package:

```bash
python -m altanalyze3.components.visualization.flow_viewer.run --bundle /path/to/bundle --port 8085
```

The current lab bundle is `/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/FlowData/rna2flow_viewer`, served at http://127.0.0.1:8085.

Choose a CITE-seq space, then an annotation and one or several populations. Search accepts `ML1a`, `ML-1a`, `MultiLin1a` or `CLP1a`. Use **Highlight cells** to inspect the population on different embeddings. Change the color menu to an ADT or a `RNA:` gene; turn off **selection overlay** to see expression colors. Use **dot size** in the toolbar to choose Auto or 1, 3, 5 or 7 pixels; the setting applies to biaxial and embedding plots and stays selected when changing spaces. Highlighted cells retain a minimum visible size. Actual marrow CITE-seq is in `cite_marrow_ADT195` and `cite_marrow_ADT112`; Grimes and Chinese thymus cells are separate spaces. The coordinate-only marrow reference and the reference with measured RNA are separate too.

Choose a transfer method and **Optimize flow gate**. Gates are fitted on measured flow events carrying the selected transferred annotation. Training-only marker associations rank both positive and exclusion markers; an exhaustive search over quantile-grid rectangles compares marker pairs with short decision-tree paths. Training events fit thresholds, validation events choose the strategy, and untouched test events score the winner. Search compares six ranked channels, 24 quantile bins and tree depths 2–4 with retention-constrained greedy cuts across all available surface markers; “optimized” means best among these candidates, not a proof of a globally optimal gate. Choose whether to balance capture/purity or emphasize either one.

For a specified sorting hypothesis, enter comma-separated signs such as:

```text
CD4-, CD8-, CD117+, CD25-, CD11c-, CD11b-, CD27+, Sca-1+
```

The **Use ML / CLP gate example** button fills that example. Each sign fixes a threshold's direction; five coordinate-update orders are compared. The signs do not supply instrument-defined positive/negative cutoffs. Thresholds are optimized on annotation membership and require experimental controls before interpreting low signal as negative. Edit the displayed threshold values to explore alternatives. Composition updates immediately; fitted test metrics refer to the original winner and are hidden after editing.

After optimization, blue marks selected flow events, orange marks the transferred target, green marks their overlap and gray marks the remainder. Switch between native flow embeddings and marrow projections to inspect the same selected events. All conditions apply to sequential gates; the two-axis outline shows only the current axes. **Export gate JSON** includes thresholds, fitting/test metadata and per-channel Logicle parameters. **Export event IDs** exports the selected 1-based FCS event indices. Manual FlowJo trees, polygon/rectangle/ellipse/quadrant drawing and parent percentages remain available.

**Compare measured ADT / PU.1 / IRFs** reports within-panel percentiles and native quartiles for actual marrow and thymus CITE measurements. PU.1 is `RNA:Spi1`; `RNA:Irf4` and `RNA:Irf8` are included along with other IRFs. No protein reporter identity is assigned to tdTomato. **Compare marrow / thymus RNA** reports descriptive RNA medians. Capture libraries and cells are not treated as independent biological replicates.

**Compare annotations** shows ARI and dominant label overlaps on the same events. **Review transfer validation** shows antibody holdout results. Complete population tables retain missing/too-small states. Same-marker MarkerFinder agreement is a consistency check; marker holdout is a separate cross-marker test. Neither establishes biological replication or functional stem-cell identity. The shuffled gate control permutes test membership for a fixed gate, not the entire transfer/fitting pipeline.

RNA HVG UMAPs are generated from RNA features without annotation labels. Existing RNA MarkerFinder/cellHarmony and scTriangulate embeddings retain their provenance. Approximate marrow placements call the existing `approximate_umap` implementation: assigned labels choose reference coordinates plus jitter. They are explicitly identified as label-conditioned placement and cannot independently validate those labels. Native flow UMAP uses measured surface channels. Original reference coordinates absent for a panel are reported rather than invented.

To rebuild the additional views on an existing event/gate bundle, use its `reference_config.json`:

```bash
python -m altanalyze3.components.visualization.flow_viewer.precompute_reference --config CONFIG.json --bundle BUNDLE
python -m altanalyze3.components.visualization.flow_viewer.precompute_additional --config CONFIG.json --bundle BUNDLE
python -m altanalyze3.components.visualization.flow_viewer.precompute_marrow --config CONFIG.json --bundle BUNDLE --flow-npz FLOW.npz
python -m altanalyze3.components.visualization.flow_viewer.precompute_transfers --bundle BUNDLE --flow-npz FLOW.npz
python -m altanalyze3.components.visualization.flow_viewer.precompute_cellharmony --bundle BUNDLE --flow-npz FLOW.npz
python -m altanalyze3.components.visualization.flow_viewer.validate --bundle BUNDLE
```

`FLOW.npz` contains the event-aligned RDS input matrix `X` and antibody names `channels`. The current copy is `validation_20261002/flow_rds_input.npz` in the bundle. Input conversion and general benchmarks support FCS, FlowJo RDS and event-by-channel CSV; label tables join on event identifiers. Shape agreement alone is not treated as an event join. Install `components/rna2flow/requirements.txt` for processing; pyInfinityFlow supplies the measured FCS Logicle transform. Its use here does not imply InfinityFlow imputation.

The added `kde_cellharmony` method ports the notebook's five assignment stages: Louvain communities, reciprocal community matching, Pearson-nearest reference cells, query prediction centroids and final Pearson-nearest centroids. Correlation matching uses bounded batches rather than a full event-by-reference correlation matrix. Tests compare dense notebook calculations with frozen communities. The default uses igraph multilevel Louvain and adapts PCA dimensions to the shared antibody count. Explicit published-compatible KDE integration, tie ordering and vtraag community options are also available; see the original-code audit for the matched-input parity results.

Browser verification uses actual Chrome, checks pixels for every embedding and PU.1 expression, selects an actual marrow population, fits its constrained flow gate, keeps its selection across embeddings, and checks gate/event exports:

```bash
python -m altanalyze3.components.visualization.flow_viewer.render_check --out RENDER_DIR
python -m altanalyze3.components.visualization.flow_viewer.render_population_check --out RENDER_DIR/population_workflow
```

The reproducible ML/CLP example writes measured ADT/RNA profiles, threshold hypotheses and FlowSOM placement tables:

```bash
python -m altanalyze3.components.rna2flow.run_population_query --out QUERY_DIR
```

Source-specific ML/CLP and DN queries can be reproduced with:

```bash
python -m altanalyze3.components.rna2flow.run_source_comparison \
  --bundle /path/to/rna2flow_viewer --out /path/to/comparison \
  --url http://127.0.0.1:8085
```

This keeps Chinese DSB and Grimes TotalVI sources separate, reports missing source labels, and exports complete FlowSOM distributions, source-cell DN/marrow associations, transferred DN/marrow associations on the same flow events, exact shared-antibody mappings, and fixed-polarity gate scores. Select `StJude` in `cite_grimes` for its broad DN annotation, or `Author_celltype` in `cite_chinese` for DN1–DN4 and gamma-delta labels. Both now offer `kde_cellharmony`; scTriangulate retains its separate DN subsets and alternative transfer methods. On the flow space, compare any two transfer label sets on the same events using the annotation comparison controls. The existing RDS-based transfers use 16 shared markers for Grimes and 14 for Chinese, excluding CD4/CD8. All four CITE panels contain CD8a; the profile export explicitly aliases flow `CD8` to CITE `CD8a` while keeping CD8b distinct.


Enrichment before isolation
--------------------------

The population panel now accepts an **Existing flow parent** (for example the imported
FlowJo Live Cells or a progenitor branch), followed by **Fixed parent thresholds**.
Enter one measured condition per line, using `>` or `<=` and a threshold in the
flow display scale. A Lin cocktail exclusion, CD34 enrichment, endothelial-marker
exclusion or CD117 gate can be used when its channel is measured. No marker is
silently inferred or omitted: the current J8DW panel has neither CD34 nor a dedicated
Lin channel. Do not remove endothelial cells solely by annotation and call that
an executable sorting gate; choose measured exclusions and verify their specificity.

Parent conditions are applied in sequence and retained in the exported strategy.
The audit reports retained events and target recovery after every stage. Splits
are established on the full starting population before filtering. Thresholds fit
within the parent, while validation selection and final test recovery include
all targets lost to parent exclusions. Conditional recovery within the parent is
reported separately. A parent that eliminates nearly all targets fails clearly.
Existing parent masks apply consistently to highlighting, composition and event exports.

The automatic search additionally compares intermediate greedy strategies:
each cut improves training purity while retaining at least 95% of the current
targets; strategies stop at 12 cuts or below 50% retained target cells. The marker
budget is configurable (2–20). Validation, rather than test outcomes, chooses
among greedy, rectangular and short-tree strategies. The returned
`pareto_candidate_indices` identifies validation tradeoffs between purity,
recovery and distinct sorting markers. These are predictions against selected
annotations, not experimental sorting purity.

Paper-specific validation remains essential. Ferchen et al. (local manuscript
`/Users/saljh8/Downloads/nihms-2119165.pdf`, pages 5–6 and 19–20) combined historical
parent-gate comparisons, refined conventional panels, cross-sorter checks and
flow-sort-scRNA-seq remapping. It progressed from a CD34-based capture to a
Lin−Sca1−CD117+CD27+ MultiLin gate; CD34 enrichment is therefore configurable,
not a universal requirement for qHSC or every marrow state. The paper's
Ab-MarkerFinder procedure and default 90% purity/50% yield are distinct from the
older notebook's 95%/95% greedy settings. The new generic greedy search scans
channels directly; it does not claim to be standard MarkerFinder. Any MarkerFinder
analysis must use the adjacent visualization SKILL.md public workflow.

### Rare-state isolation and serial virtual experiments

**Capture is optional.** Choose a published capture template or an existing workspace parent when appropriate for the target. The two paper templates describe the initial Lin−Kit+CD34+CD115−Ly6C− capture and the refined Lin−Sca1−Kit+CD27+ MultiLin capture; they require control/workspace thresholds. Missing or ambiguous antibodies are reported, and their gates cannot be silently substituted. Annotation exclusions are not physical sorting gates. Broad capture is held fixed while finer state gates are discovered.

Select **CITE → measured flow** to fit flow channels against transferred annotations, or **CITE-only ADT gate discovery** to fit measured ADTs against any source annotation. RNA channels are excluded from automated gate discovery. CITE-only thresholds preserve the source normalization (e.g. DSB versus TotalVI); they require calibration before application on an instrument. Unique antibody aliases resolve to measured panel names; ambiguous mappings require an explicit channel choice.

The search now includes training-only coordinate relaxation, pruning and forward expansion, seeded by rectangles, trees and retention-constrained gates. Relaxing earlier cuts after subsequent exclusions can recover targets that the first cut removed. Target-distribution quantiles supplement whole-population quantiles to resolve rare-state tails. The optional **official HyperGate** bridge calls the actual R package; no substitute HyperGate implementation is used. It is enabled by default in the viewer and optional in the Python API (`include_hypergate=True`). Rscript and the official `hypergate` package must be installed; `RNA2FLOW_RSCRIPT` can specify the executable. HyperGate gates exceeding the requested marker budget are pruned using training scores and named `HyperGate_budget_pruned`. Maximum markers applies to discovered child markers; fixed capture markers are additional.

Set minimum validation purity and recovery when a particular isolation requirement matters. Validation chooses a gate meeting both requirements if one exists; otherwise the result explicitly reports that no candidate met them. These constraints do not establish experimental sorting purity, and no global optimality is claimed.

After discovery, **Open serial isolation workspace** opens `/isolation` in a separate browser window/tab. This view provides:

- Capture and child steps in a selectable sequence; the plot shows only cells entering the selected stage.
- Target retained, target lost at this step, other cells retained, and other cells excluded, with circular dots. Choose one-marker serial steps or pair markers into biaxial steps; a two-marker solution can then be inspected as one gate. A one-marker search budget supports simpler isolations.
- Purity (target / retained cells), total recovery (retained target / starting target), step retention, target loss, enrichment relative to starting prevalence, and retained-state composition.
- Editable thresholds, step reordering/removal, and rectangles drawn on incoming cells to insert another serial virtual isolation.
- Validation-only candidate comparisons; the original test score is shown only for the original winner. Manual edits/reordering invalidate inherited test-score claims and show descriptive all-event counts instead.
- Export/import of the complete serial experiment, including source space, labels, capture parent, stage definitions and provenance.

Plots use a reproducible display sample of up to 800 points per target/pass category; counts use all cells/events. This sampling is for visibility and is not a density estimate. Cells/events are not biological replicates. Reordering AND conditions changes intermediate losses but not the final intersection. The workspace does not discover unions of disconnected gates or automate polygon search; it supports serial threshold/rectangle experiments.

For a new ADT-only AnnData, create a CITE-only bundle or append a source space to an existing flow bundle:

```bash
python -m altanalyze3.components.visualization.flow_viewer.precompute_cite \
  --adt-h5ad measured_ADTs.h5ad --bundle BUNDLE --name cite_experiment \
  --labels marrow_state author_cluster --scale "DSB normalized ADT"
python -m altanalyze3.components.visualization.flow_viewer.run --bundle BUNDLE
```

The input must contain ADT features in `X` (or the explicitly selected `--layer`), with cell-state annotations in `obs`. Existing cell-aligned two-dimensional coordinates in `obsm` are preserved; missing embeddings are not invented. Cell order, ADT values and normalization units are retained, and a duplicate space requires explicit `--overwrite`. Future prediction on a new CITE dataset using learned cross-study rules is not implemented: the present CITE-only mode discovers gates within the supplied measured ADT dataset.

Development benchmarks and interactive checks from 2026-10-03 are stored in the lab bundle at `validation_20261003/rare_state_isolation`. The original benchmark is reused for method development, so held-out events are withheld from fitting and candidate selection, but this is not a new external validation cohort. The rare-target stress test subsamples each target to 100 cells among 22,000 background events; there are only 20 target events in each test split. It tests prevalence sensitivity, not new biological samples. Functional identity, conventional sorting and sorted-cell RNA mapping remain separate validations.
