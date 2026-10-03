The viewer selects populations in measured CITE-seq data and evaluates their CITE-seq → flow assignments, FlowSOM overlap, marker expression and sorting gates. It does not impute unmeasured flow antibodies or establish functional stem-cell identity.

Run from the directory containing the `altanalyze3` package:

```bash
python -m altanalyze3.components.visualization.flow_viewer.run --bundle /path/to/bundle --port 8085
```

The current lab bundle is `/Users/saljh8/Dropbox/Collaborations/Grimes/Thymus/FlowData/rna2flow_viewer`, served at http://127.0.0.1:8085.

Choose a CITE-seq space, then an annotation and one or several populations. Search accepts `ML1a`, `ML-1a`, `MultiLin1a` or `CLP1a`. Use **Highlight cells** to inspect the population on different embeddings. Change the color menu to an ADT or a `RNA:` gene; turn off **selection overlay** to see expression colors. Use **dot size** in the toolbar to choose Auto or 1, 3, 5 or 7 pixels; the setting applies to biaxial and embedding plots and stays selected when changing spaces. Highlighted cells retain a minimum visible size. Actual marrow CITE-seq is in `cite_marrow_ADT195` and `cite_marrow_ADT112`; Grimes and Chinese thymus cells are separate spaces. The coordinate-only marrow reference and the reference with measured RNA are separate too.

Choose a transfer method and **Optimize flow gate**. Gates are fitted on measured flow events carrying the selected transferred annotation. MarkerFinder ranks both positive and exclusion markers; an exhaustive search over quantile-grid rectangles compares marker pairs with short decision-tree paths. Training events fit thresholds, validation events choose the strategy, and untouched test events score the winner. Search considers six ranked channels, 24 quantile bins and tree depths 2–4; “optimized” means best among these candidates, not a proof of a globally optimal gate. Choose whether to balance capture/purity or emphasize either one.

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

The added `kde_cellharmony` method ports the notebook's five assignment stages: Louvain communities, reciprocal community matching, Pearson-nearest reference cells, query prediction centroids and final Pearson-nearest centroids. Correlation matching uses bounded batches rather than a full event-by-reference correlation matrix. Tests compare dense notebook calculations with frozen communities. This implementation uses igraph multilevel Louvain instead of the notebook's vtraag backend and adapts PCA dimensions to the shared antibody count; it does not claim bitwise community parity.

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
