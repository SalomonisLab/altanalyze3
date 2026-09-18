# Plot semantics review — 2026-09-16

Reviewed the shared scALABLE interface, uploaded-job endpoints, viewer bundle endpoints, interaction and regulatory networks, expression UMAPs and violins, DotPlot/CombPlot, volcano plots, and PDF exports. Existing analysis outputs were preserved.

| Finding | Correction |
| --- | --- |
| Chat volcano plots treated missing FDR as effectively zero, producing an extreme significance score. | Missing or invalid folds/FDR are excluded; an explicitly reported FDR of zero remains supported. |
| The main volcano mixed raw p-values into an axis labeled FDR. | Use valid FDR when available. Use a labeled p-value axis only when no valid FDR is available. The API reports omitted statistics. Screen and PDF use the same label. |
| Uploaded DotPlot ignored the second annotation filter. | Apply both filters with AND. |
| Empty selections sometimes reverted to all cells/groups in uploaded and viewer gene-set plots. | Preserve empty matches and show an empty state. |
| Viewer detection fractions differed between cached and filtered paths for signed/explicit-zero data. | Both paths calculate the fraction of values greater than zero; the hover label states that definition. |
| Unknown gene requests silently substituted the first feature. | Keep the requested gene identity. RNA can use an explicitly labeled matching reference result; otherwise report missing. Only a blank request can select a default feature. |
| Missing interaction-network folds became zero and were then colored red. | Preserve nulls and render missing/zero folds in gray in screen and PDF. Marker folds require positive values. Genuine signed differential folds retain their sign. |
| Negative expression UMAP values were rendered and described as zero. | Preserve the signed values in colors and hover text; only actual zeros use the zero layer. DotPlot also preserves a negative mean color range. |
| TF activity was titled expression, and its server PDF used a different palette. | Label imputed TF activity explicitly. Align server PDF palettes and zero handling with the interface. |
| Empty filtered violin groups polluted group lists and could break PDF export. | Omit groups without observations; export a clear empty-state message when all groups are empty. Apply filters to the scatter payload too. |
| Regulatory edge widths used signed scores despite thresholds operating on magnitude. | Use absolute score for width and show the signed score on hover. |
| Cell communication's Download PDF button delivered SVG for network plots. | Convert the network to PDF using the shared vector exporter. |
| GRN edges and Explore Regulatory network displayed a modality selector that did not control their data source. | Hide it for these views in both panels; retain the user's modality when returning to expression views. |
| GRN edges opened as a dense 300-edge graph without explaining its source or aggregation. | Default to the strongest 25 connections, offer 50/100/300, show matching counts, selected genes/state/sample annotation, and the number of sample-by-state aggregates. Distinguish predicted connections from positive-marker networks. |
| GRN edge widths saturated and small scores rounded to zero in the API. | Preserve score precision; scale widths by magnitude within the displayed graph and expose scores on hover. Outline queried genes and use a concentric layout. |
| Large SVG networks could be clipped in PDF conversion. | Set the SVG viewport to the page dimensions while preserving its original coordinate system. Include the GRN scope and interpretation in its PDF. |

The marker heatmap is standardized per gene (a z-score), so negative heatmap values legitimately indicate expression below that gene's mean. They are distinct from the positive marker fold used for marker selection and network coloring.

For saved job `66ed63fe8a7a49a582ab0f86a68830e3`, all stored RNA, ADT, GRN-edge, TF-activity, and lipid expression values checked were finite and nonnegative. Its 2,274 RNA markers have positive folds, positive correlations, and query expression greater than reference expression. This check concerns the saved outputs; other deployed datasets were exercised with synthetic bundle fixtures rather than reprocessed.

Validation: 29 Python tests passed across plot semantics, individual-cell CombPlot, integration, GRN release, and bundle expression access. Browser checks passed for signed values, neutral missing folds, TF labels, missing volcano statistics, empty DotPlot, both cell-communication panels, and a downloaded PDF with a valid PDF signature. JavaScript syntax and whitespace checks passed.

GRN clarity follow-up (2026-09-17): 20 GRN/integration tests passed, including edge limits, reported totals, aggregation counts, exact mean scores, and invalid limits. Browser checks passed in both panels for hidden modality controls, changing limits, restoring the expression modality, and positive-marker legends. GRN, marker-regulatory, and cell-communication PDF exports passed; the GRN PDF was visually inspected and all 26 node labels were present. The current ZNF143 / Hillock basal example has 266 matching connections across two aggregates; the default graph shows 25.
