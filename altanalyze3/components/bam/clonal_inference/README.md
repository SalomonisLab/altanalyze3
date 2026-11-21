# Multi-event clone inference

This workflow discovers and reports **combinations of genetic events** from an already consolidated cell-by-event matrix. Direct calls, accepted passenger/CNV/LOH predictions, and accepted expression/isoform predictions can enter as binary observations. Their provenance and upstream conflict rules must be established before fitting. The workflow preserves missing entries in the input matrix, never creates artificial read counts, and never removes events/cells to make a model run. For the default reporting policy, cells with no accepted mutant call are displayed as WT, including entirely uncalled cells; this does not fill missing event genotypes.

**The default is now `detection-tree`, seed 222, following the user's visual review and explicit deployment approval on 2026-11-16.** The production core reproduces the approved `detection_h2_s222` fit exactly, including tree, genotypes, weights and every cell posterior. This is a user-reviewed model, not independent biological validation. Earlier engines and their outputs remain available as labeled historical comparisons.

The method discovers a shared gain-only event tree de novo. It does not encode SRSF2, RUNX1, BCOR or expected clone frequencies in its objective. Each descendant inherits ancestral mutations; specimen-specific mixture weights describe sampled genotypes on that tree. Direct ALT/total reads use beta-binomial likelihoods (mutant mean 0.5 for autosomes, 0.99 for male chrX; mutant concentration 2; reference error mean 0.01, concentration 100). Covered cells calibrate the specimen/event detection rates of accepted passenger/CNV/LOH evidence, with pooled shrinkage strength 2; direct reads are not counted again as independent surrogate evidence. Other accepted categorical calls retain fixed FP 0.01/FN 0.30 contributions as in the reviewed model.

The unchanged search initializes a full event tree, optimizes reattachments/label swaps, then uses BIC to collapse sampled genotype states without removing events or cells. Numerical simplex optimality and gain-only ancestry are gated. The approved default reporting policy (`no-mutant-call-WT-v1`) assigns every cell with no accepted direct or imputed MUT call to WT, including entirely uncalled cells and LOH-negative-only cells. For mutation-positive cells, original posterior memberships among malignant clones are divided by their malignant-clone sum, with zero WT membership. Clone frequencies sum these reporting memberships and divide by ALL cells in the specimen. Source heatmaps use maximum reporting membership assignments and can differ from fractional frequencies. Original calls, missing entries, native posteriors and the fitted tree are preserved. This policy is shared by new runs and rendering saved fits; no unresolved category is plotted. Additional starts 111 and 333 measure initialization sensitivity; they do not replace the reviewed seed. Sparse calibration and correlated surrogates remain limitations, and gain-only inheritance is an explicit assumption.

The prior `supported-combinations` workflow is a descriptive audit of exact observed MUT sets, not an ancestry model. The earlier `shared-genotype` workflow lacked ancestry and could overextend sparse accessory events. Earlier architecture checks are historical development diagnostics; they do not supersede the user's subsequent approval of the displayed detection-aware result. Development against known 5801 biology is not an independent validation set.

## Reviewed default: inference and complete figures

```bash
python -m altanalyze3.components.bam.clonal_inference \
  --matrix inputs/calls.npz --cells inputs/cells.tsv --events inputs/events.tsv \
  --contract inputs/contract.json --detector inputs/detector.npz --sex M \
  --output results/patient_seed222
python -m altanalyze3.components.bam.clonal_inference.detection_plotting \
  --result results/patient_seed222 --input inputs --output figures \
  --sensitivity-root results
```

The default requires direct `A`/`T` count matrices, an independently reconciled detector NPZ (`Z`, boolean `available`, exact ordered `cell_ids`/`features`), chromosome/event labels and sex M/F/U (U only when no chrX event exists). The contract must match `detection_workflow.specification`, detector hash and recorded authorization. No synthetic read counts or categorical fallback is permitted. Submit substantial analyses and rendering through the institutional scheduler.

The figure command additionally requires `source_evidence.npz`, `evidence_rows.tsv`, and `plot_source_contract.json` (source hashes and verified provenance). It exports editable-text, vector-only SVG, PDF and PNG for composition/frequency heatmaps, complete cell/source heatmaps, mutation-dot cartoons, actual inferred ancestry, nested fishplots, combined fishplot/cartoons, specimen bars and optional seed sensitivity. It verifies native baseline or current reporting tables against the complete roster before writing, then applies the same reporting policy to all panels and all sensitivity starts. MISTRG specimens are separate denominators, never later clinical stages. Fishplot smoothing interpolates sampling stages, not elapsed calendar time. Clone keys and full-precision tables accompany the figures. The wide fishplot canvas is 15.86 × 7.18 inches, with a cream/gold/brown WT background and explicit WT percentages at each clinical timepoint. Its Matplotlib rendering adapts Xuan’s Python plot only: no VAF values, CCF conversion or assumed hierarchy are imported. An explicitly schematic prediagnosis lead-in starts at all WT and stages sampled ancestral genotypes before reaching the exact first-sample frequencies; it is not inferred historical burden or timing. The clone_fish.py renderer preserves every clinical frequency knot and the all-cell denominator, without spline clipping or renormalization.

Clonal-tree panels use Arial (20-point clone IDs, 16-point italic mutation labels,
34-point titles), straight 0.5-point connectors with 18-point arrowhead scale, and unnumbered mutation dots of radius 4.68
points centered in up to three rows. A colored, italic mutation legend retains
every panel event. Install Arial and Arial Italic, or set `ALTANALYZE_ARIAL_DIR`
to a directory containing your licensed TTF files; the renderer refuses silent
font substitution. SVG text remains editable and PDFs embed the fonts. Other
figure types retain their existing typography and numbered mutation-dot keys.

`phylogeny_validation.py` supplies tested gain-only tree and frequency checks. It accepts multi-mutation gains, rejects cycles, recurrent gains, ancestral mutation losses and missing required founders, and distinguishes exclusive clone frequency from cumulative lineage frequency. This is a structural validator, not an inference engine or evidence that an old fit is corrected. No legacy engine has been silently routed through a new model.

## Reporting output and compatibility

`model.npz` always retains the native fitted posterior. `reporting_memberships.npz` contains the explicitly labeled, policy-adjusted memberships used for figures; `reporting_policy.json` records the rule. In `clone_frequencies.tsv`, `reported_percent_all` is the fractional reporting percentage and `percent_all` is the maximum-membership cell percentage. `posterior_expected_percent_all` remains a compatibility alias of `reported_percent_all`, not an unmodified native posterior estimate. Likewise, assignment `max_posterior` aliases `max_reporting_membership`. Native baseline run directories are preserved when rendering previously saved fits; corrected assignments, frequencies and evidence audits are written into the figure output directory.

A new fit requires a contract generated from the current `detection_workflow.specification(...)`, including this reporting policy. Old contracts are not silently rewritten. Saved legacy fits may be rendered after their original tables, hashes and identities pass validation; fitting and reporting provenance remain distinct. WT is the operational no-mutant-call category requested for these figures, not proof of genome-wide wild type.

## Engines

| Explicit method | Mechanism | Important limitation |
|---|---|---|
| `detection-tree` (reviewed default) | Calibrated direct/surrogate likelihoods, full gain-only tree search, specimen-specific genotype mixtures; seed 222. | RNA detection, source correlation and sparse calibration can affect clone identifiability; initialization ranges are not confidence intervals. |
| `supported-combinations` (descriptive audit only) | Exact accepted MUT-set grouping; full ternary calls retained; no additional imputation or prevalence extrapolation. | Different call combinations can be partial observations of the same biological clone. Does not infer founder/descendant relationships and is not a genetic clone solution. |
| `shared-genotype` (legacy comparison) | Shared discrete clone genotypes, specimen-specific weights, event-specific learned errors bounded at FP≤0.10 and FN≤0.30; nested masked-call selection of component count and multiple initialization starts. | Conditional independence and bounded errors are modeling assumptions, not measured assay calibration. Identical genotype components are merged for reporting. Exact genotype assignments can remain unstable between starts. No ancestry is inferred. |
| `bnpc-native` | Pinned BnpC 0.2.1 Dirichlet-process mixture, learned FP/FN, beta(.75,.75) cluster parameters; random or single-cluster initialization; free assignments throughout. | Preserves the released gamma-parameterization problems for baseline reproduction. Do not interpret this baseline as a validated result. |
| `bnpc-gamma-corrected` | Same pinned class and sampler, with exactly three corrected gamma expressions in an isolated source copy: prior `loc=0, scale=1/b`; both concentration draws use `scale=1/(b-log eta)`. | This is a separately labeled correction, not unmodified BnpC. Native concentration clipping at >=1 remains. MAP genotype parameters are not cell posterior probabilities; multiple seeds and adequate sampling are required. |
| `exhaustive` | Enumerates all 2^J binary genotypes. For each sample, maximizes `sum_i log(sum_g w_g product_j P(x_ij|g_j))`, omitting unknown entries. Fixed explicit FP/FN; no supplied exome clone configurations. | Custom method, **not native scBayes**. Error rates and conditional missingness assumptions are not automatically calibrated. Unlinked events can make complete configurations unidentifiable. No ancestry is inferred. Limited to 16 events; larger panels fail explicitly rather than being reduced. |

In the exhaustive method, identical observed patterns are compressed with multiplicities, preserving the exact likelihood and every cell. Optimization uses L-BFGS-B followed, when necessary, by likelihood-monotone Frank-Wolfe/EM polishing. A simplex directional-gap bound must be <=1e-5 per cell. Family bounds solve linear programs preserving every observed-pattern likelihood; original constraints are checked after removing numerical row redundancy. These are **likelihood-equivalence ranges under a fixed error model, not confidence intervals**.

BnpC source: commit `edb1ddb1beb1bf977d44aede4c50e6bf985e06a8`, supplied with `--bnpc-source`. Every pinned file is hashed before execution. The source remains unchanged; an output-local copy retains the original MIT license and an exact patch. Do not subclass `CRP_errors_learning`: native code dispatches error learning by `__module__`. This adapter binds the initialization method to the unchanged native class and tests class identity, deepcopy binding, free assignments and error-learning dispatch.

## Shared-genotype algorithm and model selection

For observed calls, the likelihood is
`P(x_i | specimen s) = sum_k w_sk product_(j observed) Bernoulli(x_ij; p_kj)`,
where `p_kj = FP_j + (1-FN_j-FP_j) G_kj` and `G_kj` is binary. Missing entries have likelihood contribution one. Direct and accepted imputed observations enter the consolidated target-event matrix once; passenger or expression evidence is not counted as another independent copy of that target event. The algorithm does not use treatment labels, expected drivers, manually supplied clone combinations, or expected frequencies to initialize or select fits. Specimen IDs index mixture proportions only.

Identical specimen/call patterns are compressed with multiplicities, preserving every cell and event. Initialization uses 25 soft-EM steps, then coordinate EM updates binary genotypes, weights and event errors. Weights use Dirichlet(1.5) regularization; error parameters use Beta(1.5,1.5), truncated to `[.001, FP_cap]` and `[.001, FN_cap]`. The penalized objective must be monotonic and improve by less than `1e-7` per cell to converge. A nonconverged required fit stops the workflow rather than being dropped.

Default component candidates are `2,4,6,8,12,16,24,32,48,64,96,128`, limited by the mathematical bound `2**J`; no input event is removed. For each of three starts, 20% of observed calls are held out externally, then 20% of remaining calls internally. Features sharing `validation_group` are masked together within cells (e.g. RUNX1 mutation and LOH), reducing direct cross-event leakage. The minimum mean inner log loss selects K; exact ties favor smaller K. The outer mask reports predictive performance, not genetic ground truth. The representative full-data fit has the largest training objective across starts. Full data and seeds, outer scores, component count curves, duplicate-genotype merging, and all alternative-start genotype frequencies are retained.

If `selection_at_upper_boundary` is true, extend `--k-grid` with an updated explicit algorithm contract and a new output directory before treating K selection as settled. This leaves the original run intact. The MDS-U2AF1-RAS analysis required an extension through its complete 1024-genotype space. K counts mixture components before identical genotype components are merged; it must not be reported as the number of biologically distinct clones.

FP/FN caps are operational sensitivity settings. In the 5801 comparison, the 10%/30% bounds retained explicit genotypes with similar held-out loss to broad bounds while preventing large false-positive rates from discounting accepted mutation calls. This does not prove the caps are calibrated. The specimen/event marginal predictor still slightly outperformed the chosen model on combined withheld direct and imputed calls. Corrected LACE had substantially worse direct-call prediction and extensive equal-likelihood placements. Full genotype stability, source-specific performance and contradictions are therefore reported alongside point estimates.

## Inputs and contract

- NPZ keys `X`, `cell_ids`, `features`. `X` is cells by events: `1=MUT`, `0=WT`, `-1=unknown`.
- Cell TSV: unique `cell_id`, `sample`; optional `stage`, `cell_state` and source metadata.
- Event TSV: unique `feature`; other event/provenance columns are preserved.
- JSON contract from the independently reconciled primary roster. It must contain `schema_version:1`, the explicit `method`, exact ordered `expected_cell_ids`, `expected_features`, and per-cell `expected_samples`; `input_sha256` for `matrix/cells/events`; nonempty `provenance_sources` records with `path/sha256`; and a description of `source_precedence`.
- Contract `algorithm` must exactly match the settings returned by `algorithm_spec` in `__main__.py`. This includes the seed and all sampling/error settings. Establish this specification from the intended baseline before execution; do not generate an expected roster by intersecting a candidate matrix with metadata.

Any identity, order, source hash, method, parameter, or missing-source mismatch stops before fitting. The CLI refuses existing output directories. Failed runs retain diagnostic artifacts without a `COMPLETE.json` marker. There is no skip/force gate option.

Optional `--families` JSON maps a reporting label to `{exact feature ID: 0 or 1}`. All features still enter the model. These definitions are used only **after fitting**, never as initialization, clone-count constraints or expected evolutionary history.

## Run

Legacy shared-genotype comparison (explicit method required; run each patient separately):

```bash
python -m altanalyze3.components.bam.clonal_inference \
  --matrix calls.npz --cells cells.tsv --events events.tsv \
  --contract shared_contract.json --method shared-genotype --output shared_seed5801 --seed 5801
```

The contract must name `shared-genotype` and match all defaults or explicitly supplied overrides. `--fp-cap`, `--fn-cap`, `--starts`, `--k-grid`, `--max-iterations` and `--validation-seed` are available and recorded. Unrelated engine-specific overrides are rejected. Optional cell fields `patient`, `stage`, `stage_order`, `is_xenograft` support clinical ordering; the shared workflow rejects multiple patient IDs in one fit. Event fields `event_label` and `validation_group` support biological labels and grouped masking.


Install the normal AltAnalyze3 requirements and `Bottleneck==1.3.8` for BnpC (`pip install -e '.[clonal-inference]'`). Native source provenance is checked separately; the workflow does not silently install or download an alternative BnpC revision.

```bash
python -m altanalyze3.components.bam.clonal_inference \
  --matrix calls.npz --cells cells.tsv --events events.tsv \
  --contract corrected_bnpc_contract.json --output bnpc_seed5801 \
  --method bnpc-gamma-corrected --bnpc-source /path/to/pinned/BnpC \
  --initialization together --iterations 1000 --burn-fraction 0.33 --seed 5801
```

```bash
python -m altanalyze3.components.bam.clonal_inference \
  --matrix calls.npz --cells cells.tsv --events events.tsv \
  --contract exhaustive_contract.json --output exhaustive_seed0 \
  --method exhaustive --false-positive-rate 0.10 --false-negative-rate 0.10 \
  --seed 0 --families reporting_families.json
```

For large cohorts, submit these commands through the institution's scheduler (`bsub` on the audited system). Set a sufficient wall-time; iteration count is not a convergence criterion. The 5801 audit reproduced a NumPy/OpenBLAS matrix-product error on one node's default backend. `OPENBLAS_CORETYPE=Haswell` resolved that specific discrepancy. Both startup and per-fit numerical checks remain mandatory; do not assume a backend environment setting is sufficient without passing them.

## Outputs and interpretation

The legacy shared-genotype engine exports `clone_genotypes.tsv`, `clone_frequencies.tsv`, `cell_assignments.tsv.gz`, `cell_event_posteriors.npz`, `learned_errors.tsv`, `start_stability.tsv`, `all_start_genotype_frequencies.tsv`, nested validation masks/scores, and every selected-K raw fit. Exact genotype profiles have stable SC identifiers within a patient and explicit mutation-combination labels. An event unobserved everywhere is retained and reported unknown, not arbitrarily WT or MUT. Local clone/event coverage is also flagged. Cells with no observations are labeled `NE`; raw fits retain prior-driven probabilities, which must not be interpreted as measured evidence.

`percent_all` is 100 × MAP-assigned cells / **every cell in the specimen**, including predicted WT and unresolved cells. `posterior_expected_percent_all` provides a separate probability-weighted summary with the same denominator. Start min/max fields track the exact genotype across starts; they are optimization sensitivity ranges, not confidence intervals. Genotypes found only in alternative starts remain in the all-start table.

Editable publication figures:

```bash
python -m altanalyze3.components.bam.clonal_inference.plotting \
  --result shared_seed5801 --matrix calls.npz --output publication_figures
```

Figures include genotype/frequency heatmaps, complete profile heatmaps, initialization-sensitivity heatmaps, stacked specimen compositions, and side-by-side serial cell evidence. SVG contains editable text and vector geometry, with no raster images; PDF/PNG previews and legends accompany each figure. MISTRG stays a separate specimen beside its matched clinical stage. Top-12 display pooling never changes fits or denominators; the complete plot and tables retain all profiles. Optional NPZ `direct`, `imputed_MUT_sources`, `imputed_WT_sources`, and `source_method_order` enable source-separated evidence rows and retain imputation conflicts. Plotting checks the fitted matrix hash, exact IDs, assignment-derived counts and all-cell percentage totals.


Every method exports the full cell/event rosters and a source/algorithm manifest. BnpC exports its native results, traces, MAP state, complete genotype parameter matrix, raw clone IDs, full thresholded genotype groups, sample frequencies, cell-state counts and per-clone/event MUT/WT/unknown evidence counts. A .5 parameter threshold is a display convention, **not a confidence threshold**. An event with no observations in an assigned cluster is explicitly flagged even if the model supplies an extreme parameter.

BnpC also exports saved-state mean genotype parameters, threshold-positive fractions, predictive call probabilities and the fraction of states supporting each cell’s complete MAP threshold vector. Every retained state must pass an independent likelihood check; native parameters are already post-burn-in and are not sliced twice. These support fractions quantify sampling uncertainty, not calibrated biological confidence. Cell/event IDs accompany all arrays.

The exhaustive engine exports all configuration weights, continuous cell/event posterior probabilities, full-configuration assignment probabilities, family frequencies and equivalence bounds. Sample-specific mixture weights estimate a population distribution under the conditional observation model. Expected counts from informative cells are reported separately. Cells with no evidence keep their native/prior-driven calculation but are displayed as `Unresolved_no_evidence`; they remain in all-cell denominators.

Assess seed agreement using the **same full genotype vectors**, not just label-invariant partition ARI. Inspect source-specific contradictions, event coverage, posterior ambiguity and held-out predictions. A held-out accepted imputed call is not independent genetic truth if its upstream classifier used related expression. Cell-state associations involving expression-derived calls are descriptive and may be circular. Predicted panel wild type never establishes genome-wide wild type or a normal cell. Connecting frequencies between specimens is not evidence of cellular tracking or clone ancestry.

Tests:

```bash
BNPC_SOURCE_PATH=/path/to/pinned/BnpC python -m unittest discover \
  -s altanalyze3/components/bam/clonal_inference/tests -v
```

Tests cover known-mixture recovery, non-identifiability with disjoint coverage, missing-data invariance, the full 11-event numerical kernel, source/identity/feature gates, denominator preservation, native-class sampling, MAP likelihood reconstruction, and source preservation.

## Fish-style courses and clone cartoons

The publication plotting command also exports `<patient>_fishplot_cartoons`,
`<patient>_fishplot`, and the complete `<patient>_clone_cartoons_<page>` catalog
as SVG/PDF/PNG. Ring colors match the clone colors in the specimen bars;
numbered internal dots identify predicted mutant events. Full mutation names
and an event key accompany the cartoons. No ancestry arrows are fabricated.

Frequency bands use convex cosine joins between exact specimen percentages.
They preserve nonnegative proportions and a 100% all-cell denominator. Clinical
stages are equally spaced schematically, not placed on an inferred elapsed-time
axis. The unresolved fraction is shown as gray space split equally above and
below the centered resolved profiles. MISTRG samples are separate hatched bars,
not continuations of the clinical patient trajectory. Every rare profile is
retained in the unpooled cartoon catalog and frequency tables.

To add only these panels to an existing completed result:

```bash
python -m altanalyze3.components.bam.clonal_inference.fishplot \
  --result shared_seed5801 --output publication_figures
```

This is a visualization of existing assignments, not a new lineage fit.

## Descriptive evidence grouping, withdrawn as a clone default on 2026-11-16

The former shared-genotype default is withdrawn for biological clone reconstruction.
It remains available explicitly as a comparison. Its complete genotype assignments
could propagate a sparsely supported event to most cells via other shared events.
Neither optimization convergence nor aggregate prediction scores validate that step.

`supported-combinations` remains an explicit deterministic evidence grouping, not a default or a validated clonal reconstruction:
- Each cell's accepted MUT set defines its group. All co-mutations are retained.
- No new MUT calls, no lost accepted MUT calls, and no unknown-to-WT conversion.
- Every input event and cell remains in the contract and all-cell denominator.
- WT = panel-complete WT calls; WTobs = WT-only but incomplete panel;
  NE = no calls. These are not equivalent biological categories.
- SC identifiers label supported mutation combinations, NOT proven distinct complete
  clones. Two different groups may be the same biological clone with different
  coverage; this method does not estimate unseen genetic prevalence.
- Uncertain inherited-status markers remain represented and explicitly annotated;
  including an event does not validate somatic specificity.
- No expected patient trajectory, seed selection, fitted error cap, clone scaffold,
  cell-state grouping, or bulk-VAF multiplication enters grouping.

`--method shared-genotype` retains the old method under its original explicit
contract. Existing contracts are never silently migrated: the method and algorithm
specification must match. New supported contracts retain original source/roster
hashes and document the requested no-additional-extrapolation change.

Plot a completed supported result using the regular plotting command, or:
`python -m altanalyze3.components.bam.clonal_inference.supported_plotting --result RESULT --output FIGURES`

Exports include fish-style frequency courses, mutation-dot cartoons, complete
combination/frequency pages, and all-cell ternary evidence heatmaps as editable SVG
and PDF with PNG previews. WT and unresolved cells are included in all denominators.
Curves interpolate sampled frequencies and do not imply ancestry or tracked cells.
`event_frequency_audit.tsv` proves exact equality between accepted positive counts
and reported positive-group counts for every event and specimen; the execution
gate additionally verifies equality for every cell/event. This guards unsupported
expansion, but does not calibrate the accepted upstream imputers.


## Extended default publication outputs (2026-11-17)

The default detection_plotting command now also exports compact raster source heatmaps and clinical cell-state/clone alluvials. Pass --alluvial-config for explicit manuscript state aggregation; unknown states and the complete clinical roster remain. --legacy-vector-only preserves the previous standard-figure invocation. Core inference, source calls, WT policy and frequency denominators are unchanged. Native AltAnalyze2 UMAP and ring workflows are orchestrated by visualization.clone_figure_workflow; see visualization/CLONE_FIGURES.md.

## Non-mutating lineage cross-check (2026-11-17)

`detection_workflow` additionally writes `minimum_supported_lineage_audit.tsv.gz`.
It compares fitted assignments with the least-specific active genotype containing
all accepted MUT calls and contained in every compatible active genotype. Ties or
incomparable branches remain unresolved (`-1`); no accepted MUT remains operational
WT. Reference RNA conflicts are counted separately. This audit never modifies the
fit, its cell assignments, frequencies, input calls or feature/sample roster.

`lineage_support.supported_lineage` is a diagnostic, not a new clone estimator:
a minimum supported lineage does not prove absence of descendant mutations or
validate sparse-event mutation order. A smallest mutation count alone is not a
valid choice when the compatible branches are incomparable. Dedicated tests cover
that case, missing calls, positive conflicts, WT retention and order invariance.
The 24-fit targeted likelihood-correction experiment was not promoted; the
original 5801 baseline remains preserved. Six-series candidate figures and the
executed Figure4 notebook document the comparison separately from production.

### Named candidate: `final-revised-1` (2026-11-18)

`workflow_versions/final-revised-1.json` identifies the frozen six-series
original-tree evidence cross-check candidate and its complete output/code hash
manifest. This name does not promote the candidate to the default estimator.
It preserves fitted trees and reports the minimum supported lineage; these
counts are not exclusive complete-genotype frequencies. Its 5801 assignments
and frequencies differ from the earlier approved baseline.

The plotting default is `TREE_STYLE_NAME = "5801-approved-compact"` in
`detection_plotting.py`: compact depth spacing, 20-point Arial clone labels,
16-point italic mutation labels, 0.5-point straight shafts, arrowhead scale 18,
and centered unnumbered mutation dots of radius 4.68 points. All patients share
these defaults. Wrappers must not enlarge the canvas or reduce label fonts;
long HGVS labels wrap to avoid collisions without shrinking the typography.
