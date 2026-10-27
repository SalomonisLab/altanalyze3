# Visualization CLI

`python -m altanalyze3 visualize` exposes the existing native ISV renderer and
donor-level PSI plotting without rerunning differential analyses.

```bash
python -m altanalyze3 visualize isv-domains \
  --run_dir /dataset/isv --gene ENSG_ID \
  --render-input selected_native_response.json.gz \
  --render-input-sha256 EXPECTED_SHA256 \
  --protein-features /results/protein-impact-plots \
  --caption-file scientific_caption.txt --figure-mode \
  --out integrated_ISV_domains.pdf

python -m altanalyze3 visualize psi-boxplot \
  --counts native_donor_isoform_counts.tsv \
  --donor-roster complete_donor_roster.json \
  --inclusion-isoforms INCLUSION_ID \
  --exclusion-isoforms EXCLUSION_ID_1 EXCLUSION_ID_2 \
  --cohort-order young aged MDS-pre MDS-pre-SRSF2 \
  --caption-file scientific_caption.txt --y-min 0.5 \
  --out donor_PSI.pdf
```

Saved ISV responses already contain their native query/filter selection; replay
rejects additional selection flags. The existing native junction filter accepts
`'E4-junction + E11-junction'` for an AND clause, with multiple clauses ORed.
No alternative junction-selection algorithm is introduced by these commands.

PSI input is a long TSV with `cohort`, `donor`, `isoform_id`, `assigned_RNA`
columns (column names are configurable). The independent roster JSON maps every
condition to its full ordered donor list. Each donor/isoform count, including
observed zeros, must be present exactly once; unexpected identities and missing,
duplicate, negative or nonfinite counts are errors. Inclusion and exclusion ID
sets must be disjoint. PSI is `inclusion/(inclusion+exclusion)` within that explicit
ID universe. Zero-support PSI remains undefined. Condition groups may overlap;
their donor records are preserved without deduplication across groups.

The PSI command writes a PDF, complete donor PSI TSV and provenance manifest.
Existing outputs are protected. Axis limits affect the view, not the saved donor
roster or PSI values. The manifest records whether every defined point is in view.
Boxes use median/IQR and 1.5-IQR whiskers; all defined donor points are drawn.

For existing junction-event PSI, replace `--counts` and the isoform role flags with
`--psi-table native_event_PSI.tsv --event-uid 'GENE:inclusion|GENE:exclusion'`.
The table must contain `cohort`, `donor`, `event_uid`, and `PSI`, retaining every
roster identity. Empty PSI remains undefined. `--psi-sha256` pins the input.
No isoform selection or PSI recalculation occurs in the plotter.

With `--annotations-table`, blank/NA native P values remain unavailable and are
labelled `NA`; they are not converted to `ns` or replaced with a newly computed
test. All donor roles and saved PSI values remain in the output table.

`--reference-condition aged` adds `ns`, `*`, `**`, or `***` above each other
condition, using the existing visualization Mann–Whitney U test (two-sided).
`--stat-correction holm` is the default across reference comparisons; `bh` and
`none` are explicit alternatives. Raw and adjusted P values, sample counts,
and labels are saved in `.statistics.tsv`. The reference must have disjoint
donors from every tested condition; overlapping comparison groups may both be
tested against that independent reference, without testing them against each other.

### Legacy strongest-pattern iterative MarkerFinder

The existing `components.cellHarmony.iterative_marker_finder` module now has an
explicit legacy mode. Its broad/residual mode and the standard RNA/PSI hook
selection defaults remain available unchanged.

```bash
python -m altanalyze3.components.cellHarmony.iterative_marker_finder \
  --legacy-manifest schemes.json --outdir IterativeMarkerFinder \
  --rho 0.4 --max-missing-fraction 0.3 --top-n 100 --half-range 0.2
```

The JSON manifest contains `display_h5ad`, `display_group_key`, optional
`covariate_columns` and `label_events` paths, and `schemes`, each with `name`,
`h5ad`, `group_key` and `level_order`. H5AD matrices contain pseudobulks by events,
with NaN retained. All schemes must retain the complete ordered feature roster.
The CLI invokes the standard `generate_marker_heatmap_from_adata` hook for each
scheme and exports full signed rho, raw-P and FDR tables, including untestable
source rows. `return_statistics=True` exposes these same full tables to API users.

Legacy consolidation follows AltAnalyze2: strongest up/down pattern across
schemes, strict rho > 0.4 by default, original equal-rho tie rule, up blocks then
down blocks, and descending rho/feature ordering within a block. It applies no
FDR selection cutoff; the FDR tables remain available. This differs explicitly
from the standard hook's FDR-selected marker assignments. The consolidated
figure uses row-median-centred PSI, an explicit ±0.2 scale by default and the
standard renderer. The full assignment table is retained; `top-n` limits only
the figure. Rendering uses 600 dpi in this mode, with vector text and axes.

`assess_pattern_coverage` annotates every assignment with observed PSI counts in
every requested group. Its eligibility flag permits explicit per-group coverage
sensitivity figures without changing source data or correlations. A missing
measurement is never PSI zero. Repeated patient pseudobulks are pattern-discovery
units, not independent biological replicates.

The standard cached renderer accepts `colorbar_label` for an explicitly different
saved normalization, and retains distinct iterative row-pattern labels when they
differ from column groups. Long covariate legends wrap complete entries; dense
callouts fit the vertical extent without changing the heatmap's axis limits.

For event mirror plots, `splice_event_barcharts --strict-thresholds` makes both
cutoffs strict (`P < threshold`, `abs(dPSI) > threshold`). The historical inclusive
boundary remains the default. `--cluster-order` accepts one cell state per line;
nonfinite P/dPSI values cannot qualify under either mode.

`assess_leave_patient_out_patterns` runs the same public hook with each verified
patient held out in turn, removing that patient's repeated state/time-point
pseudobulks together. It exports every original feature for every holdout and
annotates all baseline assignments with literal pattern agreement. Holdout labels
can differ at mathematically complementary up/down floating-point near ties; this
is a sensitivity diagnostic, not independent clinical significance. Coverage
annotations also expose measured group-mean PSI from the native centroid code,
so amplitude criteria can be stated separately from correlation strength.

The native `isoform_state_switches protein-outcome` CLI accepts
`--fold_scale ratio --min_usage_shift 0.1` to select a strict absolute fractional
usage shift instead of a fold ratio. Input means must lie in [0,1]. Existing
expression-fold and usage-fold defaults are unchanged. Literal source DiffMeans
values define the threshold; rounded displayed means are not substituted.

### Immutable-fit PSI figure review

`iterative_marker_finder --legacy-manifest ... --review-existing-fit <fit-dir>
--input-sha256 <exact-h5ad-sha256> --outdir <new-dir>` reviews saved native
assignments without refitting or changing their winning categories. Supply
`--rho 0.8 --top-n 5 --min-observed-per-group 3 --min-group-present 0.9
--min-mean-psi-range 0.1` for an explicit stricter figure sensitivity. Optional
`--patient-stability-table <native-holdout-table> --min-patient-agreement 0.8`
uses recorded whole-patient stability. These are stated study choices, not new
universal defaults or independent validation. Every original assignment and
eligibility flag remains in `all_original_assignments_and_review_gates.tsv`;
small patterns are never padded to the top-N cap. The exact fitted feature and
sample rosters are checked before rendering. `--title`, `--caption-file` and
`--sample-label-key` add scientific context and verified donor labels.
Use `--strict-top-n --no-collapse-pattern-logic` to retain literal saved winning
categories and enforce an exact per-pattern cap, including listed callouts.
This mode does not append unqualified reference rows. Iterative pattern-label
exports expand their page bounds to keep long literal event identifiers visible.

Missing heatmap values render black (user instruction, 2026-10-26) and remain
NaN in statistics, tables and caches. For median-centered PSI, black can therefore
mean missing or zero deviation; the missingness tables distinguish them.

Native event mirrors accept `--additional-fdr 0.5` alongside `--pval 0.005
--pval-type raw --strict-thresholds` for the joint filter, and `--caption-file`.
Protein mirrors accept `--rawp 0.005 --fdr 0.5` and `--caption_file`; original
native FDR values and correction families are reused. FDR 0.5 is not FDR 0.05.
`pseudobulk_depth_qc.expected_detected_features` is a descriptive exact
without-replacement detection expectation, not a differential test or imputation.

Iterative MarkerFinder consolidations preserve supplied splice-event callouts beyond
`top_n` using the existing exact/opposed-arm and shared-exon matching rules. Only callouts with a native discovery assignment are displayed; unassigned
reference rows are off by default (2026-10-26 user correction). The complete source
matrix, label roster and assignment table remain retained, including callout
annotations that do not qualify for this figure. `.displayed_rows.tsv` records
the displayed native assignments. Historical reference-block replay requires an
explicit API opt-in and is not used by the current figures.

Legacy iterative category labels now collapse only when their signed high/low
sample-design partitions are identical. Down-regulation of one group equals
up-regulation of its complement only for an identical sample partition. Partial
schemes retain unconstrained sample roles. `pattern_logic_aliases.tsv` records the
first-listed category and complete high/low/unconstrained sample identities;
`source_pattern` and all native winning statistics remain preserved. Use
`--no-collapse-pattern-logic` for historical label replay. Pattern identity uses no
expression correlation threshold.

For an explicit report hierarchy, `collapse_discrete_pattern_logic(...,
pattern_precedence=[(scheme, direction, level), ...])` chooses canonical
representatives in that order, appending only still-distinct design patterns.
For the ERP-5/LMPP-2 literal reference, the instructed order is up ERP-5
(shared lineage), the four up state/condition groups, up MDS-pre/post-SRSF2
(shared condition), followed by the corresponding down patterns only when
their exact signed design partition is absent above. This changes display
labels/order only; native scores, original winners, features and samples remain
unchanged. `--legacy-manifest` accepts `pattern_precedence` for fitting/replay.

Saved heatmap annotation can use `render_heatmap_from_cache(..., label_shared_introns=True, callout_color="#d62728", sample_labels={sample: sample for sample in cache.obs_names}, show_pattern_labels=True)`. Exon/intron matching remains gene-specific and uses parent region blocks from either event arm; it does not imply exact event or directional replication. Defaults preserve exon-only black callouts. Long sample labels receive additional footer space.
