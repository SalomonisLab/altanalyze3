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
