# rna2lipid reference normalization workflow

This workflow aligns the available sorted-cell and bulk lung lipid tables, builds a combined RNA–lipid reference, and retrains an abundance predictor. A normalization method is accepted only if it preserves every observed within-dataset lipid differential. It does not require the missing full-cohort raw MS data.

The completed October 1, 2026 run selected lipid-specific PMX alignment. It combines 45 sorted profiles and 19 bulk profiles from 25 donors, with 1,303 shared RNA genes and 219 complete lipid features. All 254,259 sample-pair differentials and all 2,190 population tests were preserved within numerical precision. The 1,031 significant population tests remained significant with no changed calls.

Normalization passed that requirement. Predictive accuracy is a separate result: 43 lipids meet the workflow's provisional internal support criteria. The prediction command exports those lipids by default. Disease fold accuracy has not been established for the full panel.

The workflow now includes MS1 re-extraction from all 30 lipidomics mzML runs in [MassIVE MSV000081973](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000081973). Raw replication QC supports 47 features in the combined reference; 11 also pass internal RNA prediction support. A second candidate bundle and positive abundance exports provide those 11 features on a reference-normalized MS1 intensity scale. Full original manual identification has not been reproduced for all 286 published annotations.

The later fold-magnitude audit changes the interpretation of those 11 RNA outputs: their support screen checked direction and error, not preservation of effect magnitude. RNA predictions remain attenuated. Nested ridge and PLS trials did not resolve that limitation, and no feature passed the specified magnitude criteria for both population and within-population donor contrasts. The measured-reference remapping below preserves fold magnitudes exactly; the RNA-only estimates remain exploratory for disease effect sizes.

## Source evidence

Inputs are under `/Users/saljh8/Downloads/Lipidomics`:

- `LMEX0000003692_data_table.xlsx`: original processed lipid measurements. The left block contains 45 sorted profiles. The right block is a standardized copy; intervening average and fold columns are excluded.
- `newrna_cell_clair_filtered_symbol.csv`: gene-by-profile RNA. The 19 duplicate gene-symbol rows are averaged on their supplied numeric scale, matching the existing predictor's aggregation. Additional `.1` RNA profiles are not substituted for the original matched profiles.
- `feature_blankreduiction.csv`: profile-by-gene bulk RNA.
- `Preprocessing Process/unprocessed_bulk.csv`: feature-by-profile bulk lipid data, documented as already log2 transformed. Its 28 profiles supply calibration references; 19 have paired RNA and enter regression training.
- The preprocessing PDF and final CSVs supply an audit trail, not the new training targets. The stated sorted-cell linear-scale assumption conflicts with the workbook. The final sorted table reverses some contrasts and cannot be reconstructed from those stated steps.
- `validation/MSV000081973_reference/published_lipid_profiles.csv`: extracted Table S4 from the 2018 study, containing 286 features across three donors and five populations. Ten profiles from two donors have supplied RNA; D001 lacks sorted RNA. This is processed MS evidence, not newly extracted MS1 intensities.

The input hashes, sample exclusions, annotation policies and duplicate gene records are saved with the run. Only exact lipid structure and ion-mode matches are joined: `_POS` becomes `_P`, and `_NEG` becomes `_N`. Different modes, isomers and ambiguous lipid annotations are not silently merged. Of 232 matched features, 13 have missing measurements and are excluded from combined training. The full source extraction retains them separately.

All 45 original profiles have median zero across their original lipid panel. The published three-donor table additionally has per-lipid mean zero. Adding constant lipid offsets does not change its within-lipid contrasts, so external evaluation compares each population with PMX from the same donor.

## Normalization candidates

For lipid `l`, the selected alignment estimates

`offset[l] = median(bulk_log2[l]) - median(sorted_PMX_log2[l])`

and constructs

`combined_sorted_log2[s,l] = original_sorted_log2[s,l] + offset[l]`.

Bulk log2 measurements remain on their supplied scale. Linear reference abundance is `2 ** combined_log2`, with no clipping, extra logging or pseudocount. This is a multiplicative calibration in linear space and assigns each lipid its own positive abundance baseline. It does not force a sorted population to zero abundance or erase population contrasts.

| Candidate | Differential gate | Largest log2FC change after CSV round trip |
| --- | --- | ---: |
| Global PMX shift | Passed | approximately 1e-14 |
| Lipid-specific PMX shift | Passed and selected | approximately 1e-14 |
| Lipid-specific shift using shared donor identifiers | Passed | approximately 1e-14 |
| Additional sample median normalization | Rejected | 0.700 |
| Quantile normalization | Rejected | 3.404 |
| Separate matching of population lipid medians | Rejected | 5.027 |
| No alignment | Control only; incompatible starting scales | 0 |

Offsets are common to all sorted populations for a given lipid. Separate population offsets would change population-to-population differences and are rejected. Target z-scores are used only inside regression fitting and are inverted before prediction export.

## Mandatory preservation checks

The gate compares the original and combined log2 differences for every lipid in every sample pair within each dataset: 990 sorted pairs and 171 bulk pairs, each across 219 lipids. Its absolute tolerance is 1e-9 log2 units.

It also recomputes all ten sorted-population contrasts, pairing profiles by donor. It checks the mean log2 difference, paired t-test p-value, Benjamini–Hochberg FDR, and significance call. The call threshold is FDR < 0.05 and absolute fold >= 1.5. All accepted alignments preserve these values. No cross-dataset fold is claimed to have a measured ground truth simply because an alignment has been fitted.

The saved linear reference is checked for strict positivity and finite values. Taking its log2 must reconstruct the combined log2 reference. Reference integrity is checked again after serialization, separately from predictive-model accuracy.

## Method selection and regression validation

The predictive model is linear kernel ridge regression, with RNA and target scalers fitted only on training profiles. Its five candidate penalties are 1, 10, 100, 1000 and 10000. There is no supervised gene selection using held-out measurements.

An outer leave-one-donor-out loop covers all 25 donors. Inner donor-disjoint validation selects the penalty for each accepted alignment. A held-out donor is excluded from sorted data, bulk data and calibration references together. Method selection minimizes the inner prediction loss. Methods within 1% of the best loss are compared using training-only leave-one-shared-ID-out calibration consistency. This rule is applied inside every outer fold, not selected after inspecting outer test performance.

Four identifiers, D019, D024, D036 and D041, have both sorted PMX and bulk lipid values. Holding each out gives alignment RMSEs of 1.819 log2 for the global shift, 0.958 for lipid-specific PMX alignment, and 1.098 for shared-ID alignment. These checks assess candidate bridge consistency. Matching identifiers do not independently establish the same specimen input amount or physical assay scale.

A population-mean baseline tests whether RNA predicts more than population identity. A separate nested sorted-only model tests whether combination improves sorted prediction. Combined sorted RMSE is 0.825 log2, compared with 0.830 for sorted-only training: the improvement is small.

The selected combined workflow gives 86.9% direction agreement for population contrasts within donors, and 58.5% for donor contrasts within a population. Both figures restrict direction scoring to measured absolute log2FC >= 0.5. Thus population-level performance must not be presented as proof of disease prediction within a cell type.

External evaluation uses the MSV000081973 published contrasts. D001, D008 and D011 are excluded from both training modalities and calibration; retaining their bulk profiles would leak donor information. For the two RNA-matched donors, 109 exact-mode features provide 872 population-to-PMX contrasts. Direction agreement is 82.4%, correlation is 0.615, and RMSE is 0.973 log2FC. The three-donor table has not been used to restore unmeasured sample normalization factors or to claim physical concentrations.

## Prediction support and units

The provisional internal criteria require positive prediction skill against the population-mean baseline, sorted RMSE <= 1 log2 unit, and within-population donor-contrast direction agreement >= 0.7 for contrasts with measured absolute log2FC >= 0.5. Forty-three lipids meet all three criteria. This is a descriptive support screen on cross-validation results, not a separate confirmatory disease test.

Outputs are bulk-referenced normalized relative abundance units. They are positive and support correct within-lipid fold arithmetic. Peak intensity is not automatically a concentration, and different lipid species' signals should not be added as physical quantities without analytical calibration.

New RNA must use the same numeric representation as the supplied training RNA. The command does not infer a log base or transform arbitrary counts. It requires all 1,303 model genes and rejects missing values instead of silently filling absent genes with zeros. Duplicate gene symbols are aggregated as during training.

## Commands

Run these commands from the inner `altanalyze3` repository directory. `/usr/bin/python3` on this machine supplies numpy, pandas, scipy, scikit-learn, statsmodels, threadpoolctl and anndata.

Reconstruct, compare, validate and train:

```bash
/usr/bin/python3 components/rna2lipid/reference_normalization.py run \
  --source-dir /Users/saljh8/Downloads/Lipidomics \
  --output-dir components/rna2lipid/artifacts/reference_normalization_20261001 \
  --external-reference components/rna2lipid/validation/MSV000081973_reference/published_lipid_profiles.csv
```

Recheck every saved reference without refitting models:

```bash
/usr/bin/python3 components/rna2lipid/reference_normalization.py verify \
  --output-dir components/rna2lipid/artifacts/reference_normalization_20261001
```

Predict the supplied sorted RNA in its original representation:

```bash
/usr/bin/python3 components/rna2lipid/reference_normalization.py predict \
  --bundle components/rna2lipid/artifacts/reference_normalization_20261001/reference_calibrated_bundle.pkl \
  --input /Users/saljh8/Downloads/Lipidomics/newrna_cell_clair_filtered_symbol.csv \
  --transpose \
  --output-log2 /private/tmp/lipids_supported_log2.csv \
  --output-linear /private/tmp/lipids_supported_linear.csv
```

The default exports 43 internally supported lipid features. Add `--include-unsupported` to explicitly export all 219 exploratory outputs. For sample-by-gene input, omit `--transpose`. The full candidate bundle is compatible with the existing `rna2lipid.api.load_bundle` interface, but the API itself does not apply the support filter; use this workflow's prediction command or filter by its metadata.

Run regression checks:

```bash
/usr/bin/python3 components/rna2lipid/validation/validate_reference_normalization.py
/usr/bin/python3 components/rna2lipid/validation/validate_bulk_calibration.py
/usr/bin/python3 components/rna2lipid/validation/validate_log2_folds.py
```

For disease comparisons, use donor-level profiles within the same cell state and compare the logged outputs by subtraction. For arithmetic-mean abundance ratios, exponentiate first and give each donor equal weight. `log2_fold_changes.py` implements both contrasts. Do not take ratios of logged numbers, add arbitrary linear baselines, or treat individual cells as independent biological donors.

## Result files

All files are under `artifacts/reference_normalization_20261001`:

- `combined_training_reference_log2.csv` and `combined_training_reference_linear.csv`: the selected combined reference, containing all 219 complete features.
- `paired_RNA.csv`, `sample_metadata.csv`, `manifest.json`, `excluded_lipids.csv`: matched inputs and exclusions.
- `differential_preservation.json`, `population_differential_checks.csv`: the mandatory fold, p-value, FDR and significance checks.
- `calibration_bridge_checks.json`, `calibration_bridge_heldout_residuals.csv`, `selected_calibration_offsets.csv`: calibration evidence and parameters.
- `donor_validation.json`, `outer_fold_choices.json`, `OOF_*.csv`, `per_lipid_validation.csv`: held-out predictions, baselines and support assessments.
- `MSV000081973_validation.json`, `MSV000081973_heldout_contrasts.csv`: independent donor-excluded sorted-profile evaluation.
- `reference_calibrated_bundle.pkl`: the newly trained candidate, including fitted scalers, calibration metadata and support list.
- `internally_supported_predictions_log2.csv` and `internally_supported_predictions_linear.csv`: 150 supplied RNA profiles with 43 supported outputs.
- `supplied_sorted_RNA_predictions_*.csv`: exploratory full-panel predictions. Training-set rows are resubstitution; additional supplied RNA profiles are extrapolation.
- `legacy_preprocessing_audit`: reconstruction of the old preprocessing evidence.
- `workflow_result.json`: the overall result and scope.

The production default bundle and scALABLE reference registry have not been changed. Activation of this candidate requires preserving the support policy and checking the incoming RNA representation.

## MassIVE MS1 reanalysis and donor anchoring

All 30 original lipidomics mzML runs were downloaded: three donors, five populations, and two ionization modes. The selected files are the original `peak/peak` uploads, rather than the duplicate `ccms_peak` conversion. The 15 proteomics runs were excluded. Files contain centroid MS1 and MS2 scans. The download manifest records expected byte counts and SHA-256 hashes. A browser session cookie resolves the initial HTTP 429 download failures; cookies are held in memory by the reproducible downloader.

The mass and fragment rules come from [PNNL LIQUID](https://github.com/PNNL-Comp-Mass-Spec/LIQUID), the identification software used in the original study. Their pinned commit, file hashes, authorship and Apache 2.0 license are stored under `validation/MSV000081973_reference/LIQUID_rules`. Composition reconstruction covers all 286 published targets, including ether lipids, plasmalogens, sum compositions and composite annotations. Composite annotations retain their ambiguity.

`pyNeoQuant/quant/lipid_ms1.py` streams centroid spectra, keeps MS1 separate from MS2, respects ion mode and RT units, and rejects unsupported profile data. Extraction uses a 5 ppm mass window, a 45-second half-window around locked RT, and at least three positive MS1 scans. It reports apex intensity and XIC area. Missing peaks remain missing; no pseudocount, zero filling or abundance imputation is applied.

RT locking uses pooled MS2 fragment evidence and a flagged MS1-only fallback. Fragment matching uses 0.3 Da tolerance, accommodating the CID acquisition, and a 1% relative intensity floor. Fragment support is evidence for the published target candidate, not a new identification FDR. Accurate-mass matches alone do not establish molecular identity. Targets with competing published annotations within 10 ppm are excluded from individual-species evaluation. The original manual LIQUID/MZmine peak and isomer assignments have not all been recovered.

The re-extracted MS1 apex values are log2 transformed, centered by sample median and then by lipid mean to compare with published Table S4. RT selection never reads the published abundance values. For 150 unambiguous target annotations with MS2 support and evaluable contrasts, 1,767 population-to-PMX comparisons give correlation 0.861, RMSE 0.596 log2FC and 96.35% direction agreement where published absolute log2FC is at least 0.5. These pooled results do not excuse individual failures: `per_lipid_replication.csv` records all features, including poor matches and missing profiles.

The proposed one-donor abundance anchor is tested directly. Each donor supplies its five measured population baselines in turn; the other two donors are held out. Twelve combinations of baseline-offset method and measurement normalization are evaluated. A per-lipid offset common across populations, fitted from the anchor donor's five profiles, gives these results:

| Measured scale | Held-out RMSE, log2 | Median absolute error, log2 | Published population differentials retained |
| --- | ---: | ---: | --- |
| Unnormalized peak intensity | 1.042 | 0.716 | Yes |
| Sample-median-normalized intensity | 0.618 | 0.295 | Yes |
| Separate ion-mode median normalization | 0.580 | 0.241 | Yes |

Using only the anchor donor's PMX sample for each lipid gives sample-normalized RMSE 0.641 log2. Separate population-specific lipid offsets slightly reduce abundance error but change population-to-PMX folds by as much as 5.517 log2 units on the sample-normalized scale. They fail the mandatory differential gate. Better abundance fit alone is insufficient to accept a calibration.

This validates the anchor equation on the available subset: `A_new = A_anchor * 2 ** (F_new - F_anchor)`. The anchor and fold measurements must share a normalization convention. Applying normalized folds to unnormalized instrument intensities gives larger verification errors, as the held-out results demonstrate.

The RNA model is also compared directly against re-extracted MS1 contrasts, with the MassIVE donors excluded from both original training modalities and calibration. There are 446 evaluable contrasts across 56 exact-mode features: correlation 0.607, RMSE 0.874 log2FC and 82.65% direction agreement for re-extracted absolute log2FC >= 0.5. This remains population contrast validation, rather than independent disease validation within a population.

## MS1-scaled combined training reference

Raw replication support requires all 15 profiles observed, no competing published mass annotation, MS2 support, per-lipid replication RMSE <= 0.75 log2 and Spearman correlation >= 0.7. These descriptive criteria support 121 published targets, of which 47 match the combined RNA–lipid reference. Eleven also pass its provisional internal RNA support screen. Selection by replication QC is recorded explicitly; this is not a new independent confirmatory test.

The raw signal reference equalizes sample medians to the measured D001 PMX median and uses the three-donor median PMX value for each lipid. For each of the 47 matched lipids, the difference between that PMX baseline and the existing combined-reference PMX baseline is added to **every sorted and bulk profile**. This changes the signal units by one multiplicative factor per lipid. It preserves all within-dataset differences and assigns positive baselines without independently resetting population means.

The serialized 64-profile reference passes 54,567 differential checks with maximum change 7.11e-15 log2. Its 470 population tests are checked for unchanged fold, p-value, FDR and significance calls, using the same 47-feature panel before and after calibration. Refit kernel-ridge predictions equal the prior predictions plus the calibration offset to 2.85e-14 log2 units. These checks establish calibration integrity; they do not establish equal extraction amounts or physical molar concentrations between cohorts. Median PMX composition transfer remains the biological reference assumption.

The new bundle exports the 11 features passing both support screens by default. `--include-unsupported` exports all 47 raw-supported exploratory features. The broader 219-feature bulk-scaled reference and its 43 internally supported outputs remain available separately.

### Fold magnitude audit

Direction agreement must be separated from magnitude agreement. The saved `MSV000081973_reanalysis_20261001/fold_magnitude_audit.json` reports both. For the 47 raw-supported reference features, 564 raw-versus-published population contrasts have correlation 0.928, RMSE 0.434 log2FC, median absolute error 0.268 log2FC, and OLS slope 0.848. Approximately 79.1% are within 0.5 log2 units. Raw re-extraction therefore reproduces magnitude reasonably, with residual attenuation and feature-specific errors.

RNA imputation is less accurate. For the 11 exported features, 88 donor-held-out population contrasts against re-extracted MS1 have correlation 0.753, RMSE 0.662 log2FC, median absolute error 0.307 log2FC, and OLS slope 0.715. Only 64.8% are within 0.5 log2 units. Across all 56 evaluable features, the RNA slope is 0.453. These estimates retain meaningful variation but still attenuate measured fold magnitudes. Positive intensities and exact calibration preservation do not establish accurate disease effect sizes. No post-hoc multiplication by an inverse slope is justified by this descriptive audit.

### Supplied final bulk scale

The supplied `Final_bulk_Bulk_lipids_cleaned_normalized_median_527_log210.csv` is consistent with `z = log2(1 + 10*A)`, where `A = 2**x` and `x` is the original supplied bulk log2 value. The `10` in this transformation is a multiplier, not a logarithm base. The protocol audit reproduces 11,412 of 11,421 matched values within 0.101 log2 units; nine discrepancies remain. The filename alone must not determine the scale.

To invert that transformation, use `A = (2**z - 1)/10`, or `x = log2((2**z - 1)/10)`. Matching logarithm bases does not by itself align measurement baselines. The MS1-scaled linear prediction CSV cannot be directly concatenated with this logged bulk table. The completed combined reference instead uses the original mode-preserved bulk log2 table, aligns sorted values to its per-lipid PMX reference, and then translates both datasets by the same per-lipid offset to the normalized MS1 signal scale. The resulting `combined_reference_MS1_log2.csv` and `combined_reference_MS1_linear.csv` already contain both modalities on that chosen reference scale.

### Direct remapping of the supplied final bulk dataset

`remap_bulk_abundance.py` additionally remaps the actual supplied final bulk table, rather than substituting the original bulk table for it. All 527 bulk features are restored to their native linear and log2 scales by inverting the documented transformation. Forty-seven features have independently verified ion-mode source matches and a shared sorted/MS1 calibration. Those 47 are translated by the same per-lipid factors used for the sorted reference. The other 480 retain native bulk units and are not represented as having sorted/MS1 calibration.

The output directory `artifacts/supplied_bulk_MS1_remapping_20261001` contains:

- `supplied_bulk_all_features_native_linear.csv` and `supplied_bulk_all_features_native_log2.csv`: all 27 profiles and 527 features, with the documented transformation inverted; these remain in native bulk units.
- `supplied_bulk_MS1_linear.csv` and `supplied_bulk_MS1_log2.csv`: all 27 supplied bulk profiles for the 47 calibrated features, in the same units as the sorted MS1 reference and prediction exports.
- `combined_sorted_and_supplied_bulk_MS1_linear.csv` and `combined_sorted_and_supplied_bulk_MS1_log2.csv`: 72 profiles, combining 45 sorted profiles with the actual 27 supplied bulk profiles. Bulk observations have `_BULK` identifiers, distinct from their sorted PMX observations.
- `mapping_and_exclusions.csv`, `bulk_value_translation_audit.csv`, `sample_metadata.csv` and `remapping_result.json`: source annotation checks, individual value transformations and preservation evidence.

All 63,027 sample-pair differences in the recovered log2-abundance convention are preserved through remapping and CSV serialization, with maximum error 7.11e-15 log2. The input's `+1` pseudocount is explicitly removed before these checks. For example, the supplied D075 `PC(O-18:0/18:1)` value 25.751928 maps to 170,747.48 normalized MS1 intensity units; the corresponding RNA-derived estimate is 219,081.74 in that same reference space. The difference between those two values reflects prediction/assay variation, not a remaining logarithm-base mismatch.

```bash
/usr/bin/python3 components/rna2lipid/remap_bulk_abundance.py \
  --source '/Users/saljh8/Downloads/Lipidomics/Preprocessing Process' \
  --reference components/rna2lipid/artifacts/reference_normalization_20261001 \
  --ms1-reference components/rna2lipid/artifacts/MS1_reference_normalization_20261001 \
  --out components/rna2lipid/artifacts/supplied_bulk_MS1_remapping_20261001
```

### Attempts to improve RNA-predicted fold magnitudes

`train_fold_magnitude_model.py` compares donor-nested ridge and PLS regression. Ridge penalties span 0.01–1000; PLS latent component counts span 1–20. Parameters are chosen per lipid using inner donor-disjoint predictions and an objective combining fold error with attenuation penalty. The outer donor is excluded from both training modalities and reference alignment. Reference calibration remains a shared per-lipid offset; the measured reference is never stretched to make predictions look accurate. No inverse-slope inflation is applied to predictions.

For the previous 11-output panel, the external population-to-PMX OLS slope increases only from 0.715 to 0.723 with ridge or 0.732 with PLS. RMSE is 0.662, 0.667 and 0.651 log2FC respectively. These changes do not establish magnitude preservation. Across all 47 features, outer within-population donor-contrast slopes remain approximately 0.19–0.20, with correlations approximately 0.35.

The magnitude support criteria require correlation >= 0.7, OLS slope 0.8–1.2, RMSE <= 1 log2FC and median absolute error <= 0.5 log2FC. No feature passes both contrast tiers. Ridge has one provisional population-only feature and PLS two; neither has a feature passing the within-population donor tier. These are descriptive internal screens, not independent disease validation. Candidate bundles and exploratory predictions are retained for inspection, but neither replaces the existing default.

Results are saved in `artifacts/fold_magnitude_model_20261001` and `artifacts/PLS_fold_magnitude_model_20261001`, including nested parameters, per-lipid errors, external MS1 comparisons and candidate bundles. Run the script with `--engine ridge` or `--engine pls`. Regression tests check inversion, multiplicative fold preservation and agreement between PLS component-prefix predictions and independently fitted models.

From the AltAnalyze3 package directory, reproduce the added stages with:

```bash
RAW=/Users/saljh8/Downloads/Lipidomics/MSV000081973
PYNEO=/Users/saljh8/Documents/GitHub/pyNeoQuant/.venv/bin/python
EVIDENCE=components/rna2lipid/validation/MSV000081973_reference
REANALYSIS=components/rna2lipid/artifacts/MSV000081973_reanalysis_20261001
REFERENCE=components/rna2lipid/artifacts/reference_normalization_20261001
MS1_REFERENCE=components/rna2lipid/artifacts/MS1_reference_normalization_20261001

"$PYNEO" components/rna2lipid/download_massive_lipids.py --out "$RAW"
"$PYNEO" components/rna2lipid/massive_lipid_reanalysis.py \
  --mzml-dir "$RAW/mzML" --published "$EVIDENCE/published_lipid_profiles.csv" \
  --rules-dir "$EVIDENCE/LIQUID_rules" --out "$REANALYSIS"
"$PYNEO" components/rna2lipid/validate_massive_anchors.py \
  --raw-dir "$REANALYSIS" --published "$EVIDENCE/published_lipid_profiles.csv" \
  --reference-dir "$REFERENCE"
/usr/bin/python3 components/rna2lipid/build_raw_scale_reference.py \
  --reference-dir "$REFERENCE" --raw-dir "$REANALYSIS" --out "$MS1_REFERENCE"
/usr/bin/python3 components/rna2lipid/reference_normalization.py predict \
  --bundle "$MS1_REFERENCE/MS1_reference_bundle.pkl" \
  --input /Users/saljh8/Downloads/Lipidomics/newrna_cell_clair_filtered_symbol.csv \
  --transpose --output-log2 /private/tmp/lipids_MS1_log2.csv \
  --output-linear /private/tmp/lipids_MS1_linear.csv
```

The pyNeoQuant environment supplies numpy, pandas and scipy for raw reanalysis; `/usr/bin/python3` supplies the existing sklearn/statsmodels reference-training environment. The core pyNeoQuant MS1 extractor itself uses only the Python standard library and is also available through `pyneoquant lipid-ms1` for externally annotated target TSVs.

The completed raw stages produce `replication_result.json`, `anchor_validation.json`, spectrum/fragment/peak tables, per-lipid QC, and `RNA_prediction_vs_reextracted_MS1.csv`. The MS1 reference directory contains the combined log2 and positive linear training references, fitted candidate bundle, offsets, support list, statistical preservation tables and 150-profile prediction exports. Direct prediction through the serialized bundle was checked against the stored 11-feature exports, with strictly positive values and log2 recovery at floating-point precision.

## Corrected all-reference IPF lung-tissue validation

The active IPF lipid source is `/Users/saljh8/Downloads/Lipidomics/IPF/10_results_with_statistics.csv`, supplied as lung tissue. The incorrectly supplied plasma analysis, source evidence and plasma-only code were moved to `delete/IPF_plasma_20261001`; its conclusions are superseded.

`validate_ipf_candidate.py` refits the same isolated candidates on all reference profiles, holding out no samples and fitting no IPF lipid outcomes. The primary reference has 63 matched profiles (45 sorted +18 bulk), 24 donors and 1295 available RNA genes. Earlier-reference candidates use 64 profiles. Official model selection is unchanged.

The corrected lipid table contains 544 ion-mode features and 40 tissue profiles. Its identifiers imply 10 IPF donors with three tissues each and 10 control donors with one tissue each. Primary tests average repeated IPF tissues within donors and compare 10 versus 10 donor units. Supplied pooled-tissue raw/adjusted statistics remain distinct from the recomputed donor-level Welch/BH tests. Source log base, paper citation and stage sample assignments remain pending protocol confirmation; exponentiated fold magnitudes are explicitly conditional on log2. Direction comparisons do not require that log-base assumption.

The independent RNA source remains whole-lung GSE213001: 20 IPF and 14 control donors, with 101 known-region profiles averaged within donor/region. Both tissues are now lung, but the cohorts are unmatched. Exact lipid structure and ion mode are required. Positive/negative signed folds, full-panel raw/BH values, all specific lipids and sensitivity analyses appear in [the corrected report](artifacts/IPF_candidate_validation_20261001/IPF_CANDIDATE_REPORT.md).

All-match, global-top25, top25 among matched features and donor-p ranked comparisons are retained. The supplied table has 110 raw p-values and 75 adjusted p-values rounded to zero, so its top25 is not uniquely ranked; source p ties and selected identities are disclosed rather than choosing a concordant subset.

```bash
/usr/bin/python3 components/rna2lipid/validate_ipf_candidate.py
/usr/bin/python3 components/rna2lipid/validation/validate_ipf_candidate.py
```

## Study sources

The lipid MS methods and published normalization are described in [Kyle et al. 2018](https://www.nature.com/articles/s41598-018-31640-x). The additional measurements come from its [supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-018-31640-x/MediaObjects/41598_2018_31640_MOESM2_ESM.xlsx), associated with [MSV000081973](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=1c494cc549464750b8b13358241dae21). The workflow uses the existing local processed measurements; access to the full raw-MS inventory remains unverified.
