# rna2lipid — what I changed on 2026-08-30, and what I proved

I replaced the lung lipid model with the sparse lipid-by-lipid ElasticNet
delivered in `/Users/saljh8/Downloads/rna2lipid_update/`, ported its trainer into
the module so any tissue can train the same architecture from a JSON config, and
pointed scALABLE at the new bundle.

Every number below comes from a script in
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/validation/`
that I ran on 2026-08-30 with `/usr/bin/python3` (scikit-learn 1.3.0, numpy
1.24.3, pandas 2.0.3, anndata 0.9.2).

## Scope

- **One tissue.** The shipped bundle covers human lung only.
- **Four scALABLE references.** `hs_lung_cellref_reference`,
  `hs_lung_hlca_reference`, `hs_lung_natri_reference`, `hs_lung_bpd_sun_reference`.
  The bone-marrow `lipid` modality is a different model, `rna2lipid/aml/`, and I
  did not touch it.
- **I did not retrain the lung model.** The two log2 lipid tables the delivered
  run used are absent from this Mac. I searched `/Users/saljh8/Dropbox`,
  `/Users/saljh8/Downloads` and `/Users/saljh8/Documents` and found only the
  non-log2 versions. So I could not reproduce the shipped bundle, and I say so
  again under "What I did not prove".

## Check A — the new api reproduces the prior api

Script: `validation/validate_api_lipidwise.py`. Report:
`validation/api_lipidwise_validation.json`.

I extracted the prior `api.py` from git commit `f57f8b3` and ran both versions on
the prior `MultiTaskElasticNetCV` bundle over the repository's two RNA tables.

| Rows | Lipids | Values | Max absolute difference | Bitwise identical |
| --- | --- | --- | --- | --- |
| 169 | 202 | 34,138 | 0.0 | yes |

The new api reads the prior bundle through the same `single_model` path and
returns the same numbers. The prior bundle carries an unfitted `scaler_y`, so
`target_scaling_mode` stays `none` for it, as before.

## Check B — the lipid-wise prediction math

Same script. I rebuilt the prediction independently: for each of the 202 lipids I
took its stored `coef_`, its stored `intercept_` and its stored gene list,
computed `design @ coef_ + intercept_` on the scaled RNA matrix, stacked the 202
columns in `Y_columns` order, and applied `scaler_y.inverse_transform`.

| Rows | Lipids | Values | Max absolute difference |
| --- | --- | --- | --- |
| 169 | 202 | 34,138 | 0.0 |

So `api._predict_aligned_matrix` selects the right genes per lipid, assembles the
columns in the right order, and inverts the target scaling correctly.

## Check C — the AnnData path equals the dataframe path

Same script. I built an AnnData from the same aligned matrix and predicted with
`chunk_size=7`, which forces 25 chunks.

| Values | Max absolute difference | Cause of the residual |
| --- | --- | --- |
| 34,138 | 1.64e-06 | the AnnData path casts the design matrix to float32 |

## Check D — what the shipped bundle actually covers

Same script.

| Property | Value |
| --- | --- |
| Training profiles | 50 |
| Training donors | 10 |
| Groups | END, EPI, MES, MIC, PMX, 10 profiles each |
| Bulk profiles in training | **0** |
| Lipids with an all-zero model | 1 of 202, `TG(51:1)` |
| Nonzero genes per lipid | minimum 0, median 21, maximum 57 |

`TG(51:1)` keeps no gene, so the model returns that lipid's training mean for
every sample. Treat that one column as a constant, not a prediction.

**Bulk RNA input is extrapolation.** I predicted the repository's own bulk and
cell-type RNA tables and counted how many predicted values fall outside the
training mean ± 3 SD of their lipid:

| Input rows | Fraction outside train mean ± 3 SD |
| --- | --- |
| Bulk profiles | 16.0% |
| Cell-type profiles | 0.15% |

scALABLE feeds cell-state profiles, which is the regime the model saw.

## Check E — the ported trainer equals the delivered trainer

Script: `validation/validate_trainer_equivalence.py`. Report:
`validation/trainer_equivalence.json`. Log:
`validation/trainer_equivalence_full.log`.

The delivered script reads cluster paths at import time, so I could not import
it. The script instead slices two functions out of
`provenance/train_sparse_lipidwise_delivered_2026-08-19.py` by line range
(`safe_pearson` at 648–689, `fit_sparse_lipidwise_elasticnet` at 838–1313) and
executes only those. I retyped nothing, so nothing can drift.

I fitted both implementations on the same 68 × 1,303 RNA matrix and 68 × 202
lipid matrix with the delivered hyperparameter grid.

| Compared | Count | Max absolute difference |
| --- | --- | --- |
| Lipid models | 202 | — |
| Coefficients | 30,625 | 0.0 |
| Intercepts | 202 | 0.0 |
| Predictions | 13,736 | 0.0 |

Selected gene lists, chosen top-N, chosen alpha and chosen l1 ratio match for all
202 lipids. `scaler_x`, `scaler_y`, `X_columns`, `Y_columns` and the whole
summary table are identical.

## Check F — scALABLE resolves to the new bundle

Script: `validation/validate_scalable_wiring.py`. Report:
`validation/scalable_wiring.json`.

scALABLE builds its application from `cellHarmony/webapp/app.create_app` and
computes through `cellHarmony/flask/pipeline.py`. The lipid modality enters at
`_build_imputed_lipid_adata`. The script drives that function with each lung
reference entry, loaded through the same `_lookup_reference` the job runner uses.

All four lung references resolve to
`/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/rna2lipid_hs_lung_lipidwise_bundle.pkl`,
report `architecture: lipidwise`, `target_scaling_mode: standard`, and return 202
lipids for 40 cells, with the `counts` layer and the UMAP carried over.

The input matrix is synthetic: 900 real model gene symbols plus 100 decoys. This
check proves which model runs and what shape it returns. It proves nothing about
accuracy.

## Check G — end-to-end training for a new tissue

Config: `validation/config_smoke_generic_tissue.json`. Outputs:
`validation/smoke_run/`.

I ran `cli train` through the tissue-agnostic paired-table path with a small
hyperparameter grid and the donor-disjoint holdout. It completed, wrote the
bundle and all six report files, and the manifest counted the losses:

| Manifest field | Value |
| --- | --- |
| Matched samples before the group filter | 68 |
| Matched samples after the group filter | 50 |
| Sample retained fraction | **73.5%** |
| Samples dropped by the group filter | 18, all group `UNKNOWN` |

**Retention below 90% is a finding.** Those 18 samples are the bulk profiles. A
bulk sample identifier such as `D001` carries no compartment suffix, so
`keyword_map` labels it `UNKNOWN`, and `keep_groups` drops it. The delivered
script does the same thing, which is why the shipped bundle holds 0 bulk
profiles. The old multi-task bundle used all 68.

## The delivered holdout number is not donor-held-out

The bundle carries `validation_global_metrics`: Pearson r 0.871, Spearman r
0.864, R² 0.754, RMSE 1.384, over 28 held-out profiles.

The split rule caps each donor at 3 profiles in the holdout set. Each donor
contributes exactly 5 profiles. So the cap cannot exclude a donor; it only limits
how many of that donor's profiles move. I measured it directly:
`donors_in_both_sets` is all 10 donors. The reported number therefore describes
prediction for a donor the model has already seen.

### ADDITION — how much the leakage is worth

Script: `validation/compare_holdout_leakage.py`. Report:
`validation/holdout_leakage_comparison.json`.

I ran the same trainer, on the same matrices, with the same delivered
hyperparameter grid, and changed only the holdout rule.

| Holdout rule | Train | Holdout | Donors in both | Pearson r | R² | RMSE |
| --- | --- | --- | --- | --- | --- | --- |
| Delivered constrained | 22 profiles | 28 profiles | 10 | 0.484 | 0.167 | 1.178 |
| ADDITION, donor-disjoint | 35 profiles | 15 profiles | 0 | 0.423 | −0.042 | 1.135 |

Removing the donor leakage costs 0.061 Pearson r and 0.209 R², and drives R²
below zero, which means the model beats the per-lipid training mean on none of
the held-out variance. The two rows also differ in training-set size, so the gap
mixes leakage with sample count. I claim the direction, not the exact
magnitude.

I add this comparison. The delivered constrained rule stays the default in
`configs/lung_lipidwise.json`, because that is the protocol the shipped bundle
used.

## What I did not prove

1. **I did not reproduce the shipped bundle.** Its two log2 lipid tables are not
   on this Mac. `configs/lung_lipidwise.json` names them and the README says
   where to get them.
2. **I did not reproduce Pearson r 0.871.** Both runs above, on the repository's
   lipid tables with the delivered grid, land near r 0.42–0.48. The gap tracks
   the targets, not the model: the delivered holdout reports `True_SD` 2.788,
   the repository tables give `True_SD` 1.29. The same absolute error against
   half the between-sample spread gives a much lower r and R². I have not
   verified the 0.871 figure and I cannot until the log2 tables are here.
3. **No real lung h5ad ran through scALABLE.** Check F used a synthetic matrix.
4. **Missing genes become zero before scaling.** `_align_dataframe` fills an
   absent model gene with 0.0 and then standardizes, so an absent gene enters the
   model as a strongly negative z-score, not as a neutral value. The prior api
   did the same, so this is not a regression, but it matters more now: a lipid
   whose model keeps 25 genes can rest on genes the query lacks. Every prediction
   summary reports `matched_genes` against `model_gene_count`; read it.
5. **`cli evaluate` has not run on the new bundle.** I made it architecture-aware
   and it imports, but the lung config it needs points at the absent log2 tables.

## Files I changed

| File | Change |
| --- | --- |
| `api.py` | replaced with the delivered dual-architecture version; `DEFAULT_BUNDLE_PATH` now names the lung lipid-wise bundle; added `LEGACY_MULTITASK_BUNDLE_PATH` |
| `__init__.py` | exports `LEGACY_MULTITASK_BUNDLE_PATH` |
| `pipeline.py` | added `build_paired_training_dataset`, `build_sample_metadata`, `build_dataset`; the legacy four-table builder is untouched |
| `training.py` | added `fit_sparse_lipidwise_elasticnet`, `train_lipidwise_bundle`, `build_holdout`, `holdout_global_metrics`, `resolve_alpha_grid`, `normalize_architecture`; `run_training` dispatches on architecture |
| `evaluation.py` | `_fit_and_predict_main` dispatches on architecture; the fold loop and summary name the architecture that ran |
| `configs/lung_lipidwise.json` | new, reproduces the delivered lung run |
| `configs/template_new_tissue.json` | new, the template for a new tissue |
| `configs/default_training.json` | unchanged |
| `cellHarmony/flask/pipeline.py` | `_build_imputed_lipid_adata` takes `reference_entry` and honours a per-reference `lipids` bundle path, expression scale and log base; the job log records the resolved bundle and architecture |
| `cellHarmony/flask/reference_config.json` | 20 inserted lines, an explicit `lipids` bundle block for the 4 lung references; nothing else changed |

## Files I added

| File | Purpose |
| --- | --- |
| `rna2lipid_hs_lung_lipidwise_bundle.pkl` | the new lung bundle, md5 `4c586793f2bc1d0d449ad6ea2ea35ccd`, identical to the delivered file |
| `provenance/train_sparse_lipidwise_delivered_2026-08-19.py` | the delivered training script, verbatim |
| `provenance/api_delivered_2026-08-19.py` | the delivered api, verbatim |
| `data/Bulk_lipids_cleaned_normalized_median_527.csv` | copied from `/Users/saljh8/Downloads/rna_lipid-main/data/`; `configs/default_training.json` already named it and it was absent |
| `data/cell_lipids_cleaned_norm_median_286.csv` | same |
| `validation/*.py`, `validation/*.json`, `validation/*.log` | the checks above |

I deleted nothing and I overwrote no input. The prior bundle
`newnormelastic_multitask_try.pkl` is still on disk, unchanged.

## Deployment note — git does not carry the bundle

`.gitignore` line 75 ignores `altanalyze3/components/rna2lipid/*.pkl`. That rule
predates this work and already excluded `newnormelastic_multitask_try.pkl`, so
neither lung bundle has ever entered version control. A commit of these changes
therefore ships the code and the config, not the model bytes.

Copy the file by hand to every machine that runs scALABLE:

```
/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/rna2lipid/rna2lipid_hs_lung_lipidwise_bundle.pkl
```

md5 `4c586793f2bc1d0d449ad6ea2ea35ccd`, 1,726,546 bytes.

Two alternatives, if you want git to carry it. Add a negation line to
`.gitignore`, or gzip the bundle to `.pkl.gz`, which the sibling
`rna2lipid/aml/artifacts/rna2lipid_aml_bundle.pkl.gz` already does and which the
`*.pkl` rule does not match. Both change repository policy, so I did neither.
