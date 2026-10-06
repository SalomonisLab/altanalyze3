"""Report same-ElasticNet candidate, matched-panel comparison and fold estimates."""

try:
    from .candidate_integrity import reject_retired_workflow
except ImportError:
    from candidate_integrity import reject_retired_workflow

if __name__ == "__main__":
    reject_retired_workflow('report_lungmap_elasticnet_ms1.py')

import itertools
import json
import pickle
import numpy as np
import pandas as pd
from evaluate_lungmap_elasticnet_ms1 import OUT,LABEL,HERE,PRODUCTION


def signed(delta):
    return np.where(delta>0,np.exp2(delta),np.where(delta<0,-np.exp2(-delta),0))


def report(out=OUT):
    reject_retired_workflow('report_lungmap_elasticnet_ms1.py:report')
    training=json.loads((out/'training_method_audit.json').read_text())
    inference=json.loads((out/'inference_method_audit.json').read_text())
    method=json.loads((out/'differential_method_audit.json').read_text())
    equivalence=json.loads((out/'full_delivered_trainer_equivalence.json').read_text())
    settings=json.loads((out/'production_training_settings.json').read_text())
    summary=pd.read_csv(out/'cohort_concordant_union_summary.csv')
    support=pd.read_csv(out/'cohort_concordant_lipid_celltypes.csv')
    tests=pd.read_csv(out/'all_measured_vs_imputed_tests.csv')
    best=pd.read_csv(out/'best_concordant_celltypes.csv')
    primary=summary[summary.overlap_gate.eq('both_BH010') & summary.minimum_imputed_fold.eq(1)]
    audit=pd.read_csv(HERE/'artifacts/MS1_reference_normalization_20261001/population_differential_checks.csv')
    # Check retained target differences and descriptive training prediction folds.
    directory=out/LABEL
    y=pd.read_csv(directory/'candidate_training_lipids_log2.csv',index_col=0)
    pred=pd.read_csv(directory/'training_predictions_log2.csv',index_col=0)
    populations=y.index.str.split('_').str[-1]
    rows=[]
    for a,b in itertools.combinations(sorted(set(populations)),2):
        raw=np.log2(np.exp2(y[populations==a]).mean())-np.log2(np.exp2(y[populations==b]).mean())
        inferred=np.log2(np.exp2(pred[populations==a]).mean())-np.log2(np.exp2(pred[populations==b]).mean())
        for lipid in y.columns:
            rows.append(dict(case=a,control=b,lipid=lipid,target_signed_fold=float(signed(raw[lipid])),
                             training_prediction_signed_fold=float(signed(inferred[lipid])),
                             target_log2_effect=float(raw[lipid]),prediction_log2_effect=float(inferred[lipid]),
                             effect_log2_error=float(inferred[lipid]-raw[lipid]),direction_agrees=bool(np.sign(raw[lipid])==np.sign(inferred[lipid]))))
    fold=pd.DataFrame(rows);fold.to_csv(out/'training_reference_fold_magnitude_audit.csv',index=False)
    replay=pd.DataFrame(inference['inputs_and_replay']);replay=replay[replay.production_replay_rows.notna()]
    replay.to_csv(out/'production_inference_replay_by_study.csv',index=False)
    exact=primary[primary.model.isin([LABEL,'current_matched47_explicit_log2'])]
    deployed=primary[primary.model.isin([LABEL,'current_as_deployed'])]
    intersections=[]
    selected=support[support.overlap_gate.eq('both_BH010') & support.minimum_imputed_fold.eq(1)]
    for (model,unit),g in selected.groupby(['model','unit']):
        a=g[g.contrast.eq('Adams2020__IPF_vs_Healthy')].set_index('measured_feature')
        n=g[g.contrast.eq('Natri2024__IPF_vs_Healthy')].set_index('measured_feature')
        both=set(a.index)&set(n.index)
        same=sum(bool(set(a.loc[f,'concordant_cell_types'].split('; ')) & set(n.loc[f,'concordant_cell_types'].split('; '))) for f in both)
        intersections.append({'model':model,'unit':unit,'Adams_concordant':len(a),'Natri_concordant':len(n),
                              'intersection':len(both),'same_celltype_support_in_both':same})
    intersection=pd.DataFrame(intersections)
    intersection.to_csv(out/'Adams_Natri_concordant_intersections.csv',index=False)
    lines=['# MASSIVE-corrected LungMAP candidate using the production ElasticNetCV method','',
           'This supersedes the earlier kernel-ridge candidate as the requested same-algorithm evaluation. This candidate uses the deployed lipid-wise ElasticNetCV trainer and settings; production bundles, predictions, databases and differentials were not replaced.','',
           '## What is identical and what differs','',
           'Training uses one ElasticNetCV per lipid, training-set StandardScaler for X and Y, per-lipid absolute Pearson gene ranking, top-N candidate fits, and training R² minus the saved sparsity penalty to select among gene sets. Every saved hyperparameter is inherited directly from the production bundle. Cyclic coordinate selection, default tolerance and intercept settings are unchanged. The environment uses scikit-learn 1.6.1, matching the production estimator pickle metadata.','',
           f'All {training["input_genes"]:,} production training genes and their order are retained. Original RNA preprocessing reconstructs every production X-scaler mean and scale within 1e-12. Corrected target values are read directly from combined_reference_MS1_log2.csv, restricted to the original sorted-cell groups and mapped to the existing lipid names. No bulk profiles, extra RNA study alignment, shared-gene refit or kernel ridge model is introduced. The production missing-input policy (zero fill) is retained.','',
           f'The algorithm is identical, but corrected reference coverage differs: {training["candidate_training_samples"]} of the original {training["production_training_samples"]} training profiles have corrected measurements, and this candidate is deliberately restricted to 47 of the 202 production lipid outputs by the prior raw-MS QC screen. This does not establish that the other 155 outputs cannot be reconstructed. The broader 219-feature reconstructed reference contains names matching 163 production outputs (47 selected here and 116 outside this panel); 39 production names remain unmatched. D071 is one of the ten original training donors, each contributing five sorted-cell profiles, not a patient selected from the atlas evaluation cohorts. Its five profiles are absent from this reconstructed target table; the reason has not been established. The production sample order is preserved among available profiles. This is a same-algorithm corrected-reference comparison, not an isolated target-scaling experiment on an identical 50×202 training matrix. No external holdout metric is claimed; final all-reference fitting and internal three-fold ElasticNetCV are used.','',
           'Settings:','',pd.DataFrame([{'parameter':k,'value':str(v)} for k,v in settings.items() if k!='alpha_grid']).to_markdown(index=False),'',
           f'Alpha grid: the exact 30 production values spanning 0.001–10. All {equivalence["tested_lipids"]} fitted models match the original delivered trainer on the real corrected data: maximum coefficient error {equivalence["maximum_coefficient_error"]:.3g}, prediction error {equivalence["maximum_prediction_error"]:.3g}; every selected gene set, alpha and L1 ratio also matches. This check complements the repository’s prior full trainer-equivalence audit.','',
           '## Full atlas inference replay','',
           'All 26,639 sample/state pseudobulks and 230,057 metacells were regenerated. Inputs are the stored, deployed normalized RNA H5ADs, aligned by study/sample/canonical state and metacell identity, using the same values for both candidate inference and a production replay. The COPD shard retains its deployed ln1p RNA values, while the other stored inputs use log2 CP10k; this mixed input history is preserved for parity rather than silently changed. Candidate matrix multiplication is an exact affine representation of the fitted ElasticNet models and is checked against the standard API.','',
           f'Maximum affine-versus-API prediction error: {inference["max_affine_vs_API_error"]:.4g}. Full production replay errors by study are below; the saved production matrices have float32 storage.','',
           replay[['unit','study','production_replay_rows','max_abs_production_replay_error','mean_abs_production_replay_error']].to_markdown(index=False,floatfmt='.5g'),' ',
           '## Primary comparison: same 47-lipid differential panel','',
           'Both predictors are tested on the same 47 lipid columns with explicit log2 abundance handling. This matches the moderated-test variance pool and BH family, rather than comparing the 202-panel BH correction to the 47-panel correction. A lipid counts once per cohort when any ordinary cell state has a significant change in the measured direction. Both the tissue lipid and prediction must have BH < 0.10. Opposing states do not subtract concordant support. The matched experimental BH-significant subset contains 29 features.','',
           exact.pivot(index=['contrast','unit'],columns='model',values='concordant_union').reset_index().to_markdown(index=False),'',
           '## Adams–Natri intersections','',intersection.to_markdown(index=False),'',
           '## Comparison against the existing deployed analysis','',
           'This comparison retains the deployed 202-panel testing, including its scale heuristic and full-panel BH denominator, as a separate practical baseline. It does not isolate only model-target changes. The prior replay and saved-result parity audit are retained without rerunning the production deployment.','',
           deployed.pivot(index=['contrast','unit'],columns='model',values='concordant_union').reset_index().to_markdown(index=False),'',
           '## Fold thresholds and raw-p intersections','',summary.to_markdown(index=False),'',
           '## Best individual cell types by concordant count','',best.to_markdown(index=False),'',
           '## Fold preservation and magnitude checks','',
           f'The corrected reference calibration preserves the earlier native population log2 effects to {audit.difference.abs().max():.3g}; maximum raw-p and BH changes are {abs(audit.native_pvalue-audit.combined_pvalue).max():.3g} and {abs(audit.native_FDR-audit.combined_FDR).max():.3g}, respectively. This validates reference differential preservation, not out-of-sample disease fold magnitude.','',
           f'Training-reference arithmetic fold estimates are tabulated for all 470 lipid/population comparisons. Mean absolute prediction error is {fold.effect_log2_error.abs().mean():.4g} log2-effect units. These are in-sample checks and must not be described as independent validation.','',
           'Candidate fold differences are ratios of arithmetic means after exact 2**prediction inversion. Positive signed folds indicate increases; negative folds indicate reciprocal decreases. No abundance pseudocount is added. Relative MS1 intensity units do not establish molar concentrations. Experimental log base remains unconfirmed; experimental effect signs and supplied statistics support directional testing, not confirmed magnitude agreement.','',
           '## Statistical methods and limits','',
           'All seven saved IPF/control contrasts and override comparisons are processed with the frozen deployed cellHarmony functions: pseudobulk moderated t-test (shrinkage=0.2, minimum three rows per arm), metacell Wilcoxon (tie_correct=False, minimum five rows per arm), no metacell cap and no fold cutoff in the primary test. Candidate and matched-panel comparator use explicit log2 metadata to prevent another normalization of log abundance. For differential calculations, candidate output is cast to the deployed float32 storage precision. Statistical tests and thresholds are unchanged. Underpowered pseudobulk comparisons cannot yield significance counts. Natri samples can repeat donors, fibrosis strata overlap donors, and metacell tests use metacells rather than independent donors. The union over states and best-state selection remain exploratory.','',
           '## Every matched lipid, cell-state fold and significance','']
    for (contrast,unit,model),g in tests[tests.model.isin([LABEL,'current_matched47_explicit_log2'])].groupby(['contrast','unit','model']):
        lines += [f'### {contrast} / {unit} / {model}', '',g[['population','gene','measured_feature','signed_fold','pval','fdr','supplied_raw_p','supplied_adjusted_p','direction_concordant_supplied']].to_markdown(index=False,floatfmt='.5g'),'']
    lines += ['## Reproduction','',
              '```bash','/tmp/rna2lipid_elasticnet161/bin/python components/rna2lipid/evaluate_lungmap_elasticnet_ms1.py --train-only',
              '/tmp/rna2lipid_elasticnet161/bin/python components/rna2lipid/infer_lungmap_elasticnet_ms1.py',
              '/opt/homebrew/bin/python3.11 components/rna2lipid/differential_lungmap_elasticnet_ms1.py',
              '/opt/homebrew/bin/python3.11 components/rna2lipid/report_lungmap_elasticnet_ms1.py','```','']
    integrity_notice = '> **INCOMPLETE — does not satisfy the requested production-model comparison.** The user did not authorize reducing the production panel or training roster. The required comparison must retain all 202 lipid outputs and all 50 original training profiles with the instructed deployed ElasticNetCV procedure. Results below are retained as diagnostic history; they must not be used to recommend this candidate or claim completion. Feature/sample mismatches remain unresolved failures, not grounds for exclusion.\n\n'
    (out/'ELASTICNET_MS1_EVALUATION_REPORT.md').write_text(integrity_notice + '\n'.join(lines))
    short=['# Same-ElasticNet MASSIVE-corrected evaluation','',
           'The new candidate uses the deployed lipid-wise ElasticNetCV procedure and hyperparameters, all 1,303 training genes, sorted-cell-only fitting, unchanged deployed RNA inputs and zero-fill inference. No production artifact was replaced. The earlier kernel-ridge candidates do not constitute this comparison.','',
           'The 47-lipid candidate is a restricted evaluation panel selected by raw-MS quality criteria; 155 production lipids were excluded from this candidate, not deleted from production. The broader reconstructed reference has 219 lipid features, including names matching 163 of the 202 production outputs: 47 selected here, 116 additional matching outputs outside this raw-MS-supported panel, and 39 production outputs without a matching name in that reference. Name matching alone does not establish equivalent calibration or resolve the 39 unmatched identities. All 1,303 RNA input genes remain. The original training reference contains 50 profiles from ten donors, each contributing five sorted populations. D071 is one of those training donors, separate from the atlas evaluation cohorts. Its five profiles are absent from the reconstructed target table used here, leaving 45 profiles from nine donors. The reason for that absence has not been established. This comparison changes target coverage and training donor coverage as well as target values, so it does not isolate scaling alone.','',
           '## Same-panel concordant unions, both BH < 0.10','',exact.pivot(index=['contrast','unit'],columns='model',values='concordant_union').reset_index().to_markdown(index=False),'',
           '## Versus deployed baseline','',deployed.pivot(index=['contrast','unit'],columns='model',values='concordant_union').reset_index().to_markdown(index=False),'',
           '[Full report, all lipids, fold-threshold checks and methods](ELASTICNET_MS1_EVALUATION_REPORT.md)','',
           '[Evaluation bundle](elasticnet_MS1_47/candidate_bundle.pkl)','']
    (out/'EVALUATION_SUMMARY.md').write_text(integrity_notice + '\n'.join(short))

if __name__=='__main__':report()
