"""Transfer a differential-preserving combined reference to an MS1 signal gauge.

This is a common per-lipid multiplicative calibration, not concentration
calibration. Raw replication QC and model support remain separate requirements.
"""

try:
    from .candidate_integrity import reject_retired_workflow
except ImportError:
    from candidate_integrity import reject_retired_workflow

if __name__ == "__main__":
    reject_retired_workflow('build_raw_scale_reference.py')

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from reference_normalization import fit_bundle, preservation_gate, population_differentials


def run(reference_dir,raw_dir,out):
    reject_retired_workflow('build_raw_scale_reference.py:run')
    reference=Path(reference_dir);rawdir=Path(raw_dir);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    quality=pd.read_csv(rawdir/'per_lipid_replication.csv')
    quality['raw_replication_support']=(quality.isobaric_targets.isna()&(quality.ms2_support_scans>0)&
        (quality.observed_profiles==15)&(quality.rmse_log2<=.75)&(quality.spearman>=.7))
    quality['lipid']=quality.feature_id.str.split('|').str[0]+quality.feature_id.map(lambda s:'_P' if s.endswith('+') else '_N')
    original=pd.read_csv(reference/'combined_training_reference_log2.csv',index_col=0)
    matched=quality[quality.raw_replication_support&quality.lipid.isin(original.columns)].copy()
    if matched.lipid.duplicated().any():raise ValueError('Ambiguous annotation bridge')
    matched=matched.set_index('lipid').sort_index();features=matched.index.tolist()
    measured=pd.read_csv(rawdir/'raw_apex_log2.csv',index_col=0)
    sample_medians=measured.median(axis=1)
    # Fix the overall signal gauge using one actual PMX sample. Equalize the
    # other measured sample medians before estimating per-lipid PMX baselines.
    signal_gauge=float(sample_medians.loc['D001_PMX'])
    normalized=measured.sub(sample_medians,axis=0)+signal_gauge
    metadata=pd.read_csv(reference/'sample_metadata.csv',index_col=0)
    pmx=(metadata.dataset=='sorted')&(metadata.population=='PMX')
    source_baseline=original.loc[pmx,features].median(axis=0)
    raw_baseline=pd.Series({lipid:normalized.loc[['D001_PMX','D008_PMX','D011_PMX'],row.feature_id].median()
                            for lipid,row in matched.iterrows()})
    shifts=raw_baseline-source_baseline
    combined=original[features]+shifts
    combined.to_csv(out/'combined_reference_MS1_log2.csv');linear=np.exp2(combined)
    linear.to_csv(out/'combined_reference_MS1_linear.csv')
    if not np.all(np.isfinite(linear.to_numpy())) or not np.all(linear.to_numpy()>0):raise ValueError('Invalid linear abundance')
    reloaded=pd.read_csv(out/'combined_reference_MS1_log2.csv',index_col=0)
    gate=preservation_gate(original[features],reloaded,metadata)
    if not gate['pass']:raise ValueError('Raw signal calibration changed lipid differentials')
    tests=population_differentials(original[features],reloaded,metadata,'MS1_common_lipid_shift')
    tests.to_csv(out/'population_differential_checks.csv',index=False)
    p_error=float((tests.native_pvalue-tests.combined_pvalue).abs().max())
    fdr_error=float((tests.native_FDR-tests.combined_FDR).abs().max())
    native_calls=(tests.native_FDR<.05)&(tests.native_log2fc.abs()>=np.log2(1.5))
    combined_calls=(tests.combined_FDR<.05)&(tests.combined_log2fc.abs()>=np.log2(1.5))
    changed_calls=int((native_calls!=combined_calls).sum())
    if p_error>1e-9 or fdr_error>1e-9 or changed_calls:
        raise ValueError('MS1 signal calibration changed a statistical differential')
    pd.DataFrame({'combined_reference_PMX_log2':source_baseline,'MS1_reference_PMX_log2':raw_baseline,
                  'added_common_log2_offset':shifts}).to_csv(out/'calibration_offsets.csv')
    validation=pd.read_csv(reference/'per_lipid_validation.csv',index_col=0)
    support=[f for f in features if bool(validation.loc[f,'internal_support'])]
    X=pd.read_csv(reference/'paired_RNA.csv',index_col=0)
    result=json.loads((reference/'workflow_result.json').read_text())
    alpha=float(result.get('final_alpha',100))
    model,sx,sy=fit_bundle(X,combined,'MS1_common_lipid_shift',shifts,alpha,out/'MS1_reference_bundle.pkl',support)
    predictions=sy.inverse_transform(model.predict(sx.transform(X)))
    prediction=pd.DataFrame(predictions,index=X.index,columns=features)
    prediction.to_csv(out/'training_predictions_MS1_log2.csv');np.exp2(prediction).to_csv(out/'training_predictions_MS1_linear.csv')
    supplied=pd.read_csv(reference/'supplied_sorted_RNA_predictions_log2.csv',index_col=0)
    shifted=supplied[features]+shifts
    shifted.to_csv(out/'all_raw_supported_predictions_log2.csv');np.exp2(shifted).to_csv(out/'all_raw_supported_predictions_linear.csv')
    shifted[support].to_csv(out/'internally_and_raw_supported_predictions_log2.csv')
    np.exp2(shifted[support]).to_csv(out/'internally_and_raw_supported_predictions_linear.csv')
    # Target calibration is only an output translation; refitting must retain
    # the established kernel-ridge predictions to numerical precision.
    old=pd.read_csv(reference/'training_predictions_log2.csv',index_col=0)
    prediction_error=float(np.max(np.abs((prediction-old[features]-shifts).to_numpy())))
    if prediction_error>1e-8:raise ValueError('Refitted calibrated predictor differs from original regression')
    import pickle
    path=out/'MS1_reference_bundle.pkl'
    with path.open('rb') as handle:bundle=pickle.load(handle)
    bundle['metadata'].update({'calibration_interpretation':'Reference-normalized MS1 peak intensity; median signal gauge D001_PMX; molar concentration not established',
                              'raw_replication_supported_lipids':features,'internally_supported_lipids':support,
                              'raw_replication_QC':'15 profiles observed, unambiguous published mass, MS2 support, RMSE<=0.75 log2 and Spearman>=0.7',
                              'raw_reference_donors':['D001','D008','D011'],
                              'biological_reference_assumption':'Median PMX lipid composition is transferable between the two cohorts; offsets do not establish equal specimen amounts'})
    with path.open('wb') as handle:pickle.dump(bundle,handle)
    matched.to_csv(out/'raw_replication_support.csv')
    summary={'paired_profiles':len(combined),'raw_replication_supported_features':len(features),
             'raw_and_internal_model_supported_features':len(support),'supported_lipids':support,
             'signal_gauge_reference':'D001_PMX sample median MS1 log2 peak intensity','signal_gauge_log2':signal_gauge,
             'PMX_baseline':'median of three raw PMX profiles after sample-median equalization',
             'calibration_rule':'Same per-lipid offset added to every sorted and bulk reference profile',
             'differential_gate':gate,'maximum_refitted_prediction_change_after_offset':prediction_error,
             'population_tests':len(tests),'max_pvalue_change':p_error,'max_FDR_change':fdr_error,
             'changed_significance_calls':changed_calls,'significant_tests':int(native_calls.sum()),
             'default_prediction_policy':'Only features passing both raw replication QC and internal RNA model support',
             'units':'reference-normalized MS1 peak intensity; no pseudocount',
             'QC_scope':'Descriptive replication screen; does not establish independent disease fold accuracy',
             'production_default_changed':False}
    (out/'raw_scale_reference_result.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-dir',required=True);parser.add_argument('--raw-dir',required=True);parser.add_argument('--out',required=True)
    a=parser.parse_args();run(a.reference_dir,a.raw_dir,a.out)
