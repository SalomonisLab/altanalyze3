"""Donor-nested ridge/PLS tuning for fold magnitudes without post-hoc inflation."""
import argparse
import json
from pathlib import Path
import pickle
import sys

import numpy as np
import pandas as pd
from sklearn.kernel_ridge import KernelRidge
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import LinearRegression
from threadpoolctl import threadpool_limits

from reference_normalization import fast_predict,fit_offsets,transform,effect_metrics

ALPHAS=np.array([.01,.1,1.,10.,30.,100.,300.,1000.])
ENGINE='ridge'


def candidate_predictions(xtrain,ytrain,xtest,parameters):
    if ENGINE=='ridge':return fast_predict(xtrain,ytrain,xtest,parameters)
    sx,sy=StandardScaler().fit(xtrain),StandardScaler().fit(ytrain)
    model=PLSRegression(n_components=int(max(parameters)),scale=False,max_iter=1000)
    model.fit(sx.transform(xtrain),sy.transform(ytrain))
    scores=model.transform(sx.transform(xtest))
    return [sy.inverse_transform(scores[:,:int(k)]@model.y_loadings_[:,:int(k)].T+model._y_mean) for k in parameters]


def pairs(frame,meta):
    groups=[]
    sorted_meta=meta[meta.dataset=='sorted']
    for kind in ['population','donor']:
        arrays=[]
        for _,md in sorted_meta.groupby(kind):
            a=frame.loc[md.index].to_numpy();i,j=np.triu_indices(len(a),1)
            if len(i):arrays.append(a[i]-a[j])
        groups.append(np.concatenate(arrays))
    return groups


def magnitude_scores(truth,predicted,meta):
    actual=pairs(truth,meta);estimated=pairs(predicted,meta);scores=[]
    for a,b in zip(actual,estimated):
        energy=np.maximum(np.mean(a*a,axis=0),.1)
        slope=np.sum(a*b,axis=0)/np.maximum(np.sum(a*a,axis=0),1e-10)
        # Penalize both prediction error and attenuation. The fixed slope term
        # is chosen before external evaluation; outer held-out donors test it.
        scores.append(np.mean((a-b)**2,axis=0)/energy+.5*(slope-1)**2)
    return np.mean(scores,axis=0)


def summary(a,b):
    flat_a,flat_b=a.ravel(),b.ravel();error=flat_b-flat_a
    return {'contrasts':int(a.size),'correlation':float(np.corrcoef(flat_a,flat_b)[0,1]),
            'RMSE_log2FC':float(np.sqrt(np.mean(error**2))),
            'median_absolute_error_log2FC':float(np.median(abs(error))),
            'OLS_slope':float(np.polyfit(flat_a,flat_b,1)[0]),
            'fraction_within_0_5_log2':float(np.mean(abs(error)<=.5))}


def run(reference,ms1_reference,raw,out):
    reference,ms1_reference,raw,out=map(Path,[reference,ms1_reference,raw,out]);out.mkdir(parents=True,exist_ok=True)
    truth=pd.read_csv(ms1_reference/'combined_reference_MS1_log2.csv',index_col=0)
    X=pd.read_csv(reference/'paired_RNA.csv',index_col=0)
    meta=pd.read_csv(reference/'sample_metadata.csv',index_col=0)
    native=pd.read_csv(reference/'native_lipid_targets.csv',index_col=0)[truth.columns]
    cell=native.loc[meta.dataset=='sorted']
    bulk=pd.read_csv(reference/'legacy_preprocessing_audit/bulk_native_log2_mode_preserved.csv',index_col=0)[truth.columns]
    full_offsets=fit_offsets('lipid_pmx_shift',cell,bulk)
    signal_shift=pd.read_csv(ms1_reference/'calibration_offsets.csv',index_col=0).added_common_log2_offset
    def predict(train,test,excluded,alphas):
        offsets=fit_offsets('lipid_pmx_shift',cell,bulk,excluded)
        y=transform(native,meta,offsets,'lipid_pmx_shift')
        predictions=candidate_predictions(X.loc[train].to_numpy(),y.loc[train].to_numpy(),X.loc[test].to_numpy(),alphas)
        answer=[]
        for p in predictions:
            frame=pd.DataFrame(p,index=meta.index[test],columns=truth.columns)
            frame.loc[meta.loc[frame.index].dataset=='sorted']+=full_offsets-offsets
            frame+=signal_shift
            answer.append(frame)
        return answer
    def tune(train,excluded):
        predictions=[pd.DataFrame(np.nan,index=meta.index[train],columns=truth.columns) for _ in ALPHAS]
        for held in sorted(meta.loc[train].donor.unique()):
            inner_train=train&(meta.donor!=held);inner_test=train&(meta.donor==held)
            frames=predict(inner_train,inner_test,set(excluded)|{held},ALPHAS)
            for target,frame in zip(predictions,frames):target.loc[frame.index]=frame
        scores=np.stack([magnitude_scores(truth.loc[train],p,meta.loc[train]) for p in predictions])
        return np.argmin(scores,axis=0),scores
    oof=pd.DataFrame(np.nan,index=truth.index,columns=truth.columns);baseline=oof.copy();records=[]
    for held in sorted(meta.donor.unique()):
        train=meta.donor!=held;test=~train;chosen,scores=tune(train,{held})
        predictions=predict(train,test,{held},ALPHAS)
        matrix=np.stack([p.to_numpy() for p in predictions])
        oof.loc[test]=matrix[chosen,:,np.arange(len(chosen))].T
        baseline.loc[test]=predictions[list(ALPHAS).index(100. if ENGINE=='ridge' else 5.)]
        records.append({'held_out_donor':held,'parameters_per_lipid':dict(zip(truth.columns,ALPHAS[chosen].tolist()))})
        print('nested donor',held,flush=True)
    oof.to_csv(out/'OOF_magnitude_tuned_log2.csv')
    baseline.to_csv(out/('OOF_alpha100_log2.csv' if ENGINE=='ridge' else 'OOF_PLS_fixed5_log2.csv'))
    reports={};tables=[]
    for label,p in [('magnitude_tuned',oof),('alpha100_baseline' if ENGINE=='ridge' else 'PLS_fixed5_baseline',baseline)]:
        _,arrays=effect_metrics(truth,p,meta);reports[label]={}
        for kind,(a,b) in arrays.items():
            reports[label][kind]=summary(a,b)
            for j,lipid in enumerate(truth.columns):
                tables.append({'model':label,'comparison':kind,'lipid':lipid,**summary(a[:,j],b[:,j])})
    per_lipid=pd.DataFrame(tables);per_lipid.to_csv(out/'per_lipid_fold_magnitude_validation.csv',index=False)
    chosen,scores=tune(pd.Series(True,index=meta.index),set())
    sx,sy=StandardScaler().fit(X),StandardScaler().fit(truth)
    if ENGINE=='ridge':
        model=KernelRidge(alpha=ALPHAS[chosen],kernel='linear').fit(sx.transform(X),sy.transform(truth))
    else:
        pls=PLSRegression(n_components=int(max(ALPHAS[chosen])),scale=False,max_iter=1000).fit(sx.transform(X),sy.transform(truth))
        coefs=np.stack([pls.x_rotations_[:,:int(k)]@pls.y_loadings_[j,:int(k)] for j,k in enumerate(ALPHAS[chosen])])
        model=LinearRegression();model.coef_=coefs
        model.intercept_=pls._y_mean-pls._x_mean@coefs.T;model.n_features_in_=X.shape[1]
    # Export descriptive support tiers; magnitude support is distinct from the
    # older direction-based support list and does not prove disease accuracy.
    selected=per_lipid[per_lipid.model=='magnitude_tuned']
    tiers={}
    for kind,d in selected.groupby('comparison'):
        tiers[kind]=d.loc[(d.correlation>=.7)&d.OLS_slope.between(.8,1.2)&
                         (d.RMSE_log2FC<=1)&(d.median_absolute_error_log2FC<=.5),'lipid'].tolist()
    both=sorted(set.intersection(*(set(v) for v in tiers.values())))
    bundle={'model':model,'scaler_x':sx,'scaler_y':sy,'X_columns':X.columns.tolist(),'Y_columns':truth.columns.tolist(),
            'metadata':{'expression_scale':'log2','target_scaling':{'mode':'standard'},
                        'model_name':'NestedFoldMagnitudeKernelRidge' if ENGINE=='ridge' else 'NestedFoldMagnitudePLS',
                        'normalization_method':'MS1_common_lipid_shift',
                        'engine':ENGINE,'parameters_per_lipid':dict(zip(truth.columns,ALPHAS[chosen].tolist())),
                        'internally_supported_lipids':both,'magnitude_support_tiers':tiers,
                        'calibration_interpretation':'Reference-normalized MS1 peak intensity; measured reference folds unchanged',
                        'validation':'Donor-nested parameter selection using RNA training donors only; no post-hoc inverse-slope multiplication',
                        'disease_fold_accuracy':'Not established','production_default_changed':False}}
    with (out/'magnitude_candidate_bundle.pkl').open('wb') as handle:pickle.dump(bundle,handle)
    # MassIVE donors' bulk profiles must be excluded from fitting and tuning.
    excluded={'D001','D008','D011'};train=~meta.donor.isin(excluded);ext_chosen,_=tune(train,excluded)
    offsets=fit_offsets('lipid_pmx_shift',cell,bulk,excluded)
    ext_y=transform(native,meta,offsets,'lipid_pmx_shift')
    with (ms1_reference/'MS1_reference_bundle.pkl').open('rb') as handle:old=pickle.load(handle)
    original_rna=pd.read_csv('/Users/saljh8/Downloads/Lipidomics/newrna_cell_clair_filtered_symbol.csv',index_col=0).T
    original_rna=original_rna.T.groupby(level=0).mean().T
    old_input=original_rna.loc[:,X.columns]
    ext_indices=[s for s in old_input.index if s.startswith(('D008_','D011_')) and '.' not in s]
    p=candidate_predictions(X.loc[train].to_numpy(),ext_y.loc[train].to_numpy(),old_input.loc[ext_indices].to_numpy(),ALPHAS)
    matrix=np.stack(p);external=pd.DataFrame(matrix[ext_chosen,:,np.arange(len(ext_chosen))].T,index=ext_indices,columns=truth.columns)
    measured=pd.read_csv(raw/'RNA_prediction_vs_reextracted_MS1.csv')
    rows=[]
    for r in measured.to_dict('records'):
        if r['lipid'] not in truth:continue
        r['candidate_log2FC']=float(external.loc[r['sample'],r['lipid']]-external.loc[r['reference'],r['lipid']]);rows.append(r)
    ext=pd.DataFrame(rows);ext.to_csv(out/'external_MS1_fold_comparison.csv',index=False)
    old_supported=old['metadata']['internally_supported_lipids']
    for name,d in [('all47',ext),('previous11',ext[ext.lipid.isin(old_supported)])]:
        reports['external_'+name]={model_name:summary(d.MS1_reextracted_log2FC.to_numpy(),d[column].to_numpy())
                                  for model_name,column in [('magnitude_tuned','candidate_log2FC'),('previous_model','predicted_log2fc')]}
    full=sy.inverse_transform(model.predict(sx.transform(old_input)))
    predictions=pd.DataFrame(full,index=old_input.index,columns=truth.columns)
    predictions.to_csv(out/'exploratory_predictions_log2.csv');np.exp2(predictions).to_csv(out/'exploratory_predictions_linear.csv')
    predictions[both].to_csv(out/'both_contrast_tiers_supported_predictions_log2.csv')
    np.exp2(predictions[both]).to_csv(out/'both_contrast_tiers_supported_predictions_linear.csv')
    reports['support_tiers']=tiers;reports['both_tiers_supported']=both
    reports['engine']=ENGINE
    reports['status']='Experimental candidate; report validation before selecting a default; no post-hoc fold inflation'
    (out/'fold_magnitude_model_result.json').write_text(json.dumps(reports,indent=2)+'\n')
    (out/'outer_fold_parameters.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(reports,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',required=True);p.add_argument('--ms1-reference',required=True)
    p.add_argument('--raw',required=True);p.add_argument('--out',required=True)
    p.add_argument('--engine',choices=['ridge','pls'],default='ridge');a=p.parse_args()
    ENGINE=a.engine
    if ENGINE=='pls':ALPHAS=np.array([1.,2.,3.,5.,8.,12.,16.,20.])
    with threadpool_limits(limits=1):run(a.reference,a.ms1_reference,a.raw,a.out)
