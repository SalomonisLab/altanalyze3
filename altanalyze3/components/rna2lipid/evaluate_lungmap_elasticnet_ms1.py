"""Isolated MASSIVE-corrected candidate using the deployed lipid-wise ElasticNetCV."""

try:
    from .candidate_integrity import reject_retired_workflow
except ImportError:
    from candidate_integrity import reject_retired_workflow

if __name__ == "__main__":
    reject_retired_workflow('evaluate_lungmap_elasticnet_ms1.py')

import argparse
import json
import pickle
import shutil
import sys
import time
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[2]))
sys.path.insert(0,str(HERE))
from altanalyze3.components.rna2lipid.training import fit_sparse_lipidwise_elasticnet
from evaluate_lungmap_ipf import sha,dump,ATLAS

OUT=HERE/'artifacts/LungMAP_ElasticNet_MS1_candidate_20261005'
OLD=HERE/'artifacts/LungMAP_IPF_candidate_20261004'
LABEL='elasticnet_MS1_47'
PRODUCTION=HERE/'rna2lipid_hs_lung_lipidwise_bundle.pkl'


def train(out):
    reject_retired_workflow('evaluate_lungmap_elasticnet_ms1.py:train')
    import sklearn
    if sklearn.__version__!='1.6.1':raise ValueError('Use scikit-learn 1.6.1, matching production pickle metadata')
    directory=out/LABEL;directory.mkdir(parents=True,exist_ok=True)
    p=pickle.load(PRODUCTION.open('rb'))
    x=pd.read_csv(HERE/'data/newrna_cell_clair_filtered_symbol.csv',index_col=0).T
    x=x.T.groupby(level=0).mean().T
    x.index=x.index.str.strip();x.columns=x.columns.str.strip()
    x=x.loc[p['training_samples'],p['X_columns']].apply(pd.to_numeric,errors='coerce')
    x=x.fillna(x.median())
    np.testing.assert_allclose(x.mean(0),p['scaler_x'].mean_,atol=1e-12,rtol=0)
    scale=x.std(0,ddof=0).to_numpy();scale[scale==0]=1
    np.testing.assert_allclose(scale,p['scaler_x'].scale_,atol=1e-12,rtol=0)
    ypath=HERE/'artifacts/MS1_reference_normalization_20261001/combined_reference_MS1_log2.csv'
    y=pd.read_csv(ypath,index_col=0)
    mapping={f:f[:-2] for f in y.columns}
    if not set(mapping.values()).issubset(p['Y_columns']):raise ValueError('Corrected targets are not production lipids')
    y=y.rename(columns=mapping)
    used=[s for s in p['training_samples'] if s in y.index]
    missing=[s for s in p['training_samples'] if s not in y.index]
    # No bulk, synthetic labels, extra study alignment or RNA-panel reduction.
    x=x.loc[used];y=y.loc[used,[f for f in p['Y_columns'] if f in y.columns]]
    if x.isna().any().any() or y.isna().any().any():raise ValueError('Incomplete training measurements')
    settings={k:p[k] for k in ['top_gene_options','l1_ratio_grid','alpha_grid','cv_folds','max_iter','sparsity_penalty','random_seed']}
    settings['n_jobs']=next(iter(p['models'].values()))['model'].n_jobs
    dump(out/'production_training_settings.json',settings)
    x.to_csv(directory/'candidate_training_RNA.csv');y.to_csv(directory/'candidate_training_lipids_log2.csv')
    with threadpool_limits(limits=1):
        pred,b,seconds=fit_sparse_lipidwise_elasticnet(x,y,x,**settings)
    b['training_samples']=used;b['training_metadata']=p['training_metadata'].loc[used]
    b['metadata']={'candidate_only':True,'production_default_changed':False,'algorithm':'SparseLipidwiseElasticNetCV',
                   'target_scale':'log2 reference-normalized MS1 signal','source_corrected_targets':str(ypath),
                   'heldout_samples':0,'training_scope':'sorted cells only, matching production group filter',
                   'missing_corrected_training_profiles':missing,'input_policy':'use deployed RNA inputs unchanged; zero-fill absent input genes',
                   'no_extra_RNA_study_alignment':True}
    with (directory/'candidate_bundle.pkl').open('wb') as f:pickle.dump(b,f)
    pred.to_csv(directory/'training_predictions_log2.csv')
    for k in ['summary','coefficients','candidate_models']:b[k].to_csv(directory/f'{k}.csv',index=False)
    # Compare the unchanged port against the exact delivered function on real corrected data.
    sys.path.insert(0,str(HERE/'validation'))
    from validate_trainer_equivalence import load_delivered_reference
    reference=load_delivered_reference()
    with threadpool_limits(limits=1):
        rp,rb,_=reference.fit_sparse_lipidwise_elasticnet(x_train=x,y_train=y.iloc[:,:2],x_test=x,
                            **{k:v for k,v in settings.items() if k!='n_jobs'})
    errors=[]
    for lipid in rp.columns:
        a=b['models'][lipid];r=rb['models'][lipid]
        if a['genes']!=r['genes']:raise ValueError('Trainer selected different genes')
        errors.append(float(abs(a['model'].coef_-r['model'].coef_).max()))
        if a['model'].alpha_!=r['model'].alpha_ or a['model'].l1_ratio_!=r['model'].l1_ratio_:raise ValueError('Trainer parameters differ')
    error=float(abs(pred[rp.columns]-rp).to_numpy().max())
    if error>1e-10 or max(errors)>1e-10:raise ValueError('Delivered trainer parity failed')
    audit={'production_bundle':str(PRODUCTION),'production_bundle_sha256':sha(PRODUCTION),
           'production_training_samples':len(p['training_samples']),'candidate_training_samples':len(x),
           'missing_corrected_profiles':missing,'input_genes':len(x.columns),'corrected_lipids':len(y.columns),
           'bulk_training_profiles':0,'training_seconds':seconds,'sklearn_version':sklearn.__version__,
           'delivered_trainer_coefficient_error':max(errors),'delivered_trainer_prediction_error':error,
           'corrected_target_sha256':sha(ypath),'source_to_production_lipid_names':mapping,
           'production_X_scaler_reconstruction':'all 1303 means and scales match to 1e-12',
           'algorithm_and_settings_match':True,'training_samples_identical':False,'heldout_samples':0}
    dump(out/'training_method_audit.json',audit)
    print('ElasticNet training complete',audit,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',default=str(OUT));ap.add_argument('--train-only',action='store_true')
    args=ap.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    (out/'provenance').mkdir(exist_ok=True)
    for f in ['training.py','api.py','provenance/train_sparse_lipidwise_delivered_2026-08-19.py']:
        shutil.copyfile(HERE/f,out/'provenance'/Path(f).name)
    train(out)
