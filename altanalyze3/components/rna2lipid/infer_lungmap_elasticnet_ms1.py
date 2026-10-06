"""Use deployed RNA inputs and zero-fill policy for candidate and production replay."""

try:
    from .candidate_integrity import reject_retired_workflow
except ImportError:
    from candidate_integrity import reject_retired_workflow

if __name__ == "__main__":
    reject_retired_workflow('infer_lungmap_elasticnet_ms1.py')

import pickle
import numpy as np
import pandas as pd
import anndata as ad
import h5py
from scipy import sparse
from threadpoolctl import threadpool_limits
from evaluate_lungmap_elasticnet_ms1 import OUT,OLD,LABEL,PRODUCTION,HERE,ATLAS
from evaluate_lungmap_ipf import write_predictions,dump,load_current
from altanalyze3.components.rna2lipid.api import Rna2LipidBundle


def affine(b):
    w=np.zeros((len(b['X_columns']),len(b['Y_columns'])))
    intercept=np.zeros(len(b['Y_columns']))
    ix={g:i for i,g in enumerate(b['X_columns'])}
    for j,lipid in enumerate(b['Y_columns']):
        info=b['models'][lipid];w[[ix[g] for g in info['genes']],j]=info['model'].coef_
        intercept[j]=info['model'].intercept_
    w=w*b['scaler_y'].scale_[None,:]/b['scaler_x'].scale_[:,None]
    c=b['scaler_y'].mean_+intercept*b['scaler_y'].scale_-b['scaler_x'].mean_@w
    return w,c


def run(out=OUT):
    reject_retired_workflow('infer_lungmap_elasticnet_ms1.py:run')
    out.mkdir(exist_ok=True)
    cross=pd.read_csv(OLD/'provenance/crosswalk.tsv',sep='\t')
    state_map={}
    for field in ['cell_type','cell_state','cell_type_full']:
        state_map.update(dict(zip(cross[field],cross.cell_type)))
    def row_keys(a,unit):
        if unit=='MC':
            if not a.obs_names.is_unique:raise ValueError('Non-unique metacell identifiers')
            return a.obs_names
        state=a.obs['cell_state'] if 'cell_state' in a.obs else a.obs['short_name']
        key=a.obs.Study_internal.astype(str)+'|'+a.obs.Sample.astype(str)+'|'+state.map(state_map).astype(str)
        key=pd.Index(key)
        if not key.is_unique:raise ValueError('Non-unique study/sample/state/replicate identity')
        return key
    current=load_current(out,cross)
    candidate=pickle.load((out/LABEL/'candidate_bundle.pkl').open('rb'))
    production=pickle.load(PRODUCTION.open('rb'))
    genes=candidate['X_columns'];assert genes==production['X_columns']
    weights={LABEL:affine(candidate),'production_replay':affine(production)}
    api=Rna2LipidBundle.load(out/LABEL/'candidate_bundle.pkl')
    audits=[];engineerrors=[]
    with threadpool_limits(limits=1):
        for unit,name in [('PB','pseudobulk'),('MC','metacell')]:
            target=current[unit]
            result=np.zeros((target.n_obs,len(candidate['Y_columns'])))
            replay=np.zeros((target.n_obs,len(production['Y_columns'])))
            assigned=np.zeros(target.n_obs,dtype=int)
            paths=[ATLAS/f'inputs/cellref2_v8_{"pbsample_ALL" if unit=="PB" else "mc_persample"}_log2cp10k.h5ad',
                   ATLAS/f'inputs/copd_{"pb" if unit=="PB" else "mc"}_ln1pcp10k.h5ad']
            for k,path in enumerate(paths):
                a=ad.read_h5ad(path,backed='r')
                positions=row_keys(target,unit).get_indexer(row_keys(a,unit))
                if (positions<0).any():raise ValueError('Deployed RNA rows not present in final predictions')
                gpos=a.var_names.get_indexer(genes);valid=gpos>=0
                with h5py.File(path,'r') as f:
                    group=f['X'];indptr=group['indptr'][:]
                    for start in range(0,a.n_obs,2048):
                        stop=min(start+2048,a.n_obs);lo,hi=indptr[start],indptr[stop]
                        block=sparse.csr_matrix((group['data'][lo:hi],group['indices'][lo:hi],indptr[start:stop+1]-lo),shape=(stop-start,a.n_vars))
                        x=np.zeros((stop-start,len(genes)),dtype=float)
                        x[:,valid]=block[:,gpos[valid]].toarray()
                        dest=positions[start:stop]
                        # Replay the actual deployed mixed inputs, including the COPD ln1p shard.
                        for label,(w,c) in weights.items():
                            vals=x@w+c
                            if label==LABEL:result[dest]=vals
                            else:replay[dest]=vals
                        if start==0:
                            real=api.predict_from_dataframe(pd.DataFrame(x,columns=genes)).predictions.to_numpy()
                            err=float(abs(real-result[dest]).max());engineerrors.append(err)
                            if err>1e-10:raise ValueError('Affine inference differs from production API')
                        assigned[dest]+=1
                        if start%20480==0:print(unit,path.name,stop,'/',a.n_obs,flush=True)
                audits.append({'unit':unit,'input_path':str(path),'rows':a.n_obs,'normalization':a.uns.get('normalization'),
                               'missing_input_genes':[genes[i] for i in np.flatnonzero(~valid)]})
                a.file.close()
            if not (assigned>=1).all():raise ValueError('Full atlas inference left missing rows')
            outa=write_predictions(out,LABEL,unit,result,target.obs,candidate['Y_columns'])
            outa.uns['inference_policy']='Deployed stored RNA inputs unchanged; absent genes zero-filled; no extra RNA alignment'
            outa.write_h5ad(out/LABEL/f'LungMAP_{unit}_predictions_log2.h5ad',compression='lzf')
            # Final production values were written float32; assess every replay row by study.
            for study,g in target.obs.groupby('Study_internal',observed=True):
                i=target.obs_names.get_indexer(g.index)
                errors=abs(replay[i]-np.asarray(target.X[i],dtype=float))
                audits.append({'unit':unit,'study':str(study),'production_replay_rows':len(i),
                               'max_abs_production_replay_error':float(errors.max()),'mean_abs_production_replay_error':float(errors.mean())})
            # Store replay as additional evidence, without replacing the deployed comparator.
            ad.AnnData(replay,obs=target.obs.copy(),var=target.var.copy(),uns={'expression_scale':'production_log2_predictions'}).write_h5ad(out/LABEL/f'production_replay_{unit}.h5ad',compression='lzf')
    dump(out/'inference_method_audit.json',{'inputs_and_replay':audits,'max_affine_vs_API_error':max(engineerrors),
                                        'candidate_bundle':str(out/LABEL/'candidate_bundle.pkl'),'extra_RNA_alignment':False,
                                        'all_1303_training_genes_retained':True,'zero_fill_absent_input_genes':True})
    print('Complete full atlas inference',flush=True)

if __name__=='__main__':run()
