"""Prepare read-only sample-level Chat inputs from the existing real pseudobulks.

No differential is computed. Clinical fields are joined by study and sample;
only opaque donor labels and analysis covariates are included in the viewer store.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad


def export(root,clinical,crosswalk):
    root=Path(root);m=json.loads((root/'manifest.json').read_text());a=ad.read_h5ad(m['source_pseudobulk'],backed='r')
    obs=a.obs.copy();clinical=pd.read_csv(clinical,dtype={'Sample':str,'Study_internal':str})
    keys=['Study_internal','Sample']
    if clinical.duplicated(keys).any():raise ValueError('Clinical study/sample identifiers are not unique')
    if obs.duplicated(keys+['cell_state']).any():raise ValueError('Expected one real sample/state pseudobulk')
    covs=[c for c in clinical if c in {'age_years','bmi','weight_kg','sex','copd_status','cancer_status','copd_gold_stage','gold_ordinal','smoking_status','smoking_pack_years','sampling_method','fev1_liters','fev1_percent_predicted','fev1_fvc_ratio','fvc_percent_predicted','dlco_percent_predicted'}]
    merged=obs.reset_index(drop=True).merge(clinical[keys+covs],on=keys,how='left',validate='many_to_one',indicator=True)
    if (merged['_merge']!='both').any():raise ValueError('Some pseudobulks have no clinical sample match')
    cross=pd.read_csv(crosswalk,sep='\t');names=dict(zip(cross.cell_state,cross.cell_type));names.update(zip(cross.cell_type_full,cross.cell_type));names.update(zip(cross.cell_type,cross.cell_type))
    merged['cell_state']=merged.cell_state.astype(str).map(names)
    if merged.cell_state.isna().any():raise ValueError('Unmapped states')
    identities=merged.Study_internal.astype(str)+'|'+merged.Sample.astype(str)
    opaque={k:'sample_'+str(i+1).zfill(3) for i,k in enumerate(sorted(set(identities)))}
    merged['donor']=identities.map(opaque)
    output=root/'chat_pseudobulk';output.mkdir(exist_ok=True)
    merged[['donor','cell_state','n_cells']+covs].to_csv(output/'observations.tsv',sep='\t',index=False)
    np.save(output/'genes.npy',a.var_names.to_numpy(dtype=str),allow_pickle=False)
    matrix=a.X.to_memory();matrix=matrix.toarray() if hasattr(matrix,'toarray') else np.asarray(matrix)
    if a.uns.get('normalization')!='ln1p(counts per 10,000)':raise ValueError('Unexpected pseudobulk scale')
    np.save(output/'expression.npy',np.asarray(matrix/np.log(2),dtype=np.float32),allow_pickle=False);a.file.close()
    m['chat_pseudobulk']='chat_pseudobulk';m['chat_pseudobulk_scale']='log2(1+CP10k) of real sample count sums';m['chat_clinical_source']=str(clinical.shape)+' exact study/sample join; opaque sample labels'
    (root/'manifest.json').write_text(json.dumps(m,indent=2));print(len(opaque),'samples;',len(merged),'sample/state rows;',len(covs),'covariates')

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','clinical','crosswalk'):p.add_argument('--'+k,required=True)
    export(**vars(p.parse_args()))
