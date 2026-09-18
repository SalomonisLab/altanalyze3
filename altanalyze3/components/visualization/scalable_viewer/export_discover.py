"""Read-only export of LungMAP pseudobulk results and real sample arm means.

No differential engine is invoked. The source database and pseudobulks are opened
read-only. The output is a self-contained evidence store for integrated views.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sqlite3
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse


def export(database, pseudobulk, covariates, crosswalk, out, resume=False):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    db=Path(database)
    stamp=(db.stat().st_size,db.stat().st_mtime_ns)
    source=sqlite3.connect(f"file:{db}?mode=ro",uri=True)
    contrasts=pd.read_sql_query("SELECT * FROM contrast WHERE scope='COPD' ORDER BY contrast_id",source)
    records={}
    datafile=out/'differentials.sqlite'
    if datafile.exists() and not resume:
        raise FileExistsError(f"Refusing to replace {datafile}")
    dest=sqlite3.connect(datafile)
    dest.execute('CREATE TABLE IF NOT EXISTS differential (comparison TEXT, modality TEXT, gene TEXT, population TEXT, log2fc REAL, pval REAL, fdr REAL, n_case INTEGER, n_control INTEGER, case_mean REAL, control_mean REAL, source TEXT, analysis_version TEXT, statistical_method TEXT)')
    coverage=pd.read_sql_query("SELECT * FROM analysis_coverage WHERE contrast_key LIKE 'COPD__%' AND replicate_unit='pseudobulk'",source)
    coverage.to_csv(out/'coverage.tsv',sep='\t',index=False)
    for row in contrasts.to_dict('records'):
        name=row['name'];comparison=name.removeprefix('COPD__')
        params=(row['contrast_id'],)
        query="""SELECT ? comparison,f.modality,f.name gene,p.name population,d.log2fc,d.pval,d.fdr,d.n_case,d.n_control,d.case_mean,d.control_mean,d.source,d.analysis_version,d.statistical_method
          FROM differential d JOIN feature f USING(feature_id) JOIN population p USING(population_id)
          WHERE d.contrast_id=? AND d.replicate_unit='pseudobulk'"""
        counts={}
        for chunk in ([] if resume else pd.read_sql_query(query,source,params=(comparison,*params),chunksize=100000)):
            chunk.to_sql('differential',dest,if_exists='append',index=False)
            for mod,n in chunk.groupby('modality').size().items():counts[mod]=counts.get(mod,0)+int(n)
        if resume:counts=dict(dest.execute('SELECT modality,COUNT(*) FROM differential WHERE comparison=? GROUP BY modality',(comparison,)).fetchall())
        entries=coverage.loc[coverage.contrast_key==name]
        records[comparison]={'source_contrast':name,'case_label':row['case_label'],'control_label':row['control_label'],
                             'modalities':{r.modality:{'status':r.status,'rows':counts.get(r.modality,0),'states':json.loads(r.tested_populations)} for r in entries.itertuples()},'replicate_unit':'pseudobulk'}
        print(name,counts,flush=True)
    dest.execute('CREATE INDEX IF NOT EXISTS by_comparison ON differential(comparison,modality,population)')
    dest.commit();dest.close();source.close()
    a=ad.read_h5ad(pseudobulk,backed='r')
    if a.uns.get('unit')!='pseudobulk':raise ValueError('A real sample pseudobulk input is required')
    keys=a.obs_names.astype(str)
    if keys.duplicated().any():raise ValueError('Duplicate pseudobulk replicate ids')
    if a.obs.duplicated(['Study_internal','Sample','cell_state']).any():raise ValueError('Libraries were not collapsed to samples')
    mapping=pd.read_csv(crosswalk,sep='\t',dtype=str)
    short=dict(zip(mapping.cell_state,mapping.cell_type));short.update(zip(mapping.cell_type_full,mapping.cell_type));short.update(zip(mapping.cell_type,mapping.cell_type))
    state_names=a.obs.cell_state.astype(str).map(short)
    if state_names.isna().any():raise ValueError('Unmapped pseudobulk cell states')
    states=sorted(state_names.unique());genes=a.var_names.astype(str).to_numpy()
    assignments={}
    for comp,r in records.items():
        c=pd.read_csv(Path(covariates)/(r['source_contrast']+'__PB.txt'),sep='\t',dtype=str).set_index('replicate_id').Condition
        if c.index.duplicated().any():raise ValueError('Duplicate arm assignments')
        assignments[comp]=c.reindex(keys).fillna('').to_numpy()
    arrays={comp:{arm:np.full((len(states),len(genes)),np.nan,dtype=np.float32) for arm in ['case_mean_log2cp10k','control_mean_log2cp10k','case_mean_cp10k','control_mean_cp10k']} for comp in records}
    counts={comp:np.zeros((len(states),2),dtype=np.int32) for comp in records}
    scale=str(a.uns.get('normalization',''))
    if scale!='ln1p(counts per 10,000)':raise ValueError('Unexpected source scale: '+scale)
    matrix=a.X.to_memory() if hasattr(a.X,'to_memory') else a.X
    for i,state in enumerate(states):
        mask=state_names.to_numpy()==state
        x=matrix[mask,:]
        x=x.toarray() if sparse.issparse(x) else np.asarray(x)
        for comp,groups in assignments.items():
            for j,(arm,label) in enumerate([('case','CASE'),('control','CONTROL')]):
                v=x[groups[mask]==label]
                counts[comp][i,j]=len(v)
                if len(v):
                    arrays[comp][arm+'_mean_log2cp10k'][i]=np.mean(v/np.log(2),axis=0)
                    arrays[comp][arm+'_mean_cp10k'][i]=np.mean(np.expm1(v),axis=0)
        print('means',state,len(x),flush=True)
    a.file.close()
    for comp,r in records.items():
        filename=comp+'_arm_means.npz'
        np.savez_compressed(out/filename,genes=genes.astype(str),states=np.array(states),n_samples=counts[comp],**arrays[comp])
        r['arm_means']=filename
        r['mean_scale']='mean of sample log2(1+CP10k), and mean sample CP10k stored separately'
    # A convenient, compact table for the requested disease groups.
    disease='COPD_vs_non-COPD'
    frames=[]
    for i,state in enumerate(states):
        for j,arm in enumerate(['case','control']):
            frames.append(pd.DataFrame({'gene':genes,'cell_state':state,'disease_group':records[disease][arm+'_label'],
                           'n_samples':int(counts[disease][i,j]),'mean_log2cp10k':arrays[disease][arm+'_mean_log2cp10k'][i],
                           'mean_cp10k':arrays[disease][arm+'_mean_cp10k'][i]}))
    pd.concat(frames,ignore_index=True).to_csv(out/'COPD_disease_group_expression.tsv.gz',sep='\t',index=False)
    manifest={'database':str(db),'database_size':stamp[0],'database_mtime_ns':stamp[1],
              'source_pseudobulk':str(pseudobulk),'source_covariates':str(covariates),
              'source_crosswalk':str(crosswalk),'source_scale':scale,'n_pseudobulks':len(keys),'n_genes':len(genes),'states':states,
              'replicate_unit':'one sample x cell state; libraries summed before normalization',
              'differential_database':datafile.name,'comparisons':records,'analyses_recomputed':False}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    assert stamp==(db.stat().st_size,db.stat().st_mtime_ns),'Source database changed during export'
    return manifest


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['database','pseudobulk','covariates','crosswalk','out']:ap.add_argument('--'+name,required=True)
    ap.add_argument('--resume',action='store_true',help='Resume means after a completed differential export')
    export(**vars(ap.parse_args()))
