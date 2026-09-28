#!/usr/bin/env python3
"""Add uncensored COPD donor-pseudobulk DE for states needed by cross-disease claims."""
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import ttest_ind

OUT=Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927/extended/tables')
H5=Path('/Users/saljh8/Dropbox/Transfer/CellRef2.0/reference_v7_top500/pseudobulk/CellRef2.0_pseudobulk_library_x_cellstate_UNION.h5ad')
META=Path('/Users/saljh8/Dropbox/LungMAP/refactored_website/lungmap-data/cellref2/methods/CellRef2.0_harmonized_library_metadata.txt')
CLIN=Path('/Users/saljh8/Dropbox/Transfer/CellRef2.0/Metadata/COPD-full/sources/clinFHH4__sample_name.tsv')
PRES=Path('/Users/saljh8/Dropbox/Transfer/CellRef2.0/reference_v7_top500/pseudobulk/gene_presence_by_study.tsv')
STATES={'Interstitial macrophage':'IM','Peribronchial fibroblast':'PBFB','Subpleural fibroblast':'SPFB',
        'Interstitial fibroblast':'IntFB','Multiciliated (Trachea & Bronchus)':'Ciliated-Bronch'}
def bh(p):
    p=np.asarray(p,float); out=np.full(len(p),np.nan); ok=np.isfinite(p); x=p[ok]
    if not len(x): return out
    ii=np.argsort(x); q=np.minimum.accumulate((x[ii]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1]
    z=np.empty_like(q); z[ii]=np.minimum(q,1); out[np.flatnonzero(ok)]=z; return out

h=ad.read_h5ad(H5,backed='r'); obs=h.obs.reset_index(names='row_key'); genes=np.asarray(h.var_names.astype(str))
m=pd.read_csv(META,sep='\t',dtype=str,low_memory=False)[['Study_internal','Library','Donor','disease','disease_acronym','age_years']].drop_duplicates(['Study_internal','Library'])
o=obs.merge(m,left_on=['Study','Library'],right_on=['Study_internal','Library'],how='left')
o['group']=np.where(o.disease_acronym.eq('COPD'),'COPD',np.where(o.disease.eq('normal'),'Control',None)); o['donor']=o.Donor
c=pd.read_csv(CLIN,sep='\t',dtype=str); c['group_clin']=np.where(c.Group.str.startswith('GOLD',na=False),'COPD',np.where(c.Group.eq('Normal spirometry'),'Control',None)); c['donor_clin']=c.sample3.fillna(c.sample_name)
o=o.merge(c[['Library','group_clin','donor_clin']],on='Library',how='left'); full=o.Study.eq('COPD-full'); o.loc[full,'group']=o.loc[full,'group_clin']; o.loc[full,'donor']=o.loc[full,'donor_clin']
age=pd.to_numeric(o.age_years,errors='coerce'); o=o[(o.Study.eq('COPD-full')|age.ge(18))&o.group.isin(['COPD','Control'])]
present=pd.read_csv(PRES,sep='\t',index_col=0); out=[]
for state,short in STATES.items():
    so=o[(o.cell_state==state)&(o.n_cells.astype(int)>=10)].copy()
    for study in sorted(so.Study.unique()):
        z=so[so.Study==study].copy()
        if z.group.nunique()<2: continue
        X=h.X[z.index.to_numpy(),:]; X=sp.csr_matrix(X) if not sp.issparse(X) else X
        donors=pd.Index(z.donor.astype(str).unique()); di={v:i for i,v in enumerate(donors)}
        G=sp.csr_matrix((np.ones(len(z)),([di[v] for v in z.donor.astype(str)],np.arange(len(z)))),shape=(len(donors),len(z)))
        C=G@X; lib=np.asarray(C.sum(axis=1)).ravel(); keep=lib>0; C=C[keep]; donors=donors[keep]; lib=lib[keep]
        dg=z.groupby(z.donor.astype(str)).group.first().reindex(donors).to_numpy(); case=dg=='COPD'; ctrl=dg=='Control'
        if case.sum()<3 or ctrl.sum()<3: continue
        Y=np.log2(1+C.multiply(1e4/lib[:,None]).toarray()); fc=Y[case].mean(0)-Y[ctrl].mean(0)
        _,p=ttest_ind(Y[case],Y[ctrl],axis=0,equal_var=False,nan_policy='omit')
        measured=present[study].reindex(genes).fillna(False).to_numpy(bool) if study in present else np.ones(len(genes),bool)
        p=np.where(((Y>0).mean(0)>=.1)&measured,p,np.nan); q=bh(p); ix=np.flatnonzero(np.isfinite(p))
        out.append(pd.DataFrame({'state':state,'state_short':short,'study':study,'gene':genes[ix],'log2fc':fc[ix],
          'pval':p[ix],'fdr':q[ix],'n_case_donors':case.sum(),'n_control_donors':ctrl.sum()}))
        print(short,study,int(case.sum()),int(ctrl.sum()),len(ix),flush=True)
h.file.close(); pd.concat(out,ignore_index=True).to_csv(OUT/'copd_uncensored_additional_state_differentials.tsv.gz',sep='\t',index=False)
