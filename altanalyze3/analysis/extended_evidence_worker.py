#!/usr/bin/env python3
"""Uncensored COPD pseudobulk differential programs and comprehensive drug reversal scan."""
from __future__ import annotations
import json, os, subprocess, sys
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import ttest_ind

ROOT=Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927')
OUT=ROOT/'extended'; TAB=OUT/'tables'; ENG=OUT/'drug_engine'
for d in [OUT,TAB,ENG]: d.mkdir(parents=True,exist_ok=True)
H5=Path('/Users/saljh8/Dropbox/Transfer/CellRef2.0/reference_v7_top500/pseudobulk/CellRef2.0_pseudobulk_library_x_cellstate_UNION.h5ad')
META=Path('/Users/saljh8/Dropbox/LungMAP/refactored_website/lungmap-data/cellref2/methods/CellRef2.0_harmonized_library_metadata.txt')
CLIN=Path('/Users/saljh8/Dropbox/Transfer/CellRef2.0/Metadata/COPD-full/sources/clinFHH4__sample_name.tsv')
PRES=Path('/Users/saljh8/Dropbox/Transfer/CellRef2.0/reference_v7_top500/pseudobulk/gene_presence_by_study.tsv')
SRC=Path('/Users/saljh8/Dropbox/LungMAP/Discovery/hypotheses_20260927')
FOLD=np.log2(1.2)
STATES=['Alveolar type 2','Alveolar type 1','Alveolar epithelial transitional (AT2-AT1)',
        'Alveolar macrophage','Alveolar macrophage (lipid homeostatic)',
        'Transitioning monocyte-derived macrophage','Pulmonary vein endothelial','Capillary 1',
        'Airway basal','Club (Trachea & Bronchus)','Adventitial fibroblast','Mesothelial']
SHORT={'Alveolar type 2':'AT2','Alveolar type 1':'AT1','Alveolar epithelial transitional (AT2-AT1)':'AT2_AT1_int',
       'Alveolar macrophage':'AM','Alveolar macrophage (lipid homeostatic)':'AM_lipid',
       'Transitioning monocyte-derived macrophage':'tMDM','Pulmonary vein endothelial':'PVEC','Capillary 1':'CAP1',
       'Airway basal':'Basal','Club (Trachea & Bronchus)':'Club_TB','Adventitial fibroblast':'AdvFB','Mesothelial':'Mesothelial'}

def bh(p):
    p=np.asarray(p,float); o=np.full(len(p),np.nan); ok=np.isfinite(p); x=p[ok]
    if not len(x): return o
    ii=np.argsort(x); q=np.minimum.accumulate((x[ii]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1]
    z=np.empty_like(q); z[ii]=np.minimum(q,1); o[np.flatnonzero(ok)]=z; return o

h=ad.read_h5ad(H5,backed='r'); obs=h.obs.reset_index(names='row_key'); genes=np.asarray(h.var_names.astype(str))
m=pd.read_csv(META,sep='\t',dtype=str,low_memory=False)
m=m[['Study_internal','Library','Donor','disease','disease_acronym','age_years']].drop_duplicates(['Study_internal','Library'])
o=obs.merge(m,left_on=['Study','Library'],right_on=['Study_internal','Library'],how='left')
o['group']=np.where(o.disease_acronym.eq('COPD'),'COPD',np.where(o.disease.eq('normal'),'Control',None))
o['donor']=o.Donor
# The uncensored GSE310058 clinical table supplies real COPD phenotypes missing from harmonized metadata.
c=pd.read_csv(CLIN,sep='\t',dtype=str)
c['group_clin']=np.where(c.Group.str.startswith('GOLD',na=False),'COPD',np.where(c.Group.eq('Normal spirometry'),'Control',None))
c['donor_clin']=c.sample3.fillna(c.sample_name)
o=o.merge(c[['Library','group_clin','donor_clin','Group','Smoking Status']],on='Library',how='left')
isfull=o.Study.eq('COPD-full'); o.loc[isfull,'group']=o.loc[isfull,'group_clin']; o.loc[isfull,'donor']=o.loc[isfull,'donor_clin']
# Adult-only for legacy cohorts; GSE310058 is adult and uses smoking-normal controls.
age=pd.to_numeric(o.age_years,errors='coerce'); o=o[(o.Study.eq('COPD-full') | age.ge(18)) & o.group.isin(['COPD','Control'])]
present=pd.read_csv(PRES,sep='\t',index_col=0)

all_de=[]
for state in STATES:
    so=o[(o.cell_state==state)&(o.n_cells.astype(int)>=10)].copy()
    for study in sorted(so.Study.unique()):
        z=so[so.Study==study].copy()
        if z.group.nunique()<2: continue
        rows=z.index.to_numpy(); X=h.X[rows,:]
        if not sp.issparse(X): X=sp.csr_matrix(X)
        # Collapse multiple lobes/libraries from the same donor before inference.
        donors=pd.Index(z.donor.astype(str).unique()); di={x:i for i,x in enumerate(donors)}
        G=sp.csr_matrix((np.ones(len(z)),([di[x] for x in z.donor.astype(str)],np.arange(len(z)))),shape=(len(donors),len(z)))
        C=G@X; lib=np.asarray(C.sum(axis=1)).ravel(); keep=lib>0
        C=C[keep]; donors=donors[keep]; lib=lib[keep]
        dg=z.groupby(z.donor.astype(str)).group.first().reindex(donors).to_numpy()
        case=dg=='COPD'; ctrl=dg=='Control'
        if case.sum()<3 or ctrl.sum()<3: continue
        Y=np.log2(1+(C.multiply(1e4/lib[:,None])).toarray())
        fc=np.nanmean(Y[case],axis=0)-np.nanmean(Y[ctrl],axis=0)
        _,p=ttest_ind(Y[case],Y[ctrl],axis=0,equal_var=False,nan_policy='omit')
        measured=present[study].reindex(genes).fillna(False).to_numpy(bool) if study in present else np.ones(len(genes),bool)
        detected=((Y>0).mean(axis=0)>=.1)&measured
        p=np.where(detected,p,np.nan); q=bh(p)
        ix=np.flatnonzero(np.isfinite(p))
        all_de.append(pd.DataFrame({'state':state,'state_short':SHORT[state],'study':study,'gene':genes[ix],
            'log2fc':fc[ix],'pval':p[ix],'fdr':q[ix],'n_case_donors':case.sum(),'n_control_donors':ctrl.sum()}))
        print('DE',SHORT[state],study,int(case.sum()),int(ctrl.sum()),len(ix),flush=True)
h.file.close()
de=pd.concat(all_de,ignore_index=True); de.to_csv(TAB/'copd_uncensored_cohort_differentials.tsv.gz',sep='\t',index=False)

# Concordant replication across cohorts, separately for FDR and raw-p criteria.
rep=[]; programs=[]
for (state,gene),x in de.groupby(['state','gene']):
    for criterion in ['fdr','rawp']:
        sig=(x.fdr<=.1) if criterion=='fdr' else (x.pval<.05)
        sig &= x.log2fc.abs()>=FOLD
        up=int((sig&(x.log2fc>0)).sum()); down=int((sig&(x.log2fc<0)).sum())
        direction=1 if up>=2 and down==0 else (-1 if down>=2 and up==0 else 0)
        row={'state':state,'state_short':SHORT[state],'gene':gene,'criterion':criterion,'cohorts_tested':x.study.nunique(),
             'cohorts_up':up,'cohorts_down':down,'direction':direction,'median_log2fc':x.log2fc.median(),
             'min_p':x.pval.min(),'max_p':x.pval.max(),'min_fdr':x.fdr.min(),'max_fdr':x.fdr.max(),
             'studies':';'.join(x.study)}
        rep.append(row)
        if direction: programs.append(row)
rep=pd.DataFrame(rep); prog=pd.DataFrame(programs)
rep.to_csv(TAB/'copd_uncensored_replication_screen.tsv.gz',sep='\t',index=False)
prog.to_csv(TAB/'copd_uncensored_replicated_programs.tsv',sep='\t',index=False)

# Add replicated IPF aberrant-state and acute-injury programs from exact stored rows.
queries=[]; query_info=[]
for (state,criterion),x in prog.groupby(['state','criterion']):
    if len(x)<10: continue
    name=f"COPD_{SHORT[state]}_{criterion}"
    queries.append({'name':name,'genes':x.gene.tolist(),'effects':x.median_log2fc.astype(float).tolist()})
    query_info.append({'query':name,'disease':'COPD','state':state,'criterion':criterion,'n_genes':len(x),'source':'uncensored donor pseudobulk'})
for disease,pop,fn,name in [('IPF','Aberrant basal__vs__AT2','IPF','IPF_aberrant_basal_FDR'),
                            ('infection','AT2-AT1 int.__vs__AT2','infection','INJURY_transition_FDR')]:
    d=pd.read_csv(SRC/f'rows_{fn}.csv.gz',low_memory=False); d=d[(d.modality=='rna')&(d.population==pop)]
    rows=[]
    for gene,x in d.groupby('gene'):
        sig=(x.fdr<=.1)&(x.log2fc.abs()>=FOLD); up=(sig&(x.log2fc>0)).sum(); down=(sig&(x.log2fc<0)).sum()
        if (up>=2 and down==0) or (down>=2 and up==0): rows.append((gene,float(x.log2fc.median())))
    if len(rows)>=10:
        queries.append({'name':name,'genes':[r[0] for r in rows],'effects':[r[1] for r in rows]})
        query_info.append({'query':name,'disease':disease,'state':pop,'criterion':'fdr','n_genes':len(rows),'source':'CellRef2 stored replicated RNA'})
pd.DataFrame(query_info).to_csv(TAB/'drug_query_programs.tsv',sep='\t',index=False)
request={'queries':queries,'gene_space':'all','correlations':['spearman','pearson'],'records_per_tail':1000}
(OUT/'drug_request.json').write_text(json.dumps(request))

driver='/Users/saljh8/Documents/GitHub/LungMAP-discovery/lungmap_discovery/external/site_drug_screen.py'
site='/Users/saljh8/Dropbox/LungMAP/refactored_website/LungMAP-net-refactor/app'
py='/Users/saljh8/Dropbox/LungMAP/refactored_website/LungMAP-net-refactor/.venv/bin/python'
env=dict(os.environ,PYTHONPATH=site,LUNGMAP_DATA_ROOT='/Users/saljh8/Dropbox/LungMAP/refactored_website/lungmap-data')
subprocess.run([py,driver,str(OUT/'drug_request.json'),str(ENG)],cwd=site,env=env,check=True)

# Long score table for high-confidence tails plus compound recurrence.
meta=pd.read_csv(ENG/'signature_meta.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('lmd_key')
tails=[]
for q in queries:
    z=np.load(ENG/f"{q['name']}.npz",allow_pickle=False)
    s=pd.DataFrame({k:z[k] for k in ['rho','pearson','n_shared','n_opposite','n_same','n_evaluable','percentile']},index=z['key'].astype(str))
    s['query']=q['name']; s['opposite_fraction']=s.n_opposite/s.n_evaluable.clip(lower=1)
    s=s.join(meta,how='left'); s=s[(s.n_shared>=50)&(s.percentile>=99.5)&(s.rho<0)]
    tails.append(s.reset_index(names='lmd_key'))
t=pd.concat(tails,ignore_index=True) if tails else pd.DataFrame()
t.to_csv(TAB/'drug_reversal_high_confidence_signatures.tsv.gz',sep='\t',index=False)
if len(t):
    cs=(t.groupby('compound_name').agg(queries_opposed=('query','nunique'),signatures_opposing=('query','size'),
         cell_lines_opposing=('cell_line','nunique'),median_rho=('rho','median'),best_percentile=('percentile','max'),
         queries=('query',lambda x:';'.join(sorted(set(x)))),moa=('moa','first'),targets=('targets','first')).reset_index()
        .sort_values(['queries_opposed','cell_lines_opposing','signatures_opposing','median_rho'],ascending=[False,False,False,True]))
    cs.to_csv(TAB/'drug_reversal_compound_recurrence.tsv',sep='\t',index=False)
    existing=['nintedanib','pirfenidone','roflumilast','budesonide','fluticasone','prednisone','dexamethasone','tiotropium','salmeterol','formoterol','albuterol']
    allrows=[]
    for q in queries:
        z=np.load(ENG/f"{q['name']}.npz",allow_pickle=False); s=pd.DataFrame({k:z[k] for k in ['rho','pearson','n_shared','percentile']},index=z['key'].astype(str)).join(meta)
        s=s[s.compound_name.str.lower().isin(existing)].copy(); s['query']=q['name']; allrows.append(s.reset_index(names='lmd_key'))
    pd.concat(allrows,ignore_index=True).to_csv(TAB/'existing_ipf_copd_drug_state_scores.tsv',sep='\t',index=False)
(OUT/'DONE').write_text('complete\n')
print(json.dumps({'queries':len(queries),'de_rows':len(de),'replicated_program_rows':len(prog),'high_confidence_signatures':len(t)},indent=2),flush=True)
