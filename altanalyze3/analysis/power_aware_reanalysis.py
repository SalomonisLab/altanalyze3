#!/usr/bin/env python3
"""Power-aware discovery/confirmation analysis for unequal cohort sizes.

The largest effective donor cohort is discovery. An asterisk is assigned only
when an independent cohort confirms direction at raw P < 0.05. Weighted
summaries use effective donor N = n_case*n_control/(n_case+n_control), which
prevents a very large case arm from hiding a tiny control arm.
"""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns

ROOT=Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927')
OUT=ROOT/'extended'; TAB=OUT/'tables'; FIG=OUT/'figures'
SRC=Path('/Users/saljh8/Dropbox/LungMAP/Discovery/hypotheses_20260927')

def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); q=np.empty(len(p)); q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    return np.minimum(q,1)

def pass_fold(mod,fc):
    return True if mod in ('grn_tf','fastcomm') else abs(fc)>=np.log2(1.2)

# IPF retains all stored modalities. COPD RNA is replaced by the explicitly
# requested uncensored donor-centroid analysis; non-RNA modalities remain stored.
parts=[]
for disease in ['IPF','COPD']:
    d=pd.read_csv(SRC/f'rows_{disease}.csv.gz')
    if disease=='COPD': d=d[d.modality!='rna'].copy()
    d['disease']=disease; d['cohort']=d.contrast_id.str.split('__').str[0]
    d=d.rename(columns={'n_case':'n_case_donors','n_control':'n_control_donors'})
    parts.append(d[['disease','modality','population','gene','cohort','contrast_id','log2fc','pval','fdr','n_case_donors','n_control_donors']])

u=pd.read_csv(TAB/'copd_uncensored_cohort_differentials.tsv.gz',sep='\t')
extra_path=TAB/'copd_uncensored_additional_state_differentials.tsv.gz'
if extra_path.exists(): u=pd.concat([u,pd.read_csv(extra_path,sep='\t')],ignore_index=True)
u['disease']='COPD'; u['modality']='rna'; u['population']=u.state_short; u['cohort']=u.study; u['contrast_id']=u.study+'__COPD_vs_control'
u=u.rename(columns={'pval':'pval'})
parts.append(u[['disease','modality','population','gene','cohort','contrast_id','log2fc','pval','fdr','n_case_donors','n_control_donors']])
d=pd.concat(parts,ignore_index=True)
d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=['log2fc','pval','fdr','n_case_donors','n_control_donors'])
d['effective_n']=d.n_case_donors*d.n_control_donors/(d.n_case_donors+d.n_control_donors)

keys=['disease','modality','population','gene']
d['fold_ok']=[pass_fold(m,f) for m,f in zip(d.modality,d.log2fc)]
d['signed_z']=np.sign(d.log2fc)*norm.isf(np.maximum(d.pval,1e-300)/2)
d['weighted_fc_num']=d.effective_n*d.log2fc
d['weighted_z_num']=np.sqrt(d.effective_n)*d.signed_z

# Select discovery rows without Python-level iteration over hundreds of thousands of genes.
ds=d.sort_values(keys+['effective_n','n_case_donors','n_control_donors'],ascending=[True]*4+[False,False,False])
disc=ds.drop_duplicates(keys).copy()
disc['discovery_pass']=(disc.fdr<=.10)&disc.fold_ok
disc['discovery_sign']=np.sign(disc.log2fc)
disc=disc.rename(columns={'cohort':'discovery_cohort','n_case_donors':'discovery_n_case',
    'n_control_donors':'discovery_n_control','effective_n':'discovery_effective_n',
    'log2fc':'discovery_log2fc','pval':'discovery_p','fdr':'discovery_fdr'})

other=d.merge(disc[keys+['discovery_cohort','discovery_sign']],on=keys,how='left')
other=other[other.cohort!=other.discovery_cohort]
cmask=(np.sign(other.log2fc)==other.discovery_sign)&(other.pval<.05)&other.fold_ok
xmask=(np.sign(other.log2fc)==-other.discovery_sign)&(other.pval<.05)&other.fold_ok
def summarize_hits(x,prefix):
    if x.empty: return pd.DataFrame(columns=keys+[f'n_{prefix}_cohorts',f'{prefix}_cohorts'])
    return (x.groupby(keys).cohort.agg([('n','nunique'),('names',lambda z:';'.join(sorted(set(map(str,z)))))]).reset_index()
             .rename(columns={'n':f'n_{prefix}_cohorts','names':f'{prefix}_cohorts'}))
conf=summarize_hits(other[cmask],'confirmation'); conflict=summarize_hits(other[xmask],'conflicting')

g=(d.groupby(keys,as_index=False).agg(cohorts_tested=('cohort','nunique'),total_case_donors=('n_case_donors','sum'),
    total_control_donors=('n_control_donors','sum'),sum_effective_n=('effective_n','sum'),
    weighted_fc_num=('weighted_fc_num','sum'),weighted_z_num=('weighted_z_num','sum')))
g['weighted_log2fc']=g.weighted_fc_num/g.sum_effective_n
g['weighted_stouffer_z']=g.weighted_z_num/np.sqrt(g.sum_effective_n)
g['weighted_p']=2*norm.sf(g.weighted_stouffer_z.abs())
a=(disc[keys+['discovery_cohort','discovery_n_case','discovery_n_control','discovery_effective_n','discovery_log2fc',
      'discovery_p','discovery_fdr','discovery_pass']].merge(conf,on=keys,how='left').merge(conflict,on=keys,how='left').merge(g,on=keys,how='left'))
for c in ['n_confirmation_cohorts','n_conflicting_cohorts']: a[c]=a[c].fillna(0).astype(int)
for c in ['confirmation_cohorts','conflicting_cohorts']: a[c]=a[c].fillna('')
a['independently_confirmed']=a.discovery_pass&(a.n_confirmation_cohorts>0)&(a.n_conflicting_cohorts==0)
a['display_label']=a.gene+np.where(a.independently_confirmed,'*','')
a=a.drop(columns=['sum_effective_n','weighted_fc_num','weighted_z_num'])
a['weighted_fdr']=a.groupby(['disease','modality']).weighted_p.transform(bh)
a.to_csv(TAB/'power_aware_associations.tsv.gz',sep='\t',index=False)
find=a[a.discovery_pass].sort_values(['independently_confirmed','discovery_effective_n','discovery_fdr'],ascending=[False,False,True])
find.to_csv(TAB/'power_aware_discovery_findings.tsv',sep='\t',index=False)

# Cohort power audit, retaining state-specific sample sizes.
audit=(d.groupby(['disease','modality','population','cohort'],as_index=False)
         .agg(n_case=('n_case_donors','max'),n_control=('n_control_donors','max'),effective_n=('effective_n','max')))
audit['raw_total_n']=audit.n_case+audit.n_control
audit.to_csv(TAB/'cohort_power_summary.tsv',sep='\t',index=False)

# Exact RNA associations discovered in both diseases; stars are disease-specific confirmations.
r=find[(find.modality=='rna')].copy()
cross=r.groupby(['population','gene']).filter(lambda x:x.disease.nunique()==2)
cross.to_csv(TAB/'power_aware_cross_disease_findings.tsv',sep='\t',index=False)

# Add power-aware columns to the human-readable evidence workbook without deleting prior detail.
ev=pd.read_excel(TAB/'intuitive_evidence_matrix.xlsx')
mapping={
 'IPF FN1–αv–TEAD2/ZNF322 epithelial circuit':('IPF','Aberrant basal__vs__AT2','TEAD2'),
 'TEAD2 adhesion/basement-membrane regulon':('IPF','Aberrant basal__vs__AT2','ITGA3'),
 'COPD AT2 basal/secretory priming':('COPD','AT2','TP63'),
 'COPD pulmonary venous endothelial identity loss':('COPD','PVEC','FBLIM1'),
 'Shared IPF–COPD CCL18 macrophage program':('COPD','IM','CCL18'),
 'Shared IPF–COPD stromal/endothelial remodeling':('COPD','PBFB','PTGDS')}
def lookup(name):
    if name not in mapping: return pd.Series({'power_aware_status':'not reclassified','discovery_cohort':'—','discovery_donors':'—','replication_marker':''})
    disease,pop,gene=mapping[name]; q=a[(a.disease==disease)&(a.population==pop)&(a.gene==gene)]
    if q.empty: return pd.Series({'power_aware_status':'association not present in tested table','discovery_cohort':'—','discovery_donors':'—','replication_marker':''})
    q=q.sort_values(['discovery_pass','independently_confirmed','discovery_effective_n'],ascending=False).iloc[0]
    status='discovery + independent confirmation' if q.independently_confirmed else ('discovery only' if q.discovery_pass else 'did not pass largest-cohort FDR discovery')
    return pd.Series({'power_aware_status':status,'discovery_cohort':q.discovery_cohort,
      'discovery_donors':f"{q.discovery_n_case} case / {q.discovery_n_control} control",'replication_marker':'*' if q.independently_confirmed else ''})
extra=ev.finding.apply(lookup); ev=pd.concat([ev,extra],axis=1)
ev['finding_display']=ev.finding+ev.replication_marker
with pd.ExcelWriter(TAB/'intuitive_evidence_matrix_power_aware.xlsx') as w:
    ev.to_excel(w,index=False,sheet_name='Evidence')
    audit.to_excel(w,index=False,sheet_name='Cohort power')
    find.to_excel(w,index=False,sheet_name='Discovery findings')
ev.to_csv(TAB/'intuitive_evidence_matrix_power_aware.tsv',sep='\t',index=False)

# Figure 9A documents unequal power; 9B shows top independently confirmed discoveries.
sns.set_theme(style='whitegrid',font_scale=.8)
fig,ax=plt.subplots(1,2,figsize=(15,7),gridspec_kw={'width_ratios':[.95,1.4]})
aud=audit.groupby(['disease','cohort'],as_index=False).agg(raw_total_n=('raw_total_n','max'),effective_n=('effective_n','max'))
aud['label']=aud.disease+' | '+aud.cohort
aud=aud.sort_values('raw_total_n')
ax[0].barh(aud.label,aud.raw_total_n,color='#b9c6d8',label='raw total N')
ax[0].barh(aud.label,aud.effective_n,color='#264653',label='effective N')
ax[0].set_title('A  Cohort power is unequal',loc='left',weight='bold'); ax[0].set_xlabel('Donor count / effective donor N'); ax[0].legend(frameon=False)

top=find[find.independently_confirmed & find.modality.eq('rna') & ~find.population.str.contains('__vs__',regex=False)].copy()
top['absz']=top.weighted_stouffer_z.abs(); top=top.nlargest(18,'absz').sort_values('weighted_log2fc')
top['label']=top.disease+' | '+top.population+' | '+top.gene+'*'
ax[1].barh(top.label,top.weighted_log2fc,color=np.where(top.weighted_log2fc>0,'#d95f59','#457b9d'))
ax[1].axvline(0,color='black',lw=.8); ax[1].set_xlabel('effective-N weighted log2FC')
ax[1].set_title('B  Largest-cohort discoveries with independent confirmation*',loc='left',weight='bold')
fig.suptitle('Power-aware disease associations',fontsize=15,weight='bold')
fig.text(.5,.012,'* independent cohort: concordant direction, raw P<0.05, no significant directional conflict. Discovery requires largest effective-N cohort FDR≤0.10.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.05,1,.95])
fig.savefig(FIG/'Figure9_power_aware_discovery_confirmation.pdf',bbox_inches='tight')
fig.savefig(FIG/'Figure9_power_aware_discovery_confirmation.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# Replace the earlier equal-cohort cross-disease panel with associations that
# pass discovery and independent confirmation in both diseases.
both=cross.groupby(['population','gene']).filter(lambda z:z.disease.nunique()==2 and z.independently_confirmed.all()).copy()
chosen=[('AdvFB','RGCC'),('AdvFB','CXCL12'),('AdvFB','CFI'),('IM','SGK1'),('IM','USP53'),
        ('AT1','PLPP3'),('AM','CD163'),('AM','PLIN2')]
bb=both.set_index(['population','gene','disease']).discovery_log2fc
mat=pd.DataFrame({dis:[bb.get((pop,gene,dis),np.nan) for pop,gene in chosen] for dis in ['IPF','COPD']},
                 index=[f'{pop} | {gene}*' for pop,gene in chosen])
fig,ax=plt.subplots(figsize=(6.5,5.8)); sns.heatmap(mat,cmap='vlag',center=0,annot=True,fmt='.2f',linewidths=.7,
    cbar_kws={'label':'largest-cohort discovery log2FC'},ax=ax)
ax.set_title('Power-aware biology reproduced across IPF and COPD',weight='bold'); ax.set_xlabel(''); ax.set_ylabel('cell state | gene')
fig.text(.5,.01,'* independently confirmed in both diseases',ha='center',fontsize=9); fig.tight_layout(rect=[0,.04,1,1])
fig.savefig(FIG/'Figure7_shared_IPF_COPD_biology.pdf',bbox_inches='tight')
fig.savefig(FIG/'Figure7_shared_IPF_COPD_biology.png',dpi=300,bbox_inches='tight'); plt.close(fig)

shutil.copy2(Path(__file__),OUT/'scripts'/Path(__file__).name)
print({'associations':len(a),'discovery_findings':len(find),'independently_confirmed':int(find.independently_confirmed.sum()),
       'cross_disease_rows':len(cross)})
