#!/usr/bin/env python3
"""Run drug signatures using only power-aware discovery+confirmation programs."""
from pathlib import Path
import json,os,subprocess
import numpy as np
import pandas as pd
ROOT=Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927/extended'); TAB=ROOT/'tables'; ENG=ROOT/'drug_engine_power_aware'; ENG.mkdir(exist_ok=True)
a=pd.read_csv(TAB/'power_aware_discovery_findings.tsv',sep='\t',low_memory=False)
x=a[(a.modality=='rna')&a.independently_confirmed].copy(); queries=[]; info=[]
for disease,pop in [('IPF','Aberrant basal__vs__AT2')]+[('COPD',p) for p in sorted(x.loc[x.disease.eq('COPD'),'population'].unique())]:
    z=x[(x.disease==disease)&(x.population==pop)]
    if len(z)<10: continue
    name=('IPF_aberrant_basal' if disease=='IPF' else 'COPD_'+pop.replace(' ','_'))+'_power_confirmed'
    queries.append({'name':name,'genes':z.gene.tolist(),'effects':z.weighted_log2fc.astype(float).tolist()})
    info.append({'query':name,'disease':disease,'state':pop,'criterion':'largest effective-N cohort FDR<=0.10 + independent raw-p confirmation','n_genes':len(z),'source':'uncensored COPD centroids / stored IPF donor pseudobulk'})
pd.DataFrame(info).to_csv(TAB/'power_aware_drug_query_programs.tsv',sep='\t',index=False)
req={'queries':queries,'gene_space':'all','correlations':['spearman','pearson'],'records_per_tail':1000}; (ROOT/'power_aware_drug_request.json').write_text(json.dumps(req))
driver='/Users/saljh8/Documents/GitHub/LungMAP-discovery/lungmap_discovery/external/site_drug_screen.py'; site='/Users/saljh8/Dropbox/LungMAP/refactored_website/LungMAP-net-refactor/app'; py='/Users/saljh8/Dropbox/LungMAP/refactored_website/LungMAP-net-refactor/.venv/bin/python'
env=dict(os.environ,PYTHONPATH=site,LUNGMAP_DATA_ROOT='/Users/saljh8/Dropbox/LungMAP/refactored_website/lungmap-data'); subprocess.run([py,driver,str(ROOT/'power_aware_drug_request.json'),str(ENG)],cwd=site,env=env,check=True)
meta=pd.read_csv(ENG/'signature_meta.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('lmd_key'); tails=[]
for q in queries:
    z=np.load(ENG/f"{q['name']}.npz",allow_pickle=False); s=pd.DataFrame({k:z[k] for k in ['rho','pearson','n_shared','n_opposite','n_same','n_evaluable','percentile']},index=z['key'].astype(str)); s['query']=q['name']; s=s.join(meta,how='left'); tails.append(s[(s.n_shared>=50)&(s.percentile>=99.5)&(s.rho<0)].reset_index(names='lmd_key'))
t=pd.concat(tails,ignore_index=True); t.to_csv(TAB/'power_aware_drug_reversal_signatures.tsv.gz',sep='\t',index=False)
cs=(t.groupby('compound_name').agg(queries_opposed=('query','nunique'),signatures_opposing=('query','size'),cell_lines_opposing=('cell_line','nunique'),median_rho=('rho','median'),best_percentile=('percentile','max'),queries=('query',lambda z:';'.join(sorted(set(z)))),moa=('moa','first'),targets=('targets','first')).reset_index().sort_values(['queries_opposed','cell_lines_opposing','signatures_opposing'],ascending=False))
cs.to_csv(TAB/'power_aware_drug_compound_recurrence.tsv',sep='\t',index=False)
priority={'Y-39983':('ROCK inhibitor','1'),'fostamatinib':('SYK inhibitor','2'),'PD-0325901':('MEK1/2 inhibitor','3'),
          'selumetinib':('MEK1/2 inhibitor','4'),'Y-27632':('ROCK inhibitor','5')}
pick=cs[cs.compound_name.isin(priority)].copy(); pick['mechanism']=pick.compound_name.map(lambda z:priority[z][0]); pick['priority']=pick.compound_name.map(lambda z:priority[z][1]); pick=pick.sort_values('priority')
pick['power_aware_interpretation']='Reverses independently confirmed IPF aberrant-epithelial and COPD interstitial-macrophage programs'
pick.to_csv(TAB/'ranked_drug_predictions_power_aware.tsv',sep='\t',index=False)
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style='whitegrid'); q=pick.sort_values('cell_lines_opposing')
fig,ax=plt.subplots(figsize=(8,4.8)); ax.barh(q.compound_name+'*',q.cell_lines_opposing,color=q.mechanism.map({'ROCK inhibitor':'#2a9d8f','SYK inhibitor':'#457b9d','MEK1/2 inhibitor':'#e9c46a'}))
for y,(_,r) in enumerate(q.iterrows()): ax.text(r.cell_lines_opposing+.15,y,f"{r.signatures_opposing} signatures",va='center',fontsize=9)
ax.set_xlabel('Independent perturbation cell lines with top-0.5% reversal'); ax.set_title('Power-aware drug reversals reproduced across disease programs',weight='bold')
fig.text(.5,.01,'* reverses both independently confirmed IPF aberrant-epithelial and COPD interstitial-macrophage programs',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.05,1,1]); fig.savefig(ROOT/'figures/Figure10_power_aware_drug_reversals.pdf',bbox_inches='tight'); fig.savefig(ROOT/'figures/Figure10_power_aware_drug_reversals.png',dpi=300,bbox_inches='tight'); plt.close(fig)
print({'queries':len(queries),'high_confidence_signatures':len(t),'compounds':len(cs)})
