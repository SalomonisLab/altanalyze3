#!/usr/bin/env python3
"""Rank perturbational reversals and make an interpretable drug/off-target figure."""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927')
OUT = ROOT/'extended'; TAB = OUT/'tables'; FIG = OUT/'figures'
FIG.mkdir(parents=True, exist_ok=True)

rec = pd.read_csv(TAB/'drug_reversal_compound_recurrence.tsv', sep='\t').fillna('')
by_name = rec.set_index(rec.compound_name.str.lower())

# Selection is intentionally mechanistic: recurrent across both FDR-supported IPF/repair
# programs, supported by an extra-pulmonary experiment, and not merely a generic poison.
manual = [
 ('Y-27632','ROCK1/2 inhibitor','1','Direct FN1–integrin–actomyosin/YAP-axis intervention; strongest median reversal among prioritized compounds',
  'Kidney UUO: reduced αSMA, macrophage infiltration, TGFβ, collagen I and interstitial fibrosis','https://pubmed.ncbi.nlm.nih.gov/11967018/','Preclinical; requires lung-selective delivery and epithelial repair testing'),
 ('GSK-429286A','ROCK1/2 inhibitor','2','Independent compound-level confirmation of the ROCK class across both stringent programs',
  'Class evidence from Y-27632 in kidney fibrosis supports ROCK mechanism, not this exact molecule','https://pubmed.ncbi.nlm.nih.gov/11967018/','Tool compound; class-level rather than exact-drug validation'),
 ('fostamatinib','SYK inhibitor','3','Reverses both stringent programs and connects macrophage inflammation to TGFβ/Smad fibrosis',
  'Rat peritoneal fibrosis: SYK inhibition reduced inflammatory/fibrotic signaling through TGFβ1/Smad3','https://pmc.ncbi.nlm.nih.gov/articles/PMC6918804/','Clinically used systemic drug; monitor hypertension, liver enzymes and infection risk'),
 ('selumetinib','MEK1/2 inhibitor','4','Highly recurrent reversal; targets an epithelial–stromal EGFR/MEK repair signal',
  'Bleomycin skin fibrosis: <10% skin-thickness increase with selumetinib versus >80% in controls, with less collagen staining','https://link.springer.com/article/10.1186/s13578-021-00553-0','Antiproliferative mechanism may impair normal epithelial repair; schedule/dose critical'),
 ('PD-0325901','MEK1/2 inhibitor','5','Independent MEK-class hit across both stringent programs and 18 source cell lines',
  'Post-MI mouse heart: enhanced peri-infarct vascularization, smaller scar and improved systolic function','https://www.jci.org/articles/view/152308','Potential reparative vascular benefit, but class ocular/cardiac toxicity requires caution'),
 ('verteporfin','YAP–TEAD perturbagen','Mechanistic add-on','Directly matches the nominated TEAD2 mechanotransduction circuit despite lower signature recurrence',
  'Kidney studies show reduced myofibroblast/ECM programs, but systemic verteporfin worsened UUO fibrosis in an endothelial-context study','https://pubmed.ncbi.nlm.nih.gov/34151201/','Context-dependent warning: pursue cell-targeted YAP/TEAD inhibition, not systemic verteporfin'),
]
rows=[]
for name,mech,rank,why,orth,url,caveat in manual:
    x=by_name.loc[name.lower()]
    if isinstance(x,pd.DataFrame): x=x.iloc[0]
    q=set(str(x['queries']).split(';'))
    rows.append(dict(priority=rank,compound=name,mechanism=mech,
      both_FDR_supported_queries={'IPF_aberrant_basal_FDR','INJURY_transition_FDR'}.issubset(q),
      FDR_queries_opposed='; '.join(sorted(q & {'IPF_aberrant_basal_FDR','INJURY_transition_FDR'})),
      COPD_rawp_queries_opposed=len([z for z in q if z.startswith('COPD_')]),
      total_queries_opposed=int(x.queries_opposed),signatures_opposing=int(x.signatures_opposing),
      source_cell_lines_opposing=int(x.cell_lines_opposing),median_rho=float(x.median_rho),
      best_percentile=float(x.best_percentile),selection_reason=why,
      orthogonal_extra_pulmonary_evidence=orth,literature_url=url,safety_or_interpretation_caveat=caveat))
ranked=pd.DataFrame(rows)
ranked.to_csv(TAB/'ranked_drug_predictions.tsv',sep='\t',index=False)
ranked[['compound','mechanism','orthogonal_extra_pulmonary_evidence','literature_url','safety_or_interpretation_caveat']].to_csv(
    TAB/'drug_literature_evidence.tsv',sep='\t',index=False)

# Existing-drug state summary. Positive rho means disease-state mimicry, not a clinical AE.
ex=pd.read_csv(TAB/'existing_ipf_copd_drug_state_scores.tsv',sep='\t')
ex=ex[ex.n_shared>=50].copy(); ex['compound']=ex.compound_name.str.lower()
agg=(ex.groupby(['compound','query']).agg(n_signatures=('rho','size'),median_rho=('rho','median'),
     median_reversal_percentile=('percentile','median'),fraction_disease_mimicking=('rho',lambda x:(x>0).mean()),
     fraction_disease_reversing=('rho',lambda x:(x<0).mean()),min_rho=('rho','min'),max_rho=('rho','max')).reset_index())
def call(r):
    if r.median_rho>=.05 and r.fraction_disease_mimicking>=.60: return 'transcriptomic risk flag'
    if r.median_rho<=-.05 and r.fraction_disease_reversing>=.60: return 'reversal tendency'
    return 'mixed/context-dependent'
agg['interpretation']=agg.apply(call,axis=1)
agg['evidence_boundary']=np.where(agg['query'].str.startswith('COPD_'),
    'exploratory: COPD raw-p replicated program','FDR-supported disease/repair program')
agg['not_a_clinical_adverse_event']='Signature mimicry is a hypothesis-generating cell-state flag, not evidence of toxicity or clinical harm.'
agg.to_csv(TAB/'existing_drug_predicted_offtarget_states.tsv',sep='\t',index=False)

# Figure 8: recurrence strength plus existing-therapy state effects.
sns.set_theme(style='whitegrid',font_scale=.85)
fig=plt.figure(figsize=(15,9)); gs=fig.add_gridspec(1,2,width_ratios=[.9,1.65],wspace=.35)
ax=fig.add_subplot(gs[0,0])
rr=ranked[ranked.compound!='verteporfin'].sort_values('source_cell_lines_opposing')
colors=['#2a9d8f' if 'ROCK' in m else '#e9c46a' if 'MEK' in m else '#457b9d' for m in rr.mechanism]
ax.barh(rr.compound,rr.source_cell_lines_opposing,color=colors)
for y,(_,r) in enumerate(rr.iterrows()): ax.text(r.source_cell_lines_opposing+.5,y,f"{r.total_queries_opposed} programs",va='center',fontsize=8)
ax.set_xlabel('Independent source cell lines with top-0.5% reversal')
ax.set_title('A  Recurrent, mechanism-prioritized reversals',loc='left',weight='bold')
ax.spines[['top','right']].set_visible(False)

ax=fig.add_subplot(gs[0,1])
order=['budesonide','fluticasone','dexamethasone','prednisone','formoterol','salmeterol','roflumilast','nintedanib','pirfenidone']
qorder=['IPF_aberrant_basal_FDR','INJURY_transition_FDR','COPD_AdvFB_rawp','COPD_Basal_rawp','COPD_AT1_rawp','COPD_AT2_rawp','COPD_AM_rawp','COPD_AM_lipid_rawp','COPD_tMDM_rawp','COPD_CAP1_rawp','COPD_PVEC_rawp']
mat=agg.pivot(index='compound',columns='query',values='median_rho').reindex(index=order,columns=qorder)
labels=['IPF aberrant epithelium*','Injury transition*','COPD AdvFB','COPD basal','COPD AT1','COPD AT2','COPD AM','COPD lipid-AM','COPD tMDM','COPD CAP1','COPD PVEC']
sns.heatmap(mat,cmap='vlag',center=0,vmin=-.16,vmax=.16,linewidths=.4,xticklabels=labels,
            cbar_kws={'label':'median Spearman ρ\nblue=reversal; red=disease mimicry'},ax=ax)
ax.set_title('B  Existing-drug transcriptional effects are state-specific',loc='left',weight='bold')
ax.set_xlabel(''); ax.set_ylabel(''); ax.tick_params(axis='x',rotation=55); ax.tick_params(axis='y',rotation=0)
fig.suptitle('Drug-signature evidence: recurrent candidates and predicted state liabilities',weight='bold',fontsize=15,y=.995)
fig.subplots_adjust(bottom=.22,top=.89)
fig.text(.51,.025,'* FDR-supported query. COPD columns use concordant raw-p replication and are exploratory. Red cells are risk flags, not clinical adverse events.',ha='center',fontsize=9)
fig.savefig(FIG/'Figure8_drug_reversal_and_state_liabilities.pdf',bbox_inches='tight')
fig.savefig(FIG/'Figure8_drug_reversal_and_state_liabilities.png',dpi=300,bbox_inches='tight')
plt.close(fig)

# Preserve executable provenance with the result bundle.
prov=OUT/'scripts'; prov.mkdir(exist_ok=True)
for name in ['extended_evidence_worker.py','run_extended_evidence.sh','build_extended_summary.py','finalize_extended_drugs.py']:
    shutil.copy2(Path(__file__).parent/name,prov/name)

print({'ranked_drugs':len(ranked),'state_summaries':len(agg),'risk_flags':int((agg.interpretation=='transcriptomic risk flag').sum())})
