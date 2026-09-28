#!/usr/bin/env python3
"""Create intuitive evidence, lipid, isolation, recurrence, and drug tables/figures."""
from pathlib import Path
import json, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ROOT=Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927')
OUT=ROOT/'extended'; TAB=OUT/'tables'; FIG=OUT/'figures'; FIG.mkdir(parents=True,exist_ok=True)
SRC=Path('/Users/saljh8/Dropbox/LungMAP/Discovery/hypotheses_20260927')

evidence=pd.DataFrame([
{'finding':'IPF FN1–αv–TEAD2/ZNF322 epithelial circuit','priority':'Lead','disease':'IPF','cell_state':'KRT5−/KRT17+ aberrant epithelium','cohort_replication':'2 CellRef2 cohorts','measured_RNA':'TEAD2, ITGA3, ITGB6, MMP7 ↑; ABCA3 ↓','GRN_TF_prediction':'TEAD2/ZNF322 ↑; CREB3L1 ↓','communication_prediction':'AM FN1→αvβ6/β8/β1; MT+AM SPP1→αv ↑','ADT_prediction':'CD45− CD326+ CD66c-high PDPN+ CD55+ CD49a+','lipid_prediction':'—','composition':'CellRef2 12×; Xenium donor g=3.04, FDR=5.4e−5','spatial':'KRT↔FB z −0.18→10.32/6.33; KRT↔SPP1-Mφ −0.24→2.61/6.07','orthogonal_validation':'Xenium: SPP1-Mφ FN1 g=1.91; FB FN1 g=1.77; KRT ITGB6 g=1.70','literature_boundary':'Broad SPP1-Mφ/CTHRC1-FB/basaloid niche known; exact circuit nominated','decisive_test':'FN1 matrix/SPP1-Mφ ± αvβ6 block × TEAD2/ZNF322 CRISPRi','main_caveat':'Aberrant-basal-vs-AT2 is identity-mixed; GRN/communication/ADT are RNA-derived','figure':'Figure1'},
{'finding':'TEAD2 adhesion/basement-membrane regulon','priority':'Lead mechanism','disease':'IPF','cell_state':'Aberrant epithelium','cohort_replication':'35 TEAD2 and 18 ZNF322 targets replicated-up','measured_RNA':'ITGB4, ITGA3, CDH3, JUP, HSPG2, LAMA5, PTK7, CTNND1, CDC42 ↑','GRN_TF_prediction':'Direct predicted targets','communication_prediction':'Compatible with FN1/integrin input','ADT_prediction':'CD49a/ITGA1 and epithelial surface phenotype','lipid_prediction':'—','composition':'—','spatial':'Same KRT niche','orthogonal_validation':'Target-level RNA coherence; Xenium ITGB6','literature_boundary':'Specific regulon connection not found in targeted search','decisive_test':'Epistasis and CUT&RUN/ATAC after FN1/αv perturbation','main_caveat':'Target coherence is partly circular with TF activity score','figure':'Figure1 + target table'},
{'finding':'Acute-injury ETV5/XBP1 transitional checkpoint','priority':'Second','disease':'COVID-19 + pneumonia','cell_state':'AT2→AT1 intermediate','cohort_replication':'2 injury contrasts','measured_RNA':'SFTPA1/2, SFTPC, MFSD2A ↓ vs AT2','GRN_TF_prediction':'ETV5, XBP1, CEBPA, MLXIPL ↓','communication_prediction':'HMGB1→AGER ↓; CADM1 self-input ↑','ADT_prediction':'PDPN and CD66c ↑; CD49f ↓ vs AT2','lipid_prediction':'PC(20:4/22:6), PC(20:4/20:4), PE(P-16:0/20:4) ↓ vs AT1','composition':'—','spatial':'—','orthogonal_validation':'Same direction in COVID-19 and non-COVID pneumonia','literature_boundary':'Mechanistic coupling of ETV5/XBP1 to exact lipids is prediction','decisive_test':'ETV5 CRISPRa/XBP1s in infected alveolospheres + targeted lipidomics/surfactant biophysics','main_caveat':'Two different anchors; lipid/ADT are RNA-imputed','figure':'Figure2'},
{'finding':'COPD AT2 basal/secretory priming','priority':'Exploratory','disease':'COPD','cell_state':'AT2','cohort_replication':'3 donor cohorts','measured_RNA':'Uncensored-centroid replication table generated separately','GRN_TF_prediction':'TP63 pooled g=1.61; FOXQ1 g=1.26','communication_prediction':'—','ADT_prediction':'No robust replicated isolation marker','lipid_prediction':'—','composition':'—','spatial':'—','orthogonal_validation':'Direction concordant across 3 cohorts','literature_boundary':'Driver status unresolved','decisive_test':'TP63/FOXQ1 CRISPRi in smoke-exposed distal organoids','main_caveat':'Regulons learned in other epithelial states','figure':'Figure3'},
{'finding':'COPD pulmonary venous endothelial identity loss','priority':'Exploratory','disease':'COPD','cell_state':'PVEC','cohort_replication':'2 stored cohorts + uncensored scan','measured_RNA':'FBLIM1, CHD6, DACH1 ↓','GRN_TF_prediction':'FOXA2, SREBF2, KLF5, ELF3, NFIC/NFIX, TEAD2/3 ↓','communication_prediction':'—','ADT_prediction':'No robust replicated marker','lipid_prediction':'—','composition':'—','spatial':'—','orthogonal_validation':'Independent cohort direction','literature_boundary':'Potential vascular mechanism','decisive_test':'Endothelial chip: restore FOXA2/SREBF2 and test barrier/lipid handling','main_caveat':'Needs direct protein/function validation','figure':'Figure3'},
{'finding':'Shared IPF–COPD CCL18 macrophage program','priority':'Cross-disease','disease':'IPF + COPD','cell_state':'Interstitial macrophage','cohort_replication':'2 cohorts within each disease','measured_RNA':'CCL18 ↑ (mean log2FC 1.79 IPF; 1.41 COPD), HLA-DRA ↑','GRN_TF_prediction':'—','communication_prediction':'Basaloid APP→CD74 / ALCAM circuits recur','ADT_prediction':'—','lipid_prediction':'—','composition':'—','spatial':'—','orthogonal_validation':'Four disease/cohort combinations','literature_boundary':'Suggests shared fibrotic/inflammatory macrophage axis','decisive_test':'CCL18 blockade across IPF/COPD fibroblast cocultures','main_caveat':'Shared association does not imply identical upstream cause','figure':'Figure7'},
{'finding':'Shared IPF–COPD stromal/endothelial remodeling','priority':'Cross-disease','disease':'IPF + COPD','cell_state':'PBFB / AdvFB','cohort_replication':'2 cohorts within each disease','measured_RNA':'PBFB PTGDS ↑; AdvFB RGCC and HLA-B ↑','GRN_TF_prediction':'—','communication_prediction':'SPFB LAMA2→CD44 ↑ in IM','ADT_prediction':'SPFB CD90-high CD49a+ replicated in IPF','lipid_prediction':'—','composition':'Fibroblast states enriched','spatial':'Fibrotic FB proximity to KRT state in IPF','orthogonal_validation':'Same genes recur independently in both diseases','literature_boundary':'Potential shared repair-to-fibrosis module','decisive_test':'Perturb PTGDS/LAMA2-CD44 in disease-specific stromal cocultures','main_caveat':'Cell-state abundance can influence apparent programs','figure':'Figure7'},
])
evidence.to_csv(TAB/'intuitive_evidence_matrix.tsv',sep='\t',index=False)
evidence.to_excel(TAB/'intuitive_evidence_matrix.xlsx',index=False)

# Exact database identifiers plus biochemical expansion and database pathway classes.
wp=json.loads(Path('/Users/saljh8/Dropbox/LungMAP/refactored_website/lungmap-data/cellref2/pathways/wikipathways_lipids.json').read_text())
lipids=[
('PC(20:4/22:6)','Phosphatidylcholine 20:4/22:6','arachidonoyl + docosahexaenoyl phosphatidylcholine','PC','42:10'),
('PC(20:4/20:4)','Phosphatidylcholine 20:4/20:4','di-arachidonoyl phosphatidylcholine','PC','40:8'),
('PE(P-16:0/20:4)','Plasmenyl-phosphatidylethanolamine P-16:0/20:4','16:0 vinyl-ether plasmalogen PE + arachidonoyl','PE','36:4')]
lrows=[]
for dbid,full,chains,cls,total in lipids:
    paths=[p['name'] for p in wp['pathways'] if dbid in p.get('species',[])]
    lrows.append({'lungmap_lipid_id':dbid,'expanded_name':full,'chain_interpretation':chains,'lipid_class':cls,
                  'total_composition':total,'direction':'lower in AT2–AT1 intermediate relative to AT1',
                  'replication':'COVID-19 and pneumonia','evidence_type':'RNA-imputed lipid prediction',
                  'mapped_LungMAP_WikiPathways':'; '.join(paths),'validation':'targeted LC–MS/MS with isomer-resolved standards'})
pd.DataFrame(lrows).to_csv(TAB/'lipid_name_and_pathway_mapping.tsv',sep='\t',index=False)

adt=pd.DataFrame([
{'target_population':'IPF KRT5−/KRT17+ aberrant epithelium','enrichment_gate':'Live singlet CD45− CD31− EPCAM/CD326+; enrich CD66c-high PDPN+ CD55+ CD49a+','depletion_gate':'Exclude CD45+, CD31+, CD206+, CD169+','novel_strategy':'EPCAM+ PDPN+ CD55-high CD49a+ combination separates injured aberrant epithelial cells from immune/endothelial cells; retain CD66c as intensity axis','replication':'ADT differential only Adams; RNA ITGA3/ITGB6 and Xenium support receiver phenotype','status':'Experimental proposal—ADT is RNA-imputed'},
{'target_population':'IPF subpleural fibroblast','enrichment_gate':'Live singlet CD45− EPCAM/CD326− CD31− CD90-high CD49a+ CD140a+','depletion_gate':'Exclude CD45+, EPCAM+, HLA-DR+','novel_strategy':'CD90-high/CD49a double-positive stromal gate; add PDGFRα/CD140a for specificity','replication':'CD90 and CD49a increase in Adams and Natri predictions','status':'Predicted; protein validation required'},
{'target_population':'IPF Langerhans','enrichment_gate':'Live singlet CD45+ HLA-DR-high CD1c+ CD16-low; confirm CD207/langerin and CD1A intracellular/transcript','depletion_gate':'Exclude CD14+, CD16-high, CD66c+','novel_strategy':'CD1c+ HLA-DR-high CD16-low gate followed by CD207 confirmation','replication':'CD45/CD1c direction in Adams and Natri; HLA-DR strongest in Natri','status':'Near-term validation panel'},
{'target_population':'SPP1 macrophage sender','enrichment_gate':'Live singlet CD45+ CD68/CD64+; index-sort SPP1-high or osteopontin secretion; pair with FN1 protein staining','depletion_gate':'Exclude EPCAM+, CD3+, CD19+','novel_strategy':'Functional sender isolation should combine macrophage identity with secreted SPP1/FN1 rather than rely on imputed ADT alone','replication':'Xenium directly supports disease FN1 in SPP1 macrophages','status':'Xenium-supported, flow panel needs direct optimization'}])
adt.to_csv(TAB/'adt_isolation_strategies.tsv',sep='\t',index=False)

# Exact reproduced IPF/COPD biology.
parts=[]
for cond in ['IPF','COPD']:
    d=pd.read_csv(SRC/f'screen_{cond}.csv.gz'); d=d[d.replicated].copy(); d['condition']=cond; parts.append(d)
a=pd.concat(parts); shared=a.groupby(['modality','population','gene']).filter(lambda x:x.condition.nunique()==2)
shared.to_csv(TAB/'shared_ipf_copd_reproduced_features.tsv',sep='\t',index=False)

# Evidence coverage heatmap: 2 independent/measured, 1 model-derived/supporting, 0 absent.
cols=['Replication','Measured RNA','GRN/TF','Communication','ADT','Lipid','Composition','Spatial','Orthogonal cohort']
scores=np.array([[2,2,1,1,1,0,2,2,2],[2,2,1,1,1,0,0,2,2],[2,2,1,1,1,1,0,0,1],
                 [2,0,1,0,0,0,0,0,1],[2,2,1,0,0,0,0,0,1],[2,2,0,1,0,0,0,0,2],[2,2,0,1,1,0,1,1,2]])
fig,ax=plt.subplots(figsize=(14,6)); sns.heatmap(scores,annot=np.where(scores==2,'M',np.where(scores==1,'P','')),
 fmt='',cmap=sns.color_palette(['#f3f4f6','#f4a261','#2a9d8f'],as_cmap=True),vmin=0,vmax=2,cbar=False,
 xticklabels=cols,yticklabels=evidence.finding,linewidths=.8,linecolor='white',ax=ax)
ax.set_title('Evidence architecture  |  M = measured/independently replicated; P = RNA-derived prediction/support',weight='bold'); ax.set_xlabel(''); ax.set_ylabel(''); ax.tick_params(axis='x',rotation=30); fig.tight_layout()
fig.savefig(FIG/'Figure5_intuitive_evidence_matrix.pdf',bbox_inches='tight'); fig.savefig(FIG/'Figure5_intuitive_evidence_matrix.png',dpi=300,bbox_inches='tight'); plt.close(fig)

# Cross-disease heatmap of exact mean effects.
sel=[('rna','IM','CCL18'),('rna','IM','HLA-DRA'),('rna','PBFB','PTGDS'),('rna','AdvFB','RGCC'),('rna','AdvFB','HLA-B'),
     ('rna','AT1','CD24'),('rna','Ciliated-Bronch','SLC34A2')]
s=shared.set_index(['modality','population','gene','condition']).mean_log2fc
mat=pd.DataFrame({c:[s.get((*k,c),np.nan) for k in sel] for c in ['IPF','COPD']},index=[f'{p} | {g}' for _,p,g in sel])
fig,ax=plt.subplots(figsize=(6,5)); sns.heatmap(mat,cmap='Reds',annot=True,fmt='.2f',linewidths=.7,cbar_kws={'label':'mean replicated log2FC'},ax=ax)
ax.set_title('Surprising reproduced biology across IPF and COPD',weight='bold'); ax.set_xlabel(''); ax.set_ylabel('cell state | gene'); fig.tight_layout()
fig.savefig(FIG/'Figure7_shared_IPF_COPD_biology.pdf',bbox_inches='tight'); fig.savefig(FIG/'Figure7_shared_IPF_COPD_biology.png',dpi=300,bbox_inches='tight'); plt.close(fig)

print({'evidence_rows':len(evidence),'shared_rows':len(shared),'lipids':len(lrows),'adt_strategies':len(adt)})
