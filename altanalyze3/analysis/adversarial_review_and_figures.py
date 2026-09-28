#!/usr/bin/env python3
"""Build a complete adversarial evidence inventory and replace misleading figures."""
from pathlib import Path
import shutil, textwrap
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns

ROOT=Path('/Users/saljh8/Dropbox/LungMAP/Discovery/codex_nature_predictions_20260927')
EXT=ROOT/'extended'; TAB=EXT/'tables'; FIG=EXT/'figures'; BASEFIG=ROOT/'figures'; LEG=EXT/'legacy_figures'
LEG.mkdir(exist_ok=True); sns.set_theme(style='whitegrid',font_scale=.84)
for p in list(BASEFIG.glob('Figure*.png'))+list(FIG.glob('Figure*.png')):
    q=LEG/p.name
    if not q.exists(): shutil.copy2(p,q)

pa=pd.read_csv(TAB/'power_aware_associations.tsv.gz',sep='\t',low_memory=False)
disc=pd.read_csv(TAB/'power_aware_discovery_findings.tsv',sep='\t',low_memory=False)
cross=pd.read_csv(TAB/'power_aware_cross_disease_findings.tsv',sep='\t',low_memory=False)

claims=[]
def add(cid,domain,claim,disease,state,status,score,discovery,replication,orthogonal,directness,major_risk,decision,test,star=False):
    claims.append(dict(claim_id=cid,domain=domain,claim_display=claim+('*' if star else ''),claim=claim,disease=disease,
      cell_state=state,adversarial_status=status,rigor_score_0_10=score,discovery_evidence=discovery,
      independent_replication=replication,orthogonal_evidence=orthogonal,direct_vs_inferred=directness,
      principal_threat_to_validity=major_risk,adversarial_decision=decision,decisive_next_test=test,
      asterisk_definition='Independent cohort: concordant direction, raw P<0.05, no significant directional conflict' if star else ''))

# IPF circuit: break the composite model into falsifiable components.
add('IPF-01','Integrated mechanism','FN1–αv-integrin–TEAD2/ZNF322 aberrant-epithelial circuit','IPF','Aberrant basal / KRT5−KRT17+','High-priority hypothesis',8.0,
    'CellRef identity contrast: 30 vs 24 donors','TEAD2/ZNF322/CREB3L1 and multiple RNA/communication features confirm in 4 vs 4 donors',
    'Independent Xenium composition, sender/receiver expression, donor-aggregated proximity','Mixed direct RNA/Xenium plus inferred GRN/communication',
    'Discovery contrast is aberrant basal vs AT2, not IPF vs control within the same cell state; causality remains untested',
    'Retain as the lead experimentally testable model, not a demonstrated signaling circuit','Factorial FN1/SPP1/αvβ6/ROCK × TEAD2/ZNF322 perturbation in tri-culture',True)
add('IPF-02','Regulatory','TEAD2 activity increase','IPF','Aberrant basal vs AT2','Strong association',8.0,'Adams GRN Δ=1.274, FDR=1.0e−16; RNA log2FC=0.305, FDR=0.0014','Jaiswal confirms GRN and RNA direction','35 predicted targets replicate; Xenium ITGB6 supports downstream receiver','GRN inferred; RNA measured','Cell-identity contrast and regulon/expression circularity','Retain; do not call TEAD2 causal','TEAD2 CRISPRi plus CUT&RUN/ATAC and rescue',True)
add('IPF-03','Regulatory','ZNF322 activity increase','IPF','Aberrant basal vs AT2','Moderate association',6.5,'Adams GRN Δ=1.048, FDR=3.4e−18','Jaiswal confirms GRN; ZNF322 RNA does not independently confirm','18 predicted targets replicate','GRN inferred; RNA support incomplete','Regulon trained outside this exact state; target overlap partly circular','Retain as cofactor candidate below TEAD2','ZNF322 CRISPRi and TEAD2×ZNF322 epistasis',True)
add('IPF-04','Regulatory','CREB3L1 activity loss','IPF','Aberrant basal vs AT2','Moderate association',6.0,'Adams GRN Δ=−1.377, FDR=1.6e−15','Jaiswal confirms GRN; RNA does not confirm','Opposes TEAD adhesion-state shift','GRN inferred','Cell-identity contrast; unclear functional direction','Retain as marker, not driver','CREB3L1 rescue during FN1/αv stimulation',True)
add('IPF-05','Target module','TEAD2 adhesion/basement-membrane target module','IPF','Aberrant basal vs AT2','Strong molecular endpoint',7.0,'ITGB4/ITGA3/CDH3/JUP/HSPG2/LAMA5/PTK7/CTNND1/CDC42 increase','35 TEAD2 targets and 18 ZNF322 targets recur','Xenium ITGB6 and epithelial abundance','Measured RNA targets linked to inferred TFs','Target overlap is not statistically independent of regulon scoring','Retain as perturbation readout, not proof of TF binding','CUT&RUN plus expression rescue after TEAD2 perturbation',True)
add('IPF-06','Communication','Macrophage FN1→αvβ6/β8/β1 communication','IPF','Macrophage to aberrant epithelium','Prediction',4.5,'Fastcomm scores increase in Adams','Directional score changes recur in Jaiswal','Xenium FN1 sender and ITGB6 receiver expression','RNA-derived ligand–receptor score','No protein binding, directionality, or CellChat reanalysis; identity contrast','Retain only as a nominated edge','Matrix-defined coculture with blocking antibodies',True)
add('IPF-07','Communication','SPP1→αv communication from MT+ macrophages','IPF','MT+ AM to aberrant epithelium','Prediction',4.5,'Fastcomm increase in Adams','Recurs in Jaiswal','Xenium SPP1 and FN1 in disease macrophages','RNA-derived communication','Same limitations as IPF-06','Retain as secondary sender edge','SPP1 depletion/secretion capture with αv blockade',True)
add('IPF-08','Composition','KRT5−/KRT17+ disease-state expansion','IPF','KRT5−/KRT17+ epithelium','Strong orthogonal finding',9.0,'Xenium 26 affected vs 9 control donors; g=3.04, FDR=5.4e−5','CellRef enrichment in two cohorts','1.63M-cell independent Xenium dataset','Direct cell-state abundance','Annotation and zero inflation; observational','Retain as one of the strongest results','Histologic validation with blinded state quantification',True)
add('IPF-09','Composition','Activated fibrotic fibroblast expansion','IPF','Activated fibrotic fibroblasts','Strong orthogonal finding',8.5,'Xenium g=3.19, FDR=5.4e−5','Consistent with CellRef fibrotic-state enrichment','Independent Xenium','Direct abundance','Annotation and regional sampling','Retain','Histology and donor-matched matrix deposition assay',True)
add('IPF-10','Composition','AT2 depletion','IPF','AT2','Strong orthogonal finding',8.5,'Xenium g=−1.60, FDR=1.8e−4','Consistent with alveolar identity loss','Independent Xenium','Direct abundance','End-stage tissue and regional sampling','Retain as disease association','Longitudinal/injury-repair model',True)
add('IPF-11','Sender validation','FN1 increase in SPP1 macrophages','IPF','SPP1+ macrophages','Strong orthogonal finding',8.5,'Xenium donor pseudobulk g=1.91, FDR=7.1e−4','Independent of CellRef discovery','Direct spatial transcript measurement','Direct RNA, not secreted protein','Does not prove that macrophage FN1 reaches epithelial integrins','Retain as sender-expression evidence','Protein localization and macrophage-specific FN1 perturbation',True)
add('IPF-12','Sender validation','FN1 increase in activated fibrotic fibroblasts','IPF','Activated fibrotic fibroblasts','Strong orthogonal finding',8.5,'Xenium g=1.77, FDR=0.0066','Independent of CellRef discovery','Direct spatial transcript measurement','Direct RNA, not matrix protein','Does not prove matrix organization or signaling','Retain','FN1 matrix imaging and fibroblast-specific perturbation',True)
add('IPF-13','Receiver validation','ITGB6 increase in KRT5−/KRT17+ epithelium','IPF','KRT5−/KRT17+','Moderate orthogonal finding',7.0,'Xenium g=1.70, FDR=0.050','Independent of CellRef discovery','Spatial transcript measurement','Direct RNA','Borderline multiplicity-adjusted result; only 3 control expression pseudobulks','Retain with exact FDR and uncertainty','Protein-level αvβ6 quantification',True)
add('IPF-14','Receiver validation','ITGB6 increase in transitional AT2','IPF','Transitional AT2','Strong orthogonal finding',8.0,'Xenium g=1.64, FDR=7.1e−4','Independent of CellRef discovery','Spatial transcript measurement','Direct RNA','Observational','Retain','Protein and functional ligand-binding assay',True)
add('IPF-15','Spatial','KRT-state proximity to fibrotic fibroblasts','IPF','KRT–fibroblast niche','Moderate orthogonal association',6.5,'Donor-aggregated Squidpy neighborhood score: affected vs control MW P=0.0022','Independent Xenium donors','Independent Xenium spatial geometry','Direct proximity, not signaling','Proximity is nondirectional and abundance-sensitive','Retain as donor-level spatial association','Perturb sender cells and quantify receiver-state entry',False)
add('IPF-16','Spatial','KRT-state proximity to SPP1 macrophages','IPF','KRT–macrophage niche','Descriptive only',4.0,'Donor-aggregated affected vs control MW P=0.17','Does not reach donor-level significance','Independent Xenium spatial geometry','Direct proximity, not signaling','Original region-level view overstated independence','Downgrade; do not call spatial validation','Abundance-conditioned spatial null and larger control set',False)
add('IPF-17','Regulatory','ARNT2 increase in IPF AT2','IPF','AT2','Strong reporter; driver unproven',7.0,'Natri GRN Δ=1.733 FDR=0.020; RNA log2FC=0.279 FDR=0.0083','Adams and Jaiswal confirm both GRN/RNA direction','Four-cohort pooled g=1.31','Measured RNA plus inferred regulon','Regulon learned mainly in injury/basal contexts','Retain as replicated reporter, not driver','AT2-specific ARNT2 perturbation with direct target assay',True)

# Injury and predicted multimodal claims.
add('INJ-01','Acute injury','ETV5/XBP1/CEBPA/MLXIPL regulatory loss','COVID/pneumonia','AT2→AT1 intermediate vs AT2','Replicated state association',6.0,'COVID donor contrast FDR-significant','Pneumonia contrast concordant','Surfactant RNA loss in same state','GRN inferred; RNA markers direct','State-vs-state anchors, not longitudinal transition or disease effect','Retain as state checkpoint hypothesis','Time-resolved lineage organoid with TF rescue',True)
add('INJ-02','Acute injury','SFTPA1/SFTPA2/SFTPC/MFSD2A loss','COVID/pneumonia','AT2→AT1 intermediate vs AT2','Strong state identity association',7.0,'Large measured RNA losses in COVID','Pneumonia concordant','Biologically coherent AT2 identity loss','Direct RNA','Identity contrast cannot prove failed repair','Retain as phenotype, not mechanism','Time-course lineage and surfactant function',True)
add('INJ-03','Lipid prediction','PC(20:4/22:6) depletion','COVID/pneumonia','AT2→AT1 intermediate vs AT1','Unvalidated prediction',3.0,'RNA-imputed effect −0.47 in COVID','RNA-imputed effect −0.58 in pneumonia','Named LungMAP lipid mapping','Imputed, not lipidomics','Same expression model can reproduce correlated predictions','Do not present as measured lipid biology','Isomer-resolved LC–MS/MS',True)
add('INJ-04','Lipid prediction','PC(20:4/20:4) depletion','COVID/pneumonia','AT2→AT1 intermediate vs AT1','Unvalidated prediction',3.0,'RNA-imputed effect −0.44','RNA-imputed effect −0.48','Named LungMAP lipid mapping','Imputed','No orthogonal assay','Prediction only','Isomer-resolved LC–MS/MS',True)
add('INJ-05','Lipid prediction','PE(P-16:0/20:4) depletion','COVID/pneumonia','AT2→AT1 intermediate vs AT1','Unvalidated prediction',3.0,'RNA-imputed effect −0.30','RNA-imputed effect −0.27','Plasmalogen biochemical mapping','Imputed','No orthogonal assay','Prediction only','Isomer-resolved LC–MS/MS',True)
add('INJ-06','ADT prediction','PDPN/CD55 gain and CD49f loss','COVID/pneumonia','AT2→AT1 intermediate','Unvalidated prediction',2.5,'RNA-imputed ADT in COVID','Directional support in pneumonia','Compatible with transitional phenotype','Imputed, not protein','No CITE-seq/flow confirmation','Use only for panel design','Prospective flow/CITE-seq',False)

# COPD: preserve regulatory findings that survive, explicitly reject overclaims.
add('COPD-01','Regulatory','PVEC FOXA2/SREBF2/KLF5/TEAD2/TEAD3 identity-program loss','COPD','PVEC','Power-aware replicated GRN association',7.0,'Adams 17 vs 20; individual GRN FDR 0.004–0.026','UPenn confirms direction at raw P<0.05','Coherent vascular/lipid-regulatory module','Inferred GRN scores; corresponding RNA mostly nonsignificant','Regulon biology, not TF protein/activity measurement','Retain as regulatory-program loss; retract RNA-loss wording','Endothelial chip with TF rescue and barrier/lipid assays',True)
add('COPD-02','Regulatory','AT2 TP63 basal priming','COPD','AT2','Downgraded',2.0,'COPD-full RNA log2FC −0.257, FDR=0.812; Adams GRN FDR=0.552','UPenn gives a significant opposite RNA direction','Older pooled GRN score was positive','GRN inferred; RNA fails','Largest-cohort discovery fails and direction conflicts','Do not call replicated','Prospective cohort and protein/lineage assay',False)
add('COPD-03','Regulatory','AT2 FOXQ1 secretory priming','COPD','AT2','Downgraded',2.0,'COPD-full RNA log2FC 0.010, FDR=0.883; Adams GRN FDR=0.552','No RNA confirmation','Older pooled GRN score positive','GRN inferred; RNA fails','Largest-cohort discovery fails','Do not call replicated','Prospective cohort and perturbation',False)
add('COPD-04','Endothelial RNA','PVEC FBLIM1/CHD6/DACH1 RNA loss','COPD','PVEC','Downgraded',2.5,'UPenn FBLIM1 FDR=0.553 (7 vs 13)','No eligible independent FDR discovery','Separate GRN program survives','Direct RNA','Small discovery cohort and no multiplicity support','Retract as replicated RNA loss','Larger uncensored PVEC cohort',False)

# Exact cross-disease findings surviving discovery+confirmation in both diseases.
for i,(pop,gene,desc) in enumerate([('AdvFB','RGCC','vascular/quiescence-associated stromal program'),('AdvFB','CXCL12','chemokine-rich adventitial niche'),
 ('AdvFB','CFI','complement-regulatory stromal program'),('IM','SGK1','stress-responsive macrophage program'),('IM','USP53','interstitial-macrophage state'),
 ('AT1','PLPP3','phospholipid/phosphatase endothelial-like repair signal'),('AM','CD163','loss of resident-scavenger identity'),('AM','PLIN2','loss of lipid-droplet program')],1):
    z=cross[(cross.population==pop)&(cross.gene==gene)]
    ev='; '.join(f"{r.disease}: {r.discovery_cohort} log2FC={r.discovery_log2fc:.2f}" for _,r in z.iterrows())
    cf='; '.join(f"{r.disease}: {r.confirmation_cohorts}" for _,r in z.iterrows())
    add(f'XDIS-{i:02d}','Cross-disease RNA',f'{pop} {gene}: {desc}','IPF + COPD',pop,'Confirmed in both diseases',7.0,ev,cf,'Four disease/cohort evidence chains','Direct RNA','Shared association does not imply shared upstream cause','Retain as reproduced biology','Matched perturbation in IPF and COPD primary cells',True)

# Isolation and drug hypotheses.
add('ISO-01','Isolation strategy','EPCAM+ CD66c-high PDPN+ CD55+ CD49a+ aberrant epithelium','IPF','KRT5−/KRT17+','Experimental panel',4.0,'ADT direction from Adams','No direct protein replication','RNA ITGA3/ITGB6 and Xenium receiver phenotype','Mostly RNA-imputed ADT','Marker coexpression has not been measured in single cells','Useful proposal, not validated gate','Prospective flow/CITE-seq and post-sort RNA',False)
add('ISO-02','Isolation strategy','CD90-high CD49a+ CD140a+ subpleural fibroblasts','IPF','SPFB','Experimental panel',3.5,'RNA-imputed ADT','Two predicted cohorts','Stromal identity markers','Imputed','No direct protein assay','Proposal only','Prospective flow/CITE-seq',False)
add('ISO-03','Isolation strategy','CD1c+ HLA-DR-high CD16-low with CD207 confirmation','IPF','Langerhans','Near-term panel',5.5,'Strong CD1A/CD207 RNA identity','Adams and Natri state evidence','Established lineage markers','Mixed measured RNA/imputed ADT','Disease enrichment versus identity can be conflated','Reasonable validation panel','Flow plus post-sort transcript identity',True)
for j,(drug,mech,nline,nsig,rho) in enumerate([('Y-39983','ROCK',8,9,-.259),('fostamatinib','SYK',7,11,-.209),('PD-0325901','MEK1/2',5,12,-.202),('selumetinib','MEK1/2',5,6,-.202),('Y-27632','ROCK',3,3,-.189)],1):
    add(f'DRUG-{j:02d}','Drug reversal',f'{drug} ({mech}) reversal','IPF + COPD','IPF aberrant epithelium + COPD IM','Perturbational prioritization',5.0,
      f'Top-0.5% reversal; {nsig} signatures in {nline} cell lines; median rho={rho:.3f}','Both input programs require power-aware discovery + confirmation','Extra-pulmonary mechanism evidence reviewed separately','In vitro perturbation signatures','Cell-line context, cytostasis, exposure and target specificity','Prioritize for noncytotoxic ex vivo testing, not treatment inference','Dose-matched primary lung tri-culture',True)
add('NEG-01','Negative result','Ionocyte disease enrichment','IPF/COPD','Ionocyte','Rejected',8.0,'Did not pass replicated positive-enrichment screen','No independent positive reproduction','None','Direct negative screen','Power may still be limited for rare cells','Do not promote','Only revisit with purpose-built rare-cell sampling',False)
add('BLOCK-01','Blocked analysis','CellChat network reanalysis','IPF','Multicellular network','Blocked',1.0,'RData contains CellChat object','Compatible CellChat package unavailable','Squidpy spatial matrices analyzed instead','Not performed','Software/object compatibility','Do not imply CellChat confirmation','Restore compatible R environment and rerun',False)

claims=pd.DataFrame(claims).sort_values(['rigor_score_0_10','claim_id'],ascending=[False,True])
claims.to_csv(TAB/'adversarial_claim_inventory.tsv',sep='\t',index=False)

# Full finding atlas: every confirmed power-aware row, not just a seven-row shortlist.
atlas=disc.copy(); atlas['is_state_identity_contrast']=atlas.population.str.contains('__vs__',regex=False)
atlas['finding_display']=atlas.disease+' | '+atlas.population+' | '+atlas.display_label+' | '+atlas.modality
atlas.to_csv(TAB/'complete_power_aware_finding_atlas.tsv.gz',sep='\t',index=False)

# Power-aware existing-drug summary from the updated engine.
ENG=EXT/'drug_engine_power_aware'; meta=pd.read_csv(ENG/'signature_meta.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('lmd_key')
existing=['nintedanib','pirfenidone','roflumilast','budesonide','fluticasone','prednisone','dexamethasone','salmeterol','formoterol']
qtab=pd.read_csv(TAB/'power_aware_drug_query_programs.tsv',sep='\t'); er=[]
for q in qtab['query']:
    z=np.load(ENG/f'{q}.npz',allow_pickle=False); s=pd.DataFrame({k:z[k] for k in ['rho','n_shared','percentile']},index=z['key'].astype(str)).join(meta)
    s=s[(s.n_shared>=50)&s.compound_name.str.lower().isin(existing)].copy(); s['query']=q; er.append(s.reset_index(names='lmd_key'))
er=pd.concat(er,ignore_index=True); er['compound']=er.compound_name.str.lower()
era=er.groupby(['compound','query']).agg(n_signatures=('rho','size'),median_rho=('rho','median'),fraction_positive=('rho',lambda x:(x>0).mean()),fraction_negative=('rho',lambda x:(x<0).mean())).reset_index()
era['risk_flag']=(era.median_rho>=.05)&(era.fraction_positive>=.60)
era.to_csv(TAB/'existing_drug_state_scores_power_aware.tsv',sep='\t',index=False)

# Excel replaces the inadequate seven-row file and also gets a stable named copy.
limitations=pd.DataFrame([
 ('Cell-state identity confounding','Aberrant basal vs AT2 and injury intermediate vs AT2/AT1 contrasts do not isolate within-state disease effects.'),
 ('Imputed modalities','ADT and lipid values are RNA-imputed and are not orthogonal protein/lipid measurements.'),
 ('GRN circularity','TF activity and target overlap share expression-derived information; target coherence is not independent proof of TF binding.'),
 ('Communication direction','Ligand–receptor scores and proximity do not establish signaling direction or causality.'),
 ('Spatial nesting','Regions are nested within donors; revised figures aggregate regions to donor before visual inference.'),
 ('Cohort imbalance','Discovery uses largest effective N; small cohorts can confirm but do not carry equal evidentiary weight.'),
 ('Drug transferability','Perturbation cell lines, doses, cytostasis, and tissue context limit therapeutic inference.'),
 ('Novelty','No systematic patent/full-literature novelty review was completed.'),
 ('CellChat','Compatible CellChat environment was unavailable; no CellChat network result is claimed.')],columns=['risk','consequence'])
for out in [TAB/'intuitive_evidence_matrix.xlsx',TAB/'adversarial_evidence_review.xlsx']:
    with pd.ExcelWriter(out,engine='openpyxl') as w:
        claims.to_excel(w,index=False,sheet_name='Claim inventory')
        atlas.to_excel(w,index=False,sheet_name='All confirmed findings')
        atlas[~atlas.is_state_identity_contrast & atlas.modality.eq('rna')].to_excel(w,index=False,sheet_name='Disease RNA atlas')
        cross.to_excel(w,index=False,sheet_name='Cross-disease')
        pd.read_csv(TAB/'ranked_drug_predictions_power_aware.tsv',sep='\t').to_excel(w,index=False,sheet_name='Drug priorities')
        era.to_excel(w,index=False,sheet_name='Existing drugs')
        pd.read_csv(TAB/'adt_isolation_strategies.tsv',sep='\t').to_excel(w,index=False,sheet_name='Isolation proposals')
        pd.read_csv(TAB/'lipid_name_and_pathway_mapping.tsv',sep='\t').to_excel(w,index=False,sheet_name='Lipid mapping')
        pd.read_csv(TAB/'cohort_power_summary.tsv',sep='\t').to_excel(w,index=False,sheet_name='Cohort power')
        limitations.to_excel(w,index=False,sheet_name='Validity threats')
        for ws in w.book.worksheets:
            ws.freeze_panes='A2'; ws.auto_filter.ref=ws.dimensions
            for col in ws.columns:
                letter=col[0].column_letter; ws.column_dimensions[letter].width=min(55,max(11,max(len(str(c.value or '')) for c in col[:200])+2))

# ---------------- Revised Figure 1: evidence chain with donor-level spatial aggregation ----------------
ipf=pd.read_csv(ROOT/'tables/ipf_multimodal_evidence.tsv',sep='\t')
fig=plt.figure(figsize=(15,10)); gs=fig.add_gridspec(2,2,hspace=.42,wspace=.32)
ax=fig.add_subplot(gs[0,0]); ax.axis('off'); ax.set_title('A  Testable circuit; arrows are hypotheses',loc='left',weight='bold')
nodes=[(.08,.62,'SPP1 macrophage\nFN1/SPP1', '#457b9d'),(.08,.22,'Fibrotic fibroblast\nFN1 matrix','#6a994e'),(.53,.43,'αvβ6/β8/β1\naberrant epithelium','#e9c46a'),(.86,.43,'TEAD2 / ZNF322\nadhesion program','#d95f59')]
for x,y,t,c in nodes: ax.text(x,y,t,ha='center',va='center',bbox=dict(boxstyle='round,pad=.65',fc=c,ec='white',alpha=.9),color='white' if c!='#e9c46a' else 'black',weight='bold',transform=ax.transAxes)
for y in [.62,.22]: ax.annotate('',xy=(.405,.43),xytext=(.20,y),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=2,color='#555'))
ax.annotate('',xy=(.735,.43),xytext=(.655,.43),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=2,color='#555'))
ax.text(.5,.05,'Measured support: sender FN1, receiver ITGB6, state abundance, donor-level proximity\nInferred: ligand→receptor direction and TF causality',ha='center',transform=ax.transAxes,fontsize=9)

ax=fig.add_subplot(gs[0,1]); sel=['TEAD2','ZNF322','CREB3L1','ITGA3','MMP7','ABCA3']; z=ipf[ipf.gene.isin(sel)&ipf.modality.isin(['grn_tf','rna'])].copy(); z['cohort']=z.contrast_id.str.split('__').str[0].replace({'Adams2020':'Adams','Jaiswal2026_ILD':'Jaiswal'})
z['row']=z.modality.str.replace('grn_tf','GRN',regex=False).str.replace('rna','RNA',regex=False)+' | '+z.gene
mat=z.pivot_table(index='row',columns='cohort',values='log2fc',aggfunc='first').reindex([f'GRN | {g}' for g in ['TEAD2','ZNF322','CREB3L1']]+[f'RNA | {g}' for g in ['ITGA3','MMP7','ABCA3']])
direction=np.sign(mat); ann=mat.map(lambda v:'' if pd.isna(v) else f'{v:.2f}')
sns.heatmap(direction,cmap='vlag',center=0,vmin=-1,vmax=1,annot=ann,fmt='',cbar=False,linewidths=.8,ax=ax)
ax.set_title('B  Discovery and confirmation agree in direction*\ncolor=direction; number=native modality effect',loc='left',weight='bold'); ax.set_xlabel(''); ax.set_ylabel('')

ax=fig.add_subplot(gs[1,0]); comp=pd.read_csv(ROOT/'tables/xenium_composition_stats.tsv',sep='\t'); xe=pd.read_csv(ROOT/'tables/xenium_expression_stats.tsv',sep='\t'); tests=xe[xe.group.str.contains('test')]
forest=[('KRT state abundance',float(comp.loc[comp.cell_type.eq('KRT5-/KRT17+'),'hedges_g_logit'].iloc[0]),float(comp.loc[comp.cell_type.eq('KRT5-/KRT17+'),'fdr_selected_states'].iloc[0])),
 ('Fibrotic FB abundance',float(comp.loc[comp.cell_type.eq('Activated Fibrotic FBs'),'hedges_g_logit'].iloc[0]),float(comp.loc[comp.cell_type.eq('Activated Fibrotic FBs'),'fdr_selected_states'].iloc[0])),
 ('AT2 abundance',float(comp.loc[comp.cell_type.eq('AT2'),'hedges_g_logit'].iloc[0]),float(comp.loc[comp.cell_type.eq('AT2'),'fdr_selected_states'].iloc[0]))]
for ct,g in [('SPP1 macrophage','FN1'),('Fibrotic FB','FN1'),('KRT state','ITGB6'),('Transitional AT2','ITGB6')]:
    full={'SPP1 macrophage':'SPP1+ Macrophages','Fibrotic FB':'Activated Fibrotic FBs','KRT state':'KRT5-/KRT17+','Transitional AT2':'Transitional AT2'}[ct]; r=tests[(tests.cell_type==full)&(tests.gene==g)].iloc[0]; forest.append((f'{ct}: {g}',r.q25_or_hedges_g_for_test,r.fdr_within_selected_tests))
ff=pd.DataFrame(forest,columns=['feature','g','fdr']).sort_values('g'); ax.scatter(ff.g,ff.feature,s=75,c=np.where(ff.fdr<=.05,'#2a9d8f','#999')); ax.axvline(0,color='black',lw=.8)
for y,(_,r) in enumerate(ff.iterrows()): ax.text(r.g+.08,y,f"q={r.fdr:.2g}",va='center',fontsize=8)
ax.set_xlabel('Hedges g, affected − unaffected'); ax.set_title('C  Independent Xenium donor effects\n26 affected vs 9 control for composition',loc='left',weight='bold')

ax=fig.add_subplot(gs[1,1]); sp=pd.read_csv(ROOT/'tables/xenium_spatial_pairs_donor_region.tsv',sep='\t'); pairs=[('KRT5-_KRT17+','Activated_Fibrotic_FBs','KRT↔fibrotic FB'),('KRT5-_KRT17+','SPP1+_Macrophages','KRT↔SPP1 macrophage')]; sr=[]
for a,b,label in pairs:
    q=sp[(sp.cell_a==a)&(sp.cell_b==b)].copy(); q['status']=np.where(q.group.eq('Unaffected'),'Control','Affected'); q=q.groupby(['donor_id','status'],as_index=False).squidpy_z.median();
    c=q[q.status=='Control'].squidpy_z.dropna(); d=q[q.status=='Affected'].squidpy_z.dropna(); p=mannwhitneyu(d,c,alternative='two-sided').pvalue if len(c)>0 and len(d)>0 else np.nan
    q['pair']=f'{label}\nMW p={p:.2g}'; sr.append(q)
sr=pd.concat(sr,ignore_index=True); sns.stripplot(data=sr,x='pair',y='squidpy_z',hue='status',dodge=True,jitter=.12,palette={'Control':'#6c757d','Affected':'#d95f59'},ax=ax); ax.axhline(0,color='black',lw=.7); ax.set_xlabel(''); ax.set_ylabel('median Squidpy z per donor'); ax.set_title('D  Spatial evidence aggregated to independent donors',loc='left',weight='bold'); ax.legend(frameon=False)
fig.suptitle('IPF FN1–integrin–TEAD circuit: what is measured, replicated, and still inferred',fontsize=16,weight='bold'); fig.tight_layout(rect=[0,0,1,.95])
for ext in ['png','pdf']: fig.savefig(BASEFIG/f'Figure1_IPF_FN1_integrin_TEAD_niche.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

# ---------------- Revised Figure 2: modality-specific scales ----------------
inf=pd.read_csv(ROOT/'tables/infection_multimodal_evidence.tsv',sep='\t'); inf['cohort']=inf.contrast_id.str.extract(r'__(COVID19|Pneumonia)_')[0]
fig,axs=plt.subplots(1,3,figsize=(15,7),gridspec_kw={'width_ratios':[1,1,1.15]})
panels=[('Regulatory scores (RNA-derived)','grn_tf',['ETV5','XBP1','CEBPA','MLXIPL']),('Measured RNA','rna',['SFTPA1','SFTPA2','SFTPC','MFSD2A']),('Predicted modalities—not measured',None,['Podoplanin','CD55','CD49f','PC(20:4/22:6)','PC(20:4/20:4)','PE(P-16:0/20:4)'])]
for ax,(title,mod,genes) in zip(axs,panels):
    z=inf[inf.gene.isin(genes)&(inf.modality.eq(mod) if mod else inf.modality.isin(['adt','lipid']))]; mat=z.pivot_table(index='gene',columns='cohort',values='log2fc',aggfunc='mean').reindex(genes)
    lim=np.nanmax(np.abs(mat.to_numpy())); sns.heatmap(mat,cmap='vlag',center=0,vmin=-lim,vmax=lim,annot=True,fmt='.2f',linewidths=.6,cbar=False,ax=ax)
    ax.set_title(title,weight='bold'); ax.set_xlabel(''); ax.set_ylabel(''); ax.tick_params(axis='y',rotation=0)
axs[0].text(.5,-.16,'reference: AT2',ha='center',transform=axs[0].transAxes); axs[1].text(.5,-.16,'reference: AT2',ha='center',transform=axs[1].transAxes); axs[2].text(.5,-.16,'ADT vs AT2; lipids vs AT1',ha='center',transform=axs[2].transAxes)
fig.suptitle('Acute-injury transitional state: replicated directions, separated evidence scales',fontsize=15,weight='bold'); fig.text(.5,.01,'Agreement across COVID-19 and pneumonia is replication of a state contrast—not direct longitudinal proof of transition.',ha='center',fontsize=9); fig.tight_layout(rect=[0,.04,1,.94])
for ext in ['png','pdf']: fig.savefig(BASEFIG/f'Figure2_infection_surfactant_lipid_checkpoint.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

# ---------------- Revised Figure 3: COPD survivors and retractions ----------------
fig,axs=plt.subplots(1,3,figsize=(16,6),gridspec_kw={'width_ratios':[.85,1.25,1.1]})
cx=disc[(disc.disease=='COPD')&disc.independently_confirmed&~disc.population.str.contains('__vs__',regex=False)]
cnt=cx.groupby(['population','modality']).size().unstack(fill_value=0); cnt=cnt.loc[cnt.sum(1).nlargest(10).index].sort_values('rna'); cnt.plot.barh(stacked=True,color={'rna':'#2a9d8f','grn_tf':'#e9c46a','fastcomm':'#457b9d','adt':'#9b5de5','lipid':'#f15bb5'},ax=axs[0]); axs[0].set_title('A  Confirmed features by state*',loc='left',weight='bold'); axs[0].set_xlabel('feature count'); axs[0].set_ylabel(''); axs[0].legend(fontsize=7,frameon=False)
sel=[('AdvFB','RGCC'),('AdvFB','CXCL12'),('IM','SGK1'),('IM','USP53'),('AT1','PLPP3'),('AM','CD163'),('AM','PLIN2')]; z=[]
for pop,g in sel:
    r=pa[(pa.disease=='COPD')&(pa.population==pop)&(pa.gene==g)&(pa.modality=='rna')].iloc[0]; z.append((f'{pop} | {g}*',r.discovery_log2fc,r.discovery_fdr,r.discovery_cohort))
zz=pd.DataFrame(z,columns=['feature','fc','fdr','cohort']).sort_values('fc'); axs[1].barh(zz.feature,zz.fc,color=np.where(zz.fc>0,'#d95f59','#457b9d')); axs[1].axvline(0,color='black',lw=.7)
for y,(_,r) in enumerate(zz.iterrows()): axs[1].text(r.fc+(.04 if r.fc>0 else -.04),y,f"{r.cohort}; q={r.fdr:.2g}",ha='left' if r.fc>0 else 'right',va='center',fontsize=7)
axs[1].set_title('B  Reproduced COPD RNA examples*',loc='left',weight='bold'); axs[1].set_xlabel('largest-cohort log2FC')
down=pd.DataFrame([('AT2 TP63 RNA',-.257,.812,'fails; opposite cohort'),('AT2 FOXQ1 RNA',.010,.883,'fails'),('PVEC FBLIM1 RNA',-.233,.553,'fails'),('PVEC FOXA2 GRN*',-.382,.0098,'confirmed'),('PVEC SREBF2 GRN*',-.677,.0188,'confirmed'),('PVEC KLF5 GRN*',-.570,.0043,'confirmed')],columns=['claim','effect','fdr','note']).sort_values('effect'); down['label']=down.apply(lambda r:f"{r['claim']}\nq={r.fdr:.2g}; {r.note}",axis=1)
axs[2].barh(down.label,down.effect,color=np.where(down.claim.str.endswith('*'),'#2a9d8f','#b8b8b8')); axs[2].axvline(0,color='black',lw=.7)
axs[2].set_title('C  Regulatory program survives;\nRNA overclaims do not',loc='left',weight='bold'); axs[2].set_xlabel('discovery effect (native units)')
fig.suptitle('COPD after power-aware discovery and independent confirmation',fontsize=15,weight='bold'); fig.text(.5,.01,'* concordant raw P<0.05 in an independent cohort, with no significant opposing cohort.',ha='center'); fig.tight_layout(rect=[0,.04,1,.94])
for ext in ['png','pdf']: fig.savefig(BASEFIG/f'Figure3_TF_candidates_and_COPD_endothelium.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

# ---------------- Revised Figure 4 and 5: adversarial evidence audit ----------------
show=claims[~claims.adversarial_status.isin(['Experimental panel','Near-term panel'])].head(30).sort_values('rigor_score_0_10').copy(); show.loc[show.claim_id.eq('NEG-01'),'claim_display']='Ionocyte enrichment [rejected negative result]'
palette={'Strong orthogonal finding':'#1b9e77','Strong association':'#1b9e77','High-priority hypothesis':'#66a61e','Confirmed in both diseases':'#2a9d8f','Power-aware replicated GRN association':'#2a9d8f','Moderate orthogonal finding':'#e6ab02','Moderate association':'#e6ab02','Strong molecular endpoint':'#e6ab02','Moderate orthogonal association':'#e6ab02','Replicated state association':'#e6ab02','Strong state identity association':'#e6ab02','Strong reporter; driver unproven':'#e6ab02','Prediction':'#7570b3','Perturbational prioritization':'#7570b3','Unvalidated prediction':'#d95f59','Downgraded':'#d95f59','Rejected':'#555','Blocked':'#555'}
fig,ax=plt.subplots(figsize=(11,10)); ax.barh(show.claim_display,show.rigor_score_0_10,color=[palette.get(x,'#999') for x in show.adversarial_status]); ax.axvline(7,color='#2a9d8f',ls='--',lw=1); ax.axvline(4,color='#d95f59',ls='--',lw=1); ax.set_xlim(0,10); ax.set_xlabel('adversarial rigor score (0–10)'); ax.set_title('Adversarial claim audit: direct orthogonal evidence outranks model-derived agreement',weight='bold')
fig.text(.5,.01,'Scores summarize independence, power, directness, multiplicity, and confounding; see workbook for the explicit rationale.',ha='center',fontsize=9); fig.tight_layout(rect=[0,.04,1,1])
for ext in ['png','pdf']: fig.savefig(BASEFIG/f'Figure4_replicated_disease_enriched_states.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

top=claims.head(32).copy(); dims=['powered discovery','independent cohort','orthogonal modality','direct measurement','causal evidence']
M=[]
for _,r in top.iterrows():
    M.append([2 if r.rigor_score_0_10>=6 else 1, 2 if '*' in r.claim_display else 0, 2 if any(k in r.orthogonal_evidence.lower() for k in ['xenium','direct','four disease']) else (1 if r.orthogonal_evidence not in ['None',''] else 0), 2 if r.direct_vs_inferred.lower().startswith('direct') or 'measured rna' in r.direct_vs_inferred.lower() else 1, 0])
M=np.array(M); fig,ax=plt.subplots(figsize=(10,12)); sns.heatmap(M,annot=np.where(M==2,'strong',np.where(M==1,'partial','—')),fmt='',cmap=sns.color_palette(['#f2f2f2','#f4a261','#2a9d8f'],as_cmap=True),vmin=0,vmax=2,cbar=False,xticklabels=dims,yticklabels=top.claim_display,linewidths=.5,ax=ax)
ax.set_title('Expanded evidence matrix: 32 principal claims (not a seven-row shortlist)',weight='bold'); ax.set_xlabel(''); ax.set_ylabel(''); ax.tick_params(axis='x',rotation=25); fig.tight_layout()
for ext in ['png','pdf']: fig.savefig(FIG/f'Figure5_intuitive_evidence_matrix.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

# ---------------- Revised Figure 8: only power-aware existing-drug programs ----------------
order=[x for x in existing if x in era.compound.unique()]; qorder=qtab['query'].tolist(); mat=era.pivot(index='compound',columns='query',values='median_rho').reindex(index=order,columns=qorder); mat=mat.loc[:,mat.notna().any(axis=0)]; qorder=mat.columns.tolist()
short=[q.replace('_power_confirmed','').replace('IPF_aberrant_basal','IPF aberrant epi').replace('COPD_','COPD ') for q in qorder]
fig,ax=plt.subplots(figsize=(12,6)); sns.heatmap(mat,cmap='vlag',center=0,vmin=-.16,vmax=.16,linewidths=.5,xticklabels=short,cbar_kws={'label':'median Spearman ρ\nblue reversal; red mimicry'},ax=ax)
ax.set_title('Existing therapies: state effects against power-aware confirmed programs',weight='bold'); ax.set_xlabel(''); ax.set_ylabel(''); ax.tick_params(axis='x',rotation=45); ax.tick_params(axis='y',rotation=0)
fig.text(.5,.01,'Transcriptomic mimicry is a screening flag—not a clinical adverse event or off-target binding claim.',ha='center',fontsize=9); fig.tight_layout(rect=[0,.04,1,1])
for ext in ['png','pdf']: fig.savefig(FIG/f'Figure8_drug_reversal_and_state_liabilities.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

# ---------------- Revised Figure 9: cohort power plus balanced state yield ----------------
audit=pd.read_csv(TAB/'cohort_power_summary.tsv',sep='\t'); aud=audit.groupby(['disease','cohort'],as_index=False).agg(raw_total_n=('raw_total_n','max'),effective_n=('effective_n','max')); aud['label']=aud.disease+' | '+aud.cohort; aud=aud.sort_values('raw_total_n')
fig,axs=plt.subplots(1,2,figsize=(14,6),gridspec_kw={'width_ratios':[1,1.25]}); axs[0].barh(aud.label,aud.raw_total_n,color='#b9c6d8',label='raw total N'); axs[0].barh(aud.label,aud.effective_n,color='#264653',label='effective N'); axs[0].set_xlabel('donors'); axs[0].set_title('A  Effective N exposes group imbalance',loc='left',weight='bold'); axs[0].legend(frameon=False)
z=disc[disc.independently_confirmed&disc.modality.eq('rna')&~disc.population.str.contains('__vs__',regex=False)].groupby(['disease','population']).size().reset_index(name='n'); keep=pd.concat([z[z.disease==d].nlargest(10,'n') for d in z.disease.unique()],ignore_index=True); keep['label']=keep.disease+' | '+keep.population; keep=keep.sort_values('n'); axs[1].barh(keep.label,keep.n,color=np.where(keep.label.str.startswith('COPD'),'#d95f59','#457b9d')); axs[1].set_xlabel('confirmed RNA features'); axs[1].set_title('B  Confirmation yield by disease and state*',loc='left',weight='bold')
fig.suptitle('Power-aware discovery: cohort contribution and replicated state programs',weight='bold',fontsize=15); fig.text(.5,.01,'* largest-effective-N discovery FDR≤0.10 plus independent concordant raw P<0.05.',ha='center',fontsize=9); fig.tight_layout(rect=[0,.04,1,.94])
for ext in ['png','pdf']: fig.savefig(FIG/f'Figure9_power_aware_discovery_confirmation.{ext}',dpi=300 if ext=='png' else None,bbox_inches='tight')
plt.close(fig)

shutil.copy2(Path(__file__),EXT/'scripts'/Path(__file__).name)
print({'claims':len(claims),'all_confirmed_findings':len(atlas),'direct_disease_rna_findings':int((~atlas.is_state_identity_contrast&atlas.modality.eq('rna')).sum()),'existing_drug_rows':len(era)})
