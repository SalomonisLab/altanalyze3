"""Exploratory, auditable ranking of dataset/state directional lipid support."""
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
OUT=HERE/'artifacts/LungMAP_IPF_candidate_20261004'


def rank(out=OUT):
    data=pd.read_csv(out/'all_measured_vs_imputed_tests.csv')
    census=pd.read_csv(out/'population_replication_census.csv')
    mapping=pd.read_csv(out/'lipid_matching_audit.csv')
    common=set(mapping.loc[mapping.model.eq('candidate_MS1_47') & mapping.status.eq('matched'),'measured_feature'])
    common &= set(mapping.loc[mapping.model.eq('current_as_deployed') & mapping.status.eq('matched'),'measured_feature'])
    rows=[]
    for panel in ['all_matched','shared_MS1_47']:
        selected=data[data.supplied_adjusted_p.lt(.1) & ~data.population.str.contains('__vs__')].copy()
        if panel=='shared_MS1_47':selected=selected[selected.measured_feature.isin(common)]
        for (contrast,unit,model,state),g in selected.groupby(['contrast','unit','model','population']):
            for gate,field,cut in [('rawp005','pval',.05),('BH010','fdr',.1)]:
                significant=g[g[field].lt(cut)]
                agree=significant[significant.direction_concordant_supplied]
                against=significant[~significant.direction_concordant_supplied]
                up=agree.supplied_log_effect.gt(0).sum();down=agree.supplied_log_effect.lt(0).sum()
                rows.append(dict(panel=panel,contrast=contrast,unit=unit,model=model,population=state,
                                 prediction_gate=gate,matched_experimental_differentials=g.measured_feature.nunique(),
                                 concordant_up=int(up),concordant_down=int(down),concordant_total=len(agree),
                                 discordant_total=len(against),net_concordant=len(agree)-len(against),
                                 significant_total=len(significant),direction_precision=len(agree)/len(significant) if len(significant) else np.nan))
    ranks=pd.DataFrame(rows).merge(census.drop_duplicates(['contrast','unit','model','population'])[['contrast','unit','model','population','n_case','n_control','n_case_donors','n_control_donors']],on=['contrast','unit','model','population'])
    ranks['donor_supported']=ranks.n_case_donors.ge(3)&ranks.n_control_donors.ge(3)
    ranks=ranks.sort_values(['net_concordant','concordant_total','direction_precision'],ascending=False,kind='stable')
    ranks['exploratory_rank']=ranks.groupby(['panel','unit','model','prediction_gate']).cumcount()+1
    ranks.to_csv(out/'ranked_dataset_celltype_comparisons.csv',index=False)
    # Paired comparison: identical experimental coverage; keep each model's original BH family.
    keys=['contrast','unit','population','prediction_gate']
    a=ranks[ranks.panel.eq('shared_MS1_47') & ranks.model.eq('candidate_MS1_47')]
    b=ranks[ranks.panel.eq('shared_MS1_47') & ranks.model.eq('current_as_deployed')]
    paired=a.merge(b,on=keys,suffixes=('_candidate','_current'))
    paired['concordant_gain']=paired.concordant_total_candidate-paired.concordant_total_current
    paired['net_gain']=paired.net_concordant_candidate-paired.net_concordant_current
    paired.to_csv(out/'paired_dataset_celltype_comparisons.csv',index=False)
    top=ranks[ranks.panel.eq('shared_MS1_47') & ranks.model.eq('candidate_MS1_47') & ranks.prediction_gate.eq('BH010') & ranks.donor_supported].groupby('unit',sort=False).head(12)
    cols=['unit','contrast','population','n_case_donors','n_control_donors','concordant_up','concordant_down','concordant_total','discordant_total','net_concordant','matched_experimental_differentials']
    lines=['# Exploratory dataset and cell-type selection','',
           'Selection uses the supplied experimental BH < 0.10 lipids and imputed raw p < 0.05 or BH < 0.10. Rankings prioritize concordant minus discordant lipids, then concordant count and direction precision. Both full coverage and the identical current/candidate MS1 panel are exported. The tables below use the shared panel and require at least three actual donors per arm. BH families retain each model’s full lipid panel. Metacell p-values still test metacells rather than donors.','',
           'These rankings explain where prediction agreement occurs; they do not by themselves make the model explainable AI or identify causal RNA drivers. Dataset/cell-state choices made using these outcomes are exploratory. Natri strata share donors and count as one study. Check a selected state in another independent study before claiming replicated support.','',
           '## Strongest supported MS1 candidate comparisons, BH < 0.10','',top[cols].to_markdown(index=False),'',
           '## All ranked comparisons','',ranks.to_markdown(index=False,floatfmt='.4g'),'']
    (out/'DATASET_CELLTYPE_SELECTION.md').write_text('\n'.join(lines))
    return top[cols]

if __name__=='__main__':print(rank().to_string(index=False))
