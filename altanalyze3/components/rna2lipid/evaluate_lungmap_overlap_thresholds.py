"""Both-significant overlap across fold thresholds and cell-state hypotheses."""
import numpy as np
import pandas as pd
from rank_lungmap_ipf_comparisons import OUT


def evaluate(out=OUT):
    d=pd.read_csv(out/'all_measured_vs_imputed_tests.csv')
    d=d[~d.population.str.contains('__vs__')].copy()
    census=pd.read_csv(out/'population_replication_census.csv')
    d=d.merge(census.drop_duplicates(['contrast','unit','model','population'])[['contrast','unit','model','population','n_case_donors','n_control_donors']],on=['contrast','unit','model','population'])
    mapping=pd.read_csv(out/'lipid_matching_audit.csv')
    shared=set(mapping.loc[mapping.model.eq('candidate_MS1_47') & mapping.status.eq('matched'),'measured_feature']) & set(mapping.loc[mapping.model.eq('current_as_deployed') & mapping.status.eq('matched'),'measured_feature'])
    allrows=[];summaries=[]
    for gate,expfield,predfield,cut in [('both_rawp005','supplied_raw_p','pval',.05),('both_BH005','supplied_adjusted_p','fdr',.05),('both_BH010','supplied_adjusted_p','fdr',.1)]:
        overlap=d[d[expfield].lt(cut)&d[predfield].lt(cut)].copy()
        for fold in [1.,1.1,1.2,1.5,2.]:
            selected=overlap[overlap.signed_fold.abs().ge(fold)&overlap.signed_fold.ne(0)].copy()
            selected['overlap_gate']=gate;selected['minimum_imputed_fold']=fold
            allrows.append(selected)
            for panel in ['all_matched','shared_MS1_47']:
                s=selected if panel=='all_matched' else selected[selected.measured_feature.isin(shared)]
                keys=['contrast','unit','model','population']
                base=d.groupby(keys)[['n_case_donors','n_control_donors']].first()
                counted=s.assign(both_significant_lipids=1,
                                 concordant_up=(s.direction_concordant_supplied & s.supplied_log_effect.gt(0)).astype(int),
                                 concordant_down=(s.direction_concordant_supplied & s.supplied_log_effect.lt(0)).astype(int),
                                 concordant_total=s.direction_concordant_supplied.astype(int),
                                 discordant_total=(~s.direction_concordant_supplied).astype(int))
                cols=['both_significant_lipids','concordant_up','concordant_down','concordant_total','discordant_total']
                result=base.join(counted.groupby(keys)[cols].sum()).fillna(0).reset_index()
                result[cols]=result[cols].astype(int)
                result['net_concordant']=result.concordant_total-result.discordant_total
                result['panel']=panel;result['overlap_gate']=gate;result['minimum_imputed_fold']=fold
                summaries.extend(result.to_dict('records'))
    summary=pd.DataFrame(summaries)
    summary.to_csv(out/'both_significant_fold_threshold_summary.csv',index=False)
    allrows=pd.concat(allrows,ignore_index=True)
    allrows.to_csv(out/'both_significant_fold_threshold_lipids.csv',index=False)
    hypotheses=[]
    for (contrast,unit,model,feature,gate,fold),g in allrows.groupby(['contrast','unit','model','measured_feature','overlap_gate','minimum_imputed_fold']):
        a=g[g.direction_concordant_supplied];b=g[~g.direction_concordant_supplied]
        if len(a) and len(b):
            def desc(r):return '; '.join(f'{x.population}: {x.signed_fold:+.3g} fold, p={x.pval:.3g}, BH={x.fdr:.3g}' for x in r.itertuples())
            hypotheses.append(dict(contrast=contrast,unit=unit,model=model,measured_feature=feature,overlap_gate=gate,
                                   minimum_imputed_fold=fold,concordant_states=desc(a),opposing_states=desc(b),
                                   interpretation='Hypothesized cell-state-specific opposing effects; cellular origin not experimentally validated'))
    hyp=pd.DataFrame(hypotheses);hyp.to_csv(out/'hypothesized_opposing_cellstate_effects.csv',index=False)
    main=summary[(summary.model=='candidate_MS1_47')&(summary.panel=='shared_MS1_47')&(summary.overlap_gate=='both_BH010')&summary.n_case_donors.ge(3)&summary.n_control_donors.ge(3)]
    top=main.sort_values(['net_concordant','concordant_total'],ascending=False,kind='stable').groupby(['unit','minimum_imputed_fold'],sort=False).head(5)
    cols=['unit','minimum_imputed_fold','contrast','population','both_significant_lipids','concordant_up','concordant_down','discordant_total']
    primary_hyp=hyp[(hyp.model=='candidate_MS1_47')&(hyp.overlap_gate=='both_BH010')&hyp.minimum_imputed_fold.eq(1.)]
    lines=['# Both-significant overlap and fold thresholds','',
           'Only the intersection of experimental and imputed significance is counted. Both raw p < 0.05, both BH < 0.05, and both BH < 0.10 are evaluated separately. No lipids significant in only one measurement enter these counts. Imputed absolute signed-fold thresholds are 1, 1.1, 1.2, 1.5 and 2. Experimental fold magnitudes have an unconfirmed log base, so no experimental magnitude threshold is inferred from that column. Each imputation retains its own full-panel BH correction.','',
           'Rankings are exploratory: changing cell states, datasets or fold thresholds using the same outcomes does not supply independent validation. Opposing states are recorded as potential cell-state-specific effects. Whole-lung measurements validate a lipid direction at tissue scope; they do not validate the assignment or the opposing cellular effects.','',
           '## Strongest MS1 candidate intersections at both BH < 0.10','',top[cols].to_markdown(index=False),'',
           '## Hypothesized opposing effects: candidate, both BH < 0.10, no fold cutoff','',primary_hyp.to_markdown(index=False),'',
           '## All intersection counts','',summary.to_markdown(index=False),'']
    (out/'BOTH_SIGNIFICANT_OVERLAP_REPORT.md').write_text('\n'.join(lines))
    print(top[cols].to_string(index=False))

if __name__=='__main__':evaluate()
