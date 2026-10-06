"""Per-cohort union of concordant significant lipids over any cell-state combination."""
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from rank_lungmap_ipf_comparisons import OUT


def covering_states(g,exact=False):
    sets={state:set(x.measured_feature) for state,x in g.groupby('population')}
    target=set().union(*sets.values()) if sets else set()
    covered=set();greedy=[]
    while covered!=target:
        state=max(sorted(sets),key=lambda s:len(sets[s]-covered))
        gained=sets[state]-covered;covered|=gained
        greedy.append({'step':len(greedy)+1,'population':state,'new_lipids':len(gained),'cumulative_lipids':len(covered),
                       'new_lipid_names':'; '.join(sorted(gained))})
    chosen=[x['population'] for x in greedy];status='greedy cover; minimum size not established'
    if exact and target:
        states=sorted(sets);features=sorted(target)
        a=np.array([[f in sets[s] for s in states] for f in features],dtype=float)
        res=milp(np.ones(len(states)),integrality=np.ones(len(states)),bounds=Bounds(0,1),
                 constraints=LinearConstraint(a,1,np.inf),options={'time_limit':10})
        if res.success:
            chosen=[states[i] for i in np.flatnonzero(res.x>.5)]
            status='proven minimum number of cell states; alternative optimal combinations may exist'
    return target,greedy,chosen,status


def summarize(out=OUT):
    alltests=pd.read_csv(out/'all_measured_vs_imputed_tests.csv')
    alltests=alltests[~alltests.population.str.contains('__vs__')]
    overlap=pd.read_csv(out/'both_significant_fold_threshold_lipids.csv')
    good=overlap[overlap.direction_concordant_supplied].copy()
    measured=pd.read_csv(out/'experimental_lung_IPF_lipids_all_544.csv')
    maps=pd.read_csv(out/'lipid_matching_audit.csv')
    census=pd.read_csv(out/'population_replication_census.csv')
    common=set(maps.loc[maps.model.eq('candidate_MS1_47')&maps.status.eq('matched'),'measured_feature']) & set(maps.loc[maps.model.eq('current_as_deployed')&maps.status.eq('matched'),'measured_feature'])
    summary=[];supports=[];covers=[]
    for gate,source,cut in [('both_rawp005','supplied_raw_p',.05),('both_BH005','supplied_adjusted_p',.05),('both_BH010','supplied_adjusted_p',.1)]:
        experiment=set(measured.loc[measured[source].lt(cut),'feature'])
        for (contrast,unit,model),base in census.groupby(['contrast','unit','model']):
            predicted=set(maps.loc[maps.model.eq(model)&maps.status.eq('matched'),'measured_feature'])
            for panel in ['all_matched','shared_MS1_47']:
                eligible=experiment & predicted
                if panel=='shared_MS1_47':eligible &= common
                for fold in [1.,1.1,1.2,1.5,2.]:
                    g=good[(good.contrast==contrast)&(good.unit==unit)&(good.model==model)&(good.overlap_gate==gate)&good.minimum_imputed_fold.eq(fold)&good.measured_feature.isin(eligible)]
                    exact=(fold==1 and gate in ['both_rawp005','both_BH010'])
                    union,greedy,chosen,status=covering_states(g,exact)
                    keys=dict(contrast=contrast,unit=unit,model=model,panel=panel,overlap_gate=gate,minimum_imputed_fold=fold)
                    up=set(g.loc[g.supplied_log_effect.gt(0),'measured_feature']);down=set(g.loc[g.supplied_log_effect.lt(0),'measured_feature'])
                    summary.append(dict(**keys,experimental_differentials_total=len(experiment),matched_experimental_differentials=len(eligible),
                                        tested_populations=int((base.status.eq('tested') & ~base.population.str.contains('__vs__')).sum()),
                                        concordant_union=len(union),concordant_up=len(up),concordant_down=len(down),
                                        matched_not_recapitulated=len(eligible-union),supporting_cell_states=g.population.nunique(),
                                        cover_state_count=len(chosen),cover_states='; '.join(chosen),cover_status=status))
                    for r in greedy:covers.append(dict(**keys,**r))
                    for feature,x in g.groupby('measured_feature'):
                        best=x.loc[x.fdr.idxmin()] if gate!='both_rawp005' else x.loc[x.pval.idxmin()]
                        details='; '.join(f'{r.population}: {r.signed_fold:+.4g} fold, rawp={r.pval:.4g}, BH={r.fdr:.4g}' for r in x.sort_values('population').itertuples())
                        supports.append(dict(**keys,measured_feature=feature,experimental_direction='up' if best.supplied_log_effect>0 else 'down',
                                             experimental_raw_p=best.supplied_raw_p,experimental_BH=best.supplied_adjusted_p,
                                             concordant_cell_types='; '.join(sorted(x.population.unique())),state_fold_p_details=details))
    result=pd.DataFrame(summary);result.to_csv(out/'cohort_concordant_union_summary.csv',index=False)
    support=pd.DataFrame(supports);support.to_csv(out/'cohort_concordant_lipid_celltypes.csv',index=False)
    pd.DataFrame(covers).to_csv(out/'cohort_concordant_greedy_celltype_combinations.csv',index=False)
    primary=result[(result.panel=='all_matched')&result.overlap_gate.eq('both_BH010')&result.minimum_imputed_fold.eq(1)]
    cols=['contrast','unit','model','matched_experimental_differentials','tested_populations','concordant_up','concordant_down','concordant_union','cover_state_count','cover_states']
    lines=['# Concordant lipid coverage across cell types, separately by cohort','',
           'Question: how many experimentally changed IPF lipids can be recapitulated by any combination of cell types within one RNA cohort? Each lipid counts once when at least one ordinary cell state has a significant change in the experimental direction. Opposing states do not subtract or disqualify that concordant support. Cohorts and pseudobulk/metacell units are kept separate. Natri strata are shown individually and its combined contrast is already present; their counts must not be added.','',
           'Only results significant in both measurements enter the unions. Both raw p < 0.05, both BH < 0.05, and both BH < 0.10 are tabulated at five imputed fold thresholds. Experimental fold base remains unconfirmed. The 47-lipid MS1 model and 219-lipid exploratory model have different experimental coverage. Shared-panel comparisons are supplied. The union asks whether a lipid can be explained by any state, not whether its cellular origin has been experimentally validated.','',
           'Metacell p-values remain nominal tests of metacells; UPenn’s IPF arms have only one/two donors. Searching across states is exploratory and the reported within-state BH values do not correct that search.','',
           '## Both BH < 0.10: any concordant state, no fold cutoff','',primary[cols].to_markdown(index=False),'',
           '## Shared MS1 panel: current versus candidate, raw-p and BH unions','']
    commonrows=result[result.panel.eq('shared_MS1_47')&result.minimum_imputed_fold.eq(1)&result.overlap_gate.isin(['both_rawp005','both_BH010'])]
    lines += [commonrows[cols+['overlap_gate']].to_markdown(index=False),'',
              '## Cell-type combinations','',
              'For the primary raw-p and BH analyses without a fold cutoff, integer optimization finds the minimum number of cell states covering every concordant lipid where an optimum is certified. Alternative equally small combinations can exist. Other thresholds export a greedy coverage sequence, without asserting minimum size.','']
    for row in primary[primary.model.str.startswith('candidate')].itertuples():
        lines += [f'### {row.contrast} / {row.unit} / {row.model}', '',
                  f'{row.concordant_union} distinct concordant lipids ({row.concordant_up} up, {row.concordant_down} down), out of {row.matched_experimental_differentials} matched experimental differentials. Cover: {row.cover_states or "none"}. {row.cover_status}.','']
        s=support[(support.contrast==row.contrast)&(support.unit==row.unit)&(support.model==row.model)&support.panel.eq('all_matched')&support.overlap_gate.eq('both_BH010')&support.minimum_imputed_fold.eq(1)]
        lines += [s[['measured_feature','experimental_direction','concordant_cell_types','state_fold_p_details']].to_markdown(index=False),'']
    lines += ['## All fold thresholds and union counts','',result.to_markdown(index=False),'']
    (out/'COHORT_CONCORDANT_COVERAGE.md').write_text('\n'.join(lines))
    if (out/'EVALUATION_SUMMARY.md').exists():
        from render_lungmap_ipf_evaluation import render_concordant_header
        render_concordant_header(out)
    print(primary[cols].to_string(index=False))

if __name__=='__main__':summarize()
