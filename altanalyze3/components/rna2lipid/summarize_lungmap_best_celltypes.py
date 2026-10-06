"""Rank cell types solely by distinct concordant, both-significant lipids."""
import pandas as pd
from rank_lungmap_ipf_comparisons import OUT


def best_celltypes(out=OUT):
    counts=pd.read_csv(out/'both_significant_fold_threshold_summary.csv')
    union=pd.read_csv(out/'cohort_concordant_union_summary.csv')
    cross=pd.read_csv(out/'provenance/crosswalk.tsv',sep='\t')
    names=dict(zip(cross.cell_type,cross.cell_type_full))
    selected=counts[counts.panel.eq('all_matched') & counts.minimum_imputed_fold.eq(1) & counts.overlap_gate.isin(['both_rawp005','both_BH010']) & counts.model.str.startswith('candidate')].copy()
    keys=['contrast','unit','model','overlap_gate']
    best=selected[selected.concordant_total.eq(selected.groupby(keys).concordant_total.transform('max')) & selected.concordant_total.gt(0)].copy()
    best['cell_type_full']=best.population.map(names)
    best=best.merge(union[union.panel.eq('all_matched') & union.minimum_imputed_fold.eq(1)][keys+['matched_experimental_differentials']],on=keys)
    best['donor_supported']=best.n_case_donors.ge(3)&best.n_control_donors.ge(3)
    best.to_csv(out/'best_concordant_celltypes_by_cohort.csv',index=False)
    data=pd.read_csv(out/'both_significant_fold_threshold_lipids.csv')
    data=data[data.direction_concordant_supplied & data.minimum_imputed_fold.eq(1) & data.overlap_gate.isin(['both_rawp005','both_BH010']) & data.model.str.startswith('candidate')]
    intersections=[];lipids=[]
    for (model,unit,gate),g in data.groupby(['model','unit','overlap_gate']):
        for state,s in g.groupby('population'):
            a=s[s.contrast.eq('Adams2020__IPF_vs_Healthy')].set_index('measured_feature')
            n=s[s.contrast.eq('Natri2024__IPF_vs_Healthy')].set_index('measured_feature')
            shared=set(a.index)&set(n.index)
            if not shared:continue
            up=sum(a.loc[f,'supplied_log_effect']>0 for f in shared)
            intersections.append(dict(model=model,unit=unit,overlap_gate=gate,population=state,cell_type_full=names.get(state,state),
                                      shared_concordant_total=len(shared),shared_up=up,shared_down=len(shared)-up))
            for f in sorted(shared):
                lipids.append(dict(model=model,unit=unit,overlap_gate=gate,population=state,measured_feature=f,
                                   Adams_signed_fold=a.loc[f,'signed_fold'],Adams_rawp=a.loc[f,'pval'],Adams_BH=a.loc[f,'fdr'],
                                   Natri_signed_fold=n.loc[f,'signed_fold'],Natri_rawp=n.loc[f,'pval'],Natri_BH=n.loc[f,'fdr']))
    shared=pd.DataFrame(intersections)
    shared.to_csv(out/'Adams_Natri_same_celltype_concordant_counts.csv',index=False)
    pd.DataFrame(lipids).to_csv(out/'Adams_Natri_same_celltype_concordant_lipids.csv',index=False)
    winners=shared[shared.shared_concordant_total.eq(shared.groupby(['model','unit','overlap_gate']).shared_concordant_total.transform('max'))]
    bh=best[best.overlap_gate.eq('both_BH010')]
    primary=bh[bh.contrast.isin(['Adams2020__IPF_vs_Healthy','Natri2024__IPF_vs_Healthy'])]
    cols=['contrast','unit','model','population','cell_type_full','concordant_up','concordant_down','concordant_total','matched_experimental_differentials','n_case_donors','n_control_donors']
    lines=['# Best concordant cell types, ranked by concordant count','',
           'Best means the largest number of distinct experimental lipids with a significant imputed effect in the same direction. Opposing results do not subtract from the score. All ties are retained. The primary tables require BH < 0.10 in both measurements and no fold cutoff. Raw-p < 0.05 rankings are also exported. Cohorts, model panels and pseudobulk/metacell units remain separate. The full cell-type names follow the frozen LungMAP crosswalk; iMON is annotated Classical monocyte.','',
           '## Adams and combined Natri','',primary[cols].to_markdown(index=False),'',
           '## Best shared cell type across Adams and Natri','',
           winners[winners.overlap_gate.eq('both_BH010')].to_markdown(index=False),'',
           'A shared-state count requires the same lipid to be concordant and significant in that same named cell type in both cohorts. This is stricter than the any-cell-type cross-cohort intersection.','',
           '## All cohorts and strata','',bh[cols+['donor_supported']].to_markdown(index=False),'',
           'UPenn pseudobulk arms cannot meet the three-case replication minimum. Its metacell winners have only one/two IPF donors. Metacell tests use metacells rather than independent donors. Natri strata overlap donors. These data-driven rankings are exploratory, and primary effects have no magnitude threshold; significance does not establish large effects or causal cellular origins.','',
           '## Concordant lipids in each Adams/Natri winning cell type','']
    for r in primary.itertuples():
        g=data[data.contrast.eq(r.contrast)&data.unit.eq(r.unit)&data.model.eq(r.model)&data.population.eq(r.population)&data.overlap_gate.eq('both_BH010')]
        lines += [f'### {r.contrast} / {r.unit} / {r.model} / {r.population}', '',
                  g[['measured_feature','signed_fold','pval','fdr','supplied_raw_p','supplied_adjusted_p']].to_markdown(index=False,floatfmt='.5g'),'']
    (out/'BEST_CONCORDANT_CELLTYPES.md').write_text('\n'.join(lines))
    return primary,winners[winners.overlap_gate.eq('both_BH010')]


def add_best_sections(out=OUT):
    primary,winners=best_celltypes(out)
    cols=['contrast','unit','model','population','cell_type_full','concordant_up','concordant_down','concordant_total']
    marker='<!-- best-concordant-celltypes -->'
    block=marker+'\n## Best individual cell types by concordant count\n\n'+primary[cols].to_markdown(index=False)+'\n\n'+(
        'In the MS1 candidate, Adams transitioning monocyte-derived macrophages led pseudobulks with 12 concordant increases; classical monocytes (iMON) led metacells with 18 lipids (14 up, four down). In combined Natri, proliferating AT2 and bronchiolar multiciliated cells tied at 10 pseudobulk lipids (seven up/three down and 10 up/zero down, respectively). Bronchiolar and trachea/bronchus multiciliated metacells tied at 20 (16 up/four down and 14 up/six down). Classical monocytes were the strongest same-state metacell intersection across Adams and Natri, with 12 shared concordant lipids (11 up, one down); PBFB led the pseudobulk same-state intersection with four (two up, two down).\n\n'
        'The broader exploratory candidate favored transitioning macrophages in Adams (21 pseudobulk and 50 metacell lipids), proliferating AT2 in Natri pseudobulks (26), and trachea/bronchus multiciliated cells in Natri metacells (53). These are within-panel concordant-count maxima, not rankings by opposing results.\n\n'
        'Associated method: rank each cell state separately by distinct lipids significant at BH < 0.10 in both measurements and in the experimental direction, without a fold cutoff. Retain all ties. For the same-state cross-cohort ranking, intersect lipid identities within the identical canonical state across Adams and Natri. Counts are exploratory and metacell p-values are nominal.\n\n'
        '[Full cohort rankings, donor counts and the specific lipids](BEST_CONCORDANT_CELLTYPES.md).\n\n')+marker+'\n\n'
    for filename in ['EVALUATION_SUMMARY.md','ADAMS_NATRI_RESULTS_AND_METHODS.md']:
        path=out/filename
        if not path.exists():continue
        s=path.read_text()
        if marker in s:
            before,_,rest=s.partition(marker);_,_,after=rest.partition(marker);s=before+after.lstrip('\n')
        first,rest=s.split('\n',1);path.write_text(first+'\n\n'+block+rest.lstrip('\n'))
    print(primary[cols].to_string(index=False));print(winners.to_string(index=False))

if __name__=='__main__':add_best_sections()
