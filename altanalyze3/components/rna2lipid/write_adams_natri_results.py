"""Write cohort-specific Results/Methods drafts and explicit lipid intersections."""
from pathlib import Path
import pandas as pd
from rank_lungmap_ipf_comparisons import OUT

ADAMS='Adams2020__IPF_vs_Healthy'
NATRI='Natri2024__IPF_vs_Healthy'


def write(out=OUT):
    source=pd.read_csv(out/'cohort_concordant_lipid_celltypes.csv')
    counts=pd.read_csv(out/'cohort_concordant_union_summary.csv')
    primary=source[source.panel.eq('all_matched')&source.overlap_gate.eq('both_BH010')&source.minimum_imputed_fold.eq(1)]
    allsummary=[];details=[]
    for (model,unit,panel,gate,fold),g in source.groupby(['model','unit','panel','overlap_gate','minimum_imputed_fold']):
        a=g[g.contrast.eq(ADAMS)].set_index('measured_feature');n=g[g.contrast.eq(NATRI)].set_index('measured_feature')
        ai,ni=set(a.index),set(n.index);both=ai&ni;same=0
        for lipid in sorted(ai|ni):
            aa=a.loc[lipid] if lipid in ai else None;nn=n.loc[lipid] if lipid in ni else None
            ac=set(aa.concordant_cell_types.split('; ')) if aa is not None else set()
            nc=set(nn.concordant_cell_types.split('; ')) if nn is not None else set()
            shared=ac&nc
            if shared:same+=1
            direction=aa.experimental_direction if aa is not None else nn.experimental_direction
            if aa is not None and nn is not None and aa.experimental_direction!=nn.experimental_direction:raise ValueError('Experimental direction mismatch')
            details.append(dict(model=model,unit=unit,panel=panel,overlap_gate=gate,minimum_imputed_fold=fold,
                                measured_feature=lipid,experimental_direction=direction,
                                membership='both' if lipid in both else 'Adams only' if lipid in ai else 'Natri only',
                                Adams_cell_types='; '.join(sorted(ac)),Natri_cell_types='; '.join(sorted(nc)),
                                shared_cell_types='; '.join(sorted(shared)),
                                Adams_fold_p_details=aa.state_fold_p_details if aa is not None else '',
                                Natri_fold_p_details=nn.state_fold_p_details if nn is not None else ''))
        denom=counts[(counts.model==model)&(counts.unit==unit)&(counts.panel==panel)&(counts.overlap_gate==gate)&counts.minimum_imputed_fold.eq(fold)&counts.contrast.eq(ADAMS)].matched_experimental_differentials.iloc[0]
        up=sum(a.loc[f,'experimental_direction']=='up' for f in both)
        allsummary.append(dict(model=model,unit=unit,panel=panel,overlap_gate=gate,minimum_imputed_fold=fold,
                               matched_experimental_differentials=int(denom),Adams_concordant=len(ai),Natri_concordant=len(ni),
                               intersection=len(both),intersection_up=up,intersection_down=len(both)-up,
                               Adams_only=len(ai-ni),Natri_only=len(ni-ai),same_celltype_intersection=same))
    summary=pd.DataFrame(allsummary);detail=pd.DataFrame(details)
    summary.to_csv(out/'Adams_Natri_intersection_summary.csv',index=False)
    detail.to_csv(out/'Adams_Natri_lipid_intersection_celltypes.csv',index=False)
    main=summary[summary.panel.eq('all_matched')&summary.overlap_gate.eq('both_BH010')&summary.minimum_imputed_fold.eq(1)&summary.model.str.startswith('candidate')]
    ms=main[main.model.eq('candidate_MS1_47')].set_index('unit')
    c=counts[counts.panel.eq('all_matched')&counts.overlap_gate.eq('both_BH010')&counts.minimum_imputed_fold.eq(1)&counts.model.eq('candidate_MS1_47')].set_index(['contrast','unit'])
    lines=['# Adams and Natri: Results and associated Methods drafts','',
           '## Are these the same results?', '',
           'Adams2020 and Natri2024 are distinct RNA study labels with separate IPF/control comparisons. The same fitted candidate and the same Dr. Clair lung-tissue lipidomics comparator are used for both. The recovered lipid sets overlap substantially, but their supporting cell types need not be identical. This is cross-RNA-cohort directional agreement against one shared experimental lipidomics source; it is not replication in two independent lipidomics cohorts. Natri’s combined IPF contrast is used below; its anatomical fibrosis strata are not counted as separate replications.','',
           main[['model','unit','matched_experimental_differentials','Adams_concordant','Natri_concordant','intersection','intersection_up','intersection_down','same_celltype_intersection','Adams_only','Natri_only']].to_markdown(index=False),'',
           '## Version 1 — Adams Results', '',
           f'In the Adams IPF/control RNA cohort, the repaired MS1 candidate recapitulated the direction of 27 of 29 matched experimentally differential lipids (93.1%) in sample/state pseudobulks, comprising {c.loc[(ADAMS,"PB"),"concordant_up"]} increases and {c.loc[(ADAMS,"PB"),"concordant_down"]} decreases. Metacell analysis recapitulated 28 of 29 lipids (96.6%; 20 increases and eight decreases). A lipid was counted once when at least one cell state showed an experimentally concordant change with BH FDR < 0.10 in both the measured tissue lipids and the imputed cell-state comparison. A minimum covering combination of AT2, bronchiolar multiciliated cells, classical monocytes and transitioning monocyte-derived macrophages accounted for all 28 metacell-supported lipids. This combination describes where concordant predictions occur; cellular origins were not measured by the tissue lipidomics experiment.', '',
           '### Associated Methods — Adams', '',
           'The saved Adams2020 IPF-versus-healthy registry contrast includes 32 IPF and 28 control samples/donors. Real library/state RNA count matrices were aggregated by study, sample and canonical cell state before normalization to log2(1 + counts per 10,000), using full-transcriptome library totals. The separately saved metacell input was normalized in the same RNA expression space without a metacell cap. The evaluation candidate used the repaired 47-lipid MS1 target reference and all-reference linear kernel ridge regression (alpha=100), refit on the 1,293 candidate RNA genes present in LungMAP. A gene-wise additive alignment estimated from healthy, donor-balanced whole-sample RNA was applied equally to both disease arms. Experimental IPF lipid outcomes were not used to fit the candidate, alignment or alpha. Production models were not replaced.','',
           '## Version 2 — Natri Results', '',
           'In the combined Natri IPF/control RNA contrast, the repaired MS1 candidate recapitulated 26 of 29 matched experimentally differential lipids (89.7%) in sample/state pseudobulks, including 18 increases and eight decreases. Metacell analysis recapitulated all 29 matched lipids (20 increases and nine decreases). These counts represent the union of concordant, BH-significant predictions across cell states, with each lipid counted only once. A minimum covering combination of alveolar macrophages, AT2 cells and regulatory T cells accounted for all 29 metacell-supported lipids. These cell-state assignments constitute predictive explanations of tissue-level directional agreement rather than experimental localization of the lipid changes.', '',
           '### Associated Methods — Natri', '',
           'The combined Natri2024 IPF-versus-healthy registry contrast includes 48 IPF anatomical samples from 28 donors and 43 control samples/donors. Its less- and more-fibrotic strata share donors and were not combined as independent replications. RNA aggregation, normalization, healthy-study alignment and candidate inference followed the Adams workflow with Natri-specific healthy RNA offsets. The primary pseudobulk analysis retained the deployed sample/state units; it did not collapse repeated anatomical samples within each donor. Independent-donor whole-lung RNA aggregation was evaluated separately and is not the basis of the cell-state union counts reported here.','',
           '## Version 3 — Adams–Natri intersected Results', '',
           f'To identify lipids recovered across both RNA cohorts, we intersected the cohort-specific unions of concordant cell-state predictions. The MS1 candidate recovered {ms.loc["PB","intersection"]} of 29 matched experimentally differential lipids in both pseudobulk analyses (82.8%; 18 increases and six decreases) and {ms.loc["MC","intersection"]} of 29 in both metacell analyses (96.6%; 20 increases and eight decreases). All 28 Adams-supported metacell lipids were also supported in Natri; DG(18:0/20:4/0:0)_POS was supported only in Natri at this threshold. The lipid sets were therefore nearly identical at the metacell level, while their cellular explanations were less uniform: 26 of the 28 shared metacell lipids, and nine of the 24 shared pseudobulk lipids, had concordant support in at least one identically named cell state in both cohorts. The remaining shared lipids were supported by different cell states across cohorts.', '',
           '### Associated Methods — intersected analysis', '',
           'For each cohort independently, experimentally measured lipids with supplied BH-adjusted p < 0.10 were matched unambiguously to candidate targets by lipid acyl composition and ion mode. A matched lipid entered the cohort union if any ordinary within-state IPF/control comparison had BH-adjusted p < 0.10 and the same effect direction as the experimental lipid. Significant changes in other states with an opposing direction did not remove a lipid from the concordant union. Cohorts and pseudobulk/metacell units were not pooled. The cross-cohort set intersection required the same lipid identity to qualify independently in Adams and Natri; it did not require the same cell state. A separate analysis required at least one common supporting canonical cell state. Raw-p intersections and imputed fold thresholds of 1.0, 1.1, 1.2, 1.5 and 2.0 are exported as sensitivity analyses.','',
           '## Statistical methods and interpretation applicable to all versions', '',
           'Experimental directional evidence came from Dr. Clair’s supplied lung-tissue 10_results_with_statistics.csv (544 lipid features; 30 IPF tissue measurements and 10 controls). Sample-prefix grouping provisionally identifies 10 IPF donors with three tissues each and 10 control donors. The supplied raw and BH-adjusted statistics were used for the primary comparisons; donor-recalculated sensitivity analyses are retained in the full evaluation report. The supplied log-effect sign was used without assuming an unconfirmed log base. The MS1 candidate was trained on all 63 reference profiles (45 sorted-cell and 18 bulk profiles), with no held-out training profiles in this candidate fit. The wider candidate used 64 reference profiles (45 sorted-cell and 19 bulk profiles).', '',
           'Frozen cellHarmony functions were used with the deployed settings: empirical-Bayes moderated t-tests for pseudobulks (shrinkage=0.2; at least three rows per arm), and Scanpy Wilcoxon tests for metacells (tie_correct=False; at least five metacells per arm). BH adjustment was performed within each cell-state/contrast over each model’s full independently filtered lipid panel. Candidate predictions were explicitly identified as log2 abundance to prevent repeat normalization/log transformation. Signed fold changes were computed from ratios of arithmetic arm means after exact inverse transformation (2**prediction), without an abundance pseudocount. Primary union counts used no fold-magnitude cutoff. Minimum covering state combinations were obtained by binary integer optimization, with alternative equally small combinations possible.','',
           'The experimental file contains 356 features with supplied BH < 0.10; only 29 are unambiguously matched BH-significant features in the MS1 candidate. Thus 28/29 or 29/29 describes matched-panel coverage, not recovery of all 356 experimental features. RNA and lipidomics donors are unmatched, and both RNA cohorts use the same experimental comparator and candidate model. Metacell tests treat metacells as observations, and Natri sample pseudobulks can repeat donors. The cell-state union and covering combinations are exploratory; within-state BH does not correct the search over states or cohorts. Significant directional overlap does not establish recovery of experimental effect magnitudes: most predicted changes are small, and experimental fold magnitudes cannot be confirmed until the supplied transform/log base is known.','',
           '## Broader candidate Results — optional supplementary version', '',
           'The exploratory 219-lipid candidate matched 84 experimentally BH-significant lipids. Adams recapitulated 67 in pseudobulks (45 increases, 22 decreases) and 78 in metacells (51 increases, 27 decreases). Natri recapitulated 73 in pseudobulks (44 increases, 29 decreases) and all 84 in metacells (54 increases, 30 decreases). Their intersections comprised 60 pseudobulk lipids (38 increases, 22 decreases) and 78 metacell lipids (51 increases, 27 decreases). Common-cell-state support occurred for 28 of the 60 pseudobulk intersection lipids and 73 of the 78 metacell intersection lipids. The additional targets use the earlier bulk abundance gauge and lack the 47-panel raw-MS support; this wider analysis is reported separately as exploratory.', '',
           '## Every intersected and cohort-specific lipid with cellular support', '']
    pdmain=detail[detail.panel.eq('all_matched')&detail.overlap_gate.eq('both_BH010')&detail.minimum_imputed_fold.eq(1)&detail.model.str.startswith('candidate')]
    for (model,unit),g in pdmain.groupby(['model','unit']):
        lines += [f'### {model} / {unit}', '',g[['measured_feature','experimental_direction','membership','Adams_cell_types','Natri_cell_types','shared_cell_types','Adams_fold_p_details','Natri_fold_p_details']].to_markdown(index=False),'']
    lines += ['## Files and reproduction', '',
              'Adams_Natri_intersection_summary.csv retains every model, panel, significance gate and fold threshold. Adams_Natri_lipid_intersection_celltypes.csv provides lipid membership, cohort-specific and common cell types, signed folds and raw/BH values. Source cohort unions and minimum-covering combinations are in COHORT_CONCORDANT_COVERAGE.md.', '',
              '```bash','/opt/homebrew/bin/python3.11 components/rna2lipid/write_adams_natri_results.py','```','']
    (out/'ADAMS_NATRI_RESULTS_AND_METHODS.md').write_text('\n'.join(lines))
    if (out/'best_concordant_celltypes_by_cohort.csv').exists():
        from summarize_lungmap_best_celltypes import add_best_sections
        add_best_sections(out)
    print(main.to_string(index=False))

if __name__=='__main__':write()
