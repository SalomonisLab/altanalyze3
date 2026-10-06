"""Render all numerical LungMAP IPF validation results without modifying production."""
from pathlib import Path
import argparse
import json

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LABELS = ['current_as_deployed','current_explicit_log2','candidate_MS1_47','candidate_bulk_219_exploratory']


def table(df):
    return df.to_markdown(index=False, floatfmt='.4g')


def render_concordant_header(out):
    d=pd.read_csv(out/'cohort_concordant_union_summary.csv')
    d=d[d.panel.eq('all_matched') & d.overlap_gate.eq('both_BH010') & d.minimum_imputed_fold.eq(1) & d.model.str.startswith('candidate')].copy()
    d['recapitulated']=d.concordant_union.astype(str)+'/'+d.matched_experimental_differentials.astype(str)
    wide=d.pivot(index='contrast',columns=['model','unit'],values='recapitulated').reset_index()
    wide.columns=['contrast']+[f'{model} {unit}' for model,unit in wide.columns[1:]]
    marker='<!-- concordant-coverage -->'
    block=marker+'\n## Primary question: concordant coverage across any cell-type combination\n\n'+(
        'Count each lipid once per cohort/unit when any cell state recapitulates the experimental direction with BH < 0.10 in both. Opposing states do not subtract support. The experimental file has 356 BH < 0.10 features; the MS1 candidate matches 29 and the broader exploratory candidate matches 84. Fractions below use these matched subsets, not the entire experimental lipid panel.\n\n')+table(wide)+'\n\n'+(
        'Natri strata overlap donors and must not be summed. UPenn pseudobulk zeros reflect insufficient IPF replication for testing, not evidence of absence. Metacell tests are nominal and include its one/two-donor IPF arms. This union over cell states is exploratory.\n\n'
        '[Every recapitulated lipid, supporting cell types, signed folds and p-values, plus minimum covering combinations](COHORT_CONCORDANT_COVERAGE.md).\n\n')+marker+'\n\n'
    path=out/'EVALUATION_SUMMARY.md'
    s=path.read_text()
    if marker in s:
        before,_,rest=s.partition(marker);_,_,after=rest.partition(marker);s=before+after.lstrip('\n')
    first,rest=s.split('\n',1)
    if (out/'ADAMS_NATRI_RESULTS_AND_METHODS.md').exists():
        block += '[Adams, Natri and their intersection: Results and associated Methods drafts](ADAMS_NATRI_RESULTS_AND_METHODS.md).\n\n'
    path.write_text(first+'\n\n'+block+rest.lstrip('\n'))


def render_summary(out):
    d=pd.read_csv(out/'both_significant_fold_threshold_summary.csv')
    selected=d[d.panel.eq('shared_MS1_47') & d.model.isin(['candidate_MS1_47','current_as_deployed']) &
               d.minimum_imputed_fold.eq(1) & d.overlap_gate.isin(['both_rawp005','both_BH010']) &
               (((d.contrast=='Adams2020__IPF_vs_Healthy') & d.population.isin(['tMDM','iMON'])) |
                ((d.contrast=='Natri2024__IPF_vs_Healthy') & d.population.isin(['iMON','Ciliated-Bronch'])))].copy()
    cols=['contrast','unit','population','model','overlap_gate','concordant_up','concordant_down','discordant_total']
    candidates=d[d.panel.eq('all_matched') & d.model.str.startswith('candidate') & d.overlap_gate.eq('both_BH010') &
                 d.n_case_donors.ge(3) & d.n_control_donors.ge(3)]
    maxima=candidates.groupby(['model','unit','minimum_imputed_fold']).concordant_total.max().rename('best_single_state_concordant_count').reset_index()
    lines=['# LungMAP IPF evaluation: results at a glance','',
           'Completed imputation of all 26,639 real sample/state pseudobulks and 230,057 metacells with separate 47-lipid MS1 and 219-lipid exploratory candidates. All seven saved IPF/control contrasts were evaluated with the frozen cellHarmony parameters and no metacell cap. Existing models and deployed predictions/differentials were not replaced.','',
           '## Significant in both measurements','',
           'The examples below use the identical current/MS1 experimental feature panel. Both-raw means experimental and predicted raw p < 0.05; both-BH means experimental and predicted BH < 0.10. Up/down counts are concordant with the experimental direction. Discordant results are also significant in both. No single-measurement significance is counted. Each model retains its full-panel BH denominator.','',
           table(selected[cols]),'',
           '## Fold threshold sensitivity','',
           'These maxima are the best single dataset/cell-state concordant counts at each threshold, among comparisons with at least three actual donors in each arm. Different rows may select different cell states; they are exploratory rankings, not independent validation. Metacell tests use metacells rather than donor-level replicates.','',
           table(maxima),'',
           'The strong Adams tMDM and Natri iMON pseudobulk significance overlaps both fall to zero at a 1.1-fold threshold. Thus small predicted effect sizes remain a limitation despite the corrected abundance scale. Neither candidate has a both-BH concordant cell-state result at 1.5-fold among the donor-supported comparisons. Experimental fold magnitudes cannot be filtered reliably until the source transform/log base is confirmed.','',
           '## Example of hypothesized cell-state-specific effects','',
           'In Natri pseudobulks, Cer(d18:0/24:0)_NEG agrees with the tissue increase in iMON (+1.065 fold; raw p=0.000208; BH=0.00489) and bronchiolar multiciliated cells (+1.032; p=0.00907; BH=0.0261), while PBFB is opposite (−1.118; p=0.000908; BH=0.00834) and basal cells are opposite (−1.137; p=0.00660; BH=0.0517). These are predictions of opposing state effects. The tissue experiment does not validate their cellular origin.','',
           '## Explanation and independent evidence','',
           'Exact RNA gene decompositions explain 82 selected lipid/state predictions, with maximum reconstruction error 7.5e-15. They explain geometric fold changes; differential tables separately report ratios of arithmetic abundance means. Gene coefficients are predictive associations rather than causal mechanisms. Dataset/state/fold selection using this evidence remains exploratory. Natri fibrosis strata share donors and do not provide independent cohort replications.','',
           'Whole-lung RNA summed within independent donors is separately tested. On the shared MS1 panel, whole-lung both-BH concordant counts for current → candidate are Adams 5 → 3, Jaiswal 0 → 0, Natri less-fibrotic 4 → 6, more-fibrotic 7 → 11, and combined 4 → 10. UPenn has only two/one IPF donors per technology stratum and contributes fold-only whole-lung results.','',
           'Saved pseudobulk differentials reproduce to numerical precision. All saved metacell raw-p calls are retained, but p-value residuals reach 0.00049 and the rerun adds 18 calls. This limit is audited; metacell reruns are not bit-identical.','',
           '## Detailed artifacts','',
           '- [Both-significant overlap, every threshold, and opposing-state hypotheses](BOTH_SIGNIFICANT_OVERLAP_REPORT.md)',
           '- [Dataset/cell-type rankings and paired current comparisons](DATASET_CELLTYPE_SELECTION.md)',
           '- [Exact RNA gene explanations](RNA_GENE_EXPLANATIONS.md)',
           '- [Full lipid tables, methods, donor-level validation and reproduction audit](LungMAP_IPF_CANDIDATE_REPORT.md)','']
    (out/'EVALUATION_SUMMARY.md').write_text('\n'.join(lines))
    if (out/'cohort_concordant_union_summary.csv').exists():
        render_concordant_header(out)
    if (out/'best_concordant_celltypes_by_cohort.csv').exists():
        from summarize_lungmap_best_celltypes import add_best_sections
        add_best_sections(out)


def tissue_tables(out):
    root=out/'whole_lung_donor_validation'
    summary=pd.read_csv(root/'directional_validation_summary.csv')
    matches=pd.read_csv(root/'measured_vs_imputed_all_tests.csv')
    lines=['## Donor-level whole-lung RNA validation', '',
           'This additional comparison matches the tissue experiment’s scope more closely than selecting whichever cell state agrees. Raw counts are summed across captured lung states and libraries, then within each donor in each contrast arm, before RNA normalization and imputation. It uses the same candidate fits and healthy-study offsets, and the unchanged production bundle as comparator. All methods use the frozen cellHarmony moderated t-test, minimum three independent donors per arm and no fold cutoff. These RNA sums represent captured cells, not a measured whole-tissue composition. RNA and experimental lipid donors remain unmatched.', '',
           'Natri’s full IPF contrast has 28 unique donors in 48 samples. More- and less-fibrotic strata overlap in donors and are not separate cohort replications. UPenn’s technology strata have two and one IPF donors, respectively; only folds are reported for these contrasts.', '']
    for panel in ['all_matched','shared_MS1_47','shared_bulk_219']:
        selected=summary[summary.source_gate.eq('supplied_BH010') & summary.panel.eq(panel)]
        lines += [f'### Whole lung: {panel}', '',
                  table(selected.pivot(index=['contrast','model','matched_experimental_differentials','n_case_donors','n_control_donors','status'],columns='prediction_gate',values='concordant_total').reset_index()), '']
    selected=summary[summary.source_gate.eq('supplied_BH010') & summary.panel.eq('all_matched') & summary.prediction_gate.eq('BH010')]
    lines += ['### Whole-lung BH < 0.10: concordant up/down and discordant lipids', '',
              table(selected[['contrast','model','concordant_up','concordant_down','concordant_total','discordant_total','matched_experimental_differentials']]), '',
              '### All whole-lung matched lipid results', '']
    for (contrast,label),g in matches.groupby(['contrast','model']):
        lines += [f'#### {contrast} / {label}', '',
                  table(g[['gene','measured_feature','signed_fold','pval','fdr','supplied_raw_p','supplied_adjusted_p','direction_concordant_supplied','n_case','n_control','status']]), '']
    return lines


def study_support(out, matched, census):
    c=census.drop_duplicates(['contrast','unit','model','population'])
    data=matched.merge(c[['contrast','unit','model','population','n_case_donors','n_control_donors']],on=['contrast','unit','model','population'])
    data=data[~data.population.str.contains('__vs__') & data.supplied_adjusted_p.lt(.1) & data.n_case_donors.ge(3) & data.n_control_donors.ge(3)].copy()
    data['study']=data.contrast.str.split('__').str[0]
    results,details=[],[]
    for (unit,label),g in data.groupby(['unit','model']):
        for gate,field,cut in [('rawp005','pval',.05),('BH005','fdr',.05),('BH010','fdr',.1)]:
            selected=g[g[field].lt(cut)]
            support=selected[selected.direction_concordant_supplied]
            opposing=selected[~selected.direction_concordant_supplied]
            multiple=support.groupby('measured_feature').study.nunique()
            supported_pairs=[];consistent_pairs=[]
            for (feature,state),p in support.groupby(['measured_feature','population']):
                studies=set(p.study)
                against=set(opposing.loc[opposing.measured_feature.eq(feature)&opposing.population.eq(state),'study'])
                if len(studies)>=2:
                    supported_pairs.append(feature)
                    if not against:consistent_pairs.append(feature)
                    details.append({'unit':unit,'model':label,'prediction_gate':gate,'measured_feature':feature,
                                    'population':state,'supporting_studies':';'.join(sorted(studies)),
                                    'opposing_studies':';'.join(sorted(against))})
            results.append({'unit':unit,'model':label,'prediction_gate':gate,
                            'same_lipid_two_studies_any_state':int(multiple.ge(2).sum()),
                            'same_lipid_same_state_two_studies':len(set(supported_pairs)),
                            'same_lipid_same_state_two_studies_no_opposing_study':len(set(consistent_pairs))})
    summary=pd.DataFrame(results)
    summary.to_csv(out/'independent_study_directional_support.csv',index=False)
    pd.DataFrame(details).to_csv(out/'independent_study_directional_support_lipids.csv',index=False)
    return ['## Support across independent studies', '',
            'This descriptive agreement check counts Natri as one study despite its three contrasts. Populations require at least three actual donors in each arm; UPenn’s single/two-case strata cannot contribute. No meta-analysis p-value or correction across states is implied. A metacell significance test remains nominal even when its contributing donor count is adequate.', '',
            table(summary), '',table(pd.DataFrame(details)), '']


def export_folds(out, registry):
    from evaluate_lungmap_ipf import arithmetic_stats, pairs_for, ATLAS
    cross = pd.read_csv(out / 'provenance/crosswalk.tsv', sep='\t')
    mapping = dict(zip(cross.cell_state, cross.cell_type)) | dict(zip(cross.cell_type_full, cross.cell_type)) | dict(zip(cross.cell_type,cross.cell_type))
    all_folds, representative = [], []
    for label in LABELS:
        for unit in ('PB','MC'):
            path = (out / 'current_inputs' / f'cellref2_v8_{"pseudobulk" if unit == "PB" else "metacell"}_lipid_forDE.h5ad'
                    if label.startswith('current') else out / label / f'LungMAP_{unit}_predictions_log2.h5ad')
            a = ad.read_h5ad(path, backed='r')
            obs = a.obs.copy()
            obs['cell_state'] = obs.get('cell_state', obs.get('short_name')).map(mapping)
            empty = pd.DataFrame(columns=['cell2_disease_state','cell1_reference_state'])
            for row in registry.itertuples():
                _, counts = pairs_for(obs,row,empty,3 if unit=='PB' else 5)
                uid=obs.Study_internal.astype(str)+'|'+obs.Sample.astype(str)
                def keys(value):return {'|'.join(s.split('|')[-2:]) for s in str(value).split(';') if s and s!='nan'}
                ca,co=uid.isin(keys(row.case_keys)),uid.isin(keys(row.control_keys))
                for count in counts.itertuples():
                    mask=obs.cell_state.eq(count.population)
                    ncase,nctrl=int((ca&mask).sum()),int((co&mask).sum())
                    if min(ncase,nctrl)==0:continue
                    if count.status=='tested' and count.population!='AT2':continue
                    inds=np.r_[np.flatnonzero(ca&mask),np.flatnonzero(co&mask)]
                    vals=a[inds].to_memory().X
                    vals=vals.toarray() if hasattr(vals,'toarray') else np.asarray(vals)
                    stats=arithmetic_stats(vals,ncase,a.var_names)
                    stats['signed_fold']=np.where(stats.log2fc>0,np.exp2(stats.log2fc),-np.exp2(-stats.log2fc))
                    stats.loc[stats.log2fc.eq(0),'signed_fold']=0
                    stats['case_mean_relative_linear']=np.exp2(stats.case_mean_expr)
                    stats['control_mean_relative_linear']=np.exp2(stats.control_mean_expr)
                    stats['contrast'],stats['unit'],stats['model']=row.contrast_id,unit,label
                    stats['population'],stats['n_case'],stats['n_control']=count.population,ncase,nctrl
                    stats['n_case_donors'],stats['n_control_donors']=count.n_case_donors,count.n_control_donors
                    stats.index.name='gene'
                    if count.status!='tested':
                        stats['reason']='Below saved replication gate; folds only, no p or FDR'
                        all_folds.append(stats.reset_index())
                    if count.population=='AT2':representative.append(stats.reset_index())
            a.file.close()
    folds=pd.concat(all_folds,ignore_index=True) if all_folds else pd.DataFrame()
    folds.to_csv(out/'underpowered_fold_only_results.csv',index=False)
    rep=pd.concat(representative,ignore_index=True)
    rep.to_csv(out/'representative_AT2_abundance_values.csv',index=False)
    return folds,rep


def render(out):
    audit=json.loads((out/'training_and_input_audit.json').read_text())
    summary=pd.read_csv(out/'directional_validation_summary.csv')
    measured=pd.read_csv(out/'experimental_lung_IPF_lipids_all_544.csv',index_col=0)
    matched=pd.read_csv(out/'all_measured_vs_imputed_tests.csv')
    mapping=pd.read_csv(out/'lipid_matching_audit.csv')
    registry=pd.read_csv(out/'provenance/IPF_vs_healthy_registry.tsv',sep='\t')
    census=pd.read_csv(out/'population_replication_census.csv')
    parity=pd.read_csv(out/'current_saved_vs_rerun_parity.csv')
    folds,representative=export_folds(out,registry)
    summary=summary[summary.include_state_overrides.eq(False)].copy()
    source_counts=[]
    for name,mask,eff in [('supplied raw p < 0.05',measured.supplied_raw_p.lt(.05),measured.supplied_log_effect),('supplied BH < 0.10',measured.supplied_adjusted_p.lt(.1),measured.supplied_log_effect),('donor raw p < 0.05',measured.pvalue.lt(.05),measured.geometric_log2FC),('donor BH < 0.10',measured.FDR.lt(.1),measured.geometric_log2FC)]:
        source_counts.append({'source_gate':name,'up':int(((eff>0)&mask).sum()),'down':int(((eff<0)&mask).sum()),'total':int(mask.sum())})
    lines=['# LungMAP IPF lipid candidate evaluation', '',
           'Evaluation candidates only. Production models, LungMAP prediction matrices, deployed differentials and databases were not replaced.', '',
           '## Inputs and methods', '',
           f'The full Transfer source contains {audit["PB_raw_library_rows"]:,} real library/state pseudobulks, summed before RNA normalization into {audit["PB_sample_state_rows"]:,} sample/state pseudobulks, and {audit["MC_rows"]:,} per-sample metacells. All rows were imputed. IPF/control testing uses all seven saved registry contrasts, with no metacell cap.', '',
           'Experimental validation uses Dr. Clair’s lung-tissue `10_results_with_statistics.csv`: 544 lipid features, 30 IPF tissue measurements and 10 control measurements. L-prefix grouping provisionally identifies 10 IPF donors with three tissues each and 10 control donors. Provided IPF-versus-control raw and adjusted p-values are primary source evidence. Recalculated donor-level Welch/BH statistics are separately reported. Source log base and stage metadata remain unconfirmed; source-effect signs do not require assuming the log base. Matching requires unambiguous acyl composition and ion mode. For current outputs lacking mode suffixes, the mode is recovered only from an exact single annotated training target.', '',
           'RNA inputs are log2(1 + counts per 10,000), with the denominator computed across the full transcriptome. Candidates use healthy donor-balanced whole-sample RNA to compute one additive gene offset per study, applied identically to both disease arms. Studies without normals use all healthy reference donors. IPF lipid outcomes are excluded from training, RNA calibration and alpha selection.', '',
           table(pd.DataFrame([dict(candidate=k,**v) for k,v in audit['training'].items()])), '',
           'The candidates retain the repaired target reference and fixed alpha=100. GPX1 and JMJD7-PLA2G4B are omitted because the LungMAP feature reference lacks both; models are refit on 1,293 shared genes. Target values and measured reference differentials remain unchanged. Candidate MS1 values are relative instrument-signal estimates, not molar concentrations. Additional 219-panel lipids use the earlier relative bulk gauge and lack the 47-panel raw-MS support.', '',
           'The production comparator has 202 lipid-wise ElasticNet outputs trained on sorted references. The candidates use linear kernel ridge with combined sorted/bulk references and corrected targets. Changes in performance cannot be attributed solely to output scaling; reference composition, model architecture and the shared RNA panel also change. The explicit-log2 current comparator isolates the differential routing correction on unchanged predictions.', '',
           'The frozen deployed cellHarmony functions supply both statistical tests: pseudobulks use its empirical-Bayes moderated t-test with shrinkage=0.2 and at least three rows per arm; metacells use Scanpy Wilcoxon with tie_correct=False and at least five metacells per arm. No fold cutoff is applied (fc=1). Raw p < 0.05, BH < 0.05 and BH < 0.10 are tabulated. BH is calculated within each population/contrast across each model’s full independently filtered panel; the panel sizes therefore differ. No correction across all cell states or contrasts is implied.', '',
           'Current-as-deployed preserves the saved heuristic scale handling. Current-explicit-log2 uses the same stored predictions with explicit log metadata. Candidate tests use explicit log2 abundance, preventing another normalization of logged values above 20. Candidate and explicit-current folds are ratios of arithmetic means after exact 2**prediction inversion, with no abundance pseudocount. Signed folds are positive ratios for increases and negative reciprocal ratios for decreases.', '',
           'Metacell p-values test metacells, not independent donors. Sample pseudobulks can include more than one anatomical sample from a donor, especially in Natri. The replication census includes actual donor counts, and the report does not label these nominal p-values as independent-donor replication. Comparisons against Dr. Clair’s whole-lung experiment assess directional compatibility of cell-state predictions, not measured cell-state validation.', '',
           'State-to-reference override contrasts are computed and available in the CSV summary but excluded from the primary tables below. A lipid can agree in one state and disagree in another. Unique agreement and disagreement counts can therefore overlap; counts of lipid/state results are separate.', '',
           '## Experimental differential counts', '',
           table(pd.DataFrame(source_counts)), '',
           '## Primary per-contrast counts', '',
           'Each numerator below is a distinct experimentally BH-significant lipid with the same direction and the indicated imputation significance in at least one ordinary cell state. The denominator is the experimentally significant subset matched to that model. “Opposing” means significant in the opposite direction in at least one state; it is not an exclusive category.', '']
    lines += tissue_tables(out)
    primary=summary[summary.source_gate.eq('supplied_BH010') & summary.panel.eq('all_matched')]
    for unit in ('PB','MC'):
        lines += [f'### {unit}: {"real sample pseudobulks" if unit=="PB" else "all metacells"}', '']
        selected=primary[primary.unit.eq(unit)]
        pivot=selected.pivot(index=['contrast','model','eligible_matched_lipids','matched_tested_lipids'],columns='prediction_gate',values='concordant_unique_lipids').reset_index()
        lines += [table(pivot), '']
        detail=selected[selected.prediction_gate.eq('BH010')][['contrast','model','concordant_up_lipids','concordant_down_lipids','discordant_unique_lipids','both_direction_lipids','concordant_lipid_state_results','discordant_lipid_state_results']]
        lines += ['BH < 0.10 directions and repeated state results:', '', table(detail), '']
    lines += ['## Comparisons on identical lipid panels', '',
              'Both methods are restricted to the same matched experimental features here. Each retains its original within-model BH denominator, so these are coverage-matched comparisons, not identical multiple-testing families.', '']
    for panel in ('shared_MS1_47','shared_bulk_219'):
        selected=summary[summary.source_gate.eq('supplied_BH010') & summary.panel.eq(panel)]
        lines += [f'### {panel}', '',table(selected.pivot(index=['contrast','unit','model','eligible_matched_lipids','matched_tested_lipids'],columns='prediction_gate',values='concordant_unique_lipids').reset_index()), '']
    lines += study_support(out,matched,census)
    lines += ['## Selecting datasets and cell types and explaining predictions', '',
              'The complete exploratory rankings are in [DATASET_CELLTYPE_SELECTION.md](DATASET_CELLTYPE_SELECTION.md). They rank each dataset/state by concordant minus discordant lipid counts, with separate raw-p and BH thresholds, donor counts, and identical-panel current/candidate comparisons. Selection using these outcomes is exploratory. The [both-significant overlap report](BOTH_SIGNIFICANT_OVERLAP_REPORT.md) evaluates intersections at raw/BH significance and five imputed fold thresholds, recording potential opposing cell-state effects. The strongest selected pseudobulk predictions also have exact RNA gene contribution decompositions in [RNA_GENE_EXPLANATIONS.md](RNA_GENE_EXPLANATIONS.md). These explain model inputs rather than asserting causal mechanisms.', '']
    lines += ['## Recalculated donor-level experimental evidence', '',
              'Same analysis using experimental donor-level BH < 0.10 and geometric-effect directions. This avoids relying only on the supplied tissue-level statistics.', '',
              table(summary[summary.source_gate.eq('donor_BH010') & summary.panel.eq('all_matched')].pivot(index=['contrast','unit','model','eligible_matched_lipids','matched_tested_lipids'],columns='prediction_gate',values='concordant_unique_lipids').reset_index()), '',
              '## Saved-current reproduction audit', '', table(parity), '',
              'Pseudobulk saved p-values reproduce within 1.6e-16. All saved metacell raw-p calls are retained, but metacell p-values have residual differences up to 0.00049 and the rerun adds 18 calls across seven contrasts. Saved fold differences are below 3e-8 and replication counts agree. The same frozen parameters and original Scanpy environment are used; metacell reproduction is not bit-identical. The saved-versus-rerun call sets are separately audited in current_saved_call_set_audit.csv.', '',
              '## Replication census', '',table(census[census.model.eq('candidate_MS1_47') & census.status.eq('tested')][['contrast','unit','population','n_case','n_control','n_case_donors','n_control_donors','override']]), '',
              '## Underpowered comparisons: folds only', '',
              'Populations with measurements in both arms but fewer than the saved replication minimum have no reported p-value or FDR. They cannot enter the significance-based validation totals.', '']
    if len(folds):
        lines += [table(folds[['contrast','unit','model','population','gene','n_case','n_control','signed_fold']]), '']
    else:lines += ['None with both arms measured.', '']
    lines += ['## Representative abundance estimates', '',
              'AT2, Adams IPF/control: arithmetic arm means in each model’s relative linear output space. Current values retain their original abundance gauge; candidate MS1 values use the repaired instrument-signal gauge. Sizes of values are not a cross-lipid concentration comparison.', '']
    rep=representative[representative.contrast.eq('Adams2020__IPF_vs_Healthy') & representative.unit.eq('PB')]
    rep=rep.merge(mapping[mapping.status.eq('matched')][['model','gene','measured_feature']],on=['model','gene'])
    chosen=['Cer(d18:1/16:0)_NEG','Cer(d18:1/24:0)_NEG','PE(18:0/22:6)_NEG','PI(18:0/18:1)_NEG']
    lines += [table(rep[rep.measured_feature.isin(chosen)][['model','gene','measured_feature','case_mean_relative_linear','control_mean_relative_linear','signed_fold']]), '',
              '## All experimental lipids', '',
              'Native transformed donor means are shown in their supplied measurement space. Experimental fold magnitudes are not relabeled as confirmed log2 ratios. All provided and recalculated raw/BH values are retained.', '']
    donor=pd.read_csv(HERE/'artifacts/IPF_candidate_validation_20261001/measured_lung_donor_native.csv',index_col=0)
    md=pd.read_csv(HERE/'artifacts/IPF_candidate_validation_20261001/measured_lung_sample_metadata.csv',index_col=0)
    group=md.groupby('donor').group.first()
    full=measured.copy();full['control_native_mean']=donor.loc[group.eq('NDC')].mean();full['IPF_native_mean']=donor.loc[group.eq('IPF')].mean()
    full['supplied_direction']=np.where(full.supplied_log_effect.gt(0),'up',np.where(full.supplied_log_effect.lt(0),'down','zero'))
    cols=['original_annotation','ion_mode','control_native_mean','IPF_native_mean','supplied_direction','supplied_raw_p','supplied_adjusted_p','pvalue','FDR']
    lines += [table(full.rename_axis('measured_feature').reset_index()[['measured_feature']+cols]), '',
              '## Every matched lipid: directional support by contrast and cell state', '',
              'For each matched feature, report the full model panel’s nominal p/BH tests. Positive folds indicate increases and negative folds decreases. All unmatched targets and ambiguous identities remain in lipid_matching_audit.csv.', '']
    for (contrast,unit,label),g in matched.groupby(['contrast','unit','model']):
        lines += [f'### {contrast} / {unit} / {label}', '',
                  table(g[['population','gene','measured_feature','signed_fold','pval','fdr','supplied_raw_p','supplied_adjusted_p','direction_concordant_supplied']]), '']
    lines += ['## Files and reproduction', '',
              'Complete all-feature tests: all_cellHarmony_tests.csv. Validation matches: all_measured_vs_imputed_tests.csv. Threshold and identical-panel counts: directional_validation_summary.csv. Candidate full pseudobulk/metacell predictions are stored separately as log2 and positive linear H5ADs under each candidate directory. Source matching, donor counts, frozen engine and model fits are retained.', '',
              '```bash', '/opt/homebrew/bin/python3.11 components/rna2lipid/evaluate_lungmap_ipf.py', '/opt/homebrew/bin/python3.11 components/rna2lipid/evaluate_lungmap_whole_lung.py', '/opt/homebrew/bin/python3.11 components/rna2lipid/render_lungmap_ipf_evaluation.py', '```', '']
    (out/'LungMAP_IPF_CANDIDATE_REPORT.md').write_text('\n'.join(lines))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default=str(HERE/'artifacts/LungMAP_IPF_candidate_20261004'))
    output=Path(p.parse_args().out)
    render(output)
    if (output/'both_significant_fold_threshold_summary.csv').exists():
        render_summary(output)
