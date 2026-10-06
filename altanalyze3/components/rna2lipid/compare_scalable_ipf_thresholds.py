"""Compare saved supported production predictions with the corrected candidate."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'artifacts/LungMAP_full202_BH_corrected_20261006'
NATIVE = HERE / 'artifacts/LungMAP_full202_native_log2_20261005'
OUT = SOURCE / 'scALABLE_model_comparison'
KEYS = ['contrast', 'unit', 'model', 'population']
CANDIDATE = 'candidate_native_log2_202'
CURRENT = 'current_explicit_log2'
FOLDS = [('No cutoff', None), ('>1.1', 1.1), ('>1.2', 1.2), ('>1.5', 1.5)]
GATES = [('rawp005', 'supplied_raw_p', 'pval', .05),
         ('BH005', 'supplied_adjusted_p', 'fdr', .05),
         ('BH010', 'supplied_adjusted_p', 'fdr', .1)]
MACROPHAGES = {'AM', 'AM-lipid', 'AM-prolif', 'MT+ AM', 'tMDM', 'IM'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    paths = [SOURCE / n for n in ['all_measured_vs_imputed_tests.csv',
             'all_cellHarmony_tests.csv', 'population_replication_census.csv', 'lipid_matching_audit.csv']]
    hashes = {str(p): sha(p) for p in paths}
    contract = json.loads((HERE / 'integrity/baseline_contract.json').read_text())
    bundle = Path(contract['production_bundle']['path'])
    assert sha(bundle) == contract['production_bundle']['sha256']
    bindings = list((NATIVE / 'differentials').glob('*/' + '*/' + CURRENT + '/source_binding.json'))
    assert len(bindings) == 14
    assert all(json.loads(p.read_text())['bundle'] == sha(bundle) for p in bindings)
    copied_inputs = HERE / 'artifacts/LungMAP_IPF_candidate_20261004/current_inputs'
    current_input_hashes = {str(p): sha(p) for p in [
        copied_inputs / 'cellref2_v8_pseudobulk_lipid_forDE.h5ad',
        copied_inputs / 'cellref2_v8_metacell_lipid_forDE.h5ad']}
    d = pd.read_csv(paths[0], float_precision='round_trip')
    full = pd.read_csv(paths[1], float_precision='round_trip')
    census = pd.read_csv(paths[2]); mapping = pd.read_csv(paths[3])
    expected = set(contract['lipids'])
    assert len(expected) == 202 and set(full.model) == set(d.model) == {CANDIDATE, CURRENT}
    groups = set()
    for key, g in full.groupby(KEYS):
        assert len(g) == 202 and set(g.gene) == expected, key
        assert np.isfinite(g[['pval','fdr','log2fc']].to_numpy()).all()
        assert np.max(np.abs(false_discovery_control(g.pval.to_numpy(), method='bh') - g.fdr.to_numpy())) < 1e-12
        groups.add(key)
    assert not d.duplicated(KEYS + ['measured_feature']).any()
    assert d.BH_panel_size.eq(202).all() and d.BH_policy.eq('full_panel_no_abundance_filter').all()
    assert set(d.groupby(KEYS).groups) == groups
    assert np.isfinite(d[['pval','fdr','signed_fold','log2fc','supplied_raw_p','supplied_adjusted_p','supplied_log_effect']].to_numpy()).all()
    assert np.max(np.abs(np.log2(d.signed_fold.abs()) - d.log2fc.abs())) < 1e-12
    features = {model: set(g.loc[g.status.eq('matched'),'measured_feature']) for model,g in mapping.groupby('model')}
    assert features[CURRENT] == features[CANDIDATE]
    reference = d.drop_duplicates('measured_feature').set_index('measured_feature')
    for col in ['supplied_raw_p','supplied_adjusted_p','supplied_log_effect']:
        assert np.array_equal(d[col].to_numpy(), reference.loc[d.measured_feature,col].to_numpy())
    f = full.set_index(KEYS + ['gene']); ix = d.set_index(KEYS + ['gene'])
    for col in ['pval','fdr','log2fc','signed_fold']:
        assert np.array_equal(ix[col].to_numpy(),f.loc[ix.index,col].to_numpy())
    rep = census[census.status.eq('tested')].set_index(KEYS)
    assert not rep.index.duplicated().any() and set(rep.index) == groups
    names_df = pd.read_csv(NATIVE / 'provenance/crosswalk.tsv',sep='\t')
    names = dict(zip(names_df.cell_type,names_df.cell_state))
    records, detail = [], []
    for key, g in d.groupby(KEYS,sort=True):
        assert set(g.measured_feature) == features[key[2]], key
        for gate, expcol, predcol, cutoff in GATES:
            significance = g[expcol].lt(cutoff) & g[predcol].lt(cutoff)
            for fold_label, fold in FOLDS:
                mask = significance if fold is None else significance & g.signed_fold.abs().gt(fold)
                s = g[mask].copy()
                relation = np.sign(s.log2fc) * np.sign(s.supplied_log_effect)
                c = s[relation.gt(0)]; o = s[relation.lt(0)]; z = s[relation.eq(0)]
                nc,no = len(c),len(o)
                r = dict(zip(KEYS,key)) | {
                    'cell_state_full':names.get(key[3],key[3]),
                    'comparison_kind':'state_override' if '__vs__' in key[3] else 'same_state',
                    'gate':gate,'p_cutoff_both':cutoff,'fold_rule':fold_label,'fold_cutoff':fold,
                    'concordant':nc,'opposite':no,'zero_or_ambiguous':len(z),'both_significant_after_fold':len(s),
                    'concordant_up':int(c.supplied_log_effect.gt(0).sum()),'concordant_down':int(c.supplied_log_effect.lt(0).sum()),
                    'agreement_percent':100*nc/(nc+no) if nc+no else np.nan,
                    'directional_denominator':nc+no,'above75':4*nc > 3*(nc+no),
                    'concordant_lipids':'; '.join(sorted(c.measured_feature)),
                    'opposite_lipids':'; '.join(sorted(o.measured_feature))}
                for col in ['n_case','n_control','n_case_donors','n_control_donors']:
                    r[col] = rep.loc[key,col]
                assert nc+no+len(z) == len(s)
                records.append(r)
                s['direction_relation'] = np.where(relation.gt(0),'concordant',np.where(relation.lt(0),'opposite','zero_or_ambiguous'))
                s['gate']=gate; s['fold_rule']=fold_label; s['above75_group']=r['above75']
                detail.append(s)
    all_rows = pd.DataFrame(records)
    assert len(all_rows) == len(groups)*12
    pairkeys = ['contrast','unit','population','gate','fold_rule']
    c = all_rows[all_rows.model.eq(CANDIDATE)].set_index(pairkeys)
    p = all_rows[all_rows.model.eq(CURRENT)].set_index(pairkeys)
    assert set(c.index) == set(p.index); p=p.loc[c.index]
    paired = c.copy()
    for col in ['concordant','opposite','agreement_percent','directional_denominator','above75','concordant_lipids','opposite_lipids']:
        paired[col+'_current'] = p[col]
    paired['concordant_difference'] = paired.concordant - paired.concordant_current
    paired['agreement_percentage_point_difference'] = paired.agreement_percent - paired.agreement_percent_current
    paired['concordant_count_comparison'] = np.where(paired.concordant_difference.gt(0),'candidate_higher',np.where(paired.concordant_difference.lt(0),'current_higher','equal'))
    paired = paired.reset_index()
    # Preserve all rows and separately display same-state comparisons.
    primary = all_rows[all_rows.comparison_kind.eq('same_state')]
    summaries=[]
    for key,g in primary.groupby(['model','gate','fold_rule']):
        q=g[g.above75]
        summaries.append(dict(zip(['model','gate','fold_rule'],key)) | {
            'tested_states':len(g),'states_with_directional_overlap':int(g.directional_denominator.gt(0).sum()),
            'states_above75':len(q),'above75_overlap_1_to_4':int(q.directional_denominator.lt(5).sum()),
            'above75_overlap_5_to_9':int(q.directional_denominator.between(5,9).sum()),
            'above75_overlap_at_least10':int(q.directional_denominator.ge(10).sum()),
            'concordant_lipid_state_pairs':int(g.concordant.sum()),'opposite_lipid_state_pairs':int(g.opposite.sum())})
    summary=pd.DataFrame(summaries)
    detail=pd.concat(detail,ignore_index=True)
    unions=[]
    primary_detail=detail[~detail.population.str.contains('__vs__')]
    for key,g in primary.groupby(['contrast','unit','model','gate','fold_rule']):
        feature_c=set();feature_o=set()
        for v in g.concordant_lipids:
            if v:feature_c.update(v.split('; '))
        for v in g.opposite_lipids:
            if v:feature_o.update(v.split('; '))
        gate=key[3];expcol='supplied_raw_p' if gate=='rawp005' else 'supplied_adjusted_p'
        cutoff=.1 if gate=='BH010' else .05
        eligible=int(reference[expcol].lt(cutoff).sum())
        unions.append(dict(zip(['contrast','unit','model','gate','fold_rule'],key)) | {
            'experimental_significant_matched_lipids':eligible,'concordant_any_state':len(feature_c),
            'opposite_any_state':len(feature_o),'both_directions_in_different_states':len(feature_c & feature_o),
            'concordant_coverage_percent':100*len(feature_c)/eligible,
            'concordant_lipids':'; '.join(sorted(feature_c))})
    union=pd.DataFrame(unions)
    uk=['contrast','unit','gate','fold_rule']
    cu=union[union.model.eq(CANDIDATE)].set_index(uk);pu=union[union.model.eq(CURRENT)].set_index(uk)
    assert set(cu.index)==set(pu.index);pu=pu.loc[cu.index]
    union_pair=cu.copy();union_pair['concordant_any_state_current']=pu.concordant_any_state
    union_pair['candidate_minus_current']=cu.concordant_any_state-pu.concordant_any_state
    union_pair=union_pair.reset_index()
    winners=all_rows[all_rows.above75].sort_values(['concordant','agreement_percent','population'],ascending=[False,False,True])
    best=winners[winners.comparison_kind.eq('same_state')].drop_duplicates(['contrast','unit','model','gate','fold_rule'])
    macro=winners[winners.comparison_kind.eq('same_state') & winners.population.isin(MACROPHAGES)].drop_duplicates(['contrast','unit','model','gate','fold_rule'])
    prior_scan=pd.read_csv(SOURCE/'requested_thresholds/all_requested_threshold_results.csv',float_precision='round_trip')
    old_keys=KEYS+['gate','fold_rule']
    prior_scan['fold_rule']=prior_scan.predicted_fold_strictly_greater_than.map(lambda x:f'>{x:g}')
    old=prior_scan.set_index(old_keys); new=all_rows[all_rows.fold_rule.ne('No cutoff')].set_index(old_keys)
    assert set(old.index)==set(new.index)
    for col in ['concordant','opposite','zero_or_ambiguous']:
        assert np.array_equal(new[col].to_numpy(),old.loc[new.index,col].to_numpy())
    previous_no_fold=pd.read_csv(SOURCE/'cell_state_FDR_rankings.csv')
    for gate,cut in [('BH005',.05),('BH010',.1)]:
        e=previous_no_fold[previous_no_fold.bulk_BH_cutoff.eq(cut)&previous_no_fold.predicted_BH_cutoff.eq(cut)].set_index(KEYS)
        n=all_rows[all_rows.gate.eq(gate)&all_rows.fold_rule.eq('No cutoff')].set_index(KEYS)
        assert set(e.index)==set(n.index)
        for col in ['concordant','opposite']:assert np.array_equal(n[col].to_numpy(),e.loc[n.index,col].to_numpy())
    assert all(sha(path)==hashes[str(path)] for path in paths)
    OUT.mkdir(exist_ok=True)
    for filename,frame in [('all_model_threshold_results.csv',all_rows),('candidate_vs_current_same_states_and_thresholds.csv',paired),
        ('model_summary_by_threshold.csv',summary),('cohort_concordant_coverage.csv',union),
        ('candidate_vs_current_cohort_coverage.csv',union_pair),('all_above75_comparisons.csv',winners),
        ('best_state_by_cohort_model_threshold.csv',best),('best_macrophage_by_cohort_model_threshold.csv',macro),
        ('all_significant_overlap_lipid_rows.csv',detail),('complete_replication_census.csv',census)]:
        frame.to_csv(OUT/filename,index=False)
    def md(frame):return frame.to_markdown(index=False,floatfmt='.5g')
    s=summary.pivot(index=['gate','fold_rule'],columns='model',values=['states_above75','above75_overlap_at_least10'])
    s.columns=['__'.join(col) for col in s.columns];s=s.reset_index()
    fields=['contrast','unit','model','gate','fold_rule','population','concordant','opposite','agreement_percent','n_case_donors','n_control_donors']
    lines=['# Supported scALABLE prediction model versus corrected candidate: IPF threshold comparison','',
        'Current comparator: saved 202-lipid production bundle `rna2lipid_hs_lung_lipidwise_bundle.pkl`, SHA256 '+sha(bundle)+'. '
        'The saved v8 pseudobulk/metacell `lipid_forDE.h5ad` files provide its predictions. All 14 source bindings '
        'identify this production bundle, and its saved replay audit reproduces predictions to approximately 2e-6. '
        'This is a comparison with saved supported-model predictions; the running website has not been inspected.','',
        'Both models retain all 202 outputs and the same evaluation row roster. The candidate has the authorized 45 training profiles '
        'after D071 removal; production retains its original 50 profiles. Both use lipid-wise ElasticNetCV. '
        'Saved differential results use explicit log2 fold handling and the user-authorized full-202-lipid BH correction for both models. '
        'These are not necessarily the website’s currently displayed filtered BH results. No model, prediction, raw test, or BH value is changed here.','',
        'Twelve settings: no fold cutoff, or strict predicted fold magnitude >1.1, >1.2, >1.5; '
        'raw p<0.05, BH<0.05, or BH<0.1 in both the predicted and experimental results. '
        'Experimental fold magnitudes are not thresholded because their effect log base is unconfirmed. '
        'Directional agreement is C/(C+opposite), with zero effects recorded separately; >75% is strict. '
        'No overlap-count restriction is imposed. Counts are also displayed by overlap size.','',
        'Cell-state comparisons are IPF versus healthy within the same cohort and annotated cell state. '
        'Overrides comparing different states are labeled separately. The RNA and tissue lipidomics subjects are unmatched. '
        'Natri fibrotic subgroups overlap donors. Metacell tests use metacell observations and are nominal relative to donor replication. '
        'Selecting thresholds or states is descriptive, not independent validation.','',
        '## Number of same-state comparisons exceeding 75% agreement','',md(s),'',
        'Each model has 560 tested same-state comparisons at every setting. The at-least-10 columns are a display facet; '
        'all smaller overlaps remain in the full results. More qualifying comparisons does not by itself establish better prediction.','',
        '## Counts and overlap-size distributions for both models','',md(summary),'',
        '## Distinct experimental lipids recapitulated across any state, separately by cohort','',
        'Coverage counts each experimental lipid once if at least one same-state comparison is significant and concordant. '
        'A lipid may also be opposite in another state. This union coverage is not single-state agreement.', '',
        md(union_pair[['contrast','unit','gate','fold_rule','experimental_significant_matched_lipids','concordant_any_state_current','concordant_any_state','candidate_minus_current']]),'',
        '## Best qualifying state per cohort/input/model/threshold','',md(best[fields]),'',
        '## Best qualifying macrophage state per cohort/input/model/threshold','',md(macro[fields]),'',
        '## Identical-state comparison at every reported candidate best setting','']
    b=best[best.model.eq(CANDIDATE)]
    display=paired.merge(b[['contrast','unit','population','gate','fold_rule']],on=pairkeys,how='inner',validate='one_to_one')
    assert len(display)==len(b)
    lines += [md(display[['contrast','unit','population','gate','fold_rule','concordant','opposite','agreement_percent',
        'concordant_current','opposite_current','agreement_percent_current']]),'',
        '## Every lipid in every comparison above 75% for either model','']
    for key,g in detail[detail.above75_group].groupby(KEYS+['gate','fold_rule']):
        lines += ['', '### '+' / '.join(str(x) for x in key),'',md(g[['measured_feature','gene','signed_fold','pval','fdr',
            'supplied_raw_p','supplied_adjusted_p','direction_relation']])]
    lines += ['', '## Complete outputs','',
        '- [Every identical-state/threshold candidate versus current comparison](candidate_vs_current_same_states_and_thresholds.csv)',
        '- [Cohort-specific distinct concordant lipid coverage](candidate_vs_current_cohort_coverage.csv)',
        '- [Every qualifying comparison, including small overlaps](all_above75_comparisons.csv)',
        '- [Every significant-overlap lipid row](all_significant_overlap_lipid_rows.csv)',
        '- [All tested and untested comparison identities](complete_replication_census.csv)']
    (OUT/'SCALABLE_VS_CANDIDATE_REPORT.md').write_text('\n'.join(lines)+'\n')
    audit={'source_sha256':hashes,'production_bundle_sha256':sha(bundle),'current_saved_prediction_sha256':current_input_hashes,
        'production_source_bindings_verified':len(bindings),'full202_groups_verified':len(groups),
        'matched_experimental_features':len(features[CURRENT]),'settings_per_comparison':12,'all_model_threshold_rows':len(all_rows),
        'paired_state_threshold_rows':len(paired),'previous_strict_threshold_counts_reproduced':len(old),
        'previous_no_fold_BH_counts_reproduced':len(groups)*2,'source_hashes_unchanged':all(sha(p)==hashes[str(p)] for p in paths),
        'new_training_inference_raw_tests_or_BH_recalculation':False,'live_website_verified':False,
        'script_sha256':sha(Path(__file__))}
    (OUT/'comparison_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2));print('\nSUMMARY:\n'+s.to_string(index=False))
    print('\nADAMS/NATRI COVERAGE:\n'+union_pair[union_pair.contrast.isin(['Adams2020__IPF_vs_Healthy','Natri2024__IPF_vs_Healthy'])][['contrast','unit','gate','fold_rule','concordant_any_state_current','concordant_any_state','candidate_minus_current']].to_string(index=False))


if __name__=='__main__':main()
