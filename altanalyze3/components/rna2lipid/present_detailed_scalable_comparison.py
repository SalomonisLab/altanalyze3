"""Present existing same-state model comparisons with cohort and lipid detail."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BASE = HERE / 'artifacts/LungMAP_full202_BH_corrected_20261006'
OUT = BASE / 'scALABLE_model_comparison'
CURRENT = 'current_explicit_log2'
CANDIDATE = 'candidate_native_log2_202'
KEYS = ['model', 'contrast', 'unit', 'population', 'gate', 'fold_rule']
GATES = ['rawp005', 'BH005', 'BH010']
FOLDS = ['No cutoff', '>1.1', '>1.2', '>1.5']


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source = OUT / 'all_model_threshold_results.csv'
    lipid_source = OUT / 'all_significant_overlap_lipid_rows.csv'
    hashes = {str(p): sha(p) for p in [source, lipid_source]}
    audit = json.loads((OUT / 'comparison_audit.json').read_text())
    for path, digest in audit['source_sha256'].items():
        assert sha(Path(path)) == digest
    a = pd.read_csv(source, float_precision='round_trip')
    assert len(a) == 13872 and not a.duplicated(KEYS).any()
    assert set(a.model) == {CURRENT, CANDIDATE}
    assert set(a.gate) == set(GATES) and set(a.fold_rule) == set(FOLDS)
    census = pd.read_csv(OUT / 'complete_replication_census.csv')
    primary = a[a.comparison_kind.eq('same_state')].copy()
    eligible = primary[primary.above75 & primary.n_case_donors.ge(2) & primary.n_control_donors.ge(2)]
    ranked = eligible.sort_values(['concordant', 'agreement_percent', 'opposite', 'population'], ascending=[False, False, True, True])
    global_best = ranked.drop_duplicates(['model', 'gate', 'fold_rule']).copy()
    assert len(global_best) == 24
    global_best['gate_order'] = global_best.gate.map({g:i for i,g in enumerate(GATES)})
    global_best['fold_order'] = global_best.fold_rule.map({f:i for i,f in enumerate(FOLDS)})
    global_best = global_best.sort_values(['model', 'gate_order', 'fold_order'])
    per_cohort = ranked.drop_duplicates(['model', 'contrast', 'unit', 'gate', 'fold_rule'])
    # A complete report grid also preserves groups without a qualifying state.
    roster = sorted(set(census[['contrast','unit','model']].itertuples(index=False,name=None)))
    lookup = per_cohort.set_index(['contrast','unit','model','gate','fold_rule'])
    complete = []
    for contrast,unit,model in roster:
        tested = primary[(primary.contrast==contrast)&(primary.unit==unit)&(primary.model==model)]
        for gate in GATES:
            for fold in FOLDS:
                key = (contrast,unit,model,gate,fold)
                if key in lookup.index:
                    r = lookup.loc[key].to_dict()
                    r.update(contrast=contrast,unit=unit,model=model,gate=gate,fold_rule=fold,status='above75_multiple_donor_arms')
                else:
                    r = {'contrast':contrast,'unit':unit,'model':model,'gate':gate,'fold_rule':fold,
                         'population':'—','cell_state_full':'—',
                         'status':'not_tested_insufficient_replicates' if tested.empty else 'no_above75_state_with_multiple_donors'}
                complete.append(r)
    complete = pd.DataFrame(complete)
    assert len(complete)==len(roster)*12
    ix = a.set_index(KEYS)
    # Display both directions of selection, then compare the identical state.
    pairrows = []
    for r in global_best.to_dict('records'):
        cur = ix.loc[(CURRENT,r['contrast'],r['unit'],r['population'],r['gate'],r['fold_rule'])]
        cand = ix.loc[(CANDIDATE,r['contrast'],r['unit'],r['population'],r['gate'],r['fold_rule'])]
        record = {'selected_by_model':r['model'],'contrast':r['contrast'],'unit':r['unit'],
                  'population':r['population'],'cell_state_full':r['cell_state_full'],'gate':r['gate'],'fold_rule':r['fold_rule']}
        for label,g in [('current',cur),('candidate',cand)]:
            for col in ['concordant_up','concordant_down','concordant','opposite','directional_denominator',
                        'agreement_percent','n_case_donors','n_control_donors','concordant_lipids','opposite_lipids']:
                record[label+'_'+col]=g[col]
        pairrows.append(record)
    paired = pd.DataFrame(pairrows)
    # Verify all displayed up/down/opposite counts against actual lipid rows.
    lipids = pd.read_csv(lipid_source, float_precision='round_trip')
    grouped = {key:g for key,g in lipids.groupby(KEYS)}
    display_keys = set(global_best[KEYS].itertuples(index=False,name=None))
    for r in paired.itertuples():
        for model in [CURRENT,CANDIDATE]:
            display_keys.add((model,r.contrast,r.unit,r.population,r.gate,r.fold_rule))
    verified=0
    for key in display_keys:
        g = grouped.get(key,lipids.iloc[:0]); r = ix.loc[key]
        c = g[g.direction_relation.eq('concordant')];o = g[g.direction_relation.eq('opposite')]
        assert len(c)==r.concordant and len(o)==r.opposite
        assert int(c.supplied_log_effect.gt(0).sum())==r.concordant_up
        assert int(c.supplied_log_effect.lt(0).sum())==r.concordant_down
        assert r.concordant_up+r.concordant_down==r.concordant
        verified+=1
    global_best.drop(columns=['gate_order','fold_order']).to_csv(OUT/'detailed_global_best_states.csv',index=False)
    complete.to_csv(OUT/'detailed_all_cohort_threshold_tables.csv',index=False)
    paired.to_csv(OUT/'detailed_identical_state_model_comparisons.csv',index=False)
    fields=['gate','fold_rule','contrast','unit','population','cell_state_full','concordant_up','concordant_down',
            'concordant','opposite','agreement_percent','n_case_donors','n_control_donors']
    def md(frame):return frame.to_markdown(index=False,floatfmt='.5g')
    lines=['# Detailed IPF comparison: supported scALABLE model and corrected candidate','',
           'All comparisons are IPF versus healthy within the same RNA cohort and annotated cell state. '
           'The tissue lipidomics subjects are unmatched to the RNA cohorts. The p cutoff is required in both '
           'the experimental lipidomics and that model’s predicted differential results. '
           'Raw p<0.05, BH<0.05 and BH<0.1 are displayed separately, at no fold cutoff or strict predicted fold >1.1, >1.2 and >1.5.', '',
           'Concordant up/down = lipids significant in both with matching increase/decrease. Opposite = significant in both '
           'with opposing directions. Agreement = concordant/(concordant+opposite). Best means the largest concordant count '
           'among states exceeding 75%, followed by higher agreement, fewer opposite lipids, and state name. '
           'The following display tables require at least two source donors in each arm, consistent with the prior candidate table. '
           'Single-donor results remain in the original full comparison report and datasets; no model feature or sample is removed. '
           'Source donor coverage does not turn metacell tests into donor-level replication tests.', '',
           'The current comparator is the saved supported 202-lipid production ElasticNetCV predictions. '
           'Both models use the same approved full-202 BH family. The running website has not been inspected. '
           'No predictions, model parameters, raw tests or BH values were changed. '
           'Threshold/state selection is descriptive and Natri subgroups share donors.', '',
           '## Current scALABLE: best qualifying state for each setting', '',
           md(global_best[global_best.model.eq(CURRENT)][fields]), '',
           '## Corrected candidate: best qualifying state for each setting', '',
           md(global_best[global_best.model.eq(CANDIDATE)][fields]), '',
           'Independent best states may differ between models. Compare the same state below to avoid interpreting '
           'different cell states as a direct model comparison.', '',
           '## Identical-state comparison at candidate-selected settings', '',
           md(paired[paired.selected_by_model.eq(CANDIDATE)].drop(columns=['current_concordant_lipids','current_opposite_lipids','candidate_concordant_lipids','candidate_opposite_lipids'])), '',
           '## Identical-state comparison at current-model-selected settings', '',
           md(paired[paired.selected_by_model.eq(CURRENT)].drop(columns=['current_concordant_lipids','current_opposite_lipids','candidate_concordant_lipids','candidate_opposite_lipids'])), '',
           '## Complete tables separately for every cohort, input and model', '']
    for key,g in complete.groupby(['contrast','unit','model'],sort=True):
        order = pd.MultiIndex.from_product([GATES,FOLDS],names=['gate','fold_rule'])
        g=g.set_index(['gate','fold_rule']).reindex(order).reset_index()
        lines += ['', '### '+' / '.join(key), '',md(g[['gate','fold_rule','population','cell_state_full','concordant_up',
                 'concordant_down','concordant','opposite','agreement_percent','n_case_donors','n_control_donors','status']])]
    lines += ['', '## Specific lipids for every state used in the main paired tables', '',
              'Signed +2 means twofold higher and -2 means twofold lower. Experimental log-effect base remains '
              'unconfirmed; its sign supplies the reference direction. Supplied p-values retain source rounding.']
    for key in sorted(display_keys):
        g=grouped.get(key,lipids.iloc[:0])
        lines += ['', '### '+' / '.join(key), '',md(g[['measured_feature','gene','signed_fold','pval','fdr',
                   'supplied_raw_p','supplied_adjusted_p','direction_relation']])]
    lines += ['', '## Complete results and downloads', '',
              '- [All cohort/model/threshold tables, including no-qualifying-state entries](detailed_all_cohort_threshold_tables.csv)',
              '- [Both directions of selection, compared at identical states and thresholds](detailed_identical_state_model_comparisons.csv)',
              '- [Original full comparison, including every qualifying state and lipid for both models](SCALABLE_VS_CANDIDATE_REPORT.md)',
              '- [Every state and threshold, including nonqualifying results](all_model_threshold_results.csv)']
    report=OUT/'DETAILED_COHORT_CELL_STATE_COMPARISON.md'
    report.write_text('\n'.join(lines)+'\n')
    assert all(sha(Path(path))==digest for path,digest in hashes.items())
    result={'source_sha256':hashes,'global_best_rows':len(global_best),'identical_state_display_rows':len(paired),
            'complete_cohort_model_threshold_slots':len(complete),'up_down_and_opposite_lipid_counts_verified':verified,
            'sources_unchanged':True,'new_training_inference_tests_or_BH':False,'script_sha256':sha(Path(__file__))}
    (OUT/'detailed_presentation_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print('Report:',report)


if __name__=='__main__':main()
