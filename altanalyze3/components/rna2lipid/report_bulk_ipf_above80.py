"""Extend the approved saved bulk-IPF threshold scan without refitting any tests."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

if __package__:
    from .scan_independent_bulk_lipid_thresholds import approved_bases, experimental_magnitude
else:
    from scan_independent_bulk_lipid_thresholds import approved_bases, experimental_magnitude


HERE = Path(__file__).resolve().parent
ROOT = HERE / 'artifacts/bulk_IPF_three_version_comparison_20261008'
PRIOR = ROOT / 'independent_measured_and_imputed_fold_thresholds'
OUT = ROOT / 'measured_1p2_above80_with_regulated_counts'
AUTHORIZATION = ('Regenerate this table with 1.2 for the measured as a possiblity. '
                 'Only show results above 80% and indicate the # regulated in the '
                 'table at that threshold, not just overlapping')


def main():
    state = json.loads((HERE / 'integrity/decision_state.json').read_text())
    question = next(q for q in state['questions']
                    if q['id'] == 'bulk_IPF_experimental_log_effect_base_20261008')
    bases = approved_bases(question)
    verification = json.loads((PRIOR / 'independent_verification.json').read_text())
    assert verification['all_checks_passed']
    audit = json.loads((ROOT / 'full_panel_rerun/comparison_audit.json').read_text())
    model = audit['model_checks']['candidate']
    assert hashlib.sha256(Path(model['path']).read_bytes()).hexdigest() == model['sha256']
    workbook_audit = json.loads((HERE / 'artifacts/Geremy_original_lipid_statistics_20261008/workbook_audit.json').read_text())
    source = workbook_audit['sources']['IPF']
    assert hashlib.sha256(Path(source['path']).read_bytes()).hexdigest() == source['sha256']
    measured = pd.read_excel(source['path'])
    assert len(measured) == 544 and measured.iloc[:, 0].is_unique
    detail_file = PRIOR / 'all_202_original_means_and_statistics.csv'
    details = pd.read_csv(detail_file, float_precision='round_trip')
    details = details[details.version.eq('New corrected model / CPM')].copy()
    contract = json.loads((HERE / 'integrity/baseline_contract.json').read_text())
    assert details.model_lipid.tolist() == contract['lipids'] and len(details) == 202
    assert details.status.eq('matched').sum() == 134
    source_index = measured.set_index(measured.columns[0])
    for r in details[details.status.eq('matched')].itertuples():
        s = source_index.loc[r.measured_feature]
        for source_col, value in [('log(IPF/Ctrl)', r.supplied_log_effect),
                                  ('IPF_vs_Ctrl_Ttest_p', r.supplied_raw_p),
                                  ('IPF_vs_Ctrl_Ttest_padj', r.supplied_adjusted_p)]:
            assert np.isclose(s[source_col], value, atol=1e-14, rtol=1e-12)
    rows = []
    for base, status in bases:
        source_magnitude = experimental_magnitude(measured['log(IPF/Ctrl)'], base)
        matched_magnitude = experimental_magnitude(details.supplied_log_effect, base)
        for statistic, p in [('BH', .1), ('BH', .05), ('raw_p', .05)]:
            sc = 'IPF_vs_Ctrl_Ttest_p' if statistic == 'raw_p' else 'IPF_vs_Ctrl_Ttest_padj'
            sp = details.supplied_raw_p if statistic == 'raw_p' else details.supplied_adjusted_p
            pp = details.predicted_raw_p if statistic == 'raw_p' else details.predicted_BH
            for measured_fc in [1.2, 1.5, 2.]:
                full_source = measured[sc].lt(p) & (source_magnitude > measured_fc)
                matched_source = details.status.eq('matched') & sp.lt(p) & (matched_magnitude > measured_fc)
                for predicted_fc in [1., 1.1, 1.2, 1.5, 2.]:
                    full_prediction = pp.lt(p) & details.predicted_signed_fold.abs().gt(predicted_fc)
                    overlap = matched_source & full_prediction
                    direction = details.supplied_log_effect * details.predicted_signed_fold
                    agree, oppose = overlap & direction.gt(0), overlap & direction.lt(0)
                    up, down = agree & details.supplied_log_effect.gt(0), agree & details.supplied_log_effect.lt(0)
                    # Independently reconcile the overlap using explicit lipid identity sets.
                    sm = set(details.loc[matched_source, 'model_lipid'])
                    pm = set(details.loc[full_prediction, 'model_lipid'])
                    assert set(details.loc[overlap, 'model_lipid']) == sm.intersection(pm)
                    assert agree.sum() + oppose.sum() == overlap.sum()
                    rows.append(dict(experimental_log_base=base, source_scale_status=status,
                        statistic=statistic, p_cutoff=p, minimum_measured_bulk_fold=measured_fc,
                        minimum_imputed_fold=predicted_fc,
                        measured_regulated_full_544=int(full_source.sum()),
                        measured_regulated_matchable_134=int(matched_source.sum()),
                        imputed_regulated_full_202=int(full_prediction.sum()),
                        measured_up=int((full_source & measured['log(IPF/Ctrl)'].gt(0)).sum()),
                        measured_down=int((full_source & measured['log(IPF/Ctrl)'].lt(0)).sum()),
                        imputed_up=int((full_prediction & details.predicted_signed_fold.gt(0)).sum()),
                        imputed_down=int((full_prediction & details.predicted_signed_fold.lt(0)).sum()),
                        concordant_up=int(up.sum()), concordant_down=int(down.sum()),
                        concordant_total=int(agree.sum()), discordant_total=int(oppose.sum()),
                        significant_in_both_and_both_fold_eligible=int(overlap.sum()),
                        concordance_percent=100*agree.sum()/overlap.sum() if overlap.sum() else np.nan,
                        concordant_lipids='; '.join(details.loc[agree, 'measured_feature']),
                        discordant_lipids='; '.join(details.loc[oppose, 'measured_feature'])))
    result = pd.DataFrame(rows)
    assert len(result) == 45 * len(bases)
    old = pd.read_csv(PRIOR / 'latest_model_all_60_settings.csv')
    keys = ['experimental_log_base', 'statistic', 'p_cutoff',
            'minimum_measured_bulk_fold', 'minimum_imputed_fold']
    checks = ['concordant_up', 'concordant_down', 'concordant_total',
              'discordant_total', 'significant_in_both_and_both_fold_eligible']
    joined = old.merge(result, on=keys, suffixes=('_old', '_new'), validate='one_to_one')
    assert len(joined) == len(old) == 60
    for col in checks:
        assert joined[col + '_old'].equals(joined[col + '_new'])
    selected = result[result.concordance_percent.gt(80)].copy()
    lines = ['# Latest bulk IPF results above 80 percent concordance', '',
        'Measured fold >1.2 is added to >1.5 and >2 with explicit user authorization. '
        'Models, predictions, source means, raw tests, BH families and matching are unchanged. '
        'All 202 model outputs are retained. Original results remain intact.', '',
        'Measured regulated counts include all 544 experimental lipids meeting the row’s '
        'significance and fold cutoffs. Imputed regulated counts include all 202 model outputs '
        'meeting that row’s cutoffs. These are separate totals, not overlap counts. '
        'Agreement is concordant / significant in both with both fold cutoffs. '
        'Only strict agreement >80% is shown; exactly 80% is excluded. '
        'Cutoffs are strict inequalities. Up/down counts denote concordant effects.', '',
        'The measured log base remains unconfirmed. Provisional log2 is the primary '
        'scenario; log10 is separately labeled sensitivity analysis. Predicted statistics '
        'remain donor-balanced Welch tests (20 IPF, 14 controls), with BH across 202 outputs. '
        'Measured statistics remain the supplied raw tests and BH across 544 features.', '']
    headers = ['Significance in both', 'Measured fold >', 'Imputed fold >',
               'Measured regulated', 'Imputed regulated', 'Concordant up',
               'Concordant down', 'Concordant / overlap', 'Discordant', 'Agreement']
    for base, group in selected.groupby('experimental_log_base', sort=False):
        lines += [f'## {"Provisional log2" if base == 2 else "Log10 sensitivity"}', '',
                  '| ' + ' | '.join(headers) + ' |',
                  '| ' + ' | '.join(['---']*len(headers)) + ' |']
        for r in group.itertuples():
            cells = [('Raw p' if r.statistic == 'raw_p' else 'BH') + f' <{r.p_cutoff:g}',
                     f'{r.minimum_measured_bulk_fold:g}',
                     'None' if r.minimum_imputed_fold == 1 else f'{r.minimum_imputed_fold:g}',
                     str(r.measured_regulated_full_544), str(r.imputed_regulated_full_202),
                     str(r.concordant_up), str(r.concordant_down),
                     f'{r.concordant_total}/{r.significant_in_both_and_both_fold_eligible}',
                     str(r.discordant_total), f'{r.concordance_percent:.1f}%']
            lines.append('| ' + ' | '.join(cells) + ' |')
        lines += ['']
    lines += ['## Lipid identities for each reported setting', '']
    for r in selected.itertuples():
        lines += [f'**Base {r.experimental_log_base:g}; {r.statistic} <{r.p_cutoff:g}; '
                  f'measured >{r.minimum_measured_bulk_fold:g}; imputed >{r.minimum_imputed_fold:g}**', '',
                  'Concordant: ' + r.concordant_lipids + '.', '',
                  'Discordant: ' + (r.discordant_lipids or 'None') + '.', '']
    OUT.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT / 'all_90_settings_with_regulated_counts.csv', index=False)
    selected.to_csv(OUT / 'settings_above80_with_regulated_counts.csv', index=False)
    (OUT / 'LATEST_BULK_IPF_ABOVE80.md').write_text('\n'.join(lines))
    receipt = dict(user_authorization=AUTHORIZATION, user_authorization_source='Current user message',
                   source_workbook=source, model=model, detail_source=str(detail_file),
                   detail_source_sha256=hashlib.sha256(detail_file.read_bytes()).hexdigest(),
                   complete_202_target_identities_verified=True, measured_features=544,
                   prior_60_settings_reproduced=True, set_based_overlap_checks_passed=True,
                   measured_log_base_verified=False, scale_policy=question['provisional_policy_authorization'],
                   training_performed=False, statistical_tests_changed=False,
                   selected_rows_by_base=selected.groupby('experimental_log_base').size().to_dict())
    (OUT / 'VERIFICATION.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print('\n'.join(lines[:lines.index('## Lipid identities for each reported setting')]))


if __name__ == '__main__':
    main()
