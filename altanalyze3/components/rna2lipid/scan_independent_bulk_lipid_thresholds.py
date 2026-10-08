"""Independent measured-bulk and RNA-imputed fold thresholds.

No default experimental log base: supplied effect scale must be verified first.
Uses all 202 outputs per version; statistical testing/BH are unchanged.
"""
from pathlib import Path
import json
import hashlib

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PHASE = 'bulk_IPF_separate_measured_and_predicted_folds_20261008'
SOURCE = HERE / 'artifacts/bulk_IPF_three_version_comparison_20261008'
OUT = SOURCE / 'independent_measured_and_imputed_fold_thresholds'
WORKBOOK_TABLES = HERE / 'artifacts/Geremy_original_lipid_statistics_20261008'
DUAL_POLICY = 'Run provisional log2 plus a separately labeled log10 sensitivity analysis'
LOG2_POLICY = 'Run provisional log2 only'


def approved_bases(question):
    """An assumption needs explicit user approval; never label it verified."""
    if (question.get('status') == 'resolved' and question.get('verified_log_base')
            and question.get('resolution_evidence')):
        return [(float(question['verified_log_base']), 'verified')]
    policy = question.get('provisional_policy_authorization', {})
    answer = policy.get('user_answer')
    expected = {DUAL_POLICY: [2., 10.], LOG2_POLICY: [2.]}.get(answer)
    if (not expected or policy.get('granted_by') != 'user'
            or not policy.get('response_evidence')
            or question.get('provisional_log_base_assumption_authorized') is not True
            or policy.get('approved_assumed_log_bases') != expected):
        raise ValueError('Experimental effect scale unresolved; no scientific threshold outputs written')
    if question.get('verified_log_base') is not None:
        raise ValueError('Do not confuse an assumed base with verified source provenance')
    return [(base, 'provisional_assumption') for base in expected]


def experimental_magnitude(effects, log_base):
    if log_base is None or not np.isfinite(log_base) or log_base <= 1:
        raise ValueError('Experimental log base is unverified; ask before applying source fold cutoffs')
    return np.power(float(log_base), np.abs(np.asarray(effects, float)))


def scan(frame, log_base):
    """Keep denominators explicit; source and prediction fold gates are independent."""
    source_fold = experimental_magnitude(frame.supplied_log_effect, log_base)
    records = []
    for name, group in frame.groupby('version', sort=False):
        magnitude = source_fold[frame.version.eq(name).to_numpy()]
        direction = group.supplied_log_effect * group.predicted_signed_fold
        for statistic, cutoff in [('raw_p', .05), ('BH', .05), ('BH', .1)]:
            sp = group.supplied_raw_p if statistic == 'raw_p' else group.supplied_adjusted_p
            pp = group.predicted_raw_p if statistic == 'raw_p' else group.predicted_BH
            for measured_fold in [1.5, 2.]:
                measured = group.status.eq('matched') & sp.lt(cutoff) & (magnitude > measured_fold)
                for predicted_fold in [1., 1.1, 1.2, 1.5, 2.]:
                    overlap = measured & pp.lt(cutoff) & group.predicted_signed_fold.abs().gt(predicted_fold)
                    concordant = overlap & direction.gt(0)
                    discordant = overlap & direction.lt(0)
                    records.append({'version': name, 'statistic': statistic, 'p_cutoff': cutoff,
                        'minimum_measured_bulk_fold': measured_fold, 'minimum_imputed_fold': predicted_fold,
                        'measured_significant_and_fold_eligible': int(measured.sum()),
                        'significant_in_both_and_both_fold_eligible': int(overlap.sum()),
                        'concordant_up': int((concordant & group.supplied_log_effect.gt(0)).sum()),
                        'concordant_down': int((concordant & group.supplied_log_effect.lt(0)).sum()),
                        'concordant_total': int(concordant.sum()), 'discordant_total': int(discordant.sum()),
                        'concordance_percent': 100 * concordant.sum() / overlap.sum() if overlap.sum() else np.nan,
                        'concordant_lipids': '; '.join(group.loc[concordant, 'measured_feature']),
                        'discordant_lipids': '; '.join(group.loc[discordant, 'measured_feature'])})
    return pd.DataFrame(records)


def comparison_report(details, results, common, assumption_summary):
    from .audit_geremy_lipid_workbooks import markdown
    columns = ['version', 'statistic', 'p_cutoff', 'minimum_measured_bulk_fold',
               'minimum_imputed_fold', 'measured_significant_and_fold_eligible',
               'significant_in_both_and_both_fold_eligible', 'concordant_up',
               'concordant_down', 'concordant_total', 'discordant_total',
               'concordance_percent']
    lines = ['# Bulk IPF: independent measured and imputed fold cutoffs', '',
             'This comparison uses the full-precision IPF workbook and the existing '
             'verified predictions. Fixed ElasticNetCV models and their 202-output '
             'panels are unchanged. No fitting, deployment or statistical-test change '
             'is performed. All prior result directories are preserved.', '',
             'Measured cutoff magnitudes are strictly >1.5 or >2. Imputed cutoff '
             'magnitudes are strictly >1, >1.1, >1.2, >1.5 or >2. The >1 condition '
             'means no additional effect-size cutoff. Significance is required in both '
             'datasets at raw p <0.05, BH <0.05 or BH <0.1. BH uses 544 measured '
             'features and 202 prediction features, independently.', '',
             'Measured means are centered values. Their reported effect is IPF mean '
             'minus control mean. Measured fold magnitude = base^abs(effect). Each '
             'unverified base is expressly labeled provisional; identical source p-values '
             'and signs are retained in every base scenario. Under a log2 assumption '
             'the inferred measured fold is a ratio of normalized geometric means, '
             'not an absolute MS concentration.', '',
             'Prediction folds are the established ratio of arithmetic donor-group '
             'means on each model’s linear output scale. Tests use balanced donor log '
             'predictions, with 20 IPF and 14 controls. These are independent disease '
             'and control donors; the RNA and lipidomics cohorts are unmatched.', '',
             'All three versions share 134 established measured-feature matches. '
             'Every one of the 202 predicted lipids remains in the detail table and '
             'prediction BH family. Neither matching nor fold cutoffs alter model '
             'coverage. The full measured panel and matchable eligible counts are '
             'both provided.', '',
             'A higher percentage after threshold selection describes the retained '
             'subset. Model-specific overlap denominators may differ. The common-lipid '
             'tables compare identical eligible lipids across all three versions.', '',
             '## Three versions at the same thresholds', '',
             'The first table uses BH <0.1 in both datasets, measured absolute fold '
             '>1.5 and imputed absolute fold >1.2. Under provisional log2 this '
             'recovers the largest number of concordant lipids among the scanned '
             'settings exceeding 75% concordance. This is a descriptive selection '
             'from the reported threshold grid, not a held-out accuracy estimate.', '',
             markdown(results[results.statistic.eq('BH') & results.p_cutoff.eq(.1)
                              & results.minimum_measured_bulk_fold.eq(1.5)
                              & results.minimum_imputed_fold.eq(1.2)]
                              [['experimental_log_base'] + columns]), '',
             'At this selected setting, only three lipids pass all three model '
             'versions under provisional log2. Every version agrees for two of '
             'those three lipids. The larger model-specific recovery counts therefore '
             'include different significant subsets; the percentage comparison '
             'alone does not establish improvement on an identical broad lipid set.', '',
             markdown(common[common.statistic.eq('BH') & common.p_cutoff.eq(.1)
                             & common.minimum_measured_bulk_fold.eq(1.5)
                             & common.minimum_imputed_fold.eq(1.2)]), '',
             'A stricter example uses BH <0.05 in both datasets, measured fold >2 '
             'and imputed fold >1.2. The overlap size must be considered alongside '
             'percentage concordance.', '',
             markdown(results[results.statistic.eq('BH') & results.p_cutoff.eq(.05)
                              & results.minimum_measured_bulk_fold.eq(2.)
                              & results.minimum_imputed_fold.eq(1.2)]
                              [['experimental_log_base'] + columns]), '',
             '## Source scale scenarios', '', markdown(assumption_summary), '']
    lines += ['All 544 IPF and 530 BPD original mean/statistics rows are documented '
              'separately in [the original-source tables]'
              '(../../Geremy_original_lipid_statistics_20261008/ORIGINAL_MEANS_AND_STATISTICS.md). '
              'BPD is not used for this IPF model comparison.', '']
    for base, group in results.groupby('experimental_log_base', sort=False):
        status = group.iloc[0].source_scale_status
        lines += [f'## Base {base:g}: {status}', '', markdown(group[columns]), '']
        comparable = common[common.experimental_log_base.eq(base)]
        lines += ['### Same lipids significant in all three versions', '',
                  markdown(comparable), '']
        for r in group.itertuples():
            lines += [f'**{r.version}; {r.statistic} <{r.p_cutoff:g}; '
                      f'measured fold >{r.minimum_measured_bulk_fold:g}; '
                      f'imputed fold >{r.minimum_imputed_fold:g}**', '',
                      'Concordant: ' + (r.concordant_lipids or 'None') + '.', '',
                      'Discordant: ' + (r.discordant_lipids or 'None') + '.', '']
    display = ['model_lipid', 'measured_feature', 'status', 'measured_control_mean',
               'measured_control_sd', 'measured_IPF_mean', 'measured_IPF_sd',
               'supplied_log_effect', 'supplied_raw_p', 'supplied_adjusted_p',
               'predicted_control_mean_linear', 'predicted_IPF_mean_linear',
               'predicted_signed_fold', 'predicted_raw_p', 'predicted_BH']
    display[8:8] = [c for c in details if c.startswith('measured_signed_fold_assumed_')]
    for name, group in details.groupby('version', sort=False):
        lines += [f'## Every original mean and result: {name}', '',
                  markdown(group[display]), '']
    return '\n'.join(lines)


def common_thresholds(frame, log_base):
    """Keep the same lipid denominator for comparisons among all three versions."""
    groups = {name: group.set_index('model_lipid')
              for name, group in frame.groupby('version', sort=False)}
    first = next(iter(groups.values()))
    rows = []
    for statistic, p in [('raw_p', .05), ('BH', .05), ('BH', .1)]:
        source_p = 'supplied_raw_p' if statistic == 'raw_p' else 'supplied_adjusted_p'
        prediction_p = 'predicted_raw_p' if statistic == 'raw_p' else 'predicted_BH'
        for measured_fc in [1.5, 2.]:
            source = (first.status.eq('matched') & first[source_p].lt(p)
                      & (experimental_magnitude(first.supplied_log_effect, log_base) > measured_fc))
            for predicted_fc in [1., 1.1, 1.2, 1.5, 2.]:
                common = source.copy()
                for group in groups.values():
                    assert group.index.tolist() == first.index.tolist()
                    common &= group[prediction_p].lt(p) & group.predicted_signed_fold.abs().gt(predicted_fc)
                for name, group in groups.items():
                    direction = group.supplied_log_effect * group.predicted_signed_fold
                    agree, oppose = common & direction.gt(0), common & direction.lt(0)
                    rows.append({'version': name, 'statistic': statistic, 'p_cutoff': p,
                        'minimum_measured_bulk_fold': measured_fc, 'minimum_imputed_fold': predicted_fc,
                        'identical_lipid_denominator': int(common.sum()),
                        'concordant_up': int((agree & group.supplied_log_effect.gt(0)).sum()),
                        'concordant_down': int((agree & group.supplied_log_effect.lt(0)).sum()),
                        'concordant_total': int(agree.sum()), 'discordant_total': int(oppose.sum()),
                        'concordance_percent': float(100 * agree.sum() / common.sum()) if common.sum() else np.nan})
    return pd.DataFrame(rows)


def main():
    state = json.loads((HERE / 'integrity/decision_state.json').read_text())
    role = next(q for q in state['questions'] if q['id'] == 'bulk_IPF_separate_fold_threshold_role_20261008')
    base = next(q for q in state['questions'] if q['id'] == 'bulk_IPF_experimental_log_effect_base_20261008')
    if role['status'] != 'resolved' or role.get('user_answer') != 'Measured bulk lipidomics: 1.5 and 2; imputed lipids: separate cutoffs':
        raise ValueError('Independent cutoff roles must be explicitly resolved first')
    bases = approved_bases(base)
    evidence = json.loads((SOURCE / 'THREE_VERSION_AUDIT.json').read_text())
    if not evidence['all_checks_passed']:
        raise ValueError('Baseline three-version comparison is unverified')
    workbook_audit = json.loads((WORKBOOK_TABLES / 'workbook_audit.json').read_text())
    source_info = workbook_audit['sources']['IPF']
    if hashlib.sha256(Path(source_info['path']).read_bytes()).hexdigest() != source_info['sha256']:
        raise ValueError('Experimental workbook changed; source reconciliation required')
    frame = pd.read_csv(WORKBOOK_TABLES / 'three_versions_all_202_with_original_means.csv',
                        float_precision='round_trip')
    contract = json.loads((HERE / 'integrity/baseline_contract.json').read_text())
    for _, group in frame.groupby('version'):
        if len(group) != 202 or group.model_lipid.tolist() != contract['lipids']:
            raise ValueError('Complete declared target panel required; no intersection')
    measured = pd.read_excel(source_info['path'])
    result, common, assumptions = [], [], []
    for log_base, status in bases:
        if status == 'provisional_assumption':
            label = f'measured_signed_fold_assumed_log{log_base:g}'
        else:
            label = f'measured_signed_fold_verified_log{log_base:g}'
        frame[label] = np.sign(frame.supplied_log_effect) * experimental_magnitude(frame.supplied_log_effect, log_base)
        item, shared = scan(frame, log_base), common_thresholds(frame, log_base)
        item.insert(0, 'source_scale_status', status)
        item.insert(0, 'experimental_log_base', log_base)
        shared.insert(0, 'experimental_log_base', log_base)
        result.append(item)
        common.append(shared)
        for statistic, cutoff in [('raw_p', .05), ('BH', .05), ('BH', .1)]:
            pcol = 'IPF_vs_Ctrl_Ttest_p' if statistic == 'raw_p' else 'IPF_vs_Ctrl_Ttest_padj'
            for fold in [1.5, 2.]:
                eligible = measured[pcol].lt(cutoff) & (experimental_magnitude(measured['log(IPF/Ctrl)'], log_base) > fold)
                assumptions.append({'experimental_log_base': log_base, 'scale_status': status,
                    'statistic': statistic, 'p_cutoff': cutoff, 'minimum_measured_bulk_fold': fold,
                    'full_544_measured_eligible': int(eligible.sum()),
                    'increased': int((eligible & measured['log(IPF/Ctrl)'].gt(0)).sum()),
                    'decreased': int((eligible & measured['log(IPF/Ctrl)'].lt(0)).sum())})
    result, common, assumptions = pd.concat(result, ignore_index=True), pd.concat(common, ignore_index=True), pd.DataFrame(assumptions)
    assert len(result) == len(common) == 90 * len(bases)
    OUT.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT / 'all_independent_threshold_results.csv', index=False)
    common.to_csv(OUT / 'common_lipids_all_three_versions.csv', index=False)
    assumptions.to_csv(OUT / 'full_measured_panel_threshold_counts.csv', index=False)
    frame.to_csv(OUT / 'all_202_original_means_and_statistics.csv', index=False)
    (OUT / 'BULK_IPF_INDEPENDENT_THRESHOLDS.md').write_text(comparison_report(frame, result, common, assumptions))
    receipt = {'experimental_source': source_info, 'source_scale_scenarios': bases,
               'user_scale_policy': base.get('provisional_policy_authorization'),
               'log_base_is_verified': base.get('verified_log_base') is not None,
               'three_version_baseline_verified': evidence['all_checks_passed'],
               'all_202_target_identities_retained': True,
               'threshold_rows': len(result), 'training_performed': False,
               'deployment_changed': False, 'prediction_tests_changed': False}
    (OUT / 'threshold_audit.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
