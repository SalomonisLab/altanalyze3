"""Independently verify threshold outputs from original workbook and saved stats."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


HERE = Path(__file__).resolve().parent
ROOT = HERE / 'artifacts/bulk_IPF_three_version_comparison_20261008'
OUT = ROOT / 'independent_measured_and_imputed_fold_thresholds'
VERSIONS = [('Abby supplied model / CPM', 'prior', 'log2_CPM'),
            ('Last corrected model / CP10k', 'candidate', 'log2_CP10k'),
            ('New corrected model / CPM', 'candidate', 'log2_CPM')]


def same_percentage(actual, numerator, denominator):
    expected = 100 * numerator / denominator if denominator else np.nan
    np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=0, equal_nan=True)


def main():
    contract = json.loads((HERE / 'integrity/baseline_contract.json').read_text())
    source = pd.read_excel('/Users/saljh8/Downloads/Geremy/10_results_with_statistics.xlsx')
    assert len(source) == 544
    assert source['Unnamed: 0'].is_unique
    measured = source.set_index('Unnamed: 0')
    np.testing.assert_allclose(false_discovery_control(source.IPF_vs_Ctrl_Ttest_p),
                               source.IPF_vs_Ctrl_Ttest_padj, atol=1e-12, rtol=0)
    mapping = pd.read_csv(ROOT / 'full_panel_rerun/all_202_lipid_matching.csv')
    # This is the harmonized 202-row map shared by all versions, labeled with
    # the candidate that established it, rather than one map per comparator.
    mapping_model = mapping.set_index('gene')
    assert mapping_model.index.tolist() == contract['lipids']
    declared = pd.read_csv(ROOT / 'three_versions_all_202_lipids.csv')
    expected = {}
    for version, model, representation in VERSIONS:
        version_declared = declared[declared.version.eq(version)].set_index('model_lipid')
        assert version_declared.index.tolist() == contract['lipids']
        for column in ['measured_feature', 'status']:
            assert version_declared[column].fillna('').tolist() == mapping_model[column].fillna('').tolist()
        stats = pd.read_csv(ROOT / 'full_panel_rerun' /
                            f'{model}_{representation}_donor_balanced_all_202_differentials.csv',
                            float_precision='round_trip').set_index('feature')
        assert stats.index.tolist() == contract['lipids']
        assert np.isfinite(stats[['IPF_mean', 'control_mean', 'pvalue', 'FDR']]).all().all()
        assert (stats[['IPF_mean', 'control_mean']] > 0).all().all()
        # Derive direction and fold magnitude directly from saved group abundances.
        ratio = stats.IPF_mean.to_numpy() / stats.control_mean.to_numpy()
        # The established display convention uses zero for an exactly unchanged
        # lipid; retain that convention rather than inventing a +1 direction.
        signed = np.where(ratio > 1, ratio, np.where(ratio < 1, -1 / ratio, 0.))
        np.testing.assert_allclose(signed, stats.signed_fold, atol=1e-12, rtol=1e-12)
        np.testing.assert_allclose(false_discovery_control(stats.pvalue), stats.FDR,
                                   atol=1e-12, rtol=0)
        mapped = mapping_model.status.eq('matched').to_numpy()
        ids = mapping_model.measured_feature.to_numpy()
        assert mapped.sum() == 134
        assert all(x in measured.index for x in ids[mapped])
        panel = measured.reindex(ids)
        expected[version] = {'mapped': mapped, 'ids': ids, 'signed': signed,
            'direction': np.sign(stats.IPF_mean.to_numpy() - stats.control_mean.to_numpy()),
            'effect': panel['log(IPF/Ctrl)'].to_numpy(),
            'source_raw_p': panel.IPF_vs_Ctrl_Ttest_p.to_numpy(),
            'source_BH': panel.IPF_vs_Ctrl_Ttest_padj.to_numpy(),
            'raw_p': stats.pvalue.to_numpy(), 'BH': stats.FDR.to_numpy()}

    result = pd.read_csv(OUT / 'all_independent_threshold_results.csv')
    common = pd.read_csv(OUT / 'common_lipids_all_three_versions.csv')
    full_source = pd.read_csv(OUT / 'full_measured_panel_threshold_counts.csv')
    detail = pd.read_csv(OUT / 'all_202_original_means_and_statistics.csv',
                         float_precision='round_trip')
    assert len(detail) == 606
    for version, group in detail.groupby('version', sort=False):
        assert group.model_lipid.tolist() == contract['lipids']
        item = expected[version]
        for base in [2., 10.]:
            ratio = np.exp(item['effect'] * np.log(base))
            signed = np.where(item['effect'] > 0, ratio,
                              np.where(item['effect'] < 0, -1 / ratio, 0.))
            signed[~item['mapped']] = np.nan
            np.testing.assert_allclose(group[f'measured_signed_fold_assumed_log{base:g}'],
                                       signed, atol=1e-11, rtol=1e-12, equal_nan=True)
    assert len(result) == len(common) == 180 and len(full_source) == 12
    assert set(result.experimental_log_base) == {2., 10.}
    assert set(result.source_scale_status) == {'provisional_assumption'}
    for row in result.itertuples():
        item = expected[row.version]
        # Equivalent criterion in the original effect coordinates avoids relying
        # on the implementation's exponentiation-based eligibility calculation.
        effect_threshold = np.log(row.minimum_measured_bulk_fold) / np.log(row.experimental_log_base)
        sp = item['source_raw_p'] if row.statistic == 'raw_p' else item['source_BH']
        pp = item['raw_p'] if row.statistic == 'raw_p' else item['BH']
        measured_pass = item['mapped'] & (sp < row.p_cutoff) & (abs(item['effect']) > effect_threshold)
        overlap = measured_pass & (pp < row.p_cutoff) & (abs(item['signed']) > row.minimum_imputed_fold)
        agree = overlap & (np.sign(item['effect']) == item['direction'])
        oppose = overlap & (np.sign(item['effect']) != item['direction'])
        assert int(measured_pass.sum()) == row.measured_significant_and_fold_eligible
        assert int(overlap.sum()) == row.significant_in_both_and_both_fold_eligible
        assert int(agree.sum()) == row.concordant_total
        assert int(oppose.sum()) == row.discordant_total
        assert int((agree & (item['effect'] > 0)).sum()) == row.concordant_up
        assert int((agree & (item['effect'] < 0)).sum()) == row.concordant_down
        same_percentage(row.concordance_percent, agree.sum(), overlap.sum())
        for text, flag in [(row.concordant_lipids, agree), (row.discordant_lipids, oppose)]:
            actual = [] if pd.isna(text) else text.split('; ')
            assert actual == item['ids'][flag].tolist()

    for row in common.itertuples():
        base_item = expected[VERSIONS[0][0]]
        threshold = np.log(row.minimum_measured_bulk_fold) / np.log(row.experimental_log_base)
        sp = base_item['source_raw_p'] if row.statistic == 'raw_p' else base_item['source_BH']
        mask = base_item['mapped'] & (sp < row.p_cutoff) & (abs(base_item['effect']) > threshold)
        for item in expected.values():
            pp = item['raw_p'] if row.statistic == 'raw_p' else item['BH']
            mask &= (pp < row.p_cutoff) & (abs(item['signed']) > row.minimum_imputed_fold)
        item = expected[row.version]
        agree = mask & (np.sign(item['effect']) == item['direction'])
        oppose = mask & (np.sign(item['effect']) != item['direction'])
        assert int(mask.sum()) == row.identical_lipid_denominator
        assert int(agree.sum()) == row.concordant_total
        assert int(oppose.sum()) == row.discordant_total
        assert int((agree & (item['effect'] > 0)).sum()) == row.concordant_up
        assert int((agree & (item['effect'] < 0)).sum()) == row.concordant_down
        same_percentage(row.concordance_percent, agree.sum(), mask.sum())
    for row in full_source.itertuples():
        threshold = np.log(row.minimum_measured_bulk_fold) / np.log(row.experimental_log_base)
        p = source.IPF_vs_Ctrl_Ttest_p if row.statistic == 'raw_p' else source.IPF_vs_Ctrl_Ttest_padj
        effect = source['log(IPF/Ctrl)']
        passing = (p < row.p_cutoff) & (abs(effect) > threshold)
        assert int(passing.sum()) == row.full_544_measured_eligible
        assert int((passing & (effect > 0)).sum()) == row.increased
        assert int((passing & (effect < 0)).sum()) == row.decreased
    receipt = {'all_checks_passed': True, 'verified_threshold_rows': len(result),
               'verified_same_lipid_comparison_rows': len(common),
               'verified_full_measured_panel_rows': len(full_source),
               'all_concordant_and_discordant_lipid_lists_verified': True,
               'original_source_BH_verified': True, 'full202_prediction_BH_verified': True,
               'prediction_fold_ratios_verified_from_group_means': True,
               'all_606_detail_rows_and_assumed_measured_folds_verified': True,
               'verified_log_base': None, 'assumed_log_bases': [2, 10],
               'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT / 'independent_verification.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
