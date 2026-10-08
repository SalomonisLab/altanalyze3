"""Checks for full inputs, authorized median fills, and paired significance counts."""
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.rna2lipid.compare_bulk_ipf_models import model_inputs, signed_fold, thresholds


def test_missing_gene_is_filled_from_training_and_full_source_library_is_used():
    raw = pd.DataFrame([[10., 20.], [90., 80.]], index=['ensA', 'unmapped'], columns=['s1', 's2'])
    published = np.log2(raw)
    ann = pd.DataFrame({'gene_id': ['ensA', 'ensMissing'], 'gene': ['A', 'ADORA3']})
    genes = ['A', 'ADORA3']
    matrices, audit = model_inputs(raw, published, ann, genes, pd.Series([4., 7.], index=genes))
    assert list(matrices['log2_CP10k'].columns) == genes
    np.testing.assert_allclose(matrices['log2_CP10k'].A, np.log2([1001., 2001.]))
    np.testing.assert_array_equal(matrices['log2_CP10k'].ADORA3, [7., 7.])
    np.testing.assert_array_equal(matrices['publisher_logRPKM'].ADORA3, [7., 7.])
    assert audit.set_index('gene').loc['ADORA3', 'filled']
    with pytest.raises(ValueError, match='New unresolved'):
        model_inputs(raw, published, ann, ['A', 'UNREVIEWED'], pd.Series({'A': 4., 'UNREVIEWED': 7.}))


def test_same_threshold_applies_to_both_p_values_and_counts_up_down():
    frame = pd.DataFrame({'status': ['matched'] * 4 + ['absent_or_ambiguous'],
        'measured_feature': ['up', 'down', 'wrong', 'one_only', None],
        'supplied_log_effect': [1., -1., -1., 1., np.nan],
        'predicted_log2FC': [1., -1., 1., 1., 1.],
        'supplied_raw_p': [.01] * 4 + [np.nan], 'supplied_adjusted_p': [.02] * 4 + [np.nan],
        'predicted_raw_p': [.01, .01, .01, .2, .01], 'predicted_BH': [.03, .03, .03, .3, .03]})
    result = pd.DataFrame(thresholds(frame, 'candidate', 'log2_CP10k', 'donor_balanced'))
    row = result[result.statistic.eq('raw_p') & result.minimum_predicted_fold.eq(1)].iloc[0]
    assert row.significant_in_both == 3
    assert row.concordant_up == row.concordant_down == 1
    assert row.discordant_total == 1
    assert row.concordance_percent == pytest.approx(200/3)
    assert result[result.minimum_predicted_fold.eq(2)].significant_in_both.eq(0).all()  # strict >2


def test_signed_folds_are_ratios_not_log_values():
    np.testing.assert_array_equal(signed_fold([1., -1., 0., 2., -2.]), [2., -2., 0., 4., -4.])


def test_CPM_uses_full_source_library_and_approved_median():
    raw = pd.DataFrame([[10., 20.], [90., 80.]], index=['ensA', 'unmapped'], columns=['s1', 's2'])
    ann = pd.DataFrame({'gene_id': ['ensA'], 'gene': ['A']})
    matrices, _ = model_inputs(raw, np.log2(raw), ann, ['A', 'ADORA3'], pd.Series({'A': 4., 'ADORA3': 7.}))
    np.testing.assert_allclose(matrices['log2_CPM'].A, np.log2([100001., 200001.]))
    np.testing.assert_array_equal(matrices['log2_CPM'].ADORA3, [7., 7.])


def test_reference_only_selection_does_not_use_predicted_significance():
    from altanalyze3.components.rna2lipid.compare_bulk_ipf_models import reference_only
    frame = pd.DataFrame({'model_lipid': ['up', 'down', 'not_significant'], 'status': ['matched'] * 3,
        'supplied_adjusted_p': [.01, .02, .2], 'supplied_log_effect': [1., -1., 1.],
        'predicted_log2FC': [1., -1., -1.], 'predicted_geometric_log2FC': [1., 1., -1.],
        'predicted_raw_p': [.9, .9, .001], 'predicted_BH': [.9, .9, .001]})
    summaries, details = reference_only(frame, 'candidate', 'log2_CPM', 'donor_balanced')
    assert all(r['measured_lipids'] == 2 for r in summaries)
    assert [r['concordant_total'] for r in summaries] == [1, 2, 1, 2]
    assert list(details[0].model_lipid) == ['up', 'down']


def test_unresolved_RNA_source_blocks_inference():
    from altanalyze3.components.rna2lipid.compare_bulk_ipf_models import verified_RNA_source
    with pytest.raises(ValueError, match='identity unresolved'):
        verified_RNA_source({'questions': [{'id': 'lung_RNA_training_matrix_preprocessing_20261007', 'status': 'pending'}]})


def test_CE_correction_preserves_previous_matching_and_rejects_ambiguity():
    from altanalyze3.components.rna2lipid.compare_bulk_ipf_models import corrected_matching
    base = pd.DataFrame({'gene': ['PC(16:0/18:1)', 'CE(18:2)', 'CE(20:3)'],
        'canonical': ['PC(16:0/18:1)', None, None], 'ion_mode': ['P'] * 3,
        'measured_feature': ['PC_POS', None, None], 'status': ['matched', 'absent_or_ambiguous', 'absent_or_ambiguous']}).set_index('gene')
    raw = pd.DataFrame({'Lipid': ['CE(18:2)', 'CE(20:3)']}, index=['CE(18:2)_POS', 'CE(20:3)_POS'])
    fixed, changes = corrected_matching(base, raw)
    pd.testing.assert_series_equal(fixed.loc['PC(16:0/18:1)'], base.loc['PC(16:0/18:1)'])
    assert len(changes) == 2 and fixed.status.eq('matched').all()
    duplicated = pd.concat([raw, raw.iloc[[0]].rename(index={'CE(18:2)_POS': 'CE(18:2)_iso_B_POS'})])
    with pytest.raises(ValueError, match='Ambiguous CE mapping'):
        corrected_matching(base, duplicated)
