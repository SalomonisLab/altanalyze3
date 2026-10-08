import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.rna2lipid.scan_independent_bulk_lipid_thresholds import (
    DUAL_POLICY, LOG2_POLICY, approved_bases, common_thresholds, experimental_magnitude, scan,
)


def test_unverified_source_log_base_cannot_be_used():
    with pytest.raises(ValueError, match='unverified'):
        experimental_magnitude([1., -1.], None)


def test_explicit_log_bases_have_different_abundance_inverses():
    np.testing.assert_array_equal(experimental_magnitude([1., -1.], 2.), [2., 2.])
    np.testing.assert_allclose(experimental_magnitude([1., -1.], np.e), [np.e, np.e])


def test_source_and_prediction_thresholds_are_applied_independently():
    frame = pd.DataFrame({'version':['A']*4, 'status':['matched']*4, 'measured_feature':['up','down','weak_measured','weak_predicted'],
        'supplied_log_effect':[np.log2(3),-np.log2(3),np.log2(1.2),np.log2(3)],
        'predicted_signed_fold':[1.2,-1.5,3.,1.01], 'supplied_raw_p':[.01]*4,
        'supplied_adjusted_p':[.02]*4,'predicted_raw_p':[.01]*4,'predicted_BH':[.02]*4})
    results=scan(frame,2.)
    row=results[(results.statistic=='BH')&(results.p_cutoff==.05)&(results.minimum_measured_bulk_fold==1.5)&(results.minimum_imputed_fold==1.1)].iloc[0]
    assert row.measured_significant_and_fold_eligible==3
    assert row.significant_in_both_and_both_fold_eligible==2
    assert row.concordant_up==row.concordant_down==1
    assert row.concordance_percent==100.
    strict=results[(results.statistic=='BH')&(results.p_cutoff==.05)&(results.minimum_measured_bulk_fold==2)&(results.minimum_imputed_fold==1.5)].iloc[0]
    assert strict.significant_in_both_and_both_fold_eligible==0
    assert strict.minimum_measured_bulk_fold==2 and results.minimum_measured_bulk_fold.max()==2


def test_generic_continue_and_diagnostic_permission_do_not_approve_scale_assumptions():
    for answer in ['continue', 'Best guess', DUAL_POLICY]:
        with pytest.raises(ValueError, match='unresolved'):
            approved_bases({'diagnostic_authorized': True, 'verified_log_base': None,
                            'provisional_policy_authorization': {'user_answer': answer}})


@pytest.mark.parametrize('answer,bases', [(LOG2_POLICY, [2.]), (DUAL_POLICY, [2., 10.])])
def test_explicit_assumption_policy_stays_unverified(answer, bases):
    # Synthetic authorization fixture; never writes approval to real decision state.
    question = {'verified_log_base': None,
                'provisional_log_base_assumption_authorized': True,
                'provisional_policy_authorization': {'granted_by': 'user',
                    'user_answer': answer, 'approved_assumed_log_bases': bases,
                    'response_evidence': 'Synthetic test response'}}
    assert approved_bases(question) == [(x, 'provisional_assumption') for x in bases]
    question['provisional_policy_authorization']['approved_assumed_log_bases'] = [np.e]
    with pytest.raises(ValueError, match='unresolved'):
        approved_bases(question)


def test_common_lipid_denominator_does_not_select_only_agreeing_predictions():
    frame = pd.DataFrame({'version': ['A', 'A', 'B', 'B', 'C', 'C'],
        'model_lipid': ['lipid1', 'lipid2'] * 3, 'status': ['matched'] * 6,
        'measured_feature': ['up', 'down'] * 3, 'supplied_log_effect': [1., -1.] * 3,
        'supplied_raw_p': [.01] * 6, 'supplied_adjusted_p': [.02] * 6,
        'predicted_raw_p': [.01] * 6, 'predicted_BH': [.02] * 6,
        'predicted_signed_fold': [1.5, -1.5, 1.5, 1.5, -1.5, -1.5]})
    result = common_thresholds(frame, 2.)
    selected = result[(result.statistic == 'BH') & (result.p_cutoff == .05)
                      & (result.minimum_measured_bulk_fold == 1.5)
                      & (result.minimum_imputed_fold == 1.1)]
    assert selected.identical_lipid_denominator.tolist() == [2, 2, 2]
    assert selected.concordant_total.tolist() == [2, 1, 1]
    assert selected.discordant_total.tolist() == [0, 1, 1]
