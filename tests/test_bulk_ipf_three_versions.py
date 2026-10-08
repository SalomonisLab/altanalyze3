"""Validate denominator separation and shared-lipid threshold comparisons."""
import numpy as np
import pandas as pd

from altanalyze3.components.rna2lipid.compare_three_bulk_ipf_versions import (
    VERSIONS, select_versions, frozen_reference_comparison, common_three_comparison)


def example():
    rows=[]
    for version, model, representation in VERSIONS:
        rows.append(pd.DataFrame({'version':version,'model':model,'RNA_representation':representation,'region':'donor_balanced',
            'model_lipid':['up','down','not_source_sig','pred_nonsig'], 'measured_feature':['UP','DOWN','OTHER','NO_PRED'],
            'status':'matched','supplied_log_effect':[1.,-1.,1.,1.],
            'supplied_raw_p':[.01,.01,.2,.01],'supplied_adjusted_p':[.02,.02,.3,.02],
            'predicted_raw_p':[.01,.01,.01,.8],'predicted_BH':[.02,.02,.02,.9],
            'predicted_signed_fold':[2.,-1.2,3.,2.],
            'predicted_log2FC':[1.,-np.log2(1.2),np.log2(3),1.],
            'predicted_geometric_log2FC':[1.,np.log2(1.2),np.log2(3),1.]}))
    return pd.concat(rows,ignore_index=True)


def test_three_version_selection_keeps_ordered_records_and_model_inputs():
    frame=example().drop(columns='version')
    # Two current versions share a model identity but have different RNA representations.
    chosen=select_versions(frame)
    assert len(chosen)==12 and list(chosen.version.drop_duplicates())==[v[0] for v in VERSIONS]
    assert chosen.groupby('version').size().eq(4).all()


def test_reference_only_counts_include_prediction_nonsignificant_lipid():
    result=frozen_reference_comparison(example())
    assert result.measured_significant.eq(3).all()
    assert result[result.effect_estimand=='arithmetic_abundance'].concordant_total.eq(3).all()
    assert result[result.effect_estimand=='mean_log_prediction'].concordant_total.eq(2).all()


def test_common_denominator_requires_every_version_to_pass_p_and_fold():
    frame=example()
    # One version's up lipid is nonsignificant; thus excluded in ALL versions' common set.
    frame.loc[(frame.version==VERSIONS[1][0]) & (frame.model_lipid=='up'),'predicted_BH']=.5
    result,detail=common_three_comparison(frame)
    rows=result[(result.statistic=='BH') & (result.p_cutoff==.05) & (result.minimum_predicted_fold==1.)]
    assert rows.identical_lipid_denominator.eq(1).all()
    assert rows.concordant_total.eq(1).all()
    strict=result[(result.statistic=='BH') & (result.p_cutoff==.05) & (result.minimum_predicted_fold==1.2)]
    assert strict.identical_lipid_denominator.eq(0).all()
    assert strict.concordance_percent.isna().all()
    assert set(detail[detail.statistic=='BH'].model_lipid)=={'down'}
