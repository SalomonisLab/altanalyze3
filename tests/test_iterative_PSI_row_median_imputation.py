"""Requested median filling must preserve observed zeros, source data and coverage."""
import anndata as ad
import numpy as np
import pandas as pd
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh
from altanalyze3.components.cellHarmony.iterative_marker_finder import (
    Scheme, run_legacy_iterative_marker_finder, assess_pattern_coverage,
)


def test_row_median_uses_all_samples_and_keeps_real_zero_and_unmeasured_event():
    raw=pd.DataFrame([[0.,np.nan,.2,.8,.9,1.], [np.nan]*6],index=['event','unmeasured'])
    before=raw.copy(deep=True)
    filled,mask=mh.fill_PSI_row_median(raw)
    pd.testing.assert_frame_equal(raw,before)
    np.testing.assert_array_equal(filled.loc['event'],[0.,.8,.2,.8,.9,1.])
    assert mask.loc['event'].tolist()==[False,True,False,False,False,False]
    assert filled.loc['unmeasured'].isna().all() and not mask.loc['unmeasured'].any()


def test_native_median_fit_keeps_complete_source_and_raw_replication(tmp_path):
    x=np.array([[0.,.1,np.nan],[np.nan,.2,np.nan],[.2,.15,np.nan],
                [.8,.8,np.nan],[.9,.9,np.nan],[1.,.85,np.nan]],dtype=np.float32)
    a=ad.AnnData(x.copy(),obs=pd.DataFrame({'group':['C']*3+['D']*3},index=list('abcdef')),
        var=pd.DataFrame(index=['gap','complete','all_missing']))
    a.uns['lineage_order']=['C','D']
    assignment,parameters=run_legacy_iterative_marker_finder(
        [Scheme('clinical',a,'group',['C','D'])],a,'group',tmp_path/'fit',
        rho_cutoff=.1,top_n=2,heatmap_dpi=100,missing_psi='row_median')
    np.testing.assert_array_equal(a.X,x)
    assert parameters['source_features']==list(a.var_names)
    assert parameters['schemes'][0]['original_missing_PSI_cells']==7
    assert parameters['schemes'][0]['imputed_PSI_cells']==1
    rho=pd.read_csv(tmp_path/'fit/clinical/rho.tsv',sep='\t',index_col=0)
    assert len(rho)==3 and rho.loc['all_missing'].isna().all()
    expected=np.corrcoef([0.,.8,.2,.8,.9,1.],[0,0,0,1,1,1])[0,1]
    assert np.isclose(abs(rho.loc['gap'].iloc[0]),expected)
    coverage=assess_pattern_coverage(a,assignment,'group',min_observed=3,min_fraction=1.)
    assert coverage.loc['gap','observed__C']==2 and not coverage.loc['gap','coverage_eligible']
    cache=ad.read_h5ad(parameters['figure']['heatmap_cache'])
    i=cache.var.gene.tolist().index('gap')
    assert np.isnan(cache.layers['source_PSI'][1,i]) and cache.layers['source_PSI'][0,i]==0.
    assert cache.layers['imputed_PSI'][1,i] and not cache.layers['imputed_PSI'][0,i]
