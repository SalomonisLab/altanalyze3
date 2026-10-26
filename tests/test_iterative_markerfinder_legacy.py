import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import pearsonr
from altanalyze3.components.visualization.marker_heatmap_h5ad import (
    generate_marker_heatmap_from_adata, render_consolidated_marker_heatmap,
    _draw_covariate_legends)
from altanalyze3.components.cellHarmony.iterative_marker_finder import legacy_strongest_patterns


def test_native_hook_missing_statistics_and_legacy_tie(tmp_path):
    x = np.array([[.1,.9,np.nan],[.2,.8,np.nan],[.15,.85,np.nan],
                  [.8,.2,np.nan],[.9,.1,np.nan],[np.nan,.15,np.nan]])
    a = ad.AnnData(x, obs=pd.DataFrame({'group':['A']*3+['B']*3},index=list('abcdef')),
                  var=pd.DataFrame(index=['up','down','unmeasured']))
    result = generate_marker_heatmap_from_adata(a, cluster_key='group', out=str(tmp_path/'test.pdf'),
        marker_method='markerfinder',validate_scaling=False,centroid_method='mean',
        min_rho=.4,pval_threshold=1.,max_missing_fraction=1.,return_statistics=True,cells_per_cluster=0)
    assert result['rho'].index.tolist() == ['up','down','unmeasured']
    assert result['rho'].loc['unmeasured'].isna().all()
    observed=np.isfinite(x[:,0]);expected=pearsonr(x[observed,0],np.array([0,0,0,1,1,1])[observed])
    assert np.isclose(result['rho'].loc['up','B'],expected.statistic)
    assert np.isclose(result['pval'].loc['up','B'],expected.pvalue)
    # Exact complementary ties prefer the up candidate, as the Python2 code does.
    r=pd.DataFrame({'A':[-.8,.9],'B':[.8,-.9]},index=['one','two'])
    tied=dict(rho=r,pval=r*0+.01,fdr=r*0+.02)
    assigned=legacy_strongest_patterns([('scheme',tied)],['A','B'])
    assert assigned.loc['one','pattern']=='up-B'
    assert assigned.loc['two','pattern']=='up-A'
    assert (tmp_path/'test.pdf').exists()


def test_down_program_and_literal_legacy_order():
    r=pd.DataFrame({'A':[-.9,.5],'B':[.2,.3],'C':[.3,-.2]},index=['negative','positive'])
    result=dict(rho=r,pval=r*0+.01,fdr=r*0+.02)
    table=legacy_strongest_patterns([('scheme',result)],['A','B','C'])
    assert table.index.tolist()==['positive','negative']
    assert table.loc['negative','pattern']=='down-A'


def test_consolidated_PSI_keeps_missing_and_measured_zeros(tmp_path):
    x=np.array([[0.,.1],[.2,np.nan],[.8,.7],[np.nan,.9]])
    a=ad.AnnData(x,obs=pd.DataFrame({'group':['control','control','case','case']},index=list('abcd')),
                 var=pd.DataFrame(index=['event1','event2']))
    a.uns['lineage_order']=['control','case']
    assignments=pd.DataFrame({'pattern':['up-case','down-control']},index=a.var_names)
    out=render_consolidated_marker_heatmap(a,assignments,cluster_key='group',
        out=tmp_path/'consolidated.pdf',dpi=100)
    cached=ad.read_h5ad(out['heatmap_cache'])
    assert np.array_equal(np.isnan(cached.X),np.isnan(x))
    assert np.allclose(cached.X,x-np.nanmedian(x,axis=0),equal_nan=True)
    assert np.array_equal(a.X,x,equal_nan=True)
    assert cached.obs_names.tolist()==list('abcd')


def test_long_covariate_legends_wrap_without_overlapping_text():
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(8.3,11))
    cats=['AML-SRSF2','MDS-post-SRSF2','MDS-pre-SRSF2','aged','young']
    _draw_covariate_legends(fig,{'condition':{'categories':cats,'colors':dict.fromkeys(cats,'red')}})
    fig.canvas.draw();texts=[t for t in fig.axes[0].texts if t.get_text() in cats]
    assert len(texts)==len(cats)
    boxes=[t.get_window_extent(fig.canvas.get_renderer()) for t in texts]
    for i,a in enumerate(boxes):
        for b in boxes[i+1:]:assert not a.overlaps(b)
    plt.close(fig)


def test_dense_callouts_preserve_axis_and_fit_vertical_extent():
    import matplotlib.pyplot as plt
    from altanalyze3.components.visualization.marker_heatmap_h5ad import _draw_row_callouts
    fig,ax=plt.subplots(figsize=(8,8));ax.imshow(np.zeros((200,5)),aspect='auto')
    bounds=ax.get_ylim()
    _draw_row_callouts(ax,['E'+str(i) for i in range(100)],list(range(100)),200)
    assert ax.get_ylim()==bounds
    fig.canvas.draw();renderer=fig.canvas.get_renderer();boxes=[t.get_window_extent(renderer) for t in ax.texts]
    assert all(b.y0>=ax.bbox.y0-.1 and b.y1<=ax.bbox.y1+.1 for b in boxes)
    plt.close(fig)


def test_patient_holdout_keeps_all_timepoints_together_and_all_features(tmp_path):
    from altanalyze3.components.cellHarmony.iterative_marker_finder import (
        Scheme,assess_leave_patient_out_patterns)
    x=np.array([[.1,1,np.nan],[.2,1,np.nan],[.3,1,np.nan],
                [.9,1,np.nan],[.8,1,np.nan],[.7,1,np.nan]])
    obs=pd.DataFrame({'group':['control']*3+['case']*3,'patient':['P1','P2','P3']*2},index=list('abcdef'))
    a=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=['program','constant','unmeasured']))
    baseline=pd.DataFrame({'pattern':['up-case']},index=['program'])
    audit=assess_leave_patient_out_patterns([Scheme('scheme',a,'group',['control','case'])],
        baseline,'patient',tmp_path,max_missing_fraction=1.)
    import json
    roster=json.loads((tmp_path/'patient_holdout_rosters.json').read_text())
    assert all(len(r['held_samples'])==2 and len(r['retained_samples'])==4 for r in roster)
    full=pd.read_csv(tmp_path/'all_source_features_holdout_patterns.tsv',sep='\t',index_col=0)
    assert full.index.tolist()==a.var_names.tolist() and full.shape==(3,3)
    assert np.array_equal(a.X,x,equal_nan=True)
    assert audit.loc['program','holdouts']==3


def test_native_coverage_means_and_display_alias_preserve_sources(tmp_path):
    from altanalyze3.components.cellHarmony.iterative_marker_finder import assess_pattern_coverage
    x=np.array([[0.],[.2],[np.nan],[.8],[.9],[1.]])
    a=ad.AnnData(x,obs=pd.DataFrame({'group':['control']*3+['other-token']*3},index=list('abcdef')),
                 var=pd.DataFrame(index=['event']))
    a.uns['lineage_order']=['control','other-token']
    assignment=pd.DataFrame({'pattern':['up-other-token']},index=a.var_names)
    coverage=assess_pattern_coverage(a,assignment,'group',min_observed=2,min_fraction=.5)
    assert np.isclose(coverage.loc['event','mean_PSI__control'],.1)
    assert np.isclose(coverage.loc['event','measured_group_mean_PSI_range'],.8)
    out=render_consolidated_marker_heatmap(a,assignment,cluster_key='group',out=tmp_path/'alias.pdf',
        dpi=100,label_replacements={'other-token':'other-or-unknown'})
    from pypdf import PdfReader
    text=PdfReader(out['pdf']).pages[0].extract_text()
    assert 'other-or-unknown' in text and 'other-token' not in text
    assert np.array_equal(a.X,x,equal_nan=True)
    assert a.obs.group.tolist()==['control']*3+['other-token']*3


def test_legacy_gates_reject_unresolved_groups_and_infinity_before_scoring(tmp_path):
    import pytest
    from altanalyze3.components.cellHarmony.iterative_marker_finder import Scheme,run_legacy_iterative_marker_finder
    a=ad.AnnData(np.array([[.1],[.2],[.8],[.9]]),
        obs=pd.DataFrame({'group':['A','A','B','B']},index=list('abcd')),
        var=pd.DataFrame(index=['event']))
    with pytest.raises(ValueError,match='complete observed roster'):
        run_legacy_iterative_marker_finder([Scheme('bad',a,'group',['A','C'])],a,'group',tmp_path/'bad')
    a.X[0,0]=np.inf
    with pytest.raises(ValueError,match='Infinite measurements'):
        run_legacy_iterative_marker_finder([Scheme('bad',a,'group',['A','B'])],a,'group',tmp_path/'infinite')
    assert a.obs_names.tolist()==list('abcd') and a.var_names.tolist()==['event']
