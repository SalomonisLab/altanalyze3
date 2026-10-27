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


def test_quality_review_retains_assignments_and_rejects_weak_sparse_patterns(tmp_path):
    from altanalyze3.components.cellHarmony.iterative_marker_finder import review_legacy_pattern_selection
    x=np.array([[.1,.1,.1],[.2,.2,np.nan],[.15,.15,.2],
                [.8,.8,.8],[.9,.9,.9],[.85,.85,1.]])
    a=ad.AnnData(x,obs=pd.DataFrame({'group':['control']*3+['case']*3,
                                   'donor':['C1','C2','C3','P1','P2','P3']},index=list('abcdef')),
                 var=pd.DataFrame(index=['strong','weak','sparse']))
    a.uns['lineage_order']=['control','case']
    assignments=pd.DataFrame({'pattern':['up-case']*3,'scheme':['clinical']*3,
                              'rho':[.9,.6,.95]},index=a.var_names)
    reviewed=review_legacy_pattern_selection(a,assignments,'group',rho_cutoff=.7,
                                            min_fraction=.9,min_observed=3)
    assert reviewed.index.tolist()==assignments.index.tolist()
    assert reviewed.figure_eligible.tolist()==[True,False,False]
    out=render_consolidated_marker_heatmap(a,reviewed.loc[reviewed.figure_eligible],cluster_key='group',
        out=tmp_path/'review.pdf',top_n=1,dpi=100,show_pattern_labels=True,
        sample_label_key='donor',figure_title='SRSF2 discovery',caption='Controls and cases. Missing PSI stays missing.')
    from pypdf import PdfReader
    text=PdfReader(out['pdf']).pages[0].extract_text()
    assert all(s in text for s in ['C1','P3','SRSF2 discovery','Missing PSI stays missing','up case'])
    assert out['plotted']==1 and np.array_equal(a.X,x,equal_nan=True)


def test_missing_PSI_is_displayed_black_without_replacing_nan(tmp_path,monkeypatch):
    from altanalyze3.components.visualization import marker_heatmap_h5ad as mh
    seen={}
    def capture(fig,*args,**kwargs):
        im=next(im for ax in fig.axes for im in ax.images if im.get_array().shape==(2,2))
        seen['missing_color']=im.cmap.get_bad()
        seen['mask']=np.ma.getmaskarray(im.get_array()).copy()
    monkeypatch.setattr(mh,'_save_figure',capture)
    matrix=pd.DataFrame([[0.,np.nan],[.1,-.1]],index=['a','b'],columns=['C','P'])
    mh._plot_heatmap(matrix,str(tmp_path/'color.pdf'),[('control',1),('case',1)],
        ['control','case'],['control','case'],['control','case'],dpi=100,label_events=[])
    assert np.array_equal(seen['missing_color'],[0,0,0,1])
    assert seen['mask'][0,1] and np.isnan(matrix.iloc[0,1])


def test_immutable_fit_review_hash_gate_full_roster_and_fixed_cap(tmp_path):
    import json,hashlib,pytest
    from altanalyze3.components.cellHarmony.iterative_marker_finder import review_legacy_markerfinder_fit
    weak='ENSG00000182511:E11.3-E12.2|ENSG00000182511:E10.2-E12.2'
    a=ad.AnnData(np.array([[.1,.2],[.15,.4],[.2,.3],[.8,.2],[.85,.4],[.9,.3]]),
        obs=pd.DataFrame({'group':['C']*3+['P']*3,'donor':list('abcdef')},index=list('abcdef')),
        var=pd.DataFrame(index=['strong',weak]))
    a.uns['lineage_order']=['C','P'];path=tmp_path/'source.h5ad';a.write_h5ad(path)
    digest=hashlib.sha256(path.read_bytes()).hexdigest();fit=tmp_path/'fit';fit.mkdir()
    (fit/'parameters.json').write_text(json.dumps(dict(source_features=a.var_names.tolist(),
        rho_strictly_greater_than=.4,schemes=[dict(source_samples=a.obs_names.tolist(),scored_samples=a.obs_names.tolist())])))
    pd.DataFrame(dict(pattern=['up-P','up-P'],scheme=['clinical']*2,
        direction=['up']*2,level=['P']*2,rho=[.95,.5]),
        index=a.var_names).to_csv(fit/'Consolidated-MarkerFinder.assignments.tsv',sep='\t')
    labels=tmp_path/'labels.tsv';labels.write_text('UID\n'+weak+'\n')
    manifest=dict(display_h5ad=str(path),display_group_key='group',label_events=str(labels),
        schemes=[dict(name='clinical',h5ad=str(path),group_key='group',level_order=['C','P'])])
    with pytest.raises(ValueError,match='SHA-256'):
        review_legacy_markerfinder_fit(manifest,fit,tmp_path/'bad',input_sha256='0'*64,dpi=100)
    assert not (tmp_path/'bad').exists()
    out=review_legacy_markerfinder_fit(manifest,fit,tmp_path/'reviewed',input_sha256=digest,
        top_n=1,sample_label_key='donor',dpi=100)
    audit=pd.read_csv(tmp_path/'reviewed/all_original_assignments_and_review_gates.tsv',sep='\t',index_col=0)
    assert audit.index.tolist()==['strong',weak] and audit.figure_eligible.tolist()==[True,False]
    assert out['plotted']==1 and out['assigned']==1 and out['unassigned_callout_references']==0
    assert ad.read_h5ad(out['heatmap_cache']).var.gene.tolist()==['strong']


def test_discrete_binary_complements_collapse_to_first_without_rho_identity():
    from altanalyze3.components.cellHarmony.iterative_marker_finder import Scheme, collapse_discrete_pattern_logic
    a=ad.AnnData(np.zeros((4,2)),obs=pd.DataFrame({'group':['healthy']*2+['SRSF2']*2},index=list('abcd')))
    source=pd.DataFrame({'scheme':['clinical']*2,'direction':['up','down'],
        'level':['healthy','SRSF2'],'pattern':['up-healthy','down-SRSF2'],
        'rho':[.41,.93],'signed_rho':[.41,-.93]},index=['one','two'])
    result,aliases=collapse_discrete_pattern_logic(source,[Scheme('clinical',a,'group',['healthy','SRSF2'])])
    assert result.pattern.tolist()==['up-healthy','up-healthy']
    pd.testing.assert_frame_equal(result.loc[source.index,source.columns.drop('pattern')],source.drop(columns='pattern'))
    assert result.loc['two','source_pattern']=='down-SRSF2'
    assert aliases.loc[aliases.source_pattern=='down-SRSF2','canonical_pattern'].iloc[0]=='up-healthy'
    reversed_result,_=collapse_discrete_pattern_logic(source,[Scheme('clinical',a,'group',['healthy','SRSF2'])],['SRSF2','healthy'])
    assert reversed_result.pattern.tolist()==['up-healthy','up-healthy']  # up categories still precede down


def test_discrete_pattern_logic_keeps_multigroup_and_partial_designs_distinct():
    from altanalyze3.components.cellHarmony.iterative_marker_finder import Scheme, collapse_discrete_pattern_logic
    a=ad.AnnData(np.zeros((6,2)),obs=pd.DataFrame({'group':['healthy']*2+['SRSF2']*2+['other']*2},index=list('abcdef')))
    source=pd.DataFrame({'scheme':['three','three','partial'],'direction':['up','down','up'],
        'level':['healthy','SRSF2','healthy'],'pattern':['up-healthy','down-SRSF2','up-healthy'],
        'rho':[.8,.8,.8]},index=['one','two','three'])
    result,_=collapse_discrete_pattern_logic(source,[Scheme('three',a,'group',['healthy','SRSF2','other']),
        Scheme('partial',a[:4].copy(),'group',['healthy','SRSF2'])])
    assert len(set(result.pattern))==3 and set(result.index)==set(source.index)


def test_consolidated_cap_retains_splicing_callout_beyond_top_n(tmp_path):
    x=np.array([[.1,.2],[.2,.3],[.8,.9],[.9,.8]])
    event='ENSG00000182511:E11.3-E12.2|ENSG00000182511:E10.2-E12.2'
    a=ad.AnnData(x,obs=pd.DataFrame({'group':['A','A','B','B']},index=list('abcd')),
                 var=pd.DataFrame(index=['first',event]))
    assignments=pd.DataFrame({'pattern':['up-B','up-B']},index=a.var_names)
    result=render_consolidated_marker_heatmap(a,assignments,cluster_key='group',
        out=tmp_path/'keep.pdf',top_n=1,label_events=[event],dpi=100)
    assert result['plotted']==2 and result['assigned']==2
    assert ad.read_h5ad(result['heatmap_cache']).var.gene.tolist()==a.var_names.tolist()


def test_unassigned_callouts_are_separate_references_without_marker_statistics(tmp_path):
    event='ENSG00000182511:E11.3-E12.2|ENSG00000182511:E10.2-E12.2'
    a=ad.AnnData(np.array([[.1,np.nan],[.2,.3],[.9,.4],[.8,.5]]),
        obs=pd.DataFrame({'group':['A','A','B','B']},index=list('abcd')),
        var=pd.DataFrame(index=['marker',event]))
    source=pd.DataFrame({'pattern':['up-B'],'rho':[.7]},index=['marker'])
    result=render_consolidated_marker_heatmap(a,source,cluster_key='group',
        out=tmp_path/'refs.pdf',top_n=1,label_events=[event],dpi=100,retain_unassigned_callouts=True)
    assert result['assigned']==1 and result['plotted']==2 and result['unassigned_callout_references']==1
    displayed=pd.read_csv(tmp_path/'refs.displayed_rows.tsv',sep='\t',index_col=0)
    assert displayed.loc[event,'display_role']=='unassigned_callout_reference'
    assert pd.isna(displayed.loc[event,'rho']) and 'not discovery-qualified' in displayed.loc[event,'pattern']
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path/'refs.assignments.tsv',sep='\t',index_col=0),source)
    cache=ad.read_h5ad(result['heatmap_cache'])
    assert cache.var.gene.tolist()==a.var_names.tolist() and np.isnan(cache.X[0,1])


def test_long_callout_fits_actual_vector_pdf(tmp_path):
    import pytest
    pdf=pytest.importorskip('pymupdf')
    import matplotlib.pyplot as plt
    from altanalyze3.components.visualization.marker_heatmap_h5ad import _draw_row_callouts
    fig,ax=plt.subplots(figsize=(8.3,10))
    fig.subplots_adjust(left=.26,right=.80)
    label='GAPDH:E6.10-E6.10_6537692|E6.10-E6.10_6537670'
    _draw_row_callouts(ax,[label],[1],3)
    path=tmp_path/'callout.pdf';fig.savefig(path);plt.close(fig)
    doc=pdf.open(path);page=doc[0]
    spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans']]
    callout=next(s for s in spans if s['text'].startswith('GAPDH:'))
    assert callout['text']==label and callout['bbox'][2] <= page.rect.width-2


def test_lineage_condition_mutation_precedence_keeps_nonredundant_down_roles():
    from altanalyze3.components.cellHarmony.iterative_marker_finder import Scheme,collapse_discrete_pattern_logic
    obs=pd.DataFrame({'clinical':['healthy','MDS','healthy','MDS'],
        'state':['ERP-5','ERP-5','LMPP-2','LMPP-2'],
        'joint':['ERP-5.healthy','ERP-5.MDS','LMPP-2.healthy','LMPP-2.MDS']},index=list('abcd'))
    a=ad.AnnData(np.zeros((4,16)),obs=obs)
    schemes=[Scheme('clinical',a,'clinical',['healthy','MDS']),
        Scheme('state_only',a,'state',['ERP-5','LMPP-2']),
        Scheme('state_condition',a,'joint',obs.joint.tolist())]
    roles=[(s.name,d,l) for s in schemes for d in ['up','down'] for l in s.levels()]
    original=pd.DataFrame([dict(scheme=s,direction=d,level=l,pattern=d+'-'+l,rho=.5+i/100.)
        for i,(s,d,l) in enumerate(roles)],index=['event'+str(i) for i in range(16)])
    priority=[(s,d,l) for d in ['up','down'] for s,l in [
        ('state_only','ERP-5'),*[('state_condition',x) for x in obs.joint],('clinical','MDS')]]
    result,aliases=collapse_discrete_pattern_logic(original,schemes,pattern_precedence=priority)
    assert list(dict.fromkeys(result.pattern))==[d+'-'+l for s,d,l in priority]
    assert len(set(result.pattern))==12 and set(result.index)==set(original.index)
    pd.testing.assert_frame_equal(result.loc[original.index,original.columns.drop('pattern')],original.drop(columns='pattern'))
    alias=aliases.set_index(['source_scheme','source_direction','source_level'])
    assert alias.loc[('state_only','up','LMPP-2'),'canonical_pattern']=='down-ERP-5'
    assert alias.loc[('clinical','up','healthy'),'canonical_pattern']=='down-MDS'
    assert alias.groupby('canonical_pattern').high_samples.nunique().eq(1).all()


def test_unassigned_callout_is_not_added_by_default(tmp_path):
    event='ENSG00000182511:E11.3-E12.2|ENSG00000182511:E10.2-E12.2'
    a=ad.AnnData(np.array([[.1,.2],[.2,.3],[.9,.4],[.8,.5]]),
        obs=pd.DataFrame({'group':['A','A','B','B']},index=list('abcd')),
        var=pd.DataFrame(index=['marker',event]))
    source=pd.DataFrame({'pattern':['up-B'],'rho':[.7]},index=['marker'])
    result=render_consolidated_marker_heatmap(a,source,cluster_key='group',
        out=tmp_path/'no_refs.pdf',top_n=1,label_events=[event],dpi=100)
    assert result['assigned']==result['plotted']==1 and result['unassigned_callout_references']==0
    assert ad.read_h5ad(result['heatmap_cache']).var.gene.tolist()==['marker']
    assert a.var_names.tolist()==['marker',event] and a.n_obs==4


def test_exact_highlight_callouts_do_not_promote_shared_exon_events(tmp_path):
    event='ENSG00000182511:E11.3-E12.2|ENSG00000182511:E10.2-E12.2'
    related='ENSG00000182511:E11.3-E13.1|ENSG00000182511:E10.2-E13.1'
    a=ad.AnnData(np.array([[.1,.2,.3],[.2,.3,.4],[.8,.7,.6],[.9,.8,.7]]),
        obs=pd.DataFrame({'group':['A','A','B','B']},index=list('abcd')),
        var=pd.DataFrame(index=['first',related,event]))
    source=pd.DataFrame({'pattern':['up-B']*3},index=a.var_names)
    result=render_consolidated_marker_heatmap(a,source,cluster_key='group',
        out=tmp_path/'exact.pdf',top_n=1,label_events=[event],dpi=100,label_shared_exons=False)
    assert ad.read_h5ad(result['heatmap_cache']).var.gene.tolist()==['first',event]
    labels=pd.read_csv(tmp_path/'exact_label_events.tsv',sep='\t')
    assert labels.event.tolist()==[event] and labels['match'].eq('listed event').all()
    assert a.n_vars==3 and result['assigned']==3 and result['plotted']==2


def test_explicit_seven_up_roles_precede_only_nonredundant_down_roles():
    from altanalyze3.components.cellHarmony.iterative_marker_finder import Scheme,collapse_discrete_pattern_logic
    obs=pd.DataFrame({'clinical':['healthy','MDS','healthy','MDS'],
        'state':['ERP-5','ERP-5','LMPP-2','LMPP-2'],
        'joint':['ERP-5.healthy','ERP-5.MDS','LMPP-2.healthy','LMPP-2.MDS']},index=list('abcd'))
    a=ad.AnnData(np.zeros((4,16)),obs=obs)
    schemes=[Scheme('clinical',a,'clinical',['healthy','MDS']),Scheme('state_only',a,'state',['ERP-5','LMPP-2']),Scheme('state_condition',a,'joint',obs.joint.tolist())]
    roles=[(s.name,d,l) for s in schemes for d in ['up','down'] for l in s.levels()]
    source=pd.DataFrame([dict(scheme=s,direction=d,level=l,pattern=d+'-'+l,rho=.5+i/100) for i,(s,d,l) in enumerate(roles)],index=['event'+str(i) for i in range(16)])
    order=[('state_only','ERP-5'),('state_condition','ERP-5.healthy'),('state_condition','ERP-5.MDS'),('state_only','LMPP-2'),('state_condition','LMPP-2.healthy'),('state_condition','LMPP-2.MDS'),('clinical','MDS')]
    priority=[(s,d,l) for d in ['up','down'] for s,l in order]
    result,aliases=collapse_discrete_pattern_logic(source,schemes,pattern_precedence=priority)
    assert result.pattern.drop_duplicates().tolist()==['up-'+l for s,l in order]+['down-'+l for s,l in order if s!='state_only']
    pd.testing.assert_frame_equal(result.loc[source.index,source.columns.drop('pattern')],source.drop(columns='pattern'))
    alias=aliases.set_index(['source_scheme','source_direction','source_level'])
    assert alias.loc[('state_only','down','ERP-5'),'canonical_pattern']=='up-LMPP-2'
    assert alias.loc[('state_only','down','LMPP-2'),'canonical_pattern']=='up-ERP-5'


def test_consolidated_zscore_is_opt_in_and_retains_roster_and_missingness(tmp_path):
    from altanalyze3.components.visualization.marker_heatmap_h5ad import _zscore_rows
    a=ad.AnnData(np.array([[.1,.5],[.2,.5],[.7,.5],[np.nan,.5]]),
        obs=pd.DataFrame({'group':['A','A','B','B']},index=list('abcd')),
        var=pd.DataFrame(index=['event','constant']))
    assignments=pd.DataFrame({'pattern':['up-B','up-A']},index=a.var_names)
    original=a.X.copy()
    median=render_consolidated_marker_heatmap(a,assignments,cluster_key='group',
        out=tmp_path/'median.pdf',dpi=100)
    z=render_consolidated_marker_heatmap(a,assignments,cluster_key='group',
        out=tmp_path/'z.pdf',normalization='row_zscore',half_range=1.5,dpi=100)
    m=ad.read_h5ad(median['heatmap_cache']);scaled=ad.read_h5ad(z['heatmap_cache'])
    assert m.obs_names.equals(scaled.obs_names) and m.var_names.equals(scaled.var_names)
    values=pd.DataFrame(original.T)
    expected=_zscore_rows(values).where(values.notna()).to_numpy().T
    np.testing.assert_allclose(scaled.X,expected,atol=1e-7,equal_nan=True)
    np.testing.assert_allclose(m.X,(values.sub(values.median(axis=1),axis=0)).to_numpy().T,atol=1e-7,equal_nan=True)
    np.testing.assert_array_equal(a.X,original)
    assert np.isnan(scaled.X[-1,0]) and np.all(scaled.X[:,1]==0)
    assert 'population SD' in z['normalization']
