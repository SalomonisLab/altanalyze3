import anndata as ad
import numpy as np
import pandas as pd
import pytest
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh


def test_saved_signed_roles_preserve_source_and_native_zscore(tmp_path, monkeypatch):
    groups=['healthy','pre','post','AML']
    x=np.array([[.9,.1,.3,.7],[.8,.2,.4,.6],[.3,.7,.9,.1],[.4,.6,.8,.2],
                [.3,.7,.3,.7],[.4,.6,.4,.6],[.3,.7,.3,.7],[np.nan,.6,.4,.6]])
    a=ad.AnnData(x,obs=pd.DataFrame({'group':np.repeat(groups,2)},index=['s'+str(i) for i in range(8)]),var=pd.DataFrame(index=['g'+str(i) for i in range(4)]))
    a.uns['lineage_order']=groups
    rho=pd.DataFrame([[.9,-.2,-.2,-.2],[-.9,.2,.2,.2],[-.2,.9,-.2,-.2],[.2,-.9,.2,.2]],index=a.var_names,columns=groups)
    fdr=pd.DataFrame([[1e-5,.2,.2,.2],[1e-5,.2,.2,.2],[.2,1e-5,.2,.2],[.2,1e-5,.2,.2]],index=a.var_names,columns=groups)
    captured={}
    monkeypatch.setattr(mh,'_plot_heatmap',lambda matrix,out,counts,order,cols,rows,**kw:captured.update(matrix=matrix,rows=rows,order=order))
    before=a.X.copy()
    r=mh.render_signed_marker_heatmap_from_saved_statistics(a,rho,fdr,fdr,cluster_key='group',out=tmp_path/'test.pdf',top_n=1)
    assert captured['rows']==['up healthy','up pre','down healthy','down pre']
    assert set(captured['rows'])<=set(captured['order'])
    np.testing.assert_array_equal(a.X,before)
    z=(before-np.nanmean(before,axis=0))/np.nanstd(before,axis=0)
    expected=np.nan_to_num(z[:,[0,2,1,3]],nan=0)
    cache=ad.read_h5ad(r['heatmap_cache'])
    np.testing.assert_allclose(cache.X,expected,rtol=1e-6)
    assert cache.obs_names.tolist()==a.obs_names.tolist()
    with pytest.raises(ValueError,match='every ordered source event'):
        mh.render_signed_marker_heatmap_from_saved_statistics(a,rho.iloc[:-1],fdr,fdr,cluster_key='group',out=tmp_path/'bad.pdf')


def test_cache_intronic_bulk_callouts_are_red_and_samples_retained(tmp_path):
    import pymupdf
    rows=['ENSG000001:I5.1-E6.1|ENSG000001:E4.1-E6.1',
          'ENSG000002:E3.5-E4.1|ENSG000002:E2.1-E4.1',
          'ENSG000003:E3.1-E4.1|ENSG000003:E2.1-E4.1']
    listed=['X:ENSG000001:I5.2-E9.1|ENSG000001:E8.1-E9.1',
            'Y:ENSG000002:E3.8-E8.1|ENSG000002:E7.1-E8.1']
    assert ('ENSG000001','I5') not in mh._event_exons(rows[0])
    assert ('ENSG000001','I5') in mh._event_exons(rows[0],include_introns=True)
    cells=['very_long_sample_identifier_for_cell_state.patient_1',
           'very_long_sample_identifier_for_cell_state.patient_2']
    matrix=pd.DataFrame([[1.,-1.],[-1.,1.],[0.,.2]],index=rows,columns=cells)
    cache=tmp_path/'cache.h5ad';mh._write_heatmap_cache(cache,matrix,['up group']*3,['group']*2,cells)
    original=cache.read_bytes();pdf=tmp_path/'red.pdf'
    mh.render_heatmap_from_cache(str(cache),str(pdf),load_go_terms=False,dpi=70,
        label_events=listed,species='human',label_shared_introns=True,callout_color='#d62728',
        sample_labels={x:x for x in cells},show_pattern_labels=True)
    audit=pd.read_csv(tmp_path/'red_label_events.tsv',sep='\t')
    assert set(audit.event)==set(rows[:2])
    assert 'ENSG000001:I5' in audit.shared_exons.iloc[0]
    assert cache.read_bytes()==original
    doc=pymupdf.open(pdf);text=doc[0].get_text()
    assert all(x in text for x in cells)
    spans=[s for b in doc[0].get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
    red=[s for s in spans if s['color']==int('d62728',16)]
    assert len(red)==2
    assert all(s['bbox'][1]>=-.8 and s['bbox'][3]<=doc[0].rect.height+.8 for s in spans)
