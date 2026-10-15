"""Display filters preserve parent coordinates/palette without changing cells/counts."""
import copy
from importlib import import_module
from types import SimpleNamespace
import anndata as ad
import numpy as np
import pandas as pd
import pytest

W=import_module('altanalyze3.components.cellHarmony.webapp.app')


def fixture(monkeypatch,mode='supervised'):
    obs=pd.DataFrame({'state':['A','B','C','D','E','F'], 'Library':['one','one','two','two','two','two'],
                      'author':['x','y','z','x','y','z'], 'pseudotime':[0.,1.,2.,3.,4.,5.]},index=[f'c{i}' for i in range(6)])
    a=ad.AnnData(np.arange(6,dtype=np.float32).reshape(-1,1),obs=obs,var=pd.DataFrame(index=['G']))
    a.obsm['alternative']=np.column_stack([np.arange(6)*100.,np.arange(6)*-100.])
    cache=dict(adata=a,cluster_key='state',obs_names=a.obs_names.to_numpy(),var_names=a.var_names.to_numpy(),
        populations=obs.state.to_numpy(),sample_labels=obs.Library.to_numpy(),sample_field='Library',
        umap_x=np.array([-3.,0.,2.,10.,15.,60.]),umap_y=np.array([5.,4.,3.,2.,1.,-30.]),
        obs_filter_values={'Library':obs.Library.to_numpy(),'state':obs.state.to_numpy()})
    ref=ad.AnnData(np.zeros((2,1)),obs=pd.DataFrame({'refstate':['R','A']},index=['r1','r2']))
    ref.obsm['X_umap']=np.array([[-40.,-100.],[80.,40.]])
    monkeypatch.setattr(W,'_get_expression_cache',lambda *a,**k:cache)
    monkeypatch.setattr(W,'_load_reference_adata',lambda *a:ref)
    return cache,dict(analysis_mode=mode,cluster_key='state',reference_cluster_key='refstate')


@pytest.mark.parametrize('compact',[False,True])
def test_filters_never_change_parent_ranges_or_population_colors(monkeypatch,compact):
    cache,meta=fixture(monkeypatch)
    full=W._build_umap_payload(None,meta,compact=compact)
    for library in ['one','two','absent']:
        subset=W._build_umap_payload(None,meta,compact=compact,display_filters=[('Library',[library])])
        for field in ['coordinate_bounds','population_order','population_colors']:
            assert subset[field]==full[field]
        assert subset['n_cells_selected']==int((cache['sample_labels']==library).sum())
    assert full['coordinate_bounds']=={'x':[-40.,80.],'y':[-100.,40.]}
    assert full['population_order']==['R','A','B','C','D','E','F']


def test_color_annotation_and_embedding_are_independent_of_filters(monkeypatch):
    _,meta=fixture(monkeypatch)
    full=W._build_umap_payload(None,meta,color_by='author',coords_key='alternative')
    filtered=W._build_umap_payload(None,meta,color_by='author',coords_key='alternative',display_filters=[('Library',['one'])])
    assert full['population_order']==['x','y','z'] and full['population_colors']==filtered['population_colors']
    assert full['coordinate_bounds']==filtered['coordinate_bounds']=={'x':[0.,500.],'y':[-500.,0.]}
    # Alternative coordinate spaces never inherit an unrelated atlas extent.
    assert not full['reference']


def test_metadata_axes_and_expression_use_the_complete_parent(monkeypatch):
    _,meta=fixture(monkeypatch)
    options=dict(x_field='pseudotime',y_field='__umap_2')
    categorical=W._build_umap_payload(None,meta,**options)
    full=W._build_expression_payload(None,meta,'G',view='umap',**options)
    subset=W._build_expression_payload(None,meta,'G',view='umap',display_filters=[('Library',['one'])],**options)
    assert full['coordinate_bounds']==subset['coordinate_bounds']==categorical['coordinate_bounds']=={'x':[0.,5.],'y':[-30.,5.]}
    assert [(p['barcode'],p['value'],p['x'],p['y']) for p in subset['umap']]==[('c0',0.,0.,5.),('c1',1.,1.,4.)]


def test_both_branches_preserve_their_native_coordinate_spaces_and_union_frequencies(monkeypatch):
    cache,meta=fixture(monkeypatch,'both')
    cache['adata'].obs['supervised_state']=cache['populations']
    cache['adata'].obs['umap_supervised_x']=cache['umap_x']
    cache['adata'].obs['umap_supervised_y']=cache['umap_y']
    cache['branch_embeddings']={'X_umap_supervised':np.column_stack([cache['umap_x'],cache['umap_y']]),
                                'X_umap_unsupervised':np.column_stack([np.arange(6),-np.arange(6)])}
    left=W._build_umap_payload(None,meta,color_by='supervised_state',x_field='umap_supervised_x',y_field='umap_supervised_y',display_filters=[('Library',['one'])])
    right=W._build_umap_payload(None,meta,coords_key='X_umap_unsupervised',display_filters=[('Library',['one'])])
    assert left['coordinate_bounds']=={'x':[-40.,80.],'y':[-100.,40.]}
    assert right['coordinate_bounds']=={'x':[0.,5.],'y':[-5.,0.]}
    assert sum(p['cells'] for p in left['frequency_counts'])==2
    full=W._build_umap_payload(None,meta,coords_key='X_umap_unsupervised')
    assert full['population_colors']==right['population_colors']


def test_nonfinite_rows_do_not_poison_parent_bounds_or_union_counts(monkeypatch):
    cache,meta=fixture(monkeypatch,'both')
    cache['umap_x'][-1]=np.nan
    full=W._build_umap_payload(None,meta)
    assert full['coordinate_bounds']=={'x':[-3.,15.],'y':[1.,5.]}
    assert full['n_points_drawn']==5 and sum(r['cells'] for r in full['frequency_counts'])==6


def test_parent_context_cached_without_serializing_all_parent_points(monkeypatch):
    cache,meta=fixture(monkeypatch)
    full=W._build_umap_payload(None,meta)
    parent=next(iter(cache['parent_plot_contexts'].values()))
    W._build_umap_payload(None,meta,display_filters=[('Library',['one'])])
    assert next(iter(cache['parent_plot_contexts'].values())) is parent
    assert len(cache['parent_plot_contexts'])==1
    assert set(parent)=={'population_order','population_colors','coordinate_bounds'}


def test_pdf_uses_parent_limits_and_the_same_population_palette(monkeypatch):
    _,meta=fixture(monkeypatch)
    full=W._build_umap_payload(None,meta)
    subset=W._build_umap_payload(None,meta,display_filters=[('Library',['one'])])
    captures=[]
    original=W.plt.close
    monkeypatch.setattr(W.plt,'close',lambda fig:captures.append(fig) if hasattr(fig,'axes') else None)
    try:
        for payload in [full,subset]:
            assert W._render_umap_pdf(payload,'cluster').read(4)==b'%PDF'
        assert captures[0].axes[0].get_xlim()==captures[1].axes[0].get_xlim()
        assert captures[0].axes[0].get_ylim()==captures[1].axes[0].get_ylim()
        assert W._render_umap_pdf(subset,'frequency').read(4)==b'%PDF'
        bars=captures[-1].axes[0].patches
        assert {W.matplotlib.colors.to_hex(b.get_facecolor()) for b in bars}=={full['population_colors']['A'],full['population_colors']['B']}
    finally:
        for fig in captures:original(fig)


@pytest.mark.parametrize('count',[0,1,4,5,17,70])
def test_plot_palette_matches_the_existing_browser_palette(count):
    import json,shutil,subprocess
    from pathlib import Path
    node=shutil.which('node')
    if not node:pytest.skip('Node is needed for browser palette parity')
    order=[f'state{i}' for i in range(count)]
    script=r'''
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const input=JSON.parse(fs.readFileSync(0,'utf8')),code=fs.readFileSync(input.path,'utf8'),ctx=vm.createContext({Map,Math});
const start=code.indexOf('const PAIRED_COLOR_STOPS = [');
vm.runInContext(code.slice(start,code.indexOf('];',start)+2),ctx);
vm.runInContext(code.slice(code.indexOf('function interpolatePairedColor('),code.indexOf('function buildStableUmapPopulationOrder(')),ctx);
const colors=Object.fromEntries(Array.from(ctx.buildReferencePreviewColorMap(input.order),([name,color])=>[name,color.startsWith('#')?color:'#'+color.match(/\d+/g).map(v=>Number(v).toString(16).padStart(2,'0')).join('')]));
assert.deepEqual(colors,input.colors);
'''
    data=dict(order=order,colors=W._plot_population_colors(order),path=str(Path(W.__file__).parent/'static/app.js'))
    r=subprocess.run([node,'-e',script],input=json.dumps(data),capture_output=True,text=True)
    assert r.returncode==0,r.stderr
