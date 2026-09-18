from types import SimpleNamespace
from importlib import import_module
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from altanalyze3.components.cellHarmony.webapp import cross_modal as cm


def cache(x, genes, reverse=False):
    obs = pd.DataFrame({'state': np.repeat(['A', 'B', 'C', 'D'], 6),
                        'Library': np.tile(['sample1', 'sample2', 'sample3'], 8)},
                       index=[f'c{i}' for i in range(24)])
    a = ad.AnnData(sparse.csr_matrix(x, dtype=float), obs=obs, var=pd.DataFrame(index=genes))
    if reverse:
        a = a[::-1].copy()
    return dict(adata=a, var_names=a.var_names.to_numpy(), obs_names=a.obs_names.to_numpy(),
                populations=a.obs.state.to_numpy())


def inputs():
    values = np.repeat([1., 2., 4., 8.], 6)
    rna = cache(np.column_stack([values, values, values, np.ones(24)]), ['TF1','MME','NCAM1','CONST'])
    tf = cache(np.column_stack([10-values, np.ones(24)]), ['TF1','CONST'], reverse=True)
    adt = cache(np.column_stack([values*3, 12-values, values]), ['CD10','CD56','unknown'], reverse=True)
    return rna, tf, adt


def test_pairs_align_by_cell_id_and_constant_omission():
    rna, tf, _ = inputs()
    data = cm.compute(rna, tf, 'grn_tf')
    assert data['constant'] == 1
    pair, = data['pairs']
    assert pair['rho'] == -1
    assert pair['x'] == pytest.approx([1,2,4,8])
    assert pair['y'] == pytest.approx([9,8,6,2])
    assert pair['n_cells'] == [6]*4


def test_adt_curated_lookup_including_multi_subunit():
    rna, _, adt = inputs()
    data = cm.compute(rna, adt, 'adt')
    assert {(p['feature'],p['gene'],p['rho']) for p in data['pairs']} == {('CD10','MME',1),('CD56','NCAM1',-1)}
    assert data['skipped'] == ['unknown']
    mapping = cm.adt_map()
    assert mapping['cd4_rpa.t4'] == ['CD4']
    assert mapping['cd8'] == ['CD8A','CD8B']
    assert mapping['cd45ra'] == ['PTPRC']


def test_rna_alias_and_ambiguous_gene():
    rna, tf, _ = inputs()
    rna['adata'].var['gene_symbols'] = ['TF1','MME','NCAM1','CONST']
    rna['var_names'] = np.array(['ENSG1','ENSG2','ENSG3','ENSG4'])
    assert cm.paired_features(rna, tf, 'grn_tf')[0][0][1] == 'TF1'
    rna['adata'].var.iloc[1,0] = 'TF1'
    assert cm.paired_features(rna, tf, 'grn_tf')[1] == ['TF1']


def test_baseline_chat_does_not_require_differential(monkeypatch):
    rna, tf, adt = inputs()
    web = import_module('altanalyze3.components.cellHarmony.webapp.app')
    monkeypatch.setattr(web, '_get_expression_cache', lambda app, meta, modality: {'rna':rna,'grn_tf':tf,'adt':adt}[modality])
    cm._CACHE.clear()
    app = SimpleNamespace()
    meta = {'job_id':'test', 'species':'human'}
    result = cm.answer_if_requested(app, meta, cm.examples(['grn_tf'])[0])
    assert result['status'] == 'ok'
    assert result['plot']['pairs'][0]['rho'] == -1
    assert result['reading']['source'] == 'expression'
    assert len(cm._CACHE) == 1
    again = cm.answer_if_requested(app, meta, 'Correlate TF activity with gene expression')
    assert again['table'] == result['table']
    assert len(cm._CACHE) == 1
    adt_result = cm.answer_if_requested(app, meta, 'Correlate CD56 ADT and gene expression')
    assert adt_result['table']['rows'][0]['gene'] == 'NCAM1'
    assert len(adt_result['table']['rows']) == 1


def test_insufficient_groups_and_duplicate_ids():
    rna, tf, _ = inputs()
    with pytest.raises(ValueError, match='three groups'):
        cm.compute(rna, tf, 'grn_tf', state='A')
    rna['obs_names'][0] = rna['obs_names'][1]
    with pytest.raises(ValueError, match='unique'):
        cm.compute(rna, tf, 'grn_tf')


def test_router_does_not_take_differential_queries():
    assert cm.requested_modality('Which TF activity changed in AT1?') is None
    assert cm.requested_modality('Show the regulatory network') is None
    assert cm.requested_modality('Which TFs have activity discordant with gene expression?') == 'grn_tf'


def test_group_controls_include_all_cell_states():
    web = import_module('altanalyze3.components.cellHarmony.webapp.app')
    obs = pd.DataFrame({'state':pd.Categorical([f's{i}' for i in range(86)])})
    a = ad.AnnData(np.zeros((86,1)),obs=obs)
    c = dict(adata=a,cluster_key='state',populations=obs.state.astype(str).to_numpy())
    variables = web._groupable_columns(c)
    assert variables[0]['field'] == 'state'
    assert variables[0]['n'] == 86


def test_repeated_float32_predictions_keep_ties():
    obs=pd.DataFrame({'state':np.repeat(['A','B','C','D'],[7,11,15,22])})
    vals=np.full((len(obs),1),np.float32(0.1234567))
    a=ad.AnnData(vals,obs=obs,var=pd.DataFrame(index=['TF1']))
    c=dict(adata=a,var_names=a.var_names.to_numpy(),obs_names=a.obs_names.to_numpy(),populations=obs.state.to_numpy())
    data=cm.compute(c,c,'grn_tf')
    assert data['pairs'] == [] and data['constant'] == 1


def test_redundant_marker_folds_are_used(tmp_path):
    web=import_module('altanalyze3.components.cellHarmony.webapp.app')
    marker=tmp_path/'test_markers.tsv'
    pd.DataFrame({'Gene':['A'],'Fold':[2.],'cluster':['state']}).to_csv(marker,sep='\t',index=False)
    pd.DataFrame({'Gene':['B','C'],'Fold':[-1.,0.],'cluster':['state','state']}).to_csv(tmp_path/'test_redundant_markers.tsv',sep='\t',index=False)
    edges=tmp_path/'edges.tsv'
    pd.DataFrame({'Symbol1':['A','B','C'],'Symbol2':['B','C','D']}).to_csv(edges,sep='\t',index=False)
    meta={'marker_analysis':{'markers_tsv':str(marker),'networks':[{'population':'state','tsv':str(edges)}]}}
    p=web._build_marker_network_payload(meta,'state')
    values={e['data']['id']:e['data']['log2fc'] for e in p['elements'] if 'source' not in e['data']}
    assert values == {'A':2.,'B':-1.,'C':0.,'D':None}


def test_named_state_uses_categorical_samples():
    rna, tf, _ = inputs()
    rna['populations'][:] = 'A'
    rna['adata'].obs['Library'] = pd.Categorical(np.repeat(['S1','S2','S3'],8))
    data=cm.compute(rna,tf,'grn_tf',state='A')
    assert data['n_units'] == 3
    assert data['unit'] == 'Library averages within A'
    assert data['pairs'][0]['rho'] == -1
    assert data['pairs'][0]['n_cells'] == [8,8,8]


def test_all_pairs_remain_available_and_weak_queries_rank_near_zero(monkeypatch):
    rna, tf, _ = inputs()
    web=import_module('altanalyze3.components.cellHarmony.webapp.app')
    monkeypatch.setattr(web,'_get_expression_cache',lambda app,meta,modality:rna if modality=='rna' else tf)
    pairs=[dict(feature=f'TF{i}',gene=f'TF{i}',rho=float(rho),n_units=4,
                x=[1,2,3,4],y=[4,2,3,1],labels=['A','B','C','D'],n_cells=[6]*4)
           for i,rho in enumerate(np.linspace(-1,1,217))]
    monkeypatch.setattr(cm,'compute',lambda *args:dict(pairs=pairs,skipped=[],constant=0,unit='cell-state averages',matched_observations=24,n_units=4))
    cm._CACHE.clear()
    weak=cm.answer_if_requested(None,{'job_id':'all-pairs'},'Which TFs have activity uncorrelated with gene expression?')
    assert len(weak['table']['rows']) == len(weak['plot']['pairs']) == 217
    assert weak['table']['rows'][0]['rho'] == 0
    assert weak['result_controls']['relationship'] == 'near_zero'
    assert weak['result_controls']['sort_direction'] == 'asc'
    negative=cm.answer_if_requested(None,{'job_id':'all-pairs'},'Which TFs have discordant activity and gene expression?')
    assert len(negative['table']['rows']) == 217
    assert negative['result_controls']['relationship'] == 'negative'
    assert negative['table']['rows'][0]['rho'] == -1
    assert any(row['rho']>0 for row in negative['table']['rows'])
    assert cm.requested_modality('Which ADTs are least correlated with gene expression?') == 'adt'
