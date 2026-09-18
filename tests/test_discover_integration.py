"""Scientific gates and integrated API behavior for both shared adapters."""
from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
from test_grn_web_release import uploaded, seed_comparisons
from altanalyze3.components.cellHarmony import discover_integration as d


class Evidence:
    current_contrast='case';states=['AT2']
    state_expression_scale='mean stored RNA expression in the cell state'
    def __init__(self):
        self.rows={'rna':[dict(gene='G',log2fc=2.,pval=.001,fdr=.01)],
                   'grn':[dict(gene='T|G',log2fc=1.,pval=.001,fdr=.02)],
                   'grn_tf':[dict(gene='T',log2fc=-1.,pval=.001,fdr=.01)],
                   'lipid':[dict(gene='PC(1)',log2fc=2.,pval=.001,fdr=.01),dict(gene='PC(2)',log2fc=4.,pval=.001,fdr=.02),dict(gene='PC(3)',log2fc=6.,pval=.03,fdr=.9)],
                   'metabolite':[dict(gene='Cholesterol',log2fc=-2.,pval=.001,fdr=.01)]}
    def comparison(self,c):return c or self.current_contrast
    def available(self,c,m,s):return self.comparison(c)=='case' and m in self.rows
    def differentials(self,c,s,m):return self.rows.get(m,[])
    def edge_scores(self,s):return {'T|G':.1,'S|G':99.}
    def arm_expression(self,c,s,fs):return {'T':.2},{'T':1.0}
    def marker_features(self,*a):return ['G']
    def marker_rows(self,s):return [dict(gene=g,log2fc=2.) for g in self.marker_features(s)]
    def state_expression(self,s,fs):return {'T':1.,'S':0.}
    def features(self,m):return [r['gene'] for r in self.rows.get(m,[])]


def test_own_edge_gate_and_either_arm_expression():
    e=Evidence();r=d.network(e,cell_state='AT2')
    assert [x['source'] for x in r['edges']]==['T'] # strongest static edge has no differential
    t=next(n for n in r['nodes'] if n['id']=='T')
    assert t['log2fc'] is None and t['activity_log2fc']==-1 and t['expression']==1
    assert not d.network(e,cell_state='AT2',min_expression=1)['edges']
    assert not d.network(e,cell_state='AT2',min_fold=3)['edges']
    e.rows['grn'][0]['fdr']=.5
    assert not d.network(e,cell_state='AT2',significance='fdr')['edges']
    assert d.network(e,cell_state='AT2')['edges']  # stored calls are not filtered again
    assert d.network(e,cell_state='AT2',significance='pval')['edges']
    assert not d.network(e,cell_state='AT2',contrast='unrelated')['available']


def diagrams():
    return {'diagrams':{'TEST':{'name':'Test lipid pathway','classes':['PC'],'width':300,'height':200,
        'nodes':[dict(id='g',label='G',type='GeneProduct',x=20,y=20,w=20,h=10),
                 dict(id='l',label='Phosphatidylcholine',type='Metabolite',x=60,y=20,w=20,h=10),
                 dict(id='m',label='Cholesterol',type='Metabolite',x=100,y=20,w=20,h=10)],'edges':[],'labels':[]}}}


def test_native_pathway_reported_class_mean_and_metabolite(monkeypatch):
    monkeypatch.setattr(d,'_diagrams',diagrams);monkeypatch.setattr(d,'_diagram_synonyms',lambda:{'PC':['phosphatidylcholine']})
    e=Evidence();r=d.pathway(e,id='TEST',cell_state='AT2')
    assert [n['log2fc'] for n in r['nodes']]==[2.,4.,None]
    assert r['nodes'][1]['link_feature']=='PC(3)'
    assert 'Test lipid pathway' in r['interpretation'] and 'Phosphatidylcholine' in r['interpretation']
    assert d.pathway(e,id='TEST',cell_state='AT2',max_fdr=0,significance='fdr')['nodes']==r['nodes']
    m=d.pathway(e,id='TEST',cell_state='AT2',modality='metabolite')
    assert [n['log2fc'] for n in m['nodes']]==[2.,None,-2.]
    ranks=d.pathways(e,cell_state='AT2',modality='metabolite')['pathways']
    assert ranks[0]['n_genes_painted']==1 and ranks[0]['n_lipids_painted']==1
    del e.rows['rna']
    assert not d.pathway(e,id='TEST',cell_state='AT2')['available']
    assert 'both differential expression and lipid imputed differentials' in d.pathways(e,cell_state='AT2')['note']


def test_network_ranks_connected_targets_before_limit_and_interprets():
    e=Evidence()
    e.rows['rna'].insert(0,dict(gene='UNCONNECTED',log2fc=100.,pval=.001,fdr=.01))
    r=d.network(e,cell_state='AT2',limit=1)
    assert {n['id'] for n in r['nodes']}=={'G','T'}
    assert {n['id'] for n in r['nodes']}=={e[k] for e in r['edges'] for k in ('source','target')}
    assert 'T → G' in r['interpretation'] and 'UNCONNECTED' not in r['interpretation']
    e.marker_features=lambda *args:['UNCONNECTED','G']
    assert d.network(e,cell_state='AT2',source='marker',limit=1)['edges']
    e.rows['grn']=[]
    assert not d.network(e,cell_state='AT2')['nodes']


def test_marker_network_is_independent_of_differentials():
    e=Evidence()
    e.available=lambda *a:pytest.fail('Marker view must not require differential runs')
    e.differentials=lambda *a:pytest.fail('Marker view must not read differential calls')
    e.arm_expression=lambda *a:pytest.fail('Marker view must not use comparison arms')
    r=d.network(e,source='marker',contrast='nonexistent')
    assert r['cell_state']=='AT2' and r['source']=='marker'
    assert [(x['source'],x['target']) for x in r['edges']]==[('T','G')]
    assert r['edges'][0]['log2fc'] is None
    assert next(n for n in r['nodes'] if n['id']=='G')['log2fc']==2
    assert not d.network(e,source='marker',min_score=.2)['edges']
    assert not d.network(e,source='marker',min_expression=1)['edges']
    assert not d.network(e,source='marker',gene_fold=8)['edges']
    assert not d.network(e,source='marker',features=['NOT_A_MARKER'])['edges']
    assert not d.network(e,source='marker',cell_state='invalid')['available']


def test_marker_network_excludes_nonpositive_markers_even_without_fold_threshold():
    e = Evidence()
    e.marker_rows = lambda state: [dict(gene=g, log2fc=fold) for g, fold in
                                  [('G', 2.), ('DOWN', -8.), ('ZERO', 0.), ('T', -3.)]]
    e.edge_scores = lambda state: {'T|G': .1, 'T|DOWN': .5, 'T|ZERO': .7}
    for threshold in (1, 1.2):
        result = d.network(e, source='marker', gene_fold=threshold)
        assert {n['id'] for n in result['nodes']} == {'T', 'G'}
        assert next(n for n in result['nodes'] if n['id'] == 'G')['log2fc'] == 2.
        assert next(n for n in result['nodes'] if n['id'] == 'T')['log2fc'] is None
        assert [(x['source'], x['target']) for x in result['edges']] == [('T', 'G')]
    assert not d.network(e, source='marker', features=['DOWN'])['nodes']


def test_uploaded_marker_network_without_comparisons(uploaded):
    u=uploaded
    path=u.root/'markers.tsv'
    pd.DataFrame([{'Gene':'TARGET1','Fold':2.,'cluster':'AT1'},
                  {'Gene':'TARGET3','Fold':-1.,'cluster':'AT2'}]).to_csv(path,sep='\t',index=False)
    u.store.update_job(u.job,marker_analysis={'enabled':True,'markers_tsv':str(path)})
    base=f'/api/jobs/{u.job}/integrated'
    r=u.client.get(base+'/network',params={'source':'marker','cell_state':'AT1'}).json()
    assert r['edges'] and r['source']=='marker',r
    assert {n['id'] for n in r['nodes']}=={'TF1','TARGET1'}
    assert all(e['log2fc'] is None for e in r['edges'])
    tsv=u.client.get(base+'/network.tsv',params={'source':'marker','cell_state':'AT1'})
    assert 'target_marker_log2fc' in tsv.text and 'TF1\tTARGET1' in tsv.text
    assert 'edge_log2fc' not in tsv.text
    negative = u.client.get(base+'/network', params={'source':'marker', 'cell_state':'AT2', 'gene_fold':1}).json()
    assert negative['status'] == 'no_markers' and not negative['nodes']
    assert not u.client.get(base+'/network',params={'cell_state':'AT1'}).json()['available']


def test_bundle_marker_network_uses_selected_state_without_integrated_manifest():
    from altanalyze3.components.cellHarmony.webapp.integration_data import BundleIntegrationData
    store=SimpleNamespace(features=['T|G'],stats_mean=np.array([[.1,.8]]))
    ds=SimpleNamespace(states=['A','B'],sv={},symbols=['T','G'],
                       stats_mean=np.array([[0.,2.],[1.,3.]]),
                       modality=lambda _:store,
                       markers=lambda:[dict(gene='G',cluster=s,fold=2.,p=.01) for s in ('A','B')])
    adapter=BundleIntegrationData(ds,{}, {})
    assert not d.network(adapter,source='marker',cell_state='A')['edges']
    result=d.network(adapter,source='marker',cell_state='B')
    assert result['edges'][0]['score']==.8
    assert next(n for n in result['nodes'] if n['id']=='T')['expression']==2.


def test_pathway_menu_excludes_maps_without_regulated_metabolites(monkeypatch):
    def maps():
        data=diagrams()
        data['diagrams']['GENE_ONLY']=dict(data['diagrams']['TEST'],name='Gene only',nodes=[data['diagrams']['TEST']['nodes'][0]])
        return data
    monkeypatch.setattr(d,'_diagrams',maps)
    monkeypatch.setattr(d,'_diagram_synonyms',lambda:{'PC':['phosphatidylcholine']})
    e=Evidence()
    assert [r['id'] for r in d.pathways(e,cell_state='AT2')['pathways']]==['TEST']
    assert not d.pathway(e,id='GENE_ONLY',cell_state='AT2')['nodes']
    e.rows['lipid']=[]
    result=d.pathways(e,cell_state='AT2')
    assert result['available'] and result['pathways']==[]
    assert 'No pathway contains a regulated' in result['note']


def test_arm_mean_cache_expires_when_evidence_is_replaced(tmp_path):
    from altanalyze3.components.cellHarmony.webapp.integration_data import BundleIntegrationData
    path=tmp_path/'arms.npz'
    def write(value):
        np.savez(path,states=np.array(['AT2']),genes=np.array(['G']),
                 case_mean_log2cp10k=np.array([[value]]),control_mean_log2cp10k=np.array([[.2]]))
    write(1.)
    ds=object.__new__(BundleIntegrationData);ds.root=tmp_path
    ds.manifest={'comparisons':{'case':{'arm_means':'arms.npz'}}};ds.current_contrast='case'
    assert ds.arm_expression('case','AT2',['G'])==({'G':1.},{'G':.2})
    write(2.)
    assert ds.arm_expression('case','AT2',['G'])==({'G':2.},{'G':.2})


def test_uploaded_integrated_endpoints_and_missing_pair(uploaded,monkeypatch):
    u=uploaded;runs=seed_comparisons(u);base=f'/api/jobs/{u.job}'
    r=u.client.get(base+'/integrated/network',params={'cell_state':'AT1','features':'TARGET1'})
    assert r.status_code==200,r.text
    assert r.json()['edges'],r.text
    assert 'TF1 → TARGET1' in r.json()['interpretation']
    assert r.json()['edges'][0]['log2fc']==1.5
    tsv=u.client.get(base+'/integrated/network.tsv',params={'cell_state':'AT1','features':'TARGET1'})
    assert 'TF_expression' in tsv.text and 'TF1\tTARGET1' in tsv.text
    for question,kind in [('Show common GRN impacts for marker features in AT1','integrated_network'),('Show common lipid pathway impacts in AT1','integrated_pathway')]:
        result=u.client.post(base+'/chat',json={'question':question}).json()
        assert result['plot']['kind']==kind,result
    runs.pop('run-rna');u.store.update_job(u.job,differential_history=runs)
    result=u.client.get(base+'/integrated/network',params={'cell_state':'AT1'}).json()
    assert result['missing_modalities']==['rna']
    assert 'both differential expression and grn imputed differentials' in result['note']
