import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from altanalyze3.components.cellHarmony.webapp import cross_pathways as C


@pytest.fixture
def library(monkeypatch):
    def node(label, kind='GeneProduct'):
        return dict(id=label,label=label,type=kind,x=10,y=10,w=20,h=15,ensembl='')
    diagrams = {
        'A':dict(id='A',name='Broad',nodes=[node('G1'),node('G1'),node('G2'),node('PC','Metabolite'),node('Known','Metabolite')]),
        'B':dict(id='B',name='RNA rich',nodes=[node('G'+str(i)) for i in range(1,101)]),
        'C':dict(id='C',name='Unmapped',nodes=[node('Unknown 0001','Metabolite')]),
    }
    monkeypatch.setattr(C,'_diagrams',lambda:dict(diagrams=diagrams))
    monkeypatch.setattr(C,'_diagram_synonyms',lambda:{'PC':['PC','phosphatidylcholine']})
    monkeypatch.setattr(C,'antibody_mapping',lambda species:{'cd-test':['G1']})
    C.pathway_index.cache_clear()
    yield diagrams
    C.pathway_index.cache_clear()


def features(*names):
    return {n:dict(statistic='test',value=.5) for n in names}


def test_unique_features_balanced_ranking_and_missing_modalities(library):
    f=dict(rna=features(*['G'+str(i) for i in range(1,101)]),adt=features('CD-test'),
           metabolite=features('Known','Unknown 0001'),lipid=features('PC 16:0','PC 18:0'),grn_tf=features('G1'),
           grn=features('G1|G2','G1|OUTSIDE'))
    labels={m:m for m in f}
    r=C.rank(f,labels)
    assert [v['id'] for v in r['rows']]==['A','B']
    a=r['rows'][0]
    assert a['rna']==2 and a['lipid']==2 and a['adt']==1 and a['grn']==1 and a['metabolite']==1
    assert a['modality_count']==6
    assert a['balanced_coverage']==pytest.approx((.02+1+1+1+1+1)/6)
    assert next(m for m in r['modalities'] if m['id']=='rna')['normalization_max']==100
    # Counts count original features, not diagram copies or mapped partner genes.
    d=C.diagram(dict(r,features=f,cell_state='HSC-1',source='marker',comparison='',species='human'),'A')
    assert d['cross_modal'] and d['available']
    assert [n['hits'] for n in d['nodes'] if n['label']=='G1'][0]==[n['hits'] for n in d['nodes'] if n['label']=='G1'][1]
    assert C.diagram(dict(r),'C')['available'] is False


def test_question_source_and_exact_state(monkeypatch,library):
    import pandas as pd
    frame=pd.DataFrame(dict(Gene=['G1','G2'],cluster=['HSC-1','HSC-10'],rho=[.5,.8]))
    monkeypatch.setattr(C,'marker_tables',lambda m:({'rna':frame},[]))
    meta={'modalities':{'available':[{'id':'rna','label':'RNA'}]}}
    r=C.answer_if_requested(None,meta,'What pathways have the best cross-modality HSC-1 representation?')
    assert r['status']=='ok' and r['plot']['source']=='marker' and r['plot']['contrast']==''
    assert r['table']['rows'][0]['rna']==1
    assert C.answer_if_requested(None,meta,'What pathways have the best cross-modality representation?')['status']=='clarify'
    assert C.answer_if_requested(None,meta,'What are markers of HSC-1?') is None
    assert C.requested('Cross-modality pathways for HSC-1 contrasts')


def test_unidentified_and_stereoisomer_mapping_is_conservative(library):
    idx=C.pathway_index()['A']
    assert not C.feature_nodes('Unknown 0001','metabolite',idx)
    assert not C.feature_nodes('D-Known','metabolite',idx)
    assert C.feature_nodes('Known','metabolite',idx)


def test_differential_uses_requested_contrast_and_no_new_threshold(monkeypatch,library):
    calls=[]
    class Dataset:
        current_contrast='run-a'
        def deg_manifest(self):return {'comparisons':[{'id':'run-a','comparison':'A vs B'}]}
        def available(self,contrast,modality,state):
            calls.append((contrast,modality,state));return modality=='rna'
        def differentials(self,contrast,state,modality):
            return [{'feature_key':'G1','gene':'Display name','log2fc':.01,'fdr':.8}]
    from altanalyze3.components.cellHarmony.webapp import integration_data
    monkeypatch.setattr(integration_data,'integration_data',lambda a,m:Dataset())
    meta={'modalities':{'available':[{'id':'rna'},{'id':'adt'}]}}
    result=C.collect(None,meta,'HSC-1','differential','run-a')
    assert result['features']['rna']['G1']['value']==.01
    assert result['rows'][0]['rna']==1 and result['rows'][0]['adt']==0
    assert calls==[('run-a','rna','HSC-1'),('run-a','adt','HSC-1')]
    with pytest.raises(ValueError): C.collect(None,meta,'HSC-1','differential','wrong')


def test_modality_scales_and_strongest_node_scores(library):
    f=dict(rna={'G1':dict(statistic='reported log2FC',value=-3.),'G2':dict(statistic='reported log2FC',value=1.)},
           lipid={'PC 16:0':dict(statistic='reported log2FC',value=-.02),'PC 18:0':dict(statistic='reported log2FC',value=.04)},
           adt={})
    ranked=C.rank(f,{m:m for m in f})
    result=dict(ranked,features=f,cell_state='HSC-1',source='differential',comparison='case vs control')
    d=C.diagram(result,'A'); scales={m['id']:m['scale'] for m in d['modalities']}
    assert scales['rna']['minimum']==-3 and scales['rna']['maximum']==3
    assert scales['lipid']['minimum']==-.04 and scales['lipid']['maximum']==.04
    assert scales['lipid']['n_features']==2  # not diagram node count
    assert scales['adt']['available'] is False
    lipid=next(n for n in d['nodes'] if n['label']=='PC')
    assert lipid['modality_values']['lipid']['value']==.04 and len(lipid['hits'])==2
    result['source']='marker'
    f['rna']['G1']['value']=.5
    d=C.diagram(result,'A'); scales={m['id']:m['scale'] for m in d['modalities']}
    assert scales['rna']['minimum']==0
    assert scales['rna']['statistic']=='MarkerFinder r'
