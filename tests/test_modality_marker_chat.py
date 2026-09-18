import numpy as np
import pandas as pd
import pytest
from altanalyze3.components.cellHarmony.webapp import modality_markers as M


@pytest.fixture
def meta(tmp_path):
    tables={
        'rna':[('IFIT2','HSC-1',.5082,4.7),('SECOND','HSC-1',.2,100),('OTHER','HSC-2',.9,3)],
        'adt':[('CD73','HSC-1',.15,1.7)],
        'metabolite':[('Unknown 0400','HSC-1',.2527,-2.854),('Known metabolite','HSC-1',.1,1),('constant','HSC-1',np.nan,10),('negative','HSC-1',-.9,10)],
    }
    m={'job_id':'fixture','modalities':{'available':[{'id':k,'label':k} for k in [*tables,'grn']]},'marker_analysis_by_modality':{}}
    for mod,values in tables.items():
        path=tmp_path/f'{mod}_markers.tsv'
        pd.DataFrame(values,columns=['Gene','cluster','rho','Fold']).to_csv(path,sep='\t',index=False)
        m['marker_analysis_by_modality'][mod]={'markers_tsv':str(path)}
    return m


def test_best_across_modalities_uses_positive_marker_correlation(meta):
    result=M.answer_if_requested(None,meta,'What is the best modality marker of HSC-1')
    assert result['status']=='ok'
    assert result['best']['feature']=='IFIT2'
    assert len(result['table']['rows'])==5
    assert [r['feature'] for r in result['best_by_modality']]==['IFIT2','Unknown 0400','CD73']
    assert result['provenance']['differential_required'] is False
    assert result['plot']['kind']=='modality_markers'
    assert next(r for r in result['coverage'] if r['modality']=='grn')['n_positive_markers']==0


def test_named_and_single_modality_questions(meta):
    result=M.answer_if_requested(None,meta,'What is the best named modality marker of HSC-1?')
    assert all(r['annotation']=='Named' for r in result['table']['rows'])
    assert next(r for r in result['best_by_modality'] if r['modality_id']=='metabolite')['feature']=='Known metabolite'
    result=M.answer_if_requested(None,meta,'Which is the best ADT marker of HSC-1?')
    assert result['best']['feature']=='CD73'
    assert len(result['table']['rows'])==1


def test_state_is_required_and_matches_whole_name(meta):
    for q in ['What is the best modality marker?', 'What is the best modality marker of HSC-10?', 'Best modality marker of HSC-1 and HSC-2?']:
        assert M.answer_if_requested(None,meta,q)['status']=='clarify'
    assert M.answer_if_requested(None,meta,'Best modality marker of hsc 1?')['best']['feature']=='IFIT2'
    assert M.answer_if_requested(None,meta,'What are the best marker genes of HSC-1?') is None


def test_redundant_markers_and_changed_files_are_respected(meta,tmp_path):
    path=tmp_path/'rna_redundant_markers.tsv'
    pd.DataFrame({'Gene':['REDUNDANT'],'cluster':['HSC-1'],'rho':[.8]}).to_csv(path,sep='\t',index=False)
    assert M.answer_if_requested(None,meta,'Best modality marker of HSC-1?')['best']['feature']=='REDUNDANT'
    pd.DataFrame({'Gene':['REDUNDANT'],'cluster':['HSC-1'],'rho':[.01]}).to_csv(path,sep='\t',index=False)
    assert M.answer_if_requested(None,meta,'Best modality marker of HSC-1?')['best']['feature']=='IFIT2'


def test_missing_rho_is_not_replaced_with_fold(meta,tmp_path):
    path=tmp_path/'adt_markers.tsv'
    pd.DataFrame({'Gene':['CD73'],'cluster':['HSC-1'],'Fold':[999]}).to_csv(path,sep='\t',index=False)
    result=M.answer_if_requested(None,meta,'Best modality marker of HSC-1?')
    assert not any(r['modality_id']=='adt' for r in result['table']['rows'])
    assert 'adt' in result['answer']


def test_marker_route_preserves_correlation_and_network_questions():
    for q in ['Correlate ADT with cell surface marker gene expression',
              'Show a lipid pathway for the markers of HSC-1',
              'Show the TF activity marker network of HSC-1']:
        assert not M.requested(q)
