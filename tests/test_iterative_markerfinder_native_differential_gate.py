import numpy as np
import pandas as pd
import pytest
from altanalyze3.components.cellHarmony.iterative_marker_finder import annotate_native_differential_eligibility


def write_table(tmp_path, name='contrast.tsv', **changes):
    frame=pd.DataFrame(dict(Feature=['a','b','c','d'],N_case=[3,3,2,3],N_reference=[4,4,4,4],DiffMeans_case_m_reference=[.2,-.3,.8,.1],ebayesPval=[.01,.02,.001,.001],ebayesAdjPval=[.2,.6,.01,.01]))
    for key,value in changes.items():frame[key]=value
    path=tmp_path/name;frame.to_csv(path,sep='\t',index=False);return path


def test_full_universe_native_gate_preserves_order_values_and_boundaries(tmp_path):
    p=write_table(tmp_path);before=p.read_bytes()
    a,proof=annotate_native_differential_eligibility(['d','b','a','c'],[p])
    assert list(a.index)==['d','b','a','c']
    assert a.native_differential_eligible.to_dict()==dict(d=False,b=True,a=True,c=False)
    assert a.loc['b','comparison_1__DiffMeans_case_m_reference']==-.3
    assert proof[0]['rows']==4 and proof[0]['eligible']==2 and p.read_bytes()==before


def test_native_fdr_is_additional_not_substituted_for_raw_p(tmp_path):
    a,_=annotate_native_differential_eligibility(['a','b','c','d'],[write_table(tmp_path)],max_fdr=.5)
    assert a.native_differential_eligible.tolist()==[True,False,False,False]


def test_union_can_rescue_only_independently_qualified_native_rows(tmp_path):
    one=write_table(tmp_path);two=write_table(tmp_path,'second.tsv',N_case=[3,3,3,3],ebayesPval=[1,1,.001,1])
    a,_=annotate_native_differential_eligibility(['a','b','c','d'],[one,two])
    assert a.native_differential_eligible.tolist()==[True,True,True,False]
    assert len(a)==4 and a.loc['c','comparison_1__N_case']==2 and a.loc['c','comparison_2__N_case']==3


@pytest.mark.parametrize('identities',[['a','b','c'],['a','b','c','missing'],['a','a','c','d']])
def test_unresolved_missing_or_duplicate_source_identities_are_rejected(tmp_path,identities):
    with pytest.raises(ValueError):annotate_native_differential_eligibility(identities,[write_table(tmp_path)])


def test_duplicate_native_source_row_is_rejected(tmp_path):
    p=write_table(tmp_path,Feature=['a','b','c','c'])
    with pytest.raises(ValueError):annotate_native_differential_eligibility(['a','b','c','d'],[p])


def test_nan_and_infinite_native_statistics_never_qualify(tmp_path):
    p=write_table(tmp_path,N_case=[3,3,3,3],DiffMeans_case_m_reference=[np.inf,-.3,.8,.2],ebayesPval=[.001,np.nan,.001,.001],ebayesAdjPval=[.01,.01,np.nan,.01])
    a,_=annotate_native_differential_eligibility(['a','b','c','d'],[p],max_fdr=.5)
    assert a.native_differential_eligible.tolist()==[False,False,False,True]
