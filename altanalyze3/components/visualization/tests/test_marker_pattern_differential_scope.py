import pandas as pd
import pytest
from altanalyze3.components.cellHarmony.iterative_marker_finder import apply_pattern_differential_gate


def test_gate_changes_only_requested_patterns_before_display_cap():
    ids=pd.Index(['a','b','c','d','e','f','g','h'])
    audit=pd.DataFrame({'figure_eligible':[True]*7+[False], 'rho':[.9,.88,.86,.84,.82,.80,.78,.76]},index=ids)
    original=audit.copy(deep=True)
    names=pd.Series(['p1','p2','p3','p4','p5','p6','p7','p3'],index=ids)
    gate=pd.Series([False,False,False,True,False,True,False,True],index=ids)
    result=apply_pattern_differential_gate(audit,gate,names,['p3','p4','p5','p6'])
    assert result.index.equals(ids)
    assert result.figure_eligible.tolist()==[True,True,False,True,False,True,True,False]
    assert result.native_differential_applies.tolist()==[False,False,True,True,True,True,False,True]
    pd.testing.assert_frame_equal(audit,original)
    pd.testing.assert_series_equal(result.rho,original.rho)
    # No scope keeps the existing global-gate behavior.
    global_result=apply_pattern_differential_gate(audit,gate,names)
    assert global_result.figure_eligible.tolist()==[False,False,False,True,False,True,False,False]


def test_gate_rejects_unresolved_source_and_pattern_matches():
    audit=pd.DataFrame({'figure_eligible':[True]},index=['a'])
    for eligibility,names in [(pd.Series(dtype=bool),pd.Series({'a':'p1'})),
                              (pd.Series({'a':True}),pd.Series(dtype=str))]:
        with pytest.raises(ValueError):apply_pattern_differential_gate(audit,eligibility,names,['p1'])
