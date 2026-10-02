"""SNAF-T skips and reports HLA alleles MHCflurry cannot score (Nathan, 2026-10-01)."""
import os
from types import SimpleNamespace

import pandas as pd
import pytest

from altanalyze3.components.snaf import binding
from altanalyze3.components.snaf import snaf as snaf_mod

S = ['s1.bed', 's2.bed']
HLAS = [['HLA-A*02:01', 'HLA-A*02:01', 'HLA-C*03:540'], ['HLA-A*01:01', 'HLA-B*35:578', 'HLA-B*40:01']]


def _fake(monkeypatch, bad, method='MHCflurry'):
    monkeypatch.setattr(snaf_mod, 'binding_method', method, raising=False)
    monkeypatch.setattr(binding, 'mhcflurry_unsupported_alleles',
                        lambda alleles: {a: 'no pseudosequence' for a in alleles if a in bad})
    return SimpleNamespace(junction_count_matrix=pd.DataFrame(columns=S))


def test_skips_and_reports(tmp_path, monkeypatch):
    jc = _fake(monkeypatch, {'HLA-C*03540', 'HLA-B*35578'})
    kept = snaf_mod.JunctionCountMatrixQuery._drop_unsupported_hlas(jc, HLAS, str(tmp_path))
    assert kept == [['HLA-A*02:01', 'HLA-A*02:01'], ['HLA-A*01:01', 'HLA-B*40:01']]
    rep = pd.read_csv(tmp_path / 'hla_alleles_skipped.txt', sep='\t')
    assert rep[['sample', 'allele']].values.tolist() == [['s1.bed', 'HLA-C*03:540'], ['s2.bed', 'HLA-B*35:578']]


def test_none_skipped_is_a_noop_with_empty_report(tmp_path, monkeypatch):
    jc = _fake(monkeypatch, set())
    assert snaf_mod.JunctionCountMatrixQuery._drop_unsupported_hlas(jc, HLAS, str(tmp_path)) == HLAS
    rep = pd.read_csv(tmp_path / 'hla_alleles_skipped.txt', sep='\t')
    assert list(rep.columns) == ['sample', 'allele', 'reason'] and rep.empty


def test_sample_left_with_no_allele_is_refused(tmp_path, monkeypatch):
    jc = _fake(monkeypatch, {'HLA-A*0201', 'HLA-C*03540'})
    with pytest.raises(ValueError, match='s1.bed'):
        snaf_mod.JunctionCountMatrixQuery._drop_unsupported_hlas(jc, HLAS, str(tmp_path))


def test_netmhcpan_untouched(tmp_path, monkeypatch):
    jc = _fake(monkeypatch, {'HLA-C*03540'}, method='netMHCpan')
    assert snaf_mod.JunctionCountMatrixQuery._drop_unsupported_hlas(jc, HLAS, str(tmp_path)) is HLAS
    assert not os.path.exists(tmp_path / 'hla_alleles_skipped.txt')


def test_real_mhcflurry_rejects_the_h358_alleles():
    pytest.importorskip('mhcflurry')
    try:
        bad = binding.mhcflurry_unsupported_alleles(['HLA-C*03540', 'HLA-B*35578', 'HLA-A*0301', 'HLA-C*0304'])
    except Exception as exc:
        pytest.skip('MHCflurry models unavailable: {}'.format(exc))
    assert sorted(bad) == ['HLA-B*35578', 'HLA-C*03540']
