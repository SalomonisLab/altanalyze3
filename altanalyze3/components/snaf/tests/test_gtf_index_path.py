"""SNAF_INDEX_DIR keys the parsed-GTF index by full path, not basename (2026-10-01).

Long-read catalogues are often all named combined.gff.gz; a basename key made them share one index.
"""
import os

from altanalyze3.components.snaf.surface.main import _gtf_index_path


def test_same_basename_different_dirs_get_distinct_indexes(tmp_path, monkeypatch):
    monkeypatch.setenv('SNAF_INDEX_DIR', str(tmp_path / 'index'))
    a = _gtf_index_path('/data/x/KINNEX-5/combined.gff.gz')
    b = _gtf_index_path('/data/y/ENCODE/combined.gff.gz')
    assert a != b
    assert os.path.dirname(a) == os.path.dirname(b) == str(tmp_path / 'index')
    assert os.path.basename(a).startswith('combined.gff.gz.')
    assert a.endswith('.snaf_gtf_index.pkl')


def test_same_path_is_stable(tmp_path, monkeypatch):
    monkeypatch.setenv('SNAF_INDEX_DIR', str(tmp_path / 'index'))
    assert _gtf_index_path('/data/x/combined.gff.gz') == _gtf_index_path('/data/x/combined.gff.gz')


def test_default_sidecar_unchanged(monkeypatch):
    monkeypatch.delenv('SNAF_INDEX_DIR', raising=False)
    assert _gtf_index_path('/data/x/combined.gff.gz') == '/data/x/combined.gff.gz.snaf_gtf_index.pkl'
