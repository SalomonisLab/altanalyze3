"""The SURFY + SurfaceGenie union is the default SNAF-B surface database (2026-10-01)."""
import os

import pytest

from altanalyze3.components.snaf.surface import surface_db as sdb


def test_none_and_default_resolve_to_bundled_union():
    assert sdb.resolve_surface_db(None) == sdb.DEFAULT_SURFACE_DB
    assert sdb.resolve_surface_db('default') == sdb.DEFAULT_SURFACE_DB
    assert os.path.isdir(sdb.DEFAULT_SURFACE_DB)


@pytest.mark.parametrize('key', ['alt91', 'ALT91', 'builtin', 'legacy'])
def test_alt91_keywords_select_builtin(key):
    assert sdb.resolve_surface_db(key) is None


def test_missing_path_raises():
    with pytest.raises(FileNotFoundError):
        sdb.resolve_surface_db('/nonexistent/surface_db')


def test_bundled_union_counts():
    info = sdb.load_surface_db(sdb.DEFAULT_SURFACE_DB)
    assert info['stats']['genes_listed'] == 4009
    assert info['stats']['genes_with_reference'] == 3584
    assert info['stats']['reference_sequences'] == 7598
