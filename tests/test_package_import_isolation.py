"""Optional pydata/sparse must not resolve to AltAnalyze3's sparse directory."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize('first_import', ['gff', 'direct_gff', 'markerfinder'])
def test_package_imports_do_not_expose_components_to_optional_dependencies(first_import):
    root = Path(__file__).resolve().parents[1]
    # Model an installation without the OPTIONAL pydata/sparse distribution,
    # while leaving the real in-repo namespace visible if a module pollutes
    # sys.path. This reproduces dask's SparseArray error even on a developer
    # machine that happens to have pydata/sparse installed.
    script = r'''
from importlib.abc import MetaPathFinder
from importlib.machinery import PathFinder
from pathlib import Path
import sys
components = Path(sys.argv[1]) / 'altanalyze3' / 'components'
class OptionalSparseAbsent(MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sparse':
            local_paths = [p for p in sys.path if Path(p).resolve() == components]
            spec = PathFinder.find_spec(fullname, local_paths)
            if spec is not None:
                return spec
            raise ModuleNotFoundError('Optional pydata/sparse is absent', name='sparse')
sys.meta_path.insert(0, OptionalSparseAbsent())
before = list(sys.path)
if sys.argv[2] == 'direct_gff':
    import runpy
    # Probe the same bootstrap used by python /path/to/gff_process.py, from
    # outside the checkout and without its root already in sys.path.
    sys.path[:] = [p for p in sys.path if Path(p).resolve() != Path(sys.argv[1])]
    before = list(sys.path)
    module = runpy.run_path(str(components / 'long_read' / 'gff_process.py'), run_name='gff_import_probe')
    assert module['_coord_store'] is not None
elif sys.argv[2] == 'gff':
    from altanalyze3.components.long_read import gff_process
    assert gff_process._coord_store is not None
    from altanalyze3.components.bam import coord_store
    assert gff_process._coord_store is coord_store
else:
    from altanalyze3.components.visualization import marker_heatmap_h5ad
    from altanalyze3.components.cellHarmony import markerFinder
    assert marker_heatmap_h5ad._load_marker_finder() is markerFinder
assert sys.path == before, (before, sys.path)
import anndata
import dask.array
assert 'sparse' not in sys.modules
print('PASS: canonical imports and optional-dependency isolation')
'''
    env = {**os.environ, 'PYTHONPATH': str(root)}
    result = subprocess.run([sys.executable, '-c', script, str(root), first_import],
                            cwd=root.parent if first_import == 'direct_gff' else root,
                            env=env, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr

