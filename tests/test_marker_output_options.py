import inspect
import zipfile
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from altanalyze3.components.visualization import marker_heatmap_h5ad as marker
from altanalyze3.components.cellHarmony.flask import pipeline


def data():
    rng = np.random.default_rng(4)
    x = rng.uniform(.01, .1, (80, 4))
    x[:40, :2] += 3
    x[40:, 2:] += 3
    return ad.AnnData(x, obs=pd.DataFrame({'state': ['A'] * 40 + ['B'] * 40},
                                        index=[str(i) for i in range(80)]),
                      var=pd.DataFrame(index=['G1', 'G2', 'G3', 'G4']))


def run(tmp_path, name, **options):
    folder = tmp_path / name
    folder.mkdir()
    return marker.generate_marker_heatmap_from_adata(
        data(), cluster_key='state', out=folder / 'markers.pdf', top_n=2,
        marker_method='markerfinder', validate_scaling=False, centroid_method='mean',
        write_expression_tsv=False, heatmap_dpi=60, **options)


def test_output_controls_preserve_results_and_bypass_expensive_steps(tmp_path, monkeypatch):
    baseline = run(tmp_path, 'baseline')
    assert Path(baseline['pdf']).exists()
    assert Path(baseline['pdf']).with_suffix('.svg').exists()
    def forbidden(*args, **kwargs):
        raise AssertionError('Disabled output step was executed')
    monkeypatch.setattr(marker, '_plot_heatmap', forbidden)
    interactive = run(tmp_path, 'interactive', render_heatmap=False)
    assert interactive['pdf'] is None
    assert Path(interactive['heatmap_cache']).exists()
    matrix, *_ = marker._read_heatmap_cache(interactive['heatmap_cache'])
    assert matrix.shape == (4, 80)
    monkeypatch.setattr(marker, '_build_heatmap', forbidden)
    monkeypatch.setattr(marker, 'downsample_cells_per_group', forbidden)
    minimal = run(tmp_path, 'minimal', render_heatmap=False, write_heatmap_cache=False)
    assert minimal['heatmap_cache'] is None
    assert {p.suffix for p in (tmp_path / 'minimal').iterdir()} == {'.tsv'}
    for key in ['markers_tsv', 'redundant_markers_tsv', 'centroids_tsv']:
        assert Path(baseline[key]).read_bytes() == Path(interactive[key]).read_bytes() == Path(minimal[key]).read_bytes()


def test_no_svg_and_sampling_do_not_change_scoring(tmp_path):
    a = run(tmp_path, 'a', write_svg=False, cells_per_cluster=5)
    b = run(tmp_path, 'b', render_heatmap=False, cells_per_cluster=20)
    assert Path(a['pdf']).exists()
    assert not Path(a['pdf']).with_suffix('.svg').exists()
    for key in ['markers_tsv', 'redundant_markers_tsv', 'centroids_tsv']:
        assert Path(a[key]).read_bytes() == Path(b[key]).read_bytes()


def test_web_output_options_keep_interactive_cache_and_omit_missing_pdf(tmp_path):
    options = pipeline._marker_output_options({})
    # scALABLE skips the static render by default and keeps the interactive cache.
    assert options == dict(render_heatmap=False, write_svg=True, heatmap_dpi=None,
                           cells_per_cluster=100, write_heatmap_cache=True,
                           write_heatmap_tsv=False, write_expression_tsv=False)
    # The MarkerFinder package keeps its own defaults; only scALABLE changed.
    params = inspect.signature(marker.generate_marker_heatmap_from_adata).parameters
    assert params['render_heatmap'].default is True
    assert params['write_svg'].default is True
    assert params['write_expression_tsv'].default is True
    meta = {'qc': {'marker_render_heatmap': False, 'marker_cells_per_cluster': 5}}
    # A same-job rerun must not repackage an old static heatmap that was skipped.
    folder = tmp_path / 'marker_heatmap_adt'
    folder.mkdir()
    (folder / 'cell_state_adt_marker_heatmap.pdf').write_bytes(b'old export')
    outputs = pipeline._emit_modality_marker_heatmap('adt', data(), tmp_path, 'state', meta)
    assert outputs['enabled']
    assert 'heatmap_pdf' not in outputs
    assert Path(outputs['heatmap_cache']).exists()
    assert 'expression_tsv' not in outputs
    with zipfile.ZipFile(outputs['archive']) as archive:
        assert not any(name.endswith('.pdf') for name in archive.namelist())


@pytest.mark.parametrize('dpi', [0, -1, float('nan'), float('inf')])
def test_invalid_dpi_rejected(tmp_path, dpi):
    with pytest.raises(ValueError, match='heatmap_dpi'):
        marker.generate_marker_heatmap_from_adata(data(), cluster_key='state',
                                                 out=tmp_path / 'a.pdf', heatmap_dpi=dpi)


def test_cli_markers_only_and_cache_rerender(tmp_path):
    source = tmp_path / 'input.h5ad'
    data().write_h5ad(source)
    output = tmp_path / 'cli' / 'markers.pdf'
    result = subprocess.run([
        sys.executable, '-m', marker.__name__, '--h5ad', str(source),
        '--cluster-key', 'state', '--out', str(output), '--marker-method', 'markerfinder',
        '--no-scaling-check', '--markers-only'], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert 'Skipping plot matrix construction' in result.stdout
    assert not output.exists()
    assert not list(output.parent.glob('*.h5ad'))
    assert len(list(output.parent.glob('*.tsv'))) == 3
    interactive = run(tmp_path, 'cache', render_heatmap=False)
    marker.render_heatmap_from_cache(interactive['heatmap_cache'], str(output),
                                    load_go_terms=False, dpi=60, write_svg=False)
    assert output.read_bytes().startswith(b'%PDF')
    assert not output.with_suffix('.svg').exists()


def test_web_qc_validates_and_preserves_defaults():
    from importlib import import_module
    model = import_module('altanalyze3.components.cellHarmony.webapp.app').QCSettings
    qc = model()
    assert qc.marker_render_heatmap is False
    assert qc.marker_write_svg is True
    assert qc.marker_heatmap_dpi is None
    assert qc.marker_cells_per_cluster == 100
    for payload in ({'marker_heatmap_dpi': 0}, {'marker_heatmap_dpi': float('nan')},
                    {'marker_cells_per_cluster': -1}):
        with pytest.raises(ValueError):
            model(**payload)
