"""Dense iterative plots must inherit the native raster resolution, not lose rows."""
import inspect
import anndata as ad
import numpy as np
import pandas as pd
from altanalyze3.components.cellHarmony import iterative_marker_finder as imf
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh


def test_dense_iterative_native_default_has_pixels_for_every_event(monkeypatch, tmp_path):
    # The failing real figure had 4,307 rows but only 3,556 raster scanlines at
    # the old iterative-only 600-DPI override. Exercise the actual renderer geometry.
    n = 4307
    monkeypatch.delenv('MARKER_HEATMAP_DPI', raising=False)
    a = ad.AnnData(np.tile(np.linspace(0, 1, n), (4, 1)).astype(np.float32),
                   obs=pd.DataFrame({'condition': ['young', 'aged', 'other', 'SRSF2']}),
                   var=pd.DataFrame(index=[f'feature{i}' for i in range(n)]))
    assignments = pd.DataFrame({'pattern': ['SRSF2'] * n}, index=a.var_names)
    resolutions = []

    def capture(fig, path, dpi=None, **kwargs):
        fig.canvas.draw()
        heatmap = next(im for ax in fig.axes for im in ax.images
                       if im.get_array().shape == (n, 4))
        height_inches = heatmap.axes.get_window_extent().height / fig.dpi
        resolutions.append(height_inches * (dpi or mh._heatmap_dpi()))

    monkeypatch.setattr(mh, '_save_figure', capture)
    out = mh.render_consolidated_marker_heatmap(a, assignments, cluster_key='condition',
            out=tmp_path/'dense.pdf', top_n=None, label_events=[],
            covariate_columns=['condition'])
    assert out['plotted'] == n
    assert resolutions[0] >= n
    cache = ad.read_h5ad(out['heatmap_cache'])
    assert cache.n_vars == n and cache.n_obs == 4


def test_iterative_entry_points_do_not_override_native_raster_default():
    assert inspect.signature(imf.review_legacy_markerfinder_fit).parameters['dpi'].default is None
    assert inspect.signature(imf.run_legacy_iterative_marker_finder).parameters['heatmap_dpi'].default is None
