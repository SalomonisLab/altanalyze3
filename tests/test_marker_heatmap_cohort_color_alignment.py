"""Regress the rendered cohort/pattern strip normalization, including its legend."""
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_rgba
from altanalyze3.components.visualization import marker_heatmap_h5ad as mh


@pytest.mark.parametrize('pattern_count', [0, 2, 8])
def test_rendered_integer_strip_ids_keep_exact_category_colors(monkeypatch, tmp_path, pattern_count):
    cohorts = ['young', 'aged', 'other', 'SRSF2']
    patterns = [f'pattern{i}' for i in range(pattern_count)] or cohorts
    order = cohorts + [p for p in patterns if p not in cohorts]
    matrix = pd.DataFrame(np.arange(len(patterns)*4).reshape(-1,4),
                          index=[f'event{i}' for i in range(len(patterns))], columns=['A','B','C','D'])
    cov = pd.DataFrame({'condition': pd.Categorical(cohorts, categories=cohorts)}, index=matrix.columns)
    saved = []
    monkeypatch.setattr(mh, '_save_figure', lambda fig, *args, **kwargs: saved.append(fig))
    mh._plot_heatmap(matrix, tmp_path/'view.pdf', [(c,1) for c in cohorts], order,
                     cohorts, patterns, covariate_df=cov, write_svg=False, label_events=[])
    fig = saved[0]
    # The only single-row image is the actual upper cohort strip. Row strip is n-by-1.
    strips = [im for ax in fig.axes for im in ax.images if im.get_array().shape == (1,4)]
    assert len(strips) == 1
    top = strips[0]
    _, _, legends = mh._build_covariate_track(cov)
    rgba = top.cmap(top.norm(top.get_array()))[0]
    for i,c in enumerate(cohorts):
        np.testing.assert_array_equal(rgba[i], to_rgba(legends['condition']['colors'][c]))
    row = [im for ax in fig.axes for im in ax.images if im.get_array().shape == (len(patterns),1)][0]
    palette, _ = mh._cluster_color_map(order)
    row_rgba = row.cmap(row.norm(row.get_array()))[:,0]
    for i,p in enumerate(patterns):
        np.testing.assert_array_equal(row_rgba[i], palette[p])
