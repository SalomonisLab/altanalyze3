"""The vectorized KDE map must reproduce the published algorithm, not merely resemble it."""
import numpy as np
import pytest
from altanalyze3.components.rna2flow.normalize import (
    scale_feature, kde_quantile_map, kde_quantile_map_reference, rank_plotting_positions,
)


def _panels(seed=0, n_ref=4000, n_q=3000):
    rng = np.random.default_rng(seed)
    ref = np.concatenate([rng.normal(0.3, 0.08, n_ref // 2), rng.normal(0.75, 0.12, n_ref // 2)])
    qry = np.concatenate([rng.lognormal(0.0, 0.5, n_q // 2), rng.normal(4.0, 1.0, n_q // 2)])
    return scale_feature(ref, 1, 99), scale_feature(qry, 1, 99)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_vectorized_matches_published(seed):
    ref, qry = _panels(seed)
    slow = kde_quantile_map_reference(qry, ref)
    fast = kde_quantile_map(qry, ref)
    assert np.allclose(slow, fast, atol=2e-3), np.abs(slow - fast).max()


def test_rank_positions_match_published_formula():
    rng = np.random.default_rng(7)
    x = rng.normal(size=500)
    order = np.argsort(x)
    import pandas as pd
    published = pd.Series(order, index=range(len(x))).sort_values().index.values
    published = (published + 1) / (len(x) + 2)
    assert np.allclose(rank_plotting_positions(x), published)


def test_map_is_monotone_in_the_query():
    """A quantile map must preserve order. If it does not, every downstream claim is void.

    The sort here must be STABLE, to match rank_plotting_positions. An unstable sort orders
    tied values differently from the implementation and reports false violations; that is a
    faulty test, not a faulty map.
    """
    ref, qry = _panels(3)
    out = kde_quantile_map(qry, ref)
    o = np.argsort(qry, kind="stable")
    assert np.all(np.diff(out[o]) >= -1e-9)


def test_published_ties_map_to_different_values_and_average_ties_do_not():
    """Percentile clipping manufactures ties; the published method splits them."""
    ref, _ = _panels(5)
    qry = np.concatenate([np.full(50, 0.4), np.linspace(0, 1, 50)])
    first = kde_quantile_map(qry, ref, ties="first")
    avg = kde_quantile_map(qry, ref, ties="average")
    assert np.unique(first[:50]).size > 1, "published behaviour: tied inputs split"
    assert np.unique(np.round(avg[:50], 12)).size == 1, "average: tied inputs stay tied"


def test_mapped_distribution_matches_the_reference():
    """The point of the map: query quantiles should land on reference quantiles."""
    ref, qry = _panels(4)
    out = kde_quantile_map(qry, ref)
    q = np.linspace(0.05, 0.95, 19)
    assert np.max(np.abs(np.quantile(out, q) - np.quantile(ref, q))) < 0.05


def test_constant_feature_does_not_emit_nan():
    assert not np.isnan(scale_feature(np.full(100, 3.0))).any()
