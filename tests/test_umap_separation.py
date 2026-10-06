import numpy as np
import pytest

from altanalyze3.components.clustering.benchmarking.benchmark_umap_separation import (
    separation, input_neighbor_reference,
)


def test_relative_separation_cannot_be_increased_by_rescaling():
    rng = np.random.default_rng(7)
    x = np.concatenate([rng.normal([0, 0], .1, (40, 2)), rng.normal([2, 0], .1, (40, 2))])
    states = np.repeat(["a", "b"], 40)
    a, b = separation(x, states), separation(x * 100 + 19, states)
    assert a['median_nearest_relative_centroid_distance'] == pytest.approx(
        b['median_nearest_relative_centroid_distance'])
    assert len(a['per_state']) == 2


def test_separation_reports_every_state_including_small_overlapping_states():
    rng = np.random.default_rng(8)
    x = rng.normal(size=(100, 2))
    states = np.array(["large"] * 95 + ["rare"] * 4 + ["singleton"])
    report = separation(x, states)
    assert {s['state']: s['cells'] for s in report['per_state']} == {
        'large': 95, 'rare': 4, 'singleton': 1}
    with pytest.raises(ValueError, match='Finite'):
        separation(x[:-1], states)


def test_input_neighbors_match_original_correlation_distances_without_excluding_cells():
    from scipy.spatial.distance import cdist
    x = np.random.default_rng(11).normal(size=(23, 8))
    states = np.array(["large"] * 20 + ["small"] * 3)
    anchors, actual = input_neighbor_reference(x, states, k=4, per_state=30)
    np.testing.assert_array_equal(anchors, np.arange(len(x)))
    distances = cdist(x, x, metric="correlation")
    np.fill_diagonal(distances, np.inf)
    expected = np.argsort(distances, axis=1)[:, :4]
    assert all(set(a) == set(b) for a, b in zip(actual, expected))
