"""
Does the shared whitelist do what it claims?

The claim under test: annotation 1 draws its full quota per cluster, and every
later annotation draws only the shortfall its own clusters still have. A cluster
of a later annotation that already holds 80 whitelisted cells draws 20 more, not
100.
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))))

from altanalyze3.components.sctriangulate.prefilter import (   # noqa: E402
    downsample_indices, _group_positions)


def reference_downsample(obs, keys, cap, seed=0):
    """
    The same rule written the obvious way: one full scan per cluster, clusters
    in sorted order. downsample_indices groups with a single argsort instead.
    Both must return the same rows.
    """
    rng = np.random.default_rng(seed)
    keep = np.zeros(obs.shape[0], dtype=bool)
    for key in keys:
        labels = obs[key].astype(str).to_numpy()
        for value in sorted(pd.unique(labels)):
            idx = np.flatnonzero(labels == value)
            if len(idx) <= cap:
                keep[idx] = True
                continue
            need = cap - int(keep[idx].sum())
            if need <= 0:
                continue
            available = idx[~keep[idx]]
            keep[rng.choice(available, size=need, replace=False)] = True
    return np.flatnonzero(keep)


def make_obs(n_cells=3000, seed=1):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        'ann1': ['c{}'.format(i) for i in rng.integers(0, 7, n_cells)],
        'ann2': ['d{}'.format(i) for i in rng.integers(0, 13, n_cells)],
        'ann3': ['e{}'.format(i) for i in rng.integers(0, 3, n_cells)],
    }, index=['cell{}'.format(i) for i in range(n_cells)])


def test_group_positions_partitions_every_row():
    labels = np.array(['b', 'a', 'b', 'c', 'a', 'a'])
    names, groups = _group_positions(labels)
    assert names == ['a', 'b', 'c']
    assert [g.tolist() for g in groups] == [[1, 4, 5], [0, 2], [3]]
    assert sorted(np.concatenate(groups).tolist()) == list(range(len(labels)))


def test_matches_the_obvious_implementation():
    obs = make_obs()
    keys = ['ann1', 'ann2', 'ann3']
    for cap in [50, 200, 1000]:
        rows, _info = downsample_indices(obs, keys, cap, seed=0, report=False)
        expected = reference_downsample(obs, keys, cap, seed=0)
        assert np.array_equal(rows, expected), 'cap={}'.format(cap)


def test_every_cluster_of_every_annotation_reaches_its_cap():
    obs = make_obs()
    keys = ['ann1', 'ann2', 'ann3']
    cap = 100
    rows, _info = downsample_indices(obs, keys, cap, seed=0, report=False)
    kept = obs.iloc[rows]
    for key in keys:
        full = obs[key].value_counts()
        got = kept[key].value_counts()
        for cluster, size in full.items():
            assert cluster in got, '{}:{} lost every cell'.format(key, cluster)
            assert got[cluster] >= min(size, cap), \
                '{}:{} kept {} of {}, cap {}'.format(key, cluster, got[cluster], size, cap)


def test_later_annotation_draws_only_the_shortfall():
    """
    Exact arithmetic, no approximation. ann1 is one cluster of 300 cells, so it
    contributes exactly 100 whitelisted cells at cap 100. ann2 is then built
    from that whitelist so each of its clusters has a known head start.
    """
    cap = 100
    n = 300
    obs1 = pd.DataFrame({'ann1': ['only'] * n},
                        index=['cell{}'.format(i) for i in range(n)])
    rows1, _ = downsample_indices(obs1, ['ann1'], cap, seed=0, report=False)
    assert len(rows1) == cap

    outside = np.setdiff1d(np.arange(n), rows1)
    assert len(outside) == 200

    # A: 80 already in the whitelist + 60 outside -> 140 cells, must draw 20
    # B: 20 already in the whitelist + 100 outside -> 120 cells, must draw 80
    # C: the remaining 40 outside cells -> under the cap, taken whole, draws 40
    ann2 = np.empty(n, dtype=object)
    ann2[rows1[:80]] = 'A'
    ann2[outside[:60]] = 'A'
    ann2[rows1[80:]] = 'B'
    ann2[outside[60:160]] = 'B'
    ann2[outside[160:]] = 'C'
    assert not pd.isna(pd.Series(ann2)).any()

    obs = obs1.copy()
    obs['ann2'] = ann2
    rows, info = downsample_indices(obs, ['ann1', 'ann2'], cap, seed=0, report=False)

    assert info['per_annotation']['ann1']['cells_added_here'] == 100
    assert info['per_annotation']['ann2']['cells_added_here'] == 20 + 80 + 40
    assert info['per_annotation']['ann2']['clusters_smaller_than_cap'] == 1     # C
    assert info['per_annotation']['ann2']['clusters_already_at_cap'] == 0
    assert len(rows) == 100 + 140
    assert info['cells_after'] == 240


def test_cluster_already_at_the_cap_draws_nothing():
    cap = 100
    n = 300
    obs1 = pd.DataFrame({'ann1': ['only'] * n},
                        index=['cell{}'.format(i) for i in range(n)])
    rows1, _ = downsample_indices(obs1, ['ann1'], cap, seed=0, report=False)

    outside = np.setdiff1d(np.arange(n), rows1)
    ann2 = np.empty(n, dtype=object)
    ann2[rows1] = 'full'                 # all 100 whitelisted cells
    ann2[outside[:50]] = 'full'          # 150 cells, 100 already in
    ann2[outside[50:]] = 'rest'
    obs = obs1.copy()
    obs['ann2'] = ann2

    _rows, info = downsample_indices(obs, ['ann1', 'ann2'], cap, seed=0, report=False)
    assert info['per_annotation']['ann2']['clusters_already_at_cap'] == 1
    # 'full' adds 0; 'rest' has 150 cells, none whitelisted, so it draws 100
    assert info['per_annotation']['ann2']['cells_added_here'] == 100


def test_seed_is_the_only_source_of_variation():
    obs = make_obs()
    keys = ['ann1', 'ann2', 'ann3']
    a, _ = downsample_indices(obs, keys, 100, seed=0, report=False)
    b, _ = downsample_indices(obs, keys, 100, seed=0, report=False)
    c, _ = downsample_indices(obs, keys, 100, seed=1, report=False)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)


def test_rows_are_sorted_and_unique():
    obs = make_obs()
    rows, _ = downsample_indices(obs, ['ann1', 'ann2', 'ann3'], 100, seed=0, report=False)
    assert np.array_equal(rows, np.sort(rows))
    assert len(set(rows.tolist())) == len(rows)


def test_cap_above_every_cluster_keeps_everything():
    obs = make_obs(n_cells=500)
    rows, info = downsample_indices(obs, ['ann1', 'ann2', 'ann3'], 10000, seed=0, report=False)
    assert len(rows) == 500
    assert info['reduction'] == 1.0


def test_none_cap_is_a_no_op_and_bad_cap_raises():
    obs = make_obs(n_cells=200)
    rows, info = downsample_indices(obs, ['ann1'], None, seed=0, report=False)
    assert np.array_equal(rows, np.arange(200))
    assert info['cells_after'] == 200
    with pytest.raises(ValueError):
        downsample_indices(obs, ['ann1'], 0, seed=0, report=False)


def test_missing_labels_are_refused_not_silently_grouped():
    obs = make_obs(n_cells=100)
    obs.loc[obs.index[0], 'ann1'] = np.nan
    with pytest.raises(ValueError):
        downsample_indices(obs, ['ann1'], 10, seed=0, report=False)
