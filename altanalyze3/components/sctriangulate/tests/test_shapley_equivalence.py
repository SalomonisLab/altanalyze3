"""Prove the batched Shapley engine returns upstream's numbers.

The reference is `_reference/upstream_shapley.py`, a verbatim copy of
scTriangulate 0.13.0's `sctriangulate/shapley.py` (commit 8b9598cf, 2026-03-30).
Nothing in this test imports the optimized implementations of those functions.

Run:
  cd /Users/saljh8/Documents/GitHub/altanalyze3
  python3.11 -m pytest altanalyze3/components/sctriangulate/tests/test_shapley_equivalence.py -q
"""

import os
import sys
import importlib.util

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
COMPONENT = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(COMPONENT, '..', '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from altanalyze3.components.sctriangulate import shapley as fast   # noqa: E402


def _load_reference():
    path = os.path.join(COMPONENT, '_reference', 'upstream_shapley.py')
    spec = importlib.util.spec_from_file_location('upstream_shapley', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


REF = _load_reference()
MODES = ['shapley_all_or_none', 'shapley', 'rank', 'rank_all_or_none']
BONUS = 0.01


def _reference_scores(data3d, mode, bonus=BONUS):
    """Upstream loop: one call per (cell, player)."""
    n_players, n_cells, _ = data3d.shape
    out = np.empty((n_cells, n_players))
    for i in range(n_cells):
        layer = data3d[:, i, :]
        for j in range(n_players):
            out[i, j] = REF.wrapper_shapley(j, layer, mode, bonus)
    return out


def _make_cases(rng, n_players, n_metrics, n_cells):
    """Score matrices covering plain values, exact ties, and within-bonus gaps."""
    plain = rng.random((n_players, n_cells, n_metrics))

    # exact ties: duplicate player 0's row into player 1
    tied = plain.copy()
    if n_players > 1:
        tied[1] = tied[0]

    # near ties strictly inside the bonus window, which is where cheat_add_bonus
    # decides the outcome
    near = plain.copy()
    if n_players > 1:
        near[1] = near[0] + BONUS * rng.uniform(0.1, 0.9, size=near[0].shape)

    # coarse grid values, so ties happen all over the matrix
    coarse = rng.integers(0, 3, size=(n_players, n_cells, n_metrics)).astype(float) / 2.0

    # zeros, the degenerate all-equal case
    zeros = np.zeros((n_players, n_cells, n_metrics))

    return {'plain': plain, 'tied': tied, 'near_bonus': near,
            'coarse': coarse, 'zeros': zeros}


@pytest.mark.parametrize('mode', MODES)
@pytest.mark.parametrize('n_players', [2, 3, 4, 5, 6])
def test_batch_matches_upstream(mode, n_players):
    rng = np.random.default_rng(0xC0FFEE + n_players)
    n_metrics, n_cells = 4, 25
    for name, data in _make_cases(rng, n_players, n_metrics, n_cells).items():
        expected = _reference_scores(data, mode)
        got = fast.shapley_batch(data, mode=mode, bonus=BONUS)
        assert got.shape == expected.shape, (mode, n_players, name)
        assert np.allclose(got, expected, rtol=0, atol=1e-9), \
            '{} / {} players / {}: max abs diff {}'.format(
                mode, n_players, name, np.abs(got - expected).max())


@pytest.mark.parametrize('mode', MODES)
def test_dedup_changes_nothing(mode):
    rng = np.random.default_rng(7)
    data = rng.integers(0, 2, size=(4, 200, 3)).astype(float)   # many repeated rows
    with_dedup = fast.shapley_batch(data, mode=mode, bonus=BONUS, dedup=True)
    without = fast.shapley_batch(data, mode=mode, bonus=BONUS, dedup=False)
    assert np.array_equal(with_dedup, without)


@pytest.mark.parametrize('n_players', [3, 5, 7])
def test_wide_metric_matrix(n_players):
    """More score columns than players, the shape real runs use."""
    rng = np.random.default_rng(11 * n_players)
    data = rng.random((n_players, 15, 9))
    for mode in MODES:
        expected = _reference_scores(data, mode)
        got = fast.shapley_batch(data, mode=mode, bonus=BONUS)
        assert np.allclose(got, expected, rtol=0, atol=1e-9), mode


def test_tables_are_exact_integers():
    """The lookup tables must stay integral, or ties stop being decidable."""
    for n in range(2, 13):
        aon = fast.build_all_or_none_table(n)
        sh = fast.build_shapley_table(n)
        assert all(float(x).is_integer() for x in aon.tolist())
        assert all(float(x).is_integer() for row in sh.tolist() for x in row)


def test_which_to_take_batch_matches_upstream():
    import pandas as pd
    rng = np.random.default_rng(3)
    query = ['a', 'b', 'c']
    n_cells = 400
    obs = pd.DataFrame({k: rng.integers(0, 4, n_cells).astype(str) for k in query})
    size_dict = {k: {str(c): int(v) for c, v in obs[k].value_counts().items()} for k in query}
    # integer scores force frequent exact ties, the branch that matters
    scores = rng.integers(0, 3, size=(n_cells, len(query))).astype(float)

    expected = [REF.which_to_take(scores[i], query, query[0],
                                  obs.iloc[i].loc[query].values, size_dict)
                for i in range(n_cells)]
    got = fast.which_to_take_batch(scores, obs, query, size_dict)
    assert list(got) == expected
