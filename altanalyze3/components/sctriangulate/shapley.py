'''
Shapley scoring for scTriangulate, as integrated into AltAnalyze3.

This module keeps the upstream scTriangulate 0.13.0 public API
(``wrapper_shapley``, ``shapley_all_or_none_value``, ``shapley_value``,
``rank_based_value``, ``rank_based_all_or_none_value``, ``which_to_take``,
``get_size``, ``size_sort``) and adds a vectorised engine that returns the
same numbers far faster.

WHY THE FAST ENGINE IS EXACT
----------------------------
Upstream computes, for one cell and one player ``i``, a sum over every
non-empty coalition ``S`` of the other ``n-1`` players::

    shapley_i = sum_S  coef(|S|) * surplus(S)
    coef(s)   = s! (n-s-1)! / n!

``surplus(S)`` is read off ``rankdata(total_matrix, method='max', axis=0)``
for the matrix formed by the rows of ``S`` stacked under the player row, after
``cheat_add_bonus``. That enumeration costs ``O(2**(n-1))`` per player per cell.

The surplus decomposes column by column. For score column ``j`` write
``v = data[i, j]`` and classify each other player ``o`` by its value ``d``:

    group A : d <= v            -> count p_j   (player already ranks top vs o)
    group B : v < d <= v + bonus-> count q_j   (bonus lets the player tie o)
    group C : d >  v + bonus    -> count r_j   (player loses to o)

with ``p_j + q_j + r_j = n - 1``. ``rankdata(..., method='max')`` gives the
player rank ``1 + #{o in S : d <= v}``, and ``cheat_add_bonus`` promotes that
rank to ``|S|+1`` exactly when every member of ``S`` lies in ``A | B``.

    mode 'shapley_all_or_none'
        column contribution = (|S|+1) if S subset of A|B else 0, so summing
        over coalitions of size s gives  (s+1) * C(p_j+q_j, s).

    mode 'shapley'
        column contribution = (|S|+1) if S subset of A|B, otherwise
        1 + #{o in S : o in A}. Splitting on how many members come from C and
        applying Vandermonde's identity
        (sum_i C(p,i)C(q,u-i) = C(p+q,u),  sum_i i*C(p,i)C(q,u-i) = p*C(p+q-1,u-1))
        gives, for coalitions of size s,
            (s+1)*C(p+q, s)
          + sum_{l>=1} C(r, l) * [ C(p+q, s-l) + p * C(p+q-1, s-l-1) ].

Both depend on the cell only through the per-column counts ``(p_j, q_j)``. The
engine precomputes a small integer lookup table indexed by those counts, so one
cell costs ``O(n * m)`` instead of ``O(n * 2**(n-1) * m)``.

The tables hold integers scaled by ``n!``; the division by ``n!`` happens once
at the end. Upstream accumulates ``coef(s) * surplus`` in floating point over
thousands of coalitions, so upstream carries rounding noise that the table does
not. Two players whose Shapley values are mathematically equal therefore come
out bitwise equal here, which makes the downstream tie-break deterministic.
``tests/test_shapley_equivalence.py`` checks the fast engine against the
verbatim upstream functions in ``_reference/upstream_shapley.py``.

A second exactness-preserving win: the score matrix of a cell is a function of
which cluster that cell falls in for each annotation, so cells sharing a
cluster tuple share a bitwise identical score matrix. The engine deduplicates
identical matrices and evaluates each distinct one once.
'''

import os
import math
import copy
import logging
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import rankdata


# ---------------------------------------------------------------------------
# cluster size bookkeeping (upstream API, vectorised body)
# ---------------------------------------------------------------------------

def single_size_query(obs, c):
    # c would be {gs:ERP4}
    key = list(c.keys())[0]
    cluster = list(c.values())[0]
    size = obs.loc[obs[key] == cluster, :].shape[0]
    return size


def get_size(obs, query):
    '''
    Return ``(size_dict, size_list)`` for every cluster of every annotation.

    Upstream ran one boolean mask over ``obs`` per cluster. This uses one
    ``value_counts`` per annotation, which is a single pass instead of
    ``n_clusters`` passes, and returns the identical mapping.
    '''
    size_dict = {}   # {gs:{ERP1:54,ERP2:100},...}
    size_list = []   # [  ({gs:ERP1},54),  (),()  ]
    for key in query:
        counts = obs[key].value_counts()
        # upstream iterated obs[key].unique(); restrict to observed categories so
        # unused pandas categories do not appear with size 0
        observed = pd.unique(obs[key])
        size_dict[key] = {}
        for cluster in observed:
            size = int(counts[cluster])
            size_dict[key][cluster] = size
            size_list.append(({key: cluster}, size))
    return size_dict, size_list


def size_sort(size_list):
    result = sorted(size_list, key=lambda x: x[1])
    c, s = zip(*result)   # c means {gs:ERP4}, s means size (int)
    return c, s


# ---------------------------------------------------------------------------
# upstream scalar implementations, kept verbatim in behaviour
# ---------------------------------------------------------------------------

def wrapper_shapley(index, data, mode='shapley_all_or_none', bonus=0.01):
    '''
    Dispatch to one of the four scoring modes for a single player ``index`` of a
    single cell. Kept for API compatibility; ``shapley_batch`` is the fast path.
    '''
    if mode == 'shapley_all_or_none':
        final_result = shapley_all_or_none_value(index, data, bonus)
    elif mode == 'shapley':
        final_result = shapley_value(index, data, bonus)
    elif mode == 'rank':
        final_result = rank_based_value(index, data, bonus)
    elif mode == 'rank_all_or_none':
        final_result = rank_based_all_or_none_value(index, data, bonus)
    else:
        raise ValueError('unknown shapley mode: {}'.format(mode))
    return final_result


def cheat_add_bonus(total_matrix, index_matrix, bonus):
    index_matrix = copy.deepcopy(index_matrix)   # do not mutate the caller's array
    for j in range(index_matrix.shape[1]):       # each score/column
        if index_matrix[-1, j] < index_matrix.shape[0]:   # player not win
            player_score = total_matrix[-1, j]            # a float
            rival_index = np.where(index_matrix[:, j] == index_matrix.shape[0])[0][0]
            rival_score = total_matrix[rival_index, j]    # a float
            if (player_score + bonus) >= rival_score:
                index_matrix[-1, j] = index_matrix.shape[0]   # still think it is best
    return index_matrix


def _all_coalitions(n_others):
    all_combine = []
    for r in range(n_others):
        all_combine.extend(list(combinations(np.arange(n_others), r=r + 1)))
    return all_combine


def shapley_all_or_none_value(index, data, bonus=0.01):
    others = np.delete(data, index, axis=0)
    shapley = 0
    for item in _all_coalitions(others.shape[0]):
        feature_dim = data.shape[1]
        others_matrix = others[np.array(item), :].reshape(-1, feature_dim)
        player_row = data[index, :].reshape(-1, feature_dim)
        total_matrix = np.concatenate([others_matrix, player_row], axis=0)
        index_matrix = rankdata(total_matrix, method='max', axis=0)
        index_matrix = cheat_add_bonus(total_matrix, index_matrix, bonus)
        player_rank = index_matrix[-1, :]
        good_or_not = (player_rank == total_matrix.shape[0])
        player_rank[np.logical_not(good_or_not)] = 0
        surplus = player_rank.sum()
        s = len(item)
        n = data.shape[0]
        normalize_constant = math.factorial(s) * math.factorial(n - s - 1) / math.factorial(n)
        shapley += normalize_constant * surplus
    return shapley


def shapley_value(index, data, bonus=0.01):
    others = np.delete(data, index, axis=0)
    shapley = 0
    for item in _all_coalitions(others.shape[0]):
        feature_dim = data.shape[1]
        others_matrix = others[np.array(item), :].reshape(-1, feature_dim)
        player_row = data[index, :].reshape(-1, feature_dim)
        total_matrix = np.concatenate([others_matrix, player_row], axis=0)
        index_matrix = rankdata(total_matrix, method='max', axis=0)
        index_matrix = cheat_add_bonus(total_matrix, index_matrix, bonus)
        player_rank = index_matrix[-1, :]
        surplus = player_rank.sum()
        s = len(item)
        n = data.shape[0]
        normalize_constant = math.factorial(s) * math.factorial(n - s - 1) / math.factorial(n)
        shapley += normalize_constant * surplus
    return shapley


def rank_based_value(index, data, bonus=0.01):
    index_matrix = rankdata(data, method='max', axis=0)
    index_matrix = cheat_add_bonus(data, index_matrix, bonus)
    value = index_matrix[index, :].sum()
    return value


def rank_based_all_or_none_value(index, data, bonus=0.01):
    index_matrix = rankdata(data, method='max', axis=0)
    index_matrix = cheat_add_bonus(data, index_matrix, bonus)
    player_rank = index_matrix[-1, :]
    good_or_not = (player_rank == index_matrix.shape[0])
    player_rank[np.logical_not(good_or_not)] = 0
    value = player_rank.sum()
    return value


def approximate_shapley_value(data, n_sample=6, n_time=1000):   # for big coalition
    total = np.zeros(shape=data.shape[0])
    counts = np.zeros(shape=data.shape[0])
    indices = np.arange(data.shape[0])
    for t in range(n_time):
        sampled = np.random.choice(a=indices, size=n_sample)
        sub_data = data[sampled, :]
        sub_data_shapley = [shapley_value(i, sub_data) for i in range(sub_data.shape[0])]
        for index, shapley in zip(sampled, sub_data_shapley):
            total[index] += shapley
            counts[index] += 1
    return total / counts


# ---------------------------------------------------------------------------
# fast engine: integer coalition tables
# ---------------------------------------------------------------------------

_TABLE_CACHE = {}
_INT64_MAX = np.iinfo(np.int64).max


def _binom(a, b):
    if b < 0 or a < 0 or b > a:
        return 0
    return math.comb(a, b)


def build_all_or_none_table(n):
    '''
    ``table[k]`` = n! * (column contribution when k of the n-1 other players are
    beaten once the bonus is applied), for mode ``shapley_all_or_none``::

        table[k] = sum_{s=1..n-1} s!(n-s-1)! * (s+1) * C(k, s)

    Integer valued, so summing columns never rounds.
    '''
    table = np.zeros(n, dtype=object)
    for k in range(n):            # k ranges 0..n-1
        total = 0
        for s in range(1, n):
            c = _binom(k, s)
            if c:
                total += math.factorial(s) * math.factorial(n - s - 1) * (s + 1) * c
        table[k] = total
    return table


def build_shapley_table(n):
    '''
    ``table[p, q]`` = n! * (column contribution) for mode ``shapley``, where p
    counts other players the cell already beats and q counts those the bonus
    rescues. ``r = n-1-p-q`` players are lost outright. Entries with
    ``p + q > n-1`` are unreachable and stay zero::

        table[p,q] = sum_{s=1..n-1} s!(n-s-1)! * [ (s+1)*C(p+q, s)
                     + sum_{l=1..min(r,s)} C(r,l)*( C(p+q,s-l) + p*C(p+q-1,s-l-1) ) ]
    '''
    table = np.zeros((n, n), dtype=object)
    for p in range(n):
        for q in range(n - p):
            r = n - 1 - p - q
            total = 0
            for s in range(1, n):
                coef = math.factorial(s) * math.factorial(n - s - 1)
                inner = (s + 1) * _binom(p + q, s)
                for l in range(1, min(r, s) + 1):
                    inner += _binom(r, l) * (_binom(p + q, s - l)
                                             + p * _binom(p + q - 1, s - l - 1))
                if inner:
                    total += coef * inner
            table[p, q] = total
    return table


def get_coalition_table(n, mode):
    '''Cached integer table for ``n`` players and the given mode.'''
    key = (n, mode)
    if key in _TABLE_CACHE:
        return _TABLE_CACHE[key]
    if mode == 'shapley_all_or_none':
        table = build_all_or_none_table(n)
    elif mode == 'shapley':
        table = build_shapley_table(n)
    else:
        raise ValueError('no coalition table for mode {}'.format(mode))
    biggest = int(max(table.ravel().tolist())) if table.size else 0
    if biggest > _INT64_MAX // 4096:
        # keep python ints; the caller falls back to object dtype arithmetic
        table = table.astype(object)
    else:
        table = table.astype(np.int64)
    _TABLE_CACHE[key] = table
    return table


def shapley_batch(data, mode='shapley_all_or_none', bonus=0.01, dedup=True):
    '''
    Score every player for every cell.

    :param data: ndarray ``(n_players, n_cells, n_metrics)``, the layout upstream
                 builds in ``ScTriangulate.compute_shapley``.
    :param mode: one of ``shapley_all_or_none``, ``shapley``, ``rank``,
                 ``rank_all_or_none``.
    :param bonus: float offset, as upstream.
    :param dedup: collapse bitwise identical per-cell score matrices before
                  scoring. Cells in the same cluster tuple share a matrix, so
                  this is lossless.
    :return: ndarray ``(n_cells, n_players)`` of scores.
    '''
    data = np.ascontiguousarray(data, dtype=np.float64)
    n_players, n_cells, n_metrics = data.shape
    if n_cells == 0:
        return np.zeros((0, n_players))

    # (n_cells, n_players, n_metrics)
    per_cell = np.transpose(data, (1, 0, 2))

    if dedup and n_cells > 1:
        flat = per_cell.reshape(n_cells, n_players * n_metrics)
        uniq, inverse = np.unique(flat, axis=0, return_inverse=True)
        work = uniq.reshape(-1, n_players, n_metrics)
    else:
        inverse = None
        work = per_cell

    if mode in ('shapley_all_or_none', 'shapley'):
        scores = _shapley_core(work, n_players, mode, bonus)
    elif mode == 'rank':
        scores = _rank_core(work, n_players, bonus, all_or_none=False)
    elif mode == 'rank_all_or_none':
        scores = _rank_core(work, n_players, bonus, all_or_none=True)
    else:
        raise ValueError('unknown shapley mode: {}'.format(mode))

    if inverse is not None:
        scores = scores[inverse]
    return scores


def _shapley_core(work, n_players, mode, bonus):
    '''work: (R, n_players, n_metrics) -> (R, n_players) float scores.'''
    table = get_coalition_table(n_players, mode)
    n_rows = work.shape[0]
    use_object = table.dtype == object
    acc_dtype = object if use_object else np.int64
    totals = np.zeros((n_rows, n_players), dtype=acc_dtype)

    for idx in range(n_players):
        others_idx = [j for j in range(n_players) if j != idx]
        v = work[:, idx, :]                                  # (R, m)
        others = work[:, others_idx, :]                      # (R, n-1, m)
        p = (others <= v[:, None, :]).sum(axis=1)            # (R, m)
        k = (others <= (v + bonus)[:, None, :]).sum(axis=1)  # (R, m)
        if mode == 'shapley_all_or_none':
            picked = table[k]
        else:
            picked = table[p, k - p]
        totals[:, idx] = picked.sum(axis=1)

    denom = float(math.factorial(n_players))
    if use_object:
        return np.array([[float(x) / denom for x in row] for row in totals])
    return totals.astype(np.float64) / denom


def _rank_core(work, n_players, bonus, all_or_none):
    '''
    Reproduce upstream ``rank_based_value`` / ``rank_based_all_or_none_value``.

    Upstream note, preserved deliberately: ``cheat_add_bonus`` only ever edits
    the LAST row of the rank matrix, and both rank modes call it on the full
    ``data`` matrix rather than on a player-last matrix. So in mode ``rank`` the
    bonus reaches only the last player, and in mode ``rank_all_or_none`` the
    returned value ignores ``index`` altogether and is identical for every
    player. Both quirks are upstream behaviour and are reproduced exactly.
    '''
    # rank[r, y, j] = #{x : work[r, x, j] <= work[r, y, j]}   (rankdata method='max')
    cmp = work[:, :, None, :] <= work[:, None, :, :]          # (R, x, y, m)
    ranks = cmp.sum(axis=1).astype(np.int64)                  # (R, y, m)

    last = n_players - 1
    col_max = work.max(axis=1)                                # (R, m)
    promote = (ranks[:, last, :] < n_players) & \
              ((work[:, last, :] + bonus) >= col_max)
    ranks[:, last, :] = np.where(promote, n_players, ranks[:, last, :])

    if not all_or_none:
        return ranks.sum(axis=2).astype(np.float64)

    won = ranks[:, last, :] == n_players
    value = (won.sum(axis=1) * n_players).astype(np.float64)  # (R,)
    return np.repeat(value[:, None], n_players, axis=1)


# ---------------------------------------------------------------------------
# winner selection
# ---------------------------------------------------------------------------

def which_to_take(result, query, reference, cluster_row, size_dict):
    '''
    Pick the annotation a single cell should adopt. Kept for API compatibility;
    ``which_to_take_batch`` is the fast path and gives the same answer.
    '''
    rank = rankdata(result, method='max')
    number_of_winner = len(np.where(rank == len(query))[0])
    if number_of_winner == 1:
        to_take = query[np.where(rank == len(query))[0][0]]
    else:
        winners = np.where(rank == len(query))[0]
        size = [size_dict[query[index]][cluster_row[index]] for index in winners]
        to_take = query[winners[size.index(min(size))]]
    return to_take


def which_to_take_batch(scores, obs, query, size_dict):
    '''
    Vectorised winner selection over all cells.

    :param scores: ndarray ``(n_cells, n_query)``.
    :param obs: DataFrame carrying one column per annotation in ``query``.
    :param query: list of annotation column names.
    :param size_dict: ``{annotation: {cluster: size}}``.
    :return: ndarray of dtype object, the chosen annotation name per cell.

    Matches upstream rule for rule: the winners are the players holding the row
    maximum, and among them the one whose own cluster is smallest wins, ties
    going to the earliest annotation in ``query`` (``list.index`` picks the first
    minimum upstream; ``argmin`` picks the first minimum here).
    '''
    scores = np.asarray(scores)
    n_cells, n_query = scores.shape
    is_winner = scores == scores.max(axis=1, keepdims=True)

    sizes = np.empty((n_cells, n_query), dtype=np.float64)
    for j, key in enumerate(query):
        sizes[:, j] = obs[key].map(size_dict[key]).to_numpy(dtype=np.float64)

    masked = np.where(is_winner, sizes, np.inf)
    choice = masked.argmin(axis=1)
    # a single winner short-circuits upstream; argmin over the masked sizes
    # returns that same column because every other column is inf
    return np.asarray(query, dtype=object)[choice]


def run_shapley_batch(obs, query, reference, size_dict, data, mode, bonus):
    '''
    Fast replacement for upstream ``main_class.run_shapley``.

    :return: ``(final, intermediate)`` where ``final`` is the chosen annotation
             per cell and ``intermediate`` is the ``(n_cells, n_query)`` score
             matrix as a list of per-cell lists, matching upstream's shape.
    '''
    scores = shapley_batch(data, mode=mode, bonus=bonus)
    final = which_to_take_batch(scores, obs, query, size_dict)
    return list(final), scores.tolist()
