'''
Reduce the gene axis and the cell axis before the heavy steps.

This module holds two independent reductions. Both are off by default. Both
change the input to every metric, so both are opt-in and both write what they
did into ``run_parameters.json``.

* ``select_marker_genes`` cuts the gene axis with AltAnalyze3's MarkerFinder.
* ``downsample_indices`` cuts the cell axis with a shared whitelist that every
  annotation contributes to.

GENE AXIS
---------
scTriangulate's cost and memory both scale with the gene count, and almost all of
those genes carry no information about which clustering is right. At
984,119 cells by 36,249 genes the dense matrix alone is 142.7 GB (132.9 GiB).
The same matrix as CSR at 5% nonzero is 13.3 GiB, and after reducing to the genes
that actually mark a cell state it is about 3 GiB.

THE PROCEDURE
-------------
1. Take at most ``max_cells_per_cluster`` cells from each cluster of each
   annotation. Marker ranking needs enough cells per state to be stable, not
   every cell. 1,000 per cluster caps the marker step at
   ``n_clusters * 1,000`` cells however large the dataset is.
2. Score genes against cluster membership with
   ``altanalyze3.components.cellHarmony.markerFinder.marker_finder``, which is
   already sparse end to end: it takes a CSR matrix and works from column sums,
   squared column sums and one sparse matmul. Nothing is densified.
3. Keep the top ``n_per_cluster`` genes per cluster, take the union over every
   cluster of every annotation, and subset the object to that union.

WHAT THIS CHANGES
-----------------
This is a real change to the input of every metric, not a reformulation. Genes
that were dropped can no longer appear in a marker list, contribute to a tf-idf
score, or carry a SCCAF coefficient. It is off by default for that reason, and
``benchmarks/prefilter_validate.py`` measures the effect on the final labels with
ARI before anyone relies on it.

The reduction is bounded below by the union size. MarkerFinder scores every gene
against every cluster independently, so one gene can be a top marker for several
clusters and the union is usually far smaller than
``n_clusters * n_per_cluster``.

CELL AXIS
---------
``downsample_indices`` keeps at most ``max_per_cluster`` cells for every cluster
of every annotation, and it counts cells that another annotation already
contributed.

The procedure walks the annotations in order and fills one shared whitelist:

1. Annotation 1, cluster by cluster. Take a random ``max_per_cluster`` cells, or
   the whole cluster when it is smaller. Put them in the whitelist.
2. Annotation 2, cluster by cluster. Count how many of that cluster's cells the
   whitelist already holds. Two annotations that agree share most of their
   cells, so this count is often already at the cap and the cluster adds
   nothing. A cluster holding 980 whitelisted cells draws only 20 more, and it
   draws them from the cells the whitelist does not hold.
3. Repeat for every remaining annotation.

The order matters. Annotation 1 spends its full budget; later annotations spend
only what the earlier ones left unspent. The whitelist is therefore far smaller
than ``sum over annotations of n_clusters * max_per_cluster``, and no cluster of
any annotation ever falls below its cap unless the cluster itself is smaller.

Every cluster of every annotation keeps at least one cell, so no annotation
loses a cluster and the Shapley step still sees the same competitors.

WHAT THIS CHANGES
-----------------
Cell counts are the denominator of every stability metric, so a downsampled run
is not the same experiment as a full run. Rare cells inside a large cluster can
be dropped. Use it to make a run possible, not to make a run faster, and report
that you used it.
'''

import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse

from .logger import logger_sctriangulate


def _load_marker_finder():
    '''Import AltAnalyze3's MarkerFinder, whatever the import root looks like.'''
    try:
        from ..cellHarmony import markerFinder as mf
        return mf
    except Exception:
        pass
    components_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if components_dir not in sys.path:
        sys.path.insert(0, components_dir)
    from cellHarmony import markerFinder as mf   # noqa: E402
    return mf


def subsample_indices(labels, max_per_group, seed=0):
    '''
    Row positions holding at most ``max_per_group`` cells from each label.

    Deterministic for a given seed and label vector. Groups smaller than the cap
    are taken whole. The returned positions are sorted, so the subsample keeps
    the original cell order and a CSR row slice stays cheap.
    '''
    labels = np.asarray(labels)
    rng = np.random.default_rng(seed)
    picked = []
    for value in pd.unique(labels):
        idx = np.flatnonzero(labels == value)
        if max_per_group is not None and len(idx) > max_per_group:
            idx = rng.choice(idx, size=max_per_group, replace=False)
        picked.append(idx)
    return np.sort(np.concatenate(picked)) if picked else np.array([], dtype=int)


def _group_positions(labels):
    '''
    Row positions of each distinct label, in lexicographic label order.

    One stable argsort instead of one full scan per cluster. At 984,119 cells and
    120 clusters the scan-per-cluster form costs 118 million comparisons per
    annotation; this costs one sort.

    :return: ``(names, groups)``, two lists of the same length.
    '''
    codes, uniques = pd.factorize(pd.Index(labels), sort=True)
    if (codes < 0).any():                      # pandas codes NaN as -1
        raise ValueError('annotation holds missing values; fill them before downsampling')
    counts = np.bincount(codes, minlength=len(uniques))
    order = np.argsort(codes, kind='stable')
    edges = np.concatenate([[0], np.cumsum(counts)])
    groups = [order[edges[i]:edges[i + 1]] for i in range(len(uniques))]
    return list(uniques), groups


def downsample_indices(obs, keys, max_per_cluster, seed=0, report=True):
    '''
    Row positions of a shared whitelist holding at most ``max_per_cluster`` cells
    per cluster, for every cluster of every annotation in ``keys``.

    Annotations are processed in the order given. Each one only draws the cells
    the earlier annotations did not already contribute. See the module docstring.

    :param obs: ``adata.obs``.
    :param keys: annotation column names, in priority order.
    :param max_per_cluster: cap per cluster. ``None`` returns every row.
    :param seed: fixes the draw. The same seed, obs and key order give the same
                 whitelist.

    :return: ``(rows, info)``. ``rows`` is a sorted int array of positions into
             ``obs``, so a CSR row slice stays cheap and the cell order is kept.
    '''
    t0 = time.perf_counter()
    n_cells = obs.shape[0]
    if max_per_cluster is None:
        return np.arange(n_cells), {'max_per_cluster': None, 'cells_before': int(n_cells),
                                    'cells_after': int(n_cells), 'reduction': 1.0,
                                    'seed': int(seed), 'per_annotation': {},
                                    'seconds': time.perf_counter() - t0}
    if max_per_cluster < 1:
        raise ValueError('max_per_cluster must be at least 1, got {}'.format(max_per_cluster))

    rng = np.random.default_rng(seed)
    keep = np.zeros(n_cells, dtype=bool)
    per_key = {}

    for key in keys:
        n_missing = int(pd.isna(obs[key]).sum())
        if n_missing:
            # astype(str) would turn these into a cluster literally named 'nan'
            # and then hand it a quota of its own
            raise ValueError('{} of {} cells have no {} label; fill or drop them '
                             'before downsampling'.format(n_missing, n_cells, key))
        names, groups = _group_positions(obs[key].astype(str).to_numpy())
        added = 0
        already_full = 0
        smaller_than_cap = 0
        for name, idx in zip(names, groups):
            if len(idx) <= max_per_cluster:
                smaller_than_cap += 1
                added += int((~keep[idx]).sum())       # only the ones not already in
                keep[idx] = True                       # take the whole cluster
                continue
            need = max_per_cluster - int(keep[idx].sum())
            if need <= 0:
                already_full += 1
                continue
            available = idx[~keep[idx]]
            take = rng.choice(available, size=need, replace=False)
            keep[take] = True
            added += need
        per_key[key] = {'clusters': len(names),
                        'clusters_smaller_than_cap': smaller_than_cap,
                        'clusters_already_at_cap': already_full,
                        'cells_added_here': int(added),
                        'whitelist_after': int(keep.sum())}
        if report:
            logger_sctriangulate.info(
                'downsample: {} has {} clusters, {} already at the cap, added {} cells, '
                'whitelist now {}'.format(key, len(names), already_full, added, int(keep.sum())))

    rows = np.flatnonzero(keep)
    info = {'max_per_cluster': int(max_per_cluster),
            'seed': int(seed),
            'keys_in_order': list(keys),
            'cells_before': int(n_cells),
            'cells_after': int(len(rows)),
            'reduction': (n_cells / len(rows)) if len(rows) else float('nan'),
            'per_annotation': per_key,
            'seconds': time.perf_counter() - t0}
    if report:
        logger_sctriangulate.info(
            'downsample: kept {} of {} cells, {:.1f}x fewer, in {:.1f} s'.format(
                info['cells_after'], info['cells_before'], info['reduction'], info['seconds']))
    return rows, info


def rank_markers(adata, key, max_cells_per_cluster=1000, seed=0):
    '''
    Gene-by-cluster MarkerFinder correlations on a per-cluster subsample.

    :return: ``(r_df, n_cells_used)``. ``r_df`` is genes x clusters.
    '''
    mf = _load_marker_finder()
    labels_all = adata.obs[key].astype(str).values
    rows = subsample_indices(labels_all, max_cells_per_cluster, seed)

    X = adata.X[rows] if issparse(adata.X) else csr_matrix(np.asarray(adata.X)[rows])
    X = X.tocsr() if issparse(X) else csr_matrix(X)
    labels = labels_all[rows]

    r_df, _p_df = mf.marker_finder(X, list(labels), gene_names=list(adata.var_names))
    return r_df, len(rows)


def select_marker_genes(adata, keys, n_per_cluster=200, max_cells_per_cluster=1000,
                        seed=0, report=True):
    '''
    Union of the top ``n_per_cluster`` MarkerFinder genes over every cluster of
    every annotation in ``keys``.

    :return: ``(genes, info)`` where ``genes`` is a list in ``adata.var_names``
             order and ``info`` is a dict of what happened, for the run record.
    '''
    t0 = time.perf_counter()
    var_names = pd.Index(adata.var_names)
    selected = set()
    per_key = {}

    for key in keys:
        r_df, n_used = rank_markers(adata, key, max_cells_per_cluster, seed)
        n_clusters = r_df.shape[1]
        take = min(int(n_per_cluster), r_df.shape[0])
        before = len(selected)
        for cluster in r_df.columns:
            col = r_df[cluster].to_numpy(dtype=np.float64)
            col = np.where(np.isfinite(col), col, -np.inf)
            # argpartition then sort only the top slice: O(G) instead of O(G log G)
            top = np.argpartition(-col, take - 1)[:take] if take < len(col) else np.arange(len(col))
            selected.update(var_names[top])
        per_key[key] = {'clusters': int(n_clusters), 'cells_scored': int(n_used),
                        'genes_added': len(selected) - before}
        logger_sctriangulate.info(
            'prefilter: {} scored {} clusters on {} cells, running union {} genes'.format(
                key, n_clusters, n_used, len(selected)))

    genes = [g for g in var_names if g in selected]     # keep var_names order
    info = {'method': 'markerfinder',
            'n_per_cluster': int(n_per_cluster),
            'max_cells_per_cluster': (None if max_cells_per_cluster is None
                                      else int(max_cells_per_cluster)),
            'seed': int(seed),
            'genes_before': int(adata.shape[1]),
            'genes_after': len(genes),
            'reduction': (adata.shape[1] / len(genes)) if genes else float('nan'),
            'per_annotation': per_key,
            'seconds': time.perf_counter() - t0}
    if report:
        logger_sctriangulate.info(
            'prefilter: {} of {} genes kept, {:.1f}x fewer, in {:.1f} s'.format(
                info['genes_after'], info['genes_before'], info['reduction'], info['seconds']))
    return genes, info
