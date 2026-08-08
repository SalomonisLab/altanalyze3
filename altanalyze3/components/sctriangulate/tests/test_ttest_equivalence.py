"""Check the chunked group-statistics patch against scanpy's own _basic_stats.

The patch reassociates the same sums, so means/vars/scores agree to float64
rounding rather than bitwise. What has to hold exactly is what the pipeline
actually consumes: the per-cluster marker gene lists, and therefore the purified
lists and the gene pool that reassign_score builds.

Run:
  cd /Users/saljh8/Documents/GitHub/altanalyze3
  python3.11 -m pytest altanalyze3/components/sctriangulate/tests/test_ttest_equivalence.py -q -s
"""

import os
import sys

import numpy as np
import pytest
import scanpy as sc

HERE = os.path.dirname(os.path.abspath(__file__))
COMPONENT = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(COMPONENT, '..', '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from altanalyze3.components.sctriangulate import ScTriangulate                  # noqa: E402
from altanalyze3.components.sctriangulate import metrics as M                   # noqa: E402

DEMO = os.path.join(COMPONENT, 'benchmarks', 'data', 'demo_pbmc3k.h5ad')
KEYS = ['sctri_rna_leiden_1', 'sctri_rna_leiden_2', 'sctri_rna_leiden_3']


@pytest.fixture(scope='module')
def dense_adata():
    if not os.path.exists(DEMO):
        pytest.skip('demo dataset not present at {}'.format(DEMO))
    adata = sc.read(DEMO)
    sctri = ScTriangulate(dir=os.path.join(COMPONENT, 'benchmarks', 'out_tests'),
                          adata=adata, query=KEYS, predict_doublet=False, verbose=0)
    sctri._to_dense()
    return sctri.adata


def _run(adata, key, fast):
    ad = M.check_filter_single_cluster(adata, key)
    if ad.uns.get('rank_genes_groups') is not None:
        del ad.uns['rank_genes_groups']
    with M._patched_basic_stats(enabled=fast):
        sc.tl.rank_genes_groups(ad, key, method='t-test', n_genes=ad.shape[1])
    uns = ad.uns['rank_genes_groups']
    return {f: np.asarray(uns[f].tolist()) for f in
            ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']}


@pytest.mark.parametrize('key', KEYS)
def test_marker_gene_lists_identical(dense_adata, key):
    """The only thing the pipeline reads downstream must be bit-identical."""
    slow = M.marker_gene(dense_adata, key, 'human', 2, '.', run_enrichment=False)
    M.FAST_TTEST, saved = False, M.FAST_TTEST
    try:
        M.FAST_TTEST = False
        ref = M.marker_gene(dense_adata, key, 'human', 2, '.', run_enrichment=False)
    finally:
        M.FAST_TTEST = saved

    n_match = 0
    for cluster in ref.index:
        if list(ref.loc[cluster, 'whole_marker_genes']) == list(slow.loc[cluster, 'whole_marker_genes']):
            n_match += 1
    print('\n{}: marker gene lists identical for {}/{} clusters'.format(
        key, n_match, len(ref.index)))
    assert n_match == len(ref.index)

    n_match = sum(1 for c in ref.index
                  if list(ref.loc[c, 'purify']) == list(slow.loc[c, 'purify']))
    assert n_match == len(ref.index)


@pytest.mark.parametrize('key', KEYS)
def test_numeric_agreement(dense_adata, key):
    """Report how far the reassociated statistics drift. Must be tiny."""
    fast = _run(dense_adata, key, fast=True)
    slow = _run(dense_adata, key, fast=False)

    assert (fast['names'] == slow['names']).all(), 'gene ordering changed'

    worst = {}
    for f in ['scores', 'pvals', 'pvals_adj', 'logfoldchanges']:
        a, b = fast[f].astype(float), slow[f].astype(float)
        finite = np.isfinite(a) & np.isfinite(b)
        denom = np.maximum(np.abs(b[finite]), 1e-30)
        worst[f] = float((np.abs(a[finite] - b[finite]) / denom).max())
    print('\n{}: max relative drift {}'.format(
        key, {k: '{:.2e}'.format(v) for k, v in worst.items()}))
    assert worst['scores'] < 1e-6
    assert worst['logfoldchanges'] < 1e-6
