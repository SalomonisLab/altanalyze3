"""Unit-level equivalence for the rewritten metric internals.

Each test pairs an optimized function against the upstream algorithm it
replaces. The upstream algorithms are reproduced verbatim in this file, copied
from `_reference/upstream_metrics.py` (scTriangulate 0.13.0, commit 8b9598cf),
so the comparison does not depend on upstream being importable.

Run:
  cd /Users/saljh8/Documents/GitHub/altanalyze3
  python3.11 -m pytest altanalyze3/components/sctriangulate/tests/test_metrics_equivalence.py -q
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
COMPONENT = os.path.dirname(HERE)
REPO = os.path.abspath(os.path.join(COMPONENT, '..', '..', '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from altanalyze3.components.sctriangulate import metrics as M          # noqa: E402
from altanalyze3.components.sctriangulate import main_class as MC      # noqa: E402


# ---------------------------------------------------------------------------
# verbatim upstream algorithms (see _reference/upstream_metrics.py)
# ---------------------------------------------------------------------------

def upstream_assign_loop(all_genes, all_clusters, pre_computed_dfs):
    """upstream_metrics.py :: marker_gene, the gene-to-cluster assignment block."""
    cluster2gene = {c: [] for c in all_clusters}
    for gene in all_genes:
        index_store = []
        for i, cluster in enumerate(all_clusters):
            df = pre_computed_dfs[i]
            try:
                index = np.nonzero(df['names'].values == gene)[0][0]
            except IndexError:
                index = len(all_genes)
            index_store.append(index)
        if np.all(np.array(index_store) == len(all_genes)):
            continue
        assign = all_clusters[np.argmin(np.array(index_store))]
        cluster2gene[assign].append((gene, np.min(index_store)))
    for key_, value in cluster2gene.items():
        gene = [item[0] for item in value]
        rank = [item[1] for item in value]
        temp = sorted(zip(gene, rank), key=lambda x: x[1])
        cluster2gene[key_] = [item[0] for item in temp]
    return cluster2gene


def upstream_tf_idf_bare_compute(df, cluster):
    """upstream_metrics.py :: tf_idf_bare_compute, verbatim."""
    tmp1 = df.loc[df['cluster'] == cluster, :].loc[:, df.columns != 'cluster'].values
    tf = np.count_nonzero(tmp1, axis=0) / tmp1.shape[0]
    tf = tf + 1e-5
    tmp2 = df.loc[:, df.columns != 'cluster'].values
    df_ = np.count_nonzero(tmp2, axis=0) / tmp2.shape[0]
    df_ = df_ + 1e-5
    idf = -np.log10(df_)
    return tf * idf


def upstream_doublet_compute(adata, key):
    """upstream_metrics.py :: doublet_compute, verbatim."""
    cluster_to_doublet = {}
    for cluster in adata.obs[key].astype('category').cat.categories:
        mean_score = adata[adata.obs[key] == cluster, :].obs['doublet_scores'].values.mean()
        cluster_to_doublet[cluster] = mean_score
    return cluster_to_doublet


def upstream_penalize_artifact_void(obs, query, stamps, metrics):
    """main_class.py :: penalize_artifact_void, verbatim."""
    for stamp in stamps:
        cols = [item2 + '@' + item1 for item1 in query for item2 in metrics]
        metrics_cols = obs.loc[:, cols]
        cluster_cols = obs.loc[:, query]
        df = cluster_cols.apply(
            func=lambda x: pd.Series(data=[x.name + '@' + str(item) for item in x], name=x.name),
            axis=0)
        df_repeat = pd.DataFrame(np.repeat(df.values, len(metrics), axis=1))
        truth = pd.DataFrame(data=(df_repeat == stamp).values,
                             index=metrics_cols.index, columns=metrics_cols.columns)
        obs.loc[:, cols] = metrics_cols.mask(truth, 0)
    return obs


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

def _fake_rank_tables(rng, all_genes, n_clusters, keep_frac=0.6):
    """One combo-score table per cluster: a random subset of genes in a random order."""
    dfs = []
    for _ in range(n_clusters):
        n_keep = int(len(all_genes) * keep_frac)
        picked = rng.choice(all_genes, size=n_keep, replace=False)
        dfs.append(pd.DataFrame({'names': picked}))
    return dfs


@pytest.mark.parametrize('n_clusters', [1, 3, 8])
@pytest.mark.parametrize('seed', [0, 1])
def test_assign_genes_matches_upstream_loop(n_clusters, seed):
    rng = np.random.default_rng(seed)
    all_genes = np.array(['G{:04d}'.format(i) for i in range(220)], dtype=object)
    all_clusters = pd.Index(['c{}'.format(i) for i in range(n_clusters)])
    dfs = _fake_rank_tables(rng, all_genes, n_clusters)

    expected = upstream_assign_loop(all_genes, all_clusters, dfs)
    got = M._assign_genes_to_clusters(all_genes, all_clusters, dfs)

    assert set(got) == set(expected)
    for c in all_clusters:
        assert list(got[c]) == list(expected[c]), 'cluster {}'.format(c)


def test_assign_genes_handles_empty_and_missing():
    all_genes = np.array(['A', 'B', 'C', 'D'], dtype=object)
    all_clusters = pd.Index(['c0', 'c1'])
    dfs = [pd.DataFrame({'names': np.array([], dtype=object)}),   # empty table
           pd.DataFrame({'names': np.array(['C', 'A'], dtype=object)})]
    expected = upstream_assign_loop(all_genes, all_clusters, dfs)
    got = M._assign_genes_to_clusters(all_genes, all_clusters, dfs)
    for c in all_clusters:
        assert list(got[c]) == list(expected[c])


@pytest.mark.parametrize('sparse_input', [False, True])
def test_tfidf_matrix_matches_upstream(sparse_input):
    import anndata as ad
    from scipy.sparse import csr_matrix

    rng = np.random.default_rng(4)
    n_cells, n_genes, n_clusters = 300, 120, 6
    X = rng.random((n_cells, n_genes)).astype(np.float32)
    X[X < 0.55] = 0.0                                   # realistic sparsity
    labels = rng.integers(0, n_clusters, n_cells).astype(str)
    obs = pd.DataFrame({'k': pd.Categorical(labels)},
                       index=['cell{}'.format(i) for i in range(n_cells)])
    var = pd.DataFrame(index=['g{}'.format(i) for i in range(n_genes)])
    adata = ad.AnnData(X=csr_matrix(X) if sparse_input else X, obs=obs, var=var)

    M._TFIDF_CACHE.update({'adata': None, 'matrix': None, 'ranked': None})
    clusters, matrix = M._tfidf_matrix(adata, 'k')

    df = pd.DataFrame(data=X, index=adata.obs_names, columns=adata.var_names)
    df['cluster'] = adata.obs['k'].astype('str').values
    for i, c in enumerate(clusters):
        expected = upstream_tf_idf_bare_compute(df, c)
        assert np.allclose(matrix[i], expected, rtol=0, atol=1e-12), \
            'cluster {}: max diff {}'.format(c, np.abs(matrix[i] - expected).max())


def test_doublet_compute_matches_upstream():
    import anndata as ad
    rng = np.random.default_rng(5)
    n = 200
    obs = pd.DataFrame({'k': pd.Categorical(rng.integers(0, 5, n).astype(str)),
                        'doublet_scores': rng.random(n)},
                       index=['c{}'.format(i) for i in range(n)])
    adata = ad.AnnData(X=np.zeros((n, 3), dtype=np.float32), obs=obs,
                       var=pd.DataFrame(index=['a', 'b', 'c']))
    expected = upstream_doublet_compute(adata, 'k')
    got = M.doublet_compute(adata, 'k')
    assert set(got) == set(expected)
    for c in expected:
        assert abs(got[c] - expected[c]) < 1e-12


def test_run_enrichr_matches_gseapy_local_mode():
    """The inlined hypergeometric must reproduce gseapy 0.10.4's local enrichr."""
    gp = pytest.importorskip('gseapy')
    artifact = M._artifact_class_dict('human')
    all_artifact = sorted({g for v in artifact.values() for g in v})
    rng = np.random.default_rng(9)

    for n_hits, n_other in [(60, 200), (5, 400), (0, 300), (300, 50)]:
        genes = list(rng.choice(all_artifact, size=n_hits, replace=False)) + \
                ['NOTAGENE{}'.format(i) for i in range(n_other)]
        got = M.run_enrichr(genes, 'k', 'c', '.', 'human', 2)

        enr = gp.enrichr(gene_list=genes, description='x', gene_sets=dict(artifact),
                         background=20000, outdir=None, cutoff=0.1,
                         no_plot=True, verbose=False)
        expected = {m: 0 for m in artifact}
        if enr.results.shape[0] > 0:
            for _, row in enr.results.iterrows():
                if row['Term'] in expected:
                    expected[row['Term']] = -np.log10(row['Adjusted P-value'])

        for term in expected:
            assert abs(float(got[term]) - float(expected[term])) < 1e-12, \
                '{} hits: term {} got {} want {}'.format(n_hits, term, got[term], expected[term])


def test_penalize_artifact_void_matches_upstream():
    rng = np.random.default_rng(12)
    n = 150
    query = ['a1', 'a2']
    metrics = ['reassign', 'tfidf10', 'SCCAF']
    data = {k: rng.integers(0, 4, n).astype(str) for k in query}
    for k in query:
        for m in metrics:
            data['{}@{}'.format(m, k)] = rng.random(n)
    obs = pd.DataFrame(data, index=['c{}'.format(i) for i in range(n)])

    stamps = ['a1@2', 'a2@0']
    expected = upstream_penalize_artifact_void(obs.copy(), query, stamps, metrics)
    got = MC.penalize_artifact_void(obs.copy(), query, stamps, metrics)

    cols = [m + '@' + k for k in query for m in metrics]
    assert np.allclose(got[cols].to_numpy(dtype=float),
                       expected[cols].to_numpy(dtype=float), rtol=0, atol=0)
    assert (got[cols].to_numpy(dtype=float) == 0).sum() > 0, 'test data stamped nothing'


def test_read_artifact_genes_cache_matches_fresh_read():
    """The cache must not change what the criteria select."""
    path = os.path.join(COMPONENT, 'artifact_genes.txt')
    raw = pd.read_csv(path, sep='\t', index_col=0)
    for species in ['human', 'mouse']:
        for criterion in range(1, 7):
            cached = M.read_artifact_genes(species, criterion)
            fresh = raw.loc[raw['species'] == species, :]
            if criterion == 2:
                fresh = fresh.loc[~(fresh['class'] == 'cellcycle'), :]
            elif criterion == 3:
                fresh = fresh.loc[~((fresh['class'] == 'ribosome') | (fresh['class'] == 'cellcycle')), :]
            elif criterion == 4:
                fresh = fresh.loc[~((fresh['class'] == 'ribosome') | (fresh['class'] == 'cellcylce') | (fresh['class'] == 'mitochondrial')), :]
            elif criterion == 5:
                fresh = fresh.loc[~((fresh['class'] == 'ribosome') | (fresh['class'] == 'cellcylce') | (fresh['class'] == 'mitochondrial') | (fresh['class'] == 'antisense')), :]
            elif criterion == 6:
                fresh = fresh.loc[~((fresh['class'] == 'ribosome') | (fresh['class'] == 'cellcylce') | (fresh['class'] == 'mitochondrial') | (fresh['class'] == 'antisense') | (fresh['class'] == 'predict_gene')), :]
            assert list(cached.index) == list(fresh.index), (species, criterion)


def test_get_size_matches_upstream():
    from altanalyze3.components.sctriangulate import shapley as S
    rng = np.random.default_rng(21)
    n = 500
    query = ['a1', 'a2', 'a3']
    obs = pd.DataFrame({k: rng.integers(0, 7, n).astype(str) for k in query})
    size_dict, size_list = S.get_size(obs, query)
    for k in query:
        for cluster, size in size_dict[k].items():
            assert size == int((obs[k] == cluster).sum())
    assert len(size_list) == sum(len(v) for v in size_dict.values())
