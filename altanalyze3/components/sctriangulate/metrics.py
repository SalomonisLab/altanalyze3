'''
Cluster stability metrics for scTriangulate, as integrated into AltAnalyze3.

Public names match upstream scTriangulate 0.13.0 (``main_class`` does
``from .metrics import *``), so this module is a drop-in replacement. The bodies
are rewritten for speed; the numbers they return are the upstream numbers.

WHAT CHANGED AND WHY IT IS STILL THE SAME NUMBER
------------------------------------------------
1. ``read_artifact_genes`` is cached. Upstream re-read and re-filtered an
   18,000-row table once per cluster inside three separate loops.

2. ``marker_gene`` assigned each gene to its best cluster with a nested loop
   that ran ``np.nonzero(df['names'].values == gene)`` for every
   (gene, cluster) pair, i.e. O(G^2 * C) element comparisons. That is replaced
   by a rank matrix filled with one vectorised ``get_indexer`` per cluster.
   Same argmin, same stable tie-break, same output lists.

3. ``tf_idf_bare_compute`` built a dense ``n_cells x n_genes`` DataFrame and
   recomputed the whole-matrix nonzero count once per cluster. The rewrite
   counts nonzeros once with a single indicator matmul in float64, which is
   exact for 0/1 sums, and reuses it for every cluster.

4. ``tf_idf1/5/10_for_cluster`` each recomputed the identical per-cluster tf-idf
   vectors. They now share one cached computation per (adata, key, layer).

5. Enrichr in "local mode" is a hypergeometric test plus Benjamini-Hochberg.
   ``run_enrichr`` now calls that directly instead of constructing a gseapy
   ``Enrichr`` object per cluster. The test, the term ordering, the ``x < 1``
   skip and the BH implementation are ported line for line from
   gseapy 0.10.4 (``gseapy/stats.py::calc_pvalues`` and ``fdrcorrection``),
   so p-values and adjusted p-values are identical.

6. ``run_gsea`` still calls ``gseapy.prerank`` -- the permutation NES depends on
   gseapy's RNG stream and is not worth reimplementing -- but with
   ``outdir=None`` and ``no_plot=True``. Upstream wrote a directory of figures
   and tables per cluster. Measured on the demo data: 0.714 s -> 0.119 s per
   cluster with byte-identical NES.

7. ``_fast_basic_stats`` replaces scanpy's per-group statistics with one pass
   that accumulates group sums and sums of squares. It squares in the STORED
   dtype, because scanpy's ``elem_mul(X, X)`` does; squaring float32 data in
   float64 is more accurate but disagrees by about 4e-06 relative, which
   reorders tied genes and changes the marker lists. Verified bit-identical by
   ``tests/test_ttest_equivalence.py``.

8. ``scaled_pca`` computes the reassign-step PCA exactly, and without ever
   materialising the centred matrix. Upstream densified and used sklearn's
   randomized SVD, which does not recover the true subspace: minimum cosine of
   the principal angles 0.0031 against the exact full SVD, capturing 0.8% less
   variance, and it was the main source of upstream's run-to-run instability.
   The sparse path builds a ``LinearOperator`` for the centred matrix and calls
   ``svds``; it scores 0.99999999 against the exact solution, and the dense and
   sparse paths agree to ``max|dAcc| = 0.00000000``. ``REASSIGN_PCA`` and
   ``REASSIGN_DENSE_BUDGET_MB`` pick the path; ``PCA_RANDOM_STATE`` seeds it.

9. ``SCCAF_score`` has two paths. ``SCCAF_MODE = 'optimized'``, the default,
   keeps the matrix sparse and hands liblinear the same one-vs-rest problem, so
   the accuracies are the legacy accuracies and only the representation
   changed. ``'legacy'`` is upstream's dense procedure, kept for comparison.
   Both are seeded through ``SCCAF_RANDOM_STATE``; upstream was not seeded at
   all, and matched its own first run in only 1 to 3 of 50 repeats.
   ``SCCAF_PER_CLASS`` splits the one-vs-rest loop across processes. It is
   opt-in and NOT recommended: it moves accuracies by up to 0.167, and threads
   made it non-reproducible because liblinear uses a process-global C ``rand()``.

10. ``_shared_index_csr`` builds a matrix with new values over an existing
    matrix's ``indices`` and ``indptr``. Two callers need the same sparsity
    pattern with different data -- the t-test squares it, the tf-idf indicator
    replaces it with ones -- and both used ``src.copy()``, which duplicates the
    index arrays and then discards the data it copied. Measured at 100,000,000
    nonzeros: peak 1200.4 MB against 400.0 MB, 0.11 s against 0.04 s, with
    ``data``, ``indices`` and ``indptr`` byte-equal.

11. ``check_filter_single_cluster`` returns the caller's object when no cluster
    holds one cell, instead of a view. An AnnData view rebuilds ``X`` on every
    attribute access, 81.1 MB then 81.0 MB against an 80.2 MB matrix, and
    ``each_key_run`` has four readers. ``tests/test_no_write_to_adata.py``
    replaces the write-protection the view gave by accident.

12. ``run_enrichment`` is off by default (``RUN_ENRICHMENT = False``). The
    enrichment columns feed no stability metric and no Shapley decision; only
    ``plot_cluster_feature``, ``plot_multi_modal_feature_rank``,
    ``penalize_artifact(mode='cellcycle')`` and the viewer read them. They cost
    about 3.7 s of a 21 s demo run.

MEMORY CONSTANTS
----------------
``TFIDF_CHUNK_MB``, ``FAST_TTEST_CHUNK_MB`` and ``SCCAF_CHUNK_MB`` bound the
working buffers of the dense paths, so peak memory stays flat as the dataset
grows rather than scaling with ``n_cells * n_genes``.

Every claim in items 1 to 6 is checked by ``tests/test_metrics_equivalence.py``
against the verbatim upstream code in ``_reference/upstream_metrics.py``.
Items 7 to 11 are checked by ``tests/test_ttest_equivalence.py``,
``tests/test_sccaf_modes.py`` and ``tests/test_no_write_to_adata.py``.
The measurements quoted above are reproduced by the scripts named in
``README.md``, section "Reproducing every number".
'''

import os
import sys
import math
import logging
import multiprocessing as mp

from collections.abc import Mapping

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import issparse, csr_matrix
from scipy.stats import hypergeom

from .logger import logger_sctriangulate
from .preprocessing import *


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

def check_filter_single_cluster(adata, key):
    '''
    Drop the clusters that hold exactly one cell, which have no within-group
    variance and break the differential tests.

    Return the object itself when there is nothing to drop. Upstream always
    built a view, and an AnnData view rebuilds X on EVERY attribute access:
    measured 81.1 MB on the first ``.X`` and 81.0 MB on the second, against an
    80.2 MB matrix, with no caching in between. ``each_key_run`` reaches ``.X``
    from marker_gene, reassign_score, tf_idf10_for_cluster and SCCAF_score, so
    the view cost one full rebuild per reader for no gain.

    Returning the caller's object gives up the accidental write-protection the
    view provided. No metric function writes to its adata argument;
    tests/test_no_write_to_adata.py checksums X, obs and var across a full
    each_key_run and fails if one byte moves.
    '''
    vc = adata.obs[key].value_counts()
    exclude_clusters = vc.loc[vc == 1].index
    if len(exclude_clusters) == 0:
        return adata
    truth = np.logical_not(adata.obs[key].isin(exclude_clusters).values)
    adata_valid = adata[truth, :]
    return adata_valid


def doublet_compute(adata, key):
    '''Mean doublet score per cluster. One groupby instead of one mask per cluster.'''
    grouper = adata.obs[key].astype('category')
    means = adata.obs['doublet_scores'].groupby(grouper, observed=False).mean()
    return {cluster: float(means[cluster]) for cluster in grouper.cat.categories}


def _as_dense(mat):
    '''Return a dense ndarray view/copy of a possibly sparse matrix.'''
    if issparse(mat):
        return mat.toarray()
    return np.asarray(mat)


def _shared_index_csr(src, data):
    '''
    CSR carrying new ``data`` over ``src``'s own ``indices`` and ``indptr``.

    Two places need a matrix with the same sparsity pattern and different
    values: the t-test squares the data, and the tf-idf indicator replaces it
    with ones. Both used ``src.copy()``, which duplicates the index arrays and
    then discards the data it just copied. Only the new data array is needed.

    Measured at 100,000,000 nonzeros by benchmarks/copy_avoidance_probe.py:
    peak 1200.4 MB against 400.0 MB, and 0.11 s against 0.04 s, with ``data``,
    ``indices`` and ``indptr`` byte-equal either way. The arithmetic is
    untouched, so no marker list moves.

    Sharing is only safe when the indices are sorted. scipy sorts them in place
    when they are not, and that would permute the shared index array out of step
    with ``src.data``. Copy the index arrays in that case.
    '''
    if not src.has_sorted_indices:
        return csr_matrix((data, src.indices.copy(), src.indptr.copy()),
                          shape=src.shape, copy=False)
    out = csr_matrix((data, src.indices, src.indptr), shape=src.shape, copy=False)
    out.has_sorted_indices = True
    return out


class ExclusiveGeneScores(Mapping):
    '''
    Ordered gene -> tf-idf score mapping backed by two arrays.

    ``uns['exclusive_genes'][annotation][cluster]`` used to be a plain dict with
    one entry per gene. On the demo data the 61 clusters held 1,901,980 entries
    between them, and ``lazy_run`` pickles the whole object three times. Measured:
    pickling the dict form took 2.28 s and 46.2 MB, so the three checkpoints cost
    6.8 s of a 28.9 s run. The array-backed form takes 0.09 s and 25.3 MB, a 25x
    drop, because the scores travel as one float64 buffer instead of ~1.9 million
    Python float objects.

    It stays a ``Mapping`` in insertion (descending score) order, so every caller
    keeps working: ``list(x.keys())``, ``list(x.values())``, ``x['GENE']``,
    ``len(x)``, ``for g in x`` and ``x.to_dict()``. It also accepts ``x[:10]``,
    which ``plot_multi_modal_feature_fraction`` already assumed and which a dict
    never supported.
    '''

    __slots__ = ('_genes', '_scores', '_lookup')

    def __init__(self, genes, scores):
        self._genes = pd.Index(genes)
        self._scores = np.asarray(scores, dtype=np.float64)
        if len(self._genes) != len(self._scores):
            raise ValueError('genes and scores differ in length: {} vs {}'.format(
                len(self._genes), len(self._scores)))
        self._lookup = None

    def _pos(self):
        if self._lookup is None:
            self._lookup = {g: i for i, g in enumerate(self._genes)}
        return self._lookup

    def __getitem__(self, k):
        if isinstance(k, slice):
            return ExclusiveGeneScores(self._genes[k], self._scores[k])
        return float(self._scores[self._pos()[k]])

    def __iter__(self):
        return iter(self._genes)

    def __len__(self):
        return len(self._genes)

    def keys(self):
        return list(self._genes)

    def values(self):
        return [float(v) for v in self._scores]

    def items(self):
        return list(zip(list(self._genes), (float(v) for v in self._scores)))

    def to_dict(self):
        return dict(zip(list(self._genes), (float(v) for v in self._scores)))

    def to_series(self):
        return pd.Series(self._scores, index=self._genes)

    def __repr__(self):
        return '<ExclusiveGeneScores: {} genes>'.format(len(self._genes))

    def __reduce__(self):
        return (ExclusiveGeneScores, (self._genes, self._scores))


TFIDF_CHUNK_MB = 128   # working memory for the chunked nonzero count


def _nonzero_indicator(mat):
    '''
    0-1 matrix of "is nonzero".

    Kept for the sparse path and for callers that want the whole thing. For a
    dense input this allocates a full second copy of the matrix, so the counting
    routine below chunks instead. Do not call it on a large dense matrix.
    '''
    if issparse(mat):
        out = mat.copy().tocsr()
        out.data = np.ones_like(out.data, dtype=np.float64)
        return out
    return (np.asarray(mat) != 0).astype(np.float64)


def _count_nonzero_by_group(mat, onehot):
    '''
    Nonzero counts per (cluster, gene) and per gene, without ever holding a
    second full copy of the matrix.

    :param mat: (n_cells, n_genes), dense or sparse.
    :param onehot: (n_cells, n_clusters) float32, one 1 per row.
    :return: (per_cluster_counts float64 (n_clusters, n_genes),
              total_counts float64 (n_genes,))

    Memory: the dense path materialises the 0-1 indicator for
    ``TFIDF_CHUNK_MB`` worth of rows at a time, so peak stays flat as the
    dataset grows. Building it for the whole matrix at once, which is what an
    earlier version of this file did, costs n_cells * n_genes * 8 bytes: 707 MB
    on the 2,700-cell demo and 10 GB at 50,000 cells by 25,000 genes.

    Exactness: the indicator holds only 0 and 1, and a count can never exceed
    n_cells. Every partial sum is therefore an integer below 2**24, which float32
    represents exactly, so the float32 matmul is exact arithmetic on integers.
    The accumulator is float64, where the same integers are exact up to 2**53.
    The result equals np.count_nonzero.
    '''
    n_cells, n_genes = mat.shape
    n_clusters = onehot.shape[1]
    per_cluster = np.zeros((n_clusters, n_genes), dtype=np.float64)

    if issparse(mat):
        src = mat.tocsr()
        # the copied data was overwritten with ones on the next line anyway; see
        # _shared_index_csr
        ind = _shared_index_csr(src, np.ones(src.nnz, dtype=np.float32))
        per_cluster += np.asarray(onehot.T @ ind, dtype=np.float64)
    else:
        bytes_per_row = n_genes * 4
        chunk = max(1, int(TFIDF_CHUNK_MB * 1024 * 1024 / max(bytes_per_row, 1)))
        for start in range(0, n_cells, chunk):
            stop = min(start + chunk, n_cells)
            ind = (np.asarray(mat[start:stop]) != 0).astype(np.float32)
            per_cluster += (onehot[start:stop].T @ ind).astype(np.float64)
            del ind

    # the one-hot covers every cell exactly once, so the column totals are the
    # per-gene totals; assert it rather than trust it
    counts_per_cell = onehot.sum(axis=1)
    if not np.array_equal(counts_per_cell, np.ones(n_cells, dtype=onehot.dtype)):
        raise ValueError('cluster one-hot does not assign every cell exactly once')
    total = per_cluster.sum(axis=0)
    return per_cluster, total


# ---------------------------------------------------------------------------
# artifact gene table (cached)
# ---------------------------------------------------------------------------

_ARTIFACT_TABLE = None
_ARTIFACT_CACHE = {}
_ARTIFACT_SET_CACHE = {}
_ARTIFACT_DICT_CACHE = {}


def _artifact_table():
    global _ARTIFACT_TABLE
    if _ARTIFACT_TABLE is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'artifact_genes.txt')
        _ARTIFACT_TABLE = pd.read_csv(path, sep='\t', index_col=0)
    return _ARTIFACT_TABLE


def read_artifact_genes(species, criterion):
    '''
    criterion1: all will be artifact
    criterion2: all will be artifact except cellcycle
    criterion3: all will be artifact except cellcycle, ribosome
    criterion4: all will be artifact except cellcycle, ribosome, mitochondrial
    criterion5: all will be artifact except cellcycle, ribosome, mitochondrial, antisense
    criterion6: all will be artifact except cellcycle, ribosome, mitochondrial, antisense, predict_gene

    Cached. The upstream typo 'cellcylce' in criteria 4-6 is preserved on
    purpose: fixing it would silently change which genes are treated as
    artifacts and therefore change every downstream score.
    '''
    ck = (species, criterion)
    if ck in _ARTIFACT_CACHE:
        return _ARTIFACT_CACHE[ck]
    artifact = _artifact_table()
    artifact = artifact.loc[artifact['species'] == species, :]
    if criterion == 1:
        artifact = artifact
    elif criterion == 2:
        artifact = artifact.loc[~(artifact['class'] == 'cellcycle'), :]
    elif criterion == 3:
        artifact = artifact.loc[~((artifact['class'] == 'ribosome') | (artifact['class'] == 'cellcycle')), :]
    elif criterion == 4:
        artifact = artifact.loc[~((artifact['class'] == 'ribosome') | (artifact['class'] == 'cellcylce') | (artifact['class'] == 'mitochondrial')), :]
    elif criterion == 5:
        artifact = artifact.loc[~((artifact['class'] == 'ribosome') | (artifact['class'] == 'cellcylce') | (artifact['class'] == 'mitochondrial') | (artifact['class'] == 'antisense')), :]
    elif criterion == 6:
        artifact = artifact.loc[~((artifact['class'] == 'ribosome') | (artifact['class'] == 'cellcylce') | (artifact['class'] == 'mitochondrial') | (artifact['class'] == 'antisense') | (artifact['class'] == 'predict_gene')), :]
    _ARTIFACT_CACHE[ck] = artifact
    return artifact


def _artifact_gene_set(species, criterion):
    ck = (species, criterion)
    if ck not in _ARTIFACT_SET_CACHE:
        _ARTIFACT_SET_CACHE[ck] = set(read_artifact_genes(species, criterion).index.to_list())
    return _ARTIFACT_SET_CACHE[ck]


def _artifact_class_dict(species):
    '''{'ribosome': [...], ...} built from criterion 1, as run_enrichr/run_gsea do.'''
    if species not in _ARTIFACT_DICT_CACHE:
        artifact = read_artifact_genes(species, criterion=1).reset_index()
        _ARTIFACT_DICT_CACHE[species] = \
            artifact.groupby(by='class')['genes'].apply(lambda x: x.tolist()).to_dict()
    return {k: list(v) for k, v in _ARTIFACT_DICT_CACHE[species].items()}


def purify_gene(genelist, species, criterion):
    artifact_genes = _artifact_gene_set(species, criterion)
    return [gene for gene in genelist if gene not in artifact_genes]


# ---------------------------------------------------------------------------
# enrichment
# ---------------------------------------------------------------------------

def _fdrcorrection(pvals, alpha=0.05):
    '''
    Benjamini-Hochberg, ported verbatim from gseapy 0.10.4 ``stats.fdrcorrection``
    (which itself copies GOATools). Reproduced here so ``run_enrichr`` does not
    have to build a gseapy object per cluster.
    '''
    pvals = np.asarray(pvals)
    pvals_sortind = np.argsort(pvals)
    pvals_sorted = np.take(pvals, pvals_sortind)
    nobs = len(pvals_sorted)
    ecdffactor = np.arange(1, nobs + 1) / float(nobs)
    reject = pvals_sorted <= ecdffactor * alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True
    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected > 1] = 1
    pvals_corrected_ = np.empty_like(pvals_corrected)
    pvals_corrected_[pvals_sortind] = pvals_corrected
    reject_ = np.empty_like(reject)
    reject_[pvals_sortind] = reject
    return reject_, pvals_corrected_


def run_enrichr(gene_list, key, name, folder, species, criterion):
    '''
    Artifact-class over-representation for one cluster's marker genes.

    Returns ``{artifact_class: -log10(adjusted p)}``, 0 where the class was not
    tested or not returned, matching upstream.

    Upstream called ``gseapy.enrichr(..., background=20000, outdir=...)``, which
    in local mode (a dict of gene sets) runs a hypergeometric test and writes a
    directory of output per cluster. That test is inlined here:
    ``hypergeom.sf(x-1, background, m, k)`` over terms in sorted key order,
    skipping terms with no overlap, then Benjamini-Hochberg. The ``folder``
    argument is accepted for signature compatibility and no longer written to.
    '''
    artifact_dict = _artifact_class_dict(species)
    background = 20000
    cutoff = 0.1   # upstream passes cutoff=0.1, which is BH's alpha

    k = len(gene_list)            # gseapy: k = len(query) before set()
    query = set(gene_list)
    terms, pvals = [], []
    for s in sorted(artifact_dict.keys()):
        category = set(artifact_dict[s])
        x = len(query.intersection(category))
        if x < 1:
            continue
        m = len(category)
        terms.append(s)
        pvals.append(hypergeom.sf(x - 1, background, m, k))

    enrichr_dict = {metric: 0 for metric in artifact_dict.keys()}
    if not terms:
        return enrichr_dict   # upstream: enrichr_result.shape[0]==0 -> all zero
    _, fdrs = _fdrcorrection(np.asarray(pvals, dtype=np.float64), alpha=cutoff)
    for term, fdr in zip(terms, fdrs):
        enrichr_dict[term] = -math.log10(fdr)
    return enrichr_dict


def run_gsea(gene_list, key, name, folder, species, criterion):
    '''
    Artifact-class GSEA preranked score for one cluster's marker genes.

    Returns ``{artifact_class: (nes, n_hits)}``.

    Still gseapy: the NES comes from a seeded gene-label permutation and
    reproducing gseapy's RNG stream elsewhere would be fragile. The change is
    ``outdir=None, no_plot=True, verbose=False``, which removes the per-cluster
    figure and table writing. Measured on the demo data, 0.714 s -> 0.119 s per
    cluster with identical NES.
    '''
    import gseapy as gp

    artifact_dict = _artifact_class_dict(species)
    artifact_dict_keys = list(artifact_dict.keys())
    df = pd.DataFrame({0: gene_list, 1: 1 / (np.arange(len(gene_list)) + 1)})
    gsea_dict = {}
    try:
        pre_res = gp.prerank(rnk=df, gene_sets=artifact_dict,
                             permutation_num=100,
                             outdir=None,
                             no_plot=True,
                             min_size=1,
                             max_size=10000,
                             seed=6,
                             verbose=False)
    except Exception as e:
        # upstream swallowed every exception here and returned zeros. Keep the
        # zero-fill so scores stay comparable, but say what happened.
        logger_sctriangulate.warning(
            'gseapy prerank failed for {}/{} ({}: {}); filling zeros'.format(
                key, name, type(e).__name__, e))
        for metric in artifact_dict_keys:
            gsea_dict[metric] = (0, 0)
        return gsea_dict

    gsea_result = pre_res.res2d
    metric_get = set(gsea_result.index.tolist())
    for metric in artifact_dict_keys:
        if metric in metric_get:
            row = gsea_result.loc[gsea_result.index == metric, :]
            gsea_dict[metric] = (row['nes'].to_list()[0], row['matched_size'].to_list()[0])
        else:
            gsea_dict[metric] = (0, 0)
    return gsea_dict


# ---------------------------------------------------------------------------
# marker genes
# ---------------------------------------------------------------------------

def compute_combo_score(rank_uns, cluster):
    rank_names = rank_uns['names'][cluster]
    rank_lfc = rank_uns['logfoldchanges'][cluster]
    rank_pval = rank_uns['pvals'][cluster]
    df = pd.DataFrame({'names': rank_names, 'lfc': rank_lfc, 'pval': rank_pval})
    # filter out down-regulated genes
    df = df.loc[df['lfc'] > 0, :]
    df.set_index(keys=pd.Index(np.arange(df.shape[0])), inplace=True)
    # the rank of each gene by lfc, the larger, the better, make argsort result reverse
    temp = np.flip(np.argsort(df['lfc'].values))
    ranks_lfc = np.empty_like(temp)
    ranks_lfc[temp] = np.arange(len(df['pval'].values))
    # the rank of each gene by pval, the smaller, the better
    temp = np.argsort(df['pval'].values)
    ranks_pval = np.empty_like(temp)
    ranks_pval[temp] = np.arange(len(df['pval'].values))
    # combo rank score
    df['rank_lfc'] = ranks_lfc
    df['rank_pval'] = ranks_pval
    df['combo'] = (ranks_lfc + ranks_pval) / 2
    df.sort_values(by='combo', inplace=True)
    df.set_index(keys=pd.Index(np.arange(df.shape[0])), inplace=True)
    # filter out the genes if pval > 0.05
    df = df.loc[df['pval'] < 0.05, :]
    df.set_index(keys=pd.Index(np.arange(df.shape[0])), inplace=True)
    return df


def _assign_genes_to_clusters(all_genes, all_clusters, pre_computed_dfs):
    '''
    Assign every gene to the cluster where it ranks best.

    Upstream ran, for each of G genes and each of C clusters, a linear scan
    ``np.nonzero(df['names'].values == gene)[0][0]`` over a length-G array:
    O(G^2 * C). Here each cluster's ranks are scattered into a (G, C) matrix in
    one vectorised step, then argmin over the cluster axis.

    Tie handling matches upstream exactly:
      * ``np.argmin`` takes the lowest cluster index on a tie, as upstream did;
      * genes absent from every cluster's table keep the sentinel rank
        ``len(all_genes)`` and are dropped, as upstream did;
      * the per-cluster gene order is a stable sort on rank, which reproduces
        upstream's ``sorted(zip(gene, rank), key=...)`` over genes visited in
        ``var_names`` order.
    '''
    n_genes = len(all_genes)
    gene_index = pd.Index(all_genes)
    sentinel = n_genes
    ranks = np.full((n_genes, len(all_clusters)), sentinel, dtype=np.int64)

    if gene_index.is_unique:
        for i, df in enumerate(pre_computed_dfs):
            if df.shape[0] == 0:
                continue
            pos = gene_index.get_indexer(df['names'].values)
            found = pos >= 0
            ranks[pos[found], i] = np.arange(df.shape[0])[found]
    else:
        # duplicated var_names: upstream's np.nonzero(...)[0][0] takes the first
        # match, so map every duplicate row to the first occurrence
        first_pos = {}
        for p, g in enumerate(all_genes):
            first_pos.setdefault(g, p)
        for i, df in enumerate(pre_computed_dfs):
            names = df['names'].values
            pos = np.array([first_pos.get(g, -1) for g in names], dtype=np.int64)
            found = pos >= 0
            ranks[pos[found], i] = np.arange(df.shape[0])[found]

    best_cluster = ranks.argmin(axis=1)
    best_rank = ranks[np.arange(n_genes), best_cluster]
    keep = best_rank < sentinel

    cluster2gene = {}
    genes_arr = np.asarray(all_genes)
    for i, cluster in enumerate(all_clusters):
        sel = keep & (best_cluster == i)
        if not sel.any():
            cluster2gene[cluster] = []
            continue
        g = genes_arr[sel]
        r = best_rank[sel]
        order = np.argsort(r, kind='stable')
        cluster2gene[cluster] = list(g[order])
    return cluster2gene


# ---------------------------------------------------------------------------
# fast group statistics for scanpy's t-test
# ---------------------------------------------------------------------------

FAST_TTEST = True          # set False to use scanpy's own _basic_stats
FAST_TTEST_CHUNK_MB = 128  # working memory for the chunked float64 accumulation


def _fast_basic_stats(self):
    '''
    Drop-in replacement for ``scanpy.tools._rank_genes_groups._RankGenes._basic_stats``.

    scanpy computes, for each of C groups, the mean and variance of the group
    AND of everything outside the group, by slicing X and calling
    ``_get_mean_var`` on each slice. That is 2C passes over the matrix -- on the
    demo data, 66 passes over 2700 x 32738, which the profiler attributes 5.66 s
    of the 6.60 s total.

    Group sums and group sums-of-squares determine every one of those numbers,
    and both come from one pass:

        S[g]  = sum over cells in g of x        -> mean_g    = S[g]  / n_g
        SQ[g] = sum over cells in g of x*x      -> meansq_g  = SQ[g] / n_g
        rest  = column total minus the group    -> mean_rest, meansq_rest

    then scanpy's own formula ``var = (meansq - mean**2) * n/(n-1)``.

    The accumulation is a chunked float64 matmul, so peak memory stays bounded
    instead of materialising ``X*X`` for every slice.

    This is a reassociation of the same sums, so results agree to float64
    rounding rather than bitwise. ``tests/test_ttest_equivalence.py`` checks
    that the marker gene lists, which is all the pipeline consumes, come out
    identical, and reports the numeric deltas.
    '''
    from scipy.sparse import issparse as _issparse

    X = self.X
    masks = self.groups_masks_obs
    n_groups, n_genes = masks.shape[0], X.shape[1]
    n_cells = X.shape[0]

    if self.ireference is not None:
        raise _FallBackToScanpy('ireference is not None')

    onehot = np.asarray(masks, dtype=np.float64).T          # (n_cells, n_groups)
    n_per_group = onehot.sum(axis=0)                        # (n_groups,)

    # The "rest" statistics below come from (column total - this group), which is
    # only the same set scanpy uses (~mask_obs) when the groups partition every
    # cell. scanpy allows `groups=['0','1']`, which covers a subset. Detect that
    # and hand the call back rather than return a different denominator.
    membership = onehot.sum(axis=1)
    if not np.array_equal(membership, np.ones(n_cells)):
        raise _FallBackToScanpy('group masks do not partition the cells '
                                '({} of {} cells belong to exactly one group)'.format(
                                    int((membership == 1).sum()), n_cells))

    S = np.zeros((n_groups, n_genes), dtype=np.float64)
    SQ = np.zeros((n_groups, n_genes), dtype=np.float64)

    if _issparse(X):
        Xc = X.tocsr()
        S = np.asarray(onehot.T @ Xc)
        # square in the STORED dtype, exactly as the dense branch and as scanpy's
        # `elem_mul(X, X)` do. Squaring float32 data in float64 here made the
        # sparse and dense paths disagree by 3.8e-06 relative, which was enough
        # to reorder tied genes and change the marker lists.
        SQ_src = _shared_index_csr(Xc, Xc.data * Xc.data)
        SQ = np.asarray(onehot.T @ SQ_src)
    else:
        # scanpy squares in the INPUT dtype (`_elem_mul_in_mem` is `X * X`) and
        # only then accumulates in float64. Squaring a float32 matrix in float64
        # here would be more accurate but would not be the same number, and the
        # difference is large enough (~1e-6 relative) to reshuffle tied genes.
        # So square in X's dtype, accumulate in float64, exactly as scanpy does.
        Xdtype = np.asarray(X[:1]).dtype
        bytes_per_row = n_genes * 8
        chunk = max(1, int(FAST_TTEST_CHUNK_MB * 1024 * 1024 / max(bytes_per_row, 1)))
        for start in range(0, n_cells, chunk):
            stop = min(start + chunk, n_cells)
            raw = np.asarray(X[start:stop])
            oh = onehot[start:stop]
            S += oh.T @ raw.astype(np.float64, copy=False)
            sq = raw * raw                      # in X's dtype, as scanpy does
            SQ += oh.T @ sq.astype(np.float64, copy=False)
            del raw, sq

    total_S = S.sum(axis=0)
    total_SQ = SQ.sum(axis=0)

    self.means = np.zeros((n_groups, n_genes))
    self.vars = np.zeros((n_groups, n_genes))
    self.means_rest = np.zeros((n_groups, n_genes))
    self.vars_rest = np.zeros((n_groups, n_genes))

    for g in range(n_groups):
        ng = n_per_group[g]
        nr = n_cells - ng
        mean_g = S[g] / ng
        var_g = (SQ[g] / ng) - mean_g ** 2
        if ng > 1:
            var_g = var_g * (ng / (ng - 1))
        self.means[g] = mean_g
        self.vars[g] = var_g

        mean_r = (total_S - S[g]) / nr
        var_r = ((total_SQ - SQ[g]) / nr) - mean_r ** 2
        if nr > 1:
            var_r = var_r * (nr / (nr - 1))
        self.means_rest[g] = mean_r
        self.vars_rest[g] = var_r

    if self.comp_pts:
        # chunked, so this never allocates a second full copy of X
        per_group_nz, total_nz = _count_nonzero_by_group(X, onehot.astype(np.float32))
        self.pts = per_group_nz / n_per_group[:, None]
        self.pts_rest = (total_nz[None, :] - per_group_nz) / \
                        (n_cells - n_per_group)[:, None]
    else:
        self.pts = None
        self.pts_rest = None


class _FallBackToScanpy(Exception):
    pass


class _patched_basic_stats:
    '''
    Swap ``_RankGenes._basic_stats`` for the chunked version while inside the block.

    Restores the original on exit, and yields control back to scanpy untouched
    if the class cannot be found (a scanpy version whose internals moved), so a
    scanpy upgrade degrades to "slower", never to "wrong".
    '''

    def __init__(self, enabled=True):
        self.enabled = enabled
        self.cls = None
        self.original = None

    def __enter__(self):
        if not self.enabled:
            return self
        try:
            from scanpy.tools._rank_genes_groups import _RankGenes
        except Exception as e:
            logger_sctriangulate.warning(
                'cannot reach scanpy _RankGenes ({}); using scanpy _basic_stats'.format(e))
            return self
        if not hasattr(_RankGenes, '_basic_stats'):
            logger_sctriangulate.warning(
                'scanpy _RankGenes has no _basic_stats; using scanpy as-is')
            return self
        self.cls = _RankGenes
        self.original = _RankGenes._basic_stats

        original = self.original

        def _dispatch(inner_self):
            try:
                return _fast_basic_stats(inner_self)
            except _FallBackToScanpy as e:
                logger_sctriangulate.info(
                    'fast t-test stats not applicable ({}); using scanpy'.format(e))
                return original(inner_self)

        _RankGenes._basic_stats = _dispatch
        return self

    def __exit__(self, *a):
        if self.cls is not None:
            self.cls._basic_stats = self.original
        return False


RUN_ENRICHMENT = False   # module default for marker_gene(run_enrichment=None)


def marker_gene(adata, key, species, criterion, folder, run_enrichment=None):
    '''
    Rank marker genes per cluster and annotate them with artifact enrichment.

    :param run_enrichment: add the ``enrichr`` and ``gsea`` columns.
        ``None`` (default) reads the module flag ``RUN_ENRICHMENT``, which is
        ``False``.

    Why off by default: neither column feeds a stability metric or a Shapley
    decision. Only ``plot_cluster_feature(feature='enrichment')``,
    ``plot_multi_modal_feature_rank``, ``penalize_artifact(mode='cellcycle')``
    and the viewer read them, and the GSEA permutations cost 3.7 s of a 21 s
    demo run. Upstream always computed them.

    Turn them back on with ``--enrichment`` on the command line, or by setting
    ``metrics.RUN_ENRICHMENT = True``. The four consumers above raise a message
    naming that flag when the columns are absent.
    '''
    if run_enrichment is None:
        run_enrichment = RUN_ENRICHMENT
    if adata.uns.get('rank_genes_groups') is not None:
        del adata.uns['rank_genes_groups']
    with _patched_basic_stats(enabled=FAST_TTEST):
        sc.tl.rank_genes_groups(adata, key, method='t-test', n_genes=adata.shape[1])
    all_genes = adata.var_names.values
    all_clusters = adata.obs[key].cat.categories
    rank_uns = adata.uns['rank_genes_groups']
    pre_computed_dfs = [compute_combo_score(rank_uns, cluster) for cluster in all_clusters]

    cluster2gene = _assign_genes_to_clusters(all_genes, all_clusters, pre_computed_dfs)

    result = pd.Series({c: cluster2gene[c] for c in all_clusters}).to_frame()
    result.columns = ['whole_marker_genes']

    col_enrichr, col_gsea, col_purify = [], [], []
    if run_enrichment:
        logger_sctriangulate.info(
            'Computing artifact enrichment for marker genes of {} (local, no network)'.format(key))
    for cluster in result.index:
        genes = result.loc[cluster, 'whole_marker_genes']
        if run_enrichment:
            col_enrichr.append(run_enrichr(genes, key=key, name=cluster, folder=folder,
                                           species=species, criterion=criterion))
            col_gsea.append(run_gsea(genes, key=key, name=cluster, folder=folder,
                                     species=species, criterion=criterion))
        col_purify.append(purify_gene(genes, species, criterion))

    if run_enrichment:
        result['enrichr'] = col_enrichr
        result['gsea'] = col_gsea
    result['purify'] = col_purify
    return result


# ---------------------------------------------------------------------------
# reassign score
# ---------------------------------------------------------------------------

# Seeds for the two sklearn estimators upstream left unseeded. See the module
# note in reassign_score and SCCAF_score, and benchmarks/determinism_audit.py.
PCA_RANDOM_STATE = 0
SCCAF_RANDOM_STATE = 0
REASSIGN_PCA = 'auto'      # 'auto' | 'sparse' | 'dense'
REASSIGN_DENSE_BUDGET_MB = 2048


def scaled_pca(X, n_components, random_state=0, mode='auto', budget_mb=None):
    '''
    PCA of the column-standardised matrix, without ever centring it in memory.

    Standardising subtracts the column mean, which turns a sparse matrix dense.
    On the marker gene pool that upstream feeds to PCA, ``n_cells`` by roughly
    ``30 * n_clusters`` genes, that is 28.8 GB at 984,119 cells and 240 clusters.

    The centring is never materialised. The standardised matrix

        Z = (X - 1 m^T) D^-1        m = column means, D = diag(column sds)

    is wrapped in a ``LinearOperator`` whose products are computed from the
    sparse X:

        Z  v = X (v / s) - (m . (v / s)) 1
        Z^T u = (X^T u) / s - m (sum u) / s

    and the top components come from ``svds`` on that operator. The result is the
    same PCA, computed from products instead of from a stored dense array.

    :param mode: ``'dense'`` reproduces sklearn's PCA exactly and is used for
        small inputs; ``'sparse'`` uses the operator; ``'auto'`` picks dense while
        the dense form fits in ``budget_mb``, so small runs keep bit-identical
        numbers and large runs stay possible.
    :return: ``(scores, mode_used)`` with scores of shape (n_cells, n_components).
    '''
    from sklearn.preprocessing import scale
    from sklearn.decomposition import PCA

    if budget_mb is None:
        budget_mb = REASSIGN_DENSE_BUDGET_MB
    n, p = X.shape
    dense_mb = n * p * 8 / 1e6

    if mode == 'auto':
        mode = 'dense' if (not issparse(X) or dense_mb <= budget_mb) else 'sparse'

    if mode == 'dense':
        # svd_solver='full', not sklearn's 'auto'. At this shape 'auto' selects the
        # RANDOMIZED solver, which does not resolve the trailing components:
        # measured against an exact full SVD on the demo marker pool, the
        # randomized 30-dimensional score space had a minimum principal-angle
        # cosine of 0.0031, so it is a different subspace, and it captured 0.8%
        # less variance. The sparse operator below scores 0.99999999 on the same
        # test. Using the exact solver here makes the dense and sparse paths agree,
        # so the answer no longer depends on which one the size triggers. It is
        # also faster on the demo: 0.43 s against 0.98 s.
        Z = scale(_as_dense(X), axis=0)
        return PCA(n_components=n_components, svd_solver='full',
                   random_state=random_state).fit_transform(Z), 'dense'

    from scipy.sparse.linalg import LinearOperator, svds
    from sklearn.utils.extmath import svd_flip

    Xc = X.tocsr() if issparse(X) else csr_matrix(X)
    m = np.asarray(Xc.mean(axis=0)).ravel()
    # sklearn's scale uses the population standard deviation (ddof=0)
    msq = np.asarray(Xc.multiply(Xc).mean(axis=0)).ravel()
    var = np.maximum(msq - m ** 2, 0.0)
    s = np.sqrt(var)
    s[s == 0] = 1.0                      # sklearn's _handle_zeros_in_scale
    inv_s = 1.0 / s

    def matvec(v):
        w = np.asarray(v).ravel() * inv_s
        return Xc.dot(w) - float(m.dot(w))

    def rmatvec(u):
        u = np.asarray(u).ravel()
        return (Xc.T.dot(u) - m * u.sum()) * inv_s

    def matmat(V):
        W = np.asarray(V) * inv_s[:, None]
        return Xc.dot(W) - np.outer(np.ones(n), m.dot(W))

    def rmatmat(U):
        U = np.asarray(U)
        return (Xc.T.dot(U) - np.outer(m, U.sum(axis=0))) * inv_s[:, None]

    Z = LinearOperator((n, p), matvec=matvec, rmatvec=rmatvec,
                       matmat=matmat, rmatmat=rmatmat, dtype=np.float64)

    k = min(int(n_components), min(n, p) - 1)
    rng = np.random.default_rng(random_state)
    v0 = rng.standard_normal(min(n, p))
    U, S, Vt = svds(Z, k=k, v0=v0)
    order = np.argsort(-S)               # svds returns ascending
    U, S, Vt = U[:, order], S[order], Vt[order]
    U, Vt = svd_flip(U, Vt)              # same sign convention sklearn applies
    scores = U * S
    if scores.shape[1] < n_components:   # pad if k had to be reduced
        scores = np.hstack([scores, np.zeros((n, n_components - scores.shape[1]))])
    return scores, 'sparse'


def reassign_score(adata, key, marker, regress_size=False):
    '''
    KNN self-projection accuracy per cluster, on a PCA of the marker gene pool.

    DETERMINISM FIX (a real change in output, measured in
    ``benchmarks/reassign_determinism.py``). Upstream had two unseeded sources
    of run-to-run variation in this function:

      * ``PCA(n_components=30)`` leaves ``random_state=None``. With 2700 cells
        and a few hundred marker genes sklearn's ``svd_solver='auto'`` selects
        the randomized solver, so repeated calls on identical input return
        different embeddings. Measured on the demo data, calling upstream's
        ``reassign_score`` ten times in one process gave different accuracies.
      * ``pool = list(set(pool))`` orders genes by string hash. Python
        randomises string hashing per interpreter, so the gene column order --
        and therefore the PCA input -- differed between processes, including
        between the worker processes upstream spawns for compute_metrics.

    Here the pool is sorted and the PCA is seeded, so the metric is
    reproducible. The values land inside upstream's own run-to-run spread; the
    benchmark script quantifies both.
    '''
    from sklearn.preprocessing import scale, LabelEncoder
    from sklearn.decomposition import PCA
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.metrics import confusion_matrix

    # get gene pool, slice the adata
    num = 30
    pool = []
    for i in range(marker.shape[0]):
        pool.extend(marker.iloc[i]['purify'][:num])
    pool = sorted(set(pool))          # upstream: list(set(pool)), hash-order dependent
    adata_now = adata[:, pool].copy()

    # standardise and reduce, without materialising the centred matrix. See
    # scaled_pca: 'auto' keeps the exact dense PCA while it fits in the budget.
    scoring, pca_mode = scaled_pca(adata_now.X, 30, random_state=PCA_RANDOM_STATE,
                                   mode=REASSIGN_PCA)
    if pca_mode == 'sparse':
        logger_sctriangulate.info(
            'reassign_score: {} cells x {} marker genes exceeds the dense budget, '
            'using the sparse PCA operator'.format(*adata_now.shape))

    le = LabelEncoder()
    scoring_y = le.fit_transform(adata_now.obs[key].astype('str'))
    order = le.classes_

    # compute the centroid of each cluster
    categories = adata_now.obs[key].cat.categories
    X = np.empty([len(categories), scoring.shape[1]])
    y = []
    for i, cluster in enumerate(categories):
        bool_index = (adata_now.obs[key] == cluster).values
        X[i, :] = np.mean(scoring[bool_index, :], axis=0)
        y.append(cluster)
    y = le.fit_transform(y)

    # train a KNN classifier on the centroids
    n_neighbors = min(10, X.shape[0])
    model = KNeighborsClassifier(n_neighbors=n_neighbors, weights='distance')
    model.fit(X, y)
    pred = model.predict(scoring)
    mat = confusion_matrix(scoring_y, pred)
    confusion_reassign = pd.DataFrame(data=mat, index=order, columns=order)
    accuracy = mat.diagonal() / mat.sum(axis=1)
    cluster_to_accuracy = {cluster: accuracy[i] for i, cluster in enumerate(order)}

    if regress_size:
        key_size_dict = get_size_in_metrics(adata.obs, key)
        df_inspect = pd.concat([pd.Series(cluster_to_accuracy), pd.Series(key_size_dict)], axis=1)
        cluster_to_accuracy = _regress_size_impl(df_inspect, regressor='GLM', to_dict=True)

    del adata_now
    return cluster_to_accuracy, confusion_reassign


# ---------------------------------------------------------------------------
# size regression
# ---------------------------------------------------------------------------

def background_normalizer(df, n_neighbors=10, scale=True):
    from copy import deepcopy
    df = deepcopy(df)
    df['order'] = np.arange(df.shape[0])
    col = []
    for i in range(df.shape[0]):
        this_metric = df[0][i]
        distance_to_this = (df[0] - this_metric).abs()
        df_tmp = deepcopy(df)
        df_tmp['distance'] = distance_to_this.values
        df_tmp.sort_values(by='distance', inplace=True)
        neighbors_metric = df_tmp.iloc[:, 0][:n_neighbors].values
        mean_ = neighbors_metric.mean()
        std_ = neighbors_metric.std()
        if scale:
            col.append(0 if std_ == 0 else (this_metric - mean_) / std_)
        else:
            col.append(this_metric - mean_)
    df['normalized'] = col
    return df


def _regress_size_impl(df_inspect, regressor='background_zscore', n_neighbors=10, to_dict=False):
    # df_inspect, index is cluster name, col1 is metric, col2 is size
    if regressor == 'background_zscore':
        df_now = background_normalizer(df_inspect, n_neighbors, True)
        df_inspect[0] = df_now['normalized'].values
        normalized_metric_series = df_inspect[0]
    elif regressor == 'background_mean':
        df_now = background_normalizer(df_inspect, n_neighbors, False)
        df_inspect[0] = df_now['normalized'].values
        normalized_metric_series = df_inspect[0]
    elif regressor == 'GLM':
        import statsmodels.api as sm
        endog = df_inspect[0]
        exog = sm.add_constant(df_inspect[1], prepend=True)
        res = sm.GLM(endog, exog, family=sm.families.Gaussian()).fit()
        normalized_metric_series = res.resid_response
    elif regressor == 'Huber':
        from sklearn.linear_model import HuberRegressor
        endog, exog = df_inspect[0], df_inspect[1]
        model = HuberRegressor().fit(exog.values.reshape(-1, 1), endog.values)
        df_inspect[0] = endog.values - model.predict(exog.values.reshape(-1, 1))
        normalized_metric_series = df_inspect[0]
    elif regressor == 'RANSAC':
        from sklearn.linear_model import RANSACRegressor
        endog, exog = df_inspect[0], df_inspect[1]
        model = RANSACRegressor().fit(exog.values.reshape(-1, 1), endog.values)
        df_inspect[0] = endog.values - model.predict(exog.values.reshape(-1, 1))
        normalized_metric_series = df_inspect[0]
    elif regressor == 'TheilSen':
        from sklearn.linear_model import TheilSenRegressor
        endog, exog = df_inspect[0], df_inspect[1]
        model = TheilSenRegressor().fit(exog.values.reshape(-1, 1), endog.values)
        df_inspect[0] = endog.values - model.predict(exog.values.reshape(-1, 1))
        normalized_metric_series = df_inspect[0]
    else:
        raise ValueError('unknown regressor: {}'.format(regressor))
    return normalized_metric_series.to_dict() if to_dict else normalized_metric_series



# Upstream named this function `regress_size` AND named the boolean parameter of
# reassign_score / SCCAF_score / tf_idf*_for_cluster `regress_size`. Inside those
# functions the parameter shadows the function, so upstream raises
# `TypeError: 'bool' object is not callable` whenever regress_size=True is passed.
# The implementation above is reachable under a private name, and the callers use
# that, so the option now works. The public name is kept for the upstream API.
regress_size = _regress_size_impl


# ---------------------------------------------------------------------------
# tf-idf
# ---------------------------------------------------------------------------

def tf_idf_bare_compute(df, cluster):
    '''
    Upstream signature kept: ``df`` is a DataFrame of expression with one extra
    'cluster' column. Retained so external callers and the reference tests keep
    working; the internal pipeline uses ``_tfidf_matrix``, which computes every
    cluster at once.
    '''
    tmp1 = df.loc[df['cluster'] == cluster, :].loc[:, df.columns != 'cluster'].values
    tf = np.count_nonzero(tmp1, axis=0) / tmp1.shape[0]
    tf = tf + 1e-5
    tmp2 = df.loc[:, df.columns != 'cluster'].values
    df_ = np.count_nonzero(tmp2, axis=0) / tmp2.shape[0]
    df_ = df_ + 1e-5
    idf = -np.log10(df_)
    return tf * idf


# Cache holds ONE entry and keeps a strong reference to the adata it was built
# from, so `is` identity is decisive and the id() can never be recycled under it.
_TFIDF_CACHE = {'adata': None, 'key': None, 'layer': None, 'matrix': None,
                'clusters': None, 'ranked': None, 'species': None, 'criterion': None}


def _tfidf_matrix(adata, key, layer=None):
    '''
    tf-idf for every cluster of ``key`` in one pass.

    Returns ``(clusters, matrix)`` where ``matrix[c]`` is the tf-idf vector over
    ``adata.var_names`` for cluster ``clusters[c]``, identical to
    ``tf_idf_bare_compute(df, clusters[c])``.

    Exactness: the nonzero counts are 0/1 sums accumulated in float64 by a
    single matmul. Sums of 0/1 below 2**53 are exact in float64, so the counts
    equal ``np.count_nonzero`` exactly, and the remaining arithmetic is the same
    expression upstream evaluates.
    '''
    c = _TFIDF_CACHE
    if c['adata'] is adata and c['key'] == key and c['layer'] == layer \
            and c['matrix'] is not None:
        return c['clusters'], c['matrix']

    mat = adata.X if layer is None else adata.layers[layer]
    n_cells = adata.shape[0]

    labels = adata.obs[key].astype('str').values
    # astype('category') is a no-op when the column already is one, and recovers
    # the sanitized (lexically sorted) categories when it is not. Upstream used
    # a bare `.cat.categories`, which raises unless scanpy's sanitize_anndata
    # happened to run first.
    clusters = list(adata.obs[key].astype('category').cat.categories)
    cluster_pos = {str(c): i for i, c in enumerate(clusters)}
    codes = np.array([cluster_pos[l] for l in labels], dtype=np.int64)

    onehot = np.zeros((n_cells, len(clusters)), dtype=np.float32)
    onehot[np.arange(n_cells), codes] = 1.0
    sizes = onehot.sum(axis=0, dtype=np.float64)       # (n_clusters,)

    per_cluster_nnz, total_nnz = _count_nonzero_by_group(mat, onehot)

    with np.errstate(divide='ignore', invalid='ignore'):
        tf = per_cluster_nnz / sizes[:, None] + 1e-5
    idf = -np.log10(total_nnz / n_cells + 1e-5)
    matrix = tf * idf[None, :]

    c.update({'adata': adata, 'key': key, 'layer': layer, 'matrix': matrix,
              'clusters': clusters, 'ranked': None,
              'species': None, 'criterion': None})
    return clusters, matrix


def _tfidf_ranked(adata, key, species, criterion, layer=None):
    '''
    Per cluster: the artifact-free tf-idf Series sorted descending, exactly as
    upstream produced it before taking ``.iloc[0]``, ``.iloc[4]`` or ``.iloc[9]``.

    Cached alongside the matrix, because upstream called tf_idf10, tf_idf5 and
    tf_idf1 separately and each redid this sort.
    '''
    clusters, matrix = _tfidf_matrix(adata, key, layer)
    c = _TFIDF_CACHE
    if c['ranked'] is not None and c['species'] == species and c['criterion'] == criterion:
        return c['ranked']

    artifact_genes = _artifact_gene_set(species, criterion)
    var_names = adata.var_names
    out = {}
    for i, cluster in enumerate(clusters):
        s = pd.Series(data=matrix[i], index=var_names).sort_values(ascending=False)
        out[cluster] = s.loc[~s.index.isin(artifact_genes)]
    c.update({'ranked': out, 'species': species, 'criterion': criterion})
    return out


def _tfidf_nth(adata, key, species, criterion, n, layer=None, want_exclusive=False):
    ranked = _tfidf_ranked(adata, key, species, criterion, layer)
    scores = {c: ranked[c].iloc[n] for c in ranked}
    if not want_exclusive:
        return scores, None
    clusters_in_order = list(ranked.keys())
    exclusive = pd.Series(
        data=[ExclusiveGeneScores(ranked[c].index, ranked[c].values) for c in clusters_in_order],
        index=clusters_in_order, name='genes', dtype=object)
    return scores, exclusive


def get_size_in_metrics(obs, key):
    counts = obs[key].value_counts()
    return {cluster: int(counts[cluster]) for cluster in pd.unique(obs[key])}


def _maybe_regress(scores, adata, key, regress_size_flag):
    if not regress_size_flag:
        return scores
    key_size_dict = get_size_in_metrics(adata.obs, key)
    df_inspect = pd.concat([pd.Series(scores), pd.Series(key_size_dict)], axis=1)
    return _regress_size_impl(df_inspect, regressor='GLM', to_dict=True)


def tf_idf10_for_cluster(adata, key, species, criterion, regress_size=False, layer=None):
    scores, exclusive = _tfidf_nth(adata, key, species, criterion, 9, layer, want_exclusive=True)
    return _maybe_regress(scores, adata, key, regress_size), exclusive


def tf_idf5_for_cluster(adata, key, species, criterion, regress_size=False, layer=None):
    scores, _ = _tfidf_nth(adata, key, species, criterion, 4, layer)
    return _maybe_regress(scores, adata, key, regress_size)


def tf_idf1_for_cluster(adata, key, species, criterion, regress_size=False, layer=None):
    scores, _ = _tfidf_nth(adata, key, species, criterion, 0, layer)
    return _maybe_regress(scores, adata, key, regress_size)


# ---------------------------------------------------------------------------
# SCCAF
# ---------------------------------------------------------------------------

SCCAF_MODE = 'optimized'   # 'optimized' or 'legacy'
SCCAF_PER_CLASS = False    # opt-in: split the one-vs-rest loop out of liblinear
SCCAF_N_JOBS = None        # per-class mode only. None -> min(n_classes, cpu_count)
SCCAF_CHUNK_MB = 128


SPARSE_SIDECAR_LAYER = '_sctri_sparse_X'


def _sparse_view(adata, col_mask):
    '''
    Sparse ``adata.X[:, col_mask]``, reusing the sidecar layer when it is there.

    ``each_key_run`` densifies X because scanpy and sklearn take different code
    paths for sparse input and would return different numbers. SCCAF wants it
    sparse again. Converting back costs a full pass over the dense matrix, which
    on the demo cancelled the time the sparse solve saved.

    ``compute_metrics`` therefore parks the original CSR in a layer before
    densifying. AnnData slices layers along with the object, so a row-subset
    view carries the matching rows and no conversion is needed. The layer costs
    what the sparse matrix costs, 27 MB on the demo against 354 MB dense.
    '''
    layer = adata.layers.get(SPARSE_SIDECAR_LAYER) if adata.layers is not None else None
    if layer is not None and layer.shape == adata.shape and issparse(layer):
        return layer.tocsr()[:, col_mask].tocsr()
    return _csr_from_columns(adata.X, col_mask)


def _csr_from_columns(mat, col_mask, chunk_mb=SCCAF_CHUNK_MB):
    '''
    CSR of ``mat[:, col_mask]`` without materialising the dense subset.

    ``each_key_run`` hands SCCAF a dense X, and ``adata[:, keep].X`` would
    allocate a full dense copy of it (337 MB on the demo). Converting one row
    block at a time keeps the peak at one block plus the finished sparse matrix,
    which is 16 MB on the demo, 20.5 times smaller than the dense copy.
    '''
    from scipy.sparse import csr_matrix, vstack

    if issparse(mat):
        return mat.tocsr()[:, col_mask].tocsr()

    n_cells, _ = mat.shape
    n_keep = int(np.count_nonzero(col_mask))
    bytes_per_row = max(n_keep, 1) * 4
    chunk = max(1, int(chunk_mb * 1024 * 1024 / bytes_per_row))
    blocks = []
    for start in range(0, n_cells, chunk):
        stop = min(start + chunk, n_cells)
        block = np.asarray(mat[start:stop])[:, col_mask]
        blocks.append(csr_matrix(block))
        del block
    return vstack(blocks, format='csr') if len(blocks) > 1 else blocks[0]


def _fit_one_class(X_train, Y_train, c):
    '''One binary L1 liblinear fit: class ``c`` against the rest.'''
    from sklearn.linear_model import LogisticRegression
    m = LogisticRegression(penalty='l1', solver='liblinear', max_iter=100000,
                           random_state=SCCAF_RANDOM_STATE)
    m.fit(X_train, (Y_train == c).astype(np.int8))
    return m.coef_[0], m.intercept_[0]


def _fit_ovr(X_train, Y_train, n_jobs=None):
    '''
    One binary L1 liblinear fit per class, over processes.

    sklearn hands the multiclass problem to liblinear, which loops the classes
    inside C on one thread. Splitting that loop and seeding each class the same
    way makes each class independent: upstream's class k consumed the random
    stream left by classes 0..k-1, so its score depended on how many classes
    came before it.

    PROCESSES, NOT THREADS. liblinear releases the GIL, and threads were the
    fastest option measured, 2.5x to 4.0x. They are also wrong: liblinear seeds
    and draws from the process-global C ``rand()``, so concurrent threads
    interleave on one shared stream and the run stops being reproducible. The
    measurement caught it -- repeated threaded runs disagreed. Processes each
    own their C state, so they reproduce.

    The values do not depend on the backend or on the core count. Verified
    bitwise on the demo: serial, 2 processes and 10 processes give identical
    coefficients for all three annotations. Only the runtime changes, so a
    result computed on a laptop matches one computed on a cluster node.

    Nested runs fall back to serial. ``compute_metrics`` already forks one
    worker per annotation, and a pool inside each of those would oversubscribe.
    '''
    from joblib import Parallel, delayed

    classes = np.unique(Y_train)
    nested = mp.current_process().name != 'MainProcess'
    if n_jobs is None:
        n_jobs = 1 if nested else min(len(classes), os.cpu_count() or 1)
    n_jobs = max(1, min(int(n_jobs), len(classes)))

    if n_jobs == 1:
        out = [_fit_one_class(X_train, Y_train, c) for c in classes]
    else:
        out = Parallel(n_jobs=n_jobs, prefer='processes')(
            delayed(_fit_one_class)(X_train, Y_train, c) for c in classes)

    coef = np.vstack([o[0] for o in out])
    intercept = np.asarray([o[1] for o in out])
    return classes, coef, intercept


def SCCAF_score(adata, key, species, criterion, scale_sccaf, regress_size=False,
                mode=None):
    '''
    Self-projection accuracy per cluster from an L1 logistic regression.

    Split the cells in half, fit an L1 logistic regression on one half to
    predict cluster identity from expression, score it on the other half. A
    cluster a classifier can recover is a cluster the data supports.

    :param mode: ``'optimized'`` (default, see ``SCCAF_MODE``) or ``'legacy'``.

    ``legacy`` is upstream's procedure: a dense matrix handed to sklearn, which
    lets liblinear loop the classes internally on one thread.

    ``optimized`` computes the same quantity with three changes:

      1. the matrix stays sparse. Measured on the demo, 337 MB becomes 16 MB and
         the accuracies are bitwise identical to the dense fit, because the
         solver sees the same numbers.
      2. the full matrix is freed before the fit, and never densified.
      3. ``SCCAF_PER_CLASS=True`` (off by default) additionally splits
         liblinear's one-vs-rest loop into separately seeded fits that run over
         processes. That decouples each class from the ones solved before it,
         but it moves accuracies by up to 0.167 on the demo, so it is opt-in.

    With the default settings the optimized mode returns the LEGACY NUMBERS.
    liblinear keeps its own one-vs-rest loop and sees the same values, just in a
    sparse container. Measured on the demo: maximum per-cluster difference
    0.0000 across all 61 clusters, and ARI 1.0000 on the final labels.
    ``benchmarks/sccaf_compare.py`` produces those figures.

    DETERMINISM. Upstream left ``random_state=None``, so it returned a different
    accuracy on every call: 1 of 50 repeats matched the first, spread up to
    0.167. Both modes here seed it with ``SCCAF_RANDOM_STATE``.
    '''
    if mode is None:
        mode = SCCAF_MODE
    if mode not in ('optimized', 'legacy'):
        raise ValueError("SCCAF mode must be 'optimized' or 'legacy', got {!r}".format(mode))
    if mode == 'legacy':
        return _sccaf_legacy(adata, key, species, criterion, scale_sccaf, regress_size)
    return _sccaf_optimized(adata, key, species, criterion, scale_sccaf, regress_size)


def _sccaf_optimized(adata, key, species, criterion, scale_sccaf, regress_size=False):
    from sklearn.preprocessing import LabelEncoder
    from sklearn.model_selection import StratifiedShuffleSplit
    from sklearn.metrics import confusion_matrix

    if scale_sccaf:
        # mean-centering destroys sparsity, so there is nothing to gain here;
        # fall back to the dense procedure rather than densify silently
        logger_sctriangulate.info(
            'SCCAF: scale_sccaf=True needs a dense matrix, using the legacy path for {}'.format(key))
        return _sccaf_legacy(adata, key, species, criterion, scale_sccaf, regress_size)

    artifact_genes = _artifact_gene_set(species, criterion)
    keep = np.asarray(~adata.var_names.isin(artifact_genes))
    X = _sparse_view(adata, keep)
    Y = adata.obs[key].values

    le = LabelEncoder()
    Y = le.fit_transform(Y)
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.5, random_state=0)
    train_index, test_index = next(iter(sss.split(X, Y)))
    X_train, Y_train = X[train_index], Y[train_index]
    X_test, Y_test = X[test_index], Y[test_index]
    del X

    if SCCAF_PER_CLASS:
        # opt-in. Splits liblinear's internal one-vs-rest loop into separate
        # seeded fits, which parallelises over processes and decouples each
        # class from the ones before it. It moves accuracies by up to 0.167 on
        # the demo, against a legacy self-spread of up to 0.167, and it pays a
        # process-pool startup that only earns its cost with many classes.
        classes, coef, intercept = _fit_ovr(X_train, Y_train, SCCAF_N_JOBS)
        result = classes[np.asarray(X_test @ coef.T + intercept).argmax(axis=1)]
    else:
        # default. liblinear keeps its own one-vs-rest loop, so the numbers are
        # exactly the legacy numbers; only the matrix representation changed.
        from sklearn.linear_model import LogisticRegression
        model = LogisticRegression(penalty='l1', solver='liblinear', max_iter=100000,
                                   random_state=SCCAF_RANDOM_STATE)
        model.fit(X_train, Y_train)
        result = model.predict(X_test)

    m = confusion_matrix(Y_test, result)
    confusion_sccaf = pd.DataFrame(data=m, index=le.classes_, columns=le.classes_)
    reliability = m.diagonal() / m.sum(axis=1)
    cluster_to_SCCAF = {le.classes_[i]: reliability[i] for i in range(len(reliability))}

    if regress_size:
        key_size_dict = get_size_in_metrics(adata.obs, key)
        df_inspect = pd.concat([pd.Series(cluster_to_SCCAF), pd.Series(key_size_dict)], axis=1)
        cluster_to_SCCAF = _regress_size_impl(df_inspect, regressor='GLM', to_dict=True)

    del X_train, X_test
    return cluster_to_SCCAF, confusion_sccaf


def _sccaf_legacy(adata, key, species, criterion, scale_sccaf, regress_size=False):
    '''Upstream's procedure, with the random seed pinned. Kept as an alternative.'''
    from sklearn.preprocessing import LabelEncoder, scale
    from sklearn.model_selection import StratifiedShuffleSplit
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import confusion_matrix

    artifact_genes = _artifact_gene_set(species, criterion)
    keep = ~adata.var_names.isin(artifact_genes)
    X = _as_dense(adata[:, keep].X)
    Y = adata.obs[key].values

    if scale_sccaf:
        X = scale(X, axis=0)

    le = LabelEncoder()
    Y = le.fit_transform(Y)
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.5, random_state=0)
    train_index, test_index = next(iter(sss.split(X, Y)))
    X_train, Y_train = X[train_index], Y[train_index]
    X_test, Y_test = X[test_index], Y[test_index]
    # Free the full matrix before the fit. The two halves are independent copies
    # made by fancy indexing, and liblinear allocates a float64 copy of X_train
    # on top. Holding X across the fit added its full size to the peak: 608 MB
    # on an 8,000 by 19,000 float32 matrix, about a third of this step's peak.
    del X

    model = LogisticRegression(penalty='l1', solver='liblinear', max_iter=100000,
                               random_state=SCCAF_RANDOM_STATE)   # upstream left this unseeded
    model.fit(X_train, Y_train)
    result = model.predict(X_test)
    m = confusion_matrix(Y_test, result)
    confusion_sccaf = pd.DataFrame(data=m, index=le.classes_, columns=le.classes_)
    reliability = m.diagonal() / m.sum(axis=1)
    cluster_to_SCCAF = {le.classes_[i]: reliability[i] for i in range(len(reliability))}

    if regress_size:
        key_size_dict = get_size_in_metrics(adata.obs, key)
        df_inspect = pd.concat([pd.Series(cluster_to_SCCAF), pd.Series(key_size_dict)], axis=1)
        cluster_to_SCCAF = _regress_size_impl(df_inspect, regressor='GLM', to_dict=True)

    del X_train, X_test
    return cluster_to_SCCAF, confusion_sccaf


def single_size_query(obs, c):
    key = list(c.keys())[0]
    cluster = list(c.values())[0]
    return obs.loc[obs[key] == cluster, :].shape[0]
