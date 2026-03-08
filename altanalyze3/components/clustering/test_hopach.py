import numpy as np
import scipy.sparse as sp

from altanalyze3.components.clustering.hopach import distancematrix, hopach


def _toy_data():
    cluster_a = np.array([[0.0, 0.0], [0.2, 0.1], [0.1, 0.2]])
    cluster_b = np.array([[5.0, 5.0], [5.1, 5.2], [5.2, 5.1]])
    return np.vstack([cluster_a, cluster_b])


def test_hopach_two_cluster_dense():
    data = _toy_data()
    result = hopach(
        data,
        d="euclid",
        kmax=4,
        kmin=2,
        mincluster=3,
        random_state=7,
        ordering="medoid",
    )
    sizes = sorted(result.clust.sizes.tolist())
    assert result.clust.k == 2
    assert sizes == [3, 3]
    assert len(result.order) == data.shape[0]


def test_hopach_sparse_equivalence():
    data = _toy_data()
    sparse_data = sp.csr_matrix(data)
    result_dense = hopach(
        data,
        d="euclid",
        kmax=4,
        kmin=2,
        mincluster=3,
        random_state=11,
        ordering="medoid",
    )
    result_sparse = hopach(
        sparse_data,
        d="euclid",
        kmax=4,
        kmin=2,
        mincluster=3,
        random_state=11,
        ordering="medoid",
    )
    assert sorted(result_dense.clust.sizes.tolist()) == sorted(
        result_sparse.clust.sizes.tolist()
    )
    assert len(result_sparse.order) == data.shape[0]


def test_hopach_with_precomputed_distance():
    data = _toy_data()
    dmat = distancematrix(data, d="euclid")
    result = hopach(
        data,
        dmat=dmat,
        kmax=4,
        kmin=2,
        mincluster=3,
        random_state=5,
    )
    assert result.clust.k == 2
