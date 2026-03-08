"""Fast, sparse-friendly Python implementation of the HOPACH algorithm.

This module ports the core ideas from the R hopach package:
hierarchical partitioning via k-medoids, median split silhouette
selection, optional collapsing, and ordered output.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering
from scipy.spatial.distance import cdist, squareform

try:  # Optional acceleration for pairwise distances with sparse inputs.
    from sklearn.metrics import pairwise_distances

    _HAVE_SKLEARN = True
except Exception:  # pragma: no cover - fallback path
    _HAVE_SKLEARN = False


_DISTANCE_ALIASES = {
    "euclid": "euclidean",
    "euclidean": "euclidean",
    "l2": "euclidean",
    "manhattan": "manhattan",
    "cityblock": "manhattan",
    "l1": "manhattan",
    "cor": "correlation",
    "correlation": "correlation",
    "abscor": "abscorrelation",
    "cosangle": "cosine",
    "cosine": "cosine",
    "abscosangle": "abscosine",
}


@dataclass
class HopachNode:
    node_id: int
    indices: np.ndarray
    depth: int
    medoid: int
    children: List[int] = field(default_factory=list)
    k: int = 1
    mss: float = 0.0


@dataclass
class HopachClust:
    labels: np.ndarray
    sizes: np.ndarray
    k: int
    mss: float
    medoids: np.ndarray


@dataclass
class HopachResult:
    clust: HopachClust
    order: np.ndarray
    tree: List[HopachNode]
    levels: Dict[int, np.ndarray]
    dmat: Optional[np.ndarray] = None


def distancematrix(
    data: np.ndarray | sp.spmatrix,
    d: str = "euclid",
    dmat: Optional[np.ndarray] = None,
    dmethod: Optional[str] = None,
    dtype: np.dtype = np.float32,
    n_jobs: Optional[int] = None,
) -> np.ndarray:
    """Compute pairwise distances supporting hopach metrics with sparse inputs."""
    if dmat is not None:
        return np.asarray(dmat, dtype=dtype)

    if dmethod is not None:
        d = dmethod
    metric = _DISTANCE_ALIASES.get(d, d).lower()
    is_sparse = sp.issparse(data)
    if is_sparse:
        data = data.tocsr()
    else:
        data = np.asarray(data)

    if metric == "euclidean":
        dist = _pairwise_euclidean(data, dtype=dtype)
    elif metric == "manhattan":
        dist = _pairwise_manhattan(data, dtype=dtype, n_jobs=n_jobs)
    elif metric in {"cosine", "abscosine"}:
        dist = _pairwise_cosine(data, absolute=(metric == "abscosine"), dtype=dtype)
    elif metric in {"correlation", "abscorrelation"}:
        dist = _pairwise_correlation(
            data, absolute=(metric == "abscorrelation"), dtype=dtype
        )
    else:
        # Fallback to scipy for custom metrics on dense arrays.
        if is_sparse:
            data = data.toarray()
        dist = cdist(data, data, metric=metric).astype(dtype, copy=False)

    np.fill_diagonal(dist, 0.0)
    return dist


def labelstomss(labels: Sequence[int], dmat: np.ndarray) -> Dict[str, np.ndarray | float]:
    """Compute median split silhouette (MSS) and per-sample silhouettes."""
    labels = np.asarray(labels)
    sil = _silhouette_samples(dmat, labels)
    return {"mss": float(np.median(sil)), "silhouette": sil}


def hopach(
    data: np.ndarray | sp.spmatrix,
    d: str = "euclid",
    dmat: Optional[np.ndarray] = None,
    dmethod: Optional[str] = None,
    kmax: int = 9,
    kmin: int = 2,
    mincluster: int = 5,
    maxlevel: Optional[int] = None,
    mss_threshold: float = 0.0,
    collapse: bool = True,
    ordering: str = "medoid",
    order_linkage: str = "average",
    max_iter: int = 100,
    init: str = "kmedoids++",
    random_state: Optional[int] = 0,
    seed: Optional[int] = None,
    return_dmat: bool = False,
    n_jobs: Optional[int] = None,
) -> HopachResult:
    """Run HOPACH clustering on dense or sparse inputs.

    Parameters follow the R hopach defaults where possible:
    - d: distance type ("euclid", "manhattan", "cor", "abscor", "cosangle", "abscosangle")
    - dmat: optional precomputed distance matrix
    - kmax/kmin: search range for k-medoids splits
    - mincluster: minimum cluster size to allow further splitting
    - maxlevel: maximum tree depth
    - mss_threshold: minimum MSS required to accept a split
    - collapse: if True, reject splits that do not improve MSS
    - ordering: "medoid", "hclust", or "none"
    """
    if seed is not None:
        random_state = seed
    rng = np.random.default_rng(random_state)
    dmat = distancematrix(data, d=d, dmat=dmat, dmethod=dmethod, n_jobs=n_jobs)
    n = dmat.shape[0]
    indices = np.arange(n)

    nodes: List[HopachNode] = []
    levels: Dict[int, np.ndarray] = {}

    def build_node(node_indices: np.ndarray, depth: int) -> int:
        node_id = len(nodes)
        medoid = _medoid_index(dmat, node_indices)
        node = HopachNode(
            node_id=node_id,
            indices=node_indices,
            depth=depth,
            medoid=medoid,
        )
        nodes.append(node)

        if maxlevel is not None and depth >= maxlevel:
            return node_id
        if node_indices.size < max(mincluster, kmin + 1):
            return node_id

        kmax_local = min(kmax, node_indices.size - 1)
        if kmax_local < kmin:
            return node_id

        sub_dmat = dmat[np.ix_(node_indices, node_indices)]
        best = _best_kmedoids_split(
            sub_dmat,
            kmin=kmin,
            kmax=kmax_local,
            rng=rng,
            max_iter=max_iter,
            init=init,
        )
        if best is None:
            return node_id

        labels, _, mss = best
        if mss < mss_threshold:
            return node_id
        if collapse and mss <= 0:
            return node_id

        node.k = int(np.max(labels)) + 1
        node.mss = float(mss)

        child_ids: List[int] = []
        for k in range(node.k):
            child_mask = labels == k
            if not np.any(child_mask):
                continue
            child_indices = node_indices[child_mask]
            child_id = build_node(child_indices, depth + 1)
            child_ids.append(child_id)

        if child_ids:
            node.children = child_ids
        return node_id

    root_id = build_node(indices, 0)

    leaf_nodes = [node for node in nodes if not node.children]
    labels = np.full(n, -1, dtype=int)
    medoids = np.zeros(len(leaf_nodes), dtype=int)
    sizes = np.zeros(len(leaf_nodes), dtype=int)
    for idx, node in enumerate(leaf_nodes):
        labels[node.indices] = idx + 1  # 1-based to mirror R output
        medoids[idx] = node.medoid
        sizes[idx] = node.indices.size

    final_mss = labelstomss(labels, dmat)["mss"]
    order = _order_tree(nodes, root_id, dmat, ordering, order_linkage)

    clust = HopachClust(labels=labels, sizes=sizes, k=len(leaf_nodes), mss=final_mss, medoids=medoids)
    levels[0] = labels.copy()

    result = HopachResult(
        clust=clust,
        order=order,
        tree=nodes,
        levels=levels,
        dmat=dmat if return_dmat else None,
    )
    return result


def _pairwise_euclidean(
    data: np.ndarray | sp.spmatrix, dtype: np.dtype
) -> np.ndarray:
    if sp.issparse(data):
        norms = np.asarray(data.power(2).sum(axis=1)).ravel()
        prod = data @ data.T
        dist2 = norms[:, None] + norms[None, :] - 2.0 * prod.toarray()
    else:
        norms = np.sum(data ** 2, axis=1)
        dist2 = norms[:, None] + norms[None, :] - 2.0 * np.dot(data, data.T)
    dist2[dist2 < 0] = 0.0
    dist = np.sqrt(dist2, out=dist2)
    return dist.astype(dtype, copy=False)


def _pairwise_manhattan(
    data: np.ndarray | sp.spmatrix, dtype: np.dtype, n_jobs: Optional[int]
) -> np.ndarray:
    if _HAVE_SKLEARN:
        return pairwise_distances(data, metric="manhattan", n_jobs=n_jobs).astype(dtype)
    if sp.issparse(data):
        data = data.toarray()
    return cdist(data, data, metric="cityblock").astype(dtype)


def _pairwise_cosine(
    data: np.ndarray | sp.spmatrix, absolute: bool, dtype: np.dtype
) -> np.ndarray:
    if sp.issparse(data):
        norms = np.sqrt(np.asarray(data.power(2).sum(axis=1)).ravel())
        norms[norms == 0] = 1.0
        inv_norms = 1.0 / norms
        data = data.multiply(inv_norms[:, None])
        sim = (data @ data.T).toarray()
    else:
        norms = np.linalg.norm(data, axis=1)
        norms[norms == 0] = 1.0
        sim = np.dot(data / norms[:, None], (data / norms[:, None]).T)
    if absolute:
        sim = np.abs(sim)
    sim = np.clip(sim, -1.0, 1.0)
    dist = np.arccos(sim) / np.pi
    return dist.astype(dtype, copy=False)


def _pairwise_correlation(
    data: np.ndarray | sp.spmatrix, absolute: bool, dtype: np.dtype
) -> np.ndarray:
    if sp.issparse(data):
        data = data.toarray()
    data = np.asarray(data)
    means = data.mean(axis=1, keepdims=True)
    centered = data - means
    norms = np.linalg.norm(centered, axis=1)
    norms[norms == 0] = 1.0
    sim = np.dot(centered / norms[:, None], (centered / norms[:, None]).T)
    if absolute:
        sim = np.abs(sim)
    dist = 1.0 - sim
    dist[dist < 0] = 0.0
    return dist.astype(dtype, copy=False)


def _best_kmedoids_split(
    dmat: np.ndarray,
    kmin: int,
    kmax: int,
    rng: np.random.Generator,
    max_iter: int,
    init: str,
) -> Optional[Tuple[np.ndarray, np.ndarray, float]]:
    best_labels = None
    best_medoids = None
    best_mss = -np.inf
    for k in range(kmin, kmax + 1):
        labels, medoids = _kmedoids(dmat, k, rng=rng, max_iter=max_iter, init=init)
        sil = _silhouette_samples(dmat, labels)
        mss = float(np.median(sil))
        if mss > best_mss:
            best_mss = mss
            best_labels = labels
            best_medoids = medoids
    if best_labels is None:
        return None
    return best_labels, best_medoids, best_mss


def _kmedoids(
    dmat: np.ndarray,
    k: int,
    rng: np.random.Generator,
    max_iter: int,
    init: str,
) -> Tuple[np.ndarray, np.ndarray]:
    n = dmat.shape[0]
    medoids = _init_medoids(dmat, k, rng, init)
    labels = np.argmin(dmat[:, medoids], axis=1)
    labels = _resolve_empty_clusters(dmat, labels, medoids, rng)

    for _ in range(max_iter):
        new_medoids = medoids.copy()
        for cluster_id in range(k):
            idx = np.where(labels == cluster_id)[0]
            if idx.size == 0:
                continue
            sub_d = dmat[np.ix_(idx, idx)]
            new_medoids[cluster_id] = idx[np.argmin(sub_d.sum(axis=1))]
        if np.array_equal(new_medoids, medoids):
            break
        medoids = new_medoids
        labels = np.argmin(dmat[:, medoids], axis=1)
        labels = _resolve_empty_clusters(dmat, labels, medoids, rng)
    return labels, medoids


def _init_medoids(
    dmat: np.ndarray, k: int, rng: np.random.Generator, init: str
) -> np.ndarray:
    n = dmat.shape[0]
    if k >= n:
        return np.arange(n)
    if init == "random":
        return rng.choice(n, size=k, replace=False)
    # k-medoids++ style seeding
    medoids = [int(rng.integers(0, n))]
    for _ in range(1, k):
        dist_to_nearest = np.min(dmat[:, medoids], axis=1)
        probs = dist_to_nearest ** 2
        probs_sum = probs.sum()
        if probs_sum == 0:
            remaining = np.setdiff1d(np.arange(n), medoids)
            medoids.append(int(rng.choice(remaining)))
        else:
            probs /= probs_sum
            medoids.append(int(rng.choice(n, p=probs)))
    return np.asarray(medoids, dtype=int)


def _resolve_empty_clusters(
    dmat: np.ndarray, labels: np.ndarray, medoids: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    k = medoids.size
    counts = np.bincount(labels, minlength=k)
    empty = np.where(counts == 0)[0]
    if empty.size == 0:
        return labels
    for cluster_id in empty:
        # Reassign the farthest point from current medoids to the empty cluster.
        dist_to_medoids = np.min(dmat[:, medoids], axis=1)
        farthest = int(np.argmax(dist_to_medoids))
        labels[farthest] = cluster_id
    return labels


def _silhouette_samples(dmat: np.ndarray, labels: np.ndarray) -> np.ndarray:
    labels = np.asarray(labels)
    n = labels.size
    sil = np.zeros(n, dtype=float)
    unique = np.unique(labels)
    clusters = {lab: np.where(labels == lab)[0] for lab in unique}

    for lab, idx in clusters.items():
        if idx.size == 1:
            sil[idx] = 0.0
            continue
        intra = dmat[np.ix_(idx, idx)]
        a = (intra.sum(axis=1) - np.diag(intra)) / (idx.size - 1)

        b = np.full(idx.size, np.inf)
        for other_lab, jdx in clusters.items():
            if other_lab == lab:
                continue
            inter = dmat[np.ix_(idx, jdx)].mean(axis=1)
            b = np.minimum(b, inter)

        denom = np.maximum(a, b)
        denom[denom == 0] = 1.0
        sil[idx] = (b - a) / denom
    sil[np.isnan(sil)] = 0.0
    return sil


def _medoid_index(dmat: np.ndarray, indices: np.ndarray) -> int:
    sub = dmat[np.ix_(indices, indices)]
    return int(indices[np.argmin(sub.sum(axis=1))])


def _order_tree(
    nodes: List[HopachNode],
    node_id: int,
    dmat: np.ndarray,
    ordering: str,
    linkage_method: str,
) -> np.ndarray:
    node = nodes[node_id]
    if not node.children:
        if ordering == "none":
            return node.indices
        if ordering == "hclust" and node.indices.size > 2:
            sub = dmat[np.ix_(node.indices, node.indices)]
            order = _seriation_order(sub, linkage_method)
            return node.indices[order]
        # Default: order by distance to medoid.
        return node.indices[np.argsort(dmat[node.medoid, node.indices])]

    child_ids = node.children
    child_medoids = np.array([nodes[c].medoid for c in child_ids], dtype=int)
    if child_medoids.size > 1:
        sub = dmat[np.ix_(child_medoids, child_medoids)]
        child_order = _seriation_order(sub, linkage_method)
        ordered_children = [child_ids[i] for i in child_order]
    else:
        ordered_children = child_ids

    ordered_indices: List[int] = []
    for child_id in ordered_children:
        ordered_indices.extend(_order_tree(nodes, child_id, dmat, ordering, linkage_method))
    return np.asarray(ordered_indices, dtype=int)


def _seriation_order(dmat: np.ndarray, linkage_method: str) -> np.ndarray:
    if dmat.shape[0] <= 2:
        return np.arange(dmat.shape[0])
    condensed = squareform(dmat, checks=False)
    z = linkage(condensed, method=linkage_method)
    try:
        z = optimal_leaf_ordering(z, condensed)
    except Exception:
        pass
    return leaves_list(z)
