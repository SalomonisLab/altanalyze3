"""Structural checks for a proposed gain-only clonal tree; not a tree estimator.

Genotypes are latent complete binary states, NOT observed RNA call vectors.
The checker never completes an unknown observation or alters a fitted tree.
Passing these checks establishes internal consistency, not biological accuracy.
"""
import numpy as np
from .integrity import require


def validate_gain_tree(genotypes, parents, clone_ids, features, founder=None,
                       malignant_clones=None):
    G = np.asarray(genotypes)
    parents = np.asarray(parents)
    require(G.ndim == 2, 'Genotypes must be a matrix')
    n, m = G.shape
    require(n > 0 and len(clone_ids) == n and len(features) == m, 'Identity dimensions differ')
    require(len(set(clone_ids)) == n and len(set(features)) == m, 'Duplicate identity')
    require(np.isin(G, [0, 1]).all(), 'Latent genotypes must be complete binary states; observations may remain unknown')
    require(parents.shape == (n,) and np.issubdtype(parents.dtype, np.integer), 'Invalid parent indices')
    roots = np.flatnonzero(parents == -1)
    require(len(roots) == 1, 'Tree needs exactly one root')
    root = int(roots[0])
    require(((parents >= -1) & (parents < n)).all(), 'Parent out of range')
    require((G[root] == 0).all(), 'Explicit normal root must be all zero')
    gains = G.copy()
    for child in range(n):
        seen = set()
        node = child
        while node != -1:
            require(node not in seen, 'Ancestry cycle')
            seen.add(node)
            node = int(parents[node])
        require(root in seen, 'Disconnected clone')
        parent = int(parents[child])
        if parent >= 0:
            require((G[child] >= G[parent]).all(), 'Descendant loses an ancestral mutation')
            gains[child] = G[child] - G[parent]
            require(gains[child].any(), 'Indistinguishable zero-gain nodes must not be separate genetic clones')
    require((gains.sum(axis=0) <= 1).all(), 'Mutation acquired independently more than once')
    if founder is not None:
        require(bool(founder) and set(founder) <= set(features), 'Unresolved founder feature')
        require(malignant_clones is not None, 'Explicit malignant clone identities required')
        require(set(malignant_clones) <= set(clone_ids), 'Unresolved malignant clone')
        columns = [list(features).index(f) for f in founder]
        for clone in malignant_clones:
            require(G[list(clone_ids).index(clone), columns].all(), 'Malignant descendant lacks required founder')
    return dict(root=root, acquired_events=gains)


def lineage_frequencies(exclusive_counts, parents, all_cells):
    """Counts include explicit WT and unresolved categories in the denominator.

    exclusive_counts: samples x tree nodes. Unresolved cells need not be attached.
    The exclusive founder population can be smaller than a daughter population.
    The cumulative founder lineage includes every descendant and obeys the sum rule.
    """
    counts = np.asarray(exclusive_counts, dtype=float)
    parents = np.asarray(parents)
    total = np.asarray(all_cells, dtype=float)
    require(counts.ndim == 2 and counts.shape[1] == len(parents), 'Count dimensions differ')
    require(total.shape == (counts.shape[0],), 'Sample denominators differ')
    require(np.isfinite(counts).all() and np.isfinite(total).all(), 'Nonfinite count')
    require((counts >= 0).all() and (total > 0).all(), 'Negative counts or empty specimen')
    require((counts.sum(axis=1) <= total + 1e-8).all(), 'Exclusive assignments exceed all cells')
    require(np.issubdtype(parents.dtype,np.integer) and ((parents >= -1)&(parents<len(parents))).all(), 'Invalid parent indices')
    require((parents == -1).sum() == 1, 'Tree needs exactly one root')
    cumulative = counts.copy()
    for child in range(len(parents)):
        node = int(parents[child]); seen = {child}
        while node != -1:
            require(node not in seen, 'Ancestry cycle')
            seen.add(node); cumulative[:, node] += counts[:, child]; node = int(parents[node])
    for parent in range(len(parents)):
        children = np.flatnonzero(parents == parent)
        require(np.allclose(cumulative[:, parent], counts[:, parent]+cumulative[:,children].sum(axis=1)), 'Lineage sum rule failed')
    return dict(exclusive_percent=100*counts/total[:,None],
                lineage_percent=100*cumulative/total[:,None],
                unresolved_percent=100*(total-counts.sum(axis=1))/total)
