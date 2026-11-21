"""Experimental cross-check on a fixed tree, never a replacement genotype fit.

Report a least-specific active lineage only if it contains all accepted MUT calls
and is a subset of every other compatible active genotype. This is compatibility,
not proof that descendant mutations are absent. Missing is neutral; reference RNA
is flagged, not interpreted as proven DNA absence. No MUT => operational WT.
"""
import numpy as np


def supported_lineage(X, G, active):
    X, G, active = np.asarray(X), np.asarray(G), np.asarray(active)
    if X.ndim != 2 or not np.isin(X, [-1, 0, 1]).all():
        raise ValueError('Invalid accepted-call matrix')
    if G.ndim != 2 or G.shape[1] != X.shape[1] or not np.isin(G, [0, 1]).all():
        raise ValueError('Invalid genotype matrix')
    if (active.ndim != 1 or not np.issubdtype(active.dtype, np.integer)
            or len(set(active)) != len(active) or 0 not in active
            or (active < 0).any() or (active >= len(G)).any() or G[0].any()):
        raise ValueError('Invalid active nodes or WT genotype')
    positive = X == 1
    has = positive.any(1)
    nodes = active[active != 0]
    omitted = positive.astype(np.int32) @ (1-G[nodes]).astype(np.int32).T
    compatible = omitted == 0
    # A smallest mutation count alone is insufficient: candidates may be on
    # incomparable branches. The selected genotype must be shared by all of them.
    subset = ((G[nodes].astype(np.int32) @ (1-G[nodes]).astype(np.int32).T) == 0)
    best = compatible & ((compatible.astype(np.int32) @ (~subset).astype(np.int32).T) == 0)
    count = best.sum(1)
    assigned = np.full(len(X), -1, int)
    unique = has & (count == 1)
    if len(nodes):
        assigned[unique] = nodes[best[unique].argmax(1)]
    assigned[~has] = 0
    reference_conflicts = np.full(len(X), -1, int)
    for node in active:
        ix = assigned == node
        reference_conflicts[ix] = ((X[ix] == 0) & (G[node] == 1)).sum(1)
        if (positive[ix] & (G[node] == 0)).any():
            raise ValueError('A reported lineage omits accepted mutations')
    return dict(node=assigned, compatible=compatible, candidate_nodes=nodes,
        ambiguous_minimum=has & compatible.any(1) & (count != 1),
        incompatible=has & ~compatible.any(1),
        negative_evidence_conflicts=reference_conflicts)
