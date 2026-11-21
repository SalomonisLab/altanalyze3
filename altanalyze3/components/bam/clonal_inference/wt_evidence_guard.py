"""Opt-in reporting sensitivity test; does not refit or change clone genotypes.

Cells with WT evidence and no accepted mutation cannot be assigned a malignant
clone by its prevalence prior. Missing loci remain unknown in the input matrix.
LOH-negative-only evidence is retained separately, not called variant-WT.
This module is deliberately not the workflow default.
"""
import numpy as np
from .integrity import require


def protect_wt_evidence(X, posterior, non_loh_features):
    """Return memberships [baseline columns, WT evidence, LOH-only, no calls].

Baseline memberships are unchanged for every cell with an accepted mutation.
The three new categories describe observed evidence, not complete genotypes.
"""
    X = np.asarray(X)
    P = np.asarray(posterior, dtype=float)
    variant = np.asarray(non_loh_features, dtype=bool)
    require(X.ndim == 2 and np.isin(X, [-1, 0, 1]).all(), 'Invalid calls')
    require(variant.shape == (X.shape[1],), 'Feature annotation mismatch')
    require(P.ndim == 2 and P.shape[0] == X.shape[0] and P.shape[1] > 0,
            'Posterior shape mismatch')
    require(np.isfinite(P).all() and (P >= 0).all() and
            np.allclose(P.sum(1), 1), 'Invalid posterior')
    mutant = (X == 1).any(1)
    wt_only = ~mutant & (X == 0).any(1)
    variant_wt = wt_only & (X[:, variant] == 0).any(1)
    loh_only = wt_only & ~variant_wt
    no_call = (X == -1).all(1)
    Q = np.zeros((len(X), P.shape[1]+3))
    Q[mutant, :P.shape[1]] = P[mutant]
    Q[variant_wt, -3] = 1
    Q[loh_only, -2] = 1
    Q[no_call, -1] = 1
    require(np.allclose(Q.sum(1), 1), 'Cell mass not conserved')
    require(np.array_equal(Q[mutant, :P.shape[1]], P[mutant]),
            'Changed mutation-positive memberships')
    require(not Q[wt_only, :P.shape[1]].any(), 'WT-only cell retained clone mass')
    return Q, np.select([variant_wt, loh_only, no_call],
                        ['WT_variant_evidence', 'LOH_negative_only', 'no_call'],
                        default='accepted_mutation')
