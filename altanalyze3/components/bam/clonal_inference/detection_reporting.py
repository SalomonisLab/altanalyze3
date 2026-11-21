"""Default detection-tree reporting approved on 2026-11-16.

WT means no accepted direct or imputed mutant call, including missing-only
cells. Original calls and native fitted posteriors are never overwritten.
"""
import numpy as np
import pandas as pd
from .integrity import require

POLICY = 'no-mutant-call-WT-v1'
DESCRIPTION = ('WT = no accepted direct or imputed MUT call, including no-call '
               'and LOH-negative-only cells; mutation-positive memberships are '
               'conditioned on the original malignant clones; all-cell denominator')


def memberships(X, posterior, active, G):
    X, P, active, G = map(np.asarray, (X, posterior, active, G))
    require(X.ndim == 2 and np.isin(X, [-1, 0, 1]).all(), 'Invalid calls')
    require(active.ndim == 1 and len(set(active)) == len(active) and
            np.issubdtype(active.dtype, np.integer) and (active >= 0).all() and
            (active < len(G)).all() and list(active).count(0) == 1,
            'Invalid active node identities')
    require(G.ndim == 2 and G.shape[1] == X.shape[1] and
            np.isin(G, [0, 1]).all() and not G[0].any(), 'Invalid genotypes')
    require(P.shape == (len(X), len(active)) and np.isfinite(P).all() and
            (P >= 0).all() and np.allclose(P.sum(1), 1), 'Invalid posterior')
    malignant = active != 0
    require(G[active[malignant]].any(1).all(), 'Non-WT component lacks events')
    positive = (X == 1).any(1)
    mass = P[:, malignant].sum(1)
    require((mass[positive] > 0).all(), 'Mutation-positive cell lacks malignant support')
    Q = np.zeros_like(P, dtype=float)
    Q[~positive, np.flatnonzero(active == 0)[0]] = 1
    Q[np.ix_(positive, malignant)] = P[np.ix_(positive, malignant)] / mass[positive, None]
    require(np.allclose(Q.sum(1), 1), 'Reporting membership mass loss')
    return Q


def tables(X, cells, events, model):
    """One policy shared by future runs, saved-fit reporting and all figures."""
    require(list(cells.cell_id.astype(str)) == list(model['cell_ids'].astype(str)),
            'Reporting cell identity/order mismatch')
    require(list(events.feature.astype(str)) == list(model['features'].astype(str)),
            'Reporting event identity/order mismatch')
    require(cells.cell_id.is_unique and events.feature.is_unique, 'Duplicate identity')
    G, active = model['G'], model['active']
    Q = memberships(X, model['posterior'], active, G)
    assigned = active[Q.argmax(1)]
    A = cells.copy()
    A['node'] = assigned
    A['max_posterior'] = Q.max(1)  # compatibility alias; this is reporting membership
    A['max_reporting_membership'] = Q.max(1)
    A['accepted_MUT_any'] = (X == 1).any(1)
    A['all_events_unknown'] = (X == -1).all(1)
    A['reporting_policy'] = POLICY
    rows, audit = [], []
    for sample in cells['sample'].drop_duplicates():
        ix = cells['sample'].eq(sample).to_numpy()
        n = int(ix.sum())
        for k, node in enumerate(active):
            expected = float(Q[ix, k].sum())
            count = int((ix & (assigned == node)).sum())
            rows.append(dict(sample=sample, node=int(node),
                composition='WT' if node == 0 else '|'.join(events.event_label[G[node].astype(bool)]),
                cells=count, all_cells=n, percent_all=100*count/n,
                expected_cells=expected, reported_percent_all=100*expected/n,
                posterior_expected_percent_all=100*expected/n, reporting_policy=POLICY))
        for j, event in enumerate(events.event_label):
            mask = ix & (G[assigned, j] == 1)
            audit.append(dict(sample=sample, event=event, all_cells=n,
                input_MUT=int((ix & (X[:, j] == 1)).sum()), assigned_MUT=int(mask.sum()),
                assigned_MUT_input_MUT=int((mask & (X[:, j] == 1)).sum()),
                assigned_MUT_input_WT=int((mask & (X[:, j] == 0)).sum()),
                assigned_MUT_input_unknown=int((mask & (X[:, j] == -1)).sum())))
    F = pd.DataFrame(rows)
    require(np.allclose(F.groupby('sample').reported_percent_all.sum(), 100),
            'All-cell frequency denominator failure')
    return Q, A, F, pd.DataFrame(audit)
