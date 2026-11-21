"""Non-mutating biological-review diagnostics for fitted clone assignments.

An uncalled mutation is not absent, and a reference RNA observation is not a
certified WT genotype. These diagnostics identify claims requiring review;
they neither cap clone frequency at detection frequency nor change a model.
"""
import numpy as np
import pandas as pd
from .integrity import require


def evidence_consistency_audit(X, cells, events, model, assignments):
    X = np.asarray(X)
    G = np.asarray(model['G'])
    active = np.asarray(model['active'])
    require(list(cells.cell_id.astype(str)) == list(model['cell_ids'].astype(str)) ==
            list(assignments.cell_id.astype(str)), 'Review cell identity/order mismatch')
    require(list(events.feature.astype(str)) == list(model['features'].astype(str)),
            'Review event identity/order mismatch')
    require(cells.cell_id.is_unique and events.feature.is_unique, 'Duplicate review identity')
    require(X.shape == (len(cells), len(events)) and np.isin(X, [-1, 0, 1]).all(),
            'Invalid review call matrix')
    require(G.ndim == 2 and G.shape[1] == len(events) and np.isin(G, [0, 1]).all(),
            'Invalid review genotype matrix')
    require(active.ndim == 1 and len(set(active)) == len(active) and
            np.issubdtype(active.dtype, np.integer) and
            ((active >= 0) & (active < len(G))).all(), 'Invalid review active nodes')
    nodes = assignments.node.to_numpy()
    require(np.isin(nodes, active).all(), 'Review assignment references inactive node')
    ll = model.get('log_evidence')
    if ll is not None:
        ll = np.asarray(ll)
        require(ll.shape == (len(cells), len(active)) and np.isfinite(ll).all(),
                'Invalid review data likelihoods')
    rows, pairs = [], []
    for sample in cells['sample'].drop_duplicates():
        ix = cells['sample'].eq(sample).to_numpy()
        for node in active:
            mask = ix & (nodes == node)
            for j, event in enumerate(events.event_label):
                mut = int((mask & (X[:, j] == 1)).sum())
                rows.append(dict(sample=sample, node=int(node), event=event,
                    feature=str(events.feature.iloc[j]), cells=int(mask.sum()),
                    all_cells=int(ix.sum()), genotype=int(G[node, j]),
                    accepted_MUT=mut,
                    accepted_reference=int((mask & (X[:, j] == 0)).sum()),
                    unknown=int((mask & (X[:, j] == -1)).sum()),
                    positive_calls_omitted_from_genotype=mut if not G[node, j] else 0))
        for a in range(len(active)):
            for b in range(a+1, len(active)):
                if 0 in (active[a], active[b]):
                    continue
                js = np.flatnonzero(G[active[a]] != G[active[b]])
                mask = ix & np.isin(nodes, [active[a], active[b]])
                record = dict(sample=sample, node_a=int(active[a]), node_b=int(active[b]),
                    distinguishing_features='|'.join(events.feature.iloc[js].astype(str)),
                    distinguishing_events='|'.join(events.event_label.iloc[js].astype(str)),
                    assigned_cells=int(mask.sum()),
                    distinguishing_loci_all_unknown=int((mask & (X[:, js] == -1).all(1)).sum()),
                    equal_data_likelihood_cells=None)
                if ll is not None:
                    record['equal_data_likelihood_cells'] = int((mask &
                        (np.abs(ll[:, a] - ll[:, b]) < 1e-8)).sum())
                pairs.append(record)
    event_table = pd.DataFrame(rows)
    summary = dict(status='requires biological review; not external validation',
        cells=len(cells), events=len(events), samples=int(cells['sample'].nunique()),
        omitted_positive_cell_event_calls=int(event_table.positive_calls_omitted_from_genotype.sum()),
        changed_assignments=False, excluded_cells=0, excluded_events=0,
        interpretation=('Positive-call omissions, absent discriminating observations and likelihood ties '
            'require explicit review. Sparse detection alone does not invalidate an inferred clone. '
            'Reference RNA is not proven mutation absence. Do not translate missing-call combinations '
            'into separate biological clones. MAP cell counts and posterior frequency sums differ.'))
    return event_table, pd.DataFrame(pairs), summary
