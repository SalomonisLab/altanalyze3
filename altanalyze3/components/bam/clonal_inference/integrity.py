"""Reject incomplete or altered inputs before writing or fitting."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_inputs(matrix_path, cells_path, events_path, contract_path, method):
    """Load the entire specified matrix; never inner-join, fill unknowns, or drop rows.

    Contract expected IDs must come from the independently reconciled source roster,
    not from an intersection selected by this workflow. Source records are verified.
    NPZ keys: X (-1/0/1), cell_ids, features. Extra keys are retained in the source file.
    """
    contract_path = Path(contract_path).resolve()
    contract = json.loads(contract_path.read_text())
    require(contract.get('schema_version') == 1, 'Unsupported contract schema')
    require(contract.get('method') == method, 'Contract method does not match requested method')
    require(bool(contract.get('source_precedence')), 'Document upstream source precedence')
    require(bool(contract.get('provenance_sources')), 'Independent provenance records are required')
    for record in contract['provenance_sources']:
        path = Path(record['path'])
        if not path.is_absolute():
            path = contract_path.parent / path
        require(path.is_file(), f'Unresolved provenance source: {path}')
        require(sha256(path) == record['sha256'], f'Changed provenance source: {path}')
    paths = {'matrix': matrix_path, 'cells': cells_path, 'events': events_path}
    for key, path in paths.items():
        require(sha256(path) == contract['input_sha256'][key], f'Changed {key} input')
    with np.load(matrix_path, allow_pickle=False) as z:
        X = z['X'].copy()
        ids = z['cell_ids'].astype(str)
        features = z['features'].astype(str)
    cells = pd.read_csv(cells_path, sep='\t', dtype={'cell_id': str, 'sample': str})
    events = pd.read_csv(events_path, sep='\t', dtype={'feature': str})
    require({'cell_id', 'sample'} <= set(cells.columns), 'Cell roster needs cell_id and sample')
    require('feature' in events, 'Event roster needs feature')
    require(not cells[['cell_id', 'sample']].isna().any().any(), 'Missing cell/sample identity')
    require(not events.feature.isna().any(), 'Missing event identity')
    require(list(ids) == contract['expected_cell_ids'] == list(cells.cell_id), 'Cell identity/order mismatch')
    require(list(features) == contract['expected_features'] == list(events.feature), 'Event identity/order mismatch')
    require(list(cells['sample']) == contract['expected_samples'], 'Sample membership mismatch')
    require(len(set(ids)) == len(ids), 'Duplicate cell identity')
    require(len(set(features)) == len(features), 'Duplicate event identity')
    require(X.shape == (len(ids), len(features)), 'Matrix dimensions do not match complete roster')
    require(len(ids) > 0 and len(features) > 0, 'Empty dataset')
    require(np.isin(X, [-1, 0, 1]).all(), 'Calls must be -1 unknown, 0 WT, or 1 MUT; nonfinite values rejected')
    return X, cells, events, contract


def numerical_gate():
    """Exercise the matrix shape that exposed a node-specific BLAS error in 2026.

    Do not bypass a failure. Repair/pin the numerical backend and rerun the gate.
    """
    bits = ((np.arange(2048)[:, None] >> np.arange(11)) & 1).astype(float)
    a = bits[:337].copy()
    b = np.log(.01 + .69 * bits).T
    expected = sum(a[:, j, None] * b[j] for j in range(11))
    require(np.allclose(a @ b, expected, atol=1e-10, rtol=0), 'BLAS matrix multiplication failed numerical gate')
    w = np.arange(2048, dtype=float)
    w /= w.sum()
    expected_vector = sum(b[:, j] * w[j] for j in range(2048))
    require(np.allclose(b @ w, expected_vector, atol=1e-10, rtol=0), 'BLAS matrix-vector product failed numerical gate')
