"""Shared helpers for the scALABLE-discover validation scripts. No path is stored in a result."""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def file_record(path):
    path = Path(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return {"name": path.name, "sha256": digest}


def matrix_digest(matrix):
    m = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    m.sort_indices()
    h = hashlib.sha256()
    for part in (m.indptr, m.indices, m.data):
        h.update(np.ascontiguousarray(part).tobytes())
    return h.hexdigest()


def adata_digest(adata):
    return {"shape": list(adata.shape), "X": matrix_digest(adata.X),
            "layers": {k: matrix_digest(v) for k, v in sorted(adata.layers.items())},
            "obs_names": hashlib.sha256("\n".join(adata.obs_names.astype(str)).encode()).hexdigest(),
            "var_names": hashlib.sha256("\n".join(adata.var_names.astype(str)).encode()).hexdigest(),
            "obs": hashlib.sha256(pd.util.hash_pandas_object(adata.obs.astype(str), index=True).values.tobytes()).hexdigest()}


def parse_h5(values):
    """--h5 PATH:SAMPLE, repeated."""
    out = []
    for value in values:
        path, _, sample = value.rpartition(":")
        out.append((path, sample))
    return out
