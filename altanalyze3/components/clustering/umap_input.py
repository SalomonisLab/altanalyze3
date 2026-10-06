"""Build exactly the float32 matrix UMAP consumes without a dense float64 copy."""
import numpy as np
import scipy.sparse as sp


def sparse_umap_input(adata, features, block_rows=2048, cells=None):
    """Read the complete ordered PCA panel without densifying all cells."""
    columns = adata.var_names.get_indexer(features)
    if np.any(columns < 0):
        raise ValueError('UMAP features must all be present; unresolved feature identities cannot be dropped.')
    rows = None if cells is None else adata.obs_names.get_indexer(cells)
    if rows is not None and np.any(rows < 0):
        raise ValueError('UMAP cell identities must all be present; no cells can be dropped.')
    n_rows = adata.n_obs if rows is None else len(rows)
    blocks = []
    for start in range(0, n_rows, block_rows):
        stop = min(start + block_rows, n_rows)
        selection = slice(start, stop) if rows is None else rows[start:stop]
        block = sp.csr_matrix(adata.X[selection, :][:, columns], copy=True)
        # Sum duplicate sparse entries in the source dtype before casting, just
        # as the historical dense UMAP input did.
        block.sum_duplicates()
        blocks.append(block.astype(np.float32, copy=False))
    return sp.vstack(blocks, format='csr') if blocks else sp.csr_matrix((0, len(columns)), dtype=np.float32)


def dense_umap_input(adata, features, block_rows=2048, cells=None):
    columns = adata.var_names.get_indexer(features)
    if np.any(columns < 0):
        raise ValueError('UMAP features must all be present; unresolved feature identities cannot be dropped.')
    rows = None if cells is None else adata.obs_names.get_indexer(cells)
    if rows is not None and np.any(rows < 0):
        raise ValueError('UMAP cell identities must all be present; no cells can be dropped.')
    n_rows = adata.n_obs if rows is None else len(rows)
    result = np.empty((n_rows, len(columns)), dtype=np.float32, order='C')
    for start in range(0, n_rows, block_rows):
        stop = min(start + block_rows, n_rows)
        selection = slice(start, stop) if rows is None else rows[start:stop]
        block = adata.X[selection, :][:, columns]
        # Densify before casting: duplicate sparse entries sum in the original
        # dtype, matching the historical toarray() followed by UMAP's cast.
        result[start:stop] = block.toarray() if sp.issparse(block) else np.asarray(block)
    return result
