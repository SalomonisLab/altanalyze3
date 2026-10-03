"""Bounded H5AD reads and shared statistics for cell-level differential inputs."""
import copy

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def is_disk_matrix(matrix):
    return isinstance(matrix, h5py.Dataset) or (not sp.issparse(matrix) and getattr(matrix, "format", None) in ("csr", "csc") and hasattr(matrix, "to_memory"))


def root_and_rows(adata):
    root = adata._adata_ref if adata.is_view else adata
    if not is_disk_matrix(root.X):
        return None
    if adata.is_view:
        rows = np.arange(root.n_obs)[adata._oidx]
    else:
        rows = slice(None)
    return root, rows


def read_rows(dataset, rows, columns=slice(None)):
    """Read sorted row selections using hyperslabs instead of HDF5 point queries."""
    store = getattr(dataset, '_analysis_row_store', None)
    if store is not None:
        return store.read(rows, columns)
    profiles = getattr(dataset, '_analysis_profiles', None)
    if profiles is not None:
        values, codes = profiles
        return values[codes[rows], columns]
    if isinstance(rows, slice):
        if isinstance(dataset, h5py.Dataset):
            return dataset[rows, columns]
        rows = np.arange(dataset.shape[0])[rows]
    rows = np.asarray(rows)
    if len(rows) and int(rows[-1]) - int(rows[0]) + 1 == len(rows) and np.all(np.diff(rows) == 1):
        return dataset[int(rows[0]):int(rows[-1]) + 1, columns]
    width = dataset.shape[1] if isinstance(columns, slice) else len(columns)
    if isinstance(columns, slice):
        width = len(range(*columns.indices(dataset.shape[1])))
    sparse_source = not isinstance(dataset, h5py.Dataset)
    if sparse_source:
        if not len(rows):
            return sp.csr_matrix((0, width), dtype=dataset.dtype)
        order = np.argsort(rows, kind='stable')
        ordered = rows[order]
        parts = []
        block_rows = max(1, min(512, 2_000_000 // max(1, dataset.shape[1])))
        first = int(ordered[0]) // block_rows * block_rows
        for start in range(first, int(ordered[-1]) + 1, block_rows):
            left, right = np.searchsorted(ordered, [start, start + block_rows])
            if right == left:
                continue
            slab = dataset[start:min(start + block_rows, dataset.shape[0]), :]
            parts.append(slab[ordered[left:right] - start, columns])
        result = sp.vstack(parts, format='csr')
        return result[np.argsort(order)]
    result = np.empty((len(rows), width), dtype=dataset.dtype)
    if not len(rows):
        return result
    # Cell-state indices are in source order. Preserve other selectors as well.
    order = np.argsort(rows, kind='stable')
    ordered = rows[order]
    block_rows = max(1, min(4096, 8_000_000 // max(1, width)))
    first = int(ordered[0]) // block_rows * block_rows
    for start in range(first, int(ordered[-1]) + 1, block_rows):
        left, right = np.searchsorted(ordered, [start, start + block_rows])
        if right == left:
            continue
        slab = dataset[start:min(start + block_rows, dataset.shape[0]), columns]
        result[order[left:right]] = slab[ordered[left:right] - start]
    return result


class FeatureReader:
    def __init__(self, dataset, rows):
        self.dataset, self.rows = dataset, rows
        self.shape = (dataset.shape[0] if isinstance(rows, slice) else len(rows), dataset.shape[1])

    def __getitem__(self, key):
        rows, cols = key
        if rows != slice(None):
            raise NotImplementedError('FeatureReader reads all selected rows')
        return read_rows(self.dataset, self.rows, cols)


def expression_reader(adata):
    found = root_and_rows(adata)
    if found is None:
        return None
    root, rows = found
    dataset = root.layers.get('counts', root.X)
    logged, base = False, None
    if 'counts' not in root.layers:
        info = root.uns.get('log1p', {})
        logged, base = bool(info), info.get('base')
        if not logged:
            positions = np.arange(root.n_obs) if isinstance(rows, slice) else rows
            maximum = -np.inf
            for start in range(0, len(positions), 512):
                maximum = np.fmax(maximum, read_rows(dataset, positions[start:start + 512]).max())
            logged = bool(maximum <= 20)
    return FeatureReader(dataset, rows), np.asarray(root.var_names.astype(str)), logged, base


def write_selection(adata, path, uns):
    import scipy.sparse as sp
    root, rows = root_and_rows(adata)
    cols = np.arange(root.n_vars)[adata._vidx] if adata.is_view else np.arange(root.n_vars)
    positions = np.arange(root.n_obs) if isinstance(rows, slice) else rows
    var = adata.var.copy().drop(columns=['_index'], errors='ignore')
    uns = dict(uns)
    uns.pop('broadcast_profiles', None)
    uns.pop('broadcast_profile_codes', None)
    skeleton = ad.AnnData(X=sp.csr_matrix((len(positions), len(cols)), dtype=root.X.dtype),
                         obs=adata.obs.copy(), var=var, uns=uns)
    skeleton.write_h5ad(path, compression='lzf')
    with h5py.File(path, 'r+') as handle:
        del handle['X']
        datasets = [('X', root.X)] + [('layers/' + key, value) for key, value in root.layers.items()]
        written = {}
        for key, source in datasets:
            if isinstance(source, h5py.Dataset) and source.id in written:
                handle[key] = handle[written[source.id]]
                continue
            if not isinstance(source, h5py.Dataset):
                group = handle.create_group(key)
                group.attrs.update({'encoding-type': 'csr_matrix', 'encoding-version': '0.1.0',
                                    'shape': np.asarray([len(positions), len(cols)], dtype=np.int64)})
                data = group.create_dataset('data', shape=(0,), maxshape=(None,), dtype=source.dtype, compression='lzf')
                indices = group.create_dataset('indices', shape=(0,), maxshape=(None,), dtype=np.int32, compression='lzf')
                pointers = group.create_dataset('indptr', shape=(len(positions) + 1,), dtype=np.int64, compression='lzf')
                pointers[0] = 0
                offset = 0
                for start in range(0, len(positions), 256):
                    selected = positions[start:start + 256]
                    block = read_rows(source, selected)[:, cols].tocsr()
                    stop = offset + block.nnz
                    data.resize((stop,))
                    indices.resize((stop,))
                    data[offset:stop], indices[offset:stop] = block.data, block.indices
                    pointers[start + 1:start + len(selected) + 1] = block.indptr[1:].astype(np.int64) + offset
                    offset = stop
                continue
            target = handle.create_dataset(key, shape=(len(positions), len(cols)), dtype=source.dtype,
                                            compression='lzf')
            target.attrs.update({'encoding-type': 'array', 'encoding-version': '0.2.0'})
            for start in range(0, len(positions), 256):
                selected = positions[start:start + 256]
                block = read_rows(source, selected)[:, cols]
                target[start:start + len(selected)] = block.toarray() if sp.issparse(block) else block
            if isinstance(source, h5py.Dataset):
                written[source.id] = key


def rank_features(adata, groupby, case, control, method, reference_rank, bh_fdr, filter_mask):
    root, rows = root_and_rows(adata)
    reader = FeatureReader(root.X, rows)
    n_rows, n_genes = reader.shape
    metadata = copy.deepcopy(dict(adata.uns))
    normalize = False
    totals = None
    if 'log1p' not in metadata:
        maximum = -np.inf
        totals = np.empty(n_rows, dtype=np.float64)
        positions = np.arange(root.n_obs) if isinstance(rows, slice) else rows
        for start in range(0, n_rows, 512):
            block = read_rows(root.X, positions[start:start + 512])
            if sp.issparse(block):
                maximum = np.maximum(maximum, block.max())
                totals[start:start + block.shape[0]] = np.asarray(block.sum(axis=1)).ravel()
            else:
                block = np.asarray(block, dtype=np.float64)
                maximum = np.maximum(maximum, np.max(block))
                totals[start:start + len(block)] = block.sum(axis=1)
        normalize = bool(maximum > 20)
    width = max(1, min(256, 4_000_000 // max(1, n_rows)))
    frame_parts, keep_parts = [], []
    for start in range(0, n_genes, width):
        stop = min(start + width, n_genes)
        original = reader[:, start:stop]
        keep_parts.append(filter_mask(original))
        block = original
        uns = copy.deepcopy(metadata)
        if normalize:
            if sp.issparse(original):
                block = original.copy()
                if block.dtype.kind in 'iu':
                    block = block.astype(np.float32)
                divisor = totals.astype(block.dtype) / 1e4
                divisor[divisor == 0] = 1
                block.data = np.log1p(block.data / np.repeat(divisor, np.diff(block.indptr)))
            else:
                block = np.asarray(original, dtype=np.float64)
                divisor = totals / 1e4
                divisor[divisor == 0] = 1
                block = np.log1p(block / divisor[:, None])
            uns['log1p'] = {'base': None}
        current = ad.AnnData(X=block, obs=adata.obs.copy(), var=root.var.iloc[start:stop].copy(), uns=uns)
        names, _, logfc, pvals = reference_rank(current, groupby, case, control, method)
        scores = current.uns['rank_genes_groups']['scores'][case]
        frame_parts.append(pd.DataFrame({'pval': pvals.values, 'logfc': logfc.values, 'score': scores}, index=names))
    frame = pd.concat(frame_parts).reindex(root.var_names.astype(str))
    keep = np.concatenate(keep_parts) & np.isfinite(frame['pval'].to_numpy())
    adjusted = np.ones(n_genes, dtype=float)
    adjusted[keep] = bh_fdr(frame['pval'].to_numpy()[keep])
    frame['fdr'] = adjusted
    order = np.argsort(frame['score'].to_numpy())[::-1]
    frame = frame.iloc[order]
    return frame.index, frame['fdr'], frame['logfc'], frame['pval']


def bounded_materialize(adata, max_bytes=256 * 1024**2):
    """Use RAM for small selections, preserving a disk view above a fixed budget."""
    found = root_and_rows(adata)
    if found is None:
        return adata
    root, rows = found
    datasets = {'X': root.X, **dict(root.layers)}
    required = 0
    positions = np.arange(root.n_obs) if isinstance(rows, slice) else rows
    seen = set()
    for matrix in datasets.values():
        if isinstance(matrix, h5py.Dataset):
            if matrix.id in seen:
                continue
            seen.add(matrix.id)
            required += adata.n_obs * adata.n_vars * matrix.dtype.itemsize
        elif getattr(matrix, 'format', None) == 'csr':
            group = matrix._group
            pointers = group['indptr'][:]
            nnz = int(np.sum(pointers[positions + 1] - pointers[positions], dtype=np.int64))
            required += nnz * (group['data'].dtype.itemsize + group['indices'].dtype.itemsize)
            required += (adata.n_obs + 1) * pointers.dtype.itemsize
        else:
            return adata
    if required > max_bytes:
        return adata
    x = read_rows(root.X, rows)
    layers = {key: x if isinstance(matrix, h5py.Dataset) and isinstance(root.X, h5py.Dataset) and matrix.id == root.X.id
              else read_rows(matrix, rows) for key, matrix in root.layers.items()}
    metadata = copy.deepcopy(dict(adata.uns))
    metadata.pop('broadcast_profiles', None)
    metadata.pop('broadcast_profile_codes', None)
    return ad.AnnData(X=x, obs=adata.obs.copy(), var=adata.var.copy(), layers=layers, uns=metadata)


def moderated_inputs(adata, case_mask, control_mask, expression_threshold, min_samples):
    """Two sequential passes, preserving float32 row reduction order and global filtering."""
    root, rows = root_and_rows(adata)
    positions = np.arange(root.n_obs) if isinstance(rows, slice) else rows
    dtype = root.X.dtype
    if np.dtype(dtype).kind in 'iu':
        dtype = np.float64
    sums = [np.zeros(root.n_vars, dtype=dtype), np.zeros(root.n_vars, dtype=dtype)]
    detected = np.zeros(root.n_vars, dtype=np.int64)
    width = max(1, min(512, 2_000_000 // max(1, root.n_vars)))
    masks = [np.asarray(case_mask), np.asarray(control_mask)]
    def blocks():
        for start in range(0, len(positions), width):
            block = read_rows(root.X, positions[start:start + width])
            yield start, block.toarray() if sp.issparse(block) else block
    for start, block in blocks():
        for i, mask in enumerate(masks):
            selected = block[mask[start:start + len(block)]]
            if len(selected):
                sums[i] = np.cumsum(np.vstack([sums[i], selected]), axis=0, dtype=dtype)[-1]
                detected += np.sum(selected > expression_threshold, axis=0)
    sizes = [int(mask.sum()) for mask in masks]
    means = [sums[i] / sizes[i] for i in range(2)]
    ss = [np.zeros(root.n_vars, dtype=dtype), np.zeros(root.n_vars, dtype=dtype)]
    for start, block in blocks():
        for i, mask in enumerate(masks):
            selected = block[mask[start:start + len(block)]]
            if len(selected):
                residuals = selected - means[i]
                residuals *= residuals
                ss[i] = np.cumsum(np.vstack([ss[i], residuals]), axis=0, dtype=dtype)[-1]
    return means[0], means[1], ss[0] / (sizes[0] - 1), ss[1] / (sizes[1] - 1), detected >= min_samples


def bind_broadcast_profiles(adata):
    """Lossless factorization of predictions broadcast from sample/state profiles."""
    values = adata.uns.get('broadcast_profiles')
    codes = adata.uns.get('broadcast_profile_codes')
    if values is None or codes is None:
        return
    if values.shape[1] != adata.n_vars or len(codes) != adata.n_obs or np.any(codes < 0) or np.any(codes >= len(values)):
        raise ValueError('Invalid broadcast prediction profile metadata')
    adata.X._analysis_profiles = values, codes
    if 'counts' in adata.layers:
        scale = str(adata.uns.get('expression_scale', 'linear'))
        info = adata.uns.get('log1p', {})
        if scale == 'log2':
            counts = np.maximum(np.exp2(values.astype(np.float64)) - 1, 0).astype(np.float32)
        elif scale == 'log1p':
            base = info.get('base')
            base = np.e if base is None else float(base)
            linear = np.expm1(values.astype(np.float64)) if np.isclose(base, np.e) else np.power(base, values.astype(np.float64)) - 1
            counts = np.maximum(linear, 0).astype(np.float32)
        else:
            counts = np.maximum(values, 0)
        adata.layers['counts']._analysis_profiles = counts, codes



class RowStore:
    """An anonymous, uncompressed temporary matrix with bounded resident pages.

    The kernel removes the file when the worker exits, including after SIGKILL.
    H5AD remains the portable saved artifact; this store exists only during DE.
    """
    def __init__(self, dataset):
        import tempfile
        self.file = tempfile.TemporaryFile()
        self.array = np.memmap(self.file, mode='w+', dtype=dataset.dtype, shape=dataset.shape)
        width = max(1, 16_000_000 // max(1, dataset.shape[1]))
        for start in range(0, dataset.shape[0], width):
            self.array[start:start + width] = dataset[start:start + width]
            self.array.flush()
            self.release_pages()

    def release_pages(self):
        import mmap
        import sys
        if sys.platform != 'darwin' and hasattr(self.array._mmap, 'madvise') and hasattr(mmap, 'MADV_DONTNEED'):
            self.array._mmap.madvise(mmap.MADV_DONTNEED)
        else:
            # macOS treats DONTNEED as an aging hint; reopening the mapping
            # actually releases resident pages without deleting the file.
            shape, dtype = self.array.shape, self.array.dtype
            self.array._mmap.close()
            self.array = np.memmap(self.file, mode='r+', dtype=dtype, shape=shape)

    def read(self, rows, columns):
        positions = np.arange(self.array.shape[0])[rows] if isinstance(rows, slice) else np.asarray(rows)
        width = len(range(*columns.indices(self.array.shape[1]))) if isinstance(columns, slice) else len(columns)
        out = np.empty((len(positions), width), dtype=self.array.dtype)
        for start in range(0, len(positions), 512):
            out[start:start + 512] = self.array[positions[start:start + 512], columns]
            self.release_pages()
        return out


def bind_row_store(adata):
    """Stage dense GRN values once, avoiding repeated HDF5 decompression."""
    if not isinstance(adata.X, h5py.Dataset) or hasattr(adata.X, '_analysis_profiles'):
        return
    cache = RowStore(adata.X)
    adata.X._analysis_row_store = cache
    for matrix in adata.layers.values():
        if isinstance(matrix, h5py.Dataset) and matrix.id == adata.X.id:
            matrix._analysis_row_store = cache


def write_unfiltered(adata, de_store, path):
    """Serialize immutable expression buffers without a second full matrix copy."""
    metadata = copy.deepcopy(dict(adata.uns))
    metadata.pop('broadcast_profiles', None)
    metadata.pop('broadcast_profile_codes', None)
    metadata['cellHarmony_DE'] = de_store
    aliases = [key for key, matrix in adata.layers.items() if matrix is adata.X]
    result = ad.AnnData(X=adata.X, obs=adata.obs.copy(),
                        var=adata.var.copy().drop(columns=['_index'], errors='ignore'),
                        uns=metadata,
                        layers={key: value for key, value in adata.layers.items() if key not in aliases},
                        obsm=dict(adata.obsm), varm=dict(adata.varm), obsp=dict(adata.obsp), varp=dict(adata.varp))
    result.write_h5ad(path, compression='lzf')
    if aliases:
        with h5py.File(path, 'r+') as handle:
            for key in aliases:
                handle['layers/' + key] = handle['X']
