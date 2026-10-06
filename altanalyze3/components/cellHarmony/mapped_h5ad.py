"""Bounded H5AD import and matrix selection using disposable disk-backed arrays.

The scientific operations remain in cellHarmony_lite. These helpers only change
where their arrays live. No feature intersection, downsampling or expression
precision reduction is used. Sparse index widths are chosen without loss. Work files belong to one worker and can be deleted after it exits.
"""
from pathlib import Path
import copy
import uuid
import weakref
import mmap
import shutil
import os

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
try:
    from anndata.io import read_elem
except ImportError:
    from anndata.experimental import read_elem


def matrix_bytes(node):
    if isinstance(node, h5py.Dataset):
        return int(node.size * node.dtype.itemsize)
    return sum(matrix_bytes(v) for v in node.values())


def inspect_h5ad(path):
    """Read dimensions and uncompressed sizes without reading expression values."""
    with h5py.File(path, 'r') as f:
        if not all(k in f for k in ('X', 'obs', 'var')):
            raise ValueError('H5AD input must contain X, obs and var.')
        x = f['X']
        shape = x.shape if isinstance(x, h5py.Dataset) else x.attrs.get('shape')
        if shape is None or len(shape) != 2:
            raise ValueError('H5AD X must be a two-dimensional expression matrix.')
        size = sum(matrix_bytes(f[k]) for k in ('X', 'layers', 'raw', 'obsm', 'obsp') if k in f)
        return {'cells': int(shape[0]), 'genes': int(shape[1]),
                'matrix_bytes': size, 'x_bytes': matrix_bytes(x),
                'encoding': str(x.attrs.get('encoding-type', 'array')),
                'layers': sorted(f.get('layers', {}).keys()), 'has_raw': 'raw' in f,
                'has_pair_graphs': bool(len(f.get('obsp', {})))}


def needs_disk_backed_import(headers):
    """Shared upload/worker decision, based on uncompressed arrays, not file size."""
    return (len(headers) > 1 or sum(h['cells'] for h in headers) > 100_000
            or sum(h['matrix_bytes'] for h in headers) >= 1024**3)


def validate_batch(headers):
    if len(headers) < 2:
        return
    from importlib.metadata import version
    if tuple(int(part) for part in version('anndata').split('.')[:2]) < (0, 12):
        raise ValueError('Multiple H5AD uploads require AnnData 0.12 or newer. '
                         'The server administrator must rebuild using requirements.docker.txt.')
    if any(h['layers'] != headers[0]['layers'] for h in headers):
        raise ValueError('The H5AD files have different expression layers. Supply the same layers '
                         'in each file so counts and other measurements are preserved for every sample.')
    if any(h['has_raw'] for h in headers) and not all(h['has_raw'] for h in headers):
        raise ValueError('Only some H5AD inputs contain raw measurements. Supply consistent raw '
                         'slots so the merge can retain the measurements for every sample.')


def validate_merge_scales(paths):
    """Do not apply one normalization decision to incompatible input scales."""
    from anndata.io import sparse_dataset
    from altanalyze3.components.clustering.ICGS import infer_expression_scale
    scales = {}
    for name, path in paths.items():
        with h5py.File(path, 'r') as handle:
            node = handle['X']
            matrix = node if isinstance(node, h5py.Dataset) else sparse_dataset(node)
            n = matrix.shape[0]
            rows = (np.sort(np.random.default_rng(0).choice(n, 20_000, replace=False))
                    if n > 20_000 else np.arange(n))
            parts = [sp.csr_matrix(matrix[rows[start:start + 256], :])
                     for start in range(0, len(rows), 256)]
            sample = sp.vstack(parts, format='csr') if parts else sp.csr_matrix((0, matrix.shape[1]))
            scales[name] = infer_expression_scale(sample)['verdict']
    if len(set(scales.values()) - {'empty'}) > 1:
        detail = ', '.join(f'{name}: {scale}' for name, scale in scales.items())
        raise ValueError('H5AD inputs use different expression scales (' + detail + '). '
                         'Upload files with the same preprocessing, or analyze them as separate jobs; '
                         'merging these matrices would apply incorrect QC or normalization.')


def _merge_feature_annotations(tables, index):
    """Fill annotations for union features; AnnData merge='first' leaves gaps."""
    result = pd.DataFrame(index=index)
    identifiers = {'gene_symbols', 'gene_ids', 'ensembl_id'}
    for table in tables:
        for column in table:
            incoming = table[column].astype(object).reindex(index)
            if column not in result:
                result[column] = incoming
                continue
            present = result[column].notna() & incoming.notna()
            if column in identifiers and (result.loc[present, column] != incoming[present]).any():
                raise ValueError(f'H5AD inputs contain conflicting {column} annotations for the same '
                                 'feature IDs. Reconcile those annotations before merging; '
                                 'no input features have been removed.')
            result[column] = result[column].where(result[column].notna(), incoming)
    result = result.infer_objects()
    for column in result:
        if result[column].dtype == object and all(isinstance(v, str) for v in result[column].dropna()):
            result[column] = pd.Categorical(result[column])
    return result


class Workspace:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self._maps = []

    def array(self, shape, dtype):
        # mmap cannot represent an empty file.
        if not np.prod(shape):
            return np.empty(shape, dtype=dtype)
        required = int(np.prod(shape)) * np.dtype(dtype).itemsize
        if shutil.disk_usage(self.directory).free < required + 512 * 1024**2:
            raise OSError('Not enough free disk space for bounded H5AD processing. '
                          f'The next work array needs {required / 1024**3:.1f} GiB; '
                          'uploaded inputs have not been changed.')
        out = np.memmap(self.directory / (uuid.uuid4().hex + '.bin'),
                        mode='w+', dtype=dtype, shape=shape)
        self._maps.append((weakref.ref(out), Path(out.filename)))
        return out

    def release_pages(self):
        """Flush work arrays so the OS can reclaim them under container pressure."""
        live = []
        for ref, path in self._maps:
            array = ref()
            if array is None:
                try:
                    path.unlink(missing_ok=True)
                except PermissionError:
                    pass
                continue
            live.append((ref, path))
            array.flush()
            if hasattr(array._mmap, 'madvise'):
                array._mmap.madvise(mmap.MADV_DONTNEED)
                self._release_file_cache(path)
        self._maps = live

    @staticmethod
    def _release_file_cache(path, offset=0, length=0):
        """Release clean scratch-file cache after unmapping resident pages.

        Linux can retain MADV_DONTNEED pages in its active file cache, which
        still counts towards the container working-set safety guard. Flush
        writes before calling this helper; it never changes file contents.
        Other platforms retain their existing mmap release behavior.
        """
        if not hasattr(os, 'posix_fadvise'):
            return
        try:
            descriptor = os.open(path, os.O_RDONLY)
            try:
                os.posix_fadvise(descriptor, offset, length, os.POSIX_FADV_DONTNEED)
            finally:
                os.close(descriptor)
        except OSError:
            # Some mounts cannot advise their cache. The flushed file remains
            # valid, and the normal operating-system reclaim path still works.
            pass

    def release_range(self, array, start, stop, *, write=False):
        """Flush/release a contiguous vector range in its owning work mapping.

        Sparse matrices expose ordinary ndarray views of memmaps. Match their
        address range to the registered owner rather than copying them. Inputs
        must be flushed before a read-only operation uses write=False.
        """
        if start >= stop or array.ndim != 1:
            return
        address = array.ctypes.data
        for ref, path in self._maps:
            owner = ref()
            if owner is None:
                continue
            origin = owner.ctypes.data
            if origin <= address and address + array.nbytes <= origin + owner.nbytes:
                lo = address - origin + int(start) * array.dtype.itemsize
                hi = address - origin + int(stop) * array.dtype.itemsize
                lo = lo // mmap.PAGESIZE * mmap.PAGESIZE
                hi = min(owner.nbytes, ((hi + mmap.PAGESIZE - 1) // mmap.PAGESIZE) * mmap.PAGESIZE)
                if write:
                    owner._mmap.flush(lo, hi - lo)
                if hasattr(owner._mmap, 'madvise'):
                    owner._mmap.madvise(mmap.MADV_DONTNEED, lo, hi - lo)
                    self._release_file_cache(path, lo, hi - lo)
                return

    def dataset(self, node, dtype=None):
        if isinstance(node, h5py.Dataset):
            if node.dtype.kind not in 'biufc' or node.ndim == 0:
                return read_elem(node)
            out = self.array(node.shape, dtype or node.dtype)
            width = int(np.prod(node.shape[1:])) or 1
            step = max(1, 4_000_000 // width)
            for start in range(0, node.shape[0], step):
                out[start:start + step] = node[start:start + step]
                if (start // step) % 8 == 7:
                    self.release_pages()
            self.release_pages()
            return np.asarray(out)
        encoding = node.attrs.get('encoding-type', '')
        if encoding in ('csr_matrix', 'csc_matrix'):
            shape = tuple(node.attrs['shape'])
            width = np.int64 if max(node['data'].shape[0], *shape) >= np.iinfo(np.int32).max else np.int32
            arrays = [self.dataset(node['data']), self.dataset(node['indices'], width),
                      self.dataset(node['indptr'], width)]
            ctor = sp.csr_matrix if encoding == 'csr_matrix' else sp.csc_matrix
            matrix = ctor(tuple(arrays), shape=tuple(node.attrs['shape']), copy=False)
            return self.csr(matrix) if encoding == 'csc_matrix' else matrix
        return read_elem(node)

    def csr(self, matrix):
        """Convert CSC in bounded column blocks, preserving all stored values."""
        if sp.isspmatrix_csr(matrix):
            return matrix
        if not sp.issparse(matrix):
            n, p = matrix.shape
            step = max(1, 4_000_000 // max(1, p))
            ptr = self.array((n + 1,), np.int64)
            ptr[0] = 0
            for start in range(0, n, step):
                counts = np.count_nonzero(matrix[start:start + step], axis=1)
                ptr[start + 1:start + len(counts) + 1] = ptr[start] + np.cumsum(counts)
            nnz = int(ptr[-1])
            data = self.array((nnz,), matrix.dtype)
            width = np.int64 if max(nnz, n, p) >= np.iinfo(np.int32).max else np.int32
            indices = self.array((nnz,), width)
            for start in range(0, n, step):
                block = sp.csr_matrix(matrix[start:start + step])
                lo, hi = int(ptr[start]), int(ptr[start + block.shape[0]])
                data[lo:hi], indices[lo:hi] = block.data, block.indices
                if start // step % 8 == 7:
                    self.release_pages()
            self.release_pages()
            return sp.csr_matrix((data, indices, ptr), shape=(n, p), copy=False)
        n, p = matrix.shape
        counts = np.zeros(n, dtype=np.int64)
        for start in range(0, matrix.nnz, 4_000_000):
            counts += np.bincount(matrix.indices[start:start + 4_000_000], minlength=n)
        ptr = self.array((n + 1,), np.int64)
        ptr[0] = 0
        np.cumsum(counts, out=ptr[1:])
        data = self.array((matrix.nnz,), matrix.dtype)
        width = np.int64 if max(matrix.nnz, n, p) >= np.iinfo(np.int32).max else np.int32
        indices = self.array((matrix.nnz,), width)
        cursor = np.asarray(ptr[:-1]).copy()
        start = 0
        while start < p:
            end = min(p, max(start + 1, int(np.searchsorted(
                matrix.indptr, int(matrix.indptr[start]) + 4_000_000, side='right') - 1)))
            block = matrix[:, start:end].tocsr()
            count = np.diff(block.indptr)
            dest = np.repeat(cursor - block.indptr[:-1], count) + np.arange(block.nnz)
            data[dest] = block.data
            indices[dest] = block.indices + start
            cursor += count
            start = end
        return sp.csr_matrix((data, indices, ptr), shape=(n, p), copy=False)

    def select(self, matrix, rows=None, cols=None):
        # Inputs are read-only throughout this operation. Flush any preceding
        # changes once, then release only ranges touched by this selection.
        self.release_pages()
        rows = np.arange(matrix.shape[0]) if rows is None else np.asarray(rows)
        all_columns = cols is None
        cols = np.arange(matrix.shape[1]) if cols is None else np.asarray(cols)
        all_columns = all_columns or np.array_equal(cols, np.arange(matrix.shape[1]))
        if not sp.issparse(matrix):
            out = self.array((len(rows), len(cols)), matrix.dtype)
            step = max(1, min(512, 4_000_000 // max(1, matrix.shape[1])))
            for start in range(0, len(rows), step):
                out[start:start + step] = matrix[rows[start:start + step], :][:, cols]
                if start // step % 8 == 7:
                    self.release_pages()
            self.release_pages()
            return np.asarray(out)
        # With all genes retained, row lengths are already in the CSR pointers.
        # Avoid a complete expression pass merely to count entries, and avoid
        # re-indexing an unchanged column panel for each bounded row block.
        ptr = self.array((len(rows) + 1,), np.int64)
        ptr[0] = 0
        step = max(1, min(512, 4_000_000 // max(1, matrix.shape[1])))
        # Flushing every few rows repeatedly scans every mapped array, even
        # those untouched by this operation. Bound the pages touched instead:
        # source + destination bytes, plus page-boundary overhead for each row.
        page_budget = 256 * 1024**2
        touched_bytes = 0
        if all_columns:
            matrix = self.csr(matrix)
            pointer_rows = np.where(rows < 0, rows + matrix.shape[0], rows)
            if np.any(pointer_rows < 0) or np.any(pointer_rows >= matrix.shape[0]):
                raise IndexError('row index out of range')
            lengths = (matrix.indptr[pointer_rows + 1].astype(np.int64)
                       - matrix.indptr[pointer_rows].astype(np.int64))
            np.cumsum(lengths, dtype=np.int64, out=ptr[1:])
        else:
            # Arbitrary gene selection retains the original bounded two-pass path.
            for start in range(0, len(rows), step):
                block = matrix[rows[start:start + step], :]
                touched_bytes += block.data.nbytes + block.indices.nbytes + 2 * mmap.PAGESIZE * block.shape[0]
                block = block[:, cols]
                block = block.tocsr() if sp.issparse(block) else sp.csr_matrix(block)
                ptr[start + 1:start + block.shape[0] + 1] = block.indptr[1:] + ptr[start]
                if touched_bytes >= page_budget:
                    self.release_pages()
                    touched_bytes = 0
        self.release_range(ptr, 0, len(ptr), write=True)
        nnz = int(ptr[-1])
        data = self.array((nnz,), matrix.dtype)
        index_dtype = np.int64 if max(nnz, len(cols)) >= np.iinfo(np.int32).max else np.int32
        indices = self.array((nnz,), index_dtype)
        first_output = 0
        lowest_row, highest_row = matrix.shape[0], -1

        def release_selection(last_output):
            nonlocal first_output, lowest_row, highest_row
            if matrix.format == 'csr':
                self.release_range(data, first_output, last_output, write=True)
                self.release_range(indices, first_output, last_output, write=True)
                if highest_row >= lowest_row:
                    begin, end = int(matrix.indptr[lowest_row]), int(matrix.indptr[highest_row + 1])
                    self.release_range(matrix.data, begin, end)
                    self.release_range(matrix.indices, begin, end)
                    self.release_range(matrix.indptr, lowest_row, highest_row + 2)
            else:
                self.release_pages()
            first_output = last_output
            lowest_row, highest_row = matrix.shape[0], -1

        for start in range(0, len(rows), step):
            selected_rows = rows[start:start + step]
            block = matrix[selected_rows, :]
            normalized_rows = np.where(selected_rows < 0, selected_rows + matrix.shape[0], selected_rows)
            lowest_row = min(lowest_row, int(normalized_rows.min()))
            highest_row = max(highest_row, int(normalized_rows.max()))
            touched_bytes += block.data.nbytes + block.indices.nbytes + 2 * mmap.PAGESIZE * block.shape[0]
            if not all_columns:
                block = block[:, cols]
            block = block.tocsr() if sp.issparse(block) else sp.csr_matrix(block)
            lo, hi = int(ptr[start]), int(ptr[start + block.shape[0]])
            data[lo:hi], indices[lo:hi] = block.data, block.indices
            touched_bytes += block.data.nbytes + block.indices.nbytes
            if touched_bytes >= page_budget:
                release_selection(hi)
                touched_bytes = 0
        release_selection(nnz)
        return sp.csr_matrix((data, indices, ptr), shape=(len(rows), len(cols)), copy=False)

    def load(self, path):
        with h5py.File(path, 'r') as f:
            kwargs = {k: read_elem(f[k]) for k in ('obs', 'var', 'uns') if k in f}
            for slot in ('layers', 'obsm', 'varm', 'obsp', 'varp'):
                if slot in f:
                    kwargs[slot] = {k: self.dataset(v) for k, v in f[slot].items()}
            out = ad.AnnData(X=self.dataset(f['X']), **kwargs)
            if 'scalable_upload' in out.obs:
                for key in ('Library', 'sample'):
                    if key not in out.obs:
                        out.obs[key] = out.obs['scalable_upload'].astype(str)
            if 'raw' in f:
                raw = f['raw']
                out.raw = ad.AnnData(X=self.dataset(raw['X']), obs=out.obs.copy(),
                                    var=read_elem(raw['var']),
                                    varm={k: self.dataset(v) for k, v in raw.get('varm', {}).items()})
        out._matrix_workspace = self
        self.release_pages()
        return out

    def subset(self, obj, rows=None, cols=None):
        rows = np.arange(obj.n_obs) if rows is None else np.asarray(rows)
        cols = np.arange(obj.n_vars) if cols is None else np.asarray(cols)
        out = ad.AnnData(X=self.select(obj.X, rows, cols), obs=obj.obs.iloc[rows].copy(),
                         var=obj.var.iloc[cols].copy(), uns=copy.deepcopy(obj.uns))
        for k, v in obj.layers.items():
            out.layers[k] = self.select(v, rows, cols)
        for slot, positions in (('obsm', rows), ('varm', cols)):
            for k, v in getattr(obj, slot).items():
                getattr(out, slot)[k] = v.iloc[positions].copy() if isinstance(v, pd.DataFrame) else v[positions].copy()
        for slot, positions in (('obsp', rows), ('varp', cols)):
            for k, v in getattr(obj, slot).items():
                getattr(out, slot)[k] = self.select(v, positions, positions)
        if obj.raw is not None:
            out.raw = ad.AnnData(X=self.select(obj.raw.X, rows), obs=out.obs.copy(),
                                var=obj.raw.var.copy(), varm=dict(obj.raw.varm))
        out._matrix_workspace = self
        self.release_pages()
        return out

    def copy_matrix(self, matrix, dtype=None, copy_structure=False):
        if not sp.issparse(matrix):
            out = self.array(matrix.shape, dtype or matrix.dtype)
            step = max(1, 4_000_000 // max(1, matrix.shape[1]))
            for start in range(0, matrix.shape[0], step):
                out[start:start + step] = matrix[start:start + step]
                if start // step % 8 == 7:
                    self.release_pages()
            self.release_pages()
            return np.asarray(out)
        matrix = self.csr(matrix)
        data = self.array(matrix.data.shape, dtype or matrix.dtype)
        for start in range(0, matrix.nnz, 4_000_000):
            data[start:start + 4_000_000] = matrix.data[start:start + 4_000_000]
            if start % 32_000_000 == 0:
                self.release_pages()
        self.release_pages()
        indices, ptr = matrix.indices, matrix.indptr
        if copy_structure:
            indices = self.array(matrix.indices.shape, matrix.indices.dtype)
            ptr = self.array(matrix.indptr.shape, matrix.indptr.dtype)
            ptr[:] = matrix.indptr
            for start in range(0, matrix.nnz, 4_000_000):
                indices[start:start + 4_000_000] = matrix.indices[start:start + 4_000_000]
                if start % 32_000_000 == 0:
                    self.release_pages()
        self.release_pages()
        return sp.csr_matrix((data, indices, ptr), shape=matrix.shape, copy=False)

    def scale_probe(self, obj, random_state=0):
        """Materialize exactly the existing scale classifier's sampled rows."""
        self.release_pages()
        rows = (np.sort(np.random.default_rng(random_state).choice(obj.n_obs, 20_000, replace=False))
                if obj.n_obs > 20_000 else np.arange(obj.n_obs))
        def sample(matrix):
            parts = []
            csr = sp.issparse(matrix) and matrix.format == 'csr'
            if csr:
                sizes = matrix.indptr[rows + 1].astype(np.int64) - matrix.indptr[rows].astype(np.int64)
                offsets = np.r_[0, np.cumsum(sizes, dtype=np.int64)]
            start = 0
            while start < len(rows):
                stop = min(len(rows), start + (4096 if csr else 256))
                if csr:
                    stop = min(stop, max(start + 1, int(np.searchsorted(
                        offsets, int(offsets[start]) + 4_000_000, side='right') - 1)))
                selected = rows[start:stop]
                parts.append(sp.csr_matrix(matrix[selected, :]))
                if csr:
                    lo, hi = int(matrix.indptr[selected[0]]), int(matrix.indptr[selected[-1] + 1])
                    self.release_range(matrix.data, lo, hi)
                    self.release_range(matrix.indices, lo, hi)
                    self.release_range(matrix.indptr, int(selected[0]), int(selected[-1]) + 2)
                else:
                    self.release_pages()
                start = stop
            return sp.vstack(parts, format='csr') if parts else sp.csr_matrix((0, obj.n_vars), dtype=matrix.dtype)
        result = ad.AnnData(X=sample(obj.X))
        if 'counts' in obj.layers:
            result.layers['counts'] = sample(obj.layers['counts'])
        return result

    def qc(self, obj, min_genes, min_cells, min_counts, mit_percent):
        """Same ordered QC predicates as cellHarmony, one bounded final selection."""
        matrix = obj.layers.get('counts', obj.X)
        n, p = matrix.shape
        rows = np.ones(n, dtype=bool)
        genes = np.zeros(p, dtype=np.int64)
        step = max(1, min(512, 4_000_000 // max(1, p)))
        for start in range(0, n, step):
            block = matrix[start:start + step]
            detected = np.diff(block.tocsr().indptr) if sp.issparse(block) else np.sum(block > 0, axis=1)
            keep = detected >= min_genes if min_genes is not None else np.ones(len(detected), bool)
            rows[start:start + len(keep)] = keep
            selected = block[keep]
            genes += (np.bincount(selected.tocsr().indices, minlength=p) if sp.issparse(selected)
                      else np.sum(selected > 0, axis=0))
        print(f'Cells remaining after min_genes {min_genes} filtering: {int(rows.sum())}')
        columns = np.flatnonzero(genes >= min_cells) if min_cells is not None else np.arange(p)
        mito = obj.var_names[columns].str.upper().str.startswith('MT-')
        percentages = np.zeros(n, dtype=matrix.dtype if matrix.dtype.kind == 'f' else np.float64)
        counts_pass = np.ones(n, dtype=bool)
        for start in range(0, n, step):
            block = matrix[start:start + step, :][:, columns]
            totals = np.asarray(block.sum(axis=1)).ravel()
            if min_counts is not None:
                counts_pass[start:start + len(totals)] = totals >= min_counts
            if mit_percent is not None:
                percentages[start:start + len(totals)] = np.asarray(block[:, mito].sum(axis=1)).ravel() / np.maximum(totals, 1e-12) * 100
        rows &= counts_pass
        print(f'Cells remaining after min_counts {min_counts} filtering: {int(rows.sum())}')
        if mit_percent is not None:
            obj.obs['pct_counts_mt'] = percentages
            rows &= percentages < mit_percent
            print(f'Cells remaining after mito-percent filtering: {int(rows.sum())}')
        if rows.all() and len(columns) == p:
            return obj
        return self.subset(obj, np.flatnonzero(rows), columns)


def merge_h5ads(entries, directory):
    """Outer disk merge; retain source annotations and a unique upload identity."""
    from anndata.experimental import concat_on_disk
    from anndata.io import write_elem
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = {str(name): str(path) for path, name in entries}
    if len(paths) != len(entries):
        raise ValueError('Each H5AD upload needs a unique sample name.')
    # A partial counts layer would silently change QC for some samples.
    headers = [inspect_h5ad(p) for p in paths.values()]
    validate_batch(headers)
    # concat_on_disk does not implement pairwise concatenation. Expose other
    # slots through temporary links, then write the graphs as block diagonals.
    shadows = {}
    graph_keys = set()
    source_metadata = {}
    existing_columns = set()
    encodings = {}
    for path in paths.values():
        with h5py.File(path, 'r') as source:
            existing_columns.update(source['obs'].attrs.get('column-order', []))
            matrix_keys = ['X', 'raw/X'] + [f'{slot}/{key}' for slot in ('layers', 'obsm')
                                          for key in source.get(slot, {})]
            for key in matrix_keys:
                if key in source:
                    encodings.setdefault(key, set()).add(source[key].attrs.get('encoding-type', 'array'))
    convert_keys = {key for key, kinds in encodings.items()
                    if 'csc_matrix' in kinds or (len(kinds) > 1 and 'csr_matrix' in kinds)}
    conversion = Workspace(directory / 'format_work') if convert_keys else None

    def link_or_convert(target, key, node, full_key, source_path):
        if full_key in convert_keys and node.attrs.get('encoding-type') != 'csr_matrix':
            matrix = conversion.csr(conversion.dataset(node))
            write_elem(target, key, matrix, dataset_kwargs={'compression': 'lzf'})
            del matrix
            conversion.release_pages()
        elif any(path.startswith(full_key + '/') for path in convert_keys):
            group = target.create_group(key)
            group.attrs.update(dict(node.attrs))
            for child_key, child in node.items():
                link_or_convert(group, child_key, child, full_key + '/' + child_key, source_path)
        else:
            target[key] = h5py.ExternalLink(str(Path(source_path).resolve()), '/' + full_key)

    preserved_label = 'source_scalable_upload'
    while preserved_label in existing_columns:
        preserved_label = 'source_' + preserved_label
    for i, (name, path) in enumerate(paths.items()):
        shadow = directory / f'input-{i}.h5ad'
        with h5py.File(path, 'r') as source, h5py.File(shadow, 'w') as target:
            target.attrs.update(dict(source.attrs))
            for key in source:
                if key != 'obsp':
                    if key == 'var':
                        # Assemble every annotation from the source tables below.
                        # AnnData's first-column merge both leaves union gaps and
                        # cannot serialize string columns containing those gaps.
                        write_elem(target, key, pd.DataFrame(index=read_elem(source[key]).index))
                    elif key == 'obs' and 'scalable_upload' in source['obs']:
                        annotations = read_elem(source['obs']).rename(columns={'scalable_upload': preserved_label})
                        write_elem(target, key, annotations)
                    else:
                        link_or_convert(target, key, source[key], key, path)
            graph_keys.update(source.get('obsp', {}).keys())
            source_metadata[name] = {key: read_elem(source[key]) for key in ('uns', 'var') if key in source}
            if 'raw' in source:
                source_metadata[name]['raw_var'] = read_elem(source['raw']['var'])
        shadows[name] = str(shadow)
    # Validate after storage-format conversion so CSC row sampling cannot trigger
    # an unbounded sparse read inside AnnData.
    validate_merge_scales(shadows)
    output = directory / 'merged.h5ad'
    concat_on_disk(shadows, output, join='outer', merge='first', uns_merge='first',
                   label='scalable_upload', index_unique='::', fill_value=0,
                   max_loaded_elems=4_000_000)
    with h5py.File(output, 'r+') as target:
        annotations = _merge_feature_annotations(
            [record['var'] for record in source_metadata.values()], read_elem(target['var']).index)
        write_elem(target, 'var', annotations)
        if headers[0]['has_raw']:
            raw_shadows = {}
            for i, (name, path) in enumerate(shadows.items()):
                shadow = directory / f'raw-{i}.h5ad'
                with h5py.File(path, 'r') as source, h5py.File(shadow, 'w') as raw_file:
                    raw_file.attrs.update({'encoding-type': 'anndata', 'encoding-version': '0.1.0'})
                    raw_file['obs'] = h5py.ExternalLink(str(Path(path).resolve()), '/obs')
                    for key in ('X', 'var', 'varm'):
                        if key in source['raw']:
                            if key == 'var':
                                write_elem(raw_file, key, pd.DataFrame(index=read_elem(source['raw'][key]).index))
                            else:
                                raw_file[key] = h5py.ExternalLink(str(Path(path).resolve()), '/raw/' + key)
                    for key in ('varm', 'varp', 'layers', 'obsm', 'obsp', 'uns'):
                        if key not in raw_file:
                            write_elem(raw_file, key, {})
                raw_shadows[name] = str(shadow)
            raw_group = target.create_group('raw')
            concat_on_disk(raw_shadows, raw_group, join='outer', merge='first',
                           index_unique='::', max_loaded_elems=4_000_000, fill_value=0)
            for key in list(raw_group):
                if key not in ('X', 'var', 'varm'):
                    del raw_group[key]
            raw_group.attrs.update({'encoding-type': 'raw', 'encoding-version': '0.1.0'})
            annotations = _merge_feature_annotations(
                [record['raw_var'] for record in source_metadata.values()], read_elem(raw_group['var']).index)
            write_elem(raw_group, 'var', annotations)
            for shadow in raw_shadows.values():
                Path(shadow).unlink()
        uns = target.require_group('uns')
        uns.attrs.update({'encoding-type': 'dict', 'encoding-version': '0.1.0'})
        write_elem(uns, 'scalable_sources', source_metadata)
        n = sum(inspect_h5ad(p)['cells'] for p in paths.values())
        for key in sorted(graph_keys):
            group = target.require_group('obsp').create_group(key)
            group.attrs.update({'encoding-type': 'csr_matrix', 'encoding-version': '0.1.0',
                                'shape': np.asarray([n, n], dtype=np.int64)})
            data = None
            # Result dtype includes every supplied graph's precision.
            dtypes = []
            for path in paths.values():
                with h5py.File(path, 'r') as f:
                    if key in f.get('obsp', {}):
                        node = f['obsp'][key]
                        dtypes.append(node.dtype if isinstance(node, h5py.Dataset) else node['data'].dtype)
            data = group.create_dataset('data', (0,), maxshape=(None,), dtype=np.result_type(*dtypes), compression='lzf')
            indices = group.create_dataset('indices', (0,), maxshape=(None,), dtype=np.int64, compression='lzf')
            ptr = group.create_dataset('indptr', (n + 1,), dtype=np.int64)
            ptr[0] = 0
            row_offset = nnz = 0
            for path in paths.values():
                size = inspect_h5ad(path)['cells']
                with h5py.File(path, 'r') as f:
                    node = f.get('obsp', {}).get(key)
                    if node is None:
                        ptr[row_offset + 1:row_offset + size + 1] = nnz
                    else:
                        try:
                            from anndata.io import sparse_dataset
                        except ImportError:
                            from anndata.experimental import sparse_dataset
                        source = node if isinstance(node, h5py.Dataset) else sparse_dataset(node)
                        for start in range(0, size, 256):
                            block = sp.csr_matrix(source[start:start + 256, :])
                            end = nnz + block.nnz
                            data.resize((end,)); indices.resize((end,))
                            data[nnz:end] = block.data
                            indices[nnz:end] = block.indices + row_offset
                            ptr[row_offset + start + 1:row_offset + start + block.shape[0] + 1] = block.indptr[1:] + nnz
                            nnz = end
                row_offset += size
    for shadow in shadows.values():
        Path(shadow).unlink()
    return output
