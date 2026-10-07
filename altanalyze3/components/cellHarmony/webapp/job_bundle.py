"""Serve a finished upload job's matrices from its bundle (LARGE_DATASET_DESIGN.md, step B).

A job of CELLHARMONY_BUNDLE_MIN_CELLS cells or more (default 10,000) writes ``outputs/bundle/`` at the end of
the pipeline (flask/pipeline.py ``_build_job_bundle``): the RNA store and one store per
imputed modality, in the scalable_viewer's gene-major memory-mapped format. For such a
job, the three places app.py used to read a result h5ad whole get a ``JobBundleAnnData``
instead:

  * ``_get_expression_cache``   the RNA or modality h5ad behind Explore and the figures
  * ``_open_gene_detail_adata`` the per-cell GRN edge h5ad behind the edge violin
  * ``_grn_edges_adata``        the per-cell GRN edge h5ad, when no pseudobulk network file

Only the serving matrix moves. ``obs``, ``var``, ``uns`` and the first two
columns of each ``obsm`` entry are still read from the same result h5ad, with anndata's own
element reader, so every label, category order, sample field and embedding is the object
``ad.read_h5ad`` would have built. The matrix values come from a store that
``precompute.py`` wrote from that same h5ad. Log-transformed RNA uses float16 storage
and returns float32 slices with small quantization differences. Other stores retain
float32 and exact values. Canonical analytical h5ads and raw counts are untouched.
What the h5ad path never loaded here: ``X``, ``layers`` and the wide ``obsm`` modality copies.

A bundle that does not line up with its h5ad (cell or feature order, counts) is refused,
the refusal is logged, and the caller falls back to reading the h5ad.
"""
from __future__ import annotations

import logging
import os
import threading
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:                                     # anndata >= 0.11
    from anndata.io import read_elem as _read_elem
except ImportError:                      # anndata 0.10
    from anndata.experimental import read_elem as _read_elem

from .memory_cache import BoundedCache, CacheBudget

LOG = logging.getLogger(__name__)
_BUDGET = CacheBudget(max_bytes=512 * 1024**2)

_LOCK = threading.Lock()
_DATASETS: Dict[Tuple[str, float], Any] = BoundedCache(_BUDGET)
_VIEWS: Dict[Tuple[str, str, str, float], "JobBundleAnnData"] = BoundedCache(_BUDGET)
_REFUSED: Dict[Tuple[str, str], str] = BoundedCache(_BUDGET)


# ------------------------------------------------------------------ locating a bundle


def bundle_info(meta: Dict) -> Optional[Dict[str, Any]]:
    """The job's bundle record when it was built and its files are present, else None."""
    info = (meta or {}).get("bundle") or {}
    if str(info.get("status") or "") != "completed":
        return None
    bundle_dir = str(info.get("dir") or "")
    prefix = str(info.get("prefix") or "")
    if not bundle_dir or not prefix:
        return None
    if not os.path.isfile(os.path.join(bundle_dir, f"{prefix}_metadata.json")):
        return None
    return info


def _dataset(info: Dict[str, Any]):
    from altanalyze3.components.visualization.scalable_viewer import bundle as B
    from altanalyze3.components.visualization.scalable_viewer import data_api as da

    bundle_dir, prefix = str(info["dir"]), str(info["prefix"])
    stamp = os.path.getmtime(os.path.join(bundle_dir, f"{prefix}_metadata.json"))
    key = (bundle_dir, stamp)
    ds = _DATASETS.get(key)
    if ds is None:
        ds = da.Dataset(B.BundlePaths(bundle_dir, prefix))
        _DATASETS[key] = ds
    return ds


def _same_file(a: str, b: str) -> bool:
    try:
        return bool(a) and bool(b) and os.path.samefile(a, b)
    except OSError:
        return False


def store_for(meta: Dict, h5ad_path) -> Optional[str]:
    """Which bundle store holds the matrix of ``h5ad_path``: 'rna', a modality id, or None."""
    info = bundle_info(meta)
    if info is None:
        return None
    path = str(h5ad_path or "")
    for store_id, source in (info.get("sources") or {}).items():
        if _same_file(path, str(source)):
            return str(store_id)
    return None


def view(meta: Dict, h5ad_path) -> Optional["JobBundleAnnData"]:
    """A bundle-backed AnnData stand-in for ``h5ad_path``, or None to read the h5ad.

    None when the job has no completed bundle, when the h5ad is not one of the bundle's
    sources, or when the bundle does not line up with the h5ad (logged once)."""
    store_id = store_for(meta, h5ad_path)
    if store_id is None:
        return None
    info = bundle_info(meta)
    path = str(h5ad_path)
    stamp = os.path.getmtime(path)
    key = (str(info["dir"]), path, store_id, stamp)
    with _LOCK:
        cached = _VIEWS.get(key)
        if cached is not None:
            return cached
        refused = _REFUSED.get((str(info["dir"]), path))
        if refused is not None:
            return None
        try:
            built = JobBundleAnnData(_dataset(info), store_id, path)
        except (ValueError, KeyError, FileNotFoundError, OSError) as exc:
            _REFUSED[(str(info["dir"]), path)] = f"{type(exc).__name__}: {exc}"
            LOG.warning("job bundle %s refused for %s: %s; reading the h5ad instead",
                        info["dir"], path, exc)
            return None
        _VIEWS[key] = built
        return built


def refusals() -> Dict[Tuple[str, str], str]:
    """Every (bundle, h5ad) pair this process refused, with the reason."""
    return dict(_REFUSED)


# ------------------------------------------------------------------ the AnnData stand-in


class _Column:
    """What ``adata[rows, gene]`` returns. app.py reads only ``.X``."""

    __slots__ = ("X",)

    def __init__(self, x):
        self.X = x


class _LazyObsm:
    """``adata.obsm``: the h5ad's keys in the h5ad's order, first two columns only.

    app.py reads an obsm entry only as a 2-D map (``_obsm_embedding_keys``,
    ``_coordinates_for_key``), and only its first two columns. Reading those two columns
    keeps the wide modality copies (``X_metabolite`` is 42,796 x 2,533 on the marrow job)
    out of memory. An entry with one column comes back with one, so it is excluded exactly
    as before."""

    def __init__(self, path: str):
        self._path = path
        self._cache: Dict[str, np.ndarray] = {}
        with h5py.File(path, "r") as fh:
            self._keys = list(fh["obsm"].keys()) if "obsm" in fh else []

    def __iter__(self):
        return iter(self._keys)

    def __len__(self):
        return len(self._keys)

    def __contains__(self, key):
        return key in self._keys

    def keys(self):
        return list(self._keys)

    def __getitem__(self, key):
        if key not in self._keys:
            raise KeyError(key)
        got = self._cache.get(key)
        if got is None:
            with h5py.File(self._path, "r") as fh:
                node = fh["obsm"][key]
                if isinstance(node, h5py.Dataset) and node.ndim == 2:
                    got = np.asarray(node[:, : min(2, node.shape[1])])
                else:                                    # sparse or dataframe entry
                    got = _read_elem(node)
                    if hasattr(got, "iloc"):
                        got = got.iloc[:, :2]
                    elif sp.issparse(got):
                        got = got[:, :2].toarray()
            self._cache[key] = got
        return got


class _StoreMatrix:
    """``adata.X`` over a gene-major store: the operations app.py applies, and no others.

    app.py and its helpers apply to ``X``: ``X.shape``, ``indicator @ X`` (a sparse group
    indicator on the left), ``X.sum(axis=0)``, ``X[:, j]`` and ``X[:, cols]`` (dot and comb
    plots, cross_modal), and ``X[rows]`` (GRN network, grn_data, integration_data).
    Each is computed in the order scipy or numpy uses on the h5ad's matrix. Float32
    stores agree bit for bit; float16 RNA stores sum their quantized values in float32
    (tests/test_job_bundle_serving.py). Anything else raises, so a
    new use surfaces instead of returning a wrong number."""

    __array_priority__ = 1000            # a numpy operand defers to __rmatmul__

    def __init__(self, owner: "JobBundleAnnData"):
        self._o = owner
        self.shape = (owner.n_obs, owner.n_vars)
        self.ndim = 2
        self.dtype = np.dtype(np.float32)

    def _columns(self, start: int, stop: int):
        """(cell index in h5ad row order, value, feature position) for features [start, stop)."""
        o = self._o
        indptr = np.asarray(o._indptr[start:stop + 1], dtype=np.int64)
        a, b = int(indptr[0]), int(indptr[-1])
        cells = np.asarray(o._indices[a:b], dtype=np.int64)
        values = np.asarray(o._data[a:b], dtype=np.float32)
        feature = np.repeat(np.arange(start, stop, dtype=np.int64), np.diff(indptr))
        if o._to_h5ad is not None:
            cells = o._to_h5ad[cells]
            # keep cells ascending inside each feature, the order a CSR row walk gives
            order = np.lexsort((cells, feature))
            cells, values, feature = cells[order], values[order], feature[order]
        return cells, values, feature

    def sum(self, axis=None):
        if axis != 0:
            raise NotImplementedError(f"X.sum(axis={axis!r}) is not served from a bundle")
        # scipy's csr.sum(axis=0) and numpy's dense sum(axis=0) both add each column's
        # entries in row order, in float32. np.cumsum is sequential, so its last value per
        # feature is that sum.
        n = self.shape[1]
        out = np.zeros(n, dtype=np.float32)
        step = 2048
        for s in range(0, n, step):
            e = min(s + step, n)
            _cells, values, feature = self._columns(s, e)
            if not values.size:
                continue
            ends = np.searchsorted(feature, np.arange(s, e), side="right")
            starts = np.searchsorted(feature, np.arange(s, e), side="left")
            for j in np.flatnonzero(ends > starts):
                out[s + j] = np.cumsum(values[starts[j]:ends[j]], dtype=np.float32)[-1]
        # scipy returns a (1, n) matrix for a sparse X and numpy a (n,) array for a dense
        # one; every caller flattens it, so a (1, n) array stands in for the matrix.
        return out.reshape(1, -1) if self._o._sparse else out

    def __rmatmul__(self, left):
        if not sp.issparse(left):
            raise NotImplementedError("only a sparse left operand is served from a bundle")
        left = sp.csr_matrix(left)
        if left.shape[1] != self.shape[0]:
            raise ValueError(f"shape mismatch: {left.shape} @ {self.shape}")
        coo = left.tocoo()
        per_cell = np.bincount(coo.col, minlength=self.shape[0])
        if per_cell.max(initial=0) > 1:
            raise NotImplementedError("a left operand with two entries in one column")
        group = np.full(self.shape[0], -1, dtype=np.int64)
        weight = np.zeros(self.shape[0], dtype=np.float64)
        group[coo.col] = coo.row
        weight[coo.col] = coo.data
        k, n = left.shape[0], self.shape[1]
        out = np.zeros((k, n), dtype=np.float64)
        step = 2048
        for s in range(0, n, step):
            e = min(s + step, n)
            cells, values, feature = self._columns(s, e)
            g = group[cells]
            keep = g >= 0
            if not keep.any():
                continue
            key = g[keep] * (e - s) + (feature[keep] - s)
            # scipy's csr @ csr adds a * x for each column in row order, in float64;
            # bincount adds in index order, which is cell order inside each feature.
            part = np.bincount(key, weights=weight[cells[keep]] * values[keep].astype(np.float64),
                               minlength=k * (e - s))
            out[:, s:e] = part.reshape(k, e - s)
        return out

    def group_sums_and_total(self, left):
        """The existing group contrast, with one bounded read of each feature.

        Preserve float64 weighted group accumulation and sequential float32
        totals exactly as __rmatmul__ and sum(axis=0) do separately.
        """
        left = sp.csr_matrix(left)
        if left.shape[1] != self.shape[0]:
            raise ValueError(f"shape mismatch: {left.shape} @ {self.shape}")
        coo = left.tocoo()
        if np.bincount(coo.col, minlength=self.shape[0]).max(initial=0) > 1:
            raise NotImplementedError("a left operand with two entries in one column")
        group = np.full(self.shape[0], -1, dtype=np.int64)
        weight = np.zeros(self.shape[0], dtype=np.float64)
        group[coo.col], weight[coo.col] = coo.row, coo.data
        k, n = left.shape[0], self.shape[1]
        out = np.zeros((k, n), dtype=np.float64)
        total = np.zeros(n, dtype=np.float32)
        s = 0
        while s < n:
            # At most one million stored entries, or one exceptionally dense
            # feature. Bound temporary arrays by nnz, not by dataset cell count.
            target = int(self._o._indptr[s]) + 1_000_000
            e = min(n, s + 2048, max(s + 1, int(np.searchsorted(self._o._indptr, target, side="right")) - 1))
            cells, values, feature = self._columns(s, e)
            g = group[cells]
            keep = g >= 0
            key = g[keep] * (e - s) + (feature[keep] - s)
            out[:, s:e] = np.bincount(key, weights=weight[cells[keep]] * values[keep].astype(np.float64),
                                     minlength=k * (e - s)).reshape(k, e - s)
            offsets = np.asarray(self._o._indptr[s:e + 1], dtype=np.int64)
            offsets = offsets - offsets[0]
            for j in np.flatnonzero(np.diff(offsets)):
                total[s + j] = np.cumsum(values[offsets[j]:offsets[j + 1]], dtype=np.float32)[-1]
            s = e
        return out, total

    # ---- indexing: X[rows], X[:, j], X[:, cols], X[rows, cols] ------------------------
    #
    # The result is the object the h5ad path returns, built from the same float32 values:
    # for the sparse RNA matrix a scipy CSR matrix (X[:, j] is CSR n x 1, as scipy returns),
    # for a dense modality a numpy array (X[:, j] is 1-D, as numpy returns). Downstream
    # arithmetic is then scipy's or numpy's own, on an identical operand.

    def _row_index(self, rows):
        n = self.shape[0]
        if isinstance(rows, slice):
            return None if rows == slice(None) else np.arange(n)[rows]
        idx = np.asarray(rows)
        if idx.dtype == bool:
            if idx.size != n:
                raise IndexError(f"boolean row mask of {idx.size} for {n} rows")
            return np.flatnonzero(idx)
        idx = idx.astype(np.int64).ravel()
        return np.where(idx < 0, idx + n, idx)

    def _col_index(self, cols):
        n = self.shape[1]
        if isinstance(cols, slice):
            return (None if cols == slice(None) else np.arange(n)[cols]), False
        if isinstance(cols, (int, np.integer)):
            j = int(cols)
            return np.array([j + n if j < 0 else j], dtype=np.int64), True
        idx = np.asarray(cols)
        if idx.dtype == bool:
            return np.flatnonzero(idx), False
        if idx.dtype.kind in "OUS":
            idx = np.array([self._o._position[str(v)] for v in idx.ravel()], dtype=np.int64)
        idx = idx.astype(np.int64).ravel()
        return np.where(idx < 0, idx + n, idx), False

    def _dense_block(self, row_idx, col_idx) -> np.ndarray:
        rows_n = self.shape[0] if row_idx is None else row_idx.size
        out = np.zeros((rows_n, col_idx.size), dtype=np.float32)
        for k, j in enumerate(col_idx.tolist()):
            column = self._o._dense(int(j))
            out[:, k] = column if row_idx is None else column[row_idx]
        return out

    def _rows_all_features(self, row_idx):
        """Every feature for the selected rows, in selection order."""
        n_sel, n_feat = row_idx.size, self.shape[1]
        unique = np.unique(row_idx).size == n_sel
        pos = np.full(self.shape[0], -1, dtype=np.int64)
        if unique:
            pos[row_idx] = np.arange(n_sel)
        else:
            pos[np.unique(row_idx)] = np.arange(np.unique(row_idx).size)
        step = 2048
        if not self._o._sparse:
            out = np.zeros((int((pos >= 0).sum()), n_feat), dtype=np.float32)
            for s in range(0, n_feat, step):
                e = min(s + step, n_feat)
                cells, values, feature = self._columns(s, e)
                p = pos[cells]
                keep = p >= 0
                out[p[keep], feature[keep]] = values[keep]
            return out if unique else out[np.searchsorted(np.unique(row_idx), row_idx)]
        r_parts, c_parts, v_parts = [], [], []
        for s in range(0, n_feat, step):
            e = min(s + step, n_feat)
            cells, values, feature = self._columns(s, e)
            p = pos[cells]
            keep = p >= 0
            r_parts.append(p[keep]); c_parts.append(feature[keep]); v_parts.append(values[keep])
        rr = np.concatenate(r_parts) if r_parts else np.empty(0, np.int64)
        cc = np.concatenate(c_parts) if c_parts else np.empty(0, np.int64)
        vv = np.concatenate(v_parts) if v_parts else np.empty(0, np.float32)
        n_unique = int((pos >= 0).sum())
        out = sp.coo_matrix((vv, (rr, cc)), shape=(n_unique, n_feat), dtype=np.float32).tocsr()
        out.sort_indices()
        return out if unique else out[np.searchsorted(np.unique(row_idx), row_idx)]

    def __getitem__(self, key):
        rows, cols = (key if isinstance(key, tuple) and len(key) == 2 else (key, slice(None)))
        row_idx = self._row_index(rows)
        col_idx, scalar = self._col_index(cols)
        if col_idx is None:                                 # X[rows]: every feature
            if row_idx is None:
                raise NotImplementedError("X[:] would materialise the whole matrix")
            return self._rows_all_features(row_idx)
        block = self._dense_block(row_idx, col_idx)
        if self._o._sparse:
            return sp.csr_matrix(block)                    # scipy: X[:, j] is CSR n x 1
        return block[:, 0] if scalar else block             # numpy: X[:, j] is 1-D

    def __matmul__(self, other):
        raise NotImplementedError("X @ other is not served from a bundle")


class JobBundleAnnData:
    """The AnnData surface app.py uses on a job's matrix, backed by the job's bundle."""

    def __init__(self, ds, store_id: str, h5ad_path: str):
        self.path = str(h5ad_path)
        self.store_id = str(store_id)
        with h5py.File(self.path, "r") as fh:
            self.obs = _read_elem(fh["obs"])
            var = _read_elem(fh["var"])
            self.uns = _read_elem(fh["uns"]) if "uns" in fh else {}
            x = fh["X"]
            self._sparse = isinstance(x, h5py.Group)
            shape = (tuple(int(v) for v in x.attrs["shape"]) if self._sparse
                     else tuple(int(v) for v in x.shape))
            x_dtype = np.dtype((x["data"] if self._sparse else x).dtype)
        # RNA bundles deliberately use reduced precision for display, including
        # float64-normalized ICGS3 outputs. Rejecting those sources reloads every
        # expression value for each view when the matrix exceeds the cache budget.
        # Analytical inputs and downloaded H5ADs retain their original precision.
        if x_dtype != np.float32 and not (self.store_id == "rna" and x_dtype.kind == "f"):
            raise ValueError(f"{self.path} stores X as {x_dtype}; the bundle holds float32")
        self.var = var
        self.var_names = var.index
        self.obs_names = self.obs.index
        self.n_obs, self.n_vars = int(shape[0]), int(shape[1])
        self.shape = (self.n_obs, self.n_vars)
        self.isbacked = False
        self.obsm = _LazyObsm(self.path)
        if len(self.obs) != self.n_obs or len(var) != self.n_vars:
            raise ValueError(f"{self.path}: obs/var ({len(self.obs)}, {len(var)}) do not "
                             f"match X {shape}")

        # ---- the store, and how its rows and columns line up with this h5ad ----------
        if self.store_id == "rna":
            if str(ds.sv.get("expr_dtype") or "float32") not in {"float16", "float32"}:
                raise ValueError(f"unsupported RNA store dtype: {ds.sv.get('expr_dtype')}")
            store = None
            features = [str(v) for v in ds._gene_ids] if getattr(ds, "_gene_ids", None) else None
            if features is None:
                ds._load_genes()
                features = [str(v) for v in ds._gene_ids]
            self._indptr, self._indices, self._data = ds.indptr, ds.indices, ds.data
        else:
            store = ds.modality(self.store_id)
            if store.kind != "per_cell":
                raise ValueError(f"modality '{self.store_id}' is a {store.kind} store; "
                                 "only a per-cell store stands in for a per-cell h5ad")
            if str(store.info.get("expr_dtype") or "float32") != "float32":
                raise ValueError(f"modality '{self.store_id}' store is "
                                 f"{store.info.get('expr_dtype')}, not float32")
            features = [str(v) for v in store.features]
            self._indptr, self._indices, self._data = store.indptr, store.indices, store.data
        names = [str(v) for v in self.var_names]
        if features != names:
            bad = sum(1 for a, b in zip(features, names) if a != b) + abs(len(features) - len(names))
            raise ValueError(f"{self.store_id}: bundle features differ from {self.path} var "
                             f"at {bad} of {len(names)} positions")
        from altanalyze3.components.visualization.scalable_viewer import bundle_meta as BM

        bundle_cells = np.asarray(BM._obs_names(ds), dtype=str)
        h5ad_cells = np.asarray(self.obs_names, dtype=str)
        if bundle_cells.size != h5ad_cells.size:
            raise ValueError(f"bundle holds {bundle_cells.size} cells, {self.path} "
                             f"{h5ad_cells.size}")
        if np.array_equal(bundle_cells, h5ad_cells):
            self._to_h5ad = None                        # bundle row i is h5ad row i
        else:
            where = {c: i for i, c in enumerate(h5ad_cells.tolist())}
            if len(where) != h5ad_cells.size:
                raise ValueError(f"{self.path} repeats a cell barcode")
            mapped = np.array([where.get(c, -1) for c in bundle_cells.tolist()], dtype=np.int64)
            if (mapped < 0).any():
                raise ValueError(f"{int((mapped < 0).sum())} bundle cells are absent from "
                                 f"{self.path}")
            self._to_h5ad = mapped
        self._position = {name: i for i, name in enumerate(names)}
        self.X = _StoreMatrix(self)

    # ---- the slices app.py takes -----------------------------------------------------

    def _dense(self, j: int) -> np.ndarray:
        a, b = int(self._indptr[j]), int(self._indptr[j + 1])
        out = np.zeros(self.n_obs, dtype=np.float32)
        if b > a:
            cells = np.asarray(self._indices[a:b], dtype=np.int64)
            if self._to_h5ad is not None:
                cells = self._to_h5ad[cells]
            out[cells] = np.asarray(self._data[a:b], dtype=np.float32)
        return out

    def __getitem__(self, key):
        if not (isinstance(key, tuple) and len(key) == 2):
            raise TypeError(f"JobBundleAnnData supports adata[rows, feature] only; got {key!r}")
        rows, feature = key
        if not isinstance(feature, str):
            # adata[:, [features]] (cross_modal._group_means) and adata[:, j]: the h5ad path
            # gives a view whose .X is the (n, k) submatrix, CSR for RNA, dense otherwise.
            return _Column(self.X[rows, feature if not isinstance(feature, (int, np.integer))
                                  else [int(feature)]])
        j = self._position.get(feature)
        if j is None:
            raise KeyError(feature)
        column = self._dense(j)
        if not (isinstance(rows, slice) and rows == slice(None)):
            rows = np.asarray(rows)
            column = column[rows]
        # the h5ad path returns an (n, 1) matrix: CSR for the RNA matrix, dense otherwise
        column = column.reshape(-1, 1)
        return _Column(sp.csr_matrix(column) if self._sparse else column)

    @property
    def layers(self):
        raise AttributeError("layers are not served from a job bundle")

    def __repr__(self) -> str:
        return (f"JobBundleAnnData({self.store_id}, {self.n_obs} x {self.n_vars}, "
                f"source {Path(self.path).name})")
