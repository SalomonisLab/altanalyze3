"""Binary sparse-matrix stream that feeds cellHarmony without an intermediate file.

A producer writes one cells x genes CSR matrix, its names and optional cell annotations to a
pipe, and cellHarmony reads the pipe straight into an AnnData. No h5ad, mtx or 10x file is
written on the way. An R ``dgCMatrix`` holds genes x cells in CSC order, which is byte for byte
the cells x genes CSR matrix: ``@p`` is the indptr, ``@i`` the indices and ``@x`` the data. R
therefore streams its slots unchanged (``sparse_stream_from_seurat.R``).

Layout, little-endian:

====================  =========================================================================
magic                 8 bytes ``AH3SPS01``
n_obs, n_vars, nnz    3 x int64
codes                 4 x int32: indptr, indices, data, flags (codes: 1 int32, 2 int64,
                      3 float32, 4 float64; flags bit 1 = row totals follow the data)
obs_names             int64 byte length, then UTF-8 names joined by newline
var_names             int64 byte length, then UTF-8 names joined by newline
obs_table             int64 byte length, then a UTF-8 TSV whose first column repeats obs_names
                      (length 0 = no table)
indptr                n_obs + 1 values
indices               nnz values
data                  nnz values
row_totals            n_obs float64, present when flags bit 1 is set
trailer               8 bytes ``AH3SPEND``
====================  =========================================================================

The reader checks every structural invariant and, when the producer sent row totals, checks
that each streamed row sums to the producer's own total. A truncated or corrupt stream raises.
"""

from __future__ import annotations

import argparse
import io
import struct
import sys
from typing import BinaryIO, Optional, Sequence, Union

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

MAGIC = b"AH3SPS01"
TRAILER = b"AH3SPEND"
FLAG_ROW_TOTALS = 1
_CODES = {1: np.dtype("<i4"), 2: np.dtype("<i8"), 3: np.dtype("<f4"), 4: np.dtype("<f8")}
_HEADER = struct.Struct("<8sqqqiiii")


def _code_for(dtype: np.dtype) -> int:
    dtype = np.dtype(dtype)
    for code, known in _CODES.items():
        if dtype.kind == known.kind and dtype.itemsize == known.itemsize:
            return code
    raise TypeError(f"no stream code for dtype {dtype}")


def _readinto_exact(fh: BinaryIO, buf, what: str) -> None:
    view = memoryview(buf).cast("B")
    total = len(view)
    got = 0
    while got < total:
        n = fh.readinto(view[got:])
        if not n:
            raise EOFError(f"sparse stream ended after {got:,} of {total:,} bytes of {what}")
        got += n


def _read_bytes(fh: BinaryIO, n: int, what: str) -> bytes:
    buf = bytearray(n)
    _readinto_exact(fh, buf, what)
    return bytes(buf)


def _read_array(fh: BinaryIO, code: int, n: int, what: str) -> np.ndarray:
    if code not in _CODES:
        raise ValueError(f"unknown dtype code {code} for {what}")
    arr = np.empty(int(n), dtype=_CODES[code])
    if n:
        _readinto_exact(fh, arr, what)
    return arr.astype(arr.dtype.newbyteorder("="), copy=False)


def _read_text_block(fh: BinaryIO, what: str) -> str:
    (length,) = struct.unpack("<q", _read_bytes(fh, 8, f"{what} length"))
    if length < 0:
        raise ValueError(f"negative byte length {length} for {what}")
    return _read_bytes(fh, length, what).decode("utf-8") if length else ""


def _split_names(text: str, expected: int, what: str) -> list:
    names = text.split("\n") if expected else []
    if len(names) != expected:
        raise ValueError(f"{what}: header declares {expected:,} names, stream carries {len(names):,}")
    return names


def segment_sums(indptr: np.ndarray, data: np.ndarray, chunk: int = 65536) -> np.ndarray:
    """Exact float64 sum of each CSR row (or CSC column) segment.

    scipy's ``sum(axis, dtype=np.float64)`` still accumulates a float32 matrix in float32, so a
    sum above 2**24 comes back wrong: 16,778,217 ones summed to 16,777,216 in scipy 1.14.1. This
    casts one chunk of segments at a time to float64 and differences a cumulative sum.
    """
    n = len(indptr) - 1
    out = np.zeros(n, dtype=np.float64)
    for start in range(0, n, chunk):
        stop = min(n, start + chunk)
        lo, hi = int(indptr[start]), int(indptr[stop])
        cum = np.concatenate(([0.0], np.cumsum(data[lo:hi], dtype=np.float64)))
        ptr = np.asarray(indptr[start:stop + 1], dtype=np.int64) - lo
        out[start:stop] = cum[ptr[1:]] - cum[ptr[:-1]]
    return out


def _open_source(source) -> tuple:
    if hasattr(source, "readinto"):
        return source, False
    if str(source) == "-":
        return sys.stdin.buffer, False
    return open(source, "rb"), True


def read_sparse_stream(source: Union[str, BinaryIO], *, log=print) -> ad.AnnData:
    """Read one stream into an AnnData with CSR ``X``, ``obs`` from the stream's table.

    ``source`` is a path, ``"-"`` for standard input, or a binary file object. Values in the
    obs table stay strings; an empty field becomes NA.
    """
    fh, close = _open_source(source)
    try:
        head = _read_bytes(fh, _HEADER.size, "header")
        magic, n_obs, n_vars, nnz, c_ptr, c_idx, c_dat, flags = _HEADER.unpack(head)
        if magic != MAGIC:
            raise ValueError(f"not an altanalyze3 sparse stream: magic {magic!r}, expected {MAGIC!r}")
        if min(n_obs, n_vars, nnz) < 0:
            raise ValueError(f"negative dimension in header: n_obs={n_obs} n_vars={n_vars} nnz={nnz}")
        if _CODES.get(c_ptr, np.dtype("f8")).kind != "i" or _CODES.get(c_idx, np.dtype("f8")).kind != "i":
            raise ValueError(f"indptr and indices must be integer codes, got {c_ptr} and {c_idx}")
        log(f"[stream] header: {n_obs:,} cells x {n_vars:,} genes, {nnz:,} nonzero values, "
            f"data {_CODES[c_dat]}, row totals {'yes' if flags & FLAG_ROW_TOTALS else 'no'}")

        obs_names = _split_names(_read_text_block(fh, "obs_names"), n_obs, "obs_names")
        var_names = _split_names(_read_text_block(fh, "var_names"), n_vars, "var_names")
        obs_text = _read_text_block(fh, "obs_table")

        indptr = _read_array(fh, c_ptr, n_obs + 1, "indptr")
        indices = _read_array(fh, c_idx, nnz, "indices")
        data = _read_array(fh, c_dat, nnz, "data")
        totals = _read_array(fh, 4, n_obs, "row_totals") if flags & FLAG_ROW_TOTALS else None
        trailer = _read_bytes(fh, len(TRAILER), "trailer")
        if trailer != TRAILER:
            raise ValueError(f"sparse stream trailer is {trailer!r}, expected {TRAILER!r}")
        extra = fh.read(1)
        if extra:
            raise ValueError("sparse stream carries bytes after its trailer")
    finally:
        if close:
            fh.close()

    if indptr[0] != 0 or indptr[-1] != nnz:
        raise ValueError(f"indptr runs {indptr[0]}..{indptr[-1]}, expected 0..{nnz}")
    if n_obs and np.any(np.diff(indptr) < 0):
        raise ValueError("indptr decreases, so the row pointers are corrupt")
    if nnz and (indices.min() < 0 or indices.max() >= n_vars):
        raise ValueError(f"column index range {indices.min()}..{indices.max()} leaves 0..{n_vars - 1}")
    if data.dtype.kind == "f" and nnz and not np.isfinite(data).all():
        raise ValueError(f"{int((~np.isfinite(data)).sum()):,} data values are NaN or infinite")
    if len(set(obs_names)) != n_obs:
        raise ValueError(f"{n_obs - len(set(obs_names)):,} duplicate cell names in obs_names")
    n_dup_var = n_vars - len(set(var_names))
    if n_dup_var:
        log(f"[stream] {n_dup_var:,} duplicate gene names; cellHarmony makes them unique")

    if data.dtype.kind == "i":
        data = data.astype(np.float32 if data.dtype.itemsize == 4 else np.float64)
    X = sp.csr_matrix((data, indices, indptr), shape=(n_obs, n_vars))

    if totals is not None:
        streamed = segment_sums(indptr, data)
        integer_valued = bool(nnz == 0 or np.array_equal(data, np.rint(data)))
        if integer_valued:
            bad = np.flatnonzero(streamed != totals)
        else:
            bad = np.flatnonzero(~np.isclose(streamed, totals, rtol=1e-6, atol=1e-6))
        if len(bad):
            raise ValueError(
                f"{len(bad):,} of {n_obs:,} cells do not sum to the producer's total; first "
                f"{obs_names[bad[0]]}: streamed {streamed[bad[0]]} vs {totals[bad[0]]}")
        log(f"[stream] row totals match the producer for {n_obs:,} of {n_obs:,} cells "
            f"({'exact' if integer_valued else 'rtol 1e-6'}); matrix sum {streamed.sum():,.0f}")

    obs = pd.DataFrame(index=pd.Index(obs_names, name=None))
    if obs_text:
        table = pd.read_csv(io.StringIO(obs_text), sep="\t", dtype=str, keep_default_na=False,
                            na_values=[""])
        if len(table) != n_obs:
            raise ValueError(f"obs_table has {len(table):,} rows for {n_obs:,} cells")
        if not np.array_equal(table.iloc[:, 0].astype(str).to_numpy(), np.asarray(obs_names, dtype=object).astype(str)):
            raise ValueError("obs_table first column does not repeat obs_names in order")
        obs = table.iloc[:, 1:].set_index(pd.Index(obs_names))
        log(f"[stream] obs columns: {', '.join(obs.columns) or 'none'}")

    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(var_names)))


def write_sparse_stream(sink: Union[str, BinaryIO], X, obs_names: Sequence[str],
                        var_names: Sequence[str], obs: Optional[pd.DataFrame] = None, *,
                        row_totals: bool = True) -> None:
    """Write a cells x genes matrix in the stream layout. ``sink`` is a path, ``"-"`` or a file."""
    X = sp.csr_matrix(X)
    n_obs, n_vars = X.shape
    obs_names = [str(v) for v in obs_names]
    var_names = [str(v) for v in var_names]
    if len(obs_names) != n_obs or len(var_names) != n_vars:
        raise ValueError("name counts do not match the matrix shape")
    for what, names in (("obs_names", obs_names), ("var_names", var_names)):
        if any("\n" in n for n in names):
            raise ValueError(f"{what} contain a newline")
    int_dtype = np.int32 if max(X.nnz, n_vars) < 2**31 else np.int64
    indptr = X.indptr.astype(int_dtype, copy=False)
    indices = X.indices.astype(int_dtype, copy=False)
    data = X.data if X.data.dtype in (np.float32, np.float64, np.int32, np.int64) else X.data.astype(np.float64)

    obs_text = ""
    if obs is not None:
        table = obs.copy()
        table.insert(0, "obs_name", obs_names)
        cells = table.astype(object).where(table.notna(), "").astype(str)
        if cells.apply(lambda c: c.str.contains("[\t\n\r]", regex=True)).to_numpy().any():
            raise ValueError("an obs value contains a tab or newline")
        obs_text = cells.to_csv(sep="\t", index=False)

    if hasattr(sink, "write"):
        fh, close = sink, False
    elif str(sink) == "-":
        fh, close = sys.stdout.buffer, False
    else:
        fh, close = open(sink, "wb"), True
    try:
        fh.write(_HEADER.pack(MAGIC, n_obs, n_vars, X.nnz, _code_for(indptr.dtype),
                              _code_for(indices.dtype), _code_for(data.dtype),
                              FLAG_ROW_TOTALS if row_totals else 0))
        for text in ("\n".join(obs_names), "\n".join(var_names), obs_text):
            raw = text.encode("utf-8")
            fh.write(struct.pack("<q", len(raw)))
            fh.write(raw)
        for arr in (indptr, indices, data):
            fh.write(np.ascontiguousarray(arr, dtype=arr.dtype.newbyteorder("<")).tobytes())
        if row_totals:
            fh.write(segment_sums(X.indptr, X.data).astype("<f8").tobytes())
        fh.write(TRAILER)
        fh.flush()
    finally:
        if close:
            fh.close()


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate an altanalyze3 sparse stream and print its summary.")
    parser.add_argument("source", help="stream path, or - for standard input")
    args = parser.parse_args()
    adata = read_sparse_stream(args.source)
    print(f"[stream] valid: {adata.n_obs:,} cells x {adata.n_vars:,} genes, {adata.X.nnz:,} nonzero values")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
