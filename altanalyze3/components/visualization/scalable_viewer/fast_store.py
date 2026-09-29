"""Parallel builder for the gene-major expression store of a scalable_viewer bundle.

Writes the same bytes as ``precompute.build_expression_store`` (the reference) and returns
the same arrays, in a fraction of the time.

Design:

* Phase A, one task per row block, in parallel. Read the block once, transpose it with
  SciPy's C routine ``csr_tocsc`` (through the public ``csr_matrix.tocsc``), write the
  transposed block into the block's own two spill files, and compute the block's
  statistics. No task needs another task's result, so the source is read once.
* Phase B, one task per range of genes, in parallel. With every block's per-gene counts
  known, each task copies its genes' segments out of the spill, block after block, into
  its own contiguous range of the output.
* Spill and output move as whole slices through ordinary file reads and writes, not
  through memory maps, and no two concurrent writers ever share a file block (see the
  file I/O rules above ``_WRITE_ALIGN``).

An earlier two-pass layout (count every block's genes, then re-read each block and write
its segments straight into the output) was 18.3 s against 12.9 s for this one on the
COPD metacell layer at 8 workers (checks/fast_builder_20260928/copd_design_profile*.log),
and it read the source twice. It was removed.

What made the reference slow, measured 2026-09-28 on the COPD metacell layer (one block =
8,192 rows, 63,490,212 non-zeros): ``np.argsort(kind="stable")`` + ``np.unique`` took
9.5 s of every 10 s block. ``csr_tocsc`` takes 0.50 s on the same block and returns the
same order, because it walks rows in order, so cells stay ascending inside every gene.

Identity with the reference, and why it holds:

* Row blocks are the reference's blocks (``row_block`` rows each). Each block's sums are
  the additions the reference's ``np.bincount`` makes, in the same order (see
  ``_block_stats``), and the main process adds the block results in block order, as the
  reference does. Float addition is not associative; this order makes the float64 sums
  bit-equal.
* Stored values go through the reference's casts: source -> float32 -> ``expr_dtype``.
* Invariants fail loudly: each block's transpose holds its own non-zeros, every gene range
  receives exactly its counted entries, and the bytes written equal the store's size.

Disk: the spill takes 8 bytes per non-zero beside the store while the build runs (the
store's own size again), two files per row block, removed at the end, also when the
build fails.

Workers are processes for an h5ad (h5py serialises reads inside one process) and threads
for an in-memory CSR (``csr_tocsc`` releases the GIL).

``tests/test_fast_expression_store.py`` compares this module with the reference, byte for
byte, on synthetic layers.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import time
from concurrent.futures import Executor, ProcessPoolExecutor, ThreadPoolExecutor
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp

# Bytes one worker holds per non-zero of its block at its peak: the block's indices and
# values (8), their transpose (8) and the per-entry state (2), rounded up. It only sizes
# the worker count; the validation run reports the measured peak.
_BYTES_PER_NNZ_PEAK = 20
# Entries per statistics chunk. Any size gives the same sums; this one keeps the
# temporaries small.
_STATS_CHUNK = 1 << 22
# Phase B ranges per worker, for load balance.
_RANGES_PER_WORKER = 4


# ------------------------------------------------------------------ sources


class H5Source:
    """A CSR layer inside an h5ad file, named by its HDF5 group ('X' or 'layers/<name>')."""

    kind = "h5"

    def __init__(self, path: str, group: str):
        self.path = os.path.abspath(path)
        self.group = group

    def describe(self) -> str:
        return f"{self.path}:{self.group}"


class ArraySource:
    """A CSR layer already in memory: indptr, indices, data (a scipy csr_matrix's arrays)."""

    kind = "array"

    def __init__(self, indptr: np.ndarray, indices: np.ndarray, data: np.ndarray):
        self.indptr = indptr
        self.indices = indices
        self.data = data

    @classmethod
    def from_csr(cls, matrix) -> "ArraySource":
        if matrix.format != "csr":
            raise ValueError(f"expected a CSR matrix, got {matrix.format}")
        # Unsorted genes inside a row need no sort: the transpose keeps cells ascending
        # inside every gene either way, as the reference's stable sort does.
        return cls(matrix.indptr, matrix.indices, matrix.data)

    def describe(self) -> str:
        return f"in-memory CSR ({int(self.indptr[-1]):,} non-zeros)"


# Each worker process opens the h5ad once and keeps the three datasets.
_WORKER_H5: Dict[Tuple[str, str], Tuple[object, object, object]] = {}


def _datasets(src):
    """Return (indptr, indices, data) array-likes for ``src`` inside this process."""
    if src.kind == "array":
        return src.indptr, src.indices, src.data
    key = (src.path, src.group)
    got = _WORKER_H5.get(key)
    if got is None:
        import h5py
        grp = h5py.File(src.path, "r")[src.group]
        got = (grp["indptr"], grp["indices"], grp["data"])
        _WORKER_H5[key] = got
    return got


# ------------------------------------------------------------------ shared block steps


def _read_block(src, local_indptr: np.ndarray):
    """The reference's reads and casts: gene indices as stored, values as float32."""
    _, indices, data = _datasets(src)
    s = int(local_indptr[0])
    e = int(local_indptr[-1])
    gi = np.asarray(indices[s:e])
    dv = np.asarray(data[s:e], dtype=np.float32)
    return s, e, gi, dv


def _transpose_block(block: int, r0: int, r1: int, rel_indptr: np.ndarray, gi, dv,
                     n_genes: int):
    """Gene-major copy of one row block: (bp int64, global cell ids uint32, values float32).

    tocsc keeps duplicates and never sorts; it walks the rows in order, so cells stay
    ascending inside every gene, the order the reference's stable sort gives."""
    nnz_b = int(dv.size)
    blk = sp.csr_matrix((dv, gi, rel_indptr), shape=(int(r1 - r0), int(n_genes)), copy=False)
    csc = blk.tocsc()
    bp, bi, bx = csc.indptr, csc.indices, csc.data
    del blk, csc
    if bx.dtype != np.float32 or bx.size != nnz_b or int(bp[-1]) != nnz_b:
        raise RuntimeError(f"block {block}: transpose returned {bx.size} {bx.dtype} values, "
                           f"expected {nnz_b} float32")
    # The reference stores ci.astype(uint32), ci being the global cell index.
    if bi.dtype == np.int32:                        # our own array: shift it in place
        cells = bi.view(np.uint32)
    else:
        cells = bi.astype(np.uint32)
    cells += np.uint32(r0)
    return bp.astype(np.int64), cells, bx


def _block_stats(gi, dv, rel_indptr: np.ndarray, state_block: np.ndarray, n_genes: int,
                 n_states: int):
    """The reference's per-block sums, bit for bit.

    The reference calls np.bincount(key, weights=w) once per block, which adds w[i] into
    out[key[i]] for i = 0, 1, 2 ... starting from 0.0. np.add.at over consecutive chunks
    makes the same additions in the same order, so the float64 sums are bit-equal
    (checks/fast_builder_20260928/bench_addat.log: 0 of 1,810,600 entries differ). The
    chunks keep each worker's temporaries near 100 MB instead of about 2.5 GB, which
    pushed 8 workers into memory compression (stats 0.78 s -> 5.63 s a block)."""
    nnz_b = int(dv.size)
    gs_len = n_genes * n_states
    sum_b = np.zeros(gs_len, dtype=np.float64)
    cnt_b = np.zeros(gs_len, dtype=np.int64)
    sq_b = np.zeros(n_genes, dtype=np.float64)
    entry_state = np.repeat(state_block, np.diff(rel_indptr.astype(np.int64)))
    for c0 in range(0, nnz_b, _STATS_CHUNK):
        c1 = min(c0 + _STATS_CHUNK, nnz_b)
        g64 = gi[c0:c1].astype(np.int64)
        key = g64 * n_states
        key += entry_state[c0:c1]
        w = dv[c0:c1].astype(np.float64)
        np.add.at(sum_b, key, w)
        cnt_b += np.bincount(key, minlength=gs_len)
        np.add.at(sq_b, g64, w ** 2)
        del g64, key, w
    return sum_b, cnt_b, sq_b


def _out_dtype(expr_dtype: str):
    return np.float16 if expr_dtype == "float16" else np.float32


# ------------------------------------------------------------------ file I/O
#
# Two rules, both measured on this Mac (APFS) on 2026-09-28:
#
# 1. Whole-slice writes and reads, not memory maps. A fresh file-backed memmap costs one
#    page fault per 16 KB page, and 8 processes faulting on one file serialise: the spill
#    copy took up to 2.9 s a block and the gather ran near 110 MB/s a worker
#    (checks/fast_builder_20260928/copd_design_profile.log in the COPD-atlas viewer folder).
# 2. No two concurrent writers ever share a file block. When worker processes wrote
#    adjacent byte ranges of ONE sparse spill file, the first 440 bytes of block 5's range
#    came back as zeros, ending exactly at the next 4 KB boundary: a lost update inside a
#    block that two processes wrote (110 of 666,852,999 values; checks/fast_builder_20260928/
#    run2_shared_spill_design/copd_returned_arrays.log). The old code lost 400 index values
#    in 1 of 8 repeat runs (negative_control_shared_spill.log). So every block spills to
#    its own two files, and each gene range writes only the _WRITE_ALIGN-aligned middle
#    of its output; the main process writes the unaligned edges alone, after every
#    worker has finished.
_WRITE_ALIGN = 1 << 20


def _npy_data_offset(path: str) -> int:
    """Byte offset of the array data inside a .npy file (its header length)."""
    mm = np.load(path, mmap_mode="r")
    try:
        return int(mm.offset)
    finally:
        del mm


def _as_bytes(arr: np.ndarray) -> memoryview:
    return memoryview(np.ascontiguousarray(arr)).cast("B")


def _write_bytes(path: str, byte_offset: int, mv: memoryview, mode: str = "r+b") -> None:
    with open(path, mode, buffering=0) as fh:
        if byte_offset:
            fh.seek(byte_offset)
        done = 0
        while done < len(mv):
            n = fh.write(mv[done:])
            if not n:
                raise OSError(f"short write to {path}")
            done += n


def _read_at(path: str, elem_offset: int, count: int, dtype, base: int = 0) -> np.ndarray:
    out = np.empty(count, dtype=dtype)
    mv = memoryview(out).cast("B")
    with open(path, "rb", buffering=0) as fh:
        fh.seek(base + elem_offset * out.dtype.itemsize)
        done = 0
        while done < len(mv):
            n = fh.readinto(mv[done:])
            if not n:
                raise EOFError(f"short read from {path}")
            done += n
    return out


def _write_owned(path: str, byte_start: int, mv: memoryview) -> List[Tuple[int, bytes]]:
    """Write the _WRITE_ALIGN-aligned middle of [byte_start, byte_start + len(mv)).

    Returns the unaligned edges as (byte offset, bytes) for the main process to write.
    Two ranges never share an aligned window, so their middles never share a block."""
    end = byte_start + len(mv)
    lo = -(-byte_start // _WRITE_ALIGN) * _WRITE_ALIGN
    hi = (end // _WRITE_ALIGN) * _WRITE_ALIGN
    if lo >= hi:
        return [(byte_start, bytes(mv))]
    _write_bytes(path, lo, mv[lo - byte_start:hi - byte_start])
    edges = []
    if lo > byte_start:
        edges.append((byte_start, bytes(mv[:lo - byte_start])))
    if end > hi:
        edges.append((hi, bytes(mv[hi - byte_start:])))
    return edges


def _spill_paths(spill_prefix: str, block: int) -> Tuple[str, str]:
    return f"{spill_prefix}.b{block:05d}.idx", f"{spill_prefix}.b{block:05d}.val"


# ------------------------------------------------------------------ the two phases


def _scan_block(src, block: int, r0: int, r1: int, local_indptr: np.ndarray,
                state_block: np.ndarray, n_genes: int, n_states: int, spill_prefix: str):
    """Phase A for rows [r0, r1): read once, transpose into the block's own spill files,
    compute the block's statistics."""
    tm = {}
    t0 = time.perf_counter()
    s, e, gi, dv = _read_block(src, local_indptr)
    rel = local_indptr - s
    tm["read"] = time.perf_counter() - t0; t0 = time.perf_counter()

    bp, cells, bx = _transpose_block(block, r0, r1, rel, gi, dv, n_genes)
    tm["transpose"] = time.perf_counter() - t0; t0 = time.perf_counter()

    idx_path, val_path = _spill_paths(spill_prefix, block)
    _write_bytes(idx_path, 0, _as_bytes(cells), mode="wb")
    _write_bytes(val_path, 0, _as_bytes(bx), mode="wb")
    del cells, bx
    tm["spill"] = time.perf_counter() - t0; t0 = time.perf_counter()

    sum_b, cnt_b, sq_b = _block_stats(gi, dv, rel, state_block, n_genes, n_states)
    tm["stats"] = time.perf_counter() - t0
    return block, bp, sum_b, cnt_b, sq_b, e - s, tm


def _gather_range(g0: int, g1: int, out_ip: np.ndarray, segs, spill_prefix: str,
                  out_idx_path: str, out_val_path: str, expr_dtype: str):
    """Phase B for genes [g0, g1): copy every block's segments into the output range.

    ``out_ip`` is out_indptr[g0:g1+1]. ``segs`` lists (block, the block's bp[g0:g1+1])
    for every non-empty block, in block order. Writes the aligned middle of the range and
    returns its edges for the main process."""
    t0 = time.perf_counter()
    o0, o1 = int(out_ip[0]), int(out_ip[-1])
    want = np.diff(out_ip)
    before = np.zeros(g1 - g0, dtype=np.int64)
    edges: Dict[str, List[Tuple[int, bytes]]] = {"idx": [], "val": []}
    middle = {"idx": 0, "val": 0}
    if o1 > o0:
        dt = _out_dtype(expr_dtype)
        reg_i = np.empty(o1 - o0, dtype=np.uint32)
        reg_v = np.empty(o1 - o0, dtype=dt)
        gene_start = out_ip[:-1] - o0                  # where each gene begins in the range
        for block, bps in segs:
            cnt = np.diff(bps)
            n = int(bps[-1] - bps[0])
            if n:
                # entry j of this slice belongs to gene g; it lands after the earlier
                # blocks' entries of g, at its own offset inside g's segment.
                dest = np.repeat(gene_start + before - (bps[:-1] - bps[0]), cnt)
                dest += np.arange(n, dtype=np.int64)
                idx_path, val_path = _spill_paths(spill_prefix, block)
                reg_i[dest] = _read_at(idx_path, int(bps[0]), n, np.uint32)
                vals = _read_at(val_path, int(bps[0]), n, np.float32)
                reg_v[dest] = vals if dt is np.float32 else vals.astype(dt)
                del dest, vals
            before += cnt
        if np.array_equal(before, want):
            for key, path, arr in (("idx", out_idx_path, reg_i), ("val", out_val_path, reg_v)):
                start = _npy_data_offset(path) + o0 * arr.dtype.itemsize
                mv = _as_bytes(arr)
                edges[key] = _write_owned(path, start, mv)
                middle[key] = len(mv) - sum(len(b) for _, b in edges[key])
        del reg_i, reg_v
    if not np.array_equal(before, want):
        bad = int((before != want).sum())
        raise RuntimeError(f"genes {g0}-{g1}: {bad} genes received the wrong number of entries")
    return g0, g1, o1 - o0, edges, middle, time.perf_counter() - t0


def _gene_ranges(out_indptr: np.ndarray, n_ranges: int) -> List[Tuple[int, int]]:
    """Split the genes into up to ``n_ranges`` contiguous ranges of about equal nnz."""
    n_genes = out_indptr.size - 1
    nnz = int(out_indptr[-1])
    cuts = np.searchsorted(out_indptr, np.linspace(0, nnz, n_ranges + 1)[1:-1], side="left")
    bounds = np.unique(np.concatenate([[0], np.clip(cuts, 0, n_genes), [n_genes]]))
    return [(int(a), int(b)) for a, b in zip(bounds[:-1], bounds[1:]) if b > a]


def _run_single_scan(ex, src, blocks, block_nnz, full_indptr, state16, n_genes, n_states,
                     paths, expr_dtype, nnz, n_workers, spill_dir, log):
    prefix = getattr(paths, "prefix", None) or "bundle"
    spill_prefix = os.path.join(spill_dir, f".{prefix}_spill.{os.getpid()}")
    live = [i for i in range(len(blocks)) if block_nnz[i]]
    gs_len = n_genes * n_states
    futs: list = []
    try:
        log(f"spill: {len(live)} blocks x 2 files {spill_prefix}.b*.idx/.val "
            f"({nnz * 8 / 2**30:.2f} GiB, removed at the end)")

        # ---- phase A: read each block once --------------------------------------
        t0 = time.time()
        futs = [ex.submit(_scan_block, src, i, blocks[i][0], blocks[i][1],
                          full_indptr[blocks[i][0]:blocks[i][1] + 1].copy(),
                          state16[blocks[i][0]:blocks[i][1]].copy(), n_genes, n_states,
                          spill_prefix)
                for i in live]
        bp_by_block: Dict[int, np.ndarray] = {}
        sum_gs = np.zeros(gs_len, dtype=np.float64)
        cnt_gs = np.zeros(gs_len, dtype=np.int64)
        gene_sumsq = np.zeros(n_genes, dtype=np.float64)
        scanned = 0
        for fu in futs:                                # block order: the reference's order
            bi, bp, sum_b, cnt_b, sq_b, nnz_b, tm = fu.result()
            sum_gs += sum_b
            cnt_gs += cnt_b
            gene_sumsq += sq_b
            scanned += nnz_b
            bp_by_block[bi] = bp
            r0, r1 = blocks[bi]
            log(f"phase A block {bi + 1}/{len(blocks)} cells {r0:,}-{r1:,} "
                + " ".join(f"{k} {v:.1f}s" for k, v in tm.items()))
        if scanned != nnz:
            raise RuntimeError(f"phase A read {scanned} entries, expected {nnz}")
        log(f"phase A done in {time.time() - t0:.1f}s (source read once)")

        counts = np.zeros(n_genes, dtype=np.int64)
        for bp in bp_by_block.values():
            counts += np.diff(bp)
        if counts.sum() != nnz:
            raise RuntimeError(f"block counts sum to {counts.sum()}, expected {nnz}")
        log(f"genes with zero nnz: {int((counts == 0).sum()):,} of {n_genes:,}")
        out_indptr = np.zeros(n_genes + 1, dtype=np.int64)
        np.cumsum(counts, out=out_indptr[1:])

        # ---- allocate the gene-major store (same calls as the reference) ---------
        np.save(paths.expr_indptr, out_indptr)
        dt = _out_dtype(expr_dtype)
        out_idx = np.lib.format.open_memmap(paths.expr_indices, mode="w+", dtype=np.uint32, shape=(nnz,))
        out_val = np.lib.format.open_memmap(paths.expr_data, mode="w+", dtype=dt, shape=(nnz,))
        out_idx.flush(); out_val.flush()
        del out_idx, out_val
        log(f"allocated {paths.expr_indices} ({nnz * 4 / 2**30:.2f} GiB) and "
            f"{paths.expr_data} ({nnz * np.dtype(dt).itemsize / 2**30:.2f} GiB)")

        # ---- phase B: each gene range fills its own contiguous output -----------
        t0 = time.time()
        ranges = _gene_ranges(out_indptr, n_workers * _RANGES_PER_WORKER)
        futs = [ex.submit(_gather_range, g0, g1, out_indptr[g0:g1 + 1].copy(),
                          [(i, bp_by_block[i][g0:g1 + 1].copy()) for i in live],
                          spill_prefix, paths.expr_indices, paths.expr_data, expr_dtype)
                for g0, g1 in ranges]
        written = 0
        slowest = 0.0
        pending: Dict[str, List[Tuple[int, bytes]]] = {"idx": [], "val": []}
        wrote_bytes = {"idx": 0, "val": 0}
        for fu in futs:
            _g0, _g1, n, edges, middle, secs = fu.result()
            written += n
            slowest = max(slowest, secs)
            for key in pending:
                pending[key].extend(edges[key])
                wrote_bytes[key] += middle[key]
        # the edges, written by this process alone after every worker has finished
        n_edges = 0
        for key, path in (("idx", paths.expr_indices), ("val", paths.expr_data)):
            for off, buf in pending[key]:
                _write_bytes(path, off, memoryview(buf))
                wrote_bytes[key] += len(buf)
                n_edges += 1
        want_bytes = {"idx": nnz * 4, "val": nnz * np.dtype(dt).itemsize}
        if written != nnz or wrote_bytes != want_bytes:
            raise RuntimeError(f"phase B wrote {written} entries and {wrote_bytes} bytes, "
                               f"expected {nnz} and {want_bytes}")
        log(f"phase B done in {time.time() - t0:.1f}s: {len(ranges)} gene ranges, "
            f"slowest {slowest:.1f}s, {n_edges} unaligned edges written by the main process")
    except BaseException:
        for fu in futs:                                # queued tasks never start
            fu.cancel()
        raise
    finally:
        for i in live:
            for pth in _spill_paths(spill_prefix, i):
                if os.path.exists(pth):
                    os.remove(pth)
    return out_indptr, sum_gs, cnt_gs, gene_sumsq, written


# ------------------------------------------------------------------ driver


def _physical_ram() -> int:
    try:
        return int(os.sysconf("SC_PHYS_PAGES")) * int(os.sysconf("SC_PAGE_SIZE"))
    except (ValueError, OSError, AttributeError):
        return 0


def choose_workers(requested: Optional[int], max_block_nnz: int,
                   log: Callable[[str], None] = print) -> int:
    """Worker count: ``requested`` if given, else min(8, CPUs), capped so the workers'
    peak working set stays under half the physical RAM. Logs the reason for any cap."""
    cpus = os.cpu_count() or 1
    n = int(requested) if requested and requested > 0 else min(8, cpus)
    ram = _physical_ram()
    per = max(1, max_block_nnz * _BYTES_PER_NNZ_PEAK)
    if ram:
        cap = max(1, (ram // 2) // per)
        if cap < n:
            log(f"workers capped {n} -> {cap}: each needs about {per / 2**30:.1f} GiB, "
                f"RAM is {ram / 2**30:.0f} GiB")
            n = cap
    return n


def _executor(src, workers: int) -> Executor:
    if src.kind == "h5":
        return ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn"))
    return ThreadPoolExecutor(max_workers=workers)


def build_expression_store_parallel(
    src,
    n_cells: int,
    n_genes: int,
    state_code: np.ndarray,
    n_states: int,
    paths,
    expr_dtype: str,
    row_block: int,
    workers: Optional[int] = None,
    log: Callable[[str], None] = print,
    spill_dir: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Same contract and return value as ``precompute.build_expression_store``.

    ``src`` is an ``H5Source`` or an ``ArraySource``. ``spill_dir`` holds the temporary
    transposed blocks, 8 bytes per non-zero; it defaults to the bundle directory.
    Returns (sum_gs, cnt_gs, gene_sum, gene_sumsq, nnz).
    """
    t_all = time.time()
    if src.kind == "h5":
        import h5py
        with h5py.File(src.path, "r") as fh:
            full_indptr = fh[src.group]["indptr"][:].astype(np.int64)
    else:
        full_indptr = np.asarray(src.indptr).astype(np.int64)
    if full_indptr.size != n_cells + 1:
        raise RuntimeError(f"indptr has {full_indptr.size} entries, expected {n_cells + 1}")
    nnz = int(full_indptr[-1])
    log(f"layer nnz = {nnz:,}  density = {nnz / (n_cells * n_genes):.4f}")

    blocks: List[Tuple[int, int]] = [(r0, min(r0 + row_block, n_cells))
                                     for r0 in range(0, n_cells, row_block)]
    block_nnz = [int(full_indptr[r1] - full_indptr[r0]) for r0, r1 in blocks]
    n_workers = choose_workers(workers, max(block_nnz) if block_nnz else 1, log)
    log(f"parallel builder: {len(blocks)} blocks of {row_block} rows, "
        f"{n_workers} {'processes' if src.kind == 'h5' else 'threads'}, source {src.describe()}")
    state16 = np.asarray(state_code)

    with _executor(src, n_workers) as ex:
        out_indptr, sum_gs, cnt_gs, gene_sumsq, written = _run_single_scan(
            ex, src, blocks, block_nnz, full_indptr, state16, n_genes, n_states, paths,
            expr_dtype, nnz, n_workers, spill_dir or paths.bundle_dir, log)

    if out_indptr[-1] != nnz:
        raise RuntimeError("indptr tail does not equal nnz")
    if written != nnz:
        raise RuntimeError(f"wrote {written} entries, expected {nnz}")
    if cnt_gs.sum() != nnz:
        raise RuntimeError(f"statistics count {cnt_gs.sum()} != nnz {nnz}")
    log(f"store flushed, invariants held, {time.time() - t_all:.1f}s in the parallel builder")

    gene_sum = sum_gs.reshape(n_genes, n_states).sum(axis=1)
    return (sum_gs.reshape(n_genes, n_states), cnt_gs.reshape(n_genes, n_states),
            gene_sum, gene_sumsq, nnz)


def build_expression_store_from_csr(matrix, state_code: np.ndarray, n_states: int, paths,
                                    expr_dtype: str = "float32", row_block: int = 8192,
                                    workers: Optional[int] = None,
                                    log: Callable[[str], None] = print,
                                    spill_dir: Optional[str] = None):
    """Entry point for a caller that holds the layer in memory as a scipy CSR matrix.

    Writes the expression store into ``paths`` (a ``bundle.BundlePaths``) and returns the
    same tuple as ``build_expression_store_parallel``."""
    n_cells, n_genes = matrix.shape
    if len(state_code) != n_cells:
        raise ValueError(f"state_code has {len(state_code)} entries for {n_cells} cells")
    return build_expression_store_parallel(ArraySource.from_csr(matrix), n_cells, n_genes,
                                           state_code, n_states, paths, expr_dtype,
                                           row_block, workers, log, spill_dir)
