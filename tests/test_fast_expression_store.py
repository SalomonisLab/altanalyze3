"""The parallel expression-store builder writes the reference builder's bytes.

``precompute.build_expression_store`` is the reference. ``fast_store`` must give the same
files and the same returned arrays for every source it accepts. The synthetic layers
cover: empty rows, a block with no non-zeros, genes with no non-zeros, unsorted genes
inside a row, a float64 source, integer counts, float16 output and several blocks.
"""
import hashlib
import os

import h5py
import numpy as np
import pytest
import scipy.sparse as sp

from altanalyze3.components.visualization.scalable_viewer import bundle as B
from altanalyze3.components.visualization.scalable_viewer import fast_store as F
from altanalyze3.components.visualization.scalable_viewer import precompute as P

FILES = ("expr_indptr", "expr_indices", "expr_data")


def _layer(n_cells, n_genes, density, dtype, seed, unsorted=False):
    rng = np.random.default_rng(seed)
    m = sp.random(n_cells, n_genes, density=density, format="csr", random_state=rng,
                  dtype=np.float64)
    m.data = (rng.gamma(2.0, 1.5, size=m.nnz) * 3).astype(dtype)
    dense = m.toarray()
    dense[5:9] = 0                                # empty rows
    dense[40:60] = 0                              # a whole block of empty rows (row_block 20)
    dense[:, 3] = 0                               # genes with no non-zeros
    dense[:, n_genes - 1] = 0
    m = sp.csr_matrix(dense.astype(dtype))
    if unsorted:                                  # reverse the gene order inside each row
        for r in range(n_cells):
            a, b = m.indptr[r], m.indptr[r + 1]
            m.indices[a:b] = m.indices[a:b][::-1].copy()
            m.data[a:b] = m.data[a:b][::-1].copy()
        m.has_sorted_indices = False
    return m


def _write_h5(path, m):
    with h5py.File(path, "w") as f:
        g = f.create_group("X")
        g.attrs["encoding-type"] = "csr_matrix"
        g.attrs["shape"] = np.asarray(m.shape)
        g.create_dataset("indptr", data=m.indptr.astype(np.int64))
        g.create_dataset("indices", data=m.indices.astype(np.int32), compression="gzip",
                         chunks=(97,))
        g.create_dataset("data", data=m.data, compression="gzip", chunks=(97,))


def _sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def _compare(tmp_path, m, expr_dtype, row_block, n_states=7):
    n_cells, n_genes = m.shape
    rng = np.random.default_rng(1)
    state_code = rng.integers(0, n_states, size=n_cells).astype(np.int16)
    h5 = str(tmp_path / "layer.h5")
    _write_h5(h5, m)

    ref = B.BundlePaths(str(tmp_path / "ref"), "T")
    os.makedirs(ref.bundle_dir)
    with h5py.File(h5, "r") as f:
        want = P.build_expression_store(f["X"], n_cells, n_genes, state_code, n_states,
                                        ref, expr_dtype, row_block)

    runs = {
        "h5_processes": (F.H5Source(h5, "X"), 2),
        "array_threads": (F.ArraySource(m.indptr, m.indices, m.data), 3),
        "array_serial": (F.ArraySource(m.indptr, m.indices, m.data), 1),
    }
    for name, (src, workers) in runs.items():
        out = B.BundlePaths(str(tmp_path / name), "T")
        os.makedirs(out.bundle_dir)
        got = F.build_expression_store_parallel(src, n_cells, n_genes, state_code, n_states,
                                                out, expr_dtype, row_block, workers=workers,
                                                log=lambda _m: None)
        assert sorted(os.listdir(out.bundle_dir)) == sorted(
            os.path.basename(getattr(out, k)) for k in FILES), "spill files left behind"
        for key in FILES:
            assert _sha(getattr(out, key)) == _sha(getattr(ref, key)), (name, key)
        for i, label in enumerate(("sum_gs", "cnt_gs", "gene_sum", "gene_sumsq")):
            assert got[i].dtype == want[i].dtype, (name, label)
            assert np.array_equal(got[i], want[i]), (name, label)
            assert got[i].tobytes() == want[i].tobytes(), (name, label)
        assert got[4] == want[4]


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int64])
def test_matches_reference_float32_store(tmp_path, dtype):
    _compare(tmp_path, _layer(157, 61, 0.2, dtype, seed=3), "float32", row_block=20)


def test_matches_reference_float16_store(tmp_path):
    _compare(tmp_path, _layer(157, 61, 0.2, np.float32, seed=4), "float16", row_block=20)


def test_matches_reference_unsorted_rows(tmp_path):
    _compare(tmp_path, _layer(157, 61, 0.2, np.float32, seed=5, unsorted=True), "float32",
             row_block=20)


def test_matches_reference_single_block(tmp_path):
    _compare(tmp_path, _layer(157, 61, 0.2, np.float32, seed=6), "float32", row_block=8192)


def test_matches_reference_many_stats_chunks(tmp_path, monkeypatch):
    # 13-entry chunks split every block into many statistics chunks. Threads share the
    # module, so the patched size reaches the workers; spawned processes would not see it.
    monkeypatch.setattr(F, "_STATS_CHUNK", 13)
    m = _layer(157, 61, 0.2, np.float64, seed=9)
    n_cells, n_genes = m.shape
    state_code = np.random.default_rng(3).integers(0, 6, size=n_cells).astype(np.int16)
    h5 = str(tmp_path / "layer.h5")
    _write_h5(h5, m)
    ref = B.BundlePaths(str(tmp_path / "ref"), "T"); os.makedirs(ref.bundle_dir)
    with h5py.File(h5, "r") as f:
        want = P.build_expression_store(f["X"], n_cells, n_genes, state_code, 6, ref,
                                        "float32", 40)
    out = B.BundlePaths(str(tmp_path / "par"), "T"); os.makedirs(out.bundle_dir)
    got = F.build_expression_store_parallel(F.ArraySource(m.indptr, m.indices, m.data),
                                            n_cells, n_genes, state_code, 6, out, "float32",
                                            40, workers=3, log=lambda _m: None)
    for i in range(4):
        assert got[i].tobytes() == want[i].tobytes(), i
    for key in FILES:
        assert _sha(getattr(out, key)) == _sha(getattr(ref, key)), key


def test_csr_entry_point(tmp_path):
    m = _layer(157, 61, 0.2, np.float32, seed=7)
    state_code = np.random.default_rng(2).integers(0, 5, size=m.shape[0]).astype(np.int16)
    a = B.BundlePaths(str(tmp_path / "a"), "T"); os.makedirs(a.bundle_dir)
    b = B.BundlePaths(str(tmp_path / "b"), "T"); os.makedirs(b.bundle_dir)
    F.build_expression_store_from_csr(m, state_code, 5, a, row_block=20, workers=2,
                                      log=lambda _m: None)
    F.build_expression_store_parallel(F.ArraySource(m.indptr, m.indices, m.data),
                                      m.shape[0], m.shape[1], state_code, 5, b, "float32",
                                      20, workers=1, log=lambda _m: None)
    for key in FILES:
        assert _sha(getattr(a, key)) == _sha(getattr(b, key)), key
    # the store is the transpose of the input
    ip = np.load(a.expr_indptr); ix = np.load(a.expr_indices); dx = np.load(a.expr_data)
    csc = sp.csc_matrix((dx, ix.astype(np.int64), ip), shape=m.shape)
    assert (csc != m.astype(np.float32)).nnz == 0


def test_spill_removed_when_build_fails(tmp_path, monkeypatch):
    m = _layer(60, 20, 0.3, np.float32, seed=8)
    out = B.BundlePaths(str(tmp_path / "x"), "T"); os.makedirs(out.bundle_dir)

    def boom(*_a, **_k):
        raise RuntimeError("injected phase B failure")
    monkeypatch.setattr(F, "_gather_range", boom)
    with pytest.raises(RuntimeError, match="injected phase B failure"):
        F.build_expression_store_parallel(F.ArraySource(m.indptr, m.indices, m.data), 60, 20,
                                          np.zeros(60, np.int16), 1, out, "float32", 20,
                                          workers=2, log=lambda _m: None)
    assert not [n for n in os.listdir(out.bundle_dir) if "spill" in n]


def test_gather_range_count_mismatch_raises(tmp_path):
    # A gene range whose blocks hold fewer entries than out_indptr promises must fail.
    prefix = str(tmp_path / "spill")
    si, sv = F._spill_paths(prefix, 0)
    np.arange(9, dtype=np.uint32).tofile(si)
    np.ones(9, dtype=np.float32).tofile(sv)
    out_i = str(tmp_path / "oi.npy"); out_v = str(tmp_path / "ov.npy")
    for pth, dt in ((out_i, np.uint32), (out_v, np.float32)):
        np.lib.format.open_memmap(pth, mode="w+", dtype=dt, shape=(10,)).flush()
    out_ip = np.array([0, 4, 10], dtype=np.int64)              # genes 0..1 promise 10 entries
    segs = [(0, np.array([0, 4, 9], dtype=np.int64))]           # the one block holds 9
    with pytest.raises(RuntimeError, match="wrong number of entries"):
        F._gather_range(0, 2, out_ip, segs, prefix, out_i, out_v, "float32")
    # nothing reached the output
    assert not np.load(out_i).any() and not np.load(out_v).any()


def test_write_owned_covers_each_byte_once(tmp_path, monkeypatch):
    monkeypatch.setattr(F, "_WRITE_ALIGN", 64)
    path = str(tmp_path / "f.bin")
    for start, length in ((0, 10), (5, 64), (60, 200), (128, 128), (100, 1000)):
        with open(path, "wb") as fh:
            fh.write(b"\0" * 2000)
        payload = bytes((i % 251) + 1 for i in range(length))
        edges = F._write_owned(path, start, memoryview(payload))
        with open(path, "rb") as fh:
            got = bytearray(fh.read())
        mid = [i for i in range(2000) if got[i]]
        for off, buf in edges:                   # edges are NOT written by the worker
            assert all(got[off + j] == 0 for j in range(len(buf)))
            got[off:off + len(buf)] = buf
        assert bytes(got[start:start + length]) == payload
        assert not any(got[:start]) and not any(got[start + length:])
        for i in mid:                            # the worker's bytes are aligned windows
            assert (i // 64) * 64 >= start and (i // 64 + 1) * 64 <= start + length


def test_matches_reference_with_many_unaligned_edges(tmp_path, monkeypatch):
    # 64-byte windows make nearly every range end in an edge the main process writes. The
    # patch reaches the thread runs; spawned processes keep the 1 MiB default, where these
    # small files make every range one edge.
    monkeypatch.setattr(F, "_WRITE_ALIGN", 64)
    _compare(tmp_path, _layer(157, 61, 0.2, np.float32, seed=10), "float32", row_block=20)


def test_gene_ranges_partition_all_genes():
    ip = np.concatenate([[0], np.cumsum(np.random.default_rng(0).integers(0, 50, size=997))])
    for n in (1, 3, 16, 5000):
        r = F._gene_ranges(ip.astype(np.int64), n)
        assert r[0][0] == 0 and r[-1][1] == 997
        assert all(a < b for a, b in r)
        assert all(r[i][1] == r[i + 1][0] for i in range(len(r) - 1))
