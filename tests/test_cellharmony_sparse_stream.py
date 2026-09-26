import io
import struct

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from altanalyze3.components.aggregate.pseudobulk_h5ad import (
    build_pseudobulk_from_adata,
    build_pseudobulk_h5ad,
)
from altanalyze3.components.cellHarmony import cellHarmony_lite
from altanalyze3.components.cellHarmony.sparse_stream import (
    read_sparse_stream,
    write_sparse_stream,
)


def _counts(n_obs=300, n_vars=120, seed=0):
    rng = np.random.default_rng(seed)
    dense = rng.poisson(0.4, size=(n_obs, n_vars)).astype(np.float32)
    dense[5] = 0  # one empty cell exercises an empty CSR row
    return sp.csr_matrix(dense)


def _stream_bytes(X, obs_names, var_names, obs=None):
    buf = io.BytesIO()
    write_sparse_stream(buf, X, obs_names, var_names, obs)
    return buf.getvalue()


def test_round_trip_is_exact():
    X = _counts()
    obs_names = [f"cell{i}" for i in range(X.shape[0])]
    var_names = [f"G{j}" for j in range(X.shape[1])]
    var_names[7] = var_names[8]  # duplicates pass through; cellHarmony makes them unique
    obs = pd.DataFrame({"Library": [f"L{i % 3}" for i in range(X.shape[0])],
                        "note": [None if i % 50 == 0 else f"n{i}" for i in range(X.shape[0])]})
    got = read_sparse_stream(io.BytesIO(_stream_bytes(X, obs_names, var_names, obs)), log=lambda *_: None)
    assert got.X.dtype == np.float32
    assert np.array_equal(got.X.indptr, X.indptr)
    assert np.array_equal(got.X.indices, X.indices)
    assert np.array_equal(got.X.data, X.data)
    assert list(got.obs_names) == obs_names
    assert list(got.var_names) == var_names
    assert list(got.obs["Library"]) == list(obs["Library"])
    assert got.obs["note"].isna().sum() == 6
    assert list(got.obs["note"].dropna()) == list(obs["note"].dropna())


def test_row_total_above_float32_exact_range():
    n = 2**24 + 1001  # float32 accumulation stalls at 2**24; the check must stay exact above it
    X = sp.csr_matrix((np.ones(n, dtype=np.float32), (np.zeros(n, dtype=np.int64), np.arange(n))), shape=(1, n))
    got = read_sparse_stream(io.BytesIO(_stream_bytes(X, ["deep"], [f"g{j}" for j in range(n)])), log=lambda *_: None)
    assert got.X.nnz == n


def test_truncated_stream_raises():
    X = _counts()
    raw = _stream_bytes(X, [f"c{i}" for i in range(X.shape[0])], [f"g{j}" for j in range(X.shape[1])])
    with pytest.raises(EOFError):
        read_sparse_stream(io.BytesIO(raw[:-100]), log=lambda *_: None)


def test_wrong_row_total_raises():
    X = _counts()
    raw = bytearray(_stream_bytes(X, [f"c{i}" for i in range(X.shape[0])], [f"g{j}" for j in range(X.shape[1])]))
    last_total = len(raw) - 8 - 8  # the final float64 row total sits before the trailer
    raw[last_total:last_total + 8] = struct.pack("<d", 1e9)
    with pytest.raises(ValueError, match="do not sum"):
        read_sparse_stream(io.BytesIO(bytes(raw)), log=lambda *_: None)


def test_bad_trailer_and_extra_bytes_raise():
    X = _counts()
    raw = _stream_bytes(X, [f"c{i}" for i in range(X.shape[0])], [f"g{j}" for j in range(X.shape[1])])
    with pytest.raises(ValueError, match="trailer"):
        read_sparse_stream(io.BytesIO(raw[:-8] + b"XXXXXXXX"), log=lambda *_: None)
    with pytest.raises(ValueError, match="after its trailer"):
        read_sparse_stream(io.BytesIO(raw + b"\0"), log=lambda *_: None)


def test_pseudobulk_extra_layer_and_file_path_agree(tmp_path):
    X = _counts()
    rng = np.random.default_rng(1)
    obs = pd.DataFrame({"state": rng.choice(["A", "B", "C"], X.shape[0]),
                        "sample": rng.choice(["s1", "s2"], X.shape[0])},
                       index=[f"c{i}" for i in range(X.shape[0])])
    adata = ad.AnnData(X=X.copy(), obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(X.shape[1])]))
    adata.layers["counts"] = X.copy()
    adata.layers["raw"] = X.multiply(2).tocsr()
    in_memory = build_pseudobulk_from_adata(adata, cluster_col="state", sample_col="sample",
                                            output_h5ad=tmp_path / "mem.h5ad", min_cells=10,
                                            extra_layers=["raw"])
    path = tmp_path / "cells.h5ad"
    adata.write_h5ad(path)
    from_file = build_pseudobulk_h5ad(path, cluster_col="state", sample_col="sample",
                                      output_h5ad=tmp_path / "file.h5ad", min_cells=10,
                                      extra_layers=["raw"])
    assert np.array_equal(in_memory.layers["counts"].toarray(), from_file.layers["counts"].toarray())
    assert np.allclose(in_memory.layers["raw"].toarray(), 2 * in_memory.layers["counts"].toarray())
    assert in_memory.layers["counts"].sum() == X.sum()
    assert in_memory.uns["pseudobulk"]["extra_layers"] == ["raw"]


def _synthetic_reference(tmp_path, n_genes=60, n_states=3, seed=2):
    rng = np.random.default_rng(seed)
    genes = [f"G{j}" for j in range(n_genes)]
    profiles = rng.gamma(0.5, 2.0, size=(n_states, n_genes))
    for s in range(n_states):
        profiles[s, s * 20:(s + 1) * 20] += 20.0
    ref = pd.DataFrame(np.log1p(profiles).T, index=genes, columns=[f"State{s}" for s in range(n_states)])
    ref_path = tmp_path / "SynthRef.txt"
    ref.to_csv(ref_path, sep="\t", index_label="UID")
    cells, labels = [], []
    for i in range(240):
        s = i % n_states
        cells.append(rng.poisson(profiles[s] * 2.0))
        labels.append(f"State{s}")
    X = sp.csr_matrix(np.asarray(cells, dtype=np.float32))
    return ref_path, X, genes, labels


def test_in_memory_adata_matches_h5ad_input(tmp_path):
    ref_path, X, genes, _ = _synthetic_reference(tmp_path)
    obs = pd.DataFrame({"Library": ["L1"] * X.shape[0]}, index=[f"c{i}" for i in range(X.shape[0])])
    h5ad_path = tmp_path / "query.h5ad"
    ad.AnnData(X=X.copy(), obs=obs.copy(), var=pd.DataFrame(index=genes)).write_h5ad(h5ad_path)
    kwargs = dict(h5_files=[], cellharmony_ref=str(ref_path), min_genes=1, min_counts=1,
                  mit_percent=None, alignment_mode="cosine")
    from_file = cellHarmony_lite.combine_and_align_h5(h5ad_file=str(h5ad_path),
                                                      output_dir=str(tmp_path / "file"), **kwargs)
    streamed = read_sparse_stream(io.BytesIO(_stream_bytes(X, obs.index, genes, obs)), log=lambda *_: None)
    from_stream = cellHarmony_lite.combine_and_align_h5(adata=streamed,
                                                        output_dir=str(tmp_path / "stream"), **kwargs)
    pd.testing.assert_frame_equal(from_file.reset_index(drop=True), from_stream.reset_index(drop=True))
    with pytest.raises(ValueError):
        cellHarmony_lite.combine_and_align_h5(adata=streamed, h5ad_file=str(h5ad_path),
                                              output_dir=str(tmp_path / "both"), **kwargs)
