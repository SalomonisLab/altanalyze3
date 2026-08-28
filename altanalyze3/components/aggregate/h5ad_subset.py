#!/usr/bin/env python3
"""Subset an h5ad to the cells whose obs column matches one or more values.

The reader never loads the full matrix. It reads obs and var with the AnnData element
reader, then copies only the selected CSR rows of ``X`` and of every layer, in chunks.
A 5 GB, 233,035-cell object therefore yields a 7,578-cell object inside a few hundred MB.

Entry point::

    python -m altanalyze3.components.aggregate.h5ad_subset \
        --h5ad <input.h5ad> --obs-key Library --obs-value Thymus --out <output.h5ad>

Every run prints the retention arithmetic and raises when an invariant fails.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def _read_elem():
    """Return the AnnData element reader available in this environment."""
    try:
        from anndata.experimental import read_elem
        return read_elem
    except Exception:
        from anndata._io.specs import read_elem
        return read_elem


def _csr_rows(group, rows: np.ndarray, n_vars: int, chunk_rows: int = 20000):
    """Copy the named CSR rows out of an on-disk h5ad matrix group."""
    from scipy.sparse import csr_matrix

    indptr = group["indptr"][:]
    lengths = (indptr[rows + 1] - indptr[rows]).astype(np.int64)
    total = int(lengths.sum())
    data = np.empty(total, dtype=group["data"].dtype)
    indices = np.empty(total, dtype=group["indices"].dtype)
    out = 0
    for start in range(0, len(rows), chunk_rows):
        for row in rows[start:start + chunk_rows]:
            lo, hi = int(indptr[row]), int(indptr[row + 1])
            if hi > lo:
                data[out:out + (hi - lo)] = group["data"][lo:hi]
                indices[out:out + (hi - lo)] = group["indices"][lo:hi]
            out += hi - lo
    if out != total:
        raise ValueError(f"copied {out} entries, expected {total}")
    new_indptr = np.zeros(len(rows) + 1, dtype=np.int64)
    np.cumsum(lengths, out=new_indptr[1:])
    matrix = csr_matrix((data, indices, new_indptr), shape=(len(rows), n_vars))
    if matrix.nnz != total:
        raise ValueError(f"subset nnz {matrix.nnz} differs from the summed row lengths {total}")
    return matrix


def _dense_rows(dataset, rows: np.ndarray) -> np.ndarray:
    """Copy the named rows out of an on-disk dense matrix."""
    return dataset[rows.tolist(), :]


def _matrix_rows(handle, name: Optional[str], rows: np.ndarray, n_vars: int):
    node = handle["X"] if name is None else handle["layers"][name]
    encoding = str(node.attrs.get("encoding-type", ""))
    if encoding == "csr_matrix":
        return _csr_rows(node, rows, n_vars)
    if encoding == "csc_matrix":
        raise ValueError(
            f"{'X' if name is None else name} is CSC. Rewrite it as CSR before subsetting; "
            "a column-oriented matrix has no cheap row slice."
        )
    return _dense_rows(node, rows)


def subset_h5ad(
    h5ad_path: str,
    out_path: str,
    obs_key: str,
    obs_values: Sequence[str],
    *,
    keep_layers: bool = True,
    compression: Optional[str] = "gzip",
    log=print,
) -> Tuple[int, int]:
    """Write the cells whose ``obs[obs_key]`` is in ``obs_values``. Return (kept, total)."""
    import anndata as ad
    import h5py

    read_elem = _read_elem()
    wanted = [str(v) for v in obs_values]

    with h5py.File(h5ad_path, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
        n_obs, n_vars = int(obs.shape[0]), int(var.shape[0])
        if obs_key not in obs.columns:
            raise KeyError(
                f"obs column '{obs_key}' absent from {h5ad_path}. "
                f"Available: {', '.join(map(str, obs.columns))}"
            )
        values = obs[obs_key].astype(str)
        present = sorted(values.unique().tolist())
        missing = [v for v in wanted if v not in present]
        if missing:
            raise ValueError(
                f"obs['{obs_key}'] holds no value {missing}. Available: {present}"
            )
        mask = values.isin(wanted).to_numpy()
        rows = np.flatnonzero(mask)
        if rows.size == 0:
            raise ValueError(f"obs['{obs_key}'] matched no cell for {wanted}")
        log(f"[subset] {h5ad_path}")
        log(f"[subset] obs['{obs_key}'] in {wanted}: {rows.size} of {n_obs} cells "
            f"({100.0 * rows.size / n_obs:.2f}%), {n_vars} features kept")

        X = _matrix_rows(handle, None, rows, n_vars)
        layers = {}
        if keep_layers and "layers" in handle:
            for name in list(handle["layers"].keys()):
                layers[name] = _matrix_rows(handle, name, rows, n_vars)
                log(f"[subset] layer '{name}' copied")
        uns = read_elem(handle["uns"]) if "uns" in handle else {}
        obsm = {}
        if "obsm" in handle:
            for name in list(handle["obsm"].keys()):
                obsm[name] = read_elem(handle["obsm"][name])[rows]

    sub_obs = obs.iloc[rows].copy()
    for column in sub_obs.columns:
        if isinstance(sub_obs[column].dtype, pd.CategoricalDtype):
            sub_obs[column] = sub_obs[column].cat.remove_unused_categories()
    adata = ad.AnnData(X=X, obs=sub_obs, var=var, layers=layers, uns=uns, obsm=obsm)

    if adata.n_obs != rows.size:
        raise ValueError(f"wrote {adata.n_obs} cells, expected {rows.size}")
    if adata.n_vars != n_vars:
        raise ValueError(f"wrote {adata.n_vars} features, expected {n_vars}")
    if not adata.obs_names.equals(obs.index[rows]):
        raise ValueError("output barcode order differs from the selected input rows")
    for name, matrix in layers.items():
        if matrix.shape != (rows.size, n_vars):
            raise ValueError(f"layer '{name}' has shape {matrix.shape}, expected {(rows.size, n_vars)}")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    adata.write(out_path, compression=compression)
    log(f"[subset] wrote {adata.n_obs} x {adata.n_vars} -> {os.path.abspath(out_path)}")
    return rows.size, n_obs


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Subset an h5ad to the cells whose obs column matches one or more values."
    )
    parser.add_argument("--h5ad", required=True, help="Input h5ad.")
    parser.add_argument("--out", required=True, help="Output h5ad.")
    parser.add_argument("--obs-key", required=True, help="obs column to filter on, e.g. Library.")
    parser.add_argument("--obs-value", required=True, nargs="+",
                        help="One or more values of --obs-key to keep.")
    parser.add_argument("--drop-layers", action="store_true",
                        help="Do not copy layers. By default every layer is copied.")
    parser.add_argument("--compression", default="gzip", choices=["gzip", "lzf", "none"])
    args = parser.parse_args(argv)

    compression = None if args.compression == "none" else args.compression
    subset_h5ad(
        args.h5ad, args.out, args.obs_key, args.obs_value,
        keep_layers=not args.drop_layers, compression=compression,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
