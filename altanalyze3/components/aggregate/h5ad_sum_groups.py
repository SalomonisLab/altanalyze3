#!/usr/bin/env python3
"""Sum chosen rows of an h5ad into named groups, and write the group x feature counts.

A combined pseudobulk h5ad keyed ``<state>.<library>-junction`` carries no obs columns, so the
modules that group on an obs column (``pseudobulk_h5ad``, ``h5ad_subset``) cannot select or combine
its rows. This module reads a two-column mapping, sums the listed rows of ``X`` into one row per
group, and leaves every unlisted row out::

    obs_name                    group
    HSC-1.Y1258-junction        AM75M
    MPP-1.Y5116-E1-junction     AM75M

Every mapped ``obs_name`` must exist in the source, and the summed total must equal the total of
the selected rows. The h5ad keeps every source feature. The optional text matrix keeps the features
with a non-zero count in at least one group, with ``UID`` as the first header field, which is the
junction-count layout SNAF reads through ``--juncounts``.

Entry point::

    python -m altanalyze3.components.aggregate.h5ad_sum_groups \\
        --h5ad <source.h5ad> --mapping <obs_to_group.tsv> --out-h5ad <groups.h5ad> \\
        [--out-txt <groups_counts.txt>]
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional

import numpy as np
import pandas as pd
from scipy import sparse

MAPPING_COLUMNS = ("obs_name", "group")


def _strings(dataset) -> np.ndarray:
    values = dataset[:]
    if values.dtype.kind in ("S", "O"):
        values = np.array([v.decode() if isinstance(v, bytes) else str(v) for v in values], dtype=object)
    return values


def _index(h5, key: str) -> np.ndarray:
    group = h5[key]
    name = group.attrs.get("_index", "_index")
    name = name.decode() if isinstance(name, bytes) else name
    return _strings(group[name])


def read_mapping(path: str) -> pd.DataFrame:
    """Read the obs_name -> group mapping and reject a blank, duplicated or incomplete row."""
    table = pd.read_csv(path, sep="\t", dtype=str, comment="#").dropna(how="all")
    if list(table.columns[:2]) != list(MAPPING_COLUMNS):
        raise ValueError(f"{path}: the first two columns must be {MAPPING_COLUMNS}, found {list(table.columns[:2])}")
    table = table[list(MAPPING_COLUMNS)].fillna("").apply(lambda s: s.astype(str).str.strip())
    blank = table.eq("").any(axis=1)
    if blank.any():
        raise ValueError(f"{path} holds {int(blank.sum())} row(s) with an empty field")
    duplicated = table["obs_name"].duplicated(keep=False)
    if duplicated.any():
        raise ValueError(f"{path} lists an obs_name more than once:\n"
                         f"{table.loc[duplicated].to_string(index=False)}")
    return table


def _read_rows(h5, rows: np.ndarray, n_var: int) -> sparse.csr_matrix:
    """Read the given obs rows of X (CSR group or dense dataset) into a float64 CSR matrix."""
    X = h5["X"]
    if hasattr(X, "keys"):
        encoding = X.attrs.get("encoding-type", b"")
        encoding = encoding.decode() if isinstance(encoding, bytes) else encoding
        if encoding != "csr_matrix":
            raise ValueError(f"X is stored as {encoding!r}; this module reads csr_matrix or a dense array")
        indptr = X["indptr"][:]
        data_parts, index_parts, ptr = [], [], [0]
        for r in rows:
            lo, hi = int(indptr[r]), int(indptr[r + 1])
            data_parts.append(X["data"][lo:hi].astype(np.float64))
            index_parts.append(X["indices"][lo:hi])
            ptr.append(ptr[-1] + hi - lo)
        data = np.concatenate(data_parts) if data_parts else np.zeros(0)
        indices = np.concatenate(index_parts) if index_parts else np.zeros(0, dtype=np.int64)
        return sparse.csr_matrix((data, indices, np.asarray(ptr)), shape=(len(rows), n_var))
    return sparse.csr_matrix(np.vstack([X[int(r), :] for r in rows]).astype(np.float64))


def sum_groups(
    h5ad_path: str,
    mapping_path: str,
    out_h5ad: str,
    *,
    out_txt: Optional[str] = None,
    compression: Optional[str] = "gzip",
    log=print,
) -> pd.DataFrame:
    """Write the group x feature sums and return the per-group table."""
    import anndata as ad
    import h5py

    mapping = read_mapping(mapping_path)
    with h5py.File(h5ad_path, "r") as h5:
        obs_names = _index(h5, "obs")
        var_names = _index(h5, "var")
        position = pd.Series(np.arange(len(obs_names)), index=pd.Index(obs_names))
        if position.index.duplicated().any():
            raise ValueError(f"{h5ad_path}: obs names are not unique")
        absent = [o for o in mapping["obs_name"] if o not in position.index]
        if absent:
            raise ValueError(f"{len(absent)} of {len(mapping)} mapped obs_name(s) are absent from "
                             f"{h5ad_path}; the first are {absent[:20]}")
        rows = position.loc[mapping["obs_name"]].to_numpy()
        log(f"[sum] {h5ad_path}: {len(obs_names)} obs x {len(var_names)} features; "
            f"summing {len(rows)} of {len(obs_names)} obs into {mapping['group'].nunique()} groups")
        selected = _read_rows(h5, rows, len(var_names))

    groups = list(dict.fromkeys(mapping["group"]))
    group_of = mapping["group"].map({g: i for i, g in enumerate(groups)}).to_numpy()
    indicator = sparse.csr_matrix((np.ones(len(rows)), (group_of, np.arange(len(rows)))),
                                  shape=(len(groups), len(rows)))
    summed = (indicator @ selected).tocsr()
    summed.sum_duplicates()
    summed.eliminate_zeros()

    selected_total, summed_total = float(selected.sum()), float(summed.sum())
    if not np.isclose(summed_total, selected_total, rtol=1e-12, atol=0.5):
        raise ValueError(f"conservation failed: groups sum to {summed_total}, selected rows to {selected_total}")
    row_totals = np.asarray(selected.sum(axis=1)).ravel()
    group_totals = np.asarray(summed.sum(axis=1)).ravel()
    for i, g in enumerate(groups):
        expected = float(row_totals[group_of == i].sum())
        if not np.isclose(group_totals[i], expected, rtol=1e-12, atol=0.5):
            raise ValueError(f"conservation failed for group {g}: {group_totals[i]} vs {expected}")
    integral = bool(np.all(np.mod(summed.data, 1) == 0))
    exact_float32 = integral and (summed.data.size == 0 or summed.data.max() < 2 ** 24)
    log(f"[sum] conservation: {summed_total:.0f} counts in {len(groups)} groups = "
        f"{selected_total:.0f} counts in the {len(rows)} selected obs; integer counts: {integral}")

    table = pd.DataFrame({
        "group": groups,
        "n_obs_summed": [int((group_of == i).sum()) for i in range(len(groups))],
        "total_counts": group_totals,
        "nonzero_features": np.diff(summed.indptr),
        "obs_names": [";".join(mapping["obs_name"][group_of == i]) for i in range(len(groups))],
    })
    for r in table.itertuples():
        log(f"[sum]     {r.group:<24s} {r.n_obs_summed:>4d} obs  {r.total_counts:>14.0f} counts  "
            f"{r.nonzero_features:>10d} non-zero features")

    obs = table.drop(columns="group").set_index(pd.Index(groups, name="group"))
    out = ad.AnnData(
        X=summed.astype(np.float32 if exact_float32 else np.float64),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(var_names, name="feature")),
    )
    out.uns["h5ad_sum_groups"] = {"source_h5ad": os.path.abspath(h5ad_path),
                                  "mapping": os.path.abspath(mapping_path),
                                  "selected_obs": int(len(rows)), "source_obs": int(len(obs_names))}
    os.makedirs(os.path.dirname(os.path.abspath(out_h5ad)) or ".", exist_ok=True)
    out.write(out_h5ad, compression=compression)
    log(f"[sum] wrote {len(groups)} groups x {len(var_names)} features -> {os.path.abspath(out_h5ad)}")

    table_path = os.path.splitext(out_h5ad)[0] + "_groups.tsv"
    table.to_csv(table_path, sep="\t", index=False, float_format="%.6g")
    log(f"[sum] wrote the group table -> {os.path.abspath(table_path)}")

    if out_txt:
        keep = np.unique(summed.indices)
        dense = summed[:, keep].toarray().T
        frame = pd.DataFrame(dense.astype(np.int64) if integral else dense,
                             index=pd.Index(var_names[keep], name="UID"), columns=groups)
        os.makedirs(os.path.dirname(os.path.abspath(out_txt)) or ".", exist_ok=True)
        frame.to_csv(out_txt, sep="\t")
        log(f"[sum] wrote {len(keep)} of {len(var_names)} features (non-zero in at least one group; "
            f"{len(var_names) - len(keep)} all-zero features left out) x {len(groups)} groups "
            f"-> {os.path.abspath(out_txt)}")
    return table


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Sum the mapped obs rows of an h5ad into named groups.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--h5ad", required=True, help="Source h5ad; X must hold counts (CSR or dense).")
    parser.add_argument("--mapping", required=True, help="TSV with columns obs_name and group.")
    parser.add_argument("--out-h5ad", required=True, help="Output h5ad, one row per group, every source feature.")
    parser.add_argument("--out-txt", default=None,
                        help="Optional tab-delimited matrix (UID x groups) of the features non-zero in a group.")
    parser.add_argument("--compression", default="gzip", choices=["gzip", "lzf", "none"])
    args = parser.parse_args(argv)
    sum_groups(args.h5ad, args.mapping, args.out_h5ad, out_txt=args.out_txt,
               compression=None if args.compression == "none" else args.compression)
    return 0


if __name__ == "__main__":
    sys.exit(main())
