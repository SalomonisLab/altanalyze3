"""Random-partition metacells built by summing raw counts across many h5ad files.

A metacell is the sum of raw counts over a disjoint block of cells. The module shuffles
the cells of one (group, population) stratum, cuts consecutive blocks of `size` cells, and
keeps a final short block only when it holds at least `min_tail` cells. Every cell enters at
most one metacell, and the membership table records which cells each metacell holds.

Block rule, for a stratum of n cells:

    n = 6   -> 0 metacells      (below min_tail)
    n = 7   -> 1 metacell  of 7 cells
    n = 19  -> 1 metacell  of 19 cells
    n = 26  -> 1 metacell  of 20 cells, 6 cells dropped
    n = 27  -> 2 metacells of 20 and 7 cells

The module sums raw counts and refuses to sum anything else. `verify_raw_counts` rejects a
matrix that holds negative or fractional values, and rejects a depth-normalized matrix whose
per-cell totals are constant. A caller cannot silence that gate.

Different studies carry different feature sets, so the caller supplies one reference feature
list. Every study column maps onto its own reference column; the module never folds two study
columns onto one reference column and never drops a column in silence.

Memory stays bounded. The module scans each h5ad once in cell order and writes a metacell row
as soon as the scan passes its last member cell. It streams the result into the output h5ad
one block at a time, so it never concatenates the whole matrix.

Entry point:

    python -m altanalyze3.components.metacells.metacell_random --cells CELLS.tsv \\
        --reference-features FEATURES.tsv --output OUT.h5ad --membership MEMBERS.tsv.gz

`build_metacells` is the same code path for a caller that already holds the cell table in
memory.
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import re
import sys
import time
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

LOGGER = logging.getLogger("metacell_random")

DEFAULT_SIZE = 20
DEFAULT_MIN_TAIL = 7
DEFAULT_SEED = 0
CHUNK_CELLS = 40000

#: Layers searched for raw counts, in order, before `raw/X` and `X`.
COUNT_LAYER_CANDIDATES = ("counts", "raw_counts", "raw", "umi", "UMI", "X_raw")

REQUIRED_CELL_COLUMNS = ("dataset", "h5ad", "CellBarcode", "group", "population")


###############################################################################
# Block rule
###############################################################################


def count_blocks(n: int, size: int = DEFAULT_SIZE, min_tail: int = DEFAULT_MIN_TAIL) -> int:
    """Return how many metacells a stratum of `n` cells produces."""
    if size <= 0:
        raise ValueError("size must be a positive integer")
    if not 0 < min_tail <= size:
        raise ValueError("min_tail must satisfy 0 < min_tail <= size")
    if n < min_tail:
        return 0
    full, tail = divmod(n, size)
    return full + (1 if tail >= min_tail else 0)


def partition_indices(
    indices: np.ndarray,
    rng: np.random.Generator,
    size: int = DEFAULT_SIZE,
    min_tail: int = DEFAULT_MIN_TAIL,
) -> List[np.ndarray]:
    """Shuffle `indices` and cut them into consecutive blocks under the block rule.

    Each returned block holds `size` cells, except a final block that holds between
    `min_tail` and `size - 1` cells. No index appears in two blocks.
    """
    indices = np.asarray(indices)
    n_blocks = count_blocks(indices.size, size=size, min_tail=min_tail)
    if n_blocks == 0:
        return []
    shuffled = rng.permutation(indices)
    blocks = [shuffled[start:start + size] for start in range(0, n_blocks * size, size)]
    kept = int(sum(block.size for block in blocks))
    if kept > indices.size:
        raise AssertionError("a metacell block ran past the end of the stratum")
    return blocks


def stratum_seed(base_seed: int, *keys: object) -> int:
    """Return a reproducible seed for one stratum, independent of iteration order."""
    text = "\x1f".join(str(key) for key in keys)
    digest = 0
    for character in text.encode("utf-8"):
        digest = (digest * 1000003 + character) & 0xFFFFFFFFFFFF
    return (int(base_seed) * 0x9E3779B1 + digest) & 0x7FFFFFFF


###############################################################################
# Raw-count gate
###############################################################################


def _matrix_paths(handle: h5py.File, layer: Optional[str]) -> List[str]:
    names: List[str] = []
    if layer:
        names.append(f"layers/{layer}")
    if "layers" in handle:
        for candidate in COUNT_LAYER_CANDIDATES:
            path = f"layers/{candidate}"
            if path not in names and candidate in handle["layers"]:
                names.append(path)
    if "raw" in handle and "raw/X" in handle:
        names.append("raw/X")
    names.append("X")
    return [name for name in names if name in handle]


def _probe_rows(handle: h5py.File, path: str, n_probe: int) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Return probed per-cell totals and the probed value vectors."""
    node = handle[path]
    if not isinstance(node, h5py.Group):
        raise TypeError(f"{path} is dense; this module reads CSR sparse matrices")
    indptr = node["indptr"][:]
    n_cells = len(indptr) - 1
    rows = np.unique(np.linspace(0, n_cells - 1, min(n_probe, n_cells)).astype(np.int64))
    totals, vectors = [], []
    for row in rows:
        start, stop = indptr[row], indptr[row + 1]
        values = np.asarray(node["data"][start:stop], dtype=np.float64)
        totals.append(values.sum())
        vectors.append(values)
    return np.asarray(totals), vectors


def verify_raw_counts(
    path: str,
    layer: Optional[str] = None,
    n_probe: int = 400,
) -> Dict[str, object]:
    """Find the raw-count matrix in `path` and prove it holds raw counts.

    The function walks the candidate matrices in order and returns the first that passes
    every test. It raises when none passes, and the message names each candidate and the
    test it failed.

    Tests, all of which must pass:

    1. every probed value is a non-negative whole number;
    2. probed per-cell totals are not constant, which rules out a CP10k or CPM matrix;
    3. probed per-cell totals are not all 1, which rules out a fraction matrix.

    The result also reports whether `X` reconstructs as ``log1p(counts / total * 1e4)``.
    A match proves the candidate is the matrix `X` was derived from. A mismatch is
    reported, not raised, because a file may hold a differently normalized `X`.
    """
    failures: List[str] = []
    with h5py.File(path, "r") as handle:
        candidates = _matrix_paths(handle, layer)
        if not candidates:
            raise ValueError(f"{path} holds no readable matrix")
        for candidate in candidates:
            try:
                totals, vectors = _probe_rows(handle, candidate, n_probe)
            except TypeError as error:
                failures.append(f"{candidate}: {error}")
                continue
            values = np.concatenate(vectors) if vectors else np.zeros(0)
            finite = values[np.isfinite(values)]
            if finite.size == 0:
                failures.append(f"{candidate}: no finite values in the probe")
                continue
            if finite.min() < 0:
                failures.append(f"{candidate}: holds negative values (min {finite.min():.4g})")
                continue
            fraction_whole = float(np.mean(finite == np.rint(finite)))
            if fraction_whole < 1.0:
                failures.append(
                    f"{candidate}: {(1 - fraction_whole) * 100:.2f}% of probed values are "
                    f"fractional (max {finite.max():.4g}); this is a scaled matrix")
                continue
            positive = totals[totals > 0]
            if positive.size == 0:
                failures.append(f"{candidate}: every probed cell totals zero")
                continue
            spread = float(positive.std() / positive.mean())
            if spread < 0.01:
                failures.append(
                    f"{candidate}: per-cell totals are constant at {positive.mean():.1f} "
                    f"(CV {spread:.4f}); this is a depth-normalized matrix")
                continue
            deviation = None
            if "X" in handle and isinstance(handle["X"], h5py.Group) and candidate != "X":
                deviation = _reconstruction_deviation(handle, candidate, n_probe)
            report = dict(
                path=path, matrix=candidate, n_probed_cells=int(totals.size),
                total_min=float(positive.min()), total_median=float(np.median(positive)),
                total_max=float(positive.max()), total_cv=spread,
                value_max=float(finite.max()), fraction_whole=fraction_whole,
                x_is_log1p_cp10k_of_matrix=deviation,
                rejected=failures,
            )
            LOGGER.info(
                "[counts] %s -> %s | totals median %.0f CV %.2f | max value %.0f%s",
                os.path.basename(path), candidate, report["total_median"], spread,
                report["value_max"],
                "" if deviation is None else f" | X == log1p(CP10k) to {deviation:.3g}")
            return report
    raise ValueError(
        f"{path} holds no raw-count matrix. Candidates tried:\n  " + "\n  ".join(failures))


def _reconstruction_deviation(handle: h5py.File, candidate: str, n_probe: int) -> Optional[float]:
    """Return max |X - log1p(candidate / rowsum * 1e4)| over probed cells, or None."""
    node, x_node = handle[candidate], handle["X"]
    if node["indptr"].shape != x_node["indptr"].shape:
        return None
    indptr = node["indptr"][:]
    if not np.array_equal(indptr, x_node["indptr"][:]):
        return None
    n_cells = len(indptr) - 1
    rows = np.unique(np.linspace(0, n_cells - 1, min(n_probe, n_cells)).astype(np.int64))
    worst = 0.0
    for row in rows:
        start, stop = indptr[row], indptr[row + 1]
        if stop <= start:
            continue
        counts = np.asarray(node["data"][start:stop], dtype=np.float64)
        observed = np.asarray(x_node["data"][start:stop], dtype=np.float64)
        total = counts.sum()
        if total <= 0:
            continue
        worst = max(worst, float(np.abs(observed - np.log1p(counts / total * 1e4)).max()))
    return worst


###############################################################################
# Feature mapping
###############################################################################


def build_feature_map(
    study_genes: Sequence[str],
    reference_genes: Sequence[str],
    ensembl_by_reference: Optional[Dict[str, str]] = None,
    ensembl_by_study: Optional[Sequence[str]] = None,
    dataset: str = "",
) -> Tuple[np.ndarray, List[Tuple[str, str]]]:
    """Map every study column onto its own reference column.

    A study name that is absent from the reference is retried once with the scanpy
    ``var_names_make_unique`` suffix ``NAME-N`` rewritten to the R ``make.unique`` form
    ``NAME.N``. When both files carry Ensembl IDs the rename must agree on the ID, or the
    function raises. The function raises when a column cannot be mapped and when two study
    columns would land on one reference column.
    """
    reference_index = {gene: position for position, gene in enumerate(reference_genes)}
    mapped: List[Optional[str]] = []
    renamed: List[Tuple[str, str]] = []
    for position, gene in enumerate(study_genes):
        if gene in reference_index:
            mapped.append(gene)
            continue
        match = re.match(r"^(.*)-(\d+)$", gene)
        alternate = f"{match.group(1)}.{match.group(2)}" if match else None
        if alternate and alternate in reference_index:
            if ensembl_by_study is not None and ensembl_by_reference is not None:
                source = str(ensembl_by_study[position])
                target = ensembl_by_reference.get(alternate, "")
                if source and target and source != target:
                    raise ValueError(
                        f"{dataset}: refusing to rename {gene} to {alternate}; Ensembl IDs "
                        f"disagree ({source} against {target})")
            mapped.append(alternate)
            renamed.append((gene, alternate))
            continue
        mapped.append(None)
    unmapped = [gene for gene, target in zip(study_genes, mapped) if target is None]
    if unmapped:
        raise ValueError(
            f"{dataset}: {len(unmapped)} of {len(study_genes)} features are absent from the "
            f"reference and cannot be renamed: {unmapped[:10]}")
    columns = np.array([reference_index[target] for target in mapped], dtype=np.int64)
    if np.unique(columns).size != columns.size:
        collisions = pd.Series(columns).value_counts()
        worst = collisions[collisions > 1].index[:5]
        names = [reference_genes[position] for position in worst]
        raise ValueError(
            f"{dataset}: two study columns map onto one reference column ({names}); "
            f"summing them would corrupt the counts")
    return columns, renamed


###############################################################################
# Per-dataset streaming
###############################################################################


def _read_index(handle: h5py.File, group: str) -> List[str]:
    key = handle[group].attrs.get("_index", "_index")
    key = key.decode() if isinstance(key, bytes) else key
    return [value.decode() if isinstance(value, bytes) else str(value)
            for value in handle[group][key][:]]


def _read_var_column(handle: h5py.File, name: str) -> Optional[np.ndarray]:
    if name not in handle["var"]:
        return None
    node = handle["var"][name]
    if isinstance(node, h5py.Group) and "categories" in node:
        categories = [value.decode() if isinstance(value, bytes) else str(value)
                      for value in node["categories"][:]]
        return np.asarray(categories, dtype=object)[node["codes"][:]]
    return np.asarray([value.decode() if isinstance(value, bytes) else str(value)
                       for value in node[:]], dtype=object)


def plan_dataset(
    cells: pd.DataFrame,
    barcodes: Sequence[str],
    size: int,
    min_tail: int,
    seed: int,
) -> Tuple[List[np.ndarray], pd.DataFrame]:
    """Return the metacell blocks for one dataset, ordered by their last member cell.

    Ordering the blocks by their last member cell lets the scan write a metacell row as soon
    as it passes that cell, which is what keeps memory bounded.
    """
    position_of = {barcode: position for position, barcode in enumerate(barcodes)}
    missing = [barcode for barcode in cells["CellBarcode"] if barcode not in position_of]
    if missing:
        raise ValueError(
            f"{len(missing)} cell barcodes are absent from the h5ad, for example {missing[:5]}")
    frame = cells.copy()
    frame["_position"] = frame["CellBarcode"].map(position_of).to_numpy(dtype=np.int64)

    blocks: List[np.ndarray] = []
    records: List[Dict[str, object]] = []
    for (group, population), stratum in frame.groupby(["group", "population"], sort=True):
        positions = np.sort(stratum["_position"].to_numpy())
        rng = np.random.default_rng(stratum_seed(seed, group, population))
        for ordinal, block in enumerate(
                partition_indices(positions, rng, size=size, min_tail=min_tail), start=1):
            blocks.append(np.sort(block))
            records.append(dict(group=group, population=population, n_cells=int(block.size),
                                mc_ordinal=ordinal, mc_last=int(block.max())))
    if not blocks:
        return [], pd.DataFrame(columns=["group", "population", "n_cells"])

    obs = pd.DataFrame(records)
    order = np.argsort(obs["mc_last"].to_numpy(), kind="stable")
    blocks = [blocks[position] for position in order]
    obs = obs.iloc[order].reset_index(drop=True)
    used = np.concatenate(blocks)
    if np.unique(used).size != used.size:
        raise AssertionError("a cell was placed in two metacells")
    return blocks, obs


def stream_dataset_blocks(
    source: "str | h5py.File",
    blocks: Sequence[np.ndarray],
    matrix_path: str,
    column_map: np.ndarray,
    n_reference: int,
    chunk_cells: int = CHUNK_CELLS,
) -> Iterable[sp.csr_matrix]:
    """Scan one h5ad once and yield the metacell rows as CSR blocks, in row order.

    The scan reads cells in file order. It holds only the metacell rows whose member cells it
    has begun but not finished, so the resident buffer stays near one library's worth of rows.

    `source` accepts a path or an already open `h5py.File`. Pass an open handle when the file
    may be renamed or re-synced while the scan runs, as happens under a syncing file store; the
    handle follows the inode, so the scan finishes on the bytes it opened.
    """
    if not blocks:
        return
    n_rows = len(blocks)
    row_of_cell: Dict[int, int] = {}
    for row, block in enumerate(blocks):
        for position in block:
            row_of_cell[int(position)] = row
    last_cell = np.array([int(block.max()) for block in blocks], dtype=np.int64)
    if np.any(np.diff(last_cell) < 0):
        raise AssertionError("metacell rows are not ordered by their last member cell")

    cell_positions = np.fromiter(row_of_cell.keys(), dtype=np.int64, count=len(row_of_cell))
    cell_rows = np.fromiter(row_of_cell.values(), dtype=np.int64, count=len(row_of_cell))
    order = np.argsort(cell_positions)
    cell_positions, cell_rows = cell_positions[order], cell_rows[order]

    if n_reference > np.iinfo(np.int32).max or n_rows > np.iinfo(np.int32).max:
        raise ValueError("this module indexes rows and columns with int32")
    narrow_map = column_map.astype(np.int32)

    emitted = 0
    buffer_rows: List[np.ndarray] = []
    buffer_columns: List[np.ndarray] = []
    buffer_values: List[np.ndarray] = []

    opened = h5py.File(source, "r") if isinstance(source, str) else None
    handle = opened if opened is not None else source
    try:
        node = handle[matrix_path]
        indptr = node["indptr"][:]
        n_cells = len(indptr) - 1
        data_set, index_set = node["data"], node["indices"]
        for start in range(0, n_cells, chunk_cells):
            stop = min(start + chunk_cells, n_cells)
            low, high = int(indptr[start]), int(indptr[stop])
            if high > low:
                selection = slice(
                    np.searchsorted(cell_positions, start, "left"),
                    np.searchsorted(cell_positions, stop, "left"))
                chunk_cell_positions = cell_positions[selection]
                if chunk_cell_positions.size:
                    raw = data_set[low:high]
                    if np.issubdtype(raw.dtype, np.floating):
                        if not np.all(raw == np.rint(raw)):
                            raise AssertionError(
                                "the source matrix holds a fractional count")
                        values = np.rint(raw).astype(np.int32)
                    else:
                        values = np.asarray(raw, dtype=np.int32)
                    del raw
                    columns = narrow_map[np.asarray(index_set[low:high], dtype=np.int64)]
                    local = indptr[start:stop + 1] - low
                    row_lengths = np.diff(local)
                    row_for_value = np.full(high - low, -1, dtype=np.int32)
                    starts = local[:-1]
                    for position, row in zip(chunk_cell_positions, cell_rows[selection]):
                        offset = int(position) - start
                        begin = int(starts[offset])
                        row_for_value[begin:begin + int(row_lengths[offset])] = row
                    keep = row_for_value >= 0
                    if keep.any():
                        buffer_rows.append(row_for_value[keep])
                        buffer_columns.append(columns[keep])
                        buffer_values.append(values[keep])
                    del values, columns, row_for_value, keep
            ready = int(np.searchsorted(last_cell, stop - 1, "right"))
            if ready > emitted and buffer_rows:
                rows = np.concatenate(buffer_rows)
                cols = np.concatenate(buffer_columns)
                vals = np.concatenate(buffer_values)
                flush = rows < ready
                block = sp.coo_matrix(
                    (vals[flush], (rows[flush] - emitted, cols[flush])),
                    shape=(ready - emitted, n_reference)).tocsr()
                block.sum_duplicates()
                yield block
                held = ~flush
                buffer_rows = [rows[held]]
                buffer_columns = [cols[held]]
                buffer_values = [vals[held]]
                emitted = ready
    finally:
        if opened is not None:
            opened.close()
    if emitted != n_rows:
        raise AssertionError(f"emitted {emitted} of {n_rows} metacell rows")
    leftover = sum(part.size for part in buffer_rows)
    if leftover:
        raise AssertionError(f"{leftover} counts were never written to a metacell")


###############################################################################
# Reading an existing metacell matrix
###############################################################################


def stream_h5ad_row_blocks(
    source: "str | h5py.File",
    matrix_path: str,
    column_map: np.ndarray,
    n_reference: int,
    rows_per_block: int = 4000,
) -> Iterable[sp.csr_matrix]:
    """Yield the rows of an existing CSR matrix, remapped onto the reference features.

    The function reads whole rows, so it never changes a value. It moves each column to the
    position `column_map` names, which lets a matrix built on one feature list join a matrix
    built on another. `write_h5ad_csr` sorts the indices of each block it writes.
    """
    opened = h5py.File(source, "r") if isinstance(source, str) else None
    handle = opened if opened is not None else source
    try:
        node = handle[matrix_path]
        indptr = np.asarray(node["indptr"][:], dtype=np.int64)
        n_rows = len(indptr) - 1
        narrow = column_map.astype(np.int32)
        for start in range(0, n_rows, rows_per_block):
            stop = min(start + rows_per_block, n_rows)
            low, high = int(indptr[start]), int(indptr[stop])
            local = indptr[start:stop + 1] - low
            if high == low:
                yield sp.csr_matrix((stop - start, n_reference), dtype=np.int32)
                continue
            raw = node["data"][low:high]
            if np.issubdtype(raw.dtype, np.floating):
                if not np.all(raw == np.rint(raw)):
                    raise AssertionError(f"{matrix_path} holds a fractional count")
                values = np.rint(raw).astype(np.int32)
            else:
                values = np.asarray(raw, dtype=np.int32)
            columns = narrow[np.asarray(node["indices"][low:high], dtype=np.int64)]
            yield sp.csr_matrix((values, columns, local), shape=(stop - start, n_reference))
    finally:
        if opened is not None:
            opened.close()


###############################################################################
# Incremental h5ad writer
###############################################################################


#: Filters applied to every dataset written into the output h5ad.
STORE_FILTERS = dict(compression="gzip", compression_opts=9, shuffle=True)


def _write_dataframe(parent: h5py.Group, name: str, frame: pd.DataFrame, index_name: str) -> None:
    group = parent.create_group(name)
    group.attrs["encoding-type"] = "dataframe"
    group.attrs["encoding-version"] = "0.2.0"
    group.attrs["_index"] = index_name
    group.attrs["column-order"] = np.array(list(frame.columns), dtype=h5py.special_dtype(vlen=str))
    index = group.create_dataset(
        index_name, data=np.array([str(value) for value in frame.index], dtype=object),
        dtype=h5py.special_dtype(vlen=str), compression="gzip", compression_opts=9)
    index.attrs["encoding-type"] = "string-array"
    index.attrs["encoding-version"] = "0.2.0"
    for column in frame.columns:
        values = frame[column]
        if pd.api.types.is_numeric_dtype(values) and not pd.api.types.is_bool_dtype(values):
            node = group.create_dataset(column, data=values.to_numpy(), **STORE_FILTERS)
            node.attrs["encoding-type"] = "array"
            node.attrs["encoding-version"] = "0.2.0"
            continue
        categories = pd.Categorical(values.astype(str))
        node = group.create_group(column)
        node.attrs["encoding-type"] = "categorical"
        node.attrs["encoding-version"] = "0.2.0"
        node.attrs["ordered"] = False
        levels = node.create_dataset(
            "categories", data=np.array(list(categories.categories), dtype=object),
            dtype=h5py.special_dtype(vlen=str), compression="gzip", compression_opts=9)
        levels.attrs["encoding-type"] = "string-array"
        levels.attrs["encoding-version"] = "0.2.0"
        codes = node.create_dataset("codes", data=categories.codes.astype(np.int32),
                                    **STORE_FILTERS)
        codes.attrs["encoding-type"] = "array"
        codes.attrs["encoding-version"] = "0.2.0"


def _write_uns(parent: h5py.Group, payload: Dict[str, object]) -> None:
    group = parent.create_group("uns")
    group.attrs["encoding-type"] = "dict"
    group.attrs["encoding-version"] = "0.1.0"
    for key, value in payload.items():
        if isinstance(value, dict):
            # A nested mapping must stay structured. Stringifying it writes the literal text
            # "{'base': 2.0}", and a reader that calls .get on that string raises. The one
            # mapping this project stores is uns['log1p'], which cellHarmony_differential.py
            # reads at line 733 to learn the log base before it back-transforms X.
            child = group.create_group(key)
            child.attrs["encoding-type"] = "dict"
            child.attrs["encoding-version"] = "0.1.0"
            for inner_key, inner_value in value.items():
                if isinstance(inner_value, (int, float, np.integer, np.floating)):
                    leaf = child.create_dataset(inner_key, data=np.float64(inner_value))
                    leaf.attrs["encoding-type"] = "numeric-scalar"
                    leaf.attrs["encoding-version"] = "0.2.0"
                else:
                    leaf = child.create_dataset(inner_key, data=str(inner_value))
                    leaf.attrs["encoding-type"] = "string"
                    leaf.attrs["encoding-version"] = "0.2.0"
            continue
        if isinstance(value, (list, tuple)):
            node = group.create_dataset(
                key, data=np.array([str(item) for item in value], dtype=object),
                dtype=h5py.special_dtype(vlen=str), compression="gzip", compression_opts=9)
            node.attrs["encoding-type"] = "string-array"
            node.attrs["encoding-version"] = "0.2.0"
            continue
        node = group.create_dataset(key, data=str(value))
        node.attrs["encoding-type"] = "string"
        node.attrs["encoding-version"] = "0.2.0"


def write_h5ad_csr(
    path: str,
    blocks: Iterable[sp.csr_matrix],
    obs: pd.DataFrame,
    var: pd.DataFrame,
    uns: Dict[str, object],
    dtype: str = "int32",
    compression: str = "gzip",
    compression_level: int = 9,
    index_dtype: str = "int32",
) -> Dict[str, int]:
    """Write a CSR h5ad by appending each block, so the whole matrix never sits in memory.

    X stays sparse (CSR) and holds whole counts in `dtype`. Every dataset carries the HDF5
    byte-shuffle filter and gzip at `compression_level`, which shrinks integer count arrays
    well below the gzip default.
    """
    store = dict(compression=compression, compression_opts=compression_level, shuffle=True)
    n_rows, n_columns = obs.shape[0], var.shape[0]
    written_rows, written_values = 0, 0
    with h5py.File(path, "w") as handle:
        handle.attrs["encoding-type"] = "anndata"
        handle.attrs["encoding-version"] = "0.1.0"
        matrix = handle.create_group("X")
        matrix.attrs["encoding-type"] = "csr_matrix"
        matrix.attrs["encoding-version"] = "0.1.0"
        matrix.attrs["shape"] = np.array([n_rows, n_columns], dtype=np.int64)
        data = matrix.create_dataset("data", shape=(0,), maxshape=(None,), dtype=dtype,
                                     chunks=(1 << 20,), **store)
        # scipy's csr_row_index requires indices and indptr to share a dtype. An int32
        # indices with an int64 indptr raises "Output dtype not compatible with inputs" the
        # moment anything fancy-indexes rows, which is how rna2grn failed. So the two are
        # created from ONE parameter and can never diverge. Pass index_dtype="int64" for a
        # matrix holding more than 2,147,483,647 values: the metacell network reaches 5.46
        # billion, and int32 overflows at its row 57,168 of 145,411.
        if index_dtype not in ("int32", "int64"):
            raise AssertionError(f"index_dtype must be int32 or int64, not {index_dtype!r}")
        indices = matrix.create_dataset("indices", shape=(0,), maxshape=(None,),
                                        dtype=index_dtype, chunks=(1 << 20,), **store)
        indptr = matrix.create_dataset("indptr", shape=(n_rows + 1,), dtype=index_dtype,
                                       chunks=(min(n_rows + 1, 1 << 20),), **store)
        indptr[0] = 0
        for block in blocks:
            block = block.tocsr()
            # Several source files store explicit zeros, and a block of cells that are all
            # zero for a gene sums to a stored zero. A stored zero carries no information and
            # an absent entry means the same count, so drop it and keep the matrix sparse.
            block.eliminate_zeros()
            block.sort_indices()
            if block.shape[1] != n_columns:
                raise ValueError(f"block has {block.shape[1]} columns, expected {n_columns}")
            source = block.data
            if source.size:
                target = np.dtype(dtype)
                if np.issubdtype(target, np.integer):
                    # An integer store must receive whole numbers, and they must fit.
                    if not np.all(source == np.rint(source)):
                        raise AssertionError("a value is not a whole number")
                    limits = np.iinfo(target)
                    if int(source.min()) < limits.min or int(source.max()) > limits.max:
                        raise AssertionError(
                            f"a value reached {source.max()}, outside the {dtype} range")
                elif np.all(source == np.rint(source)) and float(source.max()) >= 2 ** 24:
                    # Whole numbers above 2**24 round in float32, so refuse them. A
                    # log-transformed value is not whole and needs no such check.
                    raise AssertionError(
                        f"a whole value reached {source.max():.0f}, at or above the 2**24 "
                        f"limit where {dtype} stops holding whole numbers exactly")
            values = source.astype(dtype)
            if (values.size and np.issubdtype(np.dtype(dtype), np.integer)
                    and not np.array_equal(values.astype(np.float64),
                                           source.astype(np.float64))):
                raise AssertionError(f"casting the values to {dtype} lost a value")
            data.resize((written_values + values.size,))
            data[written_values:] = values
            indices.resize((written_values + values.size,))
            indices[written_values:] = block.indices.astype(index_dtype)
            if (index_dtype == "int32"
                    and written_values + values.size > np.iinfo(np.int32).max):
                raise AssertionError(
                    f"the matrix reached {written_values + values.size} values, above the "
                    f"int32 limit of {np.iinfo(np.int32).max}. Pass index_dtype='int64', "
                    "which widens indices AND indptr together; widening one alone breaks "
                    "scipy's csr_row_index")
            indptr[written_rows + 1:written_rows + 1 + block.shape[0]] = (
                block.indptr[1:] + written_values)
            written_rows += block.shape[0]
            written_values += values.size
        if written_rows != n_rows:
            raise AssertionError(f"wrote {written_rows} rows, obs holds {n_rows}")
        _write_dataframe(handle, "obs", obs, "metacell_id")
        _write_dataframe(handle, "var", var, var.index.name or "gene")
        _write_uns(handle, uns)
        for name in ("obsm", "varm", "obsp", "varp", "layers"):
            group = handle.create_group(name)
            group.attrs["encoding-type"] = "dict"
            group.attrs["encoding-version"] = "0.1.0"
    return dict(n_obs=written_rows, n_vars=n_columns, nnz=written_values)


###############################################################################
# Driver
###############################################################################


def build_metacells(
    cells: pd.DataFrame,
    reference_features: Sequence[str],
    output: str,
    membership: Optional[str] = None,
    size: int = DEFAULT_SIZE,
    min_tail: int = DEFAULT_MIN_TAIL,
    seed: int = DEFAULT_SEED,
    layer: Optional[str] = None,
    obs_columns: Sequence[str] = (),
    population_column: str = "population",
    keep_group_columns: bool = True,
    var_extra: Optional[pd.DataFrame] = None,
    uns_extra: Optional[Dict[str, object]] = None,
    ensembl_by_reference: Optional[Dict[str, str]] = None,
    chunk_cells: int = CHUNK_CELLS,
) -> Dict[str, object]:
    """Build metacells across every dataset in `cells` and write one combined h5ad.

    `cells` needs the columns `dataset`, `h5ad`, `CellBarcode`, `group` and `population`.
    Columns named in `obs_columns` must hold one value for each `group`, and travel to the
    output `obs`. Every other column stays out of the object.
    """
    missing = [column for column in REQUIRED_CELL_COLUMNS if column not in cells.columns]
    if missing:
        raise ValueError(f"the cell table lacks the columns {missing}")
    reference_features = [str(gene) for gene in reference_features]
    if len(set(reference_features)) != len(reference_features):
        raise ValueError("the reference feature list holds a repeated name")
    n_reference = len(reference_features)

    group_values: Dict[str, Dict[str, str]] = {}
    for column in obs_columns:
        if column not in cells.columns:
            raise ValueError(f"obs column {column!r} is absent from the cell table")
        unique_per_group = cells.groupby("group")[column].nunique()
        offenders = unique_per_group[unique_per_group > 1]
        if len(offenders):
            raise ValueError(
                f"obs column {column!r} holds more than one value inside the groups "
                f"{list(offenders.index[:5])}; a metacell could not carry a single value")
        group_values[column] = cells.groupby("group")[column].first().to_dict()

    plans: List[Dict[str, object]] = []
    open_handles: List[h5py.File] = []
    obs_frames: List[pd.DataFrame] = []
    membership_rows: List[Tuple[str, str, str]] = []
    reports: List[Dict[str, object]] = []
    started = time.time()

    # Phase 1 plans every dataset without reading a count. Phase 2 then streams the counts
    # straight into the output file, so no dataset's matrix is ever held in memory.
    for dataset, table in cells.groupby("dataset", sort=True):
        paths = table["h5ad"].unique()
        if paths.size != 1:
            raise ValueError(f"dataset {dataset} names {paths.size} h5ad files")
        h5ad_path = str(paths[0])
        report = verify_raw_counts(h5ad_path, layer=layer)
        # Hold the handle open from here until the write finishes. A syncing file store may
        # rename or re-download the file mid-run; an open handle follows the inode, so the
        # scan still reads the bytes this plan was built from.
        handle = h5py.File(h5ad_path, "r")
        open_handles.append(handle)
        barcodes = _read_index(handle, "obs")
        study_genes = _read_index(handle, "var")
        study_ensembl = _read_var_column(handle, "ensg")
        if study_ensembl is None:
            study_ensembl = _read_var_column(handle, "ensembl_id")
        column_map, renamed = build_feature_map(
            study_genes, reference_features, ensembl_by_reference=ensembl_by_reference,
            ensembl_by_study=study_ensembl, dataset=str(dataset))
        if renamed:
            LOGGER.info("[%s] renamed %d suffixed features onto the reference: %s",
                        dataset, len(renamed), renamed[:3])

        blocks, obs = plan_dataset(table, barcodes, size=size, min_tail=min_tail, seed=seed)
        if not blocks:
            LOGGER.warning("[%s] no stratum reached %d cells; no metacell built",
                           dataset, min_tail)
            continue
        identifiers = [f"{row.group}|{row.population}|MC{row.mc_ordinal:02d}"
                       for row in obs.itertuples()]
        if len(set(identifiers)) != len(identifiers):
            raise AssertionError(f"{dataset}: metacell identifiers are not unique")
        n_strata = int(obs.groupby(["group", "population"]).ngroups)
        obs = obs.drop(columns=["mc_last", "mc_ordinal"])
        for column, lookup in group_values.items():
            obs[column] = obs["group"].map(lookup)
        if population_column != "population":
            obs = obs.rename(columns={"population": population_column})
        if not keep_group_columns:
            obs = obs.drop(columns=["group"])
        else:
            obs.insert(0, "dataset", str(dataset))
        obs.index = pd.Index(identifiers, name="metacell_id")

        for identifier, block in zip(identifiers, blocks):
            for position in block:
                membership_rows.append((identifier, str(dataset), barcodes[int(position)]))
        LOGGER.info("[%s] %d metacells from %d of %d cells (%.2f%%) over %d strata",
                    dataset, len(blocks), int(obs["n_cells"].sum()), len(table),
                    100 * int(obs["n_cells"].sum()) / len(table), n_strata)

        plans.append(dict(dataset=str(dataset), h5ad=h5ad_path, handle=handle,
                          matrix=str(report["matrix"]), column_map=column_map, blocks=blocks))
        obs_frames.append(obs)
        report["dataset"] = str(dataset)
        report["n_metacells"] = int(len(blocks))
        report["n_cells_used"] = int(obs["n_cells"].sum())
        report["n_cells_available"] = int(len(table))
        report["renamed_features"] = [f"{a}->{b}" for a, b in renamed]
        reports.append(report)

    if not obs_frames:
        raise ValueError("no dataset produced a metacell")
    obs = pd.concat(obs_frames)
    lead = [column for column in ("dataset", "group", population_column) if column in obs.columns]
    obs = obs[lead + [c for c in obs.columns if c not in lead + ["n_cells"]] + ["n_cells"]]
    if obs.index.has_duplicates:
        raise AssertionError("metacell identifiers repeat across datasets")
    var = pd.DataFrame(index=pd.Index(reference_features, name="gene"))
    if var_extra is not None:
        for column in var_extra.columns:
            var[column] = var_extra.reindex(var.index)[column].to_numpy()

    payload: Dict[str, object] = {
        "description": (
            f"Metacells of up to {size} cells. Each metacell is the sum of raw counts over a "
            f"disjoint block of randomly ordered cells from one group and cell population. A "
            f"final block of fewer than {min_tail} cells is dropped."),
        "metacell_size": str(size),
        "metacell_min_tail": str(min_tail),
        "metacell_seed": str(seed),
        "counts_type": "raw_sum",
        "metacell_membership": str(membership) if membership else "not written",
    }
    if uns_extra:
        payload.update({key: value for key, value in uns_extra.items()})

    os.makedirs(os.path.dirname(os.path.abspath(output)) or ".", exist_ok=True)

    def blocks_in_row_order() -> Iterable[sp.csr_matrix]:
        """Yield every metacell row block, dataset by dataset, holding one block at a time."""
        for plan in plans:
            produced = 0
            for block in stream_dataset_blocks(
                    plan["handle"], plan["blocks"], str(plan["matrix"]),
                    plan["column_map"], n_reference, chunk_cells=chunk_cells):
                produced += block.shape[0]
                yield block
            if produced != len(plan["blocks"]):
                raise AssertionError(
                    f"{plan['dataset']}: streamed {produced} of {len(plan['blocks'])} rows")
            LOGGER.info("[%s] summed %d metacell rows", plan["dataset"], produced)

    try:
        written = write_h5ad_csr(output, blocks_in_row_order(), obs, var, payload)
    finally:
        for handle in open_handles:
            handle.close()

    if membership:
        os.makedirs(os.path.dirname(os.path.abspath(membership)) or ".", exist_ok=True)
        opener = gzip.open if membership.endswith(".gz") else open
        with opener(membership, "wt") as handle:
            handle.write("metacell_id\tdataset\tCellBarcode\n")
            for identifier, dataset, barcode in membership_rows:
                handle.write(f"{identifier}\t{dataset}\t{barcode}\n")

    summary = dict(
        output=os.path.abspath(output),
        membership=os.path.abspath(membership) if membership else None,
        n_metacells=written["n_obs"], n_features=written["n_vars"], nnz=written["nnz"],
        n_member_cells=len(membership_rows), n_cells_available=int(len(cells)),
        seconds=round(time.time() - started, 1), datasets=reports)
    LOGGER.info("[write] %s | %d metacells x %d features | %d member cells | %.1f min",
                summary["output"], summary["n_metacells"], summary["n_features"],
                summary["n_member_cells"], summary["seconds"] / 60)
    return summary


###############################################################################
# Command line
###############################################################################


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build metacells by summing raw counts over random disjoint cell blocks.")
    parser.add_argument("--cells", required=True,
                        help="Tab-delimited cell table with dataset, h5ad, CellBarcode, group, "
                             "population and any obs columns.")
    parser.add_argument("--reference-features", required=True,
                        help="File whose first column lists the reference features, in order.")
    parser.add_argument("--output", required=True, help="Destination .h5ad.")
    parser.add_argument("--membership", default=None,
                        help="Destination for the metacell-to-barcode index (.tsv or .tsv.gz).")
    parser.add_argument("--parameters", default=None, help="Destination for the run summary JSON.")
    parser.add_argument("--obs-column", action="append", default=[],
                        help="Cell-table column to carry into obs; repeatable.")
    parser.add_argument("--population-column", default="population",
                        help="Name the population column takes in the output obs.")
    parser.add_argument("--drop-group-columns", action="store_true",
                        help="Keep dataset and group out of obs; they stay in the membership "
                             "index and the run summary.")
    parser.add_argument("--size", type=int, default=DEFAULT_SIZE)
    parser.add_argument("--min-tail", type=int, default=DEFAULT_MIN_TAIL)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--layer", default=None,
                        help="Layer to read first; the raw-count gate still applies.")
    parser.add_argument("--chunk-cells", type=int, default=CHUNK_CELLS)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s", stream=sys.stdout)
    cells = pd.read_csv(args.cells, sep="\t", dtype=str)
    features = [line.split("\t")[0].strip()
                for line in open(args.reference_features) if line.strip()]
    summary = build_metacells(
        cells, features, args.output, membership=args.membership, size=args.size,
        min_tail=args.min_tail, seed=args.seed, layer=args.layer,
        obs_columns=tuple(args.obs_column), population_column=args.population_column,
        keep_group_columns=not args.drop_group_columns, chunk_cells=args.chunk_cells)
    if args.parameters:
        with open(args.parameters, "w") as handle:
            json.dump(summary, handle, indent=2)
        LOGGER.info("[write] %s", os.path.abspath(args.parameters))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
