"""Load and align the human lung CITE-seq training data.

The RNA and the ADT arrive as two separate files with two different cell-name
conventions, so every consumer in this package goes through here:

* RNA  ``/Users/saljh8/Dropbox/Transfer/COVID_with_umap_and_markers.h5ad``
  302,922 cells x 35,545 genes. ``X`` is CP10k+log1p (natural log);
  ``layers['counts']`` holds raw counts.
* ADT  ``/Users/saljh8/Dropbox/Transfer/COVID-TotalVI/denoised_ADT_HTC_covid_all_samples.txt``
  400,870 cells x 56 TotalVI-denoised antibodies on a linear scale.

Cell-name bridge, validated on the released files::

    Sample_<Library>_<BC>-1-<k>   ->   <BC>-1.<Library>
    Sample_D105_CBL_MIX_2_AAACCCAAGCTGAGCA-1-0 -> AAACCCAAGCTGAGCA-1.D105_CBL_MIX_2

302,876 of 302,922 RNA cells (99.98%) carry ADT; 302,876 of 400,870 ADT rows
(75.55%) survive, the rest being cells the RNA h5ad dropped at QC.

Target scale: the bone marrow model trains on log1p ADT values and the
cellHarmony-web registry declares ``expression_scale: log1p, log_base: e``
for it, so ``load_adt`` applies ``np.log1p`` to the denoised values by default.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .adt_rna_map import atlas_var_name


RNA_H5AD_DEFAULT = Path("/Users/saljh8/Dropbox/Transfer/COVID_with_umap_and_markers.h5ad")
ADT_TXT_DEFAULT = Path(
    "/Users/saljh8/Dropbox/Transfer/COVID-TotalVI/denoised_ADT_HTC_covid_all_samples.txt"
)

_ADT_ROW_RE = re.compile(r"^Sample_(?P<lib>.+)_(?P<bc>[ACGT]+-\d+)-(?P<k>\d+)$")


def adt_row_to_atlas_obs(name: str) -> Optional[str]:
    """``Sample_<Library>_<BC>-1-<k>`` -> ``<BC>.<Library>``; None if unparsable."""
    match = _ADT_ROW_RE.match(str(name).strip())
    if match is None:
        return None
    return f"{match.group('bc')}.{match.group('lib')}"


@dataclass(frozen=True)
class AdtTable:
    obs_names: np.ndarray          # atlas-style cell names
    var_names: np.ndarray          # Hu.* ADT names, panel order
    values: np.ndarray             # (n_cells, n_adts) float32
    n_rows_read: int
    n_unparsed: int


def load_adt(path: Path | str = ADT_TXT_DEFAULT, *, log1p: bool = True) -> AdtTable:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    n_rows_read = int(frame.shape[0])
    converted = [adt_row_to_atlas_obs(name) for name in frame.index]
    keep = np.array([name is not None for name in converted])
    n_unparsed = int((~keep).sum())
    frame = frame.loc[keep]
    obs_names = np.array([name for name in converted if name is not None], dtype=object)
    values = frame.to_numpy(dtype=np.float32)
    if np.isnan(values).any():
        raise ValueError(f"ADT table {path} holds NaN values; refusing to impute them silently")
    if values.min() < 0:
        raise ValueError(f"ADT table {path} holds negative values; log1p would be undefined")
    if log1p:
        values = np.log1p(values.astype(np.float64)).astype(np.float32)
    var_names = np.array([atlas_var_name(c) for c in frame.columns], dtype=object)
    if len(set(var_names)) != len(var_names):
        raise ValueError("ADT name collision after cleaning; inspect adt_rna_map.clean_adt_name")
    return AdtTable(obs_names=obs_names, var_names=var_names, values=values,
                    n_rows_read=n_rows_read, n_unparsed=n_unparsed)


def read_h5ad_index(path: Path | str, key: str) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        raw = handle[key][:]
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in raw], dtype=object)


def rna_obs_names(path: Path | str = RNA_H5AD_DEFAULT) -> np.ndarray:
    return read_h5ad_index(path, "obs/_index")


def rna_var_names(path: Path | str = RNA_H5AD_DEFAULT) -> np.ndarray:
    return read_h5ad_index(path, "var/_index")


def rna_obs_column(path: Path | str, column: str) -> np.ndarray:
    """Read one obs column (categorical or plain) as strings."""
    with h5py.File(path, "r") as handle:
        node = handle[f"obs/{column}"]
        if isinstance(node, h5py.Group):
            cats = [c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][:]]
            codes = node["codes"][:]
            return np.array([cats[c] if c >= 0 else "" for c in codes], dtype=object)
        raw = node[:]
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in raw], dtype=object)


def read_rna_rows(
    path: Path | str,
    row_indices: np.ndarray,
    *,
    gene_indices: Optional[np.ndarray] = None,
    layer: Optional[str] = None,
    block_rows: int = 20000,
    dense: bool = False,
    verbose: bool = True,
):
    """Stream a CSR .h5ad and return the requested rows.

    ``row_indices`` must be sorted and unique. Returns a CSR matrix over all
    genes, or a dense ``float32`` array restricted to ``gene_indices`` when
    ``dense=True``. Streaming sequentially beats fancy-indexing a 1.2e9-nnz
    CSR on disk.
    """
    row_indices = np.asarray(row_indices, dtype=np.int64)
    if row_indices.size and (np.diff(row_indices) <= 0).any():
        raise ValueError("row_indices must be strictly increasing")
    group = "X" if layer is None else f"layers/{layer}"
    with h5py.File(path, "r") as handle:
        node = handle[group]
        encoding = node.attrs.get("encoding-type", b"")
        encoding = encoding.decode() if isinstance(encoding, bytes) else str(encoding)
        if encoding != "csr_matrix":
            raise ValueError(f"{group} is {encoding!r}; only csr_matrix is supported")
        n_genes = int(node.attrs["shape"][1]) if "shape" in node.attrs else None
        indptr = node["indptr"][:]
        n_rows_total = int(indptr.size - 1)
        if n_genes is None:
            n_genes = int(node["indices"][:].max()) + 1
        gene_lookup = None
        if gene_indices is not None:
            gene_indices = np.asarray(gene_indices, dtype=np.int64)
            gene_lookup = np.full(n_genes, -1, dtype=np.int64)
            gene_lookup[gene_indices] = np.arange(gene_indices.size, dtype=np.int64)

        wanted = np.zeros(n_rows_total, dtype=bool)
        wanted[row_indices] = True
        out_position = np.full(n_rows_total, -1, dtype=np.int64)
        out_position[row_indices] = np.arange(row_indices.size, dtype=np.int64)

        if dense:
            if gene_indices is None:
                raise ValueError("dense=True requires gene_indices")
            out = np.zeros((row_indices.size, gene_indices.size), dtype=np.float32)
        else:
            chunks: List[sp.csr_matrix] = []

        data_node = node["data"]
        indices_node = node["indices"]
        for start in range(0, n_rows_total, block_rows):
            stop = min(start + block_rows, n_rows_total)
            if not wanted[start:stop].any():
                continue
            lo, hi = int(indptr[start]), int(indptr[stop])
            block_data = data_node[lo:hi]
            block_indices = indices_node[lo:hi]
            block_indptr = indptr[start:stop + 1] - lo
            local_rows = np.where(wanted[start:stop])[0]
            if dense:
                for local in local_rows:
                    a, b = int(block_indptr[local]), int(block_indptr[local + 1])
                    cols = gene_lookup[block_indices[a:b]]
                    hit = cols >= 0
                    out[out_position[start + local], cols[hit]] = block_data[a:b][hit]
            else:
                sub = sp.csr_matrix(
                    (block_data, block_indices, block_indptr),
                    shape=(stop - start, n_genes),
                )[local_rows]
                chunks.append(sub.astype(np.float32))
            if verbose:
                print(f"[rna] rows {stop}/{n_rows_total}", flush=True)
    if dense:
        return out
    matrix = sp.vstack(chunks, format="csr") if chunks else sp.csr_matrix((0, n_genes), dtype=np.float32)
    if matrix.shape[0] != row_indices.size:
        raise ValueError(f"row count mismatch: got {matrix.shape[0]}, expected {row_indices.size}")
    return matrix


@dataclass(frozen=True)
class AlignedIndex:
    rna_rows: np.ndarray           # row indices into the RNA h5ad, increasing
    adt_rows: np.ndarray           # row indices into the ADT table, same order
    obs_names: np.ndarray          # shared cell names, same order
    n_rna_cells: int
    n_adt_rows: int


def align_cells(rna_h5ad: Path | str, adt: AdtTable) -> AlignedIndex:
    rna_names = rna_obs_names(rna_h5ad)
    adt_position: Dict[str, int] = {}
    for position, name in enumerate(adt.obs_names):
        adt_position.setdefault(str(name), position)
    rna_rows: List[int] = []
    adt_rows: List[int] = []
    shared: List[str] = []
    for row, name in enumerate(rna_names):
        position = adt_position.get(str(name))
        if position is not None:
            rna_rows.append(row)
            adt_rows.append(position)
            shared.append(str(name))
    return AlignedIndex(
        rna_rows=np.asarray(rna_rows, dtype=np.int64),
        adt_rows=np.asarray(adt_rows, dtype=np.int64),
        obs_names=np.asarray(shared, dtype=object),
        n_rna_cells=int(rna_names.size),
        n_adt_rows=int(adt.obs_names.size),
    )


def report_alignment(index: AlignedIndex) -> Dict[str, object]:
    matched = int(index.rna_rows.size)
    return {
        "n_rna_cells": index.n_rna_cells,
        "n_adt_rows": index.n_adt_rows,
        "n_matched_cells": matched,
        "rna_retained_fraction": matched / index.n_rna_cells if index.n_rna_cells else 0.0,
        "adt_retained_fraction": matched / index.n_adt_rows if index.n_adt_rows else 0.0,
    }


def load_protein_coding_symbols(path: Path | str) -> set:
    """HGNC symbols of protein-coding genes, column 0 of the EnsMart100 -PC file."""
    frame = pd.read_csv(path, sep="\t", usecols=[0])
    symbols = set(frame.iloc[:, 0].astype(str).str.strip())
    symbols.discard("")
    return symbols
