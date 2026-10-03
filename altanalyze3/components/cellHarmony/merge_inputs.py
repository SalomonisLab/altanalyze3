"""Stream matching 10x libraries into a single preallocated CSR matrix."""
import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

from altanalyze3.components.cellHarmony.input_validation import validate_10x_h5


def concat_matching_10x(entries, load_sample):
    """Return a merged AnnData, or None when bounded disk concat is required.

    Only modern, RNA-only 10x files with identical ordered features qualify.
    Metadata and naming come from the normal cellHarmony loader. Each expression
    matrix is loaded once, copied into the final buffers, then released.
    """
    headers = []
    feature_names = None
    for entry in entries:
        path = entry[0] if isinstance(entry, (tuple, list)) else entry
        if not str(path).endswith(".h5"):
            return None
        with h5py.File(path, "r") as handle:
            validate_10x_h5(handle, str(path))
            if "matrix/features/name" not in handle or "matrix/features/feature_type" not in handle:
                return None
            features = handle["matrix/features"]
            if not np.all(features["feature_type"][:] == b"Gene Expression"):
                return None
            names = features["name"][:]
            if feature_names is not None and not np.array_equal(feature_names, names):
                return None
            feature_names = names
            genes, cells = map(int, handle["matrix/shape"][:])
            if genes != len(names):
                return None
            headers.append((cells, genes, int(handle["matrix/data"].shape[0])))
    if not headers:
        return None
    cells = sum(row[0] for row in headers)
    nnz = sum(row[2] for row in headers)
    genes = headers[0][1]
    index_dtype = np.int64 if max(nnz, cells, genes) >= np.iinfo(np.int32).max else np.int32
    data = np.empty(nnz, dtype=np.float32)
    indices = np.empty(nnz, dtype=index_dtype)
    indptr = np.empty(cells + 1, dtype=index_dtype)
    indptr[0] = 0
    metadata = []
    row_offset = entry_offset = 0
    expected_names = None
    for number, (entry, (sample_cells, _, sample_nnz)) in enumerate(zip(entries, headers), start=1):
        if isinstance(entry, (tuple, list)) and len(entry) >= 2:
            path, name_override = entry[0], entry[1]
        else:
            path, name_override = entry, None
        print(f"[STATUS] Streaming 10x sample {number}/{len(headers)} ({sample_cells} cells, {sample_nnz} counts).", flush=True)
        sample, _ = load_sample(path, sample_name_override=name_override)
        matrix = sample.X
        if matrix.shape != (sample_cells, genes) or matrix.nnz != sample_nnz:
            raise ValueError("10x matrix dimensions changed between inspection and import.")
        if expected_names is not None and not sample.var_names.equals(expected_names):
            raise ValueError("10x feature ordering changed during sample preparation.")
        expected_names = sample.var_names
        end = entry_offset + sample_nnz
        data[entry_offset:end] = matrix.data
        indices[entry_offset:end] = matrix.indices
        indptr[row_offset + 1:row_offset + sample_cells + 1] = matrix.indptr[1:].astype(index_dtype) + entry_offset
        metadata.append(ad.AnnData(X=sp.csr_matrix((sample_cells, genes), dtype=np.float32),
                                   obs=sample.obs.copy(), var=sample.var.copy()))
        row_offset += sample_cells
        entry_offset = end
        del matrix, sample
    combined = ad.concat(metadata, label="sample", join="outer", fill_value=0)
    combined.X = sp.csr_matrix((data, indices, indptr), shape=(cells, genes), copy=False)
    return combined
