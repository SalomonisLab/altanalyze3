"""Lightweight checks for known non-matrix 10x uploads."""
from pathlib import Path

import h5py


def validate_10x_h5(source, filename=None):
    """Reject molecule-info layouts without reading molecule or count arrays.

    Other layouts retain the existing importer's validation, including legacy 10x.
    Accept a path, an open HDF5 handle, or an upload's seekable file object.
    """
    def check(handle):
        if "matrix" not in handle and "count" in handle and "barcode_idx" in handle:
            name = filename or (Path(source).name if isinstance(source, (str, Path)) else "This upload")
            raise ValueError(
                f"{name} is a 10x molecule_info.h5 file containing per-molecule records, "
                "not a cell-by-gene count matrix. Upload filtered_feature_bc_matrix.h5 "
                "from the Cell Ranger outs folder, or an H5AD count matrix instead."
            )

    if isinstance(source, h5py.File):
        check(source)
        return
    position = source.tell() if hasattr(source, "tell") else None
    try:
        try:
            handle = h5py.File(source, "r")
        except OSError:
            return  # Not a readable HDF5 file: leave diagnosis to the normal importer.
        with handle:
            check(handle)
    finally:
        if position is not None:
            source.seek(position)
