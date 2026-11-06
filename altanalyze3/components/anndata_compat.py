"""AnnData's supported HDF5 element APIs across the pinned 0.10 and newer releases.

These functions moved from experimental to io in 0.11; data and algorithms are
unchanged. Keep a single compatibility boundary rather than upgrading scientific
packages implicitly or making individual endpoints version dependent.
"""
try:
    from anndata.io import read_elem, write_elem, sparse_dataset
except ImportError:  # AnnData 0.10.x, including the repository's pinned 0.10.9
    from anndata.experimental import read_elem, write_elem, sparse_dataset

__all__ = ['read_elem', 'write_elem', 'sparse_dataset']
