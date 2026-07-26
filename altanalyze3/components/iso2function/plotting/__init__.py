"""iso2function plotting: ISV-web-faithful isoform structure tracks + stacked isoform-usage bar/line plots.

Public callables:
    StructRenderer            -- ISV-web-faithful per-isoform structure-track renderer (mergeSegments fix)
    plot_stacked_isoform_bar  -- stacked isoform-composition bar/line + structure tracks (one call = one figure)
    load_pseudobulk_slice     -- load an AnnData-like pseudobulk slice from a .npz cache or .h5ad
"""
from .isoform_struct_view import StructRenderer
from .stacked_isoform_bar import plot_stacked_isoform_bar, load_pseudobulk_slice

__all__ = ["StructRenderer", "plot_stacked_isoform_bar", "load_pseudobulk_slice"]
