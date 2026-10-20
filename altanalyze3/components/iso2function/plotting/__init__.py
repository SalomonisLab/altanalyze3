"""iso2function plotting: ISV-web-faithful isoform structure tracks + stacked isoform-usage bar/line plots.

Public callables:
    StructRenderer            -- ISV-web-faithful per-isoform structure-track renderer (mergeSegments fix)
    plot_stacked_isoform_bar  -- stacked isoform-composition bar/line + structure tracks (one call = one figure)
    load_pseudobulk_slice     -- load an AnnData-like pseudobulk slice from a .npz cache or .h5ad
"""
__all__ = ["StructRenderer", "plot_stacked_isoform_bar", "load_pseudobulk_slice"]


def __getattr__(name):
    """Load each plotting backend only when its public callable is requested.

    Protein-feature figures do not need AnnData or the genomic ISV renderer.
    Existing public imports retain the same implementations and algorithms.
    """
    if name == "StructRenderer":
        from .isoform_struct_view import StructRenderer
        return StructRenderer
    if name in ("plot_stacked_isoform_bar", "load_pseudobulk_slice"):
        from . import stacked_isoform_bar
        return getattr(stacked_isoform_bar, name)
    raise AttributeError(name)
