"""Layer 3: assemble the interaction graph, export to Cytoscape, compute isoform-switch differentials."""

from .switch_pair_plots import (
    bait_network_with_m1h,
    isoform_set_pdi_m1h,
    usage_line,
    pair_umap,
    load_pdi_m1h,
    state_usage_from_pseudobulk,
    render_switch_pairs,
)

__all__ = ["bait_network_with_m1h", "isoform_set_pdi_m1h", "usage_line", "pair_umap", "load_pdi_m1h",
           "state_usage_from_pseudobulk", "render_switch_pairs"]
