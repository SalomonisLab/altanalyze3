"""Backward-compatible import for shared regulatory analysis in cellHarmony."""
from altanalyze3.components.cellHarmony.grn_analysis import (  # noqa: F401
    DEFAULT_EDGE_PERCENTILE, DEFAULT_EXPRESSION_PERCENTILE, DEFAULT_LIMIT,
    MAX_LIMIT, TARGET_FDR, regulator_network, tf_activity_profile, tf_activity_state_chat, read_regulatory_question,
    _split_edge, _sibling_comparison, _rna_contrast_tail, _deg_rows, _state_column,
)
