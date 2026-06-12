"""altanalyze3 UDON — a superset of the original pyudon package.

Exposes the original pyudon API (so pyudon-style code such as `run_example_dataset.py` works):
  pre_process            calculate_pseudobulk_folds, calculate_pseudobulks, pseudobulk_wrapper,
                         matched_controls, collective_controls, remove_negatives, pre_pseudobulk_qc
  feature_selection      feature_selection_wrapper, ...
  clustering_wrapper     clustering_wrapper, find_udon_clusters
  enrichment_pathwaycommons  enrichment_wrapper (PathwayCommons)
  visualizations         plot_markers_df, plot_metadata_heatmap, assign_rainbow_colors_to_groups
  satay (satay_udon_core) make_mean_based_binary, fishers_clinical_feats, cmh_clinical_feats,
                         fdr_correction, satay_udon, load_metadata, render_satay
  sashimi_udon           determine_control_cells, determine_cell_type_hvg/deg, check_similarity_distribution

plus the altanalyze3 additions:
  pseudobulk_protocol    matched sex+platform controls, build_fold_matrix, filter_udon_genes (gene filter),
                         restrict_folds, predict_sample_sex, select_matched_controls
  goelite_enrichment     run_goelite_on_udon (GO-Elite, the preferred enrichment over PathwayCommons)
  markerFinder           marker_finder_wrapper (== pyudon metadata_analysis)

Flat module imports inside the package are resolved via a sys.path shim (the modules are run both as
standalone CLIs and imported here). Demo-data loaders live in the optional `datasets` subpackage
(`from altanalyze3.components.udon.datasets import fetch_geo_raw`) and are not imported here (avoid the
hard `pooch` dependency).
"""
import os as _os
import sys as _sys

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

# pyudon-compatible API
from pre_process import *            # noqa: E402,F401,F403
from feature_selection import *      # noqa: E402,F401,F403
from clustering_wrapper import *     # noqa: E402,F401,F403
from enrichment_pathwaycommons import *  # noqa: E402,F401,F403
from visualizations import *         # noqa: E402,F401,F403
from markerFinder import *           # noqa: E402,F401,F403  (pyudon metadata_analysis)
from satay_udon_core import *        # noqa: E402,F401,F403  (canonical SATAY API)
from sashimi_udon import *           # noqa: E402,F401,F403  (control-cell determination)

# altanalyze3 additions
from pseudobulk_protocol import *    # noqa: E402,F401,F403
from goelite_enrichment import *     # noqa: E402,F401,F403
