"""Cytoscape export for iso2function graphs.

Produces (a) Cytoscape.js ``elements`` JSON (for the ISV web viewer / browser embedding) and (b) a SIF
edge file (for Cytoscape Desktop), both from the (nodes_df, edges_df) pair built by ``network.build``.

Styling follows lab convention: node/edge colors are RGB hex (no named colors). A minimal style block
is emitted alongside the elements so the graph renders consistently.
"""

import os
import json
import logging

from .. import config

logger = logging.getLogger(__name__)

# RGB-hex style palette (no named colors). Continuous scores can later use the cyan->yellow ramp.
NODE_COLORS = {
    "isoform": "#1F77B4",      # blue
    "protein": "#2CA02C",      # green  (PPI partner)
    "dna_target": "#D62728",   # red    (PDI DNA bait)
}
EDGE_COLORS = {
    "ppi": "#7F7F7F",          # grey
    "pdi": "#9467BD",          # purple
}
EDGE_UNDETECTED_COLOR = "#D9D9D9"  # light grey for assayed-negative edges


def to_cytoscape_elements(nodes_df, edges_df, include_undetected=False):
    """Return a Cytoscape.js elements dict: {'nodes': [...], 'edges': [...]}. Each element is
    {'data': {...}} with all DataFrame columns carried as data attributes. Assayed-negative edges are
    excluded unless ``include_undetected`` (they remain available in the edge table either way)."""
    nodes = []
    for _, r in nodes_df.iterrows():
        data = {k: ("" if v is None else v) for k, v in r.to_dict().items()}
        data.setdefault("id", data.get("id"))
        nodes.append({"data": data})
    edges = []
    for i, r in edges_df.iterrows():
        detected = bool(r.get("detected", True))
        if not detected and not include_undetected:
            continue
        data = {k: ("" if v is None else v) for k, v in r.to_dict().items()}
        data["id"] = f"e{i}"
        edges.append({"data": data})
    return {"nodes": nodes, "edges": edges}


def default_style():
    """A minimal Cytoscape.js style array keyed on node_type / edge_type / detected."""
    style = [{
        "selector": "node",
        "style": {"label": "data(id)", "font-size": "8px", "width": 18, "height": 18,
                  "background-color": NODE_COLORS["protein"]},
    }]
    for ntype, color in NODE_COLORS.items():
        style.append({"selector": f'node[node_type = "{ntype}"]', "style": {"background-color": color}})
    style.append({"selector": "edge",
                  "style": {"width": 1.5, "curve-style": "straight", "line-color": EDGE_COLORS["ppi"]}})
    for etype, color in EDGE_COLORS.items():
        style.append({"selector": f'edge[edge_type = "{etype}"]', "style": {"line-color": color}})
    style.append({"selector": "edge[?detected = false]",
                  "style": {"line-color": EDGE_UNDETECTED_COLOR, "line-style": "dashed"}})
    return style


def write_cytoscape_json(nodes_df, edges_df, out_path=None, include_undetected=False):
    """Write a single JSON with 'elements' + 'style' for Cytoscape.js. Returns the path."""
    out_path = out_path or os.path.join(config.ARTIFACTS_DIR, "iso2function_graph.cyjs.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    payload = {"format_version": "1.0",
               "elements": to_cytoscape_elements(nodes_df, edges_df, include_undetected),
               "style": default_style()}
    with open(out_path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    logger.info("[iso2function.network] wrote Cytoscape JSON: %s", out_path)
    return out_path


def write_sif(edges_df, out_path=None, detected_only=True):
    """Write a SIF (source <edge_type> target) file for Cytoscape Desktop. Returns the path."""
    out_path = out_path or os.path.join(config.ARTIFACTS_DIR, "iso2function_graph.sif")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        for _, r in edges_df.iterrows():
            if detected_only and not bool(r.get("detected", True)):
                continue
            fh.write(f"{r['source']}\t{r.get('assay', r['edge_type'])}\t{r['target']}\n")
    return out_path
