"""Layer 3 (network) — assemble the isoform interaction graph and compute isoform-switch differentials.

Consumes the tidy edge tables (``ingest``) and the crosswalk (``crosswalk``). The graph is represented
as two plain DataFrames (nodes, edges) — no third-party graph dependency — which is exactly the shape
the Cytoscape exporter and any downstream consumer need.

Node types: ``isoform`` (an assayed clone; structure/ENST attached only where Leg A is verified),
``protein`` (a PPI partner), ``dna_target`` (a PDI DNA bait). Edge types: ``ppi`` (Y2H/N2H),
``pdi`` (eY1H). Every edge carries ``detected`` (True/False) and ``tested=True`` so the tested-vs-
untested distinction survives into the graph: a ``detected=False`` edge is an assayed non-interaction,
not a missing one.
"""

import os
import logging

import pandas as pd

from .. import config

logger = logging.getLogger(__name__)


def load_edges(data_dir=None):
    """Load the tidy CALLED interaction edge tables into a dict (binary-only policy: Y2H + eY1H, the
    authors' own Boolean calls). N2H is intentionally NOT loaded here — it is score-only with no numeric
    cutoff in the metadata and is never rendered as a called interaction. Missing tables warn."""
    data_dir = data_dir or config.DATA_DIR
    out = {}
    for modality in ("ppi_y2h", "pdi_ey1h"):   # N2H excluded by design
        path = os.path.join(data_dir, f"{modality}.tsv")
        if os.path.exists(path):
            out[modality] = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
        else:
            logger.warning("[iso2function.network] edge table missing (run ingest): %s", path)
    return out


def _as_bool(v):
    return str(v).strip().lower() in config.TRUE_TOKENS


def build_graph(clone_ids=None, data_dir=None, crosswalk_structure=None, detected_only=False):
    """Build (nodes_df, edges_df) for the interaction graph.

    clone_ids        : optional iterable of paper2 clone_ids to restrict the graph to (e.g. the
                       isoforms of one gene, or the up/down isoforms of a switch). None = all.
    crosswalk_structure : optional DataFrame from crosswalk.build_crosswalk()['structure'] used to
                       annotate isoform nodes with ENST / structure / NMD where a VERIFIED clone->GFF
                       link exists (see note below). Not required for the graph topology.
    detected_only    : if True, drop ``detected=False`` (assayed-negative) edges from the graph.

    NOTE: interaction data is keyed on paper2 clone_id. Until Leg A is sequence-verified, isoform nodes
    are keyed on clone_id and carry NO asserted ENST/structure. ``crosswalk_structure`` annotation is
    applied only through a verified clone->GFF map (passed by the caller), never via the unverified
    id-index hint.
    """
    edges = load_edges(data_dir)
    clone_filter = set(map(str, clone_ids)) if clone_ids is not None else None
    node_attr = {}   # node_id -> dict
    edge_rows = []

    def _touch(node_id, ntype, **attrs):
        rec = node_attr.setdefault(node_id, {"id": node_id, "node_type": ntype})
        for k, v in attrs.items():
            if v not in (None, "") and k not in rec:
                rec[k] = v

    # ---- PPI: Y2H (ad_clone_id -> db_gene_symbol) ----
    if "ppi_y2h" in edges:
        df = edges["ppi_y2h"]
        for _, r in df.iterrows():
            src = r.get("ad_clone_id")
            if clone_filter is not None and str(src) not in clone_filter:
                continue
            detected = _as_bool(r.get("detected", r.get("y2h_result_raw", "")))
            if detected_only and not detected:
                continue
            partner = r.get("db_gene_symbol") or r.get("db_orf_id")
            _touch(src, "isoform", gene_symbol=r.get("ad_gene_symbol"))
            _touch(partner, "protein", gene_symbol=r.get("db_gene_symbol"),
                   category=r.get("db_gene_category"))
            edge_rows.append({"source": src, "target": partner, "edge_type": "ppi", "assay": "y2h",
                              "detected": detected, "tested": True, "score": ""})

    # N2H is deliberately NOT added as edges (score-only, no binary cutoff in metadata). Its log2 NLR
    # scores remain available in data/ppi_n2h.tsv for reference / future annotation if a cutoff is set.

    # ---- PDI: eY1H (clone_id -> DNA bait) ----
    if "pdi_ey1h" in edges:
        df = edges["pdi_ey1h"]
        for _, r in df.iterrows():
            src = r.get("clone_id")
            if clone_filter is not None and str(src) not in clone_filter:
                continue
            detected = _as_bool(r.get("detected", ""))
            if detected_only and not detected:
                continue
            bait = r.get("bait_id")
            _touch(src, "isoform", gene_symbol=r.get("gene_symbol"))
            _touch(bait, "dna_target")
            edge_rows.append({"source": src, "target": bait, "edge_type": "pdi", "assay": "ey1h",
                              "detected": detected, "tested": True, "score": ""})

    nodes_df = pd.DataFrame(list(node_attr.values()))
    edges_df = pd.DataFrame(edge_rows)

    # annotate isoform nodes with the paper-defined M1H call (activator/repressor/neither), if ingested
    m1h_path = os.path.join(data_dir or config.DATA_DIR, "activation_m1h.tsv")
    if os.path.exists(m1h_path) and len(nodes_df):
        m1h = pd.read_csv(m1h_path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
        if {"clone_id", "m1h_call"} <= set(m1h.columns):
            mp = dict(zip(m1h["clone_id"].astype(str), m1h["m1h_call"]))
            nodes_df["m1h_call"] = nodes_df["id"].map(mp).fillna("")

    # annotate isoform nodes with the canonical STRUCTURE KEY (+ ENST / observed isoform / NMD) from the
    # resolved clone->structure crosswalk (clone_to_structure.tsv). This is what makes the structure key,
    # not the paper clone label, the identity carried on each isoform node.
    cts = crosswalk_structure
    if cts is None:
        cts_path = os.path.join(data_dir or config.DATA_DIR, "clone_to_structure.tsv")
        if os.path.exists(cts_path):
            cts = pd.read_csv(cts_path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    if cts is not None and len(nodes_df) and "clone_id" in cts.columns:
        ann = cts.drop_duplicates("clone_id").set_index("clone_id")
        for col in ("structure", "enst", "final_isoform_id", "nmd_status", "match_type"):
            if col in ann.columns:
                nodes_df[col] = nodes_df["id"].map(ann[col].to_dict()).fillna("")
    logger.info("[iso2function.network] graph: %d nodes, %d edges (clones=%s, detected_only=%s)",
                len(nodes_df), len(edges_df), "all" if clone_filter is None else len(clone_filter),
                detected_only)
    return nodes_df, edges_df


def _detected_partners(edges_df, clone_id, edge_type=None):
    """Set of targets an isoform interacts with (detected=True), optionally restricted to an edge type."""
    if not len(edges_df):
        return set()
    m = (edges_df["source"].astype(str) == str(clone_id)) & (edges_df["detected"] == True)  # noqa: E712
    if edge_type:
        m &= (edges_df["edge_type"] == edge_type)
    return set(edges_df.loc[m, "target"].astype(str))


def differential_interactome(clone_up, clone_down, edges_df, edge_type=None):
    """Gained/lost/shared interaction targets for an isoform switch from ``clone_down`` -> ``clone_up``.

    Uses only assayed-positive edges of BOTH clones, so a target is called 'lost' only if it was a
    detected partner of the down-isoform AND assayed (present as an edge) for the up-isoform — i.e. we
    do not call loss where the up-isoform was simply not tested against that target. Returns a dict with
    gained/lost/shared target sets and the tested-overlap used as the denominator.
    """
    up_pos = _detected_partners(edges_df, clone_up, edge_type)
    down_pos = _detected_partners(edges_df, clone_down, edge_type)
    up_tested = set(edges_df.loc[edges_df["source"].astype(str) == str(clone_up), "target"].astype(str))
    down_tested = set(edges_df.loc[edges_df["source"].astype(str) == str(clone_down), "target"].astype(str))
    co_tested = up_tested & down_tested
    return {
        "clone_up": clone_up, "clone_down": clone_down, "edge_type": edge_type or "all",
        "gained": sorted((up_pos - down_pos) & co_tested),   # positive in up, negative-but-tested in down
        "lost": sorted((down_pos - up_pos) & co_tested),     # positive in down, negative-but-tested in up
        "shared": sorted(up_pos & down_pos),
        "co_tested_targets": sorted(co_tested),
    }


def switch_consequence(ref_clone, alt_clone, data_dir=None):
    """Return the authors' high-confidence functional-consequence categories (Data_S8) for a
    reference->alternative isoform switch, as a dict (PDI_category, PPI_category,
    M1H_activation_category, localization_category, DBD_pct_lost_in_alt, classifications), or None if
    the pair is not in Data_S8. This is the paper's own call — preferred over re-deriving from edges."""
    path = os.path.join(data_dir or config.DATA_DIR, "functional_categories.tsv")
    if not os.path.exists(path):
        logger.warning("[iso2function.network] functional_categories.tsv missing (run ingest)")
        return None
    fc = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    m = fc[(fc["reference_isoform"].astype(str) == str(ref_clone)) &
           (fc["alternative_isoform"].astype(str) == str(alt_clone))]
    if not len(m):
        return None
    return m.iloc[0].to_dict()


def write_graph(nodes_df, edges_df, out_dir=None, prefix="iso2function_graph"):
    """Write node/edge attribute tables (TSV) for inspection / Cytoscape import. Returns the paths."""
    out_dir = out_dir or config.ARTIFACTS_DIR
    os.makedirs(out_dir, exist_ok=True)
    npath = os.path.join(out_dir, f"{prefix}.nodes.tsv")
    epath = os.path.join(out_dir, f"{prefix}.edges.tsv")
    nodes_df.to_csv(npath, sep="\t", index=False)
    edges_df.to_csv(epath, sep="\t", index=False)
    return {"nodes": npath, "edges": epath}
