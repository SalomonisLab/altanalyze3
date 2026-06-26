"""ISV-web network tab: load the precomputed isoform-resolved interaction graph + per-(state x covariate)
CPM, and answer fast in-memory contrast queries. Pure data layer (no drawing).

The artifact (`isvweb_network.json`, produced by iso2function.network.precompute_isvweb) holds:
  edges[{src_iso,src_enst,src_gene,tgt_gene,tgt_ensg,type,score}], iso_cpm{enst:{ctx:cpm}},
  gene_cpm{ensg:{ctx:cpm}}, iso_meta{src_iso:{gene,enst,structure}}, contexts[], groups[].
A "context" is "<cell_state>||<covariate_group>". An edge has EXPRESSION EVIDENCE in a context iff BOTH
endpoints' CPM > threshold there (source isoform CPM and target-gene max-isoform CPM).
"""
import os, json


def load_network(path):
    with open(path) as fh:
        return json.load(fh)


def find_network_artifact(run_dir, gff_output_dir=None):
    for d in (run_dir, gff_output_dir):
        if d:
            p = os.path.join(d, "isvweb_network.json")
            if os.path.exists(p):
                return p
    return None


NMD_FLAGS = ("NMD", "Potential-NMD")   # isoforms predicted to undergo nonsense-mediated decay -> no stable
                                       # protein -> excluded from the protein-coding PPI/PDI network


def _contexts(cell_types, groups):
    """The "<state>||<group>" contexts a panel aggregates over = cell states x covariate groups."""
    return [f"{s}||{g}" for s in (cell_types or []) for g in (groups or [])]


def _gene_isoforms(net, ensg, ctxs, top=8):
    """The known isoforms underlying a target/partner gene's expression, ranked by mean CPM over the panel's
    contexts: [{enst, aa, cpm}, ...] (top N, CPM > 0)."""
    gi = net.get("gene_isoforms", {}); tic = net.get("tgt_iso_cpm", {}); taa = net.get("tgt_iso_aa", {})
    n = max(1, len(ctxs)); out = []
    for enst in gi.get(ensg, []):
        d = tic.get(enst)
        cpm = (sum(d.get(c, 0.0) for c in ctxs) / n) if d else 0.0
        if cpm > 0:
            out.append({"enst": enst, "aa": taa.get(enst), "cpm": round(cpm, 2)})
    out.sort(key=lambda x: x["cpm"], reverse=True)
    return out[:top]


def _agg_edges(net, cell_types, groups, threshold, edge_type):
    """Edges with expression evidence aggregated over the CARTESIAN PRODUCT of the given cell states x
    covariate groups: BOTH endpoints' MEAN CPM across those "<state>||<group>" contexts exceeds the
    threshold (mean counts 0 where a context is absent). Source isoforms flagged NMD/Potential-NMD are
    dropped (they make no stable protein)."""
    iso_cpm, gene_cpm = net["iso_cpm"], net["gene_cpm"]
    iso_meta = net.get("iso_meta", {})
    ctxs = _contexts(cell_types, groups)
    n = max(1, len(ctxs))
    def agg(cmap, node):
        d = cmap.get(node)
        return (sum(d.get(c, 0.0) for c in ctxs) / n) if d else 0.0
    out = []
    for e in net["edges"]:
        if edge_type and e["type"] != edge_type:
            continue
        if iso_meta.get(e["src_iso"], {}).get("nmd") in NMD_FLAGS:    # NMD isoform -> not protein-coding
            continue
        si = agg(iso_cpm, e["src_enst"]); tg = agg(gene_cpm, e["tgt_ensg"])
        if si <= threshold or tg <= threshold:
            continue
        out.append({**e, "src_cpm": round(si, 2), "tgt_cpm": round(tg, 2)})
    return out


def _restrict_to_genes(edges, genes):
    """Search: keep edges incident to a searched gene (its FIRST-DEGREE edges)."""
    gs = {g.upper() for g in genes}
    return [e for e in edges if e["src_gene"].upper() in gs or e["tgt_gene"].upper() in gs]


def _nodes(edges, iso_meta=None):
    iso_meta = iso_meta or {}
    nd = {}
    for e in edges:
        if e["src_iso"] not in nd:
            m = iso_meta.get(e["src_iso"], {})
            nd[e["src_iso"]] = {"id": e["src_iso"], "label": e["src_iso"], "kind": "isoform",
                                "gene": e["src_gene"], "cpm": e["src_cpm"], "enst": e.get("src_enst") or m.get("enst"),
                                "aa": m.get("aa"), "nmd": m.get("nmd")}
        nd[e["src_iso"]]["cpm"] = max(nd[e["src_iso"]]["cpm"], e["src_cpm"])
        t = nd.setdefault(e["tgt_gene"], {"id": e["tgt_gene"], "label": e["tgt_gene"], "kind": "gene",
                                          "gene": e["tgt_gene"], "cpm": e["tgt_cpm"]})
        t["cpm"] = max(t["cpm"], e["tgt_cpm"])
    return list(nd.values())


def _igraph_fr_layout(node_ids, edges):
    """Fruchterman-Reingold positions via igraph (same layout as visualization/NetPerspective.py:
    Graph(directed=True); g.layout("fr")). edges are (src_id, tgt_id) over node_ids (incl. same-gene-isoform
    glue edges). Returns {id: [x, y]} or None if igraph is unavailable."""
    try:
        from igraph import Graph
    except Exception:
        return None
    if not node_ids:
        return None
    idx = {n: i for i, n in enumerate(node_ids)}
    g = Graph(directed=True)
    g.add_vertices(len(node_ids))
    es = [(idx[s], idx[t]) for (s, t) in edges if s in idx and t in idx]
    if es:
        g.add_edges(es)
    try:
        lay = g.layout("fr")
    except Exception:
        return None
    return {node_ids[i]: [round(c[0], 4), round(c[1], 4)] for i, c in enumerate(lay.coords)}


def _combined_layout(a, b):
    """igraph FR over the UNION of both panels' nodes + edges + same-gene-isoform glue edges (so the two
    sides share identical positions)."""
    node_ids, kind = [], {}
    for side in (a, b):
        for n in side["nodes"]:
            if n["id"] not in kind:
                kind[n["id"]] = (n["kind"], n["gene"]); node_ids.append(n["id"])
    edges, seen = [], set()
    for side in (a, b):
        for e in side["edges"]:
            k = (e["src_iso"], e["tgt_gene"])
            if k not in seen:
                seen.add(k); edges.append(k)
    by_gene = {}
    for nid in node_ids:
        kd, gene = kind[nid]
        if kd == "isoform":
            by_gene.setdefault(gene, []).append(nid)
    for arr in by_gene.values():
        for i in range(1, len(arr)):
            edges.append((arr[i - 1], arr[i]))      # invisible glue between isoforms of one gene
    return _igraph_fr_layout(node_ids, edges)


def query_network(net, cell_types, groups, threshold=1.0, genes=None, edge_type=None,
                  panel_by="covariate", max_edges=2000):
    """Two side-by-side panels (one per WINDOW), selected by `panel_by` (mirrors the Molecule view):
      * 'covariate' : window A/B = the first two selected COVARIATE groups; each panel SUMS (mean CPM) over
                      ALL selected cell types.
      * 'cell_type' : window A/B = the first two selected CELL TYPES; each panel aggregates over ALL selected
                      covariate groups.
    >2 selected on the panel axis -> the first two are used. An edge shows in a panel iff both endpoints'
    mean CPM over that panel's contexts exceeds the threshold. genes=None -> global; else search + 1st-degree."""
    cts = list(cell_types or []); grps = list(groups or [])
    def side(panel_cts, panel_grps, label):
        ed = _agg_edges(net, panel_cts, panel_grps, threshold, edge_type)
        if genes:
            ed = _restrict_to_genes(ed, genes)
        trunc = 0
        if len(ed) > max_edges:
            ed = sorted(ed, key=lambda e: (e.get("score") or 0), reverse=True)[:max_edges]; trunc = 1
        nodes = _nodes(ed, net.get("iso_meta"))
        ctxs = _contexts(panel_cts, panel_grps)              # gene nodes: list the isoforms underlying expression
        g2ensg = {e["tgt_gene"]: e["tgt_ensg"] for e in ed}
        for nd in nodes:
            if nd.get("kind") == "gene":
                nd["isos"] = _gene_isoforms(net, g2ensg.get(nd["id"]), ctxs)
        return {"ctx": label, "edges": ed, "nodes": nodes, "truncated": trunc}
    empty = {"ctx": "(pick contexts)", "edges": [], "nodes": [], "truncated": 0}
    if panel_by == "cell_type":
        lbl = lambda s: f"{s} · {len(grps)} covariate" + ("s" if len(grps) != 1 else "")
        a = side([cts[0]], grps, lbl(cts[0])) if len(cts) >= 1 else empty
        b = side([cts[1]], grps, lbl(cts[1])) if len(cts) >= 2 else (a if len(cts) == 1 else empty)
        axis = cts[:2]
    else:
        lbl = lambda g: f"{g} · {len(cts)} cell type" + ("s" if len(cts) != 1 else "")
        a = side(cts, [grps[0]], lbl(grps[0])) if len(grps) >= 1 else empty
        b = side(cts, [grps[1]], lbl(grps[1])) if len(grps) >= 2 else (a if len(grps) == 1 else empty)
        axis = grps[:2]
    return {"A": a, "B": b, "layout": _combined_layout(a, b), "mode": ("search" if genes else "global"),
            "threshold": threshold, "edge_type": edge_type or "all", "panel_by": panel_by,
            "cell_types": cts, "groups": grps[:2], "axis": axis}
