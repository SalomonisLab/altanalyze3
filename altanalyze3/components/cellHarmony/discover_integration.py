"""Integrated views ported from LungMAP Discover (2026-09-16).

Selection, class mapping, pathway ranking and geometry follow molecule/discover.py.
Storage is supplied by a dataset adapter; serving a view never computes a differential.
"""
from __future__ import annotations
import json
import re
from pathlib import Path
from functools import lru_cache
from typing import Any
import numpy as np

RESOURCE_ROOT = Path(__file__).parent / "resources" / "pathways"

@lru_cache(maxsize=1)
def _diagrams():
    return json.loads((RESOURCE_ROOT / "wikipathways_diagrams.json").read_text())

@lru_cache(maxsize=1)
def _diagram_synonyms():
    return json.loads((RESOURCE_ROOT / "wikipathways_lipids.json").read_text()).get("class_synonyms", {})

CLASS_GROUPS = {
    "glycerophospholipid": ["PC", "PE", "PI", "PS", "PG"],
    "phospholipid": ["PC", "PE", "PI", "PS", "PG"],
    "glycosphingolipid": ["HEXCER", "LACCER", "GM3"],
    "sphingolipid": ["CER", "SM", "HEXCER", "LACCER", "GM3"],
    "neutral lipid": ["DG", "TG", "CE"],
}

def required_note(modality):
    label = "grn" if modality == "grn" else "metabolite" if modality == "metabolite" else "lipid"
    view = "network" if label == "grn" else "pathway"
    return (f"Please note: both differential expression and {label} imputed differentials "
            f"are required for {view} integrated visualization.")

def _label(significance):
    return "reported calls" if significance == "reported" else "raw p" if significance == "pval" else "FDR"

def _feature_differentials(ds, contrast, cell_state, names, modality="rna"):
    want=set(names)
    return {r["gene"]:r for r in ds.differentials(contrast, cell_state, modality) if r["gene"] in want}

def _passes(row, significance, cutoff):
    return significance == "reported" or (row.get(significance) is not None and row[significance] <= cutoff)


def _gate(significance, cutoff):
    return "the stored differential selection" if significance == "reported" else f"{_label(significance)} {cutoff}"


def _edge_stats(ds, contrast, cell_state, max_fdr, min_fold, significance):
    cut=float(np.log2(min_fold)) if min_fold>1 else 0.
    return {r["gene"]:(r["log2fc"],r.get("fdr" if significance == "reported" else significance))
            for r in ds.differentials(contrast,cell_state,"grn")
            if _passes(r,significance,max_fdr)
            and r.get("log2fc") is not None and abs(r["log2fc"])>=cut}

def _pathway_rows(ds, contrast, cell_state, modality, max_fdr, significance):
    rows=[]
    for kind in ("rna",modality):
        for r in ds.differentials(contrast,cell_state,kind):
            if r.get("log2fc") is not None:
                rows.append((r["gene"], "rna" if kind=="rna" else "metabolite" if modality=="metabolite" else "lipid",r["log2fc"]))
    return rows



def _best_group(label: str) -> tuple[str, list[str]]:
    """The parent category a metabolite label names, and the classes under it.

    The longest matching term wins, for the same reason `_best_class` needs it:
    "glycosphingolipid" contains "sphingolipid".
    """
    low = (label or "").lower()
    if not low:
        return "", []
    best, best_len = "", 0
    for term in CLASS_GROUPS:
        if term in low and len(term) > best_len:
            best, best_len = term, len(term)
    return (best, CLASS_GROUPS[best]) if best else ("", [])

def _best_class(label: str, synonyms: dict[str, list[str]]) -> str:
    """The lipid class a metabolite label names, resolving nested names by specificity.

    Class names nest as substrings: "galactosylceramide" contains "ceramide" and also
    "lactosylceramide". A first-match loop therefore depended on dictionary order and
    labelled a galactosylceramide node Cer. The longest matching synonym is the most
    specific one, so it wins; ties keep the class name for stability.
    """
    low = (label or "").lower()
    if not low:
        return ""
    best, best_len = "", 0
    for name, words in synonyms.items():
        for word in words:
            w = word.lower()
            if w and w in low and (len(w) > best_len
                                   or (len(w) == best_len and name < best)):
                best, best_len = name, len(w)
    return best

def _discover_network(ds, cell_state: str, features: list[str], limit: int = 50,
                      contrast: str = "", edge_percentile: float = 50.0,
                      expression_percentile: float = 50.0,
                      max_fdr: float = 0.05,
                      min_fold: float = 1.2,
                      min_expression: float = 0.5, significance: str = "fdr") -> dict[str, Any]:
    """Which transcription factors regulate the features the reader is looking at.

    WHAT THIS BRIDGES. Nathan, 2026-09-06: "the network should restrict TF target genes
    that are transcriptionally regulated (up or down) in the selected comparison and cell
    type ... The TFs do not need to be regulated themselves since their activity but not
    expression might be changed."

    So the TARGETS come from the finder and the EDGES come from the regulatory model's
    per-cell-state edge list.

    THREE THINGS THE FIRST VERSION GOT WRONG, all reported by Nathan on 2026-09-06:

      * every edge was drawn, so 50 factors reached 12 features through 130 edges of
        indistinguishable weight. The edge score now has a floor, set as a PERCENTILE of
        the edges in this cell state rather than as a number somebody invented: in AT2 the
        median edge scores 0.337 and the strongest tenth scores above 1.191.
      * factors that the cell type does not express were drawn anyway. BATF3 scores 0.004
        and IRF5 0.006 on log2(1 + CP10k) in AT2, which is not expression. A factor now
        has to clear a percentile of the expression of every modelled factor in this cell
        state, measured over THE DONORS OF THE COMPARISON'S TWO ARMS where the comparison
        is one the atlas can address.
      * every factor was painted the same colour, which read as "no factor is
        differential". Each factor now carries its own differential expression and its own
        differential activity, so the drawing can colour by one and hover the other.
    """

    if not cell_state:
        return {"nodes": [], "edges": [],
                "note": "Choose a cell type: a regulatory edge is scored within one."}
    wanted = [f for f in features if f]
    if not wanted:
        return {"nodes": [], "edges": [], "note": "No features to place in a network."}
    try:
        per_state = ds.edge_scores(cell_state)
    except FileNotFoundError as exc:
        return {"nodes": [], "edges": [],
                "note": f"the regulatory edge index is not deployed: {exc}"}
    if not per_state:
        return {"nodes": [], "edges": [],
                "note": f"The regulatory model records no edge in {cell_state}."}

    # ---- THE EDGE GATE. The only rule now: the edge's own differential must be significant
    # in this comparison and cell state, at the same threshold the gene expression uses.
    # The percentile floors below are kept for reporting but no longer drop anything.
    sig_edges = set(_edge_stats(ds, contrast, cell_state, max_fdr, min_fold, significance))
    all_scores = np.array([float(v[0]) if isinstance(v, (list, tuple)) else float(v)
                           for v in per_state.values()], dtype=np.float64)
    pct = min(99.0, max(0.0, float(edge_percentile)))
    edge_cut = float(np.percentile(all_scores, pct)) if all_scores.size else 0.0

    target_set = set(wanted)
    hits: dict[str, dict[str, float]] = {}
    n_edges_seen = 0
    for key, value in per_state.items():
        parts = [p for p in key.split("|") if p]
        if len(parts) < 2:
            continue
        regulators, target = parts[:-1], parts[-1]
        if target not in target_set:
            continue
        score = float(value[0]) if isinstance(value, (list, tuple)) else float(value)
        n_edges_seen += 1
        # The edge gate, and the only gate. An edge whose own differential did not clear the
        # threshold is not drawn, however strong the model thinks the connection is.
        if sig_edges is not None and key not in sig_edges:
            continue
        for regulator in regulators:
            slot = hits.setdefault(regulator, {})
            if target not in slot or score > slot[target]:
                slot[target] = score

    if not hits:
        return {"nodes": [], "edges": [], "cell_state": cell_state,
                "n_regulators_total": 0, "edge_cut": round(edge_cut, 4),
                "note": (f"No edge into these {len(wanted)} features passes {_gate(significance,max_fdr)} "
                         f"and fold {min_fold} in this comparison.")}

    # ---- the expression floor, over the donors of the comparison's arms
    every_factor = sorted({t for key in per_state for t in key.split("|")[:-1] if t})
    per_case, per_ctrl = ds.arm_expression(contrast, cell_state, every_factor)
    levels = {f: max(per_case.get(f, 0.0), per_ctrl.get(f, 0.0))
              for f in set(per_case) | set(per_ctrl)}
    arm_scope = "real sample pseudobulks in the two comparison arms"
    expression_cut = float(min_expression)
    silent = [r for r in hits if r in levels and levels[r] <= expression_cut]
    unmeasured = [r for r in hits if r not in levels]
    # APPLIED, not merely reported. A factor the cell type does not express is not drawn,
    # however significant its edge: the edge model scores a connection it learned elsewhere.
    for name in silent:
        hits.pop(name, None)

    if not hits:
        return {"nodes": [], "edges": [], "cell_state": cell_state,
                "edge_cut": round(edge_cut, 4),
                "expression_cut": round(expression_cut, 4),
                "note": (f"Every factor reaching these features sits at or below "
                         f"{expression_cut:.2f} on pseudobulk log2(1 + CP10k) in "
                         f"{cell_state}, in both arms of this comparison, so none is "
                         f"expressed here. Lower the expression floor to see them.")}

    ranked = sorted(hits.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    kept = ranked[: max(1, int(limit))]
    drawn_targets = sorted({t for _r, slot in kept for t in slot})
    order = {name: i for i, name in enumerate(wanted)}
    drawn_targets.sort(key=lambda t: order.get(t, 9999))

    # ---- what each factor did, so the drawing can colour it
    factor_names = [r for r, _slot in kept]
    de = _feature_differentials(ds, contrast, cell_state, factor_names, "rna")
    activity = _feature_differentials(ds, contrast, cell_state, factor_names, "grn_tf")
    target_de = _feature_differentials(ds, contrast, cell_state, drawn_targets, "rna")

    nodes = []
    for t in drawn_targets:
        row = target_de.get(t, {})
        nodes.append({"id": t, "kind": "target", "rank": order.get(t),
                      "log2fc": row.get("log2fc"), "fdr": row.get("fdr"),
                      "expression": round(levels[t], 3) if t in levels else None})
    for r, slot in kept:
        row = de.get(r, {})
        act = activity.get(r, {})
        nodes.append({"id": r, "kind": "regulator", "n_targets": len(slot),
                      "log2fc": row.get("log2fc"), "fdr": row.get("fdr"),
                      "activity_log2fc": act.get("log2fc"),
                      "activity_fdr": act.get("fdr"),
                      "expression": round(levels[r], 3) if r in levels else None})
    # THE EDGE CARRIES ITS OWN STATISTIC. Nathan, 2026-09-11: "Hovering over the edge should
    # show the edge differentials." The gate already decided on these numbers, so shipping
    # them with the edge costs one lookup and lets the drawing say WHY an edge survived
    # rather than only that it did. A composite element names several factors before the
    # target, so fall back to any key whose last field is the target and which lists this
    # factor.
    estats = _edge_stats(ds, contrast, cell_state, max_fdr, min_fold, significance)

    def _edge_stat(tf: str, target: str) -> tuple:
        hit = estats.get(f"{tf}|{target}")
        if hit:
            return hit
        for key, val in estats.items():
            parts = key.split("|")
            if parts[-1] == target and tf in parts[:-1]:
                return val
        return (None, None)

    edges = []
    for r, slot in kept:
        for t, score in sorted(slot.items()):
            efc, ep = _edge_stat(r, t)
            edges.append({"source": r, "target": t, "score": round(score, 4),
                          "log2fc": None if efc is None else round(float(efc), 4),
                          "fdr": None if ep is None else float(ep),
                          "significance": "FDR" if significance == "reported" else _label(significance)})
    nodes=list({node["id"]:node for node in nodes}.values())
    covered = len({t for _r, slot in ranked for t in slot})
    n_de = sum(1 for r in factor_names
               if de.get(r) and _passes(de[r],significance,max_fdr))
    n_act = sum(1 for r in factor_names
                if activity.get(r) and _passes(activity[r],significance,max_fdr))
    return {
        "cell_state": cell_state, "contrast": contrast,
        "nodes": nodes, "edges": edges,
        "n_features_asked": len(wanted), "n_features_modelled": covered,
        "n_regulators_total": len(ranked), "n_regulators_drawn": len(kept),
        "n_regulators_silent": len(silent), "n_regulators_unmeasured": len(unmeasured),
        "n_edges_before_cut": n_edges_seen, "n_edges_drawn": len(edges),
        "edge_cut": round(edge_cut, 4), "edge_percentile": pct,
        "expression_cut": round(expression_cut, 4),
        "min_expression": expression_cut,
        "n_regulators_differential": n_de, "n_regulators_activity_differential": n_act,
        "arm_scope": arm_scope,
        # SAY WHICH GATE REMOVED WHAT. Two rules now act, and a reader who cannot tell
        # them apart cannot tell "this connection did not change" from "this factor is not
        # expressed here" -- different findings with different consequences.
        # The count cap was removed from the panel: the fold and expression thresholds
        # decide what is shown. The caption states the totals so a reader can see the
        # size of the network without a control that silently truncated it.
        "note": (f"Showing {len(drawn_targets)} genes and {len(kept)} transcription "
                 f"factors. {covered} of the {len(wanted)} features shown are targets in the "
                 f"{cell_state} regulatory model. An edge is drawn only when its own "
                 f"differential clears {_gate(significance,max_fdr)} and fold {min_fold} "
                 f"in this comparison, which leaves {len(edges)} of {n_edges_seen}. "
                 f"{len(silent)} factors were removed for sitting at or below "
                 f"{expression_cut:.2f} on pseudobulk log2(1 + CP10k) in both arms, "
                 f"measured over {arm_scope}. Of the {len(kept)} drawn, {n_de} change "
                 f"expression and {n_act} change activity under {_gate(significance,max_fdr)}. "
                 f"Each factor is connected by a retained differential edge."),
    }

def _discover_pathway_ranking(ds, cell_state: str = "", contrast: str = "",
                    max_fdr: float = 0.05, modality: str = "lipid", significance: str = "fdr") -> list[dict[str, Any]]:
    """Every lipid pathway, ranked by how much of it this comparison would actually paint.

    Nathan, 2026-09-07: the menu should open on the pathway that lights up, not on whichever
    names the most head groups. Scoring by class count alone put lipid AT2 on an
    oligodendrocyte myelin map, because that map names PC, PE, PI, PS and SM -- while
    Cholesterol metabolism, which would have shown 28 genes and 4 lipid classes for the same
    query, sat further down the list.

    So the score is the count of DISTINCT names the diagram would colour: genes the
    comparison called, plus lipid classes whose species it called. Both are what a reader
    sees. Only maps with a regulated lipid class or metabolite enter the menu.

    Counting distinct names matters: WP5304 draws LDLR six times, so counting nodes would
    rank a map by how often it repeats itself.
    """
    import re as _re

    store = _diagrams().get("diagrams") or {}
    synonyms = {k.upper(): [w.lower() for w in v]
                for k, v in _diagram_synonyms().items()}

    genes_called: set[str] = set()
    classes_called: set[str] = set()
    if contrast and cell_state:
        rows = _pathway_rows(ds, contrast, cell_state, modality, max_fdr, significance)
        for name, measure, _value in rows:
            if measure == "rna":
                genes_called.add(str(name).upper())
            elif measure == "metabolite":
                classes_called.add(str(name).upper())
            else:
                head = _re.match(r"^[A-Za-z0-9-]+", str(name))
                if head:
                    classes_called.add(head.group().upper())

    out: list[dict[str, Any]] = []
    for wp, diagram in store.items():
        gene_names, class_names = set(), set()
        for node in diagram["nodes"]:
            label = (node.get("label") or "").strip()
            if not label:
                continue
            if node.get("ensembl") or node.get("type") in ("GeneProduct", "Protein", "Rna"):
                if label.upper() in genes_called:
                    gene_names.add(label.upper())
            if node.get("type") == "Metabolite":
                if modality == "metabolite":
                    if label.upper() in classes_called:class_names.add(label.upper())
                    continue
                hit = _best_class(label, synonyms)
                if hit:
                    if hit in classes_called:
                        class_names.add(hit)
                else:
                    _term, members = _best_group(label)
                    for c in members:
                        if c in classes_called:
                            class_names.add(c)
        if not class_names:
            continue
        # HOW MUCH OF THE MAP IS CHEMISTRY.
        #
        # A metabolite share separates a metabolic pathway from a signalling or
        # copy-number map that happens to name a head group. Over the 68 lipid pathways it
        # runs from 0.010 (Opioid receptor pathways, 3 metabolites to 303 genes) to 1.000
        # (Metabolism overview, 165 to 0), with a median of 0.361.
        n_met = sum(1 for n in diagram["nodes"] if n.get("type") == "Metabolite")
        n_gene = sum(1 for n in diagram["nodes"]
                     if n.get("ensembl")
                     or n.get("type") in ("GeneProduct", "Protein", "Rna"))
        share = (n_met / (n_met + n_gene)) if (n_met + n_gene) else 0.0
        out.append({"id": wp, "name": diagram["name"],
                    "classes": diagram["classes"],
                    "n_nodes": len(diagram["nodes"]),
                    "n_metabolite_nodes": n_met,
                    "n_gene_nodes": n_gene,
                    "metabolite_share": round(share, 3),
                    "n_genes_painted": len(gene_names),
                    "n_lipids_painted": len(class_names),
                    "n_painted": len(gene_names) + len(class_names)})
    # BOTH MEASUREMENTS, NOT THE LARGER SUM.
    #
    # Ranking on genes plus lipids put "Orexin receptor pathway" first for AT2 in the ACDMPV
    # age comparison, on 49 genes and ONE lipid class, ahead of Cholesterol metabolism on 28
    # genes and 4. The atlas measures 202 lipid species in about a dozen head-group classes,
    # so a lipid is intrinsically scarcer than a gene and any plain sum drowns it. Nathan
    # asked for the pathway with "multiple lipids and genes", so the leading term is the
    # SCARCER of the two counts: a map has to show both before its total matters.
    # CHEMISTRY BEFORE VOLUME, ONCE BOTH MEASUREMENTS ARE PRESENT.
    #
    # Ranking on the scarcer count and then the total still opened Ciliated-axon in
    # Jaiswal IPF vs Healthy on "2q37 copy number variation syndrome": it ties Cholesterol
    # metabolism on 3 lipid classes and beats it on genes, 28 to 13. A copy-number map is
    # not a metabolic pathway, and a reader asking about lipids does not want one.
    #
    # The share is COARSENED TO ONE DECIMAL rather than cut at a threshold I would have to
    # invent. A median cut of 0.361 would have demoted Sphingolipid metabolism in
    # senescence, which sits at 0.311 and is exactly the kind of map wanted. Bucketing puts
    # 2q37 (0.173) in 0.2 and Cholesterol metabolism (0.375) in 0.4, so the metabolic map
    # wins, while inside a bucket the amount of data still decides.
    out.sort(key=lambda r: (-min(r["n_genes_painted"], r["n_lipids_painted"]),
                            -round(r["metabolite_share"], 1),
                            -r["n_painted"], r["n_nodes"], r["name"]))
    return out

def _discover_pathway_diagram(ds, wpid: str, cell_state: str = "", contrast: str = "",
                    max_fdr: float = 0.05, modality: str = "lipid", significance: str = "fdr") -> dict[str, Any]:
    """One published pathway, at its own coordinates, with the atlas laid onto it.

    Nathan, 2026-09-06: "Pathway visualization is supposed to be a pathway diagram to
    visualize, not a dot plot."

    So the geometry is WikiPathways own and the colour is ours. A gene node is coloured
    when the comparison called that gene; a metabolite node is coloured when the atlas
    measures lipids of that head group and the comparison called them. Everything the
    atlas cannot measure stays outlined and grey, so the picture reads as the published
    pathway with data on it, never as a figure the data invented.
    """
    import re as _re

    store = _diagrams().get("diagrams") or {}
    diagram = store.get(wpid)
    if diagram is None:
        return {"error": f"{wpid} is not one of the {len(store)} lipid pathways."}

    rows = []
    fold: dict[str, float] = {}
    if contrast and cell_state:
        rows = _pathway_rows(ds, contrast, cell_state, modality, max_fdr, significance)
        for name, measure, value in rows:
            key = measure + " " + str(name).upper()
            if key not in fold or abs(value) > abs(fold[key]):
                fold[key] = float(value)

    by_class: dict[str, list[float]] = {}
    # THE STRONGEST MEASURED SPECIES IN EACH CLASS, so a metabolite node can link to a page
    # that exists. A node is labelled with a head group -- "Galactosylceramide" -- and this
    # atlas measures species, never classes, so linking the label sent a reader to
    # /molecule/?gene=Galactosylceramide, which answers "not measured in either atlas".
    class_pick: dict[str, tuple[str, float]] = {}
    for name, measure, value in rows:
        if measure != "lipid":
            continue
        match = _re.match(r"^[A-Za-z0-9-]+", str(name))
        if match:
            head = match.group().upper()
            by_class.setdefault(head, []).append(float(value))
            # THE FEATURE'S OWN SPELLING, not the upper-cased lookup key. Taking the key
            # gave a link to CER(D18:1/22:0) for a feature the index stores as
            # Cer(d18:1/22:0); the page resolved it case-insensitively, but the address a
            # reader copies should be the name the atlas uses.
            if head not in class_pick or abs(value) > abs(class_pick[head][1]):
                class_pick[head] = (str(name), float(value))
    class_fold = {k: sum(v) / len(v) for k, v in by_class.items()}
    synonyms = {k.upper(): v for k, v in _diagram_synonyms().items()}
    measured_classes = set()
    for name in ds.features(modality):
        head = _re.match(r"^[A-Za-z0-9-]+", str(name))
        if head:
            measured_classes.add(head.group().upper())

    # GREY MUST NOT MEAN TWO DIFFERENT THINGS.
    #
    # Every uncoloured node was drawn grey and the key called it "not measured by this
    # atlas". On WP5304, Cholesterol metabolism, 61 of the 74 gene labels ARE measured
    # here; none of them reached FDR 0.05 in Adams2020 IPF vs Healthy in AT2. Calling
    # those 61 unmeasured is a false statement about the atlas, and it hides the useful
    # reading, which is that the pathway was covered and did not move. A node now carries
    # `in_atlas`, so the page can separate "tested, not significant" from "never measured".
    rna_features = ds.features("rna")
    upper_rna = {str(x).upper() for x in rna_features}

    metabolite_features={str(g).upper():str(g) for g in ds.features(modality)} if modality=="metabolite" else {}
    nodes, measured, testable = [], 0, 0
    for node in diagram["nodes"]:
        value, source, in_atlas, link = None, "", False, ""
        label = (node.get("label") or "").strip()
        if node.get("ensembl") or node.get("type") in ("GeneProduct", "Protein", "Rna"):
            in_atlas = label.upper() in upper_rna
            value = fold.get("rna " + label.upper())
            if value is not None:
                source = "gene"
            if in_atlas:
                link = label
        if node.get("type") == "Metabolite":
            # THE MOST SPECIFIC NAME WINS, NOT THE FIRST ONE FOUND.
            #
            # These are substring tests, and lipid class names nest: "galactosylceramide"
            # contains "ceramide" AND "lactosylceramide". Taking the first match in
            # dictionary order labelled a galactosylceramide node as Cer and linked it to
            # Cer(d18:1/22:0), when WP5224 itself declares the class HexCer. Matching on
            # the longest synonym resolves the nesting: galactosylceramide (18 characters)
            # beats lactosylceramide (16) and ceramide (8).
            hit = _best_class(label, synonyms)
            if hit:
                in_atlas = in_atlas or hit in measured_classes
                if value is None and hit in class_fold:
                    value = class_fold[hit]
                    source = "lipid class " + hit
                if hit in class_pick:
                    link = class_pick[hit][0]
            else:
                # NO CLASS NAMED, BUT PERHAPS A CATEGORY. "Phospholipid" is not one of the
                # 13 measured classes and is not nothing either.
                term, members = _best_group(label)
                present = [c for c in members if c in class_fold]
                if members:
                    in_atlas = in_atlas or any(c in measured_classes for c in members)
                if present:
                    value = sum(class_fold[c] for c in present) / len(present)
                    source = ("lipid group " + term.upper()
                              + " (" + ", ".join(present) + ")")
                    pick = [class_pick[c] for c in present if c in class_pick]
                    if pick:
                        link = max(pick, key=lambda t: abs(t[1]))[0]
        if node.get("type") == "Metabolite" and modality == "metabolite":
            # Exact common-name matching; never infer a lipid class for metabolites.
            link=metabolite_features.get(label.upper(), "")
            in_atlas=bool(link)
            value=fold.get("metabolite " + label.upper())
            source="metabolite" if value is not None else ""
        if value is not None:
            measured += 1
        if in_atlas:
            testable += 1
        nodes.append(dict(node, log2fc=value, measured_as=source, in_atlas=in_atlas,
                          link_feature=link))

    coloured = (f"{measured} of {len(nodes)} nodes carry a result from "
                f"{contrast} in "
                f"{cell_state}, using the stored differential calls without an additional significance filter: a gene by "
                f"its own fold change, a metabolite by "
                + (f"regulated species measured in its class. A further " if modality != "metabolite" else f"measured metabolite itself. A further ")
                + f"{testable - measured:,} of the {testable:,} nodes this atlas can measure "
                f"have no reported differential result, so they are outlined rather "
                f"than filled. Only the remaining {len(nodes) - testable:,} are outside "
                f"what the atlas measures at all."
                if contrast and cell_state
                else "Choose a comparison and a cell type to colour it.")
    return {
        "id": wpid, "name": diagram["name"], "classes": diagram["classes"],
        "width": diagram["width"], "height": diagram["height"],
        "nodes": nodes, "edges": diagram["edges"], "labels": diagram["labels"],
        "n_nodes": len(nodes), "n_measured": measured,
        "n_in_atlas": testable,
        "cell_state": cell_state, "contrast": contrast,
        "note": (f"{diagram['name']} drawn at WikiPathways own coordinates, "
                 f"{len(nodes)} nodes and {len(diagram['edges'])} interactions. "
                 + coloured),
    }


def availability(ds,contrast,state,modality):
    missing=[m for m in ('rna',modality) if not ds.available(contrast,m,state)]
    return {'available':not missing,'missing_modalities':missing,
            'note':required_note(modality) if missing else '', 'comparison':ds.comparison(contrast)}


def _marker_network(ds, *, cell_state, features, limit, gene_fold, min_expression, min_score):
    """Connect stored RNA markers to expressed TFs using state-specific model scores.

    Marker folds compare a state with the rest, not experimental comparison arms.
    Neither RNA nor edge differential calls are prerequisites for this view.
    """
    states = list(ds.states)
    state = cell_state or (states[0] if states else '')
    base = dict(source='marker', cell_state=state, available=True, nodes=[], edges=[],
                available_cell_states=states, comparison='', contrast='')
    if state not in states:
        return dict(base, available=False, note='Choose an available cell state.')
    rows = ds.marker_rows(state)
    cut = float(np.log2(gene_fold))
    markers = {r['gene']: r for r in rows if r.get('log2fc') is not None
               and np.isfinite(r['log2fc']) and r['log2fc'] > 0 and r['log2fc'] >= cut}
    requested = set(features or [])
    candidates = [g for g in markers if not requested or g in requested]
    candidates.sort(key=lambda g: (-markers[g]['log2fc'], g))
    if not candidates:
        return dict(base, status='no_markers', note=f'No stored upregulated RNA markers pass marker fold {gene_fold} in {state}.')
    try:
        model = ds.edge_scores(state)
    except (KeyError, FileNotFoundError) as exc:
        return dict(base, available=False, status='missing_edges', note=f'GRN edge scores are unavailable for {state}: {exc}')
    factors = {f for key in model for f in key.split('|')[:-1] if f}
    levels = ds.state_expression(state, factors)
    targets = set(candidates)
    links = {}
    for key, value in model.items():
        parts = key.split('|')
        if len(parts) < 2 or parts[-1] not in targets:
            continue
        score = float(value[0] if isinstance(value, (list, tuple)) else value)
        if not np.isfinite(score) or score == 0 or abs(score) < min_score:
            continue
        for tf in parts[:-1]:
            if not tf or tf == parts[-1] or levels.get(tf, 0) <= min_expression:
                continue
            pair = (tf, parts[-1])
            if pair not in links or abs(score) > abs(links[pair]):
                links[pair] = score
    connected = {target for _, target in links}
    selected = set([g for g in candidates if g in connected][:limit])
    links = {pair: score for pair, score in links.items() if pair[1] in selected}
    counts = {}
    for tf, _ in links:
        counts[tf] = counts.get(tf, 0) + 1
    regulators = set(sorted(counts, key=lambda tf: (-counts[tf], tf))[:limit])
    edges = [dict(source=tf, target=target, score=round(score, 4), log2fc=None, fdr=None)
             for (tf, target), score in sorted(links.items()) if tf in regulators]
    ids = {name for e in edges for name in (e['source'], e['target'])}
    nodes = [dict(id=g, kind='regulator' if g in regulators else 'target',
                  log2fc=markers.get(g, {}).get('log2fc'), fdr=markers.get(g, {}).get('fdr'),
                  expression=levels.get(g), n_targets=sum(e['source'] == g for e in edges))
             for g in sorted(ids)]
    scale = ds.state_expression_scale
    note = (f'{len(edges)} model edges connect expressed TFs to stored RNA markers in {state}. '
            f'Marker fold ≥{gene_fold}; |edge score| ≥{min_score}; TF expression >{min_expression} '
            f'({scale}). Red shows positive marker log2 fold versus other cell states; gray TFs are expressed regulators without a qualifying positive marker fold. '
            'Edge scores are predicted activity, not differential changes or activation/repression calls.')
    if not edges:
        note = f'No expressed TF connects to the selected RNA markers in {state} at these thresholds. ' + note
    return dict(base, nodes=nodes, edges=edges, status='ready' if edges else 'no_edges',
                note=note, interpretation=note, expression_scale=scale,
                n_features_asked=len(candidates), n_features_modelled=len(selected))


def network(ds,*,contrast='',cell_state='',features=None,limit=50,min_fold=1.2,
            gene_fold=1.2,max_fdr=0.05,min_expression=0.5,significance='reported',source='differential',min_score=0):
    if source == 'marker':
        return _marker_network(ds, cell_state=cell_state, features=features, limit=limit,
                               gene_fold=gene_fold, min_expression=min_expression, min_score=min_score)
    check=availability(ds,contrast,cell_state,'grn')
    if not check['available']:return dict(check,nodes=[],edges=[],status='missing_differentials')
    if not cell_state:return dict(check,nodes=[],edges=[],note='Choose a cell state.')
    targets=list(features or [])
    if not targets:
        cut=np.log2(gene_fold) if gene_fold>1 else 0
        rows=[r for r in ds.differentials(contrast,cell_state,'rna')
              if _passes(r,significance,max_fdr)
              and r.get('log2fc') is not None and abs(r['log2fc'])>=cut]
        candidates=[r['gene'] for r in sorted(rows,key=lambda r:(-abs(r['log2fc']),r['gene']))]
        # Rank only targets with a retained edge, so unconnected high-fold genes
        # or markers cannot consume Show's limit and hide a drawable network.
        eligible=_edge_stats(ds,contrast or ds.current_contrast,cell_state,max_fdr,min_fold,significance)
        model=ds.edge_scores(cell_state)
        keys=[key for key in eligible if key in model and '|' in key]
        factors={f for key in keys for f in key.split('|')[:-1]}
        case,control=ds.arm_expression(contrast or ds.current_contrast,cell_state,factors) if factors else ({},{})
        connected={key.split('|')[-1] for key in keys if any(
            (f not in case and f not in control) or max(case.get(f,0),control.get(f,0))>min_expression
            for f in key.split('|')[:-1])}
        targets=[gene for gene in candidates if gene in connected][:limit]
        if not targets:
            selection=f"gene fold {gene_fold}"
            note=(f"No selected RNA target has a retained edge in {cell_state} under "
                  f"{_gate(significance,max_fdr)}, {selection}, edge fold {min_fold}, "
                  f"and TF expression >{min_expression}. No unconnected nodes are shown.")
            return dict(check,nodes=[],edges=[],status='no_edges',note=note,
                        interpretation=note,cell_state=cell_state,contrast=contrast)
    result=_discover_network(ds,cell_state,targets,limit=limit,contrast=contrast or ds.current_contrast,
                             max_fdr=max_fdr,min_fold=min_fold,min_expression=min_expression,significance=significance)
    note=result.get('note','')
    result.update(check,source=source,significance=significance)
    result['note']=note
    # Preserve the algorithm's explanatory note after the availability check.
    result['note']=result.get('note') or f"Edges pass {_gate(significance,max_fdr)}, fold {min_fold}; TF expression exceeds {min_expression} in either arm."
    from .integration_interpretation import interpret_network
    result['interpretation']=interpret_network(result)
    return result


def pathways(ds,*,contrast='',cell_state='',modality='lipid',max_fdr=0.05,significance='fdr',source='differential',features=None):
    check=availability(ds,contrast,cell_state,modality)
    if not check['available']:return dict(check,pathways=[],status='missing_differentials')
    if source=='marker' or features:
        selected=set(features or ds.marker_features(cell_state,'rna',50))
        original=ds
        class Selected:
            def __getattr__(self,key):return getattr(original,key)
            def differentials(self,c,s,m):
                rows=original.differentials(c,s,m)
                return [r for r in rows if r['gene'] in selected] if m=='rna' else rows
        ds=Selected()
    ranks=_discover_pathway_ranking(ds,cell_state,contrast or ds.current_contrast,max_fdr,modality,significance)
    return dict(check,pathways=ranks,note='Only pathways containing regulated lipid or metabolite features are listed; stored differential calls are used without an additional significance filter.' if ranks else 'No pathway contains a regulated lipid or metabolite feature for this comparison and cell state.')


def pathway(ds,*,id,contrast='',cell_state='',modality='lipid',max_fdr=0.05,significance='fdr'):
    check=availability(ds,contrast,cell_state,modality)
    if not check['available']:return dict(check,nodes=[],edges=[],status='missing_differentials')
    result=_discover_pathway_diagram(ds,id,cell_state,contrast or ds.current_contrast,max_fdr,modality,significance)
    if result.get('error'):return result
    if not any(n.get('log2fc') is not None and str(n.get('measured_as','')).startswith(('lipid','metabolite')) for n in result.get('nodes',[])):
        return dict(check,nodes=[],edges=[],status='no_regulated_metabolites',note='This pathway has no regulated lipid or metabolite feature in the selected comparison and cell state.')
    from .integration_interpretation import interpret_pathway
    result['interpretation']=interpret_pathway(result)
    return result
