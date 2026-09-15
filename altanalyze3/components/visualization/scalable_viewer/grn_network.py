"""The regulatory network a bundle's own differentials imply, and per-factor activity.

WHAT THIS IMPLEMENTS, AND WHOSE RULES THEY ARE

Nathan asked the scALABLE viewer's Chat to draw a GRN network by "the precise rules from
the Network option in the discover LungMAP.net viewer". Those rules live in
`refactored_website/app/lungmap/molecule/discover.py:1255 regulator_network()`, which
records his own instruction of 2026-09-06:

    "the network should restrict TF target genes that are transcriptionally regulated (up
    or down) in the selected comparison and cell type, specifically those features in the
    table selection based on the current sort and the 'Show' number of features. The TFs do
    not need to be regulated themselves since their activity but not expression might be
    changed."

So the TARGETS come from the differential and the EDGES come from the regulatory model.

That docstring also records the three faults of its own first version, all of which the
viewer's existing `/api/jobs/{id}/grn/network` still has
(`cellHarmony/webapp/app.py:2409 _build_grn_network_payload`):

  1. it draws every edge above an ABSOLUTE `threshold` somebody chose. The floor here is a
     PERCENTILE of the edge scores in this cell state, so it adapts to the state.
  2. it draws factors the cell type does not express. A factor here must clear a percentile
     of the expression of every modelled factor in this cell state.
  3. it paints every factor the same colour -- literally `log2fc = 1.0` for a TF and -1.0
     for a target. Each factor here carries its OWN differential expression and its OWN
     differential activity, so the drawing can colour by one and hover the other.

WHY THIS READS THE BUNDLE AND NOT THE LUNGMAP INDEX

`discover.regulator_network` reads `grn_edges_by_state.json`, a 111 MB file under
`CELLREF2_FAST` that only the LungMAP site deploys. A bundle already carries the same two
quantities: the GRN modality's `stats_mean` is features by cell state, and the RNA store's
`stats_mean` is genes by cell state. So the viewer stays self-contained and a bundle on
another machine draws the same network.

ONE HONEST DIFFERENCE FROM THE SITE, STATED IN THE PAYLOAD

The site measures factor expression "over THE DONORS OF THE COMPARISON'S TWO ARMS". A
bundle's `stats_mean` is the mean over the metacells of the cell state, and a metacell
pools donors, so this measures over every metacell of the state. `arm_scope` in the payload
says which of the two happened, the way the site's own fallback string does. Restricting to
the arms needs the contrast's donor keys, which a bundle does not carry.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# The site's Network tab ships these two as "strongest half" and "above the median here".
DEFAULT_EDGE_PERCENTILE = 50.0
DEFAULT_EXPRESSION_PERCENTILE = 50.0
# Nathan, 2026-09-07: "default should be 200 but user should be allowed to increase or
# decrease based on chat". The site's route caps its feature list at 200; its `limit`
# default of 50 is the Show box, not the cap.
DEFAULT_LIMIT = 200
MAX_LIMIT = 2000
# The FDR the site ranks within: "Ranked by absolute log2 fold change among results at or
# below FDR 0.05."
TARGET_FDR = 0.05


def _state_column(ds, cell_state: str) -> int:
    try:
        return list(ds.states).index(cell_state)
    except ValueError as exc:
        raise KeyError(f"{cell_state!r} is not a cell state of {ds.id}") from exc


def _split_edge(name: str) -> Optional[Tuple[List[str], str]]:
    """`A|B` is one regulator on target B. `A|B|C` is a composite element of A and B on C.

    Every field but the last is a REGULATOR; the last is the target. Reading a middle field
    as a target would turn thousands of target genes into factors.
    """
    parts = [p for p in str(name).split("|") if p]
    if len(parts) < 2:
        return None
    return parts[:-1], parts[-1]


def _deg_rows(ds, comp_id: str, state: str, fdr_max: Optional[float],
              max_rows: int) -> List[Dict[str, Any]]:
    """One differential's rows for one cell state, or [] when the bundle lacks it."""
    try:
        table = ds.deg_table(comp_id, max_rows=max_rows, fdr_max=fdr_max, state=state)
    except (KeyError, FileNotFoundError, OSError):
        return []
    return list(table.get("rows") or [])


def _sibling_comparison(ds, comp_id: str, modality: str) -> str:
    """The same comparison, same kind, computed on another modality.

    A manifest id is `<comparison>::<kind>` for RNA and `<modality>::<comparison>::<kind>`
    for everything else, so the RNA id is the tail of every other id. Matching on that tail
    keeps a TF-activity contrast tied to the RNA contrast a reader chose, instead of
    guessing by position.
    """
    tail = comp_id.split("::", 1)[1] if "::" in comp_id and not comp_id.startswith(
        tuple(m + "::" for m in ds.modality_manifest())) else comp_id
    for c in ds.deg_manifest().get("comparisons", []):
        cid = str(c.get("id") or "")
        if str(c.get("modality") or "rna") != modality:
            continue
        if cid == f"{modality}::{tail}" or cid.endswith("::" + tail):
            return cid
    return ""


def _rna_contrast_tail(comp_id: str) -> str:
    """The `<contrast>::<comparison>::<kind>` tail shared by every modality's id."""
    parts = comp_id.split("::")
    if len(parts) >= 4:          # <modality>::<contrast>::<comparison>::<kind>
        return "::".join(parts[1:])
    return comp_id


def regulator_network(
    ds, cell_state: str, *,
    contrast: str = "",
    features: Optional[List[str]] = None,
    limit: int = DEFAULT_LIMIT,
    edge_percentile: float = DEFAULT_EDGE_PERCENTILE,
    expression_percentile: float = DEFAULT_EXPRESSION_PERCENTILE,
    modality: str = "grn",
) -> Dict[str, Any]:
    """The factors that regulate the features a differential put on screen.

    `features` overrides the differential, which is how Chat can name its own targets.
    Without it, the targets are the top `limit` rows of `contrast` inside `cell_state`,
    ranked by absolute log2 fold change among rows at or below FDR 0.05.
    """
    if not cell_state:
        return {"nodes": [], "edges": [],
                "note": "Choose a cell type: a regulatory edge is scored within one."}
    try:
        col = _state_column(ds, cell_state)
    except KeyError as exc:
        return {"nodes": [], "edges": [], "note": str(exc)}

    try:
        store = ds.modality(modality)
    except (KeyError, FileNotFoundError) as exc:
        return {"nodes": [], "edges": [],
                "note": f"this bundle carries no '{modality}' modality: {exc}"}

    names = list(store.features)
    mean = store.stats_mean
    if mean.ndim != 2 or mean.shape[1] <= col:
        return {"nodes": [], "edges": [],
                "note": f"the '{modality}' store holds no column for {cell_state}"}
    scores = np.asarray(mean[:, col], dtype=np.float64)

    limit = max(1, min(int(limit or DEFAULT_LIMIT), MAX_LIMIT))

    # ---- 1. THE TARGETS. From the differential unless the caller names them.
    target_source: str
    target_rows: Dict[str, Dict[str, Any]] = {}
    if features:
        wanted = [str(f).strip() for f in features if str(f).strip()][:limit]
        target_source = f"{len(wanted)} features named by the caller"
    else:
        if not contrast:
            return {"nodes": [], "edges": [],
                    "note": ("Name a comparison or a feature list: the network needs "
                             "targets, and it takes them from a differential.")}
        rows = _deg_rows(ds, contrast, cell_state, TARGET_FDR, 200000)
        # THE SITE'S OWN RANKING. deg_table sorts by FDR first; the Network tab ranks by
        # absolute fold change among the rows that already passed the FDR, so re-sort.
        rows.sort(key=lambda r: -abs(r.get("log2fc") or 0.0))
        wanted = []
        for r in rows:
            g = str(r.get("gene") or "")
            if not g or g in target_rows:
                continue
            target_rows[g] = r
            wanted.append(g)
            if len(wanted) >= limit:
                break
        target_source = (f"the top {len(wanted)} of {len(rows)} features of {contrast} in "
                         f"{cell_state} at FDR {TARGET_FDR}, by absolute log2 fold change")
        if not wanted:
            return {"nodes": [], "edges": [], "cell_state": cell_state,
                    "contrast": contrast, "limit": limit,
                    "note": (f"{contrast} holds no feature in {cell_state} at FDR "
                             f"{TARGET_FDR}, so there is no target to place.")}

    # ---- 2. THE EDGE FLOOR, as a percentile of THIS cell state's own edges.
    finite = scores[np.isfinite(scores)]
    pct = min(99.0, max(0.0, float(edge_percentile)))
    edge_cut = float(np.percentile(finite, pct)) if finite.size else 0.0

    target_set = set(wanted)
    hits: Dict[str, Dict[str, float]] = {}
    n_edges_into_targets = 0
    for j, name in enumerate(names):
        split = _split_edge(name)
        if split is None:
            continue
        regulators, target = split
        if target not in target_set:
            continue
        score = float(scores[j])
        if not np.isfinite(score):
            continue
        n_edges_into_targets += 1
        if score < edge_cut:
            continue
        for reg in regulators:
            slot = hits.setdefault(reg, {})
            if target not in slot or score > slot[target]:
                slot[target] = score

    base = {"cell_state": cell_state, "contrast": contrast, "limit": limit,
            "modality": modality, "target_source": target_source,
            "n_targets_requested": len(wanted),
            "n_edges_into_targets": n_edges_into_targets,
            "edge_percentile": pct, "edge_cut": round(edge_cut, 4),
            "expression_percentile": min(99.0, max(0.0, float(expression_percentile)))}
    if not hits:
        return {**base, "nodes": [], "edges": [], "n_regulators_total": 0,
                "note": (f"No edge into these {len(wanted)} features reaches the "
                         f"{pct:.0f}th percentile of {cell_state}'s edge scores, "
                         f"{edge_cut:.3f}. Lower the edge weight to see weaker edges.")}
    n_reg_before_expression = len(hits)

    # ---- 3. THE EXPRESSION FLOOR, over every modelled factor in this cell state.
    every_factor = sorted({r for name in names
                           for r in (_split_edge(name) or ([], ""))[0] if r})
    rna_mean = ds.stats_mean
    levels: Dict[str, float] = {}
    for f in every_factor:
        row = ds.resolve_gene(f)
        if row is None:
            continue
        v = float(rna_mean[row, col])
        if np.isfinite(v):
            levels[f] = v
    expr_pct = base["expression_percentile"]
    expression_cut = (float(np.percentile(np.array(list(levels.values())), expr_pct))
                      if levels else 0.0)
    silent = sorted(r for r in hits if r in levels and levels[r] < expression_cut)
    unmeasured = sorted(r for r in hits if r not in levels)
    for r in silent:
        hits.pop(r, None)

    base.update({
        "expression_cut": round(expression_cut, 4),
        "n_factors_modelled": len(every_factor),
        "n_factors_measured": len(levels),
        "n_regulators_before_expression_floor": n_reg_before_expression,
        "n_regulators_dropped_as_silent": len(silent),
        "regulators_dropped_as_silent": silent[:40],
        "n_regulators_unmeasured": len(unmeasured),
        "regulators_unmeasured": unmeasured[:40],
        # THE SITE MEASURES OVER THE ARMS' DONORS; A BUNDLE CANNOT. Said plainly rather
        # than left for a reader to assume.
        "arm_scope": (f"every metacell of {cell_state} in this bundle, because a bundle "
                      f"carries no donor list for a contrast's two arms"),
    })
    if not hits:
        return {**base, "nodes": [], "edges": [],
                "note": (f"Every factor reaching these features falls below the "
                         f"{expr_pct:.0f}th percentile of factor expression in "
                         f"{cell_state}, {expression_cut:.3f}. Lower the factor "
                         f"expression floor to include quieter factors.")}

    # ---- 4. PER-FACTOR DIFFERENTIALS: activity AND expression, kept apart.
    activity: Dict[str, Dict[str, Any]] = {}
    expression: Dict[str, Dict[str, Any]] = {}
    tf_contrast = ""
    if contrast:
        tf_contrast = _sibling_comparison(ds, contrast, "grn_tf")
        if tf_contrast:
            for r in _deg_rows(ds, tf_contrast, cell_state, None, 200000):
                g = str(r.get("gene") or "")
                if g:
                    activity[g] = r
        rna_contrast = (contrast if str(contrast).count("::") <= 2
                        else _sibling_comparison(ds, contrast, "rna"))
        for r in _deg_rows(ds, rna_contrast or contrast, cell_state, None, 200000):
            g = str(r.get("gene") or "")
            if g:
                expression[g] = r
    base["tf_activity_contrast"] = tf_contrast
    base["n_factors_with_tested_activity"] = len(activity)

    nodes: List[Dict[str, Any]] = []
    edges: List[Dict[str, Any]] = []
    seen_targets: set[str] = set()
    for reg in sorted(hits):
        act = activity.get(reg) or {}
        exp = expression.get(reg) or {}
        # ABSENCE OF A ROW MEANS THE FACTOR DID NOT PASS, NOT THAT NOBODY LOOKED.
        # DEG_detailed holds only the rows that cleared the run's gates.
        nodes.append({
            "id": reg, "label": reg, "role": "factor", "shape": "diamond",
            "n_targets_here": len(hits[reg]),
            "expression_level": round(levels[reg], 4) if reg in levels else None,
            "expression_log2fc": exp.get("log2fc"), "expression_fdr": exp.get("fdr"),
            "activity_log2fc": act.get("log2fc"), "activity_fdr": act.get("fdr"),
            # WHICH NUMBER COLOURS THIS NODE, IN THREE STATES, NOT TWO.
            #
            # The site's legend reads "dashed rim = coloured by activity, not expression"
            # and greys a node with "no significant result here". So a factor whose
            # ACTIVITY moved while its expression did not gets the activity colour and the
            # dashed rim -- which is the whole reason Nathan asked for activity, since a
            # factor can act without changing its own transcript. A factor with neither row
            # is grey, and calling that "expression" would paint it as a measured zero.
            "colour_by": ("expression" if exp.get("log2fc") is not None
                          else "activity" if act.get("log2fc") is not None
                          else "none"),
            "tested_activity": bool(act), "tested_expression": bool(exp),
        })
        for tgt, score in sorted(hits[reg].items(), key=lambda kv: -kv[1]):
            if tgt not in seen_targets:
                seen_targets.add(tgt)
                trow = target_rows.get(tgt) or expression.get(tgt) or {}
                nodes.append({
                    "id": tgt, "label": tgt, "role": "target", "shape": "ellipse",
                    "expression_log2fc": trow.get("log2fc"),
                    "expression_fdr": trow.get("fdr"),
                    # A target came from the differential, so it has a fold change. When it
                    # was named by the caller instead, it may have none, and grey is then
                    # the honest colour.
                    "colour_by": ("expression" if trow.get("log2fc") is not None
                                  else "none"),
                    "tested_expression": bool(trow), "tested_activity": False,
                })
            edges.append({"source": reg, "target": tgt,
                          "score": round(float(score), 4)})

    scores_drawn = [e["score"] for e in edges]
    return {**base,
            "nodes": nodes, "edges": edges,
            "n_regulators": sum(1 for n in nodes if n["role"] == "factor"),
            "n_targets_reached": len(seen_targets),
            "n_edges_drawn": len(edges),
            "edge_score_range": [round(min(scores_drawn), 4),
                                 round(max(scores_drawn), 4)] if scores_drawn else None,
            "legend": ("diamond = factor, circle = regulated feature, dashed rim = "
                       "coloured by activity not expression, line width = edge score"),
            "note": ""}


def tf_activity_profile(ds, *, cell_state: str = "", contrast: str = "",
                        factors: Optional[List[str]] = None,
                        limit: int = 25) -> Dict[str, Any]:
    """Per-factor activity for one cell state, and its tested change when one exists.

    WHY THIS IS SEPARATE FROM THE NETWORK. Nathan, 2026-09-07: "Both GRN and TF-activity
    should be able to have separate types of plots supported in Chat." A network answers
    "who regulates these genes". This answers "which factors are most active here, and
    which of them moved", which is a ranked list, not a graph.

    The activity comes from the bundle's own `grn_tf` store, which is per METACELL, so a
    cell state's number is a mean over its metacells rather than a single broadcast value.
    """
    out: Dict[str, Any] = {"cell_state": cell_state, "contrast": contrast}
    try:
        store = ds.modality("grn_tf")
    except (KeyError, FileNotFoundError) as exc:
        return {**out, "rows": [],
                "note": f"this bundle carries no 'grn_tf' modality: {exc}"}
    out["store_kind"] = store.kind
    out["modality_label"] = store.label

    names = list(store.features)
    if cell_state:
        try:
            col = _state_column(ds, cell_state)
        except KeyError as exc:
            return {**out, "rows": [], "note": str(exc)}
        level = np.asarray(store.stats_mean[:, col], dtype=np.float64)
    else:
        level = np.asarray(store.stats_mean, dtype=np.float64).mean(axis=1)
        out["cell_state"] = "every cell state, averaged"

    tested: Dict[str, Dict[str, Any]] = {}
    tf_contrast = ""
    if contrast:
        tf_contrast = (contrast if contrast.startswith("grn_tf::")
                       else _sibling_comparison(ds, contrast, "grn_tf"))
        if tf_contrast:
            for r in _deg_rows(ds, tf_contrast, cell_state or None, None, 200000):
                g = str(r.get("gene") or "")
                if g:
                    tested[g] = r
    out["tf_activity_contrast"] = tf_contrast
    out["n_factors_with_tested_change"] = len(tested)

    want = ([f for f in (factors or []) if f] or None)
    rows: List[Dict[str, Any]] = []
    for i, name in enumerate(names):
        if want is not None and name not in want:
            continue
        t = tested.get(name) or {}
        rows.append({"factor": name, "activity": round(float(level[i]), 4),
                     "log2fc": t.get("log2fc"), "fdr": t.get("fdr"),
                     "tested": bool(t)})
    # A tested factor first, biggest move first; then the quiet ones by activity.
    rows.sort(key=lambda r: (0 if r["tested"] else 1,
                             -abs(r["log2fc"] or 0.0), -r["activity"]))
    out["n_factors"] = len(rows)
    out["rows"] = rows[:max(1, int(limit or 25))]
    out["statistic"] = ("the sum of each factor's predicted edge activity, averaged over "
                        "the metacells of the cell state")
    out["absence"] = ("a factor with no row in the contrast was tested and did not pass "
                      "that run's fold and FDR gates")
    out["note"] = ""
    return out
