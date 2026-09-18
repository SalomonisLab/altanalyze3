"""Shared regulatory networks and factor profiles for bundles and uploaded jobs.

Data adapters provide state means and completed differential tables. Targets are
selected from significant differential results (FDR <= 0.05); edge and factor
expression floors are state-specific percentiles. Factor activity and RNA
expression changes stay separate. No analysis is run by these view functions.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import re

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

    Match explicit manifest identities. Bundle ids and uploaded run ids use
    different formats; neither list order nor a partial id identifies the groups.
    """
    comparisons = ds.deg_manifest().get("comparisons", [])
    selected = next((c for c in comparisons if c.get("id") == comp_id), None)
    if selected is None:
        return ""
    for candidate in comparisons:
        if str(candidate.get("modality") or "rna") != modality:
            continue
        # Compare the contrast identity, never list position or a partial suffix.
        if all(candidate.get(key) == selected.get(key)
               for key in ("comparison", "contrast", "kind")):
            return str(candidate["id"])
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
        selected = next((c for c in ds.deg_manifest().get("comparisons", []) if c.get("id") == contrast), {})
        target_contrast = contrast
        if selected.get("modality") == "grn_tf":
            target_contrast = _sibling_comparison(ds, contrast, "rna")
        rows = _deg_rows(ds, target_contrast, cell_state, TARGET_FDR, 200000) if target_contrast else []
        # THE SITE'S OWN RANKING. deg_table sorts by FDR first; the Network tab ranks by
        # absolute fold change among the rows that already passed the FDR, so re-sort.
        rows.sort(key=lambda r: -abs(r.get("log2fc") or 0.0))
        wanted = []
        for r in rows:
            g = str(r.get("gene") or "")
            split = _split_edge(g)
            if split:
                g = split[1]
            if not g or g in wanted:
                continue
            # An edge fold change is not target RNA expression.
            if not split:
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
        "arm_scope": getattr(ds, "arm_scope", f"every metacell of {cell_state} in this bundle; not restricted to comparison arms"),
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
    if contrast and cell_state:
        tf_contrast = _sibling_comparison(ds, contrast, "grn_tf")
        if tf_contrast:
            for r in _deg_rows(ds, tf_contrast, cell_state, None, 200000):
                g = str(r.get("gene") or "")
                if g:
                    activity[g] = r
        rna_contrast = _sibling_comparison(ds, contrast, "rna")
        for r in _deg_rows(ds, rna_contrast, cell_state, None, 200000) if rna_contrast else []:
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
    if factors and not cell_state:
        # A named factor without a state asks WHERE its activity is highest.
        # Do not collapse the state axis into a single global average.
        lookup = {str(name).casefold(): (i, str(name)) for i, name in enumerate(names)}
        selected = list(dict.fromkeys(str(f).casefold() for f in factors if f))
        rows = []
        for factor in selected:
            if factor not in lookup:
                continue
            i, name = lookup[factor]
            for col, state in enumerate(ds.states):
                value = float(store.stats_mean[i, col])
                if np.isfinite(value):
                    rows.append(dict(factor=name, cell_state=str(state), activity=value,
                                     log2fc=None, fdr=None, tested=False,
                                     label=f"{name} — {state}"))
        rows.sort(key=lambda r: (-r['activity'], r['cell_state'], r['factor']))
        return dict(out, by_cell_state=True, cell_state='', contrast='',
                    rows=rows[:max(1, min(int(limit or 25), MAX_LIMIT))],
                    n_states=len({r['cell_state'] for r in rows}), n_results=len(rows),
                    n_factors=len({r['factor'] for r in rows}),
                    statistic=getattr(store, 'statistic', 'mean stored TF activity per cell state'),
                    missing_factors=[str(f) for f in factors if str(f).casefold() not in lookup],
                    note='' if rows else 'No finite activity values were found for the requested TFs in the stored TF-activity modality.')

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
    if contrast and cell_state:
        tf_contrast = _sibling_comparison(ds, contrast, "grn_tf")
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
    out["rows"] = rows[:max(1, min(int(limit or 25), MAX_LIMIT))]
    out["statistic"] = getattr(store, "statistic", "sum of predicted outgoing edge activity, averaged over the metacells of the cell state")
    out["absence"] = "No reported differential row does not establish that a factor was tested or unchanged."
    out["note"] = ""
    return out


def read_regulatory_question(question, states, genes):
    """Local, data-name-constrained routing shared by both scALABLE deployments.

    Returns None for other questions, leaving the configured AI router in charge.
    No expression values or statistical results are inferred by this parser.
    """
    text = question.lower()
    if not re.search(r"\b(tf|tfs|transcription factors?|regulators?|regulatory|regulon|grn)\b", text):
        return None
    network = bool(re.search(r"network|targets?|edges?|regulatory", text))
    differential = bool(re.search(r"differential|changed|changes|significant|upregulated|downregulated", text))
    matched = [state for state in sorted(states, key=len, reverse=True)
               if re.search(r"(?<![\w-])" + re.escape(state.lower()) + r"(?![\w-])", text)]
    index = {str(g).lower(): str(g) for g in genes}
    named_genes = [index[token.lower()] for token in re.findall(r"[A-Za-z0-9_.-]{3,}", question)
                   if token.lower() in index]
    limit = re.search(r"(?:top|show)\s+(\d+)", text)
    return {"intent": "differential" if differential else ("regulatory_driver" if network else "tf_activity"),
            "modality": "grn" if network else "grn_tf", "cell_state": matched[0] if matched else "",
            "genes": list(dict.fromkeys(named_genes)),
            "limit": max(1, min(int(limit.group(1)), MAX_LIMIT)) if limit else (DEFAULT_LIMIT if network else MAX_LIMIT if named_genes and not matched else 25),
            "router": "local_regulatory"}


def tf_activity_state_chat(answer):
    """Shared answer/table/plot contract for named-factor cell-state rankings."""
    rows = answer.get('rows') or []
    columns = ['cell_state', 'factor', 'activity']
    if answer.get('n_factors', 0) > 1: columns.append('label')
    groups = {}
    for row in rows:
        groups.setdefault(row['factor'], []).append(row)
    leaders = '; '.join(f"{factor}: " + ', '.join(
        f"{r['cell_state']} ({r['activity']:.4g})" for r in values[:5]) for factor, values in groups.items())
    text = (f"Highest mean TF activity by cell type — {leaders}. "
            f"Ranked across {answer.get('n_states', 0)} cell states using {answer.get('statistic')}. "
            'These are stored activity levels; this ranking does not require a differential analysis.') if rows else answer.get('note')
    if answer.get('missing_factors'):
        text += ' TFs absent from the activity modality: ' + ', '.join(answer['missing_factors']) + '.'
    return dict(answer=text, status='ok' if rows else 'not_found',
                table=dict(columns=columns, rows=rows),
                column_labels={'cell_state':'Cell type / state', 'activity':'Mean TF activity'},
                result_controls=dict(sort_by='activity', sort_direction='desc'),
                plot=dict(kind='barchart', label_column='cell_state' if len(groups)==1 else 'label',
                          value_column='activity', filterable=True, title='TF activity by cell type / state'))
