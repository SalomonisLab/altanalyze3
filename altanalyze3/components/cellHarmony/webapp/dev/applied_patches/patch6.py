P = "/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py"
src = open(P).read()

def sub(old, new):
    global src
    assert src.count(old) == 1, f"{src.count(old)} matches for:\n{old[:160]}"
    src = src.replace(old, new)

# --- imports ------------------------------------------------------------------
sub("""import io
import json
import logging
import re
import shutil
import threading""",
    """import io
import json
import logging
import os
import re
import shutil
import threading
import time
import urllib.request""")

# --- request model ------------------------------------------------------------
sub("""class ClientLogRequest(BaseModel):""",
    """class ChatRequest(BaseModel):
    question: str = ""


class ClientLogRequest(BaseModel):""")

# --- the chat machinery -------------------------------------------------------
CHAT = '''#: Where a chat question is read into one supported query. The model runs in the
#: LungMAP site process; this app holds none of its own. Override with
#: CELLHARMONY_ASSISTANT_URL when the site is not on the same host.
CHAT_ASSISTANT_URL = os.environ.get(
    "CELLHARMONY_ASSISTANT_URL", "http://127.0.0.1:8001/api/assistant/viewer-intent")

#: The router speaks protocol names. These are answered by the executors below.
_CHAT_PROTOCOL_ALIAS = {
    "cell_identity": "markers",
    "state_comparison": "compare",
    "expression_lookup": "expression",
    "state_contrast": "differential",
    "donor_heterogeneity": "differential",
    "patient_stratification": "differential",
    "most_affected_state": "differential",
    "shared_vs_state_specific": "differential",
    "contrast_specificity": "differential",
    "pathway_program": "goelite",
    "regulatory_driver": "network",
    "communication_rewiring": "ccc",
}

#: Protocols with a specification but no executor here. Naming the gap is the
#: point: answering one of these with a neighbouring analysis would report a
#: different statistic under the question the user asked.
_CHAT_NOT_YET = {
    "severity_gradient": "a Spearman correlation of per-donor pseudobulk against a clinical variable",
    "dose_response": "a monotonic trend test across ordered stages",
    "composition_shift": "per-donor cell-count composition",
    "coexpression_module": "per-donor co-expression around a seed gene",
    "annotation_concordance": "a cross-tabulation of the two annotations",
}


def _chat_states_by_size(cache: Dict[str, Any]) -> List[tuple]:
    """(state, n cells), largest first. Used to name states in the examples."""
    values, counts = np.unique(np.asarray(cache["populations"], dtype=str), return_counts=True)
    order = np.argsort(-counts)
    return [(str(values[i]), int(counts[i])) for i in order]


def _chat_contrast(meta: Dict) -> Dict[str, str]:
    """The one comparison this job has run, or empty when none has."""
    differential = meta.get("differential") or {}
    if str(differential.get("status") or "").lower() != "completed":
        return {}
    config = differential.get("config") or {}
    case = str(config.get("case_label") or "Group 1").strip()
    control = str(config.get("control_label") or "Group 2").strip()
    return {"case": case, "control": control, "label": f"{case} versus {control}"}


def _chat_marker_table(meta: Dict, modality: str = "rna") -> pd.DataFrame:
    """The marker table cellHarmony wrote, or an empty frame."""
    marker_analysis = _modality_marker_analysis(meta, modality) or {}
    path = Path(str(marker_analysis.get("markers_tsv", "")).strip())
    if not path.exists():
        return pd.DataFrame()
    frame = pd.read_csv(path, sep="\\t")
    if "cluster" not in frame.columns or "Gene" not in frame.columns:
        return pd.DataFrame()
    frame["cluster"] = frame["cluster"].astype(str)
    frame["Gene"] = frame["Gene"].astype(str)
    return frame


def _chat_computed_markers(cache: Dict[str, Any], state: str, limit: int) -> List[Dict[str, Any]]:
    """One-versus-rest markers for a state the marker table does not cover.

    The gap between a gene's mean inside the state and its mean everywhere else.
    Labelled `computed one-vs-rest` in the answer, because a stored fold change
    and this gap are not the same statistic.
    """
    adata = cache["adata"]
    values_of = np.asarray(cache["populations"], dtype=str)
    inside = values_of == str(state)
    if not inside.any():
        return []
    X = adata.X
    rows = np.nonzero(inside)[0]
    indicator = sp.csr_matrix((np.ones(rows.size, dtype=np.float64),
                               (np.zeros(rows.size, dtype=np.int64), rows)),
                              shape=(1, X.shape[0]))
    sums = indicator @ X
    sums = np.asarray(sums.todense() if sp.issparse(sums) else sums, dtype=np.float64).ravel()
    total = np.asarray(X.sum(axis=0), dtype=np.float64).ravel()
    n_inside = float(rows.size)
    gap = (sums / max(n_inside, 1.0)) - ((total - sums) / max(float(X.shape[0]) - n_inside, 1.0))
    var_names = [str(v) for v in cache["var_names"]]
    out = []
    for row in np.argsort(-gap)[:limit]:
        if gap[int(row)] <= 0:
            break
        out.append({"gene": var_names[int(row)], "cluster": str(state),
                    "fold": round(float(gap[int(row)]), 4), "p": None,
                    "source": "computed one-vs-rest"})
    return out


def _chat_markers_for_states(app: FastAPI, meta: Dict, cache: Dict[str, Any],
                             states: List[str], limit: int = 25) -> List[Dict[str, Any]]:
    """Marker genes for these states, from the marker table where it covers them."""
    wanted = [str(s) for s in states if s]
    if not wanted:
        return []
    frame = _chat_marker_table(meta)
    rows: List[Dict[str, Any]] = []
    covered = set()
    if not frame.empty:
        subset = frame.loc[frame["cluster"].isin(wanted)].copy()
        if "FDR p-value" in subset.columns:
            subset["_p"] = pd.to_numeric(subset["FDR p-value"], errors="coerce")
        else:
            subset["_p"] = np.nan
        subset["_fold"] = pd.to_numeric(subset.get("Fold"), errors="coerce")
        subset = subset.sort_values(["_p", "_fold"], ascending=[True, False])
        for row in subset.itertuples():
            rows.append({"gene": str(row.Gene), "cluster": str(row.cluster),
                         "fold": float(row._fold) if _is_finite_number(row._fold) else 0.0,
                         "p": float(row._p) if _is_finite_number(row._p) else None,
                         "source": "marker table"})
            covered.add(str(row.cluster))
    per_state = max(1, limit // max(1, len(wanted)))
    for state in wanted:
        if state not in covered:
            rows.extend(_chat_computed_markers(cache, state, per_state))
    return rows[:limit]


def _chat_examples(app: FastAPI, meta: Dict) -> Dict[str, Any]:
    """Example questions built from this job's own reference and cell states.

    A bone-marrow job gets bone-marrow states, a lung job gets lung states, and
    a comparison is only offered when the job has actually run one. An example
    naming a cell state the dataset does not hold would fail the moment it was
    clicked.
    """
    reference_id = str(meta.get("reference") or "")
    tissue = ""
    if "lung" in reference_id.lower():
        tissue = "lung"
    elif "_bm" in reference_id.lower() or "marrow" in reference_id.lower():
        tissue = "bone marrow"
    try:
        reference_label = str(_reference_entry_for_meta(meta).get("label") or reference_id)
    except Exception:  # noqa: BLE001 - the registry may not hold this reference any more
        reference_label = reference_id

    try:
        cache = _get_expression_cache(app, meta)
    except Exception:  # noqa: BLE001 - a job whose h5ad is gone still opens the tab
        return {"tissue": tissue, "reference": reference_label, "examples": [],
                "placeholder": "e.g. What are the best marker genes of this cell state?"}

    sizes = _chat_states_by_size(cache)
    states = [state for state, _ in sizes if state and state.lower() not in ("nan", "none")]
    first = states[0] if states else ""
    second = states[1] if len(states) > 1 else ""
    markers = _chat_markers_for_states(app, meta, cache, [first] if first else [], 6)
    genes = [str(row["gene"]) for row in markers][:2]

    examples: List[str] = []
    if first:
        examples.append(f"What are the best marker genes of {first} cells?")
    if first and second:
        examples.append(f"What distinguishes {first} from {second} cells?")
    if genes:
        examples.append(f"Where is {genes[0]} expressed?")
    if len(genes) > 1:
        examples.append(f"Which cell states express {genes[1]}?")
    contrast = _chat_contrast(meta)
    if contrast and first:
        examples.append(f"Which genes are significant in {contrast['case']} versus "
                        f"{contrast['control']} in {first} cells?")
        examples.append(f"Which cell type is most affected in {contrast['case']} versus "
                        f"{contrast['control']}?")
    return {"tissue": tissue, "reference": reference_label,
            "cluster_key": str(cache["cluster_key"]),
            "n_states": len(states), "has_contrast": bool(contrast),
            "examples": examples,
            "placeholder": (f"e.g. What are the best marker genes of {first} cells?"
                            if first else "e.g. What are the best marker genes of this cell state?")}


def _chat_read_question(question: str, cache: Dict[str, Any], meta: Dict) -> Dict[str, Any]:
    """Ask the assistant which supported query this sentence means.

    The model sees the question and the names in this dataset. It never sees the
    data, so it cannot invent a number: it chooses the reading, and the executors
    below compute the answer from the job's own files.
    """
    contrast = _chat_contrast(meta)
    payload = json.dumps({
        "question": question,
        "states": [state for state, _ in _chat_states_by_size(cache)],
        "contrasts": [contrast["label"]] if contrast else [],
        "covariates": [entry["field"] for entry in _groupable_columns(cache)],
        "modalities": [str(entry.get("id")) for entry
                       in ((meta.get("modalities") or {}).get("available") or [])] or ["rna"],
    }).encode()

    last_error = None
    for attempt in (1, 2):
        try:
            request = urllib.request.Request(
                CHAT_ASSISTANT_URL, data=payload,
                headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.loads(response.read().decode())
        except Exception as exc:  # noqa: BLE001 - the site may be down
            last_error = exc
            if attempt == 1:
                time.sleep(0.25)
    raise HTTPException(
        status_code=503,
        detail=(f"the assistant at {CHAT_ASSISTANT_URL} did not answer ({last_error}). "
                "cellHarmony web holds no model of its own."))


'''
anchor = "def _build_umap_payload("
assert src.count(anchor) == 1
src = src.replace(anchor, CHAT + anchor)

open(P, "w").write(src)
print("patch6 (chat helpers) applied")
