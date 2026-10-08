"""The discover BioMarkers reader, shared by standalone and unified Explore."""
from pathlib import Path
from typing import Dict, List
import numpy as np
import pandas as pd
from altanalyze3.components.goelite.structures import compute_z_score

GOELITE_MAX_FDR = 0.05
GOELITE_MIN_Z = 2.0

def _cluster_for_state(meta: Dict, state: str) -> str:
    """A cell state of either layer -> its ICGS3 cluster id."""
    names = (meta.get("cell_state_layers") or {}).get("names") or {}
    reverse = {str(v): str(k) for k, v in names.items()}
    return reverse.get(str(state), str(state))


def _state_label(meta: Dict, cluster: str) -> str:
    """An ICGS3 cluster id as the active layer names it."""
    layers = meta.get("cell_state_layers") or {}
    if meta.get("active_cell_state_layer") in {"cluster", "unsupervised_cluster"}:
        return cluster
    return str((layers.get("names") or {}).get(cluster, cluster))


def goelite_states(meta: Dict) -> List[str]:
    """The active layer's cell states that hold at least one BioMarkers term."""
    goelite = (meta.get("icgs3_analysis") or {}).get("goelite") or {}
    path = Path(str(goelite.get("enrichment_tsv") or ""))
    if not path.is_file():
        return []
    clusters = set(pd.read_csv(path, sep="\t", usecols=["cluster"])["cluster"].astype(str))
    order = (meta.get("icgs3_analysis") or {}).get("clusters") or sorted(clusters)
    return [_state_label(meta, c) for c in order if c in clusters]


def build_goelite_payload(meta: Dict, state: str) -> Dict:
    """One cell state's BioMarkers terms, in the shape scALABLE's GO Terms view draws.

    p, FDR and overlap are ICGS3's own (ICGS.biomarker_enrichment: hypergeometric p, BH FDR
    over every cluster and term). ICGS3 stores no z-score, so z is GO-Elite's own
    `compute_z_score` (goelite/structures.py) on ICGS3's overlap, query size, term size and
    the gene universe ICGS3 tested against.
    """
    goelite = (meta.get("icgs3_analysis") or {}).get("goelite") or {}
    path = Path(str(goelite.get("enrichment_tsv") or ""))
    background = int(goelite.get("background_size") or 0)
    cluster = _cluster_for_state(meta, state)
    label = _state_label(meta, cluster)
    payload = {"population": label, "cluster": cluster, "terms": [], "labels": [],
               "statistics": {"p_value": "ICGS3 hypergeometric p", "fdr": "ICGS3 BH FDR, all clusters x terms",
                              "z_score": "GO-Elite compute_z_score on ICGS3 overlap counts",
                              "background_size": background,
                              "highlight": f"FDR <= {GOELITE_MAX_FDR} and z > {GOELITE_MIN_Z}"}}
    if not path.is_file() or background <= 0:
        payload["message"] = "This job has no GO-Elite BioMarkers table."
        return payload
    frame = pd.read_csv(path, sep="\t")
    frame = frame.loc[frame["cluster"].astype(str) == cluster]
    if frame.empty:
        payload["message"] = f"No BioMarkers term overlaps the markers of {label}."
        return payload
    # ICGS3 labels each cluster with its first row after sorting by FDR, p and name.
    frame = frame.sort_values(["fdr", "p_value", "term_name"]).reset_index(drop=True)
    terms = []
    for index, row in frame.iterrows():
        z_score = float(compute_z_score(int(row["overlap"]), int(row["query_size"]),
                                        int(row["term_size"]), background))
        fdr = float(row["fdr"])
        genes = [g for g in str(row.get("overlap_genes") or "").split(",") if g]
        significant = fdr <= GOELITE_MAX_FDR and z_score > GOELITE_MIN_Z
        terms.append({
            "term_id": "", "term_name": str(row["term_name"]), "direction": "up",
            "fdr": fdr, "p_value": float(row["p_value"]), "z_score": z_score,
            "fdr_plot": float(min(max(fdr, 1e-300), 1.0)), "score": float(-np.log10(min(max(fdr, 1e-300), 1.0))),
            "overlap": int(row["overlap"]), "query_size": int(row["query_size"]), "term_size": int(row["term_size"]),
            "selected": significant, "is_positive_sig": significant, "is_selected_positive_sig": significant,
            "is_prediction": index == 0, "overlap_genes": genes, "selected_gene": genes[0] if genes else None,
        })
    payload["terms"] = terms
    labelled = [t for t in terms if t["is_selected_positive_sig"]][:4]
    if terms[0] not in labelled:
        labelled = [terms[0]] + labelled[:3]
    payload["labels"] = [{"term_name": t["term_name"], "z_score": t["z_score"], "fdr_plot": t["fdr_plot"],
                          "selected_gene": t["selected_gene"], "overlap_genes": t["overlap_genes"],
                          "label_color": "#1f19c7" if t["is_prediction"] else "#111827",
                          "label_rank": i, "label_role": "prediction" if t["is_prediction"] else "top"}
                         for i, t in enumerate(labelled)]
    return payload


