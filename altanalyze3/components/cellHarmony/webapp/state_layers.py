"""Request-local cell-state selection shared by every Explore and Chat reader."""
from contextvars import ContextVar
from copy import deepcopy

from ..flask.job_manager import JobStore

ACTIVE_LAYER = ContextVar("scalable_analysis_layer", default="")


def apply_layer(meta, key):
    layers = meta.get("cell_state_layers") or {}
    default = layers.get("default") or meta.get("cluster_key")
    entry = next((e for e in layers.get("layers", []) if e["key"] == key), None)
    if not layers or not entry or key == default:
        return dict(meta, active_cell_state_layer=default) if layers else meta
    result = deepcopy(meta)
    result.update(cluster_key=key, active_cell_state_layer=key)
    for field in ("marker_analysis", "fastcomm_analysis"):
        result[field] = entry.get(field) or {}
    result["marker_analysis_by_modality"] = entry.get("marker_analysis_by_modality") or {"rna": result["marker_analysis"]}
    # This default affects the form only; posted comparisons carry an explicit
    # population column and saved differential outputs retain their own column.
    result.setdefault("differential_options", {})["default_population_col"] = key
    return result


class AnalysisJobStore(JobStore):
    def get_job(self, job_id):
        return apply_layer(super().get_job(job_id), ACTIVE_LAYER.get())
