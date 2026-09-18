"""Adapt uploaded job outputs to the shared regulatory analysis data interface.

No imputation or differential analysis runs while serving these views. Statistics
come from retained, completed comparisons with matching groups and cell-state fields.
"""
from functools import cached_property
from types import SimpleNamespace
import json

import numpy as np

from altanalyze3.components.cellHarmony.modalities import modality_artifacts


def completed_differentials(meta):
    runs = dict(meta.get("differential_history") or {})
    current = meta.get("differential") or {}
    if current.get("status") == "completed" and current.get("run_id"):
        runs[current["run_id"]] = current
    return {key: value for key, value in runs.items() if value.get("status") == "completed"}


def comparison_entry(run_id, run):
    config = run.get("config") or {}
    identity = {key: config.get(key) for key in
                ("population_col", "sample_field", "comparison_type")}
    for key in ("group1_samples", "group2_samples"):
        identity[key] = sorted(config.get(key) or [])
    comparison = f"{run.get('case_label', 'Group 1')} vs {run.get('control_label', 'Group 2')}"
    return {"id": run_id, "comparison": comparison, "kind": "per_cell_state",
            "contrast": json.dumps(identity, sort_keys=True),
            "modality": config.get("modality", "rna")}


class UploadedGrnData:
    def __init__(self, app, meta):
        from importlib import import_module
        web = import_module("altanalyze3.components.cellHarmony.webapp.app")
        self.web, self.app, self.meta = web, app, meta
        self.id = meta["job_id"]
        self.runs = completed_differentials(meta)
        self.current_contrast = (meta.get("differential") or {}).get("run_id", "")
        self._stores = {}
        self.arm_scope = "all uploaded cells in the selected cell state; not restricted to comparison arms"

    @cached_property
    def rna(self):
        return self.web._get_expression_cache(self.app, self.meta)

    @cached_property
    def states(self):
        return list(dict.fromkeys(map(str, self.rna["populations"])))

    def _mean(self, adata, labels):
        return np.stack([
            np.asarray(adata.X[np.asarray(labels, dtype=str) == state].mean(axis=0)).ravel()
            if np.any(np.asarray(labels, dtype=str) == state)
            else np.full(adata.n_vars, np.nan)
            for state in self.states], axis=1)

    @cached_property
    def stats_mean(self):
        return self._mean(self.rna["adata"], self.rna["populations"])

    def resolve_gene(self, gene):
        names = list(map(str, self.rna["var_names"]))
        return names.index(gene) if gene in names else None

    def modality_manifest(self):
        return modality_artifacts(self.meta)

    def modality(self, modality):
        if modality in self._stores:
            return self._stores[modality]
        if modality not in self.modality_manifest():
            raise KeyError(modality)
        if modality == "grn":
            adata = self.web._grn_edges_adata(self.meta)
            key = self.rna["cluster_key"]
            if key not in adata.obs:
                raise KeyError(f"GRN output has no {key} cell-state column")
            labels = adata.obs[key].astype(str).to_numpy()
            kind = "per_sample_cell_state"
        else:
            cache = self.web._get_expression_cache(self.app, self.meta, modality=modality)
            adata, labels, kind = cache["adata"], cache["populations"], "per_cell"
        statistic = str(adata.uns.get("activity_statistic") or
                        ("legacy standardized target-set enrichment" if
                         self.modality_manifest()[modality].get("legacy_enrichment") else
                         "sum of predicted outgoing edge activity"))
        store = SimpleNamespace(features=list(adata.var_names.astype(str)),
                                stats_mean=self._mean(adata, labels), kind=kind,
                                label=self.web._modality_definition(self.meta, modality)["label"],
                                statistic=statistic + ", averaged over uploaded cells in the state")
        self._stores[modality] = store
        return store

    def deg_manifest(self):
        # Latest matching completed run wins when the user repeats a comparison.
        return {"comparisons": [comparison_entry(key, value)
                                for key, value in reversed(list(self.runs.items()))]}

    def deg_table(self, comp_id, *, max_rows=200000, fdr_max=None, state=None):
        run = self.runs[comp_id]
        snapshot = dict(self.meta, differential=run)
        frame = self.web._get_differential_detail_table(self.app, snapshot).copy()
        if state:
            frame = frame.loc[frame["population"].astype(str) == state]
        if fdr_max is not None:
            frame = frame.loc[frame["fdr"] <= fdr_max]
        # JSON round-trip converts pandas/NumPy values and non-finite numbers.
        return {"rows": json.loads(frame.head(max_rows).to_json(orient="records"))}
