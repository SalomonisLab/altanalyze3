"""Adapt uploaded job outputs to the shared regulatory analysis data interface.

No imputation or differential analysis runs while serving these views. Statistics
come from retained, completed comparisons with matching groups and cell-state fields.
"""
from functools import cached_property
from types import SimpleNamespace
from pathlib import Path
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


def comparison_metadata(meta, contrast):
    """A completed comparison owns its population field, independent of Explore."""
    run = completed_differentials(meta).get(contrast)
    if run is None:
        raise ValueError("Unknown completed comparison.")
    column = (run.get("config") or {}).get("population_col")
    if not column:
        raise ValueError("The completed comparison has no recorded cell-state field.")
    from .state_layers import apply_layer
    return dict(apply_layer(meta, column), cluster_key=column,
                active_cell_state_layer=column, differential=run)


def aggregate_state_column(meta, record, adata, key):
    """Resolve the recorded supervised rename without regrouping native aggregates."""
    if key in adata.obs:return key
    if key != 'supervised_state':raise KeyError(f"GRN output has no {key} cell-state column")
    alias=(record.get('state_columns') or {}).get(key)
    source=Path(record.get('h5ad') or record.get('differential_h5ad') or '')
    # Older Both jobs retained the native branch record rather than an alias map.
    branch=source.parent.parent/'job.json'
    if not alias and branch.is_file():
        native=json.loads(branch.read_text())
        if native.get('job_id')==meta['job_id'] and branch.parent.parent.name=='supervised':
            alias=native.get('cluster_key')
    if not alias or alias not in adata.obs or not source.is_file():
        raise KeyError(f"GRN output has no recorded {key} cell-state column")
    import h5py
    from anndata.io import read_elem
    with h5py.File(source,'r') as f:obs=read_elem(f['obs'])
    if key not in obs or alias not in obs or obs[[key,alias]].isna().any().any() or not np.array_equal(obs[key].astype(str).to_numpy(),obs[alias].astype(str).to_numpy()):
        raise ValueError('The recorded supervised GRN annotation alias does not match every source cell; no aggregate labels were substituted.')
    return alias


class UploadedGrnData:
    def __init__(self, app, meta, contrast=""):
        from importlib import import_module
        web = import_module("altanalyze3.components.cellHarmony.webapp.app")
        if contrast:
            meta = comparison_metadata(meta, contrast)
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
            key = aggregate_state_column(self.meta,self.modality_manifest()[modality],adata,key)
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
