"""Adapter: a precomputed scalable_viewer bundle -> the cellHarmony webapp `meta` contract.

The cellHarmony webapp (scALABLE) is the reference implementation. Its payload builders
and PDF renderers all take `(app, meta, ...)` where

  * `app`  is a FastAPI app carrying `app.state.<cache>` dictionaries, and
  * `meta` is the job-metadata dict a finished job writes.

Neither reads a job object. So a precomputed bundle can drive the whole scALABLE app
unchanged, provided two things are supplied:

  1. a `meta` dict with the keys the builders read, and
  2. pre-seeded caches, so the two functions that would open an h5ad
     (`_get_expression_cache`, `_get_marker_heatmap_cache_entry`) return the bundle's
     memory-mapped arrays instead.

Nothing in this module draws a plot, computes a statistic or duplicates a builder.
Every figure, colour ramp and payload comes from
/Users/saljh8/Documents/GitHub/altanalyze3/altanalyze3/components/cellHarmony/webapp/app.py

Reference: app.py:1353 `_get_expression_cache` (cache contract), app.py:1260
`_get_marker_heatmap_cache_entry` (marker matrix contract), app.py:1524
`_get_differential_cache_entry` (differential table contract).
"""
from __future__ import annotations

import copy
import json
import os
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import data_api as da

# ---------------------------------------------------------------------------------
# The webapp module. Imported, never copied.
# ---------------------------------------------------------------------------------
from importlib import import_module  # noqa: E402

# `from ...webapp import app` would bind the FastAPI instance the package __init__
# re-exports, not the module. Import the module explicitly.
W = import_module("altanalyze3.components.cellHarmony.webapp.app")


# =================================================================================
# 1. AnnData surface over the memory-mapped bundle
# =================================================================================

class _ColumnView:
    """What `adata[rows, gene]` returns. Only `.X` is read by app.py."""

    __slots__ = ("X",)

    def __init__(self, values: np.ndarray):
        self.X = values


class BundleAnnData:
    """The minimal AnnData surface app.py uses on the expression cache entry.

    app.py touches exactly four things on `cache_entry["adata"]`:
      `adata[:, gene].X`            app.py:3531 and app.py:3479
      `adata[mask, gene].X`         app.py:3188
      `adata.obs[col]`              app.py:3179-3180
      `adata.obs.columns`           app.py:3162
    A dense gene column is one slice of the gene-major memmap: no h5ad is opened and
    no matrix is held in memory.
    """

    def __init__(self, ds: da.Dataset, obs: pd.DataFrame, var_names: pd.Index):
        self._ds = ds
        self.obs = obs
        self.var_names = var_names
        self.n_obs = int(ds.n_cells)
        self.n_vars = int(len(var_names))
        self.isbacked = False

    def _dense_column(self, gene: str) -> np.ndarray:
        row = self._ds.resolve_gene(gene)
        if row is None:
            raise KeyError(gene)
        idx, val = self._ds.gene_column(row)
        out = np.zeros(self.n_obs, dtype=np.float32)
        if idx.size:
            out[idx.astype(np.int64)] = val
        return out

    def __getitem__(self, key):
        if not (isinstance(key, tuple) and len(key) == 2):
            raise TypeError(f"BundleAnnData only supports [rows, gene]; got {key!r}")
        rows, gene = key
        if not isinstance(gene, str):
            raise TypeError(f"BundleAnnData only supports a single gene name; got {gene!r}")
        values = self._dense_column(gene)
        if isinstance(rows, slice):
            if rows != slice(None):
                values = values[rows]
        else:
            values = values[np.asarray(rows)]
        return _ColumnView(values)


# =================================================================================
# 2. meta contract
# =================================================================================

_EXCLUDED_FILTER_SUFFIXES = ("__n_obs", "__homogeneous")


def _categorical_covariates(ds: da.Dataset) -> Dict[str, List[str]]:
    """Categorical covariates offered as `Filter data to display` annotations.

    Excludes the bookkeeping columns precompute.py carries (`*__n_obs`,
    `*__homogeneous`) and single-level columns, which cannot filter anything.
    """
    out: Dict[str, List[str]] = {}
    for name, info in ds.covariate_names().items():
        if info.get("kind") != "categorical":
            continue
        if name.endswith(_EXCLUDED_FILTER_SUFFIXES):
            continue
        cats = [str(c) for c in info.get("categories", [])]
        if len(cats) < 2 or len(cats) > 120:
            continue
        out[name] = cats
    return out


def _pick_sample_field(ds: da.Dataset, categorical: Dict[str, List[str]]) -> str:
    """The field `Cell frequency` normalizes by. Same preference order the webapp uses
    (app.py:1441 `("Library", "group", "sample")`), then the bundle's own sample-like
    covariates."""
    for candidate in ("Library", "group", "sample", "Donor", "pool", "sample_id", "donor"):
        if candidate in categorical:
            return candidate
    return next(iter(categorical), "")


def build_obs(ds: da.Dataset) -> Tuple[pd.DataFrame, Dict[str, np.ndarray], str]:
    """The obs frame, the display-filter value arrays and the sample field.

    Values come from `<prefix>_cells.npz`, written by precompute.py. Nothing is
    recomputed here.
    """
    cells = ds.cells
    barcodes = _obs_names(ds)
    states = np.asarray(ds.states, dtype=object)
    codes = np.asarray(ds.state_code)
    population = states[codes.astype(np.int64)].astype(str)

    obs_filter_values: Dict[str, np.ndarray] = {ds.cluster_key: population.astype(str)}
    columns: Dict[str, Any] = {ds.cluster_key: pd.Categorical(population, categories=ds.states)}

    categorical = _categorical_covariates(ds)
    for name, cats in categorical.items():
        if name == ds.cluster_key:
            continue
        arr = np.asarray(cells["cov_cat_" + name])
        cat_arr = np.asarray(cats, dtype=object)
        text = np.where(arr >= 0, cat_arr[np.clip(arr, 0, len(cats) - 1)], "").astype(str)
        obs_filter_values[name] = text
        columns[name] = pd.Categorical(text, categories=[""] + cats)

    obs = pd.DataFrame(columns, index=pd.Index(barcodes, name="barcode"))
    sample_field = _pick_sample_field(ds, categorical)
    return obs, obs_filter_values, sample_field


def _obs_names(ds: da.Dataset) -> np.ndarray:
    """Cell barcodes, in bundle row order, read from `<prefix>_clusters.tsv`."""
    cached = getattr(ds, "_viewer_obs_names", None)
    if cached is not None:
        return cached
    names: List[str] = []
    with open(ds.paths.clusters, "r") as fh:
        fh.readline()
        for line in fh:
            names.append(line.split("\t", 1)[0])
    arr = np.asarray(names, dtype=str)
    if arr.size != ds.n_cells:
        raise ValueError(
            f"{ds.paths.clusters} holds {arr.size} barcodes but the bundle declares "
            f"{ds.n_cells} cells"
        )
    ds._viewer_obs_names = arr  # type: ignore[attr-defined]
    return arr


def _umap_xy(ds: da.Dataset) -> Tuple[np.ndarray, np.ndarray]:
    emb = np.asarray(ds.embedding, dtype=float)
    return emb[:, 0].copy(), emb[:, 1].copy()


def build_meta(
    ds: da.Dataset,
    *,
    state_dir: Path,
    marker_source_tsv: Optional[str] = None,
    marker_heatmap_cache: Optional[str] = None,
    fastcomm_analysis: Optional[Dict[str, Any]] = None,
    marker_networks: Optional[List[Dict[str, str]]] = None,
    differential_assets: Optional[Dict[str, Dict[str, Any]]] = None,
    contrast_id: Optional[str] = None,
) -> Dict[str, Any]:
    """The job-metadata dict the scALABLE builders read, derived from one bundle."""
    obs, obs_filter_values, sample_field = build_obs(ds)
    categorical = _categorical_covariates(ds)
    sample_values = sorted({v for v in obs_filter_values.get(sample_field, []) if v}) if sample_field else []

    combined_h5ad = str(ds.meta.get("source_h5ad") or "")
    if combined_h5ad and not os.path.isfile(combined_h5ad):
        combined_h5ad = ""

    marker_analysis: Dict[str, Any] = {
        "enabled": bool(ds.sv.get("markers", {}).get("available")),
        "status": "completed",
        "markers_tsv": marker_source_tsv or ds.sv.get("markers", {}).get("source") or ds.paths.markers,
        "populations": list(ds.sv.get("markers", {}).get("clusters") or []),
    }
    if marker_heatmap_cache and os.path.isfile(marker_heatmap_cache):
        marker_analysis["heatmap_cache"] = marker_heatmap_cache
    if marker_networks:
        marker_analysis["networks"] = marker_networks

    deg_comparisons = ds.deg_manifest().get("comparisons", [])
    per_state = [c for c in deg_comparisons if c.get("kind") == "per_cell_state"]
    chosen = None
    if contrast_id:
        chosen = next((c for c in per_state if c["id"] == contrast_id), None)
    if chosen is None and per_state:
        chosen = per_state[0]

    meta: Dict[str, Any] = {
        "job_id": ds.id,
        "status": "completed",
        "progress": 100,
        "message": f"Precomputed bundle {ds.paths.prefix} loaded from {ds.paths.bundle_dir}.",
        "created_at": ds.sv.get("built_utc"),
        "updated_at": ds.sv.get("built_utc"),
        "species": "Hs",
        "reference": ds.label,
        "cluster_key": ds.cluster_key,
        "reference_cluster_key": ds.cluster_key,
        "reference_coords_tsv": ds.paths.umap,
        "reference_clusters_tsv": ds.paths.clusters,
        "default_gene": _default_gene(ds),
        "files": [{"sample_name": value, "filename": f"{value}.h5ad"} for value in (sample_values or [ds.id])],
        "artifacts": {
            "combined_h5ad": combined_h5ad,
            "umap_coordinates": ds.paths.umap,
            "assignments": ds.paths.clusters,
        },
        "modalities": {"default": "rna", "available": [dict(W._DEFAULT_MODALITY_DEFINITIONS["rna"])]},
        "marker_analysis": marker_analysis,
        "marker_analysis_by_modality": {"rna": marker_analysis},
        "fastcomm_analysis": fastcomm_analysis or {"enabled": False},
        "scalable_viewer": {
            "bundle_dir": ds.paths.bundle_dir,
            "prefix": ds.paths.prefix,
            "state_dir": str(state_dir),
            "deg_comparisons": [
                {"id": c["id"], "comparison": c["comparison"], "n_rows": c.get("n_rows")}
                for c in per_state
            ],
        },
    }
    meta["artifacts"] = {k: v for k, v in meta["artifacts"].items() if v}

    meta["differential_options"] = {
        "enabled": True,
        "sample_names": sample_values,
        "sample_fields": [{"value": f, "label": f} for f in categorical],
        "sample_values": {f: cats for f, cats in categorical.items()},
        "default_sample_field": _contrast_sample_field(ds, chosen, categorical) or sample_field,
        "population_columns": [{"value": ds.cluster_key, "label": ds.cluster_key, "n_categories": len(ds.states)}],
        "default_population_col": ds.cluster_key,
        "comparison_types": ["pseudobulk", "cells"],
        "upload_profile": {
            "total_files": max(len(sample_values), 2),
            "extensions": ["h5ad"] * max(len(sample_values), 2),
            "single_h5": False,
            "single_h5ad": False,
            "multiple_h5": False,
            "differential_eligible": True,
            "allow_alternate_population_fields": False,
            "show_secondary_display_filter": True,
        },
    }
    meta["scalable_viewer"]["differential_assets"] = differential_assets or {}
    meta["differential"] = build_differential_block(
        ds, chosen, categorical, (differential_assets or {}).get(chosen["id"]) if chosen else None
    )
    return meta


def _default_gene(ds: da.Dataset) -> str:
    rows = ds.markers()
    if rows:
        return str(rows[0]["gene"])
    return ds.symbols[0] if ds.n_genes else ""


def _norm_label(value: str) -> str:
    return "".join(ch for ch in str(value).lower() if ch.isalnum())


def _contrast_group_field(ds, comparison, categorical: Dict[str, List[str]]
                          ) -> Tuple[str, str, str]:
    """(covariate field, case value, control value) for one precomputed contrast.

    The DEG table names its groups in `case_label` / `control_label`. Those strings do
    not always spell the covariate value exactly - `cancer_vs_no_cancer` writes
    `no_cancer` where obs holds `no cancer` - so matching ignores case, spaces,
    underscores and hyphens, and the covariate's own spelling is returned.
    """
    if not comparison:
        return "", "", ""
    case, control = _contrast_labels(comparison)
    case_key, control_key = _norm_label(case), _norm_label(control)
    for field, cats in categorical.items():
        keys = {_norm_label(c): c for c in cats}
        if case_key in keys and control_key in keys and case_key != control_key:
            return field, keys[case_key], keys[control_key]
    return "", "", ""


def _contrast_sample_field(ds, comparison, categorical: Dict[str, List[str]]) -> str:
    return _contrast_group_field(ds, comparison, categorical)[0]


def _contrast_labels(comparison: Dict[str, Any]) -> Tuple[str, str]:
    """Group labels for a contrast, read from the DEG table's own case/control columns."""
    path = str(comparison.get("path") or "")
    if path and os.path.isfile(path):
        try:
            head = pd.read_csv(path, sep="\t", nrows=1)
            case = str(head.get("case_label", pd.Series([""])).iloc[0]).strip()
            control = str(head.get("control_label", pd.Series([""])).iloc[0]).strip()
            if case and control and case != "nan" and control != "nan":
                return case, control
        except (OSError, ValueError, IndexError):
            pass
    name = str(comparison.get("comparison") or "")
    if "_vs_" in name:
        case, control = name.split("_vs_", 1)
        return case, control
    return name or "Group 1", "Group 2"


def build_differential_block(ds, comparison, categorical: Dict[str, List[str]],
                             assets: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """meta['differential'] for one precomputed contrast, in the shape app.py:1524 reads."""
    if not comparison:
        return {"status": "unavailable", "progress": 0, "artifacts": {}, "networks": [],
                "message": "This bundle carries no per-cell-state differential table."}
    case, control = _contrast_labels(comparison)
    field, case_value, control_value = _contrast_group_field(ds, comparison, categorical)
    assets = assets or {}
    artifacts: Dict[str, str] = {f"DEG_detailed_{comparison['comparison']}": comparison["path"]}
    go_tsv = str(assets.get("goelite_tsv") or "")
    if go_tsv and os.path.isfile(go_tsv):
        artifacts["goelite_tsv"] = go_tsv
    networks = [n for n in (assets.get("networks") or []) if os.path.isfile(str(n.get("tsv", "")))]
    return {
        "status": "completed",
        "progress": 100,
        "message": f"Precomputed differential: {comparison['comparison']} ({comparison.get('n_rows')} rows).",
        "run_id": comparison["id"],
        "case_label": case,
        "control_label": control,
        "go_terms_included": bool(go_tsv),
        "config": {
            "modality": "rna",
            "population_col": ds.cluster_key,
            "sample_field": field,
            "group1_samples": [case_value] if field else [],
            "group2_samples": [control_value] if field else [],
            "comparison_type": "pseudobulk",
        },
        "artifacts": artifacts,
        "networks": networks,
        "go_terms_available": bool(go_tsv),
    }


# =================================================================================
# 3. Cache seeding - the step that keeps every h5ad closed
# =================================================================================

def seed_expression_cache(app, ds: da.Dataset, meta: Dict[str, Any]) -> Dict[str, Any]:
    """Pre-fill `app.state.expression_cache` so app.py:1353 returns on its cache hit.

    The signature app.py:1370-1376 compares is (h5ad_path, umap_path, cluster_key,
    modality); this builds exactly that key, so `_get_expression_cache` never reaches
    its `ad.read_h5ad` at app.py:1391.
    """
    obs, obs_filter_values, sample_field = build_obs(ds)
    var_names = np.asarray(ds.symbols, dtype=str)
    adata = BundleAnnData(ds, obs, pd.Index(var_names))
    umap_x, umap_y = _umap_xy(ds)
    obs_names = _obs_names(ds)

    fields = [{"value": f, "label": f} for f in obs_filter_values]
    values = {f: sorted({v for v in vals if v}) for f, vals in obs_filter_values.items()}
    default_primary = sample_field if sample_field in obs_filter_values else (fields[0]["value"] if fields else "")
    default_secondary = ds.cluster_key if ds.cluster_key != default_primary else ""

    entry = {
        "h5ad_path": str(W._modality_h5ad_path(meta, "rna")),
        "umap_path": str(meta.get("artifacts", {}).get("umap_coordinates") or ""),
        "cluster_key": str(ds.cluster_key),
        "modality": "rna",
        "adata": adata,
        "obs_names": obs_names,
        "var_names": var_names,
        "populations": obs[ds.cluster_key].astype(str).to_numpy(),
        "sample_field": sample_field,
        "sample_labels": obs_filter_values.get(sample_field) if sample_field else None,
        "umap_x": umap_x,
        "umap_y": umap_y,
        "obs_filter_values": obs_filter_values,
        "display_filters_meta": {
            "fields": fields,
            "values": values,
            "default_primary_field": default_primary,
            "default_secondary_field": default_secondary,
            "show_secondary": True,
        },
    }
    app.state.expression_cache[f"{ds.id}:rna"] = entry
    return entry


def seed_marker_heatmap_cache(app, ds: da.Dataset, meta: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Pre-fill `app.state.marker_heatmap_cache` from the bundle's marker matrix npz.

    app.py:1284 returns the cached entry when the signature matches, so the npz is
    read once here and never again.
    """
    marker_analysis = meta.get("marker_analysis") or {}
    cache_path = str(marker_analysis.get("heatmap_cache") or "")
    if not cache_path or not os.path.isfile(cache_path):
        return None
    with np.load(cache_path, allow_pickle=False) as npz:
        matrix = np.asarray(npz["matrix"], dtype=np.float32)
        row_ids = np.asarray(npz["row_ids"], dtype=str)
        col_ids = np.asarray(npz["col_ids"], dtype=str)
        col_barcodes = np.asarray(npz["col_barcodes"], dtype=str)
    entry = {
        "signature": {
            "modality": "rna",
            "heatmap_cache": cache_path,
            "heatmap_tsv": str(marker_analysis.get("heatmap_tsv") or ""),
            "expression_tsv": str(marker_analysis.get("expression_tsv") or ""),
        },
        "modality": "rna",
        "source": "cache",
        "source_path": Path(cache_path),
        "matrix": matrix,
        "row_ids": row_ids,
        "col_ids": col_ids,
        "col_barcodes": col_barcodes,
        "tsv_path": None,
    }
    app.state.marker_heatmap_cache[f"{ds.id}:rna"] = entry
    return entry


# =================================================================================
# 4. A JobStore over bundles
# =================================================================================

class BundleJobStore:
    """The read side of `JobStore` (flask/job_manager.py:16) backed by bundles.

    Every scALABLE endpoint begins with `_job_resources(app)` (app.py:684) then
    `store.job_exists(job_id)` / `store.get_job(job_id)`. Satisfying those two calls
    is enough for the whole app to serve a bundle. Job creation, QC and alignment are
    not supported: a bundle is already aligned.
    """

    def __init__(self, catalog: da.Catalog, state_dir: Path, marker_assets: Dict[str, Dict[str, str]]):
        self._catalog = catalog
        self._state_dir = Path(state_dir)
        self._marker_assets = marker_assets
        self._meta: Dict[str, Dict[str, Any]] = {}
        self._lock = threading.Lock()

    # ---- construction -------------------------------------------------------

    def dataset(self, job_id: str) -> da.Dataset:
        return self._catalog.get(job_id)

    def ensure(self, app, job_id: str) -> Dict[str, Any]:
        """Build the meta for one bundle and seed its caches. Idempotent."""
        with self._lock:
            if job_id in self._meta:
                return self._meta[job_id]
            ds = self._catalog.get(job_id)
            assets = self._marker_assets.get(job_id, {})
            meta = build_meta(
                ds,
                state_dir=self._state_dir / job_id,
                # `marker_analysis["markers_tsv"]` has exactly one consumer in scALABLE:
                # `_build_marker_network_payload` (app.py:2217), which colours network nodes
                # by fold change. Prefer the merged lookup, which covers the genes of both
                # the unique and the redundant marker table (prepare_assets.py).
                marker_source_tsv=assets.get("marker_fold_lookup_tsv") or assets.get("markers_tsv"),
                marker_heatmap_cache=assets.get("heatmap_cache"),
                fastcomm_analysis=assets.get("fastcomm_analysis"),
                marker_networks=assets.get("networks"),
                differential_assets=assets.get("differential"),
            )
            self._meta[job_id] = meta
            seed_expression_cache(app, ds, meta)
            seed_marker_heatmap_cache(app, ds, meta)
            return meta

    # ---- JobStore surface ---------------------------------------------------

    def job_exists(self, job_id: str) -> bool:
        return str(job_id) in self._catalog.ids()

    def get_job(self, job_id: str) -> Dict[str, Any]:
        meta = self._meta.get(str(job_id))
        if meta is None:
            raise KeyError(job_id)
        return copy.deepcopy(meta)

    def update_job(self, job_id: str, **changes) -> Dict[str, Any]:
        with self._lock:
            meta = self._meta[str(job_id)]
            meta.update(changes)
            return copy.deepcopy(meta)

    def add_artifact(self, job_id: str, key: str, path: Path) -> None:
        with self._lock:
            self._meta[str(job_id)].setdefault("artifacts", {})[key] = str(path)

    def append_log(self, job_id: str, message: str) -> None:
        path = self.logs_dir(job_id) / "pipeline.log"
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {message}\n")

    def _dir(self, job_id: str, name: str) -> Path:
        path = self._state_dir / str(job_id) / name
        path.mkdir(parents=True, exist_ok=True)
        return path

    def uploads_dir(self, job_id: str) -> Path:
        return self._dir(job_id, "uploads")

    def outputs_dir(self, job_id: str) -> Path:
        return self._dir(job_id, "outputs")

    def logs_dir(self, job_id: str) -> Path:
        return self._dir(job_id, "logs")

    def list_jobs(self) -> List[Dict[str, Any]]:
        return [copy.deepcopy(m) for m in self._meta.values()]

    def purge_old_jobs(self, *, max_age_hours: float = 8.0) -> int:
        return 0


# =================================================================================
# 5. Marker heatmap matrix: GCT -> the npz scALABLE reads
# =================================================================================

def gct_to_marker_cache(gct_path: str, out_npz: str) -> Dict[str, Any]:
    """Convert a Morpheus GCT (genes x cells) to the npz app.py:1294 reads.

    The GCT is an existing project artifact; this only reshapes it into the four
    arrays scALABLE's marker heatmap cache holds. No value is changed.
    """
    with open(gct_path, "r") as fh:
        fh.readline()
        dims = fh.readline().split()
        n_row_annot = int(dims[2])
        n_col_annot = int(dims[3])
        header = fh.readline().rstrip("\n").split("\t")
        col_ids = header[1 + n_row_annot:]
        for _ in range(n_col_annot):        # GCT 1.3 column-annotation rows
            fh.readline()
        rows: List[str] = []
        values: List[List[float]] = []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 1 + n_row_annot:
                continue
            rows.append(parts[0])
            values.append([float(v) if v not in ("", "NA", "nan") else 0.0
                           for v in parts[1 + n_row_annot:]])
    matrix = np.asarray(values, dtype=np.float32)
    col_barcodes = np.asarray([c.split(":", 1)[1] if ":" in c else c for c in col_ids], dtype=str)
    os.makedirs(os.path.dirname(out_npz), exist_ok=True)
    np.savez_compressed(
        out_npz,
        matrix=matrix,
        row_ids=np.asarray(rows, dtype=str),
        col_ids=np.asarray(col_ids, dtype=str),
        col_barcodes=col_barcodes,
    )
    return {"path": out_npz, "n_rows": matrix.shape[0], "n_cols": matrix.shape[1]}


# =================================================================================
# 6. fastComm: read a completed run of the real backend
# =================================================================================

def fastcomm_analysis_from_dir(run_dir: str, state_key: str, sample_key: str = "") -> Dict[str, Any]:
    """Build `meta['fastcomm_analysis']` from a finished fastComm run.

    Mirrors the success block the pipeline writes at flask/pipeline.py:1885. Keys and
    file names are the ones `altanalyze3.components.fastComm.cli score` writes.
    """
    run = Path(run_dir)
    scores = run / "fastcomm_scores.tsv"
    if not scores.is_file():
        return {"enabled": False, "status": "unavailable",
                "message": f"No fastComm scores at {scores}", "state_key": state_key}
    populations: List[str] = []
    frame = pd.read_csv(scores, sep="\t", usecols=["sender_state", "receiver_state"])
    populations = sorted({str(v).strip() for v in
                          pd.concat([frame["sender_state"], frame["receiver_state"]]) if str(v).strip()})
    summary: Dict[str, Any] = {}
    summary_path = run / "summary.json"
    if summary_path.is_file():
        with open(summary_path, "r") as fh:
            summary = json.load(fh)
    analysis: Dict[str, Any] = {
        "enabled": True,
        "status": "completed",
        "state_key": state_key,
        "sample_key": sample_key,
        "populations": populations,
        "scores_tsv": str(scores),
        "state_pair_tsv": str(run / "state_pair_summary.tsv"),
        "state_expression_tsv": str(run / "state_expression.tsv"),
        "significance_threshold": 0.25,
        "summary": summary,
        "n_rows": int(frame.shape[0]),
    }
    split_long = run / "per_sample" / "split_scores_long.tsv"
    if split_long.is_file():
        analysis["per_sample"] = {
            "split_scores_long_tsv": str(split_long),
            "output_dir": str(run / "per_sample"),
            "split_key": sample_key,
            "state_key": state_key,
        }
    return analysis
