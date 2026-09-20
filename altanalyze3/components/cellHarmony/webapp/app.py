from __future__ import annotations

import io
import json
import math
import logging
import os
import re
import shutil
import threading
import time
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from fastapi import FastAPI, File, Form, HTTPException, Query, Request, UploadFile
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel, Field

from altanalyze3.components.cellHarmony.flask import pipeline as pipeline_mod
from altanalyze3.components.cellHarmony.flask.job_manager import JobStore
from altanalyze3.components.cellHarmony.flask.tasks import JobRunner
from altanalyze3.components.visualization import approximate_umap as approx_mod
from altanalyze3.components.rna2metabolite import annotations as metabolite_annotations

from .config import BASE_DIR, load_config
from .grn_data import UploadedGrnData, completed_differentials, comparison_entry
from altanalyze3.components.cellHarmony import grn_analysis as gnet

_POOLED_OVERALL_LABEL = "Pooled overall"
_GO_ELITE_HIGHLIGHT_KEYWORDS = (
    "cell cycle",
    "mitotic",
    "splicing",
    "mrna processing",
    "proliferation",
    "cytokine",
    "death",
    "chromatin",
    "lipid",
    "circadian",
    "tp53",
    "wnt",
    "tgf",
    "tnf",
    "granule",
)


class QCSettings(BaseModel):
    min_genes: int = 500
    min_counts: int = 1000
    min_cells: int = 0
    mit_percent: int = 15
    align_cutoff: float = 0.4
    ambient_correction: str = Field(default="no", pattern="^(no|yes)$")
    impute_modality: Optional[str] = "none"          # legacy single-select
    impute_modalities: Optional[List[str]] = None    # multi-select
    marker_render_heatmap: bool = False
    marker_write_svg: bool = True
    marker_heatmap_dpi: Optional[float] = Field(default=None, gt=0, allow_inf_nan=False)
    marker_cells_per_cluster: int = Field(default=100, ge=0)


class JobConfigSettings(BaseModel):
    species: str
    reference: str
    ambient_option: Optional[str] = None


class DifferentialSettings(BaseModel):
    modality: Optional[str] = "rna"
    population_col: str
    sample_field: Optional[str] = None
    group1_samples: List[str]
    group2_samples: List[str]
    comparison_type: str = "cells"


class ApproximateUMAPRequest(BaseModel):
    query: Optional[str] = None
    reference: Optional[str] = None
    query_clusters_tsv: Optional[str] = None
    reference_coords_tsv: Optional[str] = None
    reference_clusters_tsv: Optional[str] = None
    query_cluster_key: str
    reference_cluster_key: Optional[str] = None
    umap_key: str = "X_umap"
    jitter: float = 0.05
    num_reference_cells: int = 1
    random_state: Optional[int] = None
    custom_colors_tsv: Optional[str] = None
    restrict_obs_field: Optional[str] = None
    restrict_obs_value: Optional[str] = None
    output_prefix: Optional[str] = None
    outdir: str
    output_h5ad: Optional[str] = None
    output_pdf: Optional[str] = None
    save_updated_h5ad: bool = False
    verbose: bool = False


class ChatRequest(BaseModel):
    question: str = ""


class ClientLogRequest(BaseModel):
    message: str


def _secure_filename(filename: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", filename).strip("._")
    return cleaned or "upload"


def _require_within(candidate: Path, root: Path, label: str) -> Path:
    """Resolve `candidate` and refuse it unless it stays under `root`.

    Every request-supplied path reaches this before it is opened or written.
    `Path.resolve()` collapses ".." and follows symlinks, so the comparison is
    made on the real location rather than on the text the caller sent.
    """
    resolved = Path(candidate).resolve()
    base = Path(root).resolve()
    if resolved != base and base not in resolved.parents:
        raise HTTPException(status_code=400, detail=f"{label} must stay inside {base}.")
    return resolved


def _require_within_any(candidate: str, roots: List[Path], label: str) -> Path:
    resolved = Path(candidate).resolve()
    for root in roots:
        base = Path(root).resolve()
        if resolved == base or base in resolved.parents:
            return resolved
    raise HTTPException(status_code=400, detail=f"{label} is outside the permitted directories.")


def _tool_input_roots(app: FastAPI) -> List[Path]:
    """Where /api/tools/* may read from: job storage and the reference tree."""
    cfg = app.state.config
    return [Path(cfg["JOB_STORAGE"]), Path(cfg["REFERENCE_REGISTRY"]).parent]


def _normalize_root_path(root_path: Optional[str]) -> str:
    value = str(root_path or "").strip()
    if not value or value == "/":
        return ""
    return "/" + value.strip("/")


def _with_root_path(root_path: str, path: str) -> str:
    normalized_path = path if path.startswith("/") else f"/{path}"
    return f"{root_path}{normalized_path}" if root_path else normalized_path


def _is_api_request(request: Request) -> bool:
    request_path = str(request.url.path or "")
    scope_path = str(request.scope.get("path") or "")
    root_path = _normalize_root_path(request.scope.get("root_path") or "")
    api_prefix = _with_root_path(root_path, "/api/")
    return (
        request_path.startswith(api_prefix)
        or request_path.startswith("/api/")
        or scope_path.startswith("/api/")
    )


def _flatten_expr(values) -> np.ndarray:
    if sp.issparse(values):
        return np.asarray(values.todense()).ravel()
    return np.asarray(values).ravel()


def _is_finite_number(value: object) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def _normalize_h5ad_compression(value: Optional[str]) -> Optional[str]:
    raw = str(value or "").strip().lower()
    if not raw or raw in {"none", "off", "false", "null"}:
        return None
    if raw in {"lzf", "gzip"}:
        return raw
    return "lzf"


from altanalyze3.components.cellHarmony.modalities import (
    MODALITY_DEFINITIONS as _DEFAULT_MODALITY_DEFINITIONS,
    normalize_modality_id as _normalize_modality_id,
    modality_artifacts as _modality_artifacts,
)


def _modalities_state(meta: Dict) -> Dict[str, object]:
    stored = dict(meta.get("modalities") or {})
    available_raw = list(stored.get("available") or [])
    artifacts = _modality_artifacts(meta)
    if artifacts.get("grn_tf", {}).get("legacy_enrichment"):
        available_raw = [dict(entry, **_DEFAULT_MODALITY_DEFINITIONS["grn"])
                         if entry.get("id") == "grn" else entry for entry in available_raw]
        available_raw.append(dict(_DEFAULT_MODALITY_DEFINITIONS["grn_tf"],
                                  label="TF enrichment (legacy)", supports_differential=False))
    available: List[Dict[str, object]] = []
    seen: set[str] = set()
    for entry in available_raw:
        if not isinstance(entry, dict):
            continue
        modality_id = _normalize_modality_id(entry.get("id"), default="")
        if not modality_id or modality_id in seen:
            continue
        base = dict(_DEFAULT_MODALITY_DEFINITIONS.get(modality_id, {"id": modality_id, "label": modality_id.upper(), "feature_label": "feature"}))
        base.update(entry)
        base["id"] = modality_id
        available.append(base)
        seen.add(modality_id)
    if "rna" not in seen:
        available.insert(0, dict(_DEFAULT_MODALITY_DEFINITIONS["rna"]))
        seen.add("rna")
    default_modality = _normalize_modality_id(stored.get("default"), default="rna")
    if default_modality not in seen:
        default_modality = "rna"
    return {"default": default_modality, "available": available}


def _modality_definition(meta: Dict, modality: object) -> Dict[str, object]:
    normalized = _normalize_modality_id(modality)
    for entry in _modalities_state(meta)["available"]:
        if _normalize_modality_id(entry.get("id"), default="") == normalized:
            return dict(entry)
    # A DEG-only modality is absent from the bundle's modality list, because the bundle
    # stores no feature matrix for it. Cell communication is one: it arrives through
    # --deg-modality alone. Falling straight through to RNA gave it RNA's feature label,
    # so the panel read "gene" where it means "ligand-receptor interaction", and RNA's
    # supports_* flags with it. Use the modality's OWN definition when one exists.
    if normalized in _DEFAULT_MODALITY_DEFINITIONS:
        return dict(_DEFAULT_MODALITY_DEFINITIONS[normalized])
    return dict(_DEFAULT_MODALITY_DEFINITIONS["rna"])


def _modality_marker_analysis(meta: Dict, modality: object) -> Dict[str, object]:
    normalized = _normalize_modality_id(modality)
    by_modality = meta.get("marker_analysis_by_modality") or {}
    if isinstance(by_modality, dict):
        entry = by_modality.get(normalized)
        if normalized == "grn_tf" and _modality_artifacts(meta).get("grn_tf", {}).get("legacy_enrichment"):
            entry = by_modality.get("grn")
        if isinstance(entry, dict):
            return dict(entry)
    if normalized == "rna":
        return dict(meta.get("marker_analysis") or {})
    return {}


def _modality_h5ad_path(meta: Dict, modality: object) -> Path:
    normalized = _normalize_modality_id(modality)
    if normalized == "rna":
        raw_path = str(meta.get("artifacts", {}).get("combined_h5ad", "")).strip()
    else:
        raw_path = str((_modality_artifacts(meta).get(normalized) or {}).get("h5ad", "")).strip()
    if not raw_path:
        raise FileNotFoundError(f"AnnData output unavailable for modality '{normalized}'.")
    path = Path(raw_path)
    if not path.exists():
        raise FileNotFoundError(f"AnnData output unavailable for modality '{normalized}'.")
    return path


def _normalize_gene_token(value: object) -> str:
    return re.sub(r"[\s_.-]+", "", str(value or "").strip()).upper()


def _resolve_gene_name(candidates, requested_gene: str) -> Optional[str]:
    requested = str(requested_gene or "").strip()
    if not requested:
        return None
    if requested in candidates:
        return requested

    requested_norm = _normalize_gene_token(requested)
    normalized_map: Dict[str, str] = {}
    for candidate in candidates:
        candidate_str = str(candidate)
        normalized_map.setdefault(_normalize_gene_token(candidate_str), candidate_str)
    return normalized_map.get(requested_norm)


def _configure_matplotlib_pdf_style() -> None:
    plt.rcParams["axes.linewidth"] = 0.5
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    plt.rcParams["figure.facecolor"] = "white"


def _clear_directory_contents(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for child in path.iterdir():
        if child.is_dir():
            shutil.rmtree(child, ignore_errors=False)
        else:
            child.unlink(missing_ok=True)


def _build_marker_heatmap_viewer_html(job_id: str, root_path: str) -> str:
    job_id_json = json.dumps(str(job_id))
    root_path_json = json.dumps(str(root_path or ""))
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Marker Heatmap</title>
  <link rel="stylesheet" href="https://software.broadinstitute.org/morpheus/css/morpheus-latest.min.css" />
  <style>
    html, body, #morpheus-target {{
      width: 100%;
      height: 100%;
      margin: 0;
      padding: 0;
      overflow: hidden;
      background: #ffffff;
    }}
    body {{
      font-family: Arial, sans-serif;
    }}
    #fallback-message {{
      display: none;
      box-sizing: border-box;
      width: 100%;
      height: 100%;
      padding: 24px;
      color: #475569;
      font-size: 14px;
      line-height: 1.4;
    }}
  </style>
</head>
<body>
  <div id="morpheus-target"></div>
  <div id="fallback-message"></div>
  <script src="https://software.broadinstitute.org/morpheus/js/morpheus-external-latest.min.js"></script>
  <script src="https://software.broadinstitute.org/morpheus/js/morpheus-latest.min.js"></script>
  <script>
    (function () {{
      const jobId = {job_id_json};
      const rootPath = {root_path_json};
      function apiPath(path) {{
        const normalizedPath = path.startsWith("/") ? path : `/${{path}}`;
        return rootPath ? `${{rootPath}}${{normalizedPath}}` : normalizedPath;
      }}
      async function logClient(message) {{
        try {{
          await fetch(apiPath(`/api/jobs/${{jobId}}/client-log`), {{
            method: "POST",
            headers: {{ "Content-Type": "application/json" }},
            body: JSON.stringify({{ message: `[marker-viewer] ${{message}}` }}),
          }});
        }} catch (err) {{
          console.debug("Viewer log failed", err);
        }}
      }}
      function showFallback(message) {{
        const target = document.getElementById("morpheus-target");
        const fallback = document.getElementById("fallback-message");
        if (target) {{
          target.style.display = "none";
        }}
        fallback.style.display = "block";
        fallback.textContent = message;
      }}
      function loadScript(src) {{
        return new Promise((resolve, reject) => {{
          const existing = Array.from(document.scripts).find((script) => script.src === src);
          if (existing) {{
            if (existing.dataset.loaded === "true") {{
              resolve();
              return;
            }}
            existing.addEventListener("load", () => resolve(), {{ once: true }});
            existing.addEventListener("error", () => reject(new Error(`Failed to load ${{src}}`)), {{ once: true }});
            return;
          }}
          const script = document.createElement("script");
          script.src = src;
          script.async = true;
          script.crossOrigin = "anonymous";
          script.addEventListener("load", () => {{
            script.dataset.loaded = "true";
            resolve();
          }}, {{ once: true }});
          script.addEventListener("error", () => reject(new Error(`Failed to load ${{src}}`)), {{ once: true }});
          document.head.appendChild(script);
        }});
      }}
      async function ensurePdfLibraries() {{
        if (!(window.jspdf && window.jspdf.jsPDF)) {{
          const sources = [
            "https://cdn.jsdelivr.net/npm/jspdf@2.5.1/dist/jspdf.umd.min.js",
            "https://unpkg.com/jspdf@2.5.1/dist/jspdf.umd.min.js",
          ];
          let lastError = null;
          for (const source of sources) {{
            try {{
              await loadScript(source);
              if (window.jspdf && window.jspdf.jsPDF) {{
                break;
              }}
            }} catch (err) {{
              lastError = err;
            }}
          }}
          if (!(window.jspdf && window.jspdf.jsPDF)) {{
            throw lastError || new Error("Unable to load jsPDF.");
          }}
        }}
        if (!window.html2canvas) {{
          const sources = [
            "https://cdn.jsdelivr.net/npm/html2canvas@1.4.1/dist/html2canvas.min.js",
            "https://unpkg.com/html2canvas@1.4.1/dist/html2canvas.min.js",
          ];
          let lastError = null;
          for (const source of sources) {{
            try {{
              await loadScript(source);
              if (window.html2canvas) {{
                break;
              }}
            }} catch (err) {{
              lastError = err;
            }}
          }}
          if (!window.html2canvas) {{
            throw lastError || new Error("Unable to load html2canvas.");
          }}
        }}
      }}
      async function exportPdf(filename) {{
        await ensurePdfLibraries();
        const target = document.getElementById("morpheus-target");
        const canvas = await window.html2canvas(target, {{
          backgroundColor: "#ffffff",
          scale: 2,
          useCORS: true,
          logging: false,
        }});
        const dataUrl = canvas.toDataURL("image/png");
        const width = canvas.width || target.clientWidth || 1;
        const height = canvas.height || target.clientHeight || 1;
        const orientation = width >= height ? "landscape" : "portrait";
        const pdf = new window.jspdf.jsPDF({{
          orientation,
          unit: "pt",
          format: "a4",
          compress: true,
        }});
        const pageWidth = pdf.internal.pageSize.getWidth();
        const pageHeight = pdf.internal.pageSize.getHeight();
        const margin = 18;
        const availableWidth = pageWidth - margin * 2;
        const availableHeight = pageHeight - margin * 2;
        const scale = Math.min(availableWidth / width, availableHeight / height, 1);
        const renderWidth = width * scale;
        const renderHeight = height * scale;
        const x = (pageWidth - renderWidth) / 2;
        const y = (pageHeight - renderHeight) / 2;
        pdf.addImage(dataUrl, "PNG", x, y, renderWidth, renderHeight, undefined, "FAST");
        pdf.save(String(filename || "marker_heatmap.pdf"));
        await logClient(`marker heatmap PDF exported filename=${{filename || "marker_heatmap.pdf"}}`);
      }}
      window.addEventListener("message", (event) => {{
        if (event.origin !== window.location.origin) {{
          return;
        }}
        const payload = event.data || {{}};
        if (payload.type !== "marker-heatmap-export-pdf") {{
          return;
        }}
        exportPdf(payload.filename).catch((err) => {{
          logClient(`marker heatmap PDF export failed: ${{err && err.message ? err.message : err}}`);
          showFallback(`Unable to export the current marker heatmap PDF. ${{err && err.message ? err.message : err}}`);
        }});
      }});
      window.addEventListener("error", (event) => {{
        const details = [
          event.message || "unknown error",
          event.filename || "-",
          event.lineno || 0,
          event.colno || 0,
        ].join(" | ");
        logClient(`window.error ${{details}}`);
      }});
      window.addEventListener("unhandledrejection", (event) => {{
        const reason = event.reason && event.reason.message ? event.reason.message : String(event.reason || "unknown rejection");
        logClient(`window.unhandledrejection ${{reason}}`);
      }});
      async function fetchDatasetText(datasetUrl) {{
        try {{
          const headResp = await fetch(datasetUrl, {{ method: "HEAD", cache: "no-store" }});
          await logClient(
            `dataset HEAD status=${{headResp.status}} ok=${{headResp.ok}} content_type=${{headResp.headers.get("content-type") || "-"}} content_length=${{headResp.headers.get("content-length") || "-"}}`
          );
        }} catch (err) {{
          await logClient(`dataset HEAD failed: ${{err && err.message ? err.message : err}}`);
        }}
        try {{
          const getResp = await fetch(datasetUrl, {{ method: "GET", cache: "no-store" }});
          await logClient(
            `dataset GET status=${{getResp.status}} ok=${{getResp.ok}} content_type=${{getResp.headers.get("content-type") || "-"}} content_length=${{getResp.headers.get("content-length") || "-"}}`
          );
          if (!getResp.ok) {{
            throw new Error(`Dataset GET returned ${{getResp.status}}.`);
          }}
          const text = await getResp.text();
          await logClient(`dataset text length=${{text.length}}`);
          return text;
        }} catch (err) {{
          await logClient(`dataset GET failed: ${{err && err.message ? err.message : err}}`);
          throw err;
        }}
      }}
      async function render() {{
        const datasetUrl = `${{window.location.origin}}${{apiPath(`/api/jobs/${{jobId}}/marker/heatmap.tsv`)}}${{window.location.search || ""}}`;
        await logClient(
          `Initializing local Morpheus viewer. dataset_url=${{datasetUrl}} protocol=${{window.location.protocol}} ready_state=${{document.readyState}} ua=${{navigator.userAgent}}`
        );
        const datasetText = await fetchDatasetText(datasetUrl);
        if (!window.morpheus || !window.morpheus.HeatMap) {{
          await logClient(
            `Morpheus scripts did not load. morpheus_present=${{Boolean(window.morpheus)}} heatmap_present=${{Boolean(window.morpheus && window.morpheus.HeatMap)}}`
          );
          showFallback("Morpheus assets failed to load.");
          return;
        }}
        try {{
          const blobUrl = URL.createObjectURL(new Blob([datasetText], {{ type: "text/tab-separated-values" }}));
          await logClient(`Created blob URL for marker heatmap dataset.`);
          new window.morpheus.HeatMap({{
            el: document.getElementById("morpheus-target"),
            dataset: blobUrl,
            rowSize: 14,
            columnSize: 10,
            drawGrid: false,
            rows: [{{ field: "id", display: ["text"] }}],
            columns: [{{ field: "id", display: ["text"] }}],
            colorScheme: {{
              scalingMode: "fixed",
              stepped: false,
              values: [-2, 0, 2],
              colors: ["#00f0ff", "#000000", "#ffff00"],
            }},
          }});
          await logClient("Local Morpheus viewer initialized.");
        }} catch (err) {{
          await logClient(`Local Morpheus viewer failed: ${{err && err.message ? err.message : err}}`);
          showFallback(`Unable to initialize Morpheus viewer. ${{err && err.message ? err.message : err}}`);
        }}
      }}
      render();
    }}());
  </script>
</body>
</html>
"""


_PAIRED_COLOR_STOPS = [
    (0.0, (0.6509804129600525, 0.8078431487083435, 0.8901960849761963)),
    (0.09090909090909091, (0.12156862765550613, 0.47058823704719543, 0.7058823704719543)),
    (0.18181818181818182, (0.6980392336845398, 0.8745098114013672, 0.5411764979362488)),
    (0.2727272727272727, (0.20000000298023224, 0.6274510025978088, 0.1725490242242813)),
    (0.36363636363636365, (0.9843137264251709, 0.6039215922355652, 0.6000000238418579)),
    (0.45454545454545453, (0.8901960849761963, 0.10196078568696976, 0.10980392247438431)),
    (0.5454545454545454, (0.9921568632125854, 0.7490196228027344, 0.43529412150382996)),
    (0.6363636363636364, (1.0, 0.49803921580314636, 0.0)),
    (0.7272727272727273, (0.7921568751335144, 0.6980392336845398, 0.8392156958580017)),
    (0.8181818181818182, (0.4156862795352936, 0.239215686917305, 0.6039215922355652)),
    (0.9090909090909091, (1.0, 1.0, 0.6000000238418579)),
    (1.0, (0.6941176652908325, 0.3490196168422699, 0.1568627506494522)),
]


def _interpolate_paired_color(value: float) -> tuple[float, float, float]:
    clipped = max(0.0, min(1.0, float(value)))
    for index in range(1, len(_PAIRED_COLOR_STOPS)):
        right_pos, right_rgb = _PAIRED_COLOR_STOPS[index]
        left_pos, left_rgb = _PAIRED_COLOR_STOPS[index - 1]
        if clipped <= right_pos or index == len(_PAIRED_COLOR_STOPS) - 1:
            span = max(right_pos - left_pos, 1e-9)
            ratio = max(0.0, min(1.0, (clipped - left_pos) / span))
            return tuple(left_rgb[channel] + (right_rgb[channel] - left_rgb[channel]) * ratio for channel in range(3))
    return _PAIRED_COLOR_STOPS[-1][1]


def _custom_shuffle_indices(indices: list[int]) -> list[int]:
    shuffled: list[int] = []
    for index, value in enumerate(indices):
        if value not in shuffled:
            shuffled.append(value)
        from_end = indices[len(indices) - 1 - index] if indices else value
        if from_end not in shuffled:
            shuffled.append(from_end)
        middle_index = int((index + len(indices)) / 2)
        from_middle = indices[middle_index] if middle_index < len(indices) else indices[-1]
        if from_middle not in shuffled:
            shuffled.append(from_middle)
    return shuffled


def _seeded_shuffle(items: list[int], seed: int = 0) -> list[int]:
    rng = np.random.default_rng(seed)
    values = list(items)
    for index in range(len(values) - 1, 0, -1):
        swap_index = int(rng.integers(0, index + 1))
        values[index], values[swap_index] = values[swap_index], values[index]
    return values


def _build_preview_palette(populations: list[str]) -> dict[str, tuple[float, float, float]]:
    ordered = list(populations)
    if len(ordered) <= 4:
        base = ["#ff0000", "#0000ff", "#ffff00", "#00aa00", "#ffffff", "#000000", "#ff00ff"]
        return {population: matplotlib.colors.to_rgb(base[index % len(base)]) for index, population in enumerate(ordered)}
    indices = _seeded_shuffle(_custom_shuffle_indices(list(range(len(ordered)))), 0)
    denominator = max(len(ordered) - 1, 1)
    colors = [_interpolate_paired_color(index / denominator) for index in indices]
    return {population: colors[index] for index, population in enumerate(ordered)}


def _resolve_output_path(base_dir: Path, candidate: Optional[str]) -> Optional[Path]:
    if not candidate:
        return None
    path = Path(candidate)
    return path if path.is_absolute() else base_dir / path


def _load_reference_registry(app: FastAPI) -> Dict:
    registry_path = Path(app.state.config["REFERENCE_REGISTRY"])
    if not registry_path.exists():
        return {"species": []}
    with registry_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _allowed_file(app: FastAPI, filename: str) -> bool:
    return "." in filename and filename.rsplit(".", 1)[1].lower() in app.state.config["ALLOWED_EXTENSIONS"]


def _job_resources(app: FastAPI) -> tuple[JobStore, JobRunner]:
    return app.state.job_store, app.state.job_runner



def _job_metadata_with_recovery(store, runner, job_id):
    # Published viewers replace the upload store with a bundle-backed store.
    # Their completed results have no upload worker or on-disk job.json to recover.
    if getattr(runner, "store", None) is store:
        return runner.recover_interrupted_differential(job_id)
    return store.get_job(job_id)


def _build_reference_preview_payload(app: FastAPI, species: str, reference_id: str) -> Dict:
    registry_path = Path(app.state.config["REFERENCE_REGISTRY"])
    reference_entry = pipeline_mod._lookup_reference(species, reference_id, registry_path)
    pipeline_mod._ensure_reference_fields(reference_entry)
    cluster_key = reference_entry.get("cluster_key", Path(reference_entry["states_tsv"]).stem)

    adata = approx_mod._load_reference_from_tsv(
        reference_entry["reference_coords_tsv"],
        reference_entry["reference_clusters_tsv"],
        umap_key="X_umap",
        cluster_key=cluster_key,
    )

    coords = np.asarray(adata.obsm["X_umap"])
    labels = adata.obs[cluster_key].astype(str).tolist()
    frame = pd.DataFrame(
        {
            "x": coords[:, 0],
            "y": coords[:, 1],
            "population": labels,
        }
    )
    frame = frame[
        frame["x"].map(_is_finite_number)
        & frame["y"].map(_is_finite_number)
        & frame["population"].astype(str).str.strip().ne("")
    ].copy()

    if frame.empty:
        raise HTTPException(status_code=500, detail="Reference preview contains no finite UMAP coordinates.")

    centroids = (
        frame.groupby("population", as_index=False)[["x", "y"]]
        .median()
        .sort_values("population")
    )

    return {
        "species": species,
        "reference": reference_id,
        "reference_label": reference_entry.get("label", reference_id),
        "cluster_key": cluster_key,
        "points": [
            {
                "x": float(row.x),
                "y": float(row.y),
                "population": str(row.population),
            }
            for row in frame.itertuples(index=False)
        ],
        "labels": [
            {
                "x": float(row.x),
                "y": float(row.y),
                "population": str(row.population),
            }
            for row in centroids.itertuples(index=False)
        ],
    }


def _get_cache_lock(app: FastAPI, lock_name: str, cache_key: str) -> threading.Lock:
    lock_map = getattr(app.state, lock_name, None)
    if not isinstance(lock_map, dict):
        lock_map = {}
        setattr(app.state, lock_name, lock_map)
    global_lock = getattr(app.state, "cache_registry_lock", None)
    if global_lock is None:
        global_lock = threading.Lock()
        app.state.cache_registry_lock = global_lock
    with global_lock:
        lock = lock_map.get(cache_key)
        if lock is None:
            lock = threading.Lock()
            lock_map[cache_key] = lock
        return lock


def _load_reference_adata(app: FastAPI, meta: Dict) -> Optional[ad.AnnData]:
    coords = meta.get("reference_coords_tsv")
    clusters = meta.get("reference_clusters_tsv")
    if not coords or not clusters:
        return None
    cluster_key = meta.get("reference_cluster_key") or meta.get("cluster_key")
    cache_key = json.dumps(
        {
            "coords": str(coords),
            "clusters": str(clusters),
            "cluster_key": str(cluster_key or ""),
        },
        sort_keys=True,
    )
    cache = getattr(app.state, "reference_adata_cache", None)
    if not isinstance(cache, dict):
        cache = {}
        app.state.reference_adata_cache = cache
    cached = cache.get(cache_key)
    if isinstance(cached, ad.AnnData):
        return cached

    lock = _get_cache_lock(app, "reference_adata_cache_locks", cache_key)
    with lock:
        cached = cache.get(cache_key)
        if isinstance(cached, ad.AnnData):
            return cached
        ref_adata = approx_mod._load_reference_from_tsv(
            coords,
            clusters,
            umap_key="X_umap",
            cluster_key=cluster_key,
        )
        cache[cache_key] = ref_adata
        return ref_adata


def _filter_qc_log_lines(lines: List[str]) -> List[str]:
    qc_markers = (
        "...performing QC",
        "[qc]",
        "reimported adata shape",
        "Cells remaining after min_genes",
        "Cells remaining after min_counts",
        "Cells remaining after mito-percent",
        "Job failed:",
    )
    return [line for line in lines if any(marker in line for marker in qc_markers)]


def _derive_live_pipeline_message(status: object, log_lines: List[str], fallback: object) -> str:
    normalized_status = str(status or "").strip().lower()
    if normalized_status != "processing":
      return str(fallback or "")
    stage_markers = (
        "Running fastComm receptor-ligand communication analysis.",
        "fastComm analysis complete:",
        "Running rna2adt ADT imputation.",
        "rna2adt ADT imputation complete.",
        "Running rna2lipid lipid imputation.",
        "rna2lipid lipid imputation complete.",
        "Running rna2metabolite imputation.",
        "rna2metabolite imputation complete.",
        "Running rna2lipid (AML) imputation.",
        "rna2lipid (AML) imputation complete.",
        "Running rna2grn imputation.",
        "rna2grn imputation complete.",
        "Ambient RNA correction",
        "ambient RNA correction",
        "ambient correction",
        "Running approximate UMAP placement.",
        "Exporting NetPerspective marker networks.",
        "Identifying cell-state marker genes.",
        "Running cellHarmony_lite pipeline.",
        "Reference metadata loaded.",
    )
    for line in reversed(log_lines):
        text = str(line or "")
        for marker in stage_markers:
            if marker in text:
                return marker
        if "Aligning cells to reference" in text:
            return "Aligning cells to reference..."
        if "Normalization steps" in text:
            return "Normalizing expression values..."
    return str(fallback or "")


def _job_sample_names(meta: Dict) -> List[str]:
    sample_names: List[str] = []
    for record in meta.get("files", []):
        sample_name = str(record.get("sample_name", "")).strip()
        if sample_name:
            sample_names.append(sample_name)
    return sample_names


def _differential_options(meta: Dict) -> Dict:
    stored = dict(meta.get("differential_options") or {})
    sample_names = stored.get("sample_names") or _job_sample_names(meta)
    sample_values = stored.get("sample_values") or {}
    population_columns = list(stored.get("population_columns", []))
    sample_fields = list(stored.get("sample_fields", []))
    max_group_values = max((len(values) for values in sample_values.values() if isinstance(values, list)), default=0)
    upload_profile = dict(stored.get("upload_profile") or pipeline_mod._upload_profile(meta))
    if upload_profile.get("single_h5ad"):
        combined_h5ad_raw = str(meta.get("artifacts", {}).get("combined_h5ad", "")).strip()
        if combined_h5ad_raw:
            combined_h5ad_path = Path(combined_h5ad_raw)
        else:
            combined_h5ad_path = None
        if combined_h5ad_path is not None and combined_h5ad_path.exists():
            rebuilt_fields, rebuilt_values = pipeline_mod._candidate_group_fields(
                combined_h5ad_path,
                preferred=["Library", "group", "sample"],
                max_categories=None,
            )
            if rebuilt_fields:
                sample_fields = rebuilt_fields
                sample_values = rebuilt_values
    default_population_col = stored.get("default_population_col") or meta.get("cluster_key")
    if not upload_profile.get("allow_alternate_population_fields"):
        population_columns = [entry for entry in population_columns if entry.get("value") == default_population_col]
        if not population_columns and default_population_col:
            population_columns = [{"value": default_population_col, "label": default_population_col, "n_categories": 0}]
    default_sample_field = stored.get("default_sample_field") or ""
    if default_sample_field and default_sample_field not in {str(entry.get("value") or "").strip() for entry in sample_fields}:
        default_sample_field = sample_fields[0]["value"] if sample_fields else ""
    if default_population_col and default_population_col not in {str(entry.get("value") or "").strip() for entry in population_columns}:
        default_population_col = population_columns[0]["value"] if population_columns else default_population_col
    max_group_values = max((len(values) for values in sample_values.values() if isinstance(values, list)), default=0)
    enabled = bool(
        stored.get(
            "enabled",
            bool(upload_profile.get("differential_eligible")),
        )
    )
    comparison_types = stored.get("comparison_types")
    if not isinstance(comparison_types, list) or not comparison_types:
        pseudobulk_allowed = bool(upload_profile.get("single_h5ad") or upload_profile.get("total_files", 0) >= 4)
        comparison_types = ["cells", "pseudobulk"] if pseudobulk_allowed else ["cells"]
    modalities_state = _modalities_state(meta)
    modalities_available = [entry for entry in modalities_state["available"]
                            if entry.get("supports_differential", True)]
    # fastComm being available does NOT mean a cell-communication DIFFERENTIAL exists.
    # For a precomputed bundle fastComm scores communication per cell state for the
    # Explore tab; running cellHarmony-differential on it is a separate analysis. Adding
    # the modality regardless left a reader able to pick "Cell communication" while the
    # panel stayed on the previous modality, because no contrast carries it. Offer it
    # only when the bundle actually holds such a contrast, or when this is an
    # interactive job that will compute one.
    fastcomm_analysis = meta.get("fastcomm_analysis") or {}
    if isinstance(fastcomm_analysis, dict) and fastcomm_analysis.get("enabled"):
        precomputed = (meta.get("scalable_viewer") or {}).get("deg_comparisons")
        if precomputed is None:
            has_contrast = True                     # an interactive job computes on demand
        else:
            has_contrast = any(
                _normalize_modality_id(str(c.get("modality") or ""), default="")
                == "cell_communication"
                for c in precomputed if isinstance(c, dict))
        if has_contrast and not any(
                _normalize_modality_id(entry.get("id"), default="") == "cell_communication"
                for entry in modalities_available):
            modalities_available.append(dict(_DEFAULT_MODALITY_DEFINITIONS["cell_communication"]))
    return {
        "enabled": enabled,
        "sample_names": sample_names,
        "sample_fields": sample_fields,
        "sample_values": sample_values,
        "default_sample_field": default_sample_field,
        "population_columns": population_columns,
        "default_population_col": default_population_col,
        "modalities": modalities_available,
        "default_modality": modalities_state["default"],
        "comparison_types": comparison_types,
        "upload_profile": upload_profile,
    }


def _build_differential_payload(app: FastAPI, job_id: str, meta: Dict, root_path: str = "") -> Dict:
    options = _differential_options(meta)
    differential = dict(meta.get("differential") or {})
    artifacts = differential.get("artifacts") or {}
    config = differential.get("config") or {}
    selected_modality = _normalize_modality_id(config.get("modality") or options.get("default_modality") or "rna")
    modality_info = _modality_definition(meta, selected_modality)
    result_populations: List[str] = []
    visualization_populations: Dict[str, List[str]] = {"summary": [], "heatmap": [], "volcano": [], "network": [], "go": [], "table": []}
    if str(differential.get("status") or "") == "completed":
        try:
            result_populations = _differential_result_populations(app, meta)
        except Exception:
            result_populations = []
        try:
            visualization_populations["heatmap"] = _differential_heatmap_populations(app, meta)
        except Exception:
            visualization_populations["heatmap"] = []
        visualization_populations["volcano"] = result_populations
        # The DEG-count chart draws every cell state at once. It carries the same list
        # so the selected state survives a switch into and out of the chart.
        visualization_populations["summary"] = result_populations
        if selected_modality == "cell_communication":
            visualization_populations["network"] = result_populations
            visualization_populations["table"] = result_populations
        if selected_modality != "cell_communication" and bool(modality_info.get("supports_differential_network")):
            visualization_populations["network"] = list(
                dict.fromkeys(
                    str(entry.get("population", "")).strip()
                    for entry in differential.get("networks") or []
                    if str(entry.get("population", "")).strip()
                )
            )
        if bool(modality_info.get("supports_differential_go")):
            visualization_populations["go"] = _differential_go_populations(app, meta)
    networks = []
    for entry in differential.get("networks") or []:
        network_id = str(entry.get("id", "")).strip()
        if not network_id:
            continue
        networks.append(
            {
                "id": network_id,
                "population": str(entry.get("population", network_id)),
                "png_url": _with_root_path(root_path, f"/api/jobs/{job_id}/differential/network/{network_id}?format=png"),
                "pdf_url": _with_root_path(root_path, f"/api/jobs/{job_id}/differential/network/{network_id}?format=pdf"),
                "tsv_url": _with_root_path(root_path, f"/api/jobs/{job_id}/differential/network/{network_id}?format=tsv"),
            }
        )

    status = str(differential.get("status") or ("idle" if options["enabled"] else "unavailable"))
    if selected_modality == "cell_communication":
        visualization_modes = [
            {"value": "summary", "label": "Differential counts"},
            {"value": "volcano", "label": "Score delta"},
            {"value": "network", "label": "Cell-state network"},
            {"value": "table", "label": "Top interaction table"},
        ]
    else:
        visualization_modes = [
            {"value": "summary", "label": "Differential counts"},
            {"value": "heatmap", "label": "Heatmap"},
            {"value": "volcano", "label": "Volcano"},
        ]
        # A view with no data must not be offered. `supports_differential_*` is a
        # MODALITY capability, not a statement about this run, so RNA kept both flags and
        # the menu listed Network and GO Terms for a contrast that ships neither.
        # Selecting one printed "No network data are available for this differential run."
        # Measured on the COPD bundle, contrast `cancer_vs_no_cancer`: 0 networks and 0 GO
        # populations. A run still in progress keeps its entry, because its artifacts
        # arrive later and the menu must not flicker while it computes.
        run_completed = status == "completed"
        if bool(modality_info.get("supports_differential_network")) and (
                networks or not run_completed):
            visualization_modes.append({"value": "network", "label": "Network"})
        if bool(modality_info.get("supports_differential_go")) and (
                visualization_populations.get("go") or not run_completed):
            visualization_modes.append({"value": "go", "label": "GO Terms"})
        if selected_modality in {"rna", "grn", "grn_tf"}:
            visualization_modes.append({"value": "integrated_network", "label": "Regulatory network"})
            visualization_populations["integrated_network"] = result_populations
        if selected_modality in {"lipid", "lipids", "metabolite"}:
            visualization_modes.append({"value": "integrated_pathway", "label": "Pathway"})
            visualization_populations["integrated_pathway"] = result_populations
    return {
        **options,
        "status": status,
        "progress": int(differential.get("progress", 0) or 0),
        "message": differential.get("message") or (
            "Select a cell-state field and numerator/denominator sample groups to run cellHarmony-differential."
            if options["enabled"]
            else "Differential gene analyses between biological groups (i.e., disease versus controls) are only enabled when two or more samples (multiple h5 files or a single h5ad) are uploaded for the job."
        ),
        "config": {**config, "modality": selected_modality},
        "completed_comparisons": [comparison_entry(key, run) for key, run in completed_differentials(meta).items()],
        "selected_modality": selected_modality,
        "feature_label": str(differential.get("feature_label") or modality_info.get("feature_label") or "gene"),
        "visualization_modes": visualization_modes,
        "run_id": differential.get("run_id"),
        "case_label": differential.get("case_label"),
        "control_label": differential.get("control_label"),
        "go_terms_included": bool(differential.get("go_terms_included")),
        "archive_url": _with_root_path(root_path, f"/api/jobs/{job_id}/differential/archive") if artifacts.get("archive") else None,
        "heatmap_svg_url": _with_root_path(root_path, f"/api/jobs/{job_id}/differential/heatmap?format=svg") if artifacts.get("heatmap_svg") else None,
        "heatmap_pdf_url": _with_root_path(root_path, f"/api/jobs/{job_id}/differential/heatmap?format=pdf") if artifacts.get("heatmap_pdf") else None,
        "heatmap_png_url": _with_root_path(root_path, f"/api/jobs/{job_id}/differential/heatmap?format=png") if artifacts.get("heatmap_png") else None,
        "cell_frequency_grouped_png_url": (
            _with_root_path(root_path, f"/api/jobs/{job_id}/differential/artifact/cell_frequency_grouped_png")
            if artifacts.get("cell_frequency_grouped_png")
            else None
        ),
        "cell_frequency_grouped_pdf_url": (
            _with_root_path(root_path, f"/api/jobs/{job_id}/differential/artifact/cell_frequency_grouped_pdf")
            if artifacts.get("cell_frequency_grouped_pdf")
            else None
        ),
        "cell_frequency_stacked_pdf_url": (
            _with_root_path(root_path, f"/api/jobs/{job_id}/differential/artifact/cell_frequency_stacked_pdf")
            if artifacts.get("cell_frequency_stacked_pdf")
            else None
        ),
        "cell_frequency_tsv_url": (
            _with_root_path(root_path, f"/api/jobs/{job_id}/differential/artifact/cell_frequency_tsv")
            if artifacts.get("cell_frequency_tsv")
            else None
        ),
        "visualization_populations": visualization_populations,
        "result_populations": result_populations,
        "default_result_population": (
            visualization_populations["heatmap"][0]
            if visualization_populations["heatmap"]
            else (result_populations[0] if result_populations else None)
        ),
        "networks": networks,
    }


def _get_differential_detail_table(app: FastAPI, meta: Dict) -> pd.DataFrame:
    cache_entry = _get_differential_cache_entry(app, meta)
    if isinstance(cache_entry.get("detail_table"), pd.DataFrame):
        return cache_entry["detail_table"]
    raw_path = cache_entry["signature"]["detail_path"]
    if not raw_path:
        raise FileNotFoundError("Differential DEG detail table is unavailable.")
    path = Path(raw_path)
    if not path.exists():
        raise FileNotFoundError("Differential DEG detail table is unavailable.")
    frame = pd.read_csv(path, sep="\t")
    if "population" not in frame.columns:
        raise ValueError("Differential DEG detail table is missing the population column.")
    frame["population"] = frame["population"].astype(str)
    frame["gene"] = frame["gene"].astype(str)
    # AN OPTIONAL DISPLAY MAP RENAMES THE FEATURE FOR EVERY READER OF THIS FRAME.
    #
    # A precomputed bundle publishes `meta["feature_display"]` for the chosen contrast's
    # modality (bundle_meta.build_meta), so the Differential Explorer names lipids
    # "16:0/22:4 Phosphatidylethanolamine" rather than "PE(16:0/22:4)", which is what
    # Explore already shows. `feature_key` keeps the key on every row for anything that
    # joins against the fold matrix or the expression store.
    #
    # A cellHarmony job sets no such key, so this whole block is skipped and the frame is
    # exactly what it was.
    display = meta.get("feature_display") or {}
    frame["feature_key"] = frame["gene"]
    if isinstance(display, dict) and display:
        frame["gene"] = frame["gene"].map(lambda g: display.get(g, g))
    cache_entry["detail_table"] = frame
    return frame


def _get_differential_heatmap_table(app: FastAPI, meta: Dict) -> pd.DataFrame:
    cache_entry = _get_differential_cache_entry(app, meta)
    if isinstance(cache_entry.get("heatmap_table"), pd.DataFrame):
        return cache_entry["heatmap_table"]
    path = _get_differential_artifact(meta, "heatmap_tsv")
    frame = pd.read_csv(path, sep="\t", index_col=0)
    cache_entry["heatmap_table"] = frame
    return frame


def _get_differential_go_table(app: FastAPI, meta: Dict) -> pd.DataFrame:
    cache_entry = _get_differential_cache_entry(app, meta)
    if isinstance(cache_entry.get("go_table"), pd.DataFrame):
        return cache_entry["go_table"]
    path = _get_differential_artifact(meta, "goelite_tsv")
    frame = pd.read_csv(path, sep="\t")
    if "population" in frame.columns:
        frame["population"] = frame["population"].astype(str)
    if "term_name" in frame.columns:
        frame["term_name"] = frame["term_name"].astype(str)
    cache_entry["go_table"] = frame
    return frame


def _split_population_direction(raw_value: object) -> tuple[str, str]:
    raw_text = str(raw_value or "").strip()
    if "__" not in raw_text:
        return raw_text, ""
    population, direction = raw_text.rsplit("__", 1)
    return population, direction


def _parse_heatmap_row_key(raw_key: object) -> Dict[str, str]:
    row_text = str(raw_key or "").strip()
    cluster_group, _, gene = row_text.partition(":")
    population, direction = _split_population_direction(cluster_group)
    return {
        "row_key": row_text,
        "population": population,
        "direction": direction or "unknown",
        "gene": gene or row_text,
    }


def _cell_state_populations(ordered: List[str]) -> List[str]:
    """Cell states only. `Pooled overall` is a whole-sample test, not a population.

    It answers a different question from a per-cell-state test, and offering it in the
    same selector made a reader treat it as a cell state: the network view then had
    nothing to draw and reported a missing modality instead of a missing cell state.
    """
    return [value for value in ordered if value and value != _POOLED_OVERALL_LABEL]


def _differential_result_populations(app: FastAPI, meta: Dict) -> List[str]:
    detailed = _get_differential_detail_table(app, meta)
    return _cell_state_populations(
        list(dict.fromkeys(detailed["population"].astype(str).tolist())))


def _differential_heatmap_populations(app: FastAPI, meta: Dict) -> List[str]:
    try:
        detailed = _get_differential_detail_table(app, meta)
    except Exception:
        return []
    return _cell_state_populations(
        list(dict.fromkeys(detailed["population"].astype(str).tolist())))


def _differential_go_populations(app: FastAPI, meta: Dict) -> List[str]:
    try:
        frame = _get_differential_go_table(app, meta)
    except Exception:
        return []
    ordered: List[str] = []
    seen: set[str] = set()
    for raw_value in frame.get("population", pd.Series(dtype=str)).astype(str).tolist():
        population, _ = _split_population_direction(raw_value)
        if population and population not in seen:
            seen.add(population)
            ordered.append(population)
    return ordered


def _differential_network_entry(meta: Dict, population: str) -> Optional[Dict]:
    differential = meta.get("differential", {}) or {}
    for entry in differential.get("networks", []) or []:
        if str(entry.get("population", "")).strip() == population:
            return entry
    return None


def _differential_group_labels(meta: Dict) -> tuple[str, str]:
    differential = meta.get("differential", {}) or {}
    config = differential.get("config", {}) or {}
    group1_samples = [str(value).strip() for value in config.get("group1_samples", []) if str(value).strip()]
    group2_samples = [str(value).strip() for value in config.get("group2_samples", []) if str(value).strip()]
    case_label = str(differential.get("case_label") or pipeline_mod._group_display_label(group1_samples) or "Group 1").strip()
    control_label = str(differential.get("control_label") or pipeline_mod._group_display_label(group2_samples) or "Group 2").strip()
    return case_label, control_label


def _differential_population_col(meta: Dict) -> str:
    differential = meta.get("differential", {}) or {}
    config = differential.get("config", {}) or {}
    population_col = str(config.get("population_col", "")).strip()
    if not population_col:
        raise ValueError("Differential population column is not configured.")
    return population_col


def _differential_gene_h5ad_path(meta: Dict) -> Path:
    for key in ("differentials_only_h5ad",):
        raw_path = str(meta.get("differential", {}).get("artifacts", {}).get(key, "")).strip()
        if raw_path and Path(raw_path).exists():
            return Path(raw_path)
    modality = _normalize_modality_id((meta.get("differential", {}).get("config", {}) or {}).get("modality"), default="rna")
    # GRN differential operates on edges, served from the edge-level h5ad.
    diff_h5ad = str(((meta.get("modality_artifacts") or {}).get(modality) or {}).get("differential_h5ad", "")).strip()
    if diff_h5ad and Path(diff_h5ad).exists():
        return Path(diff_h5ad)
    try:
        return _modality_h5ad_path(meta, modality)
    except FileNotFoundError:
        pass
    raise FileNotFoundError("Aligned AnnData output is unavailable for differential gene detail.")


_GRN_RAGGED_CACHE: Dict[str, Any] = {}


def _grn_ragged_values(meta: Dict, population: str, edge: str, obs_names):
    """Per-METACELL values of one GRN edge in one cell state, or None.

    rna2grn imputed the shipped GRN sidecar per CELL STATE, so a violin of one state
    held a single value repeated across its metacells: distinct=1, sd=0. The ragged
    store gives every metacell its own value, and holds only the edges significant in
    that state, which is what keeps it near 0.6 GB rather than 17.8 GB.

    The violin must stay on metacells. The donor-level pseudobulk that produced the
    statistics may not be published one point at a time, because a point is a donor.
    """
    root = str(meta.get("grn_ragged_dir") or "").strip()
    if not root or not population or not edge:
        return None
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in str(population))
    path = os.path.join(root, "%s.h5ad" % safe)
    if not os.path.isfile(path):
        return None
    entry = _GRN_RAGGED_CACHE.get(path)
    mtime = os.path.getmtime(path)
    if not entry or entry.get("mtime") != mtime:
        adata = ad.read_h5ad(path)
        entry = {"mtime": mtime, "adata": adata,
                 "edges": {str(e): i for i, e in enumerate(adata.var_names.astype(str))},
                 "rows": {str(b): i for i, b in enumerate(adata.obs_names.astype(str))}}
        _GRN_RAGGED_CACHE[path] = entry
    col = entry["edges"].get(str(edge))
    if col is None:
        return None
    column = entry["adata"].X[:, col]
    column = column.toarray().ravel() if hasattr(column, "toarray") else np.asarray(column).ravel()
    rows = entry["rows"]
    out = np.full(len(obs_names), np.nan, dtype=float)
    for i_, barcode in enumerate(np.asarray(obs_names, dtype=str)):
        j = rows.get(barcode)
        if j is not None:
            out[i_] = float(column[j])
    return out if np.isfinite(out).any() else None



def _has_edge_level_h5ad(meta: Dict) -> bool:
    """True when a real edge-level h5ad backs this run's GRN differential.

    A cellHarmony job writes `..._grn_edges.h5ad` and registers it, and the GRN gene
    detail must read THAT rather than the per-TF activity matrix. A precomputed bundle
    has no such file: its GRN modality store already holds the edges, so the expression
    cache is the edge matrix and the h5ad branch must not run.
    """
    artifacts = (meta.get("differential", {}) or {}).get("artifacts", {}) or {}
    candidates = [str(artifacts.get("differentials_only_h5ad", "")).strip()]
    grn = (meta.get("modality_artifacts") or {}).get("grn") or {}
    candidates.append(str(grn.get("differential_h5ad", "")).strip())
    return any(path and Path(path).is_file() for path in candidates)


def _open_gene_detail_adata(app: FastAPI, meta: Dict, gene: str) -> tuple[ad.AnnData, Path]:
    cache_entry = _get_differential_cache_entry(app, meta)
    primary = _differential_gene_h5ad_path(meta)
    adata = cache_entry.get("primary_adata")
    if not isinstance(adata, ad.AnnData):
        adata = ad.read_h5ad(primary)
        cache_entry["primary_adata"] = adata
        cache_entry["primary_var_names"] = adata.var_names.astype(str).to_numpy()
    if gene in set(cache_entry["primary_var_names"]):
        return adata, primary

    combined_path = str(meta.get("artifacts", {}).get("combined_h5ad", "")).strip()
    if combined_path and Path(combined_path).exists():
        fallback = Path(combined_path)
        fallback_adata = cache_entry.get("fallback_adata")
        if not isinstance(fallback_adata, ad.AnnData):
            fallback_adata = ad.read_h5ad(fallback)
            cache_entry["fallback_adata"] = fallback_adata
            cache_entry["fallback_var_names"] = fallback_adata.var_names.astype(str).to_numpy()
        if gene in set(cache_entry["fallback_var_names"]):
            return fallback_adata, fallback
    raise KeyError(f"Gene '{gene}' not found in the aligned AnnData output.")


def _close_backed_adata(adata: Optional[ad.AnnData]) -> None:
    if adata is None:
        return
    try:
        if getattr(adata, "isbacked", False):
            adata.file.close()
    except Exception:
        pass


def _invalidate_expression_cache(app: FastAPI, job_id: str) -> None:
    cache = getattr(app.state, "expression_cache", None)
    if isinstance(cache, dict):
        prefix = f"{str(job_id)}:"
        for key in list(cache.keys()):
            if key == str(job_id) or str(key).startswith(prefix):
                cache.pop(key, None)
    locks = getattr(app.state, "expression_cache_locks", None)
    if isinstance(locks, dict):
        prefix = f"{str(job_id)}:"
        for key in list(locks.keys()):
            if key == str(job_id) or str(key).startswith(prefix):
                locks.pop(key, None)


def _invalidate_differential_cache(app: FastAPI, job_id: str) -> None:
    cache = getattr(app.state, "differential_cache", None)
    if not isinstance(cache, dict):
        return
    entry = cache.pop(str(job_id), None)
    if not isinstance(entry, dict):
        return
    for key in ("primary_adata", "fallback_adata"):
        adata_obj = entry.get(key)
        if isinstance(adata_obj, ad.AnnData):
            try:
                _close_backed_adata(adata_obj)
            except Exception:
                pass


def _invalidate_marker_heatmap_cache(app: FastAPI, job_id: str) -> None:
    cache = getattr(app.state, "marker_heatmap_cache", None)
    if isinstance(cache, dict):
        prefix = f"{str(job_id)}:"
        for key in list(cache.keys()):
            if key == str(job_id) or str(key).startswith(prefix):
                cache.pop(key, None)


def _invalidate_fastcomm_cache(app: FastAPI, job_id: str) -> None:
    cache = getattr(app.state, "fastcomm_cache", None)
    if isinstance(cache, dict):
        cache.pop(str(job_id), None)


def _get_marker_heatmap_cache_entry(app: FastAPI, meta: Dict, modality: str = "rna",
                                    compact: bool = True) -> Dict[str, Any]:
    """The marker matrix for one modality, at one of two column densities.

    `compact` picks the run that kept 10 metacells per cell population, which is the
    default because a 441-column heatmap stays readable where 1,827 columns do not. The
    all-column run is a SEPARATE MarkerFinder run and not a superset: MarkerFinder ranks
    markers from the cells it is given, so the two share 760 of 1,248 rows, 60.9%, on the
    COPD v7 build. A bundle that ships only the compact run keeps serving it whatever
    `compact` says, because no second matrix exists to switch to.
    """
    job_id = str(meta.get("job_id") or "").strip()
    if not job_id:
        raise ValueError("Job metadata is missing a job_id.")

    normalized_modality = _normalize_modality_id(modality)
    marker_analysis = _modality_marker_analysis(meta, normalized_modality) or {}
    full_path_text = str(marker_analysis.get("heatmap_cache_full", "")).strip()
    if not compact and full_path_text:
        cache_path_text = full_path_text
        density = "all"
    else:
        cache_path_text = str(marker_analysis.get("heatmap_cache", "")).strip()
        density = "compact"
    heatmap_tsv_text = str(marker_analysis.get("heatmap_tsv", "")).strip()
    expression_tsv_text = str(marker_analysis.get("expression_tsv", "")).strip()

    cache_path = Path(cache_path_text) if cache_path_text else None
    heatmap_tsv_path = Path(heatmap_tsv_text) if heatmap_tsv_text else None
    expression_tsv_path = Path(expression_tsv_text) if expression_tsv_text else None

    cache_signature = {
        "modality": normalized_modality,
        "density": density,
        "heatmap_cache": str(cache_path or ""),
        "heatmap_tsv": str(heatmap_tsv_path or ""),
        "expression_tsv": str(expression_tsv_path or ""),
    }
    cache = app.state.marker_heatmap_cache
    cache_key = f"{job_id}:{normalized_modality}:{density}"
    existing = cache.get(cache_key)
    if existing and existing.get("signature") == cache_signature:
        return existing

    tsv_path: Optional[Path] = None
    if heatmap_tsv_path and heatmap_tsv_path.exists():
        tsv_path = heatmap_tsv_path
    elif expression_tsv_path and expression_tsv_path.exists():
        tsv_path = expression_tsv_path

    if cache_path and cache_path.exists():
        # The fold-matrix cache moved from .npz to .h5ad, at
        # visualization/marker_heatmap_h5ad.py:73. np.load cannot open an HDF5
        # file, so every MarkerHeatmap request raised
        # "Cannot load file containing pickled data when allow_pickle=False".
        # `_read_heatmap_cache` is the writer's own reader and it opens both
        # formats, so this calls it rather than parsing the cache a second way.
        # The import sits here because that module loads scanpy, which would add
        # seconds to app start-up for a view many jobs never open.
        from altanalyze3.components.visualization.marker_heatmap_h5ad import _read_heatmap_cache

        heatmap_df, row_clusters, column_clusters, _, _ = _read_heatmap_cache(cache_path)
        matrix = np.asarray(heatmap_df.to_numpy(), dtype=np.float32)
        col_barcodes = heatmap_df.columns.astype(str).to_numpy()
        genes = heatmap_df.index.astype(str).to_numpy()
        # The ids the .npz held were "cluster:gene" and "cluster:barcode". They
        # are rebuilt here so every reader downstream sees what it saw before.
        row_ids = np.asarray([f"{c}:{g}" for c, g in zip(row_clusters, genes)], dtype=str)
        col_ids = np.asarray([f"{c}:{b}" for c, b in zip(column_clusters, col_barcodes)], dtype=str)
        entry = {
            "signature": cache_signature,
            "modality": normalized_modality,
            "source": "cache",
            "source_path": cache_path,
            "matrix": matrix,
            "row_ids": row_ids,
            "col_ids": col_ids,
            "col_barcodes": col_barcodes,
            "tsv_path": tsv_path,
        }
        cache[cache_key] = entry
        return entry

    if tsv_path and tsv_path.exists():
        frame = pd.read_csv(tsv_path, sep="\t", index_col=0)
        col_ids = frame.columns.astype(str).to_numpy()
        entry = {
            "signature": cache_signature,
            "modality": normalized_modality,
            "source": "tsv",
            "source_path": tsv_path,
            "matrix": np.asarray(frame.to_numpy(), dtype=np.float32),
            "row_ids": frame.index.astype(str).to_numpy(),
            "col_ids": col_ids,
            "col_barcodes": np.asarray(
                [value.split(":", 1)[1] if ":" in value else value for value in col_ids],
                dtype=str,
            ),
            "tsv_path": tsv_path,
        }
        cache[cache_key] = entry
        return entry

    raise FileNotFoundError("Marker heatmap matrix unavailable.")


def _marker_heatmap_subset_to_tsv(
    matrix: np.ndarray,
    row_ids: np.ndarray,
    col_ids: np.ndarray,
) -> str:
    frame = pd.DataFrame(matrix, index=row_ids, columns=col_ids)
    buffer = io.StringIO()
    frame.to_csv(buffer, sep="\t", float_format="%.4g")
    return buffer.getvalue()


def _get_expression_cache(app: FastAPI, meta: Dict, modality: str = "rna") -> Dict[str, Any]:
    job_id = str(meta.get("job_id") or "").strip()
    if not job_id:
        raise ValueError("Job metadata is missing a job_id.")

    normalized_modality = _normalize_modality_id(modality)
    # Precomputed viewers own their expression stores and do not require the
    # original analysis H5AD to remain available after a bundle is published.
    bundle_cache = getattr(app.state.job_store, "get_expression_cache", None)
    if callable(bundle_cache):
        return bundle_cache(app, meta, normalized_modality)
    artifacts = meta.get("artifacts", {})
    h5ad_path = _modality_h5ad_path(meta, normalized_modality)

    cluster_key = meta.get("cluster_key")
    if not cluster_key:
        raise ValueError("Cluster assignments missing from AnnData output.")

    umap_path = artifacts.get("umap_coordinates")
    is_bundle=hasattr(app.state.job_store,'dataset')
    source_stamp=None
    if not is_bundle and Path(h5ad_path).is_file():
        stat=Path(h5ad_path).stat();source_stamp=(stat.st_mtime_ns,stat.st_size)
    cache = app.state.expression_cache
    cache_key = f"{job_id}:{normalized_modality}"
    cache_entry = cache.get(cache_key)
    if (
        cache_entry
        and cache_entry.get("h5ad_path") == str(h5ad_path)
        and cache_entry.get("umap_path") == str(umap_path or "")
        and cache_entry.get("cluster_key") == str(cluster_key)
        and cache_entry.get("modality") == normalized_modality
        and (is_bundle or cache_entry.get("source_stamp") == source_stamp)
    ):
        return cache_entry

    lock = _get_cache_lock(app, "expression_cache_locks", cache_key)
    with lock:
        cache_entry = cache.get(cache_key)
        if (
            cache_entry
            and cache_entry.get("h5ad_path") == str(h5ad_path)
            and cache_entry.get("umap_path") == str(umap_path or "")
            and cache_entry.get("cluster_key") == str(cluster_key)
            and cache_entry.get("modality") == normalized_modality
            and (is_bundle or cache_entry.get("source_stamp") == source_stamp)
        ):
            return cache_entry

        adata = ad.read_h5ad(h5ad_path)
        if cluster_key not in adata.obs.columns:
            raise ValueError("Cluster assignments missing from AnnData output.")

        obs_names = adata.obs_names.astype(str).to_numpy()
        populations = adata.obs[cluster_key].astype(str).to_numpy()
        umap_x = np.full(adata.n_obs, np.nan, dtype=float)
        umap_y = np.full(adata.n_obs, np.nan, dtype=float)

        if umap_path and Path(umap_path).exists():
            coords_df = pd.read_csv(umap_path, sep="\t")
            barcode_col = coords_df.columns[0]
            coords_df = coords_df.rename(columns={barcode_col: "CellBarcode", "umap_0": "UMAP1", "umap_1": "UMAP2"})
            coords_df["CellBarcode"] = coords_df["CellBarcode"].astype(str)
            coords_df["UMAP1"] = pd.to_numeric(coords_df["UMAP1"], errors="coerce")
            coords_df["UMAP2"] = pd.to_numeric(coords_df["UMAP2"], errors="coerce")
            coords_indexed = coords_df.drop_duplicates(subset=["CellBarcode"]).set_index("CellBarcode")
            reindexed = coords_indexed.reindex(obs_names)
            umap_x = reindexed["UMAP1"].to_numpy(dtype=float)
            umap_y = reindexed["UMAP2"].to_numpy(dtype=float)

        obs_filter_values: Dict[str, np.ndarray] = {}
        filter_fields: List[Dict[str, str]] = []
        obs = adata.obs
        upload_profile = pipeline_mod._upload_profile(meta)
        excluded_columns = {"UMAP-X", "UMAP-Y", "UMAP_leiden_X", "UMAP_leiden_Y"}
        for column in obs.columns:
            column_text = str(column).strip()
            if not column_text or column_text in excluded_columns or "umap" in column_text.lower():
                continue
            series = obs[column]
            if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
                continue
            values = (
                series.astype(str)
                .str.strip()
                .replace({"nan": "", "None": ""})
            )
            unique_values = [value for value in pd.Index(values[values != ""]).unique().tolist() if value]
            if len(unique_values) < 2 or len(unique_values) > 120:
                continue
            obs_filter_values[column_text] = values.to_numpy(dtype=str)
            filter_fields.append({"value": column_text, "label": column_text})

        sample_names = _job_sample_names(meta)
        sample_filter_field = None
        if sample_names:
            try:
                sample_filter_field, _ = pipeline_mod._resolve_samples_for_adata(adata, meta, sample_names)
            except Exception:
                sample_filter_field = next((column for column in ("Library", "group", "sample") if column in obs.columns), None)
        else:
            sample_filter_field = next((column for column in ("Library", "group", "sample") if column in obs.columns), None)

        if sample_filter_field and str(sample_filter_field) not in obs_filter_values:
            sample_values = (
                obs[sample_filter_field].astype(str).str.strip().replace({"nan": "", "None": ""})
            )
            unique_values = [value for value in pd.Index(sample_values[sample_values != ""]).unique().tolist() if value]
            if len(unique_values) >= 2:
                obs_filter_values[str(sample_filter_field)] = sample_values.to_numpy(dtype=str)
                filter_fields.insert(0, {"value": str(sample_filter_field), "label": str(sample_filter_field)})

        seen_fields = set()
        deduped_fields = []
        for entry in filter_fields:
            value = str(entry["value"])
            if value in seen_fields:
                continue
            seen_fields.add(value)
            deduped_fields.append(entry)
        filter_fields = deduped_fields

        default_secondary_field = str(cluster_key) if str(cluster_key) in obs_filter_values else ""
        if upload_profile.get("single_h5"):
            preferred_primary_field = str(cluster_key) if str(cluster_key) in obs_filter_values else (filter_fields[0]["value"] if filter_fields else "")
            filter_fields = [entry for entry in filter_fields if entry["value"] == preferred_primary_field] if preferred_primary_field else filter_fields[:1]
            filter_values = {
                field: [value for value in pd.Index(values).unique().tolist() if value]
                for field, values in obs_filter_values.items()
                if field == preferred_primary_field
            }
            default_primary_field = preferred_primary_field
            default_secondary_field = ""
        else:
            if default_secondary_field and default_secondary_field == str(sample_filter_field or ""):
                default_secondary_field = next(
                    (entry["value"] for entry in filter_fields if entry["value"] != str(sample_filter_field or "")),
                    "",
                )
            filter_values = {
                field: [value for value in pd.Index(values).unique().tolist() if value]
                for field, values in obs_filter_values.items()
            }
            default_primary_field = str(sample_filter_field or "")

        sample_labels = None
        if sample_filter_field and str(sample_filter_field) in obs.columns:
            sample_labels = (
                obs[sample_filter_field]
                .astype(str)
                .str.strip()
                .replace({"nan": "", "None": ""})
                .to_numpy(dtype=str)
            )

        cache_entry = {
            "job_id": job_id,
            "modality": normalized_modality,
            "h5ad_path": str(h5ad_path),
            "source_stamp": source_stamp,
            "umap_path": str(umap_path or ""),
            "cluster_key": str(cluster_key),
            "adata": adata,
            "obs_names": obs_names,
            "var_names": adata.var_names.astype(str).to_numpy(),
            "populations": populations,
            "sample_field": str(sample_filter_field or ""),
            "sample_labels": sample_labels,
            "umap_x": umap_x,
            "umap_y": umap_y,
            "obsm_keys": _obsm_embedding_keys(adata),
            "obs_filter_values": obs_filter_values,
            "display_filters_meta": {
                "fields": filter_fields,
                "values": filter_values,
                "default_primary_field": default_primary_field,
                "default_secondary_field": default_secondary_field,
                "show_secondary": bool(upload_profile.get("show_secondary_display_filter", True)),
            },
        }
        cache[cache_key] = cache_entry
        return cache_entry


def _get_differential_cache_entry(app: FastAPI, meta: Dict) -> Dict[str, Any]:
    job_id = str(meta.get("job_id") or "").strip()
    if not job_id:
        raise ValueError("Job metadata is missing a job_id.")

    differential = meta.get("differential", {}) or {}
    artifacts = differential.get("artifacts", {}) or {}
    detail_path = next(
        (str(path).strip() for key, path in artifacts.items() if key.startswith("DEG_detailed_") and str(path).strip()),
        "",
    )
    heatmap_path = str(artifacts.get("heatmap_tsv", "")).strip()
    fold_matrix_path = str(artifacts.get("fold_matrix_tsv", "")).strip()
    go_path = str(artifacts.get("goelite_tsv", "")).strip()
    primary_h5ad_path = str(artifacts.get("differentials_only_h5ad", "")).strip()
    fallback_h5ad_path = str(meta.get("artifacts", {}).get("combined_h5ad", "")).strip()
    network_paths = {
        str(entry.get("population", "")).strip(): str(entry.get("tsv", "")).strip()
        for entry in differential.get("networks", []) or []
        if str(entry.get("population", "")).strip() and str(entry.get("tsv", "")).strip()
    }

    cache = app.state.differential_cache
    entry = cache.get(job_id)
    current_signature = {
        "detail_path": detail_path,
        "heatmap_path": heatmap_path,
        "fold_matrix_path": fold_matrix_path,
        "go_path": go_path,
        "primary_h5ad_path": primary_h5ad_path,
        "fallback_h5ad_path": fallback_h5ad_path,
        "network_paths": network_paths,
    }
    if entry and entry.get("signature") == current_signature:
        return entry

    _invalidate_differential_cache(app, job_id)
    entry = {
        "signature": current_signature,
        "detail_table": None,
        "heatmap_table": None,
        "fold_matrix": None,
        "go_table": None,
        "network_tables": {},
        "primary_adata": None,
        "fallback_adata": None,
        "primary_var_names": None,
        "fallback_var_names": None,
    }
    cache[job_id] = entry
    return entry


def _get_differential_fold_matrix(app: FastAPI, meta: Dict) -> pd.DataFrame:
    """The complete gene x cell-state log2 fold matrix the differential produced.

    Every pseudobulk contrast yields a fold in every cell state, significant or not.
    `fold_matrix_tsv` holds that matrix; flask/pipeline.py writes it from
    de_store['fold_matrix'] (cellHarmony_differential.py:1419). A run made before that
    artifact existed has no such file, and the caller falls back to the heatmap TSV.
    """
    cache_entry = _get_differential_cache_entry(app, meta)
    if isinstance(cache_entry.get("fold_matrix"), pd.DataFrame):
        return cache_entry["fold_matrix"]
    path = _get_differential_artifact(meta, "fold_matrix_tsv")
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = [str(label) for label in frame.columns]
    cache_entry["fold_matrix"] = frame
    return frame


def _differential_fold_rows(app: FastAPI, meta: Dict) -> tuple[Dict[str, List[float]], List[str], str]:
    """gene -> its log2 fold in every cell state, the column order, and the source file.

    Keyed by gene, never by (population, gene). A heatmap TSV row label names the
    BLOCK the differential placed a gene in - `coreg_global__up:CLUH`,
    `AT1__down:SFTPC` - and each gene sits in one block only. Keying by that block
    population missed every gene whose block was not the selected cell state, and the
    caller then read 0.0 from the significant-only detail table for all 45 columns.
    """
    # Cached per contrast. Rebuilding this for every heatmap request made the
    # cell-communication modality unusable: its fold matrix holds 126,726 interactions over
    # 55 cell states, and the old row-at-a-time loop below built about 7 million Python
    # floats each time, pinning one core at 99% and 3.1 GB for minutes.
    cache_entry = _get_differential_cache_entry(app, meta)
    cached = cache_entry.get("fold_rows")
    if cached is not None:
        return cached

    result: tuple[Dict[str, List[float]], List[str], str] = ({}, [], "detail_table")
    for source, loader in (
        ("fold_matrix_tsv", _get_differential_fold_matrix),
        ("heatmap_tsv", _get_differential_heatmap_table),
    ):
        try:
            frame = loader(app, meta)
        except Exception:
            continue
        if frame is None or frame.empty:
            continue
        labels = [str(label).split(":", 1)[0] for label in frame.columns]
        if len(set(labels)) != len(frame.columns):
            labels = [str(label) for label in frame.columns]
        # Vectorised: coerce every column once, then read whole rows out of one array.
        # Preserve unreported folds as missing; zero would imply a tested null change.
        numeric = frame.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        numeric = np.where(np.isfinite(numeric), numeric, np.nan)
        if source == "fold_matrix_tsv":
            keys = [str(label) for label in frame.index]
        else:
            keys = [_parse_heatmap_row_key(label)["gene"] for label in frame.index]
        rows: Dict[str, List[float]] = {}
        for position, gene in enumerate(keys):
            if not gene or gene in rows:
                continue
            rows[gene] = numeric[position].tolist()
        if rows:
            result = (rows, labels, source)
            break
    cache_entry["fold_rows"] = result
    return result


def _build_differential_heatmap_payload(app: FastAPI, meta: Dict, population: str) -> Dict:
    detailed = _get_differential_detail_table(app, meta).copy()
    detailed["population"] = detailed["population"].astype(str)
    subset = detailed.loc[detailed["population"] == population].copy()
    if subset.empty:
        raise HTTPException(status_code=404, detail=f"No heatmap rows were found for '{population}'.")

    subset["gene"] = subset["gene"].astype(str)
    subset["log2fc"] = pd.to_numeric(subset.get("log2fc"), errors="coerce")
    subset["fdr"] = pd.to_numeric(subset.get("fdr"), errors="coerce")
    subset["pval"] = pd.to_numeric(subset.get("pval"), errors="coerce")
    subset["sig_metric"] = subset["fdr"].fillna(subset["pval"])
    subset = subset.dropna(subset=["gene", "log2fc"])
    subset = subset.sort_values(
        ["log2fc", "sig_metric", "pval", "gene"],
        ascending=[False, True, True, True],
    ).drop_duplicates(subset=["gene"])

    detailed_matrix = (
        detailed.assign(
            gene=detailed["gene"].astype(str),
            log2fc=pd.to_numeric(detailed.get("log2fc"), errors="coerce"),
        )
        .dropna(subset=["gene", "log2fc"])
        .pivot_table(index="gene", columns="population", values="log2fc", aggfunc="first")
    )
    fold_rows, column_labels, fold_source = _differential_fold_rows(app, meta)

    if not fold_rows:
        column_labels = list(dict.fromkeys(detailed["population"].astype(str).tolist()))

    detailed_matrix = detailed_matrix.reindex(columns=column_labels)

    rows = []
    rows_from_fold = 0
    rows_from_detail = 0
    for row in subset.itertuples():
        gene = str(row.gene)
        values = fold_rows.get(gene)
        if values is None:
            rows_from_detail += 1
            values = (
                detailed_matrix.loc[gene].tolist()
                if gene in detailed_matrix.index
                else [None] * len(column_labels)
            )
        else:
            rows_from_fold += 1
        row_values = [float(value) if _is_finite_number(value) else None for value in values]
        rows.append(
            {
                "row_key": f"{population}:{gene}",
                "gene": gene,
                "direction": "up" if float(row.log2fc) >= 0 else "down",
                "values": row_values,
            }
        )

    return {
        "population": population,
        "columns": column_labels,
        "rows": rows,
        "default_gene": rows[0]["gene"],
        "fold_source": fold_source,
        "rows_with_measured_folds": rows_from_fold,
        "rows_without_measured_folds": rows_from_detail,
    }


def _build_differential_summary_payload(app: FastAPI, meta: Dict) -> Dict:
    """How many features the differential called up and down in each cell state.

    The detailed DEG table holds one row per feature that passed the run's own fold and
    significance filters, so these counts are the run's own calls and nothing is
    re-thresholded here: log2fc > 0 counts up, log2fc < 0 counts down.

    Cell states follow the fold matrix's column order, which is the lineage order the
    differential used. Every tested state is reported, including the states where the
    run called nothing, so the chart never presents a subset as the whole. A state in
    the detailed table but absent from the fold matrix is appended in the table's order.
    """
    detailed = _get_differential_detail_table(app, meta).copy()
    detailed["population"] = detailed["population"].astype(str)
    detailed["log2fc"] = pd.to_numeric(detailed.get("log2fc"), errors="coerce")
    scored = detailed.dropna(subset=["log2fc"])

    _, fold_columns, fold_source = _differential_fold_rows(app, meta)
    ordered: List[str] = [str(label) for label in fold_columns if str(label)]
    for value in scored["population"].tolist():
        if value and value not in ordered:
            ordered.append(value)

    up_counts = scored.loc[scored["log2fc"] > 0, "population"].value_counts()
    down_counts = scored.loc[scored["log2fc"] < 0, "population"].value_counts()
    sizes = (
        detailed.drop_duplicates(subset=["population"]).set_index("population")
        if "population" in detailed.columns
        else pd.DataFrame()
    )

    rows: List[Dict[str, object]] = []
    for population in ordered:
        up = int(up_counts.get(population, 0))
        down = int(down_counts.get(population, 0))
        entry: Dict[str, object] = {
            "population": population,
            "up": up,
            "down": down,
            "total": up + down,
        }
        if population in getattr(sizes, "index", []):
            for key, column in (("n_case", "n_case"), ("n_control", "n_control")):
                value = sizes.loc[population].get(column)
                entry[key] = int(value) if _is_finite_number(value) else None
        rows.append(entry)

    case_label, control_label = _differential_group_labels(meta)
    modality = _normalize_modality_id(
        (meta.get("differential", {}).get("config", {}) or {}).get("modality"), default="rna"
    )
    feature_label = str(_modality_definition(meta, modality).get("feature_label") or "gene")
    max_count = max([0] + [max(row["up"], row["down"]) for row in rows])
    return {
        "rows": rows,
        "max_count": int(max_count),
        "case_label": case_label,
        "control_label": control_label,
        "feature_label": feature_label,
        "population_col": str((meta.get("differential", {}).get("config", {}) or {}).get("population_col") or ""),
        "n_states_tested": len(rows),
        "n_states_with_calls": int(sum(1 for row in rows if row["total"] > 0)),
        "n_features_called": int(len(scored)),
        "state_order_source": fold_source,
        "default_gene": None,
    }


def _build_differential_volcano_payload(app: FastAPI, meta: Dict, population: str) -> Dict:
    detailed = _get_differential_detail_table(app, meta)
    subset = detailed.loc[detailed["population"] == population].copy()
    if subset.empty:
        raise HTTPException(status_code=404, detail=f"No differential genes were found for '{population}'.")
    subset["log2fc"] = pd.to_numeric(subset.get("log2fc"), errors="coerce")
    subset["fdr"] = pd.to_numeric(subset.get("fdr"), errors="coerce")
    subset["pval"] = pd.to_numeric(subset.get("pval"), errors="coerce")
    subset = subset.dropna(subset=["gene", "log2fc"])
    valid_fdr = np.isfinite(subset["fdr"]) & subset["fdr"].between(0, 1)
    valid_pval = np.isfinite(subset["pval"]) & subset["pval"].between(0, 1)
    # A run with no test has no statistic to put on the y axis. Cell communication with
    # fewer than 2 samples a side is the case: it reports an effect and no p-value.
    # Dropping every point left the panel saying "No volcano data were found", which
    # reads as a broken view rather than an untested comparison. The effect still plots.
    untested = not valid_fdr.any() and not valid_pval.any()
    reason = None
    if untested:
        effect = pd.to_numeric(subset.get("abs_delta_score"), errors="coerce")
        if not np.isfinite(effect).any():
            effect = subset["log2fc"].abs()
        subset = subset.loc[np.isfinite(effect)].copy()
        subset["score"] = effect.loc[subset.index].astype(float)
        statistic = "effect"
        statistic_label = "effect size"
        n_missing_statistic = 0
        reason = str(subset.get("no_test_reason").iloc[0]) if (
            "no_test_reason" in subset.columns and not subset.empty
            and str(subset["no_test_reason"].iloc[0]).strip()) else (
            "No statistical test was run for this comparison, so the y axis shows "
            "effect size rather than significance.")
    else:
        statistic = "fdr" if valid_fdr.any() else "pval"
        statistic_label = "FDR" if statistic == "fdr" else "p-value"
        valid = np.isfinite(subset[statistic]) & subset[statistic].between(0, 1)
        n_missing_statistic = int((~valid).sum())
        subset = subset.loc[valid].copy()
        subset["fdr_clamped"] = subset[statistic].clip(lower=1e-300)
        subset["score"] = -np.log10(subset["fdr_clamped"])
    subset["direction"] = np.where(subset["log2fc"] >= 0, "up", "down")
    if untested:
        subset = subset.sort_values(["score", "gene"], ascending=[False, True]).reset_index(drop=True)
    else:
        subset = subset.sort_values(["fdr_clamped", "pval", "gene"], ascending=[True, True, True]).reset_index(drop=True)
    points = [
        {
            "gene": str(row.gene),
            "log2fc": float(row.log2fc),
            "score": float(row.score),
            "fdr": float(row.fdr) if _is_finite_number(row.fdr) else None,
            "pval": float(row.pval) if _is_finite_number(row.pval) else None,
            "direction": str(row.direction),
        }
        for row in subset.itertuples()
        if _is_finite_number(row.log2fc) and _is_finite_number(row.score)
    ]
    return {
        "population": population,
        "statistic": statistic, "statistic_label": statistic_label,
        "n_missing_statistic": n_missing_statistic,
        "untested": bool(untested),
        "reason": reason,
        "points": points,
        "default_gene": points[0]["gene"] if points else None,
    }


def _build_differential_go_payload(app: FastAPI, meta: Dict, population: str) -> Dict:
    # GO terms describe GENES. An antibody panel, a lipid species and a regulon edge have
    # no gene-set annotation, so only the RNA contrasts run with --goelite_species and
    # only they carry a table. An empty panel with no explanation reads as a broken view,
    # so say which modality is showing and why it has no terms.
    modality = _normalize_modality_id(
        ((meta.get("differential") or {}).get("config") or {}).get("modality"), default="rna")
    empty_reason = None
    if modality != "rna":
        label = str(_modality_definition(meta, modality).get("label") or modality)
        empty_reason = ("GO enrichment describes genes, so it is not computed for "
                        "%s. Switch the modality to RNA to see enriched terms." % label)

    def _empty():
        payload = {"population": population, "terms": [], "default_gene": None}
        if empty_reason:
            payload["message"] = empty_reason
        return payload

    try:
        frame = _get_differential_go_table(app, meta)
    except HTTPException:
        return _empty()
    except FileNotFoundError:
        return _empty()

    if frame.empty or "population" not in frame.columns:
        return _empty()

    directions = frame["population"].map(_split_population_direction)
    frame = frame.assign(
        base_population=[value[0] for value in directions],
        direction=[value[1] or "unknown" for value in directions],
    )
    subset = frame.loc[frame["base_population"] == population].copy()
    if subset.empty:
        return {"population": population, "terms": [], "default_gene": None}

    subset["fdr"] = pd.to_numeric(subset.get("fdr"), errors="coerce")
    subset["p_value"] = pd.to_numeric(subset.get("p_value"), errors="coerce")
    subset["z_score"] = pd.to_numeric(subset.get("z_score"), errors="coerce")
    subset["fdr_plot"] = subset["fdr"].fillna(subset["p_value"]).clip(lower=1e-300, upper=1.0)
    subset["score"] = -np.log10(subset["fdr_plot"])
    subset["is_positive_sig"] = (
        subset["fdr_plot"].le(0.05)
        & subset["z_score"].gt(2.0)
    )
    subset = subset.sort_values(["fdr_plot", "p_value", "z_score", "term_name"], ascending=[True, True, False, True]).reset_index(drop=True)
    terms = []
    for row in subset.itertuples():
        overlap_genes = [gene.strip() for gene in str(getattr(row, "overlap_genes", "")).split(",") if gene.strip()]
        terms.append(
            {
                "term_id": str(getattr(row, "term_id", "")),
                "term_name": str(getattr(row, "term_name", "")),
                "direction": str(getattr(row, "direction", "unknown")),
                "fdr": float(row.fdr) if _is_finite_number(row.fdr) else None,
                "p_value": float(row.p_value) if _is_finite_number(row.p_value) else None,
                "z_score": float(row.z_score) if _is_finite_number(row.z_score) else None,
                "score": float(row.score) if _is_finite_number(row.score) else None,
                "fdr_plot": float(row.fdr_plot) if _is_finite_number(row.fdr_plot) else None,
                "selected": bool(getattr(row, "selected", False)),
                "is_positive_sig": bool(getattr(row, "is_positive_sig", False)),
                "is_selected_positive_sig": bool(getattr(row, "selected", False) and getattr(row, "is_positive_sig", False)),
                "overlap_genes": overlap_genes,
                "selected_gene": overlap_genes[0] if overlap_genes else None,
            }
        )
    positive_sig = [term for term in terms if term.get("is_selected_positive_sig")]
    positive_sig.sort(key=lambda term: (
        float(term.get("fdr_plot", 1.0) or 1.0),
        float(term.get("p_value", 1.0) or 1.0),
        -float(term.get("z_score", 0.0) or 0.0),
        str(term.get("term_name", "")),
    ))
    top_labels = positive_sig[:4]
    used_names = {str(term.get("term_name", "")).strip().lower() for term in top_labels}
    keyword_term = None
    for term in positive_sig[4:]:
        term_name_lc = str(term.get("term_name", "")).strip().lower()
        if term_name_lc in used_names:
            continue
        if any(keyword in term_name_lc for keyword in _GO_ELITE_HIGHLIGHT_KEYWORDS):
            keyword_term = term
            break
    if keyword_term is None:
        for term in terms:
            term_name_lc = str(term.get("term_name", "")).strip().lower()
            if term_name_lc in used_names:
                continue
            if float(term.get("z_score", 0.0) or 0.0) <= 0.0:
                continue
            if any(keyword in term_name_lc for keyword in _GO_ELITE_HIGHLIGHT_KEYWORDS):
                keyword_term = term
                break
    labels = []
    for index, term in enumerate(top_labels):
        labels.append(
            {
                "term_name": str(term.get("term_name", "")),
                "z_score": float(term.get("z_score", 0.0) or 0.0),
                "fdr_plot": float(term.get("fdr_plot", 1.0) or 1.0),
                "selected_gene": term.get("selected_gene"),
                "overlap_genes": term.get("overlap_genes") or [],
                "label_color": "#111827",
                "label_rank": index,
                "label_role": "top",
            }
        )
    if keyword_term is not None:
        labels.append(
            {
                "term_name": str(keyword_term.get("term_name", "")),
                "z_score": float(keyword_term.get("z_score", 0.0) or 0.0),
                "fdr_plot": float(keyword_term.get("fdr_plot", 1.0) or 1.0),
                "selected_gene": keyword_term.get("selected_gene"),
                "overlap_genes": keyword_term.get("overlap_genes") or [],
                "label_color": "#4f6ef7",
                "label_rank": len(labels),
                "label_role": "keyword",
            }
        )
    default_gene = next((term["selected_gene"] for term in terms if term["selected_gene"]), None)
    return {"population": population, "terms": terms, "labels": labels, "default_gene": default_gene}


def _communication_focus_positions(states: List[str], focus: str) -> Dict[str, Dict[str, float]]:
    ordered_states = [str(state).strip() for state in states if str(state).strip()]
    if not ordered_states:
        return {}
    focus_state = str(focus).strip()
    if focus_state and focus_state in ordered_states:
        ordered_states = [focus_state] + [state for state in ordered_states if state != focus_state]
    if len(ordered_states) == 1:
        return {ordered_states[0]: {"x": 0.0, "y": 0.0}}
    positions: Dict[str, Dict[str, float]] = {}
    if focus_state and focus_state in ordered_states:
        positions[focus_state] = {"x": 0.0, "y": 0.0}
        ring_states = [state for state in ordered_states if state != focus_state]
    else:
        ring_states = ordered_states
    if not ring_states:
        return positions
    radius = 260.0
    for index, state in enumerate(ring_states):
        angle = (2.0 * np.pi * float(index)) / float(max(len(ring_states), 1))
        positions[state] = {
            "x": float(np.cos(angle) * radius),
            "y": float(np.sin(angle) * radius),
        }
    return positions


def _build_differential_cell_communication_network_payload(app: FastAPI, meta: Dict, population: str) -> Dict:
    detailed = _get_differential_detail_table(app, meta).copy()
    detailed["population"] = detailed["population"].astype(str)
    subset = detailed.loc[detailed["population"] == population].copy()
    if subset.empty:
        return {"population": population, "network_type": "cell_communication_diff", "elements": [], "default_gene": None}

    for column in ("delta_score", "case_mean_score", "control_mean_score", "fdr", "pval", "lr_expression_score", "receiver_response_score"):
        subset[column] = pd.to_numeric(subset.get(column), errors="coerce")
    subset["abs_delta_score"] = pd.to_numeric(subset.get("abs_delta_score"), errors="coerce").fillna(subset["delta_score"].abs())
    subset = subset.sort_values(["abs_delta_score", "fdr", "pval"], ascending=[False, True, True]).reset_index(drop=True)

    pair_rows: List[Dict[str, object]] = []
    for (sender, receiver), pair_df in subset.groupby(["sender_state", "receiver_state"], sort=False):
        ranked = pair_df.sort_values(["abs_delta_score", "fdr", "pval"], ascending=[False, True, True]).reset_index(drop=True)
        top_rows = ranked.head(8)
        top_interactions = []
        tooltip_lines = [
            f"{sender} -> {receiver}",
            f"{ranked.shape[0]} differential ligand-receptor interaction(s)",
        ]
        for interaction in top_rows.itertuples(index=False):
            interaction_gene = str(getattr(interaction, "gene", "") or "").strip()
            interaction_label = str(getattr(interaction, "interaction", "") or "").strip()
            ligand_symbol = str(getattr(interaction, "ligand", "") or "").strip()
            receptor_symbol = str(getattr(interaction, "receptor", "") or "").strip()
            delta_score = float(getattr(interaction, "delta_score", 0.0) or 0.0)
            case_mean = float(getattr(interaction, "case_mean_score", 0.0) or 0.0)
            control_mean = float(getattr(interaction, "control_mean_score", 0.0) or 0.0)
            lr_score = float(getattr(interaction, "lr_expression_score", 0.0) or 0.0)
            response_score = float(getattr(interaction, "receiver_response_score", 0.0) or 0.0)
            top_interactions.append(
                {
                    "gene": interaction_gene,
                    "interaction": interaction_label,
                    "ligand": ligand_symbol,
                    "receptor": receptor_symbol,
                    "delta_score": delta_score,
                    "case_mean_score": case_mean,
                    "control_mean_score": control_mean,
                    "lr_expression_score": lr_score,
                    "receiver_response_score": response_score,
                    "fdr": float(getattr(interaction, "fdr", np.nan)) if _is_finite_number(getattr(interaction, "fdr", np.nan)) else None,
                    "response_support_genes": str(getattr(interaction, "response_support_genes", "") or ""),
                }
            )
            tooltip_lines.append(
                f"{interaction_label}: "
                f"delta={delta_score:.3f}; case={case_mean:.3f}; control={control_mean:.3f}; "
                f"LR={lr_score:.3f}; response={response_score:.3f}"
            )
        total_delta = float(ranked["delta_score"].sum())
        pair_rows.append(
            {
                "sender": str(sender).strip(),
                "receiver": str(receiver).strip(),
                "top_score": float(ranked["abs_delta_score"].max()),
                "total_score": abs(total_delta),
                "signed_delta": total_delta,
                "case_mean_score": float(ranked["case_mean_score"].sum()),
                "control_mean_score": float(ranked["control_mean_score"].sum()),
                "n_interactions": int(ranked.shape[0]),
                "top_interactions": top_interactions,
                "tooltip": "\n".join(tooltip_lines),
            }
        )

    pair_rows = _annotate_fastcomm_pair_visual_weights(pair_rows)
    states = sorted({str(population).strip(), *subset["sender_state"].astype(str).tolist(), *subset["receiver_state"].astype(str).tolist()})
    positions = _communication_focus_positions(states, population)
    nodes: Dict[str, Dict[str, object]] = {}
    for state in states:
        node_id = f"state::{state}"
        is_focus = state == str(population).strip()
        nodes[node_id] = {
            "data": {
                "id": node_id,
                "label": state,
                "node_type": "focus" if is_focus else "state",
                "color": "#0f766e" if is_focus else "#e0f2fe",
                "node_score": 1.0 if is_focus else 0.0,
            },
            "position": positions.get(state, {"x": 0.0, "y": 0.0}),
        }

    edges: List[Dict[str, object]] = []
    for index, pair in enumerate(pair_rows):
        top_interactions = list(pair.get("top_interactions") or [])
        top_interaction = top_interactions[0] if top_interactions else {}
        signed_delta = float(pair.get("signed_delta", 0.0) or 0.0)
        edge_color = "#dc2626" if signed_delta > 0 else ("#2563eb" if signed_delta < 0 else "#64748b")
        label = str(top_interaction.get("interaction", "") or f"{int(pair.get('n_interactions', 0) or 0)} LR")
        edges.append(
            {
                "data": {
                    "id": f"ccd{index}",
                    "source": f"state::{str(pair.get('sender', '')).strip()}",
                    "target": f"state::{str(pair.get('receiver', '')).strip()}",
                    "interaction_type": "cell_communication_diff",
                    "direction": "positive" if signed_delta >= 0 else "negative",
                    "score": float(pair.get("top_score", 0.0) or 0.0),
                    "delta_score": signed_delta,
                    "abs_delta_score": float(pair.get("total_score", 0.0) or 0.0),
                    "case_mean_score": float(pair.get("case_mean_score", 0.0) or 0.0),
                    "control_mean_score": float(pair.get("control_mean_score", 0.0) or 0.0),
                    "n_interactions": int(pair.get("n_interactions", 0) or 0),
                    "weight": float(pair.get("weight", 2.0) or 2.0),
                    "edge_opacity": float(pair.get("edge_opacity", 0.65) or 0.65),
                    "edge_color": edge_color,
                    "label": label,
                    "gene": str(top_interaction.get("gene", "") or ""),
                    "tooltip": (
                        str(pair.get("tooltip", ""))
                        + f"\nSummed delta={signed_delta:.3f}; visual weight={float(pair.get('weight', 2.0) or 2.0):.1f}"
                    ),
                    "top_interactions": top_interactions,
                }
            }
        )

    default_gene = str(subset.iloc[0].get("gene", "") or "").strip() or None
    return {
        "population": population,
        "network_type": "cell_communication_diff",
        "elements": list(nodes.values()) + edges,
        "default_gene": default_gene,
    }


def _build_differential_cell_communication_table_payload(app: FastAPI, meta: Dict, population: str, limit: int = 40) -> Dict:
    detailed = _get_differential_detail_table(app, meta).copy()
    detailed["population"] = detailed["population"].astype(str)
    subset = detailed.loc[detailed["population"] == population].copy()
    if subset.empty:
        return {"population": population, "columns": [], "rows": [], "default_gene": None}
    for column in ("delta_score", "case_mean_score", "control_mean_score", "log2fc", "fdr", "pval", "lr_expression_score", "receiver_response_score"):
        subset[column] = pd.to_numeric(subset.get(column), errors="coerce")
    subset["abs_delta_score"] = pd.to_numeric(subset.get("abs_delta_score"), errors="coerce").fillna(subset["delta_score"].abs())
    subset = subset.sort_values(["fdr", "abs_delta_score", "max_sample_score"], ascending=[True, False, False]).head(max(1, int(limit))).copy()
    columns = [
        "sender_state",
        "receiver_state",
        "ligand",
        "receptor",
        "delta_score",
        "case_mean_score",
        "control_mean_score",
        "log2fc",
        "fdr",
        "pval",
        "lr_expression_score",
        "receiver_response_score",
        "response_support_genes",
    ]
    available_columns = [column for column in columns if column in subset.columns]
    rows = []
    for record in subset.where(pd.notnull(subset), None).to_dict("records"):
        cleaned = {}
        for key, value in record.items():
            if isinstance(value, (np.floating, float)):
                cleaned[key] = float(value) if _is_finite_number(value) else None
            elif isinstance(value, (np.integer, int)):
                cleaned[key] = int(value)
            else:
                cleaned[key] = value
        rows.append(cleaned)
    default_gene = str(subset.iloc[0].get("gene", "") or "").strip() or None
    return {
        "population": population,
        "plot_type": "cell_communication_diff_table",
        "columns": available_columns,
        "rows": rows,
        "default_gene": default_gene,
    }


def _build_differential_network_payload(app: FastAPI, meta: Dict, population: str, root_path: str = "") -> Dict:
    modality = _normalize_modality_id((meta.get("differential", {}) or {}).get("config", {}).get("modality"), default="rna")
    if modality == "cell_communication":
        return _build_differential_cell_communication_network_payload(app, meta, population)

    entry = _differential_network_entry(meta, population)
    if not entry:
        return {"population": population, "elements": [], "default_gene": None}

    tsv_path = Path(str(entry.get("tsv", "")).strip())
    if not tsv_path.exists():
        return {"population": population, "elements": [], "default_gene": None}

    detailed = _get_differential_detail_table(app, meta)
    subset = detailed.loc[detailed["population"] == population, ["gene", "log2fc", "fdr", "pval"]].copy()
    subset["gene"] = subset["gene"].astype(str)
    subset["log2fc"] = pd.to_numeric(subset.get("log2fc"), errors="coerce")
    fc_map = dict(zip(subset["gene"], subset["log2fc"]))

    cache_entry = _get_differential_cache_entry(app, meta)
    network_tables = cache_entry["network_tables"]
    frame = network_tables.get(population)
    if not isinstance(frame, pd.DataFrame):
        frame = pd.read_csv(tsv_path, sep="\t")
        network_tables[population] = frame
    nodes: Dict[str, Dict[str, object]] = {}
    edges = []
    edge_index = 0
    for row in frame.itertuples():
        source = str(getattr(row, "Symbol1", "")).strip()
        targets = [value.strip() for value in str(getattr(row, "Symbol2", "")).split("|") if value.strip()]
        interaction_type = str(getattr(row, "InteractionType", "")).strip()
        direction = str(getattr(row, "Direction", "")).strip().lower() or "neutral"
        if not source or not targets:
            continue
        if source not in nodes:
            source_fc = fc_map.get(source)
            nodes[source] = {
                "data": {
                    "id": source,
                    "label": source,
                    "log2fc": float(source_fc) if _is_finite_number(source_fc) else None,
                }
            }
        for target in targets:
            if target not in nodes:
                target_fc = fc_map.get(target)
                nodes[target] = {
                    "data": {
                        "id": target,
                        "label": target,
                        "log2fc": float(target_fc) if _is_finite_number(target_fc) else None,
                    }
                }
            edges.append(
                {
                    "data": {
                        "id": f"e{edge_index}",
                        "source": source,
                        "target": target,
                        "interaction_type": interaction_type,
                        "direction": direction,
                    }
                }
            )
            edge_index += 1

    ranked_nodes = sorted(
        (
            (
                abs(float(element["data"].get("log2fc", 0.0) or 0.0)),
                str(element["data"].get("id", "")),
            )
            for element in nodes.values()
        ),
        reverse=True,
    )
    default_gene = ranked_nodes[0][1] if ranked_nodes else None
    return {
        "population": population,
        "elements": list(nodes.values()) + edges,
        "default_gene": default_gene,
        "pdf_url": _with_root_path(root_path, f"/api/jobs/{meta['job_id']}/differential/network/{entry['id']}?format=pdf"),
    }


_GRN_EDGE_CACHE: Dict[str, Any] = {}


def _grn_edges_adata(meta: Dict) -> ad.AnnData:
    # `network_h5ad` is the GRN NETWORK's own key. It must not be `differential_h5ad`,
    # because `_has_edge_level_h5ad` (app.py:1193) treats that key as "this run has a
    # real edge-level h5ad" and routes the DETAIL VIOLIN to it. For a precomputed
    # bundle that object is the donor-level pseudobulk, so the violin would draw one
    # point per donor. The COPD atlas may not publish per-donor points, and the violin
    # must stay on metacells. A cellHarmony job still sets differential_h5ad, so that
    # remains the fallback.
    grn_meta = (meta.get("modality_artifacts") or {}).get("grn") or {}
    edges_path = str(grn_meta.get("network_h5ad")
                     or grn_meta.get("differential_h5ad", "")).strip()
    if not edges_path or not Path(edges_path).exists():
        raise FileNotFoundError("GRN edge output is unavailable for this job.")
    mtime = Path(edges_path).stat().st_mtime
    cached = _GRN_EDGE_CACHE.get(edges_path)
    if not cached or cached[0] != mtime:
        _GRN_EDGE_CACHE[edges_path] = (mtime, ad.read_h5ad(edges_path))
    return _GRN_EDGE_CACHE[edges_path][1]


_LINEAGE_ORDER_CACHE: Dict[str, Any] = {}


def _as_str_list(value):   # uns['lineage_order'] may be a numpy array
    return [] if value is None else [str(v) for v in list(value)]


def _combined_lineage_order(meta: Dict) -> List[str]:
    """Cached lineage_order from the combined RNA h5ad (read once per job, not per request)."""
    cpath = str((meta.get("artifacts") or {}).get("combined_h5ad", "")).strip()
    if not cpath or not Path(cpath).exists():
        return []
    mtime = Path(cpath).stat().st_mtime
    cached = _LINEAGE_ORDER_CACHE.get(cpath)
    if not cached or cached[0] != mtime:
        try:
            order = _as_str_list(ad.read_h5ad(cpath, backed="r").uns.get("lineage_order"))
        except Exception:
            order = []
        _LINEAGE_ORDER_CACHE[cpath] = (mtime, order)
    return _LINEAGE_ORDER_CACHE[cpath][1]


def _order_by_lineage(values, adata, meta: Dict) -> List[str]:
    """Order cell-state values by uns['lineage_order'] (from the edge h5ad, or cached from the
    combined RNA h5ad for older jobs); unknowns sorted at the end."""
    vals = {str(v) for v in values}
    order = _as_str_list(adata.uns.get("lineage_order")) or _combined_lineage_order(meta)
    ordered = [c for c in order if c in vals]
    ordered += sorted(v for v in vals if v not in set(ordered))
    return ordered


def _build_grn_network_payload(meta: Dict, genes, *, sample: str = "", cell_state: str = "",
                               threshold: float = 0.0, max_edges: int = 25) -> Dict:
    """Cytoscape network of rna2grn edges (TF->target) involving the requested genes, scored
    for a selected sample/cell-state. A blank state selects the first in lineage order;
    blank genes select the TF with the largest summed absolute outgoing score in that state.
    Scores average the selected sample x state pseudobulks, never unrelated cell states."""
    adata = _grn_edges_adata(meta)
    var_names = adata.var_names.astype(str).to_numpy()
    cluster_key = str(meta.get("cluster_key") or meta.get("reference_cluster_key") or "")
    # value lists for the explorer dropdowns (condition-style sample column preferred)
    sample_pref = ("Library", "library", "group", "condition", "Sample", "sample", "dataset", "Dataset", "orig.ident")
    disp_col = next((c for c in sample_pref
                     if c in adata.obs.columns and 1 < int(adata.obs[c].astype(str).nunique()) <= 200), None)
    available_samples = sorted(set(adata.obs[disp_col].astype(str))) if disp_col else []
    cell_values = (set(adata.obs[cluster_key].astype(str)) if cluster_key and cluster_key in adata.obs.columns else set())
    available_cell_states = _order_by_lineage(cell_values, adata, meta)
    cell_state = cell_state or (available_cell_states[0] if available_cell_states else '')
    req = {_normalize_gene_token(g) for g in (genes or []) if str(g).strip()}
    empty = {"sample": sample, "cell_state": cell_state, "threshold": threshold, "elements": [], "n_edges": 0,
             "score_type": "mean predicted edge score", "sample_field": disp_col,
             "n_matching_edges": 0, "n_aggregates": 0, "max_edges": max_edges,
             "available_samples": available_samples, "available_cell_states": available_cell_states}
    if not available_cell_states:
        return {**empty, "message": "The GRN output has no cell-state annotations."}
    if cell_state not in available_cell_states:
        return {**empty, "message": f"Cell type '{cell_state}' not found."}

    mask = np.ones(adata.n_obs, dtype=bool)
    if sample:
        # match the requested value against whichever sample-like column actually holds it
        # (e.g. the user picks "AML" from Library, or "0" from sample)
        scol = next((c for c in dict.fromkeys([disp_col, *pipeline_mod._PB_SAMPLE_COLS]) if c
                     if c in adata.obs.columns and str(sample) in set(adata.obs[c].astype(str))), None)
        if scol is None:
            return {**empty, "message": f"Sample '{sample}' not found."}
        mask &= (adata.obs[scol].astype(str).to_numpy() == str(sample))
    if cell_state and cluster_key and cluster_key in adata.obs.columns:
        mask &= (adata.obs[cluster_key].astype(str).to_numpy() == str(cell_state))
    if int(mask.sum()) == 0:
        return {**empty, "message": "No pseudobulks match the selected sample / cell type."}

    sub = adata.X[mask]
    # Average predicted scores across the selected aggregates, without a contrast
    # or marker-enrichment calculation. Preserve their original scale and precision.
    scores = np.asarray(sub.mean(axis=0), dtype=float).ravel()
    factors = {}
    for key, score in zip(var_names, scores):
        if not np.isfinite(score):
            continue
        for tf in key.split('|')[:-1]:
            factors[tf] = factors.get(tf, 0.0) + abs(float(score))
    ranked_factors = sorted(factors, key=lambda tf: (-factors[tf], tf))
    selected_genes = [str(g) for g in (genes or []) if str(g).strip()]
    if not req and ranked_factors:
        selected_genes = [ranked_factors[0]]
        req = {_normalize_gene_token(ranked_factors[0])}
    empty.update(genes=selected_genes, available_factors=ranked_factors, n_aggregates=int(mask.sum()))

    # keep edges whose TF or target matches a requested gene, above the |score| threshold
    keep = []
    for j in range(len(var_names)):
        ev = var_names[j]
        if "|" not in ev:
            continue
        *tfs, tgt = ev.split("|")
        sc = float(scores[j])
        if not np.isfinite(sc) or sc == 0 or abs(sc) < float(threshold):
            continue
        for tf in tfs:
            if tf and tgt and (_normalize_gene_token(tf) in req or _normalize_gene_token(tgt) in req):
                keep.append((j, tf, tgt, sc))
    keep.sort(key=lambda t: -abs(t[3]))
    n_matching_edges = len(keep)
    keep = keep[:int(max_edges)]

    nodes: Dict[str, Dict] = {}
    edges_out = []
    for j, tf, tgt, sc in keep:
        for nid, role in ((tf, "tf"), (tgt, "target")):
            node = nodes.setdefault(nid, {"data": {"id": nid, "label": nid, "role": role,
                                                    "queried": _normalize_gene_token(nid) in req}})
            if role == "tf":
                node["data"]["role"] = "tf"
        edges_out.append({"data": {
            "id": f"{tf}__{tgt}__{j}", "source": tf, "target": tgt, "score": sc,
            "interaction_type": "transcription", "direction": "neutral"}})
    msg = "" if edges_out else "No GRN edges matched the requested genes at this threshold."
    return {**empty, "sample": sample, "cell_state": cell_state, "threshold": threshold,
            "n_matching_edges": n_matching_edges,
            "n_edges": len(edges_out), "elements": list(nodes.values()) + edges_out, "message": msg,
            "genes": selected_genes, "available_factors": ranked_factors, "node_encoding": "role",
            "available_samples": available_samples, "available_cell_states": available_cell_states}


def _build_marker_network_payload(meta: Dict, population: str, modality: str = "rna") -> Dict:
    marker_analysis = _modality_marker_analysis(meta, modality) or {}
    networks = marker_analysis.get("networks") or []
    entry = next((item for item in networks if str(item.get("population", "")).strip() == str(population).strip()), None)
    if not entry:
        return {"population": population, "elements": [], "default_gene": None}

    tsv_path = Path(str(entry.get("tsv", "")).strip())
    if not tsv_path.exists():
        return {"population": population, "elements": [], "default_gene": None}

    marker_stats_path = Path(str(marker_analysis.get("markers_tsv", "")).strip())
    fc_map: Dict[str, float] = {}
    # Networks use redundant per-state markers, not only the unique top marker list.
    redundant_path = Path(str(marker_analysis.get("redundant_markers_tsv") or
                              marker_stats_path.with_name(marker_stats_path.stem.replace("_markers", "_redundant_markers") + ".tsv")))
    paths = [path for path in (marker_stats_path, redundant_path) if path.is_file()]
    if paths:
        stats_df = pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)
        if "cluster" in stats_df.columns:
            subset = stats_df.loc[stats_df["cluster"].astype(str) == str(population)].copy()
        else:
            subset = stats_df.copy()
        if "Gene" in subset.columns and "Fold" in subset.columns:
            subset["Gene"] = subset["Gene"].astype(str)
            subset["Fold"] = pd.to_numeric(subset.get("Fold"), errors="coerce")
            fc_map = {
                str(gene): float(fold)
                for gene, fold in zip(subset["Gene"], subset["Fold"])
                if _is_finite_number(fold)
            }

    frame = pd.read_csv(tsv_path, sep="\t")
    nodes: Dict[str, Dict[str, object]] = {}
    edges = []
    edge_index = 0
    for row in frame.itertuples():
        source = str(getattr(row, "Symbol1", "")).strip()
        targets = [value.strip() for value in str(getattr(row, "Symbol2", "")).split("|") if value.strip()]
        interaction_type = str(getattr(row, "InteractionType", "")).strip()
        direction = str(getattr(row, "Direction", "")).strip().lower() or "neutral"
        if not source or not targets:
            continue
        if source not in nodes:
            source_fc = fc_map.get(source)
            nodes[source] = {"data": {"id": source, "label": source, "log2fc": source_fc}}
        for target in targets:
            if target not in nodes:
                target_fc = fc_map.get(target)
                nodes[target] = {"data": {"id": target, "label": target, "log2fc": target_fc}}
            edges.append(
                {
                    "data": {
                        "id": f"m{edge_index}",
                        "source": source,
                        "target": target,
                        "interaction_type": interaction_type,
                        "direction": direction,
                    }
                }
            )
            edge_index += 1

    ranked_nodes = sorted(
        ((abs(float(element["data"].get("log2fc", 0.0) or 0.0)), str(element["data"].get("id", ""))) for element in nodes.values()),
        reverse=True,
    )
    default_gene = ranked_nodes[0][1] if ranked_nodes else None
    return {
        "population": population,
        "elements": list(nodes.values()) + edges,
        "default_gene": default_gene,
    }


def _fastcomm_analysis(meta: Dict) -> Dict[str, object]:
    analysis = meta.get("fastcomm_analysis") or {}
    if not isinstance(analysis, dict):
        return {}
    return dict(analysis)


def _fastcomm_scores_path(meta: Dict) -> Path:
    analysis = _fastcomm_analysis(meta)
    raw_path = str(analysis.get("scores_tsv") or meta.get("artifacts", {}).get("fastcomm_scores") or "").strip()
    if not raw_path:
        raise FileNotFoundError("fastComm scores are unavailable.")
    path = Path(raw_path)
    if not path.exists():
        raise FileNotFoundError("fastComm scores are unavailable.")
    return path


def _fastcomm_scores_table(app: FastAPI, meta: Dict) -> pd.DataFrame:
    job_id = str(meta.get("job_id") or "").strip()
    if not job_id:
        raise ValueError("Job metadata is missing a job_id.")
    if not hasattr(app.state, "fastcomm_cache"):
        app.state.fastcomm_cache = {}
    cache = app.state.fastcomm_cache
    path = _fastcomm_scores_path(meta)
    signature = (str(path), path.stat().st_mtime_ns, path.stat().st_size)
    entry = cache.get(job_id)
    if isinstance(entry, dict) and entry.get("signature") == signature and isinstance(entry.get("scores"), pd.DataFrame):
        return entry["scores"]
    scores = pd.read_csv(path, sep="\t")
    cache[job_id] = {"signature": signature, "scores": scores}
    return scores


def _fastcomm_split_scores_path(meta: Dict) -> Optional[Path]:
    analysis = _fastcomm_analysis(meta)
    per_sample = analysis.get("per_sample") or {}
    raw_path = str(per_sample.get("split_scores_long_tsv") or "").strip() if isinstance(per_sample, dict) else ""
    if not raw_path and isinstance(per_sample, dict):
        output_dir = str(per_sample.get("output_dir") or "").strip()
        if output_dir:
            raw_path = str(Path(output_dir) / "split_scores_long.tsv")
    if not raw_path:
        scores_path = str(analysis.get("scores_tsv") or "").strip()
        if scores_path:
            raw_path = str(Path(scores_path).parent / "per_sample" / "split_scores_long.tsv")
    if not raw_path:
        return None
    path = Path(raw_path)
    return path if path.exists() else None


def _fastcomm_split_scores_table(app: FastAPI, meta: Dict) -> pd.DataFrame:
    job_id = str(meta.get("job_id") or "").strip()
    path = _fastcomm_split_scores_path(meta)
    if path is None:
        raise FileNotFoundError("Per-sample fastComm scores are unavailable.")
    if not hasattr(app.state, "fastcomm_cache"):
        app.state.fastcomm_cache = {}
    cache = app.state.fastcomm_cache
    signature = (str(path), path.stat().st_mtime_ns, path.stat().st_size)
    cache_key = f"{job_id}:split_scores"
    entry = cache.get(cache_key)
    if isinstance(entry, dict) and entry.get("signature") == signature and isinstance(entry.get("scores"), pd.DataFrame):
        return entry["scores"]
    scores = pd.read_csv(path, sep="\t")
    cache[cache_key] = {"signature": signature, "scores": scores}
    return scores


def _prepare_fastcomm_scores(scores: pd.DataFrame) -> pd.DataFrame:
    required = {"sender_state", "receiver_state", "ligand", "receptor", "fastcomm_score"}
    missing = required.difference(scores.columns)
    if missing:
        raise ValueError(f"fastComm scores missing required columns: {', '.join(sorted(missing))}")
    frame = scores.copy()
    for column in ("sender_state", "receiver_state", "ligand", "receptor"):
        frame[column] = frame[column].astype(str)
    for column in ("response_key", "response_support_genes", "pathway"):
        if column in frame.columns:
            frame[column] = frame[column].astype(str).replace({"nan": "", "None": ""})
    for column in ("fastcomm_score", "receiver_response_score", "lr_expression_score_scaled", "lr_expression_score"):
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce").fillna(0.0)
    if "lr_expression_score_scaled" not in frame.columns and "lr_expression_score" in frame.columns:
        frame["lr_expression_score_scaled"] = frame["lr_expression_score"]
    if "receiver_response_score" not in frame.columns:
        frame["receiver_response_score"] = 0.0
    if "lr_expression_score_scaled" not in frame.columns:
        frame["lr_expression_score_scaled"] = 0.0
    return frame


def _fastcomm_filter_focus(scores: pd.DataFrame, focus: str, direction: str) -> pd.DataFrame:
    direction_key = str(direction or "incoming").strip().lower()
    if direction_key == "outgoing":
        return scores.loc[scores["sender_state"].astype(str) == focus].copy()
    return scores.loc[scores["receiver_state"].astype(str) == focus].copy()


def _first_nonempty_value(values: pd.Series) -> str:
    for value in values:
        text = str(value or "").strip()
        if text:
            return text
    return ""


def _fastcomm_selected_splits(
    app: FastAPI,
    meta: Dict,
    display_filters: Optional[List[tuple[str, List[str]]]],
) -> Optional[List[str]]:
    if not display_filters:
        return None
    analysis = _fastcomm_analysis(meta)
    per_sample = analysis.get("per_sample") or {}
    split_key = str(per_sample.get("split_key") or "").strip() if isinstance(per_sample, dict) else ""
    sample_key = str(analysis.get("sample_key") or split_key or "").strip()
    if not sample_key:
        return None
    selected = _display_filter_values(display_filters, sample_key)
    if not selected:
        return None
    return sorted(dict.fromkeys(str(value).strip() for value in selected if str(value).strip()))


def _display_filter_values(display_filters: Optional[List[tuple[str, List[str]]]], field: str) -> List[str]:
    target = str(field or "").strip()
    if not target or not display_filters:
        return []
    values: List[str] = []
    for filter_field, filter_values in display_filters:
        if str(filter_field or "").strip() != target:
            continue
        values.extend(str(value).strip() for value in filter_values if str(value).strip())
    return list(dict.fromkeys(values))


def _fastcomm_state_filter_values(meta: Dict, display_filters: Optional[List[tuple[str, List[str]]]]) -> List[str]:
    analysis = _fastcomm_analysis(meta)
    state_key = str(analysis.get("state_key") or meta.get("cluster_key") or "").strip()
    return _display_filter_values(display_filters, state_key)


def _filter_fastcomm_scores_by_states(scores: pd.DataFrame, state_values: Sequence[str]) -> pd.DataFrame:
    selected = {str(value).strip() for value in state_values if str(value).strip()}
    if not selected or scores.empty:
        return scores
    return scores.loc[
        scores["sender_state"].astype(str).isin(selected)
        | scores["receiver_state"].astype(str).isin(selected)
    ].copy()


def _aggregate_fastcomm_scores(scores: pd.DataFrame) -> pd.DataFrame:
    if scores.empty:
        return scores.copy()
    group_cols = ["sender_state", "receiver_state", "ligand", "receptor"]
    frame = scores.copy()
    frame["response_key"] = frame.get("response_key", pd.Series("", index=frame.index)).astype(str)
    frame["response_support_genes"] = frame.get("response_support_genes", pd.Series("", index=frame.index)).astype(str)
    aggregations = {
        "fastcomm_score": "mean",
        "receiver_response_score": "mean",
        "lr_expression_score_scaled": "mean",
        "response_key": _first_nonempty_value,
        "response_support_genes": _first_nonempty_value,
    }
    if "empirical_percentile" in frame.columns:
        aggregations["empirical_percentile"] = "mean"
    if "empirical_rank" in frame.columns:
        aggregations["empirical_rank"] = "mean"
    aggregated = frame.groupby(group_cols, as_index=False).agg(aggregations)
    return aggregated.sort_values("fastcomm_score", ascending=False).reset_index(drop=True)


def _fastcomm_populations(app: FastAPI, meta: Dict) -> List[str]:
    try:
        scores = _fastcomm_scores_table(app, meta)
    except Exception:
        return []
    values = pd.concat(
        [
            scores.get("sender_state", pd.Series(dtype=str)),
            scores.get("receiver_state", pd.Series(dtype=str)),
        ],
        ignore_index=True,
    )
    return sorted({str(value).strip() for value in values if str(value).strip()})


def _fastcomm_interaction_payload(interaction: pd.Series) -> Dict[str, object]:
    return {
        "sender_state": str(interaction.get("sender_state", "")).strip(),
        "receiver_state": str(interaction.get("receiver_state", "")).strip(),
        "ligand": str(interaction.get("ligand", "")).strip(),
        "receptor": str(interaction.get("receptor", "")).strip(),
        "score": float(interaction.get("fastcomm_score", 0.0) or 0.0),
        "receiver_response_score": float(interaction.get("receiver_response_score", 0.0) or 0.0),
        "lr_expression_score": float(interaction.get("lr_expression_score_scaled", 0.0) or 0.0),
        "pathway": str(interaction.get("response_key", "") or interaction.get("pathway", "") or "").strip(),
        "supporting_response_genes": str(interaction.get("response_support_genes", "") or "").strip(),
    }


def _fastcomm_pair_records(scores: pd.DataFrame, *, limit: int = 35) -> List[Dict[str, object]]:
    if scores.empty:
        return []
    frame = scores.sort_values("fastcomm_score", ascending=False)
    records: List[Dict[str, object]] = []
    for (sender, receiver), pair_df in frame.groupby(["sender_state", "receiver_state"], sort=False):
        ranked = pair_df.sort_values("fastcomm_score", ascending=False)
        top_interactions = [_fastcomm_interaction_payload(interaction) for _, interaction in ranked.head(8).iterrows()]
        tooltip_lines = [
            f"{sender} -> {receiver}",
            f"{ranked.shape[0]} significant ligand-receptor interaction(s)",
        ]
        for item in top_interactions:
            tooltip_lines.append(
                f"{item['ligand']}->{item['receptor']}: "
                f"score={item['score']:.3f}; LR={item['lr_expression_score']:.3f}; response={item['receiver_response_score']:.3f}"
            )
        records.append(
            {
                "sender": str(sender).strip(),
                "receiver": str(receiver).strip(),
                "top_score": float(ranked["fastcomm_score"].max()),
                "total_score": float(ranked["fastcomm_score"].sum()),
                "n_interactions": int(ranked.shape[0]),
                "top_interactions": top_interactions,
                "tooltip": "\n".join(tooltip_lines),
            }
        )
    return sorted(records, key=lambda item: float(item["total_score"]), reverse=True)[: max(1, int(limit))]


def _annotate_fastcomm_pair_visual_weights(pair_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    if not pair_rows:
        return []
    totals = [float(row.get("total_score", 0.0) or 0.0) for row in pair_rows]
    min_total = min(totals)
    max_total = max(totals)
    denom = max_total - min_total
    count = len(pair_rows)
    weighted_rows: List[Dict[str, object]] = []
    for index, row in enumerate(pair_rows):
        total_score = float(row.get("total_score", 0.0) or 0.0)
        score_norm = ((total_score - min_total) / denom) if denom > 1e-12 else 0.0
        rank_norm = 1.0 if count <= 1 else 1.0 - (index / float(count - 1))
        visual_norm = score_norm if denom > 1e-12 else rank_norm
        if denom > 1e-12:
            visual_norm = 0.8 * score_norm + 0.2 * rank_norm
        weight = 2.0 + 12.0 * visual_norm
        edge_opacity = 0.3 + 0.65 * visual_norm
        weighted_rows.append(
            {
                **row,
                "visual_norm": float(visual_norm),
                "weight": float(weight),
                "edge_opacity": float(edge_opacity),
            }
        )
    return weighted_rows


def _build_fastcomm_network_elements(pair_rows: List[Dict[str, object]], *, focus: str = "") -> List[Dict[str, object]]:
    pair_rows = _annotate_fastcomm_pair_visual_weights(pair_rows)
    nodes: Dict[str, Dict[str, object]] = {}
    if focus:
        focus_id = f"state::{focus}"
        nodes[focus_id] = {"data": {"id": focus_id, "label": focus, "node_type": "focus", "color": "#0f766e"}}
    edges: List[Dict[str, object]] = []
    for index, pair in enumerate(pair_rows):
        sender = str(pair["sender"]).strip()
        receiver = str(pair["receiver"]).strip()
        for state in (sender, receiver):
            node_id = f"state::{state}"
            if node_id not in nodes:
                nodes[node_id] = {
                    "data": {
                        "id": node_id,
                        "label": state,
                        "node_type": "state",
                        "color": "#e0f2fe",
                    }
                }
        if focus and f"state::{focus}" in nodes:
            nodes[f"state::{focus}"]["data"]["node_type"] = "focus"
            nodes[f"state::{focus}"]["data"]["color"] = "#0f766e"
        top_interactions = list(pair.get("top_interactions") or [])
        top_interaction = top_interactions[0] if top_interactions else {}
        ligand = str(top_interaction.get("ligand", "")).strip()
        receptor = str(top_interaction.get("receptor", "")).strip()
        total_score = float(pair.get("total_score", 0.0) or 0.0)
        edges.append(
            {
                "data": {
                    "id": f"fc{index}",
                    "source": f"state::{sender}",
                    "target": f"state::{receiver}",
                    "interaction_type": "fastcomm",
                    "direction": "positive",
                    "ligand": ligand,
                    "receptor": receptor,
                    "pathway": str(top_interaction.get("pathway", "")).strip(),
                    "score": float(pair.get("top_score", 0.0) or 0.0),
                    "total_score": total_score,
                    "n_interactions": int(pair.get("n_interactions", 0) or 0),
                    "weight": float(pair.get("weight", 2.0) or 2.0),
                    "edge_opacity": float(pair.get("edge_opacity", 0.65) or 0.65),
                    "visual_norm": float(pair.get("visual_norm", 0.0) or 0.0),
                    "label": f"{ligand}->{receptor}" if ligand and receptor else f"{pair.get('n_interactions', 0)} LR",
                    "tooltip": (
                        str(pair.get("tooltip", ""))
                        + f"\nTotal pair score={total_score:.3f}; visual weight={float(pair.get('weight', 2.0) or 2.0):.1f}"
                    ),
                    "top_interactions": top_interactions,
                    "supporting_response_genes": str(top_interaction.get("supporting_response_genes", "") or ""),
                    "receiver_response_score": float(top_interaction.get("receiver_response_score", 0.0) or 0.0),
                    "lr_expression_score": float(top_interaction.get("lr_expression_score", 0.0) or 0.0),
                }
            }
        )
    return list(nodes.values()) + edges


def _build_fastcomm_payload(app: FastAPI, meta: Dict, population: str, direction: str = "incoming", limit: int = 35) -> Dict:
    scores = _prepare_fastcomm_scores(_fastcomm_scores_table(app, meta))

    focus = str(population or "").strip()
    if not focus:
        populations = _fastcomm_populations(app, meta)
        focus = populations[0] if populations else ""
    if not focus:
        return {"population": "", "direction": direction, "elements": [], "message": "No fastComm populations are available."}

    direction_key = str(direction or "incoming").strip().lower()
    if direction_key not in {"incoming", "outgoing"}:
        direction_key = "incoming"

    subset = _fastcomm_filter_focus(scores, focus, direction_key)

    if subset.empty:
        return {
            "population": focus,
            "direction": direction_key,
            "elements": [],
            "message": f"No fastComm {direction_key} interactions were available for {focus}.",
        }

    pair_rows = _fastcomm_pair_records(subset, limit=limit)
    elements = _build_fastcomm_network_elements(pair_rows, focus=focus)

    return {
        "population": focus,
        "direction": direction_key,
        "plot_type": "focused_outgoing" if direction_key == "outgoing" else "focused_incoming",
        "state_key": str(_fastcomm_analysis(meta).get("state_key") or meta.get("cluster_key") or ""),
        "elements": elements,
        "summary": {
            "n_edges": int(subset.shape[0]),
            "n_displayed_edges": int(sum(1 for item in elements if "source" in item.get("data", {}))),
            "top_score": max((float(item["top_score"]) for item in pair_rows), default=0.0),
        },
    }


def _build_fastcomm_focus_payload_from_scores(
    scores: pd.DataFrame,
    *,
    meta: Dict,
    population: str,
    direction: str = "incoming",
    limit: int = 35,
    selected_splits: Optional[List[str]] = None,
) -> Dict:
    focus = str(population or "").strip()
    if not focus:
        values = pd.concat(
            [
                scores.get("sender_state", pd.Series(dtype=str)),
                scores.get("receiver_state", pd.Series(dtype=str)),
            ],
            ignore_index=True,
        )
        populations = sorted({str(value).strip() for value in values if str(value).strip()})
        focus = populations[0] if populations else ""
    if not focus:
        return {"population": "", "direction": direction, "elements": [], "message": "No fastComm populations are available."}

    direction_key = str(direction or "incoming").strip().lower()
    if direction_key not in {"incoming", "outgoing"}:
        direction_key = "incoming"
    subset = _fastcomm_filter_focus(scores, focus, direction_key)
    if subset.empty:
        return {
            "population": focus,
            "direction": direction_key,
            "plot_type": "focused_outgoing" if direction_key == "outgoing" else "focused_incoming",
            "elements": [],
            "message": f"No fastComm {direction_key} interactions were available for {focus}.",
            "selected_splits": selected_splits or [],
        }
    pair_rows = _fastcomm_pair_records(subset, limit=limit)
    elements = _build_fastcomm_network_elements(pair_rows, focus=focus)
    return {
        "population": focus,
        "direction": direction_key,
        "plot_type": "focused_outgoing" if direction_key == "outgoing" else "focused_incoming",
        "state_key": str(_fastcomm_analysis(meta).get("state_key") or meta.get("cluster_key") or ""),
        "selected_splits": selected_splits or [],
        "elements": elements,
        "summary": {
            "n_edges": int(subset.shape[0]),
            "n_displayed_edges": int(sum(1 for item in elements if "source" in item.get("data", {}))),
            "top_score": max((float(item["top_score"]) for item in pair_rows), default=0.0),
        },
    }


def _build_fastcomm_plot_payload(
    app: FastAPI,
    meta: Dict,
    population: str,
    plot_type: str = "focused_incoming",
    limit: int = 60,
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
) -> Dict:
    scores = _prepare_fastcomm_scores(_fastcomm_scores_table(app, meta))
    state_filter_values = _fastcomm_state_filter_values(meta, display_filters)
    selected_splits = _fastcomm_selected_splits(app, meta, display_filters)
    if selected_splits:
        try:
            split_scores = _prepare_fastcomm_scores(_fastcomm_split_scores_table(app, meta))
        except FileNotFoundError:
            split_scores = pd.DataFrame()
        if "split" not in split_scores.columns:
            selected_splits = []
        else:
            split_scores["split"] = split_scores["split"].astype(str)
            filtered_split_scores = split_scores.loc[split_scores["split"].isin(selected_splits)].copy()
            if filtered_split_scores.empty:
                raise HTTPException(status_code=404, detail="No per-sample cell-communication scores match the current display filters.")
            scores = _aggregate_fastcomm_scores(filtered_split_scores)
    if state_filter_values:
        scores = _filter_fastcomm_scores_by_states(scores, state_filter_values)
    plot_key = str(plot_type or "focused_incoming").strip().lower()
    aliases = {
        "incoming": "focused_incoming",
        "outgoing": "focused_outgoing",
        "network": "cell_state_network",
        "global_network": "cell_state_network",
        "heatmap": "state_heatmap",
        "dotplot": "lr_dotplot",
        "table": "top_table",
        "sample": "per_sample",
    }
    plot_key = aliases.get(plot_key, plot_key)
    if plot_key not in {
        "focused_incoming",
        "focused_outgoing",
        "cell_state_network",
        "lr_dotplot",
        "state_heatmap",
        "top_table",
        "per_sample",
    }:
        plot_key = "focused_incoming"

    focus = str(population or "").strip()
    if state_filter_values and (not focus or focus not in set(state_filter_values)):
        focus = state_filter_values[0]
    if not focus:
        populations = _fastcomm_populations(app, meta)
        focus = populations[0] if populations else ""

    if plot_key in {"focused_incoming", "focused_outgoing"}:
        direction = "outgoing" if plot_key == "focused_outgoing" else "incoming"
        return _build_fastcomm_focus_payload_from_scores(
            scores,
            meta=meta,
            population=focus,
            direction=direction,
            limit=min(limit, 60),
            selected_splits=selected_splits,
        )

    if plot_key == "cell_state_network":
        pair_rows = _fastcomm_pair_records(scores, limit=limit)
        elements = _build_fastcomm_network_elements(pair_rows)
        return {
            "plot_type": plot_key,
            "population": focus,
            "direction": "global",
            "state_key": str(_fastcomm_analysis(meta).get("state_key") or meta.get("cluster_key") or ""),
            "elements": elements,
            "selected_splits": selected_splits or [],
            "summary": {
                "n_edges": int(scores.shape[0]),
                "n_displayed_edges": int(sum(1 for item in elements if "source" in item.get("data", {}))),
                "top_score": max((float(item["top_score"]) for item in pair_rows), default=0.0),
            },
        }

    if plot_key == "state_heatmap":
        grouped = (
            scores.groupby(["sender_state", "receiver_state"], as_index=False)
            .agg(total_score=("fastcomm_score", "sum"), max_score=("fastcomm_score", "max"), n_interactions=("fastcomm_score", "size"))
            .sort_values("total_score", ascending=False)
        )
        senders = sorted(grouped["sender_state"].astype(str).unique().tolist())
        receivers = sorted(grouped["receiver_state"].astype(str).unique().tolist())
        matrix = grouped.pivot_table(index="sender_state", columns="receiver_state", values="total_score", aggfunc="sum", fill_value=0.0)
        matrix = matrix.reindex(index=senders, columns=receivers, fill_value=0.0)
        hover = grouped.set_index(["sender_state", "receiver_state"]).to_dict("index")
        text = []
        for sender in senders:
            row_text = []
            for receiver in receivers:
                item = hover.get((sender, receiver), {})
                row_text.append(
                    f"{sender} -> {receiver}<br>"
                    f"total score={float(item.get('total_score', 0.0) or 0.0):.3f}<br>"
                    f"max score={float(item.get('max_score', 0.0) or 0.0):.3f}<br>"
                    f"interactions={int(item.get('n_interactions', 0) or 0)}"
                )
            text.append(row_text)
        return {
            "plot_type": plot_key,
            "population": focus,
            "selected_splits": selected_splits or [],
            "senders": senders,
            "receivers": receivers,
            "z": matrix.to_numpy(dtype=float).tolist(),
            "text": text,
            "summary": {"n_edges": int(scores.shape[0]), "n_state_pairs": int(grouped.shape[0])},
        }

    direction = "outgoing" if plot_key == "focused_outgoing" else "incoming"
    focused = _fastcomm_filter_focus(scores, focus, direction) if focus else scores.copy()
    if focused.empty:
        return {
            "plot_type": plot_key,
            "population": focus,
            "direction": direction,
            "rows": [],
            "points": [],
            "message": f"No fastComm interactions were available for {focus}.",
        }
    focused = focused.sort_values("fastcomm_score", ascending=False)

    if plot_key == "lr_dotplot":
        other_col = "receiver_state" if direction == "outgoing" else "sender_state"
        top = focused.head(max(1, int(limit))).copy()
        points = []
        for _, row in top.iterrows():
            interaction = f"{row['ligand']}->{row['receptor']}"
            other = str(row.get(other_col, "")).strip()
            points.append(
                {
                    "x": other,
                    "y": interaction,
                    "score": float(row.get("fastcomm_score", 0.0) or 0.0),
                    "receiver_response_score": float(row.get("receiver_response_score", 0.0) or 0.0),
                    "lr_expression_score": float(row.get("lr_expression_score_scaled", 0.0) or 0.0),
                    "sender_state": str(row.get("sender_state", "")).strip(),
                    "receiver_state": str(row.get("receiver_state", "")).strip(),
                    "supporting_response_genes": str(row.get("response_support_genes", "") or ""),
                }
            )
        return {
            "plot_type": plot_key,
            "population": focus,
            "direction": direction,
            "selected_splits": selected_splits or [],
            "points": points,
            "summary": {"n_edges": int(focused.shape[0]), "n_displayed": int(len(points))},
        }

    if plot_key == "top_table":
        columns = [
            "sender_state",
            "receiver_state",
            "ligand",
            "receptor",
            "fastcomm_score",
            "lr_expression_score_scaled",
            "receiver_response_score",
            "response_key",
            "response_support_genes",
            "empirical_percentile",
            "empirical_rank",
        ]
        table_source = _fastcomm_filter_focus(scores, focus, "incoming") if focus else scores.copy()
        if table_source.empty:
            table_source = scores.copy()
        table_source = table_source.sort_values("fastcomm_score", ascending=False)
        available_columns = [column for column in columns if column in table_source.columns]
        table = table_source.loc[:, available_columns].head(max(1, int(limit))).copy()
        return {
            "plot_type": plot_key,
            "population": focus,
            "direction": "incoming",
            "selected_splits": selected_splits or [],
            "columns": available_columns,
            "rows": table.to_dict("records"),
            "summary": {"n_edges": int(table_source.shape[0]), "n_displayed": int(table.shape[0])},
        }

    if plot_key == "per_sample":
        try:
            split_scores = _prepare_fastcomm_scores(_fastcomm_split_scores_table(app, meta))
        except FileNotFoundError:
            return {
                "plot_type": plot_key,
                "population": focus,
                "direction": direction,
                "sample_key": str((_fastcomm_analysis(meta).get("sample_key") or "sample")),
                "selected_splits": selected_splits or [],
                "rows": [],
                "message": "Per-sample fastComm scores are unavailable. Re-run the job with the current reduced fastComm settings.",
            }
        if "split" not in split_scores.columns:
            raise ValueError("Per-sample fastComm scores missing required column: split")
        split_scores["split"] = split_scores["split"].astype(str)
        if selected_splits:
            split_scores = split_scores.loc[split_scores["split"].isin(selected_splits)].copy()
            if split_scores.empty:
                raise HTTPException(status_code=404, detail="No per-sample cell-communication scores match the current display filters.")
        if state_filter_values:
            split_scores = _filter_fastcomm_scores_by_states(split_scores, state_filter_values)
        split_focused = _fastcomm_filter_focus(split_scores, focus, direction) if focus else split_scores.copy()
        other_col = "receiver_state" if direction == "outgoing" else "sender_state"
        grouped = (
            split_focused.groupby(["split", other_col], as_index=False)
            .agg(total_score=("fastcomm_score", "sum"), max_score=("fastcomm_score", "max"), n_interactions=("fastcomm_score", "size"))
            .sort_values("total_score", ascending=False)
        )
        top_states = (
            grouped.groupby(other_col, as_index=False)["total_score"].sum().sort_values("total_score", ascending=False).head(12)[other_col].astype(str).tolist()
        )
        grouped = grouped.loc[grouped[other_col].astype(str).isin(top_states)].copy()
        rows = [
            {
                "sample": str(row.split),
                "state": str(getattr(row, other_col)),
                "total_score": float(row.total_score),
                "max_score": float(row.max_score),
                "n_interactions": int(row.n_interactions),
            }
            for row in grouped.itertuples(index=False)
        ]
        return {
            "plot_type": plot_key,
            "population": focus,
            "direction": direction,
            "sample_key": str((_fastcomm_analysis(meta).get("sample_key") or "sample")),
            "selected_splits": selected_splits or [],
            "rows": rows,
            "summary": {"n_rows": int(len(rows)), "n_samples": int(grouped["split"].nunique()) if not grouped.empty else 0},
        }

    return _build_fastcomm_payload(app, meta, focus, direction="incoming", limit=min(limit, 60))


def _build_cell_communication_feature_expression(
    *,
    app: FastAPI,
    meta: Dict,
    feature_symbol: str,
    feature_role: str,
    sender_state: str,
    receiver_state: str,
    case_label: str,
    control_label: str,
    group1_samples: List[str],
    group2_samples: List[str],
    sample_field: str,
    population_col: str,
) -> Optional[Dict]:
    expression_cache = _get_expression_cache(app, meta, modality="rna")
    adata = expression_cache["adata"]
    resolved_gene = _resolve_gene_name(expression_cache["var_names"], feature_symbol)
    if not resolved_gene:
        raise KeyError(f"Gene '{feature_symbol}' not found in the aligned AnnData output.")
    if population_col not in adata.obs.columns:
        raise ValueError(f"'{population_col}' is not present in the aligned AnnData observations.")

    if sample_field and sample_field in adata.obs.columns:
        sample_col = sample_field
        available_values = set(adata.obs[sample_col].astype(str))
        missing = sorted((set(group1_samples) | set(group2_samples)) - available_values)
        if missing:
            raise ValueError(
                f"Selected group values were not found in obs['{sample_col}']: {', '.join(missing)}"
            )
        resolved_group1 = list(group1_samples)
        resolved_group2 = list(group2_samples)
    else:
        sample_col, resolved_samples = pipeline_mod._resolve_samples_for_adata(
            adata, meta, group1_samples + group2_samples
        )
        resolved_group1 = [resolved_samples[sample] for sample in group1_samples]
        resolved_group2 = [resolved_samples[sample] for sample in group2_samples]

    population_values = adata.obs[population_col].astype(str)
    sample_values = adata.obs[sample_col].astype(str)

    role = (feature_role or "").strip().lower()
    if role == "ligand":
        target_state = sender_state
    elif role == "receptor":
        target_state = receiver_state
    else:
        target_state = receiver_state or sender_state

    target_state = str(target_state or "").strip()
    if not target_state:
        return None

    mask = (population_values == target_state) & sample_values.isin(resolved_group1 + resolved_group2)
    if int(np.asarray(mask).sum()) == 0:
        return None

    subset = adata[mask.to_numpy(), resolved_gene]
    values = _flatten_expr(subset.X).astype(float)
    subset_samples = sample_values.loc[mask].astype(str)
    group_labels = np.where(subset_samples.isin(resolved_group1), case_label, control_label)

    groups = []
    for label in (case_label, control_label):
        group_values = values[group_labels == label]
        finite_values = group_values[np.isfinite(group_values)]
        groups.append(
            {
                "label": label,
                "values": [float(value) for value in finite_values],
                "n_cells": int(len(finite_values)),
                "mean": float(np.mean(finite_values)) if len(finite_values) else 0.0,
            }
        )

    return {
        "gene": resolved_gene,
        "feature_symbol": resolved_gene,
        "feature_role": role,
        "feature_state": target_state,
        "view_kind": "feature_expression",
        "groups": groups,
    }


def _build_differential_gene_detail_payload(
    app: FastAPI,
    meta: Dict,
    population: str,
    gene: str,
    feature: Optional[str] = None,
) -> Dict:
    differential = meta.get("differential", {}) or {}
    config = differential.get("config", {}) or {}
    modality = _normalize_modality_id(config.get("modality"), default="rna")
    modality_info = _modality_definition(meta, modality)
    group1_samples = [str(value).strip() for value in config.get("group1_samples", []) if str(value).strip()]
    group2_samples = [str(value).strip() for value in config.get("group2_samples", []) if str(value).strip()]
    sample_field = str(config.get("sample_field", "")).strip()
    if not group1_samples or not group2_samples:
        # Some contrasts cannot be grouped from metacell annotation at all: the six
        # tertile contrasts (FEV1, FVC, FEV1/FVC, DLCO, pack-years, weight) split donors
        # on a continuous covariate that obs does not carry, and current_vs_never needs a
        # "current" smoking value no metacell holds. The statistics are still valid, so
        # say that rather than answering HTTP 500.
        raise HTTPException(
            status_code=404,
            detail=("This comparison groups samples on a covariate the released "
                    "data do not carry, so the per-replicate distribution cannot "
                    "be drawn. The differential statistics above are unaffected."))
    population_col = _differential_population_col(meta)
    case_label, control_label = _differential_group_labels(meta)

    if modality == "cell_communication":
        detailed = _get_differential_detail_table(app, meta)
        population_rows = detailed.loc[detailed["population"] == population].copy()
        if population_rows.empty:
            raise KeyError(f"No differential interactions were found for '{population}'.")
        stats = population_rows.loc[population_rows["gene"] == gene].copy()
        if stats.empty:
            population_rows["abs_delta_score"] = pd.to_numeric(
                population_rows.get("abs_delta_score"), errors="coerce"
            ).fillna(pd.to_numeric(population_rows.get("delta_score"), errors="coerce").abs())
            stats = population_rows.sort_values(
                ["abs_delta_score", "fdr", "pval"], ascending=[False, True, True]
            ).head(1).copy()
        row = stats.sort_values(["fdr", "pval"], ascending=[True, True]).iloc[0]
        ligand_symbol = str(row.get("ligand", "") or "").strip()
        receptor_symbol = str(row.get("receptor", "") or "").strip()
        sender_state = str(row.get("sender_state", "") or "").strip()
        receiver_state = str(row.get("receiver_state", "") or "").strip() or population

        feature_request = str(feature or "").strip()
        feature_role = ""
        feature_symbol = ""
        if feature_request:
            if feature_request.lower() == ligand_symbol.lower() and ligand_symbol:
                feature_role = "ligand"
                feature_symbol = ligand_symbol
            elif feature_request.lower() == receptor_symbol.lower() and receptor_symbol:
                feature_role = "receptor"
                feature_symbol = receptor_symbol
            else:
                feature_symbol = feature_request

        per_cell_payload = None
        if feature_symbol:
            try:
                per_cell_payload = _build_cell_communication_feature_expression(
                    app=app,
                    meta=meta,
                    feature_symbol=feature_symbol,
                    feature_role=feature_role,
                    sender_state=sender_state,
                    receiver_state=receiver_state,
                    case_label=case_label,
                    control_label=control_label,
                    group1_samples=group1_samples,
                    group2_samples=group2_samples,
                    sample_field=sample_field,
                    population_col=population_col,
                )
            except KeyError:
                per_cell_payload = None
            except (FileNotFoundError, ValueError):
                per_cell_payload = None

        if per_cell_payload is not None:
            per_cell_payload["population"] = population
            per_cell_payload["interaction"] = str(row.get("interaction") or gene)
            per_cell_payload["interaction_key"] = gene
            per_cell_payload["ligand"] = ligand_symbol
            per_cell_payload["receptor"] = receptor_symbol
            per_cell_payload["sender_state"] = sender_state
            per_cell_payload["receiver_state"] = receiver_state
            per_cell_payload["modality"] = modality
            per_cell_payload["feature_label"] = "expression"
            per_cell_payload["stats"] = {
                "p_value": float(row["pval"]) if _is_finite_number(row.get("pval")) else None,
                "fdr": float(row["fdr"]) if _is_finite_number(row.get("fdr")) else None,
                "log2fc": float(row["log2fc"]) if _is_finite_number(row.get("log2fc")) else None,
                "n_case": int(row.get("n_case", 0) or 0),
                "n_control": int(row.get("n_control", 0) or 0),
            }
            return per_cell_payload

        case_mean = float(row.get("case_mean_score", 0.0) or 0.0)
        control_mean = float(row.get("control_mean_score", 0.0) or 0.0)
        return {
            "population": population,
            "gene": str(row.get("interaction") or gene),
            "interaction": str(row.get("interaction") or gene),
            "interaction_key": gene,
            "ligand": ligand_symbol,
            "receptor": receptor_symbol,
            "sender_state": sender_state,
            "receiver_state": receiver_state,
            "modality": modality,
            "feature_label": str(modality_info.get("feature_label") or "ligand-receptor interaction"),
            "view_kind": "interaction_summary",
            "groups": [
                {"label": case_label, "values": [case_mean], "n_cells": int(row.get("n_case", 0) or 0), "mean": case_mean},
                {"label": control_label, "values": [control_mean], "n_cells": int(row.get("n_control", 0) or 0), "mean": control_mean},
            ],
            "stats": {
                "p_value": float(row["pval"]) if _is_finite_number(row.get("pval")) else None,
                "fdr": float(row["fdr"]) if _is_finite_number(row.get("fdr")) else None,
                "log2fc": float(row["log2fc"]) if _is_finite_number(row.get("log2fc")) else None,
                "n_case": int(row.get("n_case", 0) or 0),
                "n_control": int(row.get("n_control", 0) or 0),
            },
        }

    if modality == "grn" and _has_edge_level_h5ad(meta):
        # GRN detail is an edge ("TF|target") distribution — read the edge-level h5ad, not the
        # per-TF activity matrix served by the expression cache.
        adata, _grn_detail_path = _open_gene_detail_adata(app, meta, gene)
        var_names = adata.var_names.astype(str).to_numpy()
    else:
        expression_cache = _get_expression_cache(app, meta, modality=modality)
        adata = expression_cache["adata"]
        var_names = expression_cache["var_names"]
    resolved_gene = _resolve_gene_name(var_names, gene)
    ragged_only = False
    if not resolved_gene and modality == "grn":
        # The per-cell-state GRN sidecar and the per-metacell ragged store hold DIFFERENT
        # edge sets. The sidecar scores every edge in every state; the ragged store keeps
        # only the edges significant in that state, which is what holds it near 0.6 GB.
        # A GRN differential therefore names edges the sidecar lacks. Measured on the COPD
        # bundle: `KLF6|NTN4` is 0 of 63,647 sidecar edges and 1 of the 4,035 edges in
        # `imputed_v7_grn_ragged/AT1.h5ad`, so clicking that volcano point answered
        # "Gene 'KLF6|NTN4' not found in the aligned AnnData output" while its p-value,
        # FDR and fold change were on screen beside it.
        # The ragged store supplies the violin in either branch, so resolve against it and
        # take the values from it alone. A name in neither store still raises below.
        if _grn_ragged_values(meta, population, gene, adata.obs.index) is not None:
            resolved_gene, ragged_only = gene, True
    if not resolved_gene:
        raise KeyError(f"Gene '{gene}' not found in the aligned AnnData output.")
    if population_col not in adata.obs.columns:
        raise ValueError(f"'{population_col}' is not present in the aligned AnnData observations.")
    if sample_field and sample_field in adata.obs.columns:
        sample_col = sample_field
        available_values = set(adata.obs[sample_col].astype(str))
        missing = sorted((set(group1_samples) | set(group2_samples)) - available_values)
        if missing:
            raise ValueError(
                f"Selected group values were not found in obs['{sample_col}']: {', '.join(missing)}"
            )
        resolved_group1 = group1_samples
        resolved_group2 = group2_samples
    else:
        sample_col, resolved_samples = pipeline_mod._resolve_samples_for_adata(adata, meta, group1_samples + group2_samples)
        resolved_group1 = [resolved_samples[sample] for sample in group1_samples]
        resolved_group2 = [resolved_samples[sample] for sample in group2_samples]

    population_values = adata.obs[population_col].astype(str)
    sample_values = adata.obs[sample_col].astype(str)
    if str(population).strip() == _POOLED_OVERALL_LABEL:
        mask = sample_values.isin(resolved_group1 + resolved_group2)
    else:
        mask = (population_values == population) & sample_values.isin(resolved_group1 + resolved_group2)
    if int(np.asarray(mask).sum()) == 0:
        # The volcano and the violin read DIFFERENT objects. The statistics come from
        # the sample pseudobulks, which cover every cell state the alignment found.
        # The violin draws the replicate unit this dataset releases, which may cover
        # fewer states: the COPD atlas releases 50 metacell states of 81 tested. A
        # state with results but no replicates is not a failure, so say which of the
        # two is missing rather than reporting the panel broken.
        present = sorted(set(population_values[sample_values.isin(
            resolved_group1 + resolved_group2)].unique()))
        raise HTTPException(
            status_code=404,
            detail=(f"'{population}' carries differential statistics but no replicate "
                    f"profiles in this atlas, so the distribution cannot be drawn. "
                    f"The comparison covers {len(present)} cell states with replicates."))

    if ragged_only:
        # No sidecar column exists for this edge; every value comes from the ragged store.
        values = np.full(int(mask.to_numpy().sum()), np.nan, dtype=float)
    else:
        subset = adata[mask.to_numpy(), resolved_gene]
        values = _flatten_expr(subset.X).astype(float)
    if modality == "grn":
        # Prefer the per-metacell ragged store when it covers this state and edge. The
        # bundle sidecar holds one value per CELL STATE, which draws a flat line.
        # BundleAnnData exposes obs but not obs_names, so the index is read directly.
        ragged = _grn_ragged_values(meta, population, resolved_gene, adata.obs.index)
        if ragged is not None:
            picked = ragged[mask.to_numpy()]
            if np.isfinite(picked).any():
                values = np.where(np.isfinite(picked), picked, values)
    subset_samples = sample_values.loc[mask].astype(str)
    group_labels = np.where(subset_samples.isin(resolved_group1), case_label, control_label)

    groups = []
    for label in (case_label, control_label):
        group_values = values[group_labels == label]
        finite_values = group_values[np.isfinite(group_values)]
        groups.append(
            {
                "label": label,
                "values": [float(value) for value in finite_values],
                "n_cells": int(len(finite_values)),
                "mean": float(np.mean(finite_values)) if len(finite_values) else 0.0,
            }
        )

    detailed = _get_differential_detail_table(app, meta)
    stats = detailed.loc[(detailed["population"] == population) & (detailed["gene"] == resolved_gene)].copy()
    stats_payload = {
        "p_value": None,
        "fdr": None,
        "log2fc": None,
        "n_case": groups[0]["n_cells"],
        "n_control": groups[1]["n_cells"],
    }
    if not stats.empty:
        stats = stats.sort_values(["fdr", "pval"], ascending=[True, True]).iloc[0]
        stats_payload = {
            "p_value": float(stats["pval"]) if _is_finite_number(stats.get("pval")) else None,
            "fdr": float(stats["fdr"]) if _is_finite_number(stats.get("fdr")) else None,
            "log2fc": float(stats["log2fc"]) if _is_finite_number(stats.get("log2fc")) else None,
            "n_case": int(stats["n_case"]) if _is_finite_number(stats.get("n_case")) else groups[0]["n_cells"],
            "n_control": int(stats["n_control"]) if _is_finite_number(stats.get("n_control")) else groups[1]["n_cells"],
        }

    return {
        "population": population,
        "gene": resolved_gene,
        "modality": modality,
        "value_label": ("Predicted TF activity (sum)" if modality == "grn_tf" else
                        "GRN edge score" if modality == "grn" else "Normalized expression"),
        "feature_label": str(modality_info.get("feature_label") or "gene"),
        "groups": groups,
        "stats": stats_payload,
    }


def _validate_differential_request(meta: Dict, payload: DifferentialSettings) -> Dict:
    options = _differential_options(meta)
    if not options["enabled"]:
        raise HTTPException(status_code=400, detail="Differential analysis requires two or more uploaded samples.")
    if meta.get("status") != "completed":
        raise HTTPException(status_code=400, detail="Run the cellHarmony analysis before starting differential analysis.")

    modality = _normalize_modality_id(payload.modality or options.get("default_modality") or "rna")
    modality_options = {_normalize_modality_id(entry.get("id"), default="") for entry in options.get("modalities", []) if isinstance(entry, dict)}
    if modality_options and modality not in modality_options:
        raise HTTPException(status_code=400, detail="Selected modality is not available for this job.")

    population_options = {entry["value"] for entry in options.get("population_columns", []) if entry.get("value")}
    if population_options and payload.population_col not in population_options:
        raise HTTPException(status_code=400, detail="Selected cell-state field is not available for this job.")

    sample_field = str(payload.sample_field or options.get("default_sample_field") or "").strip()
    sample_field_options = {entry["value"] for entry in options.get("sample_fields", []) if entry.get("value")}
    if sample_field_options and sample_field not in sample_field_options:
        raise HTTPException(status_code=400, detail="Selected group obs field is not available for this job.")
    if sample_field and sample_field == payload.population_col:
        raise HTTPException(status_code=400, detail="Group values must come from a different obs field than the selected cell-state field.")

    group1_samples = [str(sample).strip() for sample in payload.group1_samples if str(sample).strip()]
    group2_samples = [str(sample).strip() for sample in payload.group2_samples if str(sample).strip()]
    if not group1_samples or not group2_samples:
        raise HTTPException(status_code=400, detail="Both differential groups must include at least one sample.")

    available_samples = set((options.get("sample_values") or {}).get(sample_field, [])) if sample_field else set(options["sample_names"])
    missing = sorted((set(group1_samples) | set(group2_samples)) - available_samples)
    if missing:
        raise HTTPException(status_code=400, detail=f"Unknown samples selected: {', '.join(missing)}")

    overlap = sorted(set(group1_samples) & set(group2_samples))
    if overlap:
        raise HTTPException(status_code=400, detail=f"Samples cannot appear in both groups: {', '.join(overlap)}")

    comparison_type = str(payload.comparison_type or "cells").strip().lower()
    if comparison_type not in {"cells", "pseudobulk"}:
        raise HTTPException(status_code=400, detail="Comparison Type must be either 'cells' or 'pseudobulk'.")
    upload_profile = dict((meta.get("differential_options") or {}).get("upload_profile") or pipeline_mod._upload_profile(meta))
    pseudobulk_allowed = bool(upload_profile.get("single_h5ad") or upload_profile.get("total_files", 0) >= 4)
    if comparison_type == "pseudobulk" and not pseudobulk_allowed:
        comparison_type = "cells"

    return {
        "modality": modality,
        "population_col": payload.population_col,
        "sample_field": sample_field,
        "group1_samples": group1_samples,
        "group2_samples": group2_samples,
        "comparison_type": comparison_type,
    }


def _get_differential_artifact(meta: Dict, key: str) -> Path:
    raw_path = str(meta.get("differential", {}).get("artifacts", {}).get(key, "")).strip()
    if not raw_path:
        raise HTTPException(status_code=404, detail="Differential artifact unavailable.")
    path = Path(raw_path)
    if not path.exists():
        raise HTTPException(status_code=404, detail="Differential artifact unavailable.")
    return path


def _obsm_embedding_keys(adata) -> List[str]:
    """The obsm entries that can be drawn as a 2-D map, in the h5ad's own order.

    An entry needs at least two columns; the first two are plotted. PCA and
    scVI latent spaces qualify as much as a UMAP does, so nothing is filtered on
    the name - a dataset that stores `X_umap_harmony`, `X_tsne` or `X_scvi` gets
    all of them.
    """
    keys = []
    for key in getattr(adata, "obsm", {}) or {}:
        try:
            matrix = adata.obsm[key]
            if getattr(matrix, "ndim", 0) == 2 and matrix.shape[1] >= 2:
                keys.append(str(key))
        except Exception:  # noqa: BLE001 - an unreadable entry is simply not offered
            continue
    return keys


def _numeric_obs_columns(cache: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The numeric obs columns a panel may plot on an axis.

    ShinyCell lets a reader put any numeric cell annotation on an axis, so a
    stored `UMAP_1`/`UMAP_2` pair, a pseudotime, a module score or a QC measure
    all work. Counts live in X, not in obs, so nothing here is expression. A
    boolean column is a category, not a measurement, and is left out. So is a
    column with no finite value, which would draw an empty panel.
    """
    adata = cache["adata"]
    out = []
    for name in adata.obs.columns:
        series = adata.obs[name]
        if pd.api.types.is_bool_dtype(series) or not pd.api.types.is_numeric_dtype(series):
            continue
        values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(values)
        if not finite.any():
            continue
        # Nathan's rule, 2026-09-01, in his words: "any numerical columns
        # non-counts or non-scaled counts obs field that are 100% float values".
        # A column whose every value is a whole number is a count or an index,
        # not a measurement: n_cells, n_counts, metacell, and the per-covariate
        # `<name>__n_obs` tallies. Plotting one against a measurement draws a
        # meaningless panel, so it is not offered as an axis. The test is on the
        # VALUES, so no field name is hardcoded and a new count is caught too.
        if np.all(values[finite] == np.round(values[finite])):
            continue
        out.append({
            "field": str(name),
            "n_finite": int(finite.sum()),
            "n_missing": int(values.size - finite.sum()),
            "n_unique": int(np.unique(values[finite]).size),
            "min": float(values[finite].min()),
            "max": float(values[finite].max()),
        })
    return out


UMAP_AXIS_FIELDS = ("__umap_1", "__umap_2")


def _axis_values(cache: Dict[str, Any], field: str) -> Optional[np.ndarray]:
    """One numeric obs column as floats, or None when it cannot serve as an axis."""
    column = str(field or "").strip()
    if not column:
        return None
    # The cellHarmony embedding is not an obs column, but a reader picking X and Y
    # must be able to choose it, otherwise the default pair cannot be the UMAP.
    if column == UMAP_AXIS_FIELDS[0]:
        return np.asarray(cache["umap_x"], dtype=float)
    if column == UMAP_AXIS_FIELDS[1]:
        return np.asarray(cache["umap_y"], dtype=float)
    adata = cache["adata"]
    if column not in adata.obs.columns:
        return None
    series = adata.obs[column]
    if pd.api.types.is_bool_dtype(series) or not pd.api.types.is_numeric_dtype(series):
        return None
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    return values if np.isfinite(values).any() else None


def _axis_field_options(cache: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Every field a coordinate axis may use, the cellHarmony UMAP pair first.

    "Coordinates" is two dropdowns, not one. The first two entries are the two
    axes of the embedding the panel has always drawn, so an untouched panel picks
    them by position and draws the UMAP. Everything after is a float measurement.
    """
    x = np.asarray(cache["umap_x"], dtype=float)
    y = np.asarray(cache["umap_y"], dtype=float)
    out = []
    for field, label, values in ((UMAP_AXIS_FIELDS[0], "cellHarmony UMAP 1", x),
                                 (UMAP_AXIS_FIELDS[1], "cellHarmony UMAP 2", y)):
        finite = np.isfinite(values)
        out.append({"field": field, "label": label,
                    "n_finite": int(finite.sum()),
                    "n_missing": int(values.size - finite.sum()),
                    "n_unique": int(np.unique(values[finite]).size) if finite.any() else 0,
                    "min": float(values[finite].min()) if finite.any() else 0.0,
                    "max": float(values[finite].max()) if finite.any() else 0.0})
    # An obs column that simply repeats the embedding above is the SAME axis under
    # another name. The approximate-projection UMAP is written into obs and is also
    # what the panel draws by default, so it appeared twice. The test compares
    # VALUES, not names, so any future duplicate is caught the same way.
    for entry in _numeric_obs_columns(cache):
        values = _axis_values(cache, entry["field"])
        if values is not None and any(
                values.shape == ref.shape
                and np.allclose(values, ref, rtol=0, atol=1e-6, equal_nan=True)
                for ref in (x, y)):
            continue
        entry.setdefault("label", _axis_label(entry["field"]))
        out.append(entry)
    return out


# A stored coordinate column says what it is, not what it was called in a script.
# `umap_averaged_x` told a reader nothing about which embedding it holds.
# scanpy_x / scanpy_y keep the names Nathan chose, so they are not relabelled.
# The umap_averaged_* entries are gone: those columns held a mean of the wrong
# coordinate export and were deleted from the objects on 2026-09-01.
_AXIS_LABELS = {
    "umap_approx_x": "approximate-projection UMAP 1",
    "umap_approx_y": "approximate-projection UMAP 2",
}


def _axis_label(field: str) -> str:
    return _AXIS_LABELS.get(str(field), str(field))


def _umap_coordinate_options(cache: Dict[str, Any]) -> List[Dict[str, str]]:
    """The coordinate sets a UMAP panel may be drawn on.

    The first entry is the cellHarmony projection, which is what the panel has
    always drawn: the coordinates the alignment wrote, read from the job's
    `umap_coordinates` artifact.
    """
    options = [{"key": "", "label": "cellHarmony UMAP"}]
    for key in cache.get("obsm_keys", []) or []:
        options.append({"key": str(key), "label": str(key)})
    # A coordinate set is a PAIR, so a stored `<stem>_x` / `<stem>_y` pair is one
    # choice, not two fields a reader has to assemble by hand. A precomputed
    # bundle keeps its embeddings as obs columns rather than in obsm, so without
    # this the alternative UMAPs were reachable only by picking both axes
    # manually from a list of 30-odd covariates.
    for stem in _obs_coordinate_pairs(cache):
        options.append({"key": OBS_PAIR_PREFIX + stem,
                        "label": stem.replace("_", " ")})
    return options


OBS_PAIR_PREFIX = "obspair:"


def _obs_coordinate_pairs(cache: Dict[str, Any]) -> List[str]:
    """Stems of every `<stem>_x` / `<stem>_y` float obs pair, in a stable order."""
    fields = {c["field"] for c in _numeric_obs_columns(cache)}
    stems = sorted({f[:-2] for f in fields
                    if f.endswith("_x") and (f[:-2] + "_y") in fields})
    return stems


def _coordinates_for_key(cache: Dict[str, Any], coords_key: str = "") -> tuple:
    """(x, y, resolved key). An unknown key falls back to the cellHarmony one."""
    key = str(coords_key or "").strip()
    if not key:
        return cache["umap_x"], cache["umap_y"], ""
    if key.startswith(OBS_PAIR_PREFIX):
        stem = key[len(OBS_PAIR_PREFIX):]
        x = _axis_values(cache, stem + "_x")
        y = _axis_values(cache, stem + "_y")
        if x is not None and y is not None:
            return x, y, key
        return cache["umap_x"], cache["umap_y"], ""
    adata = cache["adata"]
    obsm = getattr(adata, "obsm", {}) or {}
    if key not in obsm:
        return cache["umap_x"], cache["umap_y"], ""
    matrix = np.asarray(obsm[key])
    if matrix.ndim != 2 or matrix.shape[1] < 2:
        return cache["umap_x"], cache["umap_y"], ""
    return (np.asarray(matrix[:, 0], dtype=float),
            np.asarray(matrix[:, 1], dtype=float), key)


def _labels_for_color_by(cache: Dict[str, Any], color_by: str = "") -> tuple:
    """(per-cell labels, resolved column). Empty means the cellHarmony states."""
    column = str(color_by or "").strip()
    if not column or column == str(cache["cluster_key"]):
        return cache["populations"], ""
    adata = cache["adata"]
    if column not in adata.obs.columns:
        return cache["populations"], ""
    values = (adata.obs[column].astype(str).str.strip()
              .replace({"nan": "", "None": ""}).to_numpy(dtype=str))
    return values, column


#: Where a chat question is read into one supported query. The model runs in the
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
    "patient_stratification": "differential",
    "shared_vs_state_specific": "differential",
    "contrast_specificity": "differential",
    "pathway_program": "goelite",
    "regulatory_driver": "network",
    "communication_rewiring": "ccc",
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
    case = str(differential.get("case_label") or config.get("case_label") or "Group 1").strip()
    control = str(differential.get("control_label") or config.get("control_label") or "Group 2").strip()
    return {"case": case, "control": control, "label": f"{case} versus {control}"}


def _chat_marker_table(meta: Dict, modality: str = "rna") -> pd.DataFrame:
    """The marker table cellHarmony wrote, or an empty frame."""
    marker_analysis = _modality_marker_analysis(meta, modality) or {}
    path = Path(str(marker_analysis.get("markers_tsv", "")).strip())
    if not path.is_file():
        return pd.DataFrame()
    frame = pd.read_csv(path, sep="\t")
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
        # Not "_p" / "_fold": itertuples renames any column starting with an
        # underscore to a positional name, and the attribute then does not exist.
        if "FDR p-value" in subset.columns:
            subset["rank_p"] = pd.to_numeric(subset["FDR p-value"], errors="coerce")
        else:
            subset["rank_p"] = np.nan
        subset["rank_fold"] = pd.to_numeric(subset.get("Fold"), errors="coerce")
        subset = subset.sort_values(["rank_p", "rank_fold"], ascending=[True, False])
        for row in subset.itertuples():
            rows.append({"gene": str(row.Gene), "cluster": str(row.cluster),
                         "fold": float(row.rank_fold) if _is_finite_number(row.rank_fold) else 0.0,
                         "p": float(row.rank_p) if _is_finite_number(row.rank_p) else None,
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

    from .cross_modal import examples as cross_examples
    examples: List[str] = cross_examples(_modality_artifacts(meta))
    if first:
        marker_state = next((s for s in states if s.casefold() == "hsc-1"), first)
        examples.extend([f"What pathways have the best cross-modality {marker_state} representation?", f"What is the best modality marker of {marker_state}?"])
    if first:
        examples.append(f"What are the best marker genes of {first} cells?")
    if first and second:
        examples.append(f"What distinguishes {first} from {second} cells?")
    if genes:
        examples.append(f"Where is {genes[0]} expressed?")
    if len(genes) > 1:
        examples.append(f"Which cell states express {genes[1]}?")
    if first and "grn_tf" in _modality_artifacts(meta):
        examples.append(f"Which TFs are most active in {first} cells?")
        examples.append(f"Show the regulatory network in {first} cells")
    contrast = _chat_contrast(meta)
    if contrast and first:
        examples.append(f"Which genes are significant in {contrast['case']} versus "
                        f"{contrast['control']} in {first} cells?")
        examples.append(f"Which cell type is most affected in {contrast['case']} versus "
                        f"{contrast['control']}?")
    covariates=cache['adata'].obs
    numeric=next((c for c in covariates if pd.api.types.is_numeric_dtype(covariates[c]) and c not in {'n_counts','n_genes','n_cells'}), '')
    if first and numeric:examples.append(f"Which genes in {first} track {numeric}?")
    if first and genes:examples.append(f"Which genes coexpress with {genes[0]} in {first}?")
    if first and contrast:examples.append(f"Is the signature in {first} present in all donors or a subset?")
    return {"tissue": tissue, "reference": reference_label,
            "cluster_key": str(cache["cluster_key"]),
            "n_states": len(states), "has_contrast": bool(contrast),
            "examples": examples,
            "placeholder": (f"e.g. What are the best marker genes of {first} cells?"
                            if first else "e.g. What are the best marker genes of this cell state?")}


def _chat_states_in_question(question: str, states: List[str]) -> List[str]:
    """The cell states this sentence actually names, in the order it names them.

    The router matches a state name anywhere in the sentence, so "pre-aceNKP"
    also reports "aceNKP" and a two-state question came back comparing a state
    with itself. A match that sits inside a longer match is dropped here, which
    leaves the states the reader wrote.
    """
    text = str(question or "").lower()
    spans = []
    for state in sorted([str(s) for s in states if s], key=len, reverse=True):
        needle = state.lower()
        start = 0
        while True:
            at = text.find(needle, start)
            if at < 0:
                break
            if not any(begin <= at and at + len(needle) <= end for begin, end, _ in spans):
                spans.append((at, at + len(needle), state))
            start = at + 1
    spans.sort()
    out: List[str] = []
    for _, _, state in spans:
        if state not in out:
            out.append(state)
    return out


#: Words that are also gene symbols in some annotations. The question scan below
#: only runs when the router found no gene, and these would turn an ordinary
#: sentence into a gene lookup.
_CHAT_GENE_STOPWORDS = {
    "and", "are", "can", "cell", "cells", "for", "gene", "genes", "has", "how",
    "impact", "many", "max", "most", "not", "rest", "set", "she", "state",
    "states", "the", "was", "what", "when", "where", "which", "who", "why",
}


def _chat_genes_in_question(question: str, cache: Dict[str, Any]) -> List[str]:
    """Genes this sentence names, matched against the dataset's own gene list.

    The router returns an empty gene list for a plain sentence such as "Where is
    Cdca3 expressed?". Nothing is invented here: a token is only taken when the
    dataset holds a gene of that name.
    """
    index = {}
    for name in cache.get("var_names", []):
        index.setdefault(str(name).lower(), str(name))
    out: List[str] = []
    for token in re.findall(r"[A-Za-z0-9_.\-]{3,}", str(question or "")):
        key = token.lower()
        if key in _CHAT_GENE_STOPWORDS:
            continue
        name = index.get(key)
        if name and name not in out:
            out.append(name)
    return out


def _chat_read_question(question: str, cache: Dict[str, Any], meta: Dict) -> Dict[str, Any]:
    """Ask the assistant which supported query this sentence means.

    The model sees the question and the names in this dataset. It never sees the
    data, so it cannot invent a number: it chooses the reading, and the executors
    below compute the answer from the job's own files.
    """
    reading = gnet.read_regulatory_question(question, [s for s, _ in _chat_states_by_size(cache)],
                                            cache.get("var_names", []))
    if reading is not None:
        return reading
    from .chat_service import read_question
    local=read_question(question,[s for s,_ in _chat_states_by_size(cache)],cache.get('var_names',[]),cache['adata'].obs.columns)
    if local is not None:return local
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


def _build_umap_payload(
    app: FastAPI,
    meta: Dict,
    modality: str = "rna",
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
    color_by: str = "",
    coords_key: str = "",
    x_field: str = "",
    y_field: str = "",
) -> Dict[str, List[Dict]]:
    cache_entry = _get_expression_cache(app, meta, modality=modality)
    obs_names = cache_entry["obs_names"]
    populations, resolved_color_by = _labels_for_color_by(cache_entry, color_by)
    sample_field = str(cache_entry.get("sample_field") or "").strip()
    sample_labels = cache_entry.get("sample_labels")
    umap_x, umap_y, resolved_coords = _coordinates_for_key(cache_entry, coords_key)
    # A pair of numeric obs columns replaces the embedding outright. Both have to
    # resolve: one axis alone would silently mix a metadata value against a UMAP
    # coordinate, which reads as a map and is not one.
    axis_x = _axis_values(cache_entry, x_field)
    axis_y = _axis_values(cache_entry, y_field)
    resolved_x, resolved_y = "", ""
    if axis_x is not None and axis_y is not None:
        umap_x, umap_y = axis_x, axis_y
        resolved_x, resolved_y = str(x_field).strip(), str(y_field).strip()
        resolved_coords = ""
    axes_source = "obs" if resolved_x else ("obsm" if resolved_coords else "cellharmony")
    display_mask = _apply_display_filter_mask(cache_entry, display_filters)

    query_points = [
        {
            "barcode": barcode,
            "population": population,
            "sample": str(sample) if sample_labels is not None else "",
            "x": float(x),
            "y": float(y),
        }
        for barcode, population, sample, x, y, keep in zip(
            obs_names,
            populations,
            sample_labels if sample_labels is not None else np.repeat("", len(obs_names)),
            umap_x,
            umap_y,
            display_mask,
        )
        if keep and _is_finite_number(x) and _is_finite_number(y)
    ]

    # The reference atlas is drawn behind the query only while both choices are
    # the default. Another obs column has no counterpart in the reference, and a
    # second embedding is a different coordinate space, so an overlay there would
    # put the reference cells in positions that mean nothing.
    reference_points = []
    ref_adata = (None if (resolved_color_by or resolved_coords or resolved_x)
                 else _load_reference_adata(app, meta))
    if ref_adata is not None:
        ref_cluster_key = meta.get("reference_cluster_key") or meta.get("cluster_key")
        coords = ref_adata.obsm["X_umap"]
        labels = ref_adata.obs[ref_cluster_key].astype(str).tolist()
        for barcode, (x, y), population in zip(ref_adata.obs_names, coords, labels):
            if not (_is_finite_number(x) and _is_finite_number(y)):
                continue
            reference_points.append(
                {
                    "barcode": str(barcode),
                    "population": population,
                    "x": float(x),
                    "y": float(y),
                }
            )
    # How many cells the panel could not place. A metadata axis is often
    # recorded for part of the dataset only, and a silently shorter plot would
    # read as a real absence of cells.
    n_kept = int(np.count_nonzero(np.asarray(display_mask, dtype=bool)))
    return {"reference": reference_points, "query": query_points,
            "sample_field": sample_field,
            "color_by": resolved_color_by,
            "color_label": resolved_color_by or str(cache_entry["cluster_key"]),
            "coords_key": resolved_coords,
            "coords_label": (resolved_coords[len(OBS_PAIR_PREFIX):].replace("_", " ")
                             if resolved_coords.startswith(OBS_PAIR_PREFIX)
                             else (resolved_coords or "cellHarmony UMAP")),
            "axes_source": axes_source,
            "x_field": resolved_x,
            "y_field": resolved_y,
            "x_label": resolved_x or "UMAP 1",
            "y_label": resolved_y or "UMAP 2",
            "n_cells_selected": n_kept,
            "n_points_drawn": len(query_points),
            "n_dropped_no_coordinate": max(0, n_kept - len(query_points)),
            "reference_hidden": bool(resolved_color_by or resolved_coords or resolved_x)}


def _reference_entry_for_meta(meta: Dict) -> Dict:
    registry_path = Path(load_config().get("REFERENCE_REGISTRY"))
    reference_entry = pipeline_mod._lookup_reference(meta["species"], meta["reference"], registry_path)
    pipeline_mod._ensure_reference_fields(reference_entry)
    return reference_entry


def _load_reference_states_table(meta: Dict) -> pd.DataFrame:
    reference_entry = _reference_entry_for_meta(meta)
    states_path = Path(reference_entry["states_tsv"])
    if not states_path.exists():
        raise FileNotFoundError("Reference states TSV is unavailable for this job.")
    return pd.read_csv(states_path, sep="\t", index_col=0)


def _build_reference_expression_payload(app: FastAPI, meta: Dict, gene: str) -> Dict:
    reference_df = _load_reference_states_table(meta)
    resolved_gene = _resolve_gene_name(reference_df.index.astype(str), gene)
    if not resolved_gene:
        return {
            "gene": gene,
            "requested_gene": gene,
            "resolved_gene": None,
            "source": "missing",
            "message": f"Gene '{gene}' was not found in the aligned data or the selected reference.",
            "scatter": [],
            "violin": [],
            "umap": [],
        }

    series = pd.to_numeric(reference_df.loc[resolved_gene], errors="coerce").dropna()
    ref_adata = _load_reference_adata(app, meta)
    ref_cluster_key = meta.get("reference_cluster_key") or meta.get("cluster_key")

    umap_points = []
    scatter_data = []
    if ref_adata is not None and ref_cluster_key and ref_cluster_key in ref_adata.obs.columns:
        coords = np.asarray(ref_adata.obsm["X_umap"])
        populations = ref_adata.obs[ref_cluster_key].astype(str).tolist()
        value_map = {str(pop): float(val) for pop, val in series.items() if _is_finite_number(val)}
        for barcode, (x, y), population in zip(ref_adata.obs_names, coords, populations):
            if not (_is_finite_number(x) and _is_finite_number(y)):
                continue
            value = float(value_map.get(str(population), 0.0))
            umap_points.append(
                {
                    "barcode": str(barcode),
                    "population": str(population),
                    "value": value,
                    "x": float(x),
                    "y": float(y),
                }
            )
        scatter_data = [
            {"population": str(pop), "value": float(val)}
            for pop, val in value_map.items()
            if _is_finite_number(val)
        ]

    violin_data = [
        {
            "population": str(pop),
            "values": [float(val)],
            "mean": float(val),
        }
        for pop, val in sorted(series.items(), key=lambda item: float(item[1]), reverse=True)[:10]
        if _is_finite_number(val)
    ]
    series_values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float, copy=False)
    finite_series = series_values[np.isfinite(series_values)]
    global_min = float(np.min(finite_series)) if finite_series.size else 0.0
    global_max = float(np.max(finite_series)) if finite_series.size else 0.0

    return {
        "gene": resolved_gene,
        "requested_gene": gene,
        "resolved_gene": resolved_gene,
        "source": "reference",
        "message": (
            f"Showing reference expression for '{resolved_gene}' because the gene was not found "
            "in the aligned AnnData output."
        ),
        "scatter": scatter_data,
        "violin": violin_data,
        "umap": umap_points,
        "global_min": global_min,
        "global_max": global_max,
    }


def split_gene_list(text: str, keep_pipe: bool = False) -> List[str]:
    """Gene symbols out of whatever a user pasted.

    A column copied out of Excel arrives newline separated, a row tab separated,
    and people also type commas, spaces and semicolons. Splitting on commas alone
    turned "SFTPC AGER" into one symbol of that name and found nothing.
    Duplicates are dropped and order is kept, so a figure reads in the order the
    genes were listed.

    `keep_pipe` stops the vertical bar being a separator. A GRN edge feature is named
    `TF|target`, so splitting on the bar turned one edge name into two names the
    modality does not hold.
    """
    pattern = r"[\s,;]+" if keep_pipe else r"[\s,;|]+"
    parts = re.split(pattern, str(text or ""))
    seen, out = set(), []
    for part in parts:
        gene = part.strip().strip('"').strip("'")
        if gene and gene not in seen:
            seen.add(gene)
            out.append(gene)
    return out


def _split_expression_features(text, cache):
    """Preserve literal metabolite names before falling back to gene separators."""
    raw = str(text or "").strip()
    names = set(map(str, cache["var_names"]))
    if raw in names:
        return [raw]
    for separator in (r"[\n\t]+", r"[,;\n\t]+"):
        parts = [part.strip().strip('"').strip("'") for part in re.split(separator, raw) if part.strip()]
        if parts and all(part in names for part in parts):
            return list(dict.fromkeys(parts))
    return split_gene_list(raw, keep_pipe=_names_have_pipe(cache))


def _names_have_pipe(cache: Dict[str, Any]) -> bool:
    """True when this modality's feature names carry `|`, as a GRN edge does."""
    cached = cache.get("names_have_pipe")
    if cached is None:
        names = cache.get("var_names")
        cached = bool(names is not None and any("|" in str(name) for name in names))
        cache["names_have_pipe"] = cached
    return bool(cached)


def _gene_rows(cache: Dict[str, Any], wanted: List[str]) -> tuple:
    """Row index of each requested gene, and the ones this dataset lacks."""
    var_names = cache["var_names"]
    index = {str(name): i for i, name in enumerate(var_names)}
    upper = {str(name).upper(): i for i, name in enumerate(var_names)}
    rows, labels, missing = [], [], []
    for gene in wanted:
        row = index.get(gene, upper.get(gene.upper()))
        if row is None:
            missing.append(gene)
        else:
            rows.append(row)
            labels.append(gene)
    return rows, labels, missing


def _dense_column(adata, row: int) -> np.ndarray:
    """One gene's values for every cell, dense."""
    column = adata.X[:, row]
    if sp.issparse(column):
        return np.asarray(column.todense()).ravel()
    return np.asarray(column).ravel()


def _default_marker_genes(cache: Dict[str, Any], group_by: str = "",
                          limit: int = 12) -> List[str]:
    """One marker gene per group, used when the user gives no gene set.

    Picked as the gene with the largest gap between its mean inside a group and
    its mean everywhere else, over a sample of genes. That is a plain
    one-versus-rest contrast, not a stored marker table, so it works for any
    grouping variable the user switches to, not only cell state.
    """
    adata = cache["adata"]
    _, groups, values_of = _group_axis(cache, group_by)
    var_names = [str(v) for v in cache["var_names"]]
    if not groups or not var_names:
        return []
    groups = list(groups[:limit])

    # One pass over the matrix, not one column slice per (gene, group). The
    # previous version sliced `adata.X[:, row]` inside a nested loop: on a
    # 2,797-cell CSR matrix that is 7.7 ms per slice, 4,764 sampled genes x 12
    # groups = 57,168 slices, measured at 7.3 minutes for one figure. The
    # endpoint held the event loop for all of it, so every other panel in the
    # app froze. An indicator matrix multiplied into X gives the same sums for
    # every gene at once, and it reads every gene rather than a sample of them.
    X = adata.X
    n_cells, n_genes = X.shape
    index = {g: i for i, g in enumerate(groups)}
    codes = pd.Series(values_of).map(index).fillna(-1).to_numpy(dtype=np.int64)
    inside_any = codes >= 0
    if not inside_any.any():
        return []
    rows = np.nonzero(inside_any)[0]
    indicator = sp.csr_matrix(
        (np.ones(rows.size, dtype=np.float64), (codes[rows], rows)),
        shape=(len(groups), n_cells))
    sums = indicator @ X
    sums = np.asarray(sums.todense() if sp.issparse(sums) else sums, dtype=np.float64)
    counts = np.bincount(codes[rows], minlength=len(groups)).astype(np.float64)
    # The contrast is against every other cell in the dataset, which is what the
    # previous `values[~inside]` measured, so the statistic is unchanged.
    total = np.asarray(X.sum(axis=0), dtype=np.float64).ravel()
    inside_mean = sums / np.maximum(counts[:, None], 1.0)
    outside_mean = (total[None, :] - sums) / np.maximum(float(n_cells) - counts[:, None], 1.0)
    gap = inside_mean - outside_mean

    chosen, seen = [], set()
    for position, group in enumerate(groups):
        if counts[position] <= 0:
            continue
        for row in np.argsort(-gap[position])[:200]:
            if gap[position][int(row)] <= 0:
                break
            name = var_names[int(row)]
            if name in seen:
                continue
            seen.add(name)
            chosen.append(name)
            break
    return chosen


def _subset_mask(cache: Dict[str, Any], subset_by: str = "",
                 subset_values: Optional[List[str]] = None) -> Optional[np.ndarray]:
    """Which cells a figure is restricted to, or None for all of them.

    This is what makes a cell-state-specific contrast possible. Grouping by
    copd_status alone contrasts COPD against control over every cell in the
    atlas, which mixes 39 cell types together. Restricting to AT2 first and then
    grouping by copd_status asks the question the user meant: within AT2, what
    differs between disease and control.
    """
    if not subset_by or not subset_values:
        return None
    adata = cache["adata"]
    if subset_by not in adata.obs.columns:
        return None
    values = adata.obs[subset_by].astype(str).to_numpy()
    keep = set(str(v) for v in subset_values)
    mask = np.isin(values, list(keep))
    return mask


def _gene_set_filter_mask(cache, subset_by="", subset_values=None, subset2_by="", subset2_values=None):
    mask = np.ones(cache["adata"].n_obs, dtype=bool)
    for field, values in ((subset_by, subset_values), (subset2_by, subset2_values)):
        selected = _subset_mask(cache, field, values)
        if selected is not None:
            mask &= selected
    return mask


def _group_axis(cache: Dict[str, Any], group_by: str = "") -> tuple:
    """The obs column the DotPlot and CombPlot group their columns by.

    Cell state is the default. Choosing another variable regroups the whole
    figure by that variable's categories, so the same genes can be read across
    disease, sex or assay instead of across cell types. The categories keep the
    dataset's own order when the column is categorical, which is the centroid
    order for cell state.
    """
    adata = cache["adata"]
    cluster_key = str(cache["cluster_key"])
    column = str(group_by or "").strip() or cluster_key
    if column not in adata.obs.columns:
        column = cluster_key
    series = adata.obs[column]
    if str(series.dtype) == "category":
        groups = [str(c) for c in series.cat.categories]
    elif column == cluster_key and len(cache.get("populations", [])):
        # `populations` is a numpy array. `array or []` calls bool() on it and
        # raises "truth value of an array ... is ambiguous", which took down
        # every DotPlot and CombPlot request.
        groups = list(dict.fromkeys(str(s) for s in cache["populations"]))
    else:
        groups = sorted({str(v) for v in series.astype(str).to_numpy()})
    return column, groups, series.astype(str).to_numpy()


def _groupable_columns(cache: Dict[str, Any], max_categories: int = 60) -> List[Dict[str, Any]]:
    """The obs columns a user may group or filter by.

    Only categorical-like columns with a workable number of levels are offered.
    A per-cell numeric column such as n_counts has one level per cell and would
    draw a column per cell, so it is left out.
    """
    adata = cache["adata"]
    cluster_key = str(cache["cluster_key"])
    _, states, _ = _group_axis(cache, cluster_key)
    out = [{"field": cluster_key, "values": list(states), "n": len(states)}]
    for name in adata.obs.columns:
        if str(name) == cluster_key:
            continue
        series = adata.obs[name]
        if str(series.dtype) == "category":
            levels = [str(c) for c in series.cat.categories]
        elif series.dtype == object or str(series.dtype).startswith(("bool", "str")):
            levels = sorted({str(v) for v in series.astype(str).to_numpy()})
        else:
            continue
        if 1 < len(levels) <= max_categories:
            out.append({"field": str(name), "values": levels, "n": len(levels)})
    return out


def _gene_state_stats(cache: Dict[str, Any], wanted: List[str],
                      group_by: str = "", keep_groups: Optional[List[str]] = None,
                      subset_by: str = "",
                      subset_values: Optional[List[str]] = None,
                      subset2_by: str = "", subset2_values: Optional[List[str]] = None) -> Dict[str, Any]:
    """Mean expression and detected fraction of each gene in each group.

    Groups are cell states by default. `group_by` regroups by any other
    categorical variable, and `keep_groups` restricts the figure to the groups
    the user selected, which is how the cell-state filter narrows the plot.
    """
    adata = cache["adata"]
    column, groups, values_of = _group_axis(cache, group_by)
    if keep_groups:
        chosen = [g for g in groups if g in set(keep_groups)]
        groups = chosen
    restrict = _gene_set_filter_mask(cache, subset_by, subset_values, subset2_by, subset2_values)
    rows, labels, missing = _gene_rows(cache, wanted)
    masks = [(values_of == group) for group in groups]
    if restrict is not None:
        masks = [m & restrict for m in masks]
    counts = [int(m.sum()) for m in masks]
    mean, frac = [], []
    for row in rows:
        column_values = _dense_column(adata, row)
        mean.append([float(column_values[m].mean()) if n else 0.0
                     for m, n in zip(masks, counts)])
        frac.append([float((column_values[m] > 0).mean()) if n else 0.0
                     for m, n in zip(masks, counts)])
    return {"genes": labels, "states": groups, "groups": groups,
            "group_by": column, "group_label": column,
            "subset_by": subset_by or "", "subset_values": list(subset_values or []),
            "state_n": counts, "mean": mean, "frac": frac,
            "colors": _state_colors(cache, groups),
            "n_requested": len(wanted), "n_returned": len(labels),
            "n_missing": len(missing), "missing": missing}

def _sample_plot_cells(cache, indices, cells_per_sample):
    """Sample each sample × cell-type stratum, preserving the supplied plot order."""
    from altanalyze3.components.visualization.cell_sampling import sample_cell_indices
    obs = cache["adata"].obs
    # Match the viewer's donor/sample identity order; a display-filter default
    # can be a condition or sex and must not silently become a sample identity.
    sample_field = next((c for c in ("meta_sample", "donor", "Donor", "Library", "sample", "sample_id", "pool") if c in obs), "")
    samples = obs[sample_field].astype(str).to_numpy() if sample_field else np.repeat("dataset", len(obs))
    ids = np.asarray(cache.get("obs_names", obs.index), dtype=str)
    states = np.asarray(cache["populations"], dtype=str)
    try:
        keep = sample_cell_indices(ids[indices], samples[indices], limit=cells_per_sample, group_labels=states[indices])
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    description = (f"Up to {cells_per_sample} individual cells per sample per cell type"
                   if cells_per_sample else "All individual cells")
    description += f"; sample annotation: {sample_field}." if sample_field else "; no sample annotation, treating the dataset as one sample."
    return indices[keep], {"cells_per_sample": cells_per_sample, "sample_field": sample_field,
                           "n_available": len(indices), "n_selected": len(keep), "description": description}


def _gene_cell_values(cache: Dict[str, Any], wanted: List[str], group_by: str = "",
                      keep_groups: Optional[List[str]] = None, subset_by: str = "",
                      subset_values: Optional[List[str]] = None,
                      subset2_by: str = "", subset2_values: Optional[List[str]] = None, cells_per_sample: int = 0) -> Dict[str, Any]:
    """Individual cell values in group order, without donor averaging or sampling."""
    adata = cache["adata"]
    column, groups, values_of = _group_axis(cache, group_by)
    if keep_groups:
        groups = [g for g in groups if g in set(keep_groups)]
    valid = np.isin(values_of, groups)
    for field, values in ((subset_by, subset_values), (subset2_by, subset2_values)):
        if field and values and field in adata.obs.columns:
            valid &= adata.obs[field].astype(str).isin(values).to_numpy()
    order = {g: i for i, g in enumerate(groups)}
    indices = np.flatnonzero(valid)
    indices = indices[np.argsort([order[values_of[i]] for i in indices], kind="stable")]
    indices, sampling = _sample_plot_cells(cache, indices, cells_per_sample)
    sample = cache.get("sample_field")
    donors = adata.obs[sample].astype(str).to_numpy() if sample in adata.obs.columns else None
    names = cache.get("obs_names", adata.obs.index)
    columns = [{"cell": str(names[i]), "group": str(values_of[i]), "state": str(values_of[i]),
                "donor": str(donors[i]) if donors is not None else "", "n_cells": 1}
               for i in indices]
    rows, labels, missing = _gene_rows(cache, wanted)
    series = [np.round(_dense_column(adata, row)[indices].astype(float), 5).tolist() for row in rows]
    return {"genes": labels, "values": series, "columns": columns,
            "colors": _state_colors(cache, [c["group"] for c in columns]),
            "unit": "cells", "observation_unit": "cells", "sampling": sampling, "groups": groups, "states": groups,
            "group_by": column, "group_label": column, "n_columns": len(columns),
            "n_cells_kept": len(columns), "n_groups_dropped": 0,
            "n_requested": len(wanted), "n_returned": len(labels),
            "n_missing": len(missing), "missing": missing}


def _gene_donor_state_means(cache: Dict[str, Any], wanted: List[str],
                            min_cells: int, group_by: str = "",
                            keep_groups: Optional[List[str]] = None,
                            subset_by: str = "",
                            subset_values: Optional[List[str]] = None,
                            subset2_by: str = "", subset2_values: Optional[List[str]] = None) -> Dict[str, Any]:
    """Per-donor pseudobulk of each gene, within each group.

    One column per (group, donor). Groups are cell states by default; `group_by`
    regroups by any other categorical variable and `keep_groups` restricts the
    figure to the selected groups.
    """
    adata = cache["adata"]
    donor_field = cache.get("sample_field") or ""
    if not donor_field or donor_field not in adata.obs.columns:
        for candidate in ("meta_sample", "donor", "Donor", "Library", "sample", "pool"):
            if candidate in adata.obs.columns:
                donor_field = candidate
                break
    if not donor_field or donor_field not in adata.obs.columns:
        return {"error": "this dataset records no donor column"}

    column, groups, values_of = _group_axis(cache, group_by)
    if keep_groups:
        chosen = [g for g in groups if g in set(keep_groups)]
        groups = chosen
    donor_of = adata.obs[donor_field].astype(str).to_numpy()
    donors = sorted(set(donor_of))
    group_index = {g: i for i, g in enumerate(groups)}
    donor_index = {d: i for i, d in enumerate(donors)}

    restrict = _subset_mask(cache, subset_by, subset_values)
    codes = np.array([group_index.get(v, -1) for v in values_of], dtype=np.int64)
    if restrict is not None:
        codes = np.where(restrict, codes, -1)
    if subset2_by and subset2_values and subset2_by in adata.obs.columns:
        codes = np.where(adata.obs[subset2_by].astype(str).isin(subset2_values), codes, -1)
    donor_codes = np.array([donor_index[d] for d in donor_of], dtype=np.int64)
    group_id = codes * len(donors) + donor_codes
    valid = codes >= 0
    n_groups = len(groups) * len(donors)
    per_group = np.bincount(group_id[valid], minlength=n_groups)
    keep = np.nonzero(per_group >= min_cells)[0]
    dropped = int(np.count_nonzero((per_group > 0) & (per_group < min_cells)))
    if not keep.size:
        return {"error": f"no donor contributes {min_cells} or more cells to a group"}

    columns = [{"group": groups[int(g) // len(donors)],
                "state": groups[int(g) // len(donors)],
                "donor": donors[int(g) % len(donors)],
                "n_cells": int(per_group[int(g)])} for g in keep]
    colors = _state_colors(cache, [c["group"] for c in columns])

    rows, labels, missing = _gene_rows(cache, wanted)
    series = []
    for row in rows:
        column_values = _dense_column(adata, row)
        sums = np.zeros(n_groups, dtype=np.float64)
        np.add.at(sums, group_id[valid], column_values[valid])
        means = np.where(per_group > 0, sums / np.maximum(per_group, 1), 0.0)
        series.append([round(float(v), 5) for v in means[keep]])

    return {"unit": "donor", "genes": labels, "values": series, "columns": columns, "colors": colors,
            "states": groups, "groups": groups,
            "subset_by": subset_by or "", "subset_values": list(subset_values or []),
            "group_by": column, "group_label": column, "donor_key": donor_field,
            "n_donors": len(donors), "n_states": len(groups),
            "n_columns": len(columns), "n_groups_dropped": dropped,
            "n_requested": len(wanted), "n_returned": len(labels),
            "n_missing": len(missing), "missing": missing}

#: The 12 Paired colours the Explore panel already uses for cell states, as hex.
#: Kept identical to PAIRED_COLOR_STOPS in static/app.js so a state is the same
#: colour whichever side of the app drew it.
_PAIRED_HEX = [
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
]


def _state_colors(cache: Dict[str, Any], states: List[str]) -> List[str]:
    """The colour each cell state is drawn in, falling back to a neutral grey.

    `adata.uns[f"{cluster_key}_colors"]` holds one colour per cell-state
    CATEGORY. Two defects lived in the previous version. It read
    `cache["populations"]`, which is one label per CELL, so the length test
    never matched and every state came back grey. And it wrote
    `cache.get("populations") or []`, which calls bool() on a numpy array and
    raises, so both gene-set figures returned HTTP 500 on every request.
    """
    adata = cache["adata"]
    cluster_key = cache["cluster_key"]
    stored = adata.uns.get(f"{cluster_key}_colors")
    order = _state_order(cache)
    lookup = {}
    if stored is not None and len(stored) == len(order):
        lookup = {s: str(c) for s, c in zip(order, stored)}
    else:
        # cellHarmony output h5ads carry no `<cluster_key>_colors`, so every bar
        # of the CombPlot came back grey and its state blocks were unreadable.
        # The fallback is the same Paired palette the front end draws cell states
        # with, assigned by position in the dataset's own state order, so a state
        # keeps one colour across every figure.
        lookup = {state: _PAIRED_HEX[index % len(_PAIRED_HEX)]
                  for index, state in enumerate(order)}
    return [lookup.get(s, "#BBBBBB") for s in states]


def _state_order(cache: Dict[str, Any]) -> List[str]:
    """The cell-state categories, in the dataset's own order.

    The categorical order is the centroid order every other figure reads, so a
    dataset that stores the column as a plain string column keeps first-seen
    order rather than being re-sorted alphabetically here.
    """
    adata = cache["adata"]
    cluster_key = str(cache["cluster_key"])
    series = adata.obs[cluster_key] if cluster_key in adata.obs.columns else None
    if series is not None and str(series.dtype) == "category":
        return [str(c) for c in series.cat.categories]
    return list(dict.fromkeys(str(s) for s in cache.get("populations", [])))



def _as_int(value: Any, fallback: int, low: int, high: int) -> int:
    """An int within bounds, whatever arrived.

    A route called directly rather than through FastAPI receives its own
    `Query(...)` default object, and int() on that raises. Anything that will
    not convert falls back rather than breaking the request.
    """
    try:
        number = int(value)
    except (TypeError, ValueError):
        number = fallback
    return max(low, min(number, high))


def _build_expression_payload(
    app: FastAPI,
    meta: Dict,
    gene: str,
    modality: str = "rna",
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
    violin_limit: int = 10,
    x_field: str = "",
    y_field: str = "",
) -> Dict:
    """`violin_limit` is how many cell states the violin plot draws.

    Ten fits the half-width panel the two-window layout gives. With one window
    the panel is twice as wide, so the viewer asks for more and the plot uses
    the space instead of leaving it empty.
    """
    def _expression_global_range(raw_values: np.ndarray) -> tuple[float, float]:
        finite = np.asarray(raw_values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return 0.0, 0.0
        return float(np.min(finite)), float(np.max(finite))

    normalized_modality = _normalize_modality_id(modality)
    cache_entry = _get_expression_cache(app, meta, modality=normalized_modality)
    adata = cache_entry["adata"]
    populations = cache_entry["populations"]
    obs_names = cache_entry["obs_names"]
    umap_x = cache_entry["umap_x"]
    umap_y = cache_entry["umap_y"]
    # Every UMAP panel takes the same coordinate control, not only the cell-type
    # one. Nathan, 2026-09-01: "all UMAP plots should have the option to change
    # umap coordinates". Both axes must resolve, or the panel would mix a
    # metadata value against a UMAP axis and read as a map that is not one.
    _ax = _axis_values(cache_entry, x_field)
    _ay = _axis_values(cache_entry, y_field)
    if _ax is not None and _ay is not None:
        umap_x, umap_y = _ax, _ay
    display_mask = _apply_display_filter_mask(cache_entry, display_filters)
    resolved_gene = _resolve_gene_name(cache_entry["var_names"], gene)
    if not resolved_gene and not str(gene or "").strip() and len(cache_entry["var_names"]):
        resolved_gene = str(cache_entry["var_names"][0])
    if not resolved_gene:
        if normalized_modality == "rna":
            try:
                return _build_reference_expression_payload(app, meta, gene)
            except (FileNotFoundError, KeyError, ValueError):
                pass
        raise KeyError(f"Feature '{gene}' not found in the selected modality output.")

    values = _flatten_expr(adata[:, resolved_gene].X)
    global_min, global_max = _expression_global_range(values)
    scatter_data = [
        {"population": pop, "value": float(val)}
        for pop, val, keep in zip(populations, values, display_mask)
        if keep and _is_finite_number(val)
    ]

    umap_points = [
        {
            "barcode": barcode,
            "population": pop,
            "value": float(val),
            "x": float(x),
            "y": float(y),
        }
        for barcode, pop, val, x, y, keep in zip(obs_names, populations, values.astype(float), umap_x, umap_y, display_mask)
        if keep and _is_finite_number(val) and _is_finite_number(x) and _is_finite_number(y)
    ]
    umap_points.sort(key=lambda point: (point["value"], point["population"], point["barcode"]))

    violin_data = []
    for pop in sorted(pd.unique(populations)):
        mask = (populations == pop) & display_mask
        pop_values = values[mask].astype(float)
        finite_values = pop_values[np.isfinite(pop_values)]
        if not len(finite_values):
            continue
        violin_data.append(
            {
                "population": pop,
                "values": [float(v) for v in finite_values],
                "mean": float(np.mean(finite_values)) if len(finite_values) else 0.0,
            }
        )
    violin_data = sorted(violin_data, key=lambda x: x["mean"], reverse=True)[:violin_limit]

    return {
        "gene": resolved_gene,
        "requested_gene": gene,
        "resolved_gene": resolved_gene,
        "source": "query",
        "modality": normalized_modality,
        "message": None if int(np.asarray(display_mask).sum()) else "No cells match the current Display only filters.",
        "scatter": scatter_data,
        "violin": violin_data,
        "umap": umap_points,
        "global_min": global_min,
        "global_max": global_max,
    }


def _build_gene_suggestions_payload(app: FastAPI, meta: Dict, modality: str = "rna") -> Dict:
    cache_entry = _get_expression_cache(app, meta, modality=modality)
    genes = [str(gene) for gene in cache_entry["var_names"].tolist()]
    modality_info = _modality_definition(meta, modality)
    return {
        "genes": genes,
        "modality": _normalize_modality_id(modality),
        "feature_label": str(modality_info.get("feature_label") or "gene"),
    }


def _build_display_filter_payload(app: FastAPI, meta: Dict, modality: str = "rna") -> Dict[str, Any]:
    cache_entry = _get_expression_cache(app, meta, modality=modality)
    return dict(cache_entry.get("display_filters_meta") or {})


def _normalize_display_filter_values(values: Optional[List[str]]) -> List[str]:
    if not values:
        return []
    ordered: List[str] = []
    seen: set[str] = set()
    for value in values:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            ordered.append(text)
    return ordered


def _display_filter_specs(
    filter1_field: Optional[str],
    filter1_values: Optional[List[str]],
    filter2_field: Optional[str],
    filter2_values: Optional[List[str]],
) -> List[tuple[str, List[str]]]:
    specs: List[tuple[str, List[str]]] = []
    for raw_field, raw_values in (
        (filter1_field, filter1_values),
        (filter2_field, filter2_values),
    ):
        field = str(raw_field or "").strip()
        values = _normalize_display_filter_values(raw_values)
        if field and values:
            specs.append((field, values))
    return specs


def _apply_display_filter_mask(cache_entry: Dict[str, Any], display_filters: Optional[List[tuple[str, List[str]]]]) -> np.ndarray:
    obs_names = cache_entry["obs_names"]
    mask = np.ones(len(obs_names), dtype=bool)
    if not display_filters:
        return mask
    obs_filter_values = cache_entry.get("obs_filter_values") or {}
    for field, values in display_filters:
        field_values = obs_filter_values.get(field)
        if field_values is None or not values:
            continue
        mask &= np.isin(field_values, values)
    return mask


def _sampled_marker_heatmap(app, meta, modality, display_filters, cells_per_sample):
    from altanalyze3.components.visualization.cell_sampling import CELL_SAMPLE_LIMITS
    if cells_per_sample not in (0, *CELL_SAMPLE_LIMITS):
        raise HTTPException(400, "Cells per sample must be 5, 10, 20, 50, or 0 for all cells.")
    try:
        expression = _get_expression_cache(app, meta, modality=modality)
        base = _get_marker_heatmap_cache_entry(app, meta, modality=modality)
    except FileNotFoundError as exc:
        raise HTTPException(404, f"Cell-level expression and marker genes are required for sampling: {exc}")
    # Reuse the same standardized values for HEAD, viewer GET, and PDF.
    key = f"{meta['job_id']}:{modality}:sampled:{cells_per_sample}:" + json.dumps(display_filters or [], sort_keys=True)
    cached = app.state.marker_heatmap_cache.get(key)
    if cached and cached.get("base_signature") == base.get("signature") and cached.get("expression_entry") is expression:
        return cached
    mask = _apply_display_filter_mask(expression, display_filters)
    indices = np.flatnonzero(mask)
    states = np.asarray(expression["populations"], dtype=str)
    row_ids = np.asarray(base["row_ids"], dtype=str)
    genes = [row.split(":", 1)[-1] for row in row_ids]
    group_order = list(dict.fromkeys(row.split(":", 1)[0] for row in row_ids))
    order = {group: i for i, group in enumerate(group_order)}
    indices = indices[np.argsort([order.get(states[i], len(order)) for i in indices], kind="stable")]
    indices, sampling = _sample_plot_cells(expression, indices, cells_per_sample)
    if not len(indices):
        raise HTTPException(404, "No cells match the current display filters for the marker heatmap.")
    adata = expression["adata"]
    names = np.asarray(expression["var_names"], dtype=str)
    lookup = {gene: i for i, gene in enumerate(names)}
    found = [(row, gene, lookup[gene]) for row, gene in zip(row_ids, genes) if gene in lookup]
    if not found:
        raise HTTPException(404, "The stored marker genes are unavailable in the expression matrix.")
    values = np.empty((len(found), len(indices)), dtype=np.float32)
    # Slice/upload matrices once and use gene-major access; the viewer already has
    # a gene-major sparse store. Avoid repeated scans of the full RNA CSR matrix.
    selected_matrix = adata.X[:, [col for _, _, col in found]] if isinstance(adata, ad.AnnData) else None
    if sp.issparse(selected_matrix):
        selected_matrix = selected_matrix.tocsc()
    for out, (_, gene, col) in enumerate(found):
        raw = selected_matrix[:, out] if selected_matrix is not None else adata[:, gene].X
        vector = _flatten_expr(raw).astype(float)
        finite = vector[np.isfinite(vector)]
        mean = float(finite.mean()) if finite.size else 0.
        std = float(finite.std()) if finite.size else 0.
        values[out] = np.nan_to_num((vector[indices] - mean) / std) if std > 0 else 0.
    ids = np.asarray(expression["obs_names"], dtype=str)
    entry = {"base_signature": base.get("signature"), "expression_entry": expression,
             "matrix": values, "row_ids": np.asarray([r for r, _, _ in found]),
             "col_ids": np.asarray([f"{states[i]}:{ids[i]}" for i in indices]),
             "col_barcodes": ids[indices], "sampling": sampling}
    # Bound interactive caches when readers try many filter/sample combinations.
    prefix = f"{meta['job_id']}:{modality}:sampled:"
    older = [k for k in app.state.marker_heatmap_cache if k.startswith(prefix)]
    total_bytes = values.nbytes + sum(app.state.marker_heatmap_cache[k]["matrix"].nbytes for k in older)
    while older and (len(older) >= 3 or total_bytes > 128 * 1024 * 1024):
        old = older.pop(0)
        removed = app.state.marker_heatmap_cache.pop(old)
        total_bytes -= removed["matrix"].nbytes
    app.state.marker_heatmap_cache[key] = entry
    return entry


def _filter_marker_heatmap_matrix(
    app: FastAPI,
    meta: Dict,
    modality: str,
    display_filters: Optional[List[tuple[str, List[str]]]],
) -> tuple[str, int]:
    cache_entry = _get_expression_cache(app, meta, modality=modality)
    marker_matrix = _get_marker_heatmap_cache_entry(app, meta, modality=modality)
    display_mask = _apply_display_filter_mask(cache_entry, display_filters)
    allowed_barcodes = np.asarray(
        [str(barcode) for barcode, keep in zip(cache_entry["obs_names"], display_mask) if bool(keep)],
        dtype=str,
    )
    if allowed_barcodes.size == 0:
        raise HTTPException(status_code=404, detail="No cells match the current display filters for the marker heatmap.")

    keep_mask = np.isin(marker_matrix["col_barcodes"], allowed_barcodes)
    keep_count = int(np.count_nonzero(keep_mask))
    if keep_count == 0:
        raise HTTPException(status_code=404, detail="No heatmap columns match the current display filters.")

    return (
        _marker_heatmap_subset_to_tsv(
            marker_matrix["matrix"][:, keep_mask],
            marker_matrix["row_ids"],
            marker_matrix["col_ids"][keep_mask],
        ),
        keep_count,
    )


def _square_umap_axes(ax: plt.Axes, point_groups: List[List[Dict]]) -> None:
    coords = [
        (float(point["x"]), float(point["y"]))
        for group in point_groups
        for point in group
        if _is_finite_number(point.get("x")) and _is_finite_number(point.get("y"))
    ]
    if not coords:
        return

    points = np.asarray(coords, dtype=float)
    min_x, max_x = float(points[:, 0].min()), float(points[:, 0].max())
    min_y, max_y = float(points[:, 1].min()), float(points[:, 1].max())
    span_x = max(max_x - min_x, 1e-6)
    span_y = max(max_y - min_y, 1e-6)
    half_span = max(span_x, span_y) / 2.0
    padding = max(half_span * 0.08, 0.5)
    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    extent = half_span + padding

    ax.set_xlim(center_x - extent, center_x + extent)
    ax.set_ylim(center_y - extent, center_y + extent)
    ax.set_aspect("equal", adjustable="box")
    ax.set_box_aspect(1)


def _render_umap_pdf(payload: Dict[str, List[Dict]], mode: str) -> io.BytesIO:
    _configure_matplotlib_pdf_style()
    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    if mode == "frequency":
        query_points = payload.get("query", []) or []
        sample_field = str(payload.get("sample_field") or "sample").strip() or "sample"
        frame = pd.DataFrame(query_points)
        if frame.empty or "sample" not in frame.columns or not frame["sample"].astype(str).str.strip().any():
            ax.text(0.5, 0.5, "Sample labels were not available for this job.", ha="center", va="center", transform=ax.transAxes, color="#64748b")
            ax.axis("off")
        else:
            frame["population"] = frame["population"].astype(str)
            frame["sample"] = frame["sample"].astype(str)
            counts = (
                frame.groupby(["sample", "population"], observed=False)
                .size()
                .rename("count")
                .reset_index()
            )
            sample_totals = counts.groupby("sample", observed=False)["count"].sum().rename("sample_total")
            counts = counts.merge(sample_totals, on="sample", how="left")
            counts["fraction"] = counts["count"] / counts["sample_total"].replace(0, np.nan)
            ranked = (
                counts.groupby("population", observed=False)["fraction"]
                .mean()
                .sort_values(ascending=False)
                .index.tolist()
            )
            pivot = (
                counts.pivot_table(index="population", columns="sample", values="fraction", aggfunc="first")
                .reindex(ranked)
                .fillna(0.0)
            )
            positions = np.arange(len(pivot.index), dtype=float)
            left = np.zeros(len(pivot.index), dtype=float)
            samples = list(pivot.columns)
            palette = _build_preview_palette(samples)
            for sample in samples:
                values = pivot[sample].to_numpy(dtype=float)
                ax.barh(
                    positions,
                    values,
                    left=left,
                    color=palette.get(sample, matplotlib.colors.to_rgb("#94a3b8")),
                    edgecolor="none",
                    label=sample,
                    height=0.8,
                )
                left = left + values
            ax.set_yticks(positions)
            ax.set_yticklabels(pivot.index.tolist())
            ax.invert_yaxis()
            ax.set_xlim(0, 1)
            ax.set_xlabel(f"Fraction of filtered cells per {sample_field}")
            ax.set_title("Cell frequency")
            if samples:
                ax.legend(frameon=False, fontsize=8, loc="lower right")
    elif mode == "cluster":
        if payload["reference"]:
            ax.scatter(
                [p["x"] for p in payload["reference"]],
                [p["y"] for p in payload["reference"]],
                s=1,
                c="#e5e7eb",
                alpha=0.3,
                linewidths=0,
            )
        populations = list(dict.fromkeys(p["population"] for p in payload["query"] if p.get("population")))
        palette = _build_preview_palette(populations)
        for population in populations:
            subset = [p for p in payload["query"] if p["population"] == population]
            ax.scatter(
                [p["x"] for p in subset],
                [p["y"] for p in subset],
                s=4,
                color=[palette.get(population, matplotlib.colors.to_rgb("#64748b"))],
                alpha=0.5,
                linewidths=0,
            )
            if subset:
                xs = [float(point["x"]) for point in subset]
                ys = [float(point["y"]) for point in subset]
                ax.text(
                    float(np.median(xs)),
                    float(np.median(ys)),
                    population,
                    fontsize=9,
                    color="#0f172a",
                    ha="center",
                    va="center",
                )
        ax.set_title("UMAP cell types")
    else:
        if payload["reference"]:
            ax.scatter(
                [p["x"] for p in payload["reference"]],
                [p["y"] for p in payload["reference"]],
                s=4,
                c="#94a3b8",
                linewidths=0,
                label="Reference",
            )
        if payload["query"]:
            ax.scatter(
                [p["x"] for p in payload["query"]],
                [p["y"] for p in payload["query"]],
                s=4,
                c="#f97316",
                linewidths=0,
                label="Query",
            )
        ax.legend(frameon=False)
        ax.set_title("UMAP broad")
    if mode != "frequency":
        _square_umap_axes(ax, [payload["reference"], payload["query"]])
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def _render_expression_pdf(payload: Dict, mode: str) -> io.BytesIO:
    _configure_matplotlib_pdf_style()
    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    gene = payload["gene"]
    modality = _normalize_modality_id(payload.get("modality"), default="rna")
    measurement = {"rna": "expression", "grn_tf": "imputed TF activity", "grn": "predicted edge score"}.get(modality, "abundance")
    violin_data = [entry for entry in payload.get("violin", []) if len(entry.get("values", []))]
    if (mode == "violin" and not violin_data) or (mode != "violin" and not payload.get("umap")):
        ax.text(.5, .5, payload.get("message") or "No observations match the selected filters.",
                ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        buf = io.BytesIO()
        fig.savefig(buf, format="pdf", bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return buf
    global_min = float(payload.get("global_min", 0.0) or 0.0)
    global_max = float(payload.get("global_max", 0.0) or 0.0)
    if global_max <= global_min:
        global_max = global_min + 1e-9
    if mode == "violin":
        positions = np.arange(1, len(violin_data) + 1)
        parts = ax.violinplot(
            [entry["values"] for entry in violin_data],
            positions=positions,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )
        for body in parts["bodies"]:
            body.set_facecolor("#94a3b8")
            body.set_edgecolor("#475569")
            body.set_alpha(0.55)
        for idx, entry in enumerate(violin_data, start=1):
            vals = np.asarray(entry["values"], dtype=float)
            jitter = np.random.default_rng(0).normal(0, 0.035, size=len(vals))
            ax.scatter(np.full(len(vals), idx) + jitter, vals, s=4, c="#0f172a", alpha=0.35, linewidths=0)
        ax.set_xticks(positions)
        ax.set_xticklabels([entry["population"] for entry in violin_data], rotation=45, ha="right")
        ax.set_title(f"{gene} {measurement} (top {len(violin_data)} states by mean)")
        ax.set_ylabel(measurement.capitalize())
        pad = max((global_max - global_min) * 0.04, 0.05)
        ax.set_ylim(global_min - pad, global_max + pad)
    else:
        umap_points = payload["umap"]
        if _normalize_modality_id(payload.get("modality"), default="rna") in {"lipids", "adt", "metabolite", "lipid", "grn", "grn_tf"}:
            expression_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "expression_blue_yellow_red",
                [
                    (0.0, "#2563eb"),
                    (0.5, "#fde047"),
                    (1.0, "#dc2626"),
                ],
            )
        else:
            expression_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                "expression_grey_red",
                [
                    (0.0, "#f3f4f6"),
                    (0.15, "#fecaca"),
                    (0.35, "#fca5a5"),
                    (0.6, "#ef4444"),
                    (1.0, "#b91c1c"),
                ],
            )
        zero_points = [p for p in umap_points if p["value"] == 0]
        measured_points = [p for p in umap_points if p["value"] != 0]
        if zero_points:
            ax.scatter([p["x"] for p in zero_points], [p["y"] for p in zero_points],
                       s=4, c="#e5e7eb", linewidths=0)
        if measured_points:
            sc = ax.scatter(
                [p["x"] for p in measured_points], [p["y"] for p in measured_points],
                s=4, c=[p["value"] for p in measured_points], cmap=expression_cmap,
                vmin=global_min, vmax=global_max, linewidths=0)
            cbar = fig.colorbar(sc, ax=ax)
            cbar.set_label(gene)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title(f"{gene} {measurement}")
        _square_umap_axes(ax, [umap_points])
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def _render_differential_gene_pdf(payload: Dict) -> io.BytesIO:
    _configure_matplotlib_pdf_style()
    fig, ax = plt.subplots(figsize=(5.0, 5.4))

    groups = payload.get("groups", [])
    positions = np.arange(1, len(groups) + 1)
    if len(groups):
        parts = ax.violinplot(
            [np.asarray(group.get("values", []), dtype=float) for group in groups],
            positions=positions,
            showmeans=False,
            showmedians=True,
            showextrema=False,
        )
        colors = ["#dc2626", "#2563eb"]
        for idx, body in enumerate(parts["bodies"]):
            color = colors[idx % len(colors)]
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.35)
        for idx, group in enumerate(groups, start=1):
            vals = np.asarray(group.get("values", []), dtype=float)
            if vals.size == 0:
                continue
            color = colors[(idx - 1) % len(colors)]
            jitter = np.random.default_rng(0).normal(0, 0.035, size=len(vals))
            ax.scatter(np.full(len(vals), idx) + jitter, vals, s=8, c=color, alpha=0.3, linewidths=0)

    ax.set_xticks(positions)
    ax.set_xticklabels([str(group.get("label", f"Group {idx}")) for idx, group in enumerate(groups, start=1)])
    ax.set_ylabel(payload.get("value_label") or "Normalized expression")
    ax.set_title(f"{payload.get('gene', 'Gene')} in {payload.get('population', 'population')}")
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def _render_differential_heatmap_pdf(payload: Dict) -> io.BytesIO:
    _configure_matplotlib_pdf_style()
    rows = payload.get("rows", []) or []
    columns = payload.get("columns", []) or []
    if not rows or not columns:
        raise HTTPException(status_code=404, detail="No differential heatmap rows were found.")

    matrix = np.asarray([row.get("values", []) for row in rows], dtype=float)
    finite = matrix[np.isfinite(matrix)]
    max_abs = float(np.max(np.abs(finite))) if finite.size else 1.0
    color_extent = max(max_abs / 3.0, 1.0)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "differential_heatmap",
        ["#00f0ff", "#000000", "#ffff00"],
    )
    cmap.set_bad("#000000")
    figure_height = max(6.0, min(28.0, len(rows) * 0.18 + 2.0))
    figure_width = max(7.0, min(18.0, len(columns) * 0.65 + 2.5))
    fig, ax = plt.subplots(figsize=(figure_width, figure_height))
    image = ax.imshow(np.ma.masked_invalid(matrix), aspect="auto", cmap=cmap, vmin=-color_extent, vmax=color_extent)
    ax.set_title(f"Heatmap: {payload.get('population', '')}")
    ax.set_xticks(np.arange(len(columns), dtype=float))
    ax.set_xticklabels([str(value) for value in columns], rotation=40, ha="right")
    y_labels = [str(row.get("gene", "")) for row in rows]
    tick_step = max(1, int(np.ceil(len(y_labels) / 80)))
    tick_positions = np.arange(0, len(y_labels), tick_step, dtype=int)
    ax.set_yticks(tick_positions)
    ax.set_yticklabels([y_labels[index] for index in tick_positions], fontsize=8)
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("log2FC")
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


#: The two bar colours of the DEG-count chart, read from the reference figure.
DEG_COUNT_UP_COLOR = "#C75252"
DEG_COUNT_DOWN_COLOR = "#7CC7E9"


def _symmetric_count_ticks(limit: int) -> List[float]:
    """Ticks either side of zero for a count axis: about four per side, 1/2/5 x 10^k."""
    target = max(1.0, float(limit) / 4.0)
    exponent = math.floor(math.log10(target))
    for multiple in (1.0, 2.0, 5.0, 10.0):
        step = multiple * (10.0 ** exponent)
        if step >= target:
            break
    count = int(math.ceil(limit / step))
    return [index * step for index in range(-count, count + 1)]


def _render_differential_summary_pdf(payload: Dict) -> io.BytesIO:
    """The DEG-count chart as vector PDF: one row per cell state, down left, up right."""
    _configure_matplotlib_pdf_style()
    rows = payload.get("rows", []) or []
    if not rows:
        raise HTTPException(status_code=404, detail="No differential counts were found.")

    labels = [str(row.get("population", "")) for row in rows]
    up = [int(row.get("up", 0) or 0) for row in rows]
    down = [int(row.get("down", 0) or 0) for row in rows]
    # Top of the axis is the first cell state of the lineage order, as the chart reads.
    positions = np.arange(len(rows), dtype=float)[::-1]

    height = max(3.0, min(24.0, len(rows) * 0.30 + 1.6))
    fig, ax = plt.subplots(figsize=(6.5, height))
    ax.barh(positions, up, height=0.62, color=DEG_COUNT_UP_COLOR, label="Upregulated",
            edgecolor="none")
    ax.barh(positions, [-value for value in down], height=0.62, color=DEG_COUNT_DOWN_COLOR,
            label="Downregulated", edgecolor="none")
    ax.axvline(0.0, color="#000000", linewidth=0.8)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    limit = max(1, int(payload.get("max_count", 0) or 0))
    # Ticks first, limits last: set_xticks widens the axis to hold every tick it is
    # given, so setting the limit first and the ticks after left the axis at +/-800
    # when the largest count was 557.
    ticks = _symmetric_count_ticks(limit)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(int(abs(value))) for value in ticks])
    ax.set_xlim(-limit * 1.08, limit * 1.08)
    feature_label = str(payload.get("feature_label") or "gene")
    ax.set_xlabel(f"Number of differential {feature_label}s")
    case_label = str(payload.get("case_label") or "case")
    control_label = str(payload.get("control_label") or "control")
    ax.set_title(f"{case_label} versus {control_label}")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def _render_differential_volcano_pdf(payload: Dict) -> io.BytesIO:
    _configure_matplotlib_pdf_style()
    points = payload.get("points", []) or []
    if not points:
        raise HTTPException(status_code=404, detail="No differential volcano points were found.")

    up = [point for point in points if str(point.get("direction", "")).lower() == "up"]
    down = [point for point in points if str(point.get("direction", "")).lower() == "down"]

    fig, ax = plt.subplots(figsize=(8.0, 6.5))
    if down:
        ax.scatter(
            [float(point.get("log2fc", 0.0)) for point in down],
            [float(point.get("score", 0.0)) for point in down],
            s=26,
            c="#2563eb",
            alpha=0.72,
            linewidths=0,
            label="Down",
        )
    if up:
        ax.scatter(
            [float(point.get("log2fc", 0.0)) for point in up],
            [float(point.get("score", 0.0)) for point in up],
            s=26,
            c="#dc2626",
            alpha=0.72,
            linewidths=0,
            label="Up",
        )
    ax.set_title(f"Volcano: {payload.get('population', '')}")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel(f"-log10({payload.get('statistic_label', 'FDR')})")
    ax.axvline(0.0, color="#94a3b8", linewidth=0.8, alpha=0.6)
    if up or down:
        ax.legend(frameon=False)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def _render_differential_go_pdf(payload: Dict) -> io.BytesIO:
    _configure_matplotlib_pdf_style()
    terms = payload.get("terms", []) or []
    if not terms:
        raise HTTPException(status_code=404, detail="No differential GO terms were available.")

    positive_terms = []
    background_terms = []
    x_values = []
    y_values = []
    for term in terms:
        z_score = float(term.get("z_score", np.nan) or np.nan)
        fdr_plot = float(term.get("fdr_plot", np.nan) or np.nan)
        if not (_is_finite_number(z_score) and _is_finite_number(fdr_plot) and fdr_plot > 0):
            continue
        x_values.append(z_score)
        y_values.append(fdr_plot)
        (positive_terms if bool(term.get("is_selected_positive_sig")) else background_terms).append((z_score, fdr_plot, term))
    if not x_values or not y_values:
        raise HTTPException(status_code=404, detail="No differential GO terms were available.")

    fig, ax = plt.subplots(figsize=(9.2, 7.8))
    if background_terms:
        ax.scatter(
            [entry[0] for entry in background_terms],
            [entry[1] for entry in background_terms],
            s=42,
            c="#d1d5db",
            alpha=0.95,
            linewidths=0,
            zorder=2,
        )
    if positive_terms:
        ax.scatter(
            [entry[0] for entry in positive_terms],
            [entry[1] for entry in positive_terms],
            s=48,
            c="#1f19c7",
            alpha=0.98,
            linewidths=0,
            zorder=3,
        )

    labels = payload.get("labels", []) or []
    annotation_offsets = [(-6, 38), (-2, 10), (2, -10), (6, -34), (18, -60)]
    for index, label in enumerate(labels):
        z_score = float(label.get("z_score", np.nan) or np.nan)
        fdr_plot = float(label.get("fdr_plot", np.nan) or np.nan)
        if not (_is_finite_number(z_score) and _is_finite_number(fdr_plot) and fdr_plot > 0):
            continue
        dx, dy = annotation_offsets[min(index, len(annotation_offsets) - 1)]
        ax.annotate(
            str(label.get("term_name", "")),
            xy=(z_score, fdr_plot),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=9,
            color=str(label.get("label_color", "#111827")),
            arrowprops={
                "arrowstyle": "-",
                "color": str(label.get("label_color", "#111827")),
                "linewidth": 1.0,
                "alpha": 0.9,
                "shrinkA": 0,
                "shrinkB": 0,
            },
            zorder=4,
        )

    x_min = min(-10.0, float(np.floor(min(x_values) - 0.5)))
    x_max = max(20.0, float(np.ceil(max(x_values) + 2.5)))
    y_min = max(min(y_values) * 0.5, 1e-300)
    ax.set_yscale("log")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, 1.0)
    ax.set_xlabel("Z-Score")
    ax.set_ylabel("Fishers FDR p")
    ax.set_title(f"GO terms: {payload.get('population', '')}")
    ax.axvline(0.0, color="#111827", linewidth=1.2, alpha=0.95, zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(False)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def _network_positions(node_ids: List[str]) -> Dict[str, tuple[float, float]]:
    ordered = list(node_ids)
    positions: Dict[str, tuple[float, float]] = {}
    if not ordered:
        return positions
    if len(ordered) == 1:
        positions[ordered[0]] = (0.0, 0.0)
        return positions

    positions[ordered[0]] = (0.0, 0.0)
    remaining = ordered[1:]
    placed = 0
    ring_index = 1
    while placed < len(remaining):
        ring_capacity = max(8, ring_index * 10)
        ring_nodes = remaining[placed:placed + ring_capacity]
        radius = 1.3 * ring_index
        for index, node_id in enumerate(ring_nodes):
            angle = (2.0 * np.pi * index) / max(len(ring_nodes), 1)
            positions[node_id] = (radius * float(np.cos(angle)), radius * float(np.sin(angle)))
        placed += len(ring_nodes)
        ring_index += 1
    return positions


def _render_network_pdf(payload: Dict, title: str) -> io.BytesIO:
    _configure_matplotlib_pdf_style()
    elements = payload.get("elements", []) or []
    node_map: Dict[str, Dict[str, Any]] = {}
    preset_positions: Dict[str, tuple[float, float]] = {}
    edges: List[Dict[str, Any]] = []
    for element in elements:
        data = element.get("data") or {}
        if not data:
            continue
        if data.get("source") and data.get("target"):
            edges.append(data)
            continue
        node_id = str(data.get("id", "")).strip()
        if node_id:
            node_map[node_id] = data
            position = element.get("position") or {}
            x = position.get("x")
            y = position.get("y")
            if _is_finite_number(x) and _is_finite_number(y):
                preset_positions[node_id] = (float(x), float(y))
    if not node_map:
        raise HTTPException(status_code=404, detail="No interaction network was available.")

    ordered_nodes = [
        node_id
        for _, node_id in sorted(
            (
                (abs(float((data.get("log2fc") or 0.0))), str(node_id))
                for node_id, data in node_map.items()
            ),
            reverse=True,
        )
    ]
    positions = preset_positions if len(preset_positions) == len(node_map) else _network_positions(ordered_nodes)
    fig, ax = plt.subplots(figsize=(8.2, 7.4))
    for edge in edges:
        source = str(edge.get("source", "")).strip()
        target = str(edge.get("target", "")).strip()
        if source not in positions or target not in positions:
            continue
        interaction_type = str(edge.get("interaction_type", "")).lower()
        color = str(edge.get("edge_color", "")).strip() or "#9ca3af"
        if "transcription" in interaction_type:
            color = "#ef4444"
        elif "tbar" in interaction_type:
            color = "#60a5fa"
        elif interaction_type == "cell_communication_diff" and str(edge.get("edge_color", "")).strip():
            color = str(edge.get("edge_color", "")).strip()
        start = positions[source]
        end = positions[target]
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops={
                "arrowstyle": "-|>",
                "color": color,
                "linewidth": max(1.0, float(edge.get("weight", 3.0) or 3.0) / 3.2),
                "shrinkA": 14,
                "shrinkB": 14,
                "alpha": min(1.0, max(0.2, float(edge.get("edge_opacity", 0.8) or 0.8))),
            },
        )

    for node_id in ordered_nodes:
        data = node_map[node_id]
        x, y = positions[node_id]
        log2fc = float(data.get("log2fc", 0.0) or 0.0)
        color = str(data.get("color", "")).strip() or ("#cbd5e1" if not _is_finite_number(data.get("log2fc")) or log2fc == 0 else "#fca5a5" if log2fc > 0 else "#7dd3fc")
        node_size = 460 if str(data.get("node_type", "")).strip() == "focus" else 320
        ax.scatter([x], [y], s=node_size, c=[color], edgecolors="white", linewidths=1.0, zorder=3)
        ax.text(x, y, str(data.get("label") or node_id), ha="center", va="center", fontsize=9, zorder=4)

    coords = np.asarray(list(positions.values()), dtype=float)
    if coords.size:
        span = max(float(coords[:, 0].max() - coords[:, 0].min()), float(coords[:, 1].max() - coords[:, 1].min()), 1.0)
        padding = span * 0.2
        ax.set_xlim(float(coords[:, 0].min()) - padding, float(coords[:, 0].max()) + padding)
        ax.set_ylim(float(coords[:, 1].min()) - padding, float(coords[:, 1].max()) + padding)
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def _render_marker_heatmap_pdf(frame: pd.DataFrame, title: str = "Marker heatmap") -> io.BytesIO:
    _configure_matplotlib_pdf_style()
    if frame.empty:
        raise HTTPException(status_code=404, detail="Marker heatmap matrix unavailable.")
    values = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    matrix = values.to_numpy(dtype=float)
    finite = matrix[np.isfinite(matrix)]
    max_abs = float(np.max(np.abs(finite))) if finite.size else 2.0
    color_extent = max(2.0, max_abs)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "marker_heatmap",
        ["#00f0ff", "#000000", "#ffff00"],
    )

    fig_height = max(5.5, min(24.0, matrix.shape[0] * 0.12 + 1.8))
    fig_width = max(7.5, min(24.0, matrix.shape[1] * 0.08 + 2.6))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-color_extent, vmax=color_extent)
    ax.set_title(title)

    row_labels = [str(index) for index in values.index]
    row_step = max(1, int(np.ceil(len(row_labels) / 80)))
    row_positions = np.arange(0, len(row_labels), row_step, dtype=int)
    ax.set_yticks(row_positions)
    ax.set_yticklabels([row_labels[index] for index in row_positions], fontsize=8)

    col_labels = [str(label) for label in values.columns]
    if len(col_labels) <= 60:
        ax.set_xticks(np.arange(len(col_labels), dtype=float))
        ax.set_xticklabels(col_labels, rotation=90, fontsize=6)
    else:
        col_step = max(1, int(np.ceil(len(col_labels) / 40)))
        col_positions = np.arange(0, len(col_labels), col_step, dtype=int)
        ax.set_xticks(col_positions)
        ax.set_xticklabels([col_labels[index] for index in col_positions], rotation=90, fontsize=6)

    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("Fold")
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf


def create_app(test_config: dict | None = None) -> FastAPI:
    cfg = load_config(test_config)
    cfg["ROOT_PATH"] = _normalize_root_path(cfg.get("ROOT_PATH"))
    template_dir = Path(str(cfg.get("TEMPLATE_DIR") or (BASE_DIR / "templates")))
    static_dir = Path(str(cfg.get("STATIC_DIR") or (BASE_DIR / "static")))
    index_template = str(cfg.get("INDEX_TEMPLATE") or "index.html")
    templates = Jinja2Templates(directory=str(template_dir))
    app = FastAPI(title=cfg["APP_TITLE"], root_path=cfg["ROOT_PATH"])
    app.state.config = cfg
    app.state.root_path = cfg["ROOT_PATH"]
    app.state.job_store = JobStore(Path(cfg["JOB_STORAGE"]))
    app.state.expression_cache = {}
    app.state.expression_cache_locks = {}
    app.state.differential_cache = {}
    app.state.marker_heatmap_cache = {}
    app.state.fastcomm_cache = {}
    app.state.reference_adata_cache = {}
    app.state.reference_adata_cache_locks = {}
    app.state.cache_registry_lock = threading.Lock()
    app.state.job_runner = JobRunner(
        app.state.job_store,
        Path(cfg["REFERENCE_REGISTRY"]),
        max_workers=cfg["JOB_WORKERS"],
        export_approx_pdfs=cfg.get("EXPORT_APPROX_PDFS", False),
        h5ad_compression=cfg.get("H5AD_COMPRESSION", "lzf"),
    )

    @app.exception_handler(RequestValidationError)
    async def handle_validation_error(request: Request, exc: RequestValidationError):
        if _is_api_request(request):
            return JSONResponse(status_code=422, content={"detail": str(exc)})
        raise exc

    @app.exception_handler(Exception)
    async def handle_unexpected_error(request: Request, exc: Exception):
        if _is_api_request(request):
            logging.exception("Unhandled API error for %s", request.url.path, exc_info=exc)
            return JSONResponse(status_code=500, content={"detail": f"{type(exc).__name__}: {exc}" if str(exc) else type(exc).__name__})
        raise exc

    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

    @app.get("/", response_class=HTMLResponse)
    async def index(request: Request):
        registry = _load_reference_registry(app)
        css_version = int((static_dir / "styles.css").stat().st_mtime) if (static_dir / "styles.css").exists() else 0
        js_version = max((path.stat().st_mtime_ns for path in
                          (static_dir / "app.js", static_dir / "integrated.js")
                          if path.exists()), default=0)
        return templates.TemplateResponse(
            request,
            index_template,
            {
                "request": request,
                "registry_json": json.dumps(registry),
                "app_title": cfg["APP_TITLE"],
                "app_root_path": cfg["ROOT_PATH"],
                "styles_version": css_version,
                "app_js_version": js_version,
            },
        )

    @app.get("/api/meta/species")
    async def meta_species():
        return JSONResponse(_load_reference_registry(app))

    @app.get("/api/meta/reference-preview")
    def reference_preview(species: str = Query(...), reference: str = Query(...)):
        return JSONResponse(_build_reference_preview_payload(app, species, reference))

    @app.post("/api/jobs")
    async def create_job(
        species: str = Form(...),
        reference: str = Form(...),
        ambient_option: str | None = Form(None),
        soupx_option: str | None = Form(None),
        sample_names: List[str] = Form(...),
        files: List[UploadFile] = File(...),
    ):
        store, _ = _job_resources(app)
        if len(files) > cfg["MAX_FILES_PER_JOB"]:
            raise HTTPException(status_code=400, detail=f"Maximum {cfg['MAX_FILES_PER_JOB']} files are allowed per job.")
        if len(sample_names) != len(files):
            raise HTTPException(status_code=400, detail="Each uploaded file must include a matching sample name.")

        normalized_samples: List[str] = []
        for name in sample_names:
            clean_name = name.strip()
            if not clean_name:
                raise HTTPException(status_code=400, detail="Sample names cannot be empty.")
            if clean_name in normalized_samples:
                raise HTTPException(status_code=400, detail=f"Duplicate sample name '{clean_name}' detected.")
            normalized_samples.append(clean_name)

        effective_ambient_option = ambient_option if ambient_option is not None else soupx_option
        metadata = store.create_job(species, reference, effective_ambient_option, files=[])
        job_id = metadata["job_id"]
        _invalidate_expression_cache(app, job_id)
        _invalidate_differential_cache(app, job_id)
        _invalidate_marker_heatmap_cache(app, job_id)
        _invalidate_fastcomm_cache(app, job_id)
        uploads_dir = store.uploads_dir(job_id)
        records = []
        used_names: set[str] = set()
        for sample, upload in zip(normalized_samples, files):
            if not upload.filename:
                raise HTTPException(status_code=400, detail="One uploaded file is missing a filename.")
            if not _allowed_file(app, upload.filename):
                raise HTTPException(status_code=400, detail=f"Unsupported file extension for {upload.filename}.")
            # The sample name reaches the filesystem, so it is sanitised the same
            # way the uploaded filename is. Without this a name like "../../x"
            # writes outside the job's uploads directory. The unsanitised value is
            # still kept as sample_name, which is only ever displayed.
            dest_name = f"{_secure_filename(sample)}_{_secure_filename(upload.filename)}"
            if dest_name in used_names:
                raise HTTPException(status_code=400, detail=f"Duplicate sample name '{sample}' detected.")
            used_names.add(dest_name)
            dest_path = _require_within(uploads_dir / dest_name, uploads_dir, "Upload destination")
            content = await upload.read()
            dest_path.write_bytes(content)
            records.append({"sample_name": sample, "filename": dest_name, "size": dest_path.stat().st_size})

        store.update_job(job_id, files=records, message="Upload complete. Configure QC to proceed.")
        return JSONResponse({"job_id": job_id, "status": "uploaded"})

    @app.post("/api/jobs/{job_id}/qc")
    async def update_qc(job_id: str, qc: QCSettings):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        registry_path = Path(app.state.config["REFERENCE_REGISTRY"])
        reference_entry = pipeline_mod._lookup_reference(meta["species"], meta["reference"], registry_path)
        supported_modalities = {
            _normalize_modality_id(value, default="")
            for value in (reference_entry.get("impute_modalities") or [])
        }
        if qc.impute_modalities is not None:
            requested_raw = qc.impute_modalities
        elif qc.impute_modality and qc.impute_modality != "none":
            requested_raw = [qc.impute_modality]
        else:
            requested_raw = []
        requested = [m for m in (_normalize_modality_id(v, default="") for v in requested_raw) if m]
        # "all" selects every supported modality (expanded at run time); not a literal modality
        unsupported = [m for m in requested if m != "all" and m not in supported_modalities]
        if unsupported:
            raise HTTPException(status_code=400,
                                detail=f"Impute modalities not supported for this reference: {', '.join(unsupported)}")
        meta = store.update_job(job_id, qc=qc.model_dump(), message="QC parameters saved.")
        return JSONResponse({"job_id": job_id, "qc": meta["qc"]})

    @app.post("/api/jobs/{job_id}/configure")
    async def configure_job(job_id: str, config: JobConfigSettings):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")

        registry_path = Path(app.state.config["REFERENCE_REGISTRY"])
        reference_entry = pipeline_mod._lookup_reference(config.species, config.reference, registry_path)
        pipeline_mod._ensure_reference_fields(reference_entry)

        meta = store.get_job(job_id)
        requested_ambient = config.ambient_option if config.ambient_option is not None else meta.get("ambient_option", meta.get("soupx_option"))
        changed = (
            str(meta.get("species") or "") != config.species
            or str(meta.get("reference") or "") != config.reference
            or str(meta.get("ambient_option", meta.get("soupx_option")) or "") != str(requested_ambient or "")
        )
        if changed:
            _invalidate_expression_cache(app, job_id)
            _invalidate_differential_cache(app, job_id)
            _invalidate_marker_heatmap_cache(app, job_id)
            _clear_directory_contents(store.outputs_dir(job_id))
            _clear_directory_contents(store.logs_dir(job_id))
            meta = store.update_job(
                job_id,
                species=config.species,
                reference=config.reference,
                ambient_option=requested_ambient,
                status="uploaded",
                progress=0,
                message="Reference updated. Configure QC and rerun alignment.",
                artifacts={},
                marker_analysis={},
                marker_analysis_by_modality={},
                fastcomm_analysis={},
                modality_artifacts={},
                modalities={"default": "rna", "available": [dict(_DEFAULT_MODALITY_DEFINITIONS["rna"])]},
                differential={},
            )
        return JSONResponse(
            {
                "job_id": job_id,
                "status": meta.get("status"),
                "species": meta.get("species"),
                "reference": meta.get("reference"),
                "ambient_option": meta.get("ambient_option", meta.get("soupx_option")),
                "changed": changed,
            }
        )

    @app.post("/api/jobs/{job_id}/run")
    async def run_job(job_id: str):
        store, runner = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        _invalidate_expression_cache(app, job_id)
        _invalidate_differential_cache(app, job_id)
        _invalidate_marker_heatmap_cache(app, job_id)
        _invalidate_fastcomm_cache(app, job_id)
        store.update_job(
            job_id,
            marker_analysis={},
            marker_analysis_by_modality={},
            fastcomm_analysis={},
            modality_artifacts={},
            modalities={"default": "rna", "available": [dict(_DEFAULT_MODALITY_DEFINITIONS["rna"])]},
            differential={},
        )
        runner.submit(job_id)
        store.append_log(job_id, "Job queued by user request.")
        store.update_job(job_id, message="Job submitted to worker.", status="queued", progress=15)
        return JSONResponse({"job_id": job_id, "status": "queued"})

    @app.get("/api/jobs/{job_id}/status")
    async def status(job_id: str):
        store, runner = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = _job_metadata_with_recovery(store, runner, job_id)
        log_path = store.logs_dir(job_id) / "pipeline.log"
        log_head: List[str] = []
        log_tail: List[str] = []
        if log_path.exists():
            all_lines = log_path.read_text(encoding="utf-8").splitlines(True)
            log_head = all_lines[:80]
            log_tail = all_lines[-200:]
        meta["log_head"] = log_head
        meta["log_tail"] = log_tail
        meta["qc_log_tail"] = log_tail if meta.get("status") == "failed" else _filter_qc_log_lines(log_tail)
        meta["message"] = _derive_live_pipeline_message(meta.get("status"), log_tail, meta.get("message"))
        if isinstance(meta.get("fastcomm_analysis"), dict) and meta["fastcomm_analysis"].get("enabled"):
            if not meta["fastcomm_analysis"].get("populations"):
                meta["fastcomm_analysis"] = {
                    **meta["fastcomm_analysis"],
                    "populations": _fastcomm_populations(app, meta),
                }
        meta["differential_ui"] = _build_differential_payload(app, job_id, meta, root_path=app.state.root_path)
        return JSONResponse(meta, headers={"Cache-Control": "no-store"})

    @app.get("/api/jobs/{job_id}/differential/status")
    async def differential_status(job_id: str):
        store, runner = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = _job_metadata_with_recovery(store, runner, job_id)
        return JSONResponse(
            _build_differential_payload(app, job_id, meta, root_path=app.state.root_path),
            headers={"Cache-Control": "no-store"},
        )

    @app.post("/api/jobs/{job_id}/differential")
    async def run_differential(job_id: str, payload: DifferentialSettings):
        store, runner = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = _job_metadata_with_recovery(store, runner, job_id)
        if (meta.get("differential") or {}).get("status") in {"queued", "processing"}:
            raise HTTPException(409, "A differential analysis is already running.")
        history = completed_differentials(meta)
        if history:
            store.update_job(job_id, differential_history=history)
        config = _validate_differential_request(meta, payload)
        options = _differential_options(meta)
        _invalidate_differential_cache(app, job_id)
        store.update_job(
            job_id,
            differential={
                "status": "queued",
                "worker_pid": os.getpid(),
                "progress": 5,
                "message": "Differential analysis queued.",
                "config": config,
                "artifacts": {},
                "networks": [],
                "go_terms_included": False,
                "default_population_col": options.get("default_population_col"),
            },
        )
        store.append_log(job_id, "Differential analysis queued by user request.")
        runner.submit_differential(job_id)
        updated = store.get_job(job_id)
        return JSONResponse(_build_differential_payload(app, job_id, updated, root_path=app.state.root_path))

    @app.get("/api/jobs/{job_id}/umap")
    def umap(
        job_id: str,
        modality: str = Query("rna"),
        filter1_field: Optional[str] = Query(None),
        filter1_values: List[str] = Query([]),
        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
        color_by: str = Query(""),
        coords: str = Query(""),
        x_field: str = Query(""),
        y_field: str = Query(""),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
        try:
            return JSONResponse(_build_umap_payload(
                app, meta, modality=modality, display_filters=display_filters,
                color_by=color_by, coords_key=coords,
                x_field=x_field, y_field=y_field))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/expression")
    async def expression(
        job_id: str,
        gene: str = Query(...),
        modality: str = Query("rna"),
        filter1_field: Optional[str] = Query(None),
        filter1_values: List[str] = Query([]),
        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
        violin_limit: int = Query(10),
        x_field: str = Query(""),
        y_field: str = Query(""),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
        try:
            return JSONResponse(_build_expression_payload(
                app, meta, gene, modality=modality, display_filters=display_filters,
                x_field=str(x_field or ""), y_field=str(y_field or ""),
                # Coerced defensively: this route is also called directly by
                # the scALABLE viewer's wrapper, where an unpassed argument
                # arrives as a FastAPI Query object rather than an int.
                violin_limit=_as_int(violin_limit, 10, 1, 80)))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/marker/network")
    def marker_network(job_id: str, population: str = Query(...), modality: str = Query("rna")):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        try:
            return JSONResponse(_build_marker_network_payload(meta, population, modality=modality))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/grn/regulator-network")
    def regulator_network(job_id: str, cell_state: str = Query(""), contrast: str = Query(""),
                          features: str = Query(""), limit: int = Query(200, ge=1, le=2000),
                          edge_percentile: float = Query(50, ge=0, le=100),
                          expression_percentile: float = Query(50, ge=0, le=100)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(404, "Job not found.")
        ds = UploadedGrnData(app, store.get_job(job_id))
        if contrast and contrast not in ds.runs:
            raise HTTPException(404, "Unknown completed comparison.")
        return gnet.regulator_network(ds, cell_state, contrast=contrast or ds.current_contrast,
                                      features=features.replace(",", " ").split() or None,
                                      limit=limit, edge_percentile=edge_percentile,
                                      expression_percentile=expression_percentile)

    @app.get("/api/jobs/{job_id}/grn/tf-activity")
    def tf_activity(job_id: str, cell_state: str = Query(""), contrast: str = Query(""),
                    factors: str = Query(""), limit: int = Query(25, ge=1, le=2000)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(404, "Job not found.")
        ds = UploadedGrnData(app, store.get_job(job_id))
        if contrast and contrast not in ds.runs:
            raise HTTPException(404, "Unknown completed comparison.")
        return gnet.tf_activity_profile(ds, cell_state=cell_state, contrast=contrast or ds.current_contrast,
                                       factors=factors.replace(",", " ").split() or None, limit=limit)

    @app.post("/api/jobs/{job_id}/differential/select")
    def select_completed_differential(job_id: str, contrast: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(404, "Job not found.")
        meta = store.get_job(job_id)
        if (meta.get("differential") or {}).get("status") in {"queued", "processing"}:
            raise HTTPException(409, "Wait for the running comparison to finish.")
        run = completed_differentials(meta).get(contrast)
        if run is None:
            raise HTTPException(404, "Unknown completed comparison.")
        _invalidate_differential_cache(app, job_id)
        store.update_job(job_id, differential=run)
        return JSONResponse(_build_differential_payload(app, job_id, store.get_job(job_id), root_path=app.state.root_path))

    @app.get("/api/jobs/{job_id}/grn/network")
    def grn_network(
        job_id: str,
        genes: List[str] = Query(default=[]),
        sample: str = Query(default=""),
        cell_state: str = Query(default=""),
        threshold: float = Query(default=0.0, ge=0),
        max_edges: int = Query(default=25, ge=1, le=1000),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        gene_list = []
        for entry in genes:
            gene_list.extend(str(entry).replace(",", " ").split())
        try:
            return JSONResponse(_build_grn_network_payload(
                meta, gene_list, sample=sample, cell_state=cell_state,
                threshold=threshold, max_edges=max_edges))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/fastcomm/network")
    def fastcomm_network(
        job_id: str,
        population: str = Query(""),
        direction: str = Query("incoming"),
        limit: int = Query(35),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        try:
            return JSONResponse(_build_fastcomm_payload(app, meta, population, direction=direction, limit=limit))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/fastcomm/plot")
    def fastcomm_plot(
        job_id: str,
        population: str = Query(""),
        plot_type: str = Query("focused_incoming"),
        limit: int = Query(60),
        filter1_field: str = Query(""),
        filter1_values: List[str] = Query(default=[]),
        filter2_field: str = Query(""),
        filter2_values: List[str] = Query(default=[]),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        try:
            display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
            return JSONResponse(
                _build_fastcomm_plot_payload(
                    app,
                    meta,
                    population,
                    plot_type=plot_type,
                    limit=limit,
                    display_filters=display_filters,
                )
            )
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except HTTPException:
            raise
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/feature-annotations")
    def feature_annotations(job_id: str, modality: str = Query("metabolite")):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        return JSONResponse(metabolite_annotations.for_dataset(store.get_job(job_id), modality))

    @app.get("/api/jobs/{job_id}/genes")
    def job_genes(job_id: str, modality: str = Query("rna")):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        try:
            return JSONResponse(_build_gene_suggestions_payload(app, meta, modality=modality))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/dotplot")
    def dotplot(job_id: str, genes: str = Query(""), modality: str = Query("rna"),
                      group_by: str = Query(""), groups: List[str] = Query([]),
                      subset_by: str = Query(""), subset_values: List[str] = Query([]),
                      subset2_by: str = Query(""), subset2_values: List[str] = Query([])):
        """Mean expression and detected fraction per (gene, cell state).

        Colour and dot size for the DotPlot. Cell states are returned in the
        dataset's own order, which is the centroid ordering, so the figure reads
        the same way as every other plot in the tool.
        """
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        cache = _get_expression_cache(app, store.get_job(job_id), modality=modality)
        wanted = (_split_expression_features(genes, cache)
                  or _default_marker_genes(cache, group_by))
        if not wanted:
            raise HTTPException(status_code=400, detail="Give at least one gene.")
        payload = _gene_state_stats(cache, wanted, group_by, list(groups),
                                    subset_by, list(subset_values), subset2_by, list(subset2_values))
        if not payload["genes"]:
            raise HTTPException(
                status_code=404,
                detail=f"none of the {len(wanted)} requested genes are in this dataset")
        payload["modality"] = _normalize_modality_id(modality)
        return JSONResponse(payload)

    @app.get("/api/jobs/{job_id}/combplot")
    def combplot(job_id: str, genes: str = Query(""), modality: str = Query("rna"),
                       min_cells: int = Query(5), cells_per_sample: int = Query(10), unit: str = Query("cells", pattern="^(cells|donor)$"),
                       group_by: str = Query(""), groups: List[str] = Query([]),
                      subset_by: str = Query(""), subset_values: List[str] = Query([]),
                      subset2_by: str = Query(""), subset2_values: List[str] = Query([])):
        """Individual cells by default, with opt-in per-donor means."""
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        cache = _get_expression_cache(app, store.get_job(job_id), modality=modality)
        # Blank means the marker gene of every group, the same default the
        # DotPlot uses, so switching between the two keeps the same gene set.
        wanted = (_split_expression_features(genes, cache)
                  or _default_marker_genes(cache, group_by))
        if not wanted:
            raise HTTPException(status_code=400, detail="Give at least one gene.")
        if unit == "cells":
            payload = _gene_cell_values(cache, wanted, group_by, list(groups),
                                        subset_by, list(subset_values), subset2_by, list(subset2_values), cells_per_sample)
        else:
            payload = _gene_donor_state_means(cache, wanted, max(1, int(min_cells)),
                                            group_by, list(groups), subset_by, list(subset_values),
                                            subset2_by, list(subset2_values))
        if payload.get("error"):
            raise HTTPException(status_code=404, detail=payload["error"])
        payload["modality"] = _normalize_modality_id(modality)
        return JSONResponse(payload)


    @app.get("/api/jobs/{job_id}/plot-variables")
    def plot_variables(job_id: str, modality: str = Query("rna")):
        """The variables the DotPlot and CombPlot may group or filter by.

        Only categorical columns with a workable number of levels; a per-cell
        numeric column would draw one column per cell.
        """
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        try:
            cache = _get_expression_cache(app, store.get_job(job_id), modality=modality)
        except FileNotFoundError:
            # A DEG-only modality stores no feature matrix. Cell communication is one: it
            # reaches a bundle through --deg-modality alone. The variables returned here
            # are obs columns, which every modality shares, so read them from RNA instead
            # of answering 500. Before this, plot-variables?modality=cell_communication
            # raised FileNotFoundError and the server returned 500.
            cache = _get_expression_cache(app, store.get_job(job_id), modality="rna")
        variables = _groupable_columns(cache)
        return JSONResponse({"cluster_key": str(cache["cluster_key"]),
                             "variables": variables,
                             # The UMAP panel colours by any of the same columns,
                             # draws on any 2-D embedding the h5ad carries, and
                             # takes any pair of numeric obs columns as axes.
                             "color_variables": variables,
                             "coords": _umap_coordinate_options(cache),
                             "numeric_variables": _axis_field_options(cache)})

    @app.get("/api/jobs/{job_id}/display-filters")
    def job_display_filters(job_id: str, modality: str = Query("rna")):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        try:
            return JSONResponse(_build_display_filter_payload(app, meta, modality=modality))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/umap/pdf")
    def umap_pdf(
        job_id: str,
        mode: str = Query("relative"),
        modality: str = Query("rna"),
        filter1_field: Optional[str] = Query(None),
        filter1_values: List[str] = Query([]),
        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
        try:
            payload = _build_umap_payload(app, meta, modality=modality, display_filters=display_filters)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))
        pdf = _render_umap_pdf(payload, mode)
        return StreamingResponse(
            pdf,
            media_type="application/pdf",
            headers={"Content-Disposition": f'attachment; filename="{job_id}_umap.pdf"'},
        )

    @app.get("/api/jobs/{job_id}/expression/pdf")
    def expression_pdf(
        job_id: str,
        gene: str = Query(...),
        mode: str = Query("umap"),
        modality: str = Query("rna"),
        filter1_field: Optional[str] = Query(None),
        filter1_values: List[str] = Query([]),
        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
        try:
            payload = _build_expression_payload(app, meta, gene, modality=modality, display_filters=display_filters)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))
        pdf = _render_expression_pdf(metabolite_annotations.pdf_payload(payload, meta, modality), mode)
        return StreamingResponse(
            pdf,
            media_type="application/pdf",
            headers={"Content-Disposition": f'attachment; filename="{job_id}_{gene}_expression.pdf"'},
        )

    @app.post("/api/jobs/{job_id}/client-log")
    async def client_log(job_id: str, payload: ClientLogRequest):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        message = str(payload.message or "").strip()
        if message:
            store.append_log(job_id, f"[client] {message}")
        return JSONResponse({"ok": True})

    @app.get("/api/jobs/{job_id}/chat-examples")
    def chat_examples(job_id: str):
        """Example questions for THIS job, named after its own cell states.

        The Chat tab shipped eight fixed lung sentences, so a bone-marrow job
        offered AT2 and COPD examples that its data cannot answer.
        """
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        return JSONResponse(_chat_examples(app, store.get_job(job_id)))

    @app.post("/api/jobs/{job_id}/chat")
    def chat(job_id: str, payload: ChatRequest):
        """Answer one question about this job, with numbers from this job.

        Two steps, deliberately separated. The assistant reads the sentence into
        one supported query and sees no data. This route then runs that query
        against the job's own files. So the model chooses the question and the
        data answers it; no number here comes from the model.
        """
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        question = str(payload.question or "").strip()
        if not question:
            raise HTTPException(status_code=400, detail="question is required")

        from .cross_pathways import answer_if_requested as pathway_answer
        pathways = pathway_answer(app, meta, question)
        if pathways is not None:
            return JSONResponse(pathways)
        from .modality_markers import answer_if_requested as marker_answer
        markers = marker_answer(app, meta, question)
        if markers is not None:
            return JSONResponse(markers)
        from .cross_modal import answer_if_requested
        cross_answer = answer_if_requested(app, meta, question)
        if cross_answer is not None:
            return JSONResponse(cross_answer)
        cache = _get_expression_cache(app, meta)
        from .integration_chat import read_question as read_integrated, answer as integrated_answer
        integrated = read_integrated(question, [s for s, _ in _chat_states_by_size(cache)], cache.get("var_names", []))
        if integrated is not None:
            return JSONResponse(integrated_answer(app, meta, question, integrated))
        reading = _chat_read_question(question, cache, meta)
        intent = str(reading.get("intent") or "")
        state = str(reading.get("cell_state") or "")
        state2 = str(reading.get("cell_state_2") or "")
        genes = [str(g) for g in (reading.get("genes") or []) if str(g).strip()]
        states = [s for s, _ in _chat_states_by_size(cache)]
        contrast = _chat_contrast(meta)

        # The router reads the sentence with a keyword matcher that overlaps
        # state names and drops gene names. Both are repaired against the names
        # this job actually holds, and the repaired reading is returned so the
        # answer says which state and which gene it used.
        named = _chat_states_in_question(question, states)
        if named:
            if state not in named:
                state = named[0]
            if state2 and state2 not in named:
                state2 = next((s for s in named if s != state), "")
            if not state2 and len(named) > 1:
                state2 = next((s for s in named if s != state), "")
        if not genes:
            genes = _chat_genes_in_question(question, cache)
        reading = dict(reading)
        reading["cell_state"] = state
        reading["cell_state_2"] = state2
        reading["genes"] = genes
        result: Dict[str, Any] = {"question": question, "reading": reading, "intent": intent}

        if intent == "clarify" and genes and not state:
            # The router's keyword matcher returns `clarify` for a plain sentence
            # such as "Which cell states express Mki67?". The question names a
            # gene this dataset holds, so the reading is stated and answered
            # rather than bounced back as a question.
            intent = "expression"
            result["intent"] = intent
            result["status"] = "read_from_question"

        if intent == "clarify":
            result["answer"] = (
                f"I need to know {reading.get('missing') or 'a little more'}. "
                f"This dataset has {len(states)} cell states and "
                f"{1 if contrast else 0} completed comparison(s).")
            result["choices"] = {"states": states,
                                 "contrasts": [contrast["label"]] if contrast else []}
            return JSONResponse(result)

        if intent == "unsupported":
            result["answer"] = (
                "I can answer four things about this dataset: the marker genes of a "
                "cell state, genes differing between two groups within a cell state, "
                "where named genes are expressed, and what separates two cell states.")
            return JSONResponse(result)

        from .chat_service import execute as execute_protocol
        protocol_result=execute_protocol(app,meta,question,reading)
        if protocol_result is not None:return JSONResponse(protocol_result)

        if intent in {"tf_activity", "regulator_activity", "regulatory_driver"}:
            ds = UploadedGrnData(app, meta)
            selected = ds.current_contrast
            if intent != "regulatory_driver" or reading.get("modality") == "grn_tf":
                answer = gnet.tf_activity_profile(ds, cell_state=state, contrast=selected,
                                                 factors=genes or None, limit=int(reading.get("limit") or 25))
                result.update(answer)
                if answer.get("by_cell_state"):
                    result.update(gnet.tf_activity_state_chat(answer))
                    return JSONResponse(result)
                rows = answer.get("rows") or []
                result["table"] = {"columns": ["factor", "activity", "log2fc", "fdr"],
                                   "rows": [[r[k] for k in ("factor", "activity", "log2fc", "fdr")] for r in rows]}
                result["plot"] = {"kind": "barchart", "label_column": "factor",
                                  "value_column": "activity", "sign_column": "log2fc"}
                result["answer"] = (f"{len(rows)} factors in {answer.get('cell_state')}. "
                                    f"Activity is {answer.get('statistic')}. "
                                    "Reported differential changes are listed separately; missing statistics are not zero."
                                    if rows else answer.get("note", "No TF activity available."))
            else:
                answer = gnet.regulator_network(ds, state, contrast=selected,
                                                features=genes or None, limit=int(reading.get("limit") or 200))
                result.update(answer)
                result["plot"] = {"kind": "network"}
                result["answer"] = (answer.get("note") or
                                    f"{answer.get('n_regulators', 0)} factors and {answer.get('n_edges_drawn', 0)} regulatory edges in {state}. "
                                    "Factor expression and activity changes are reported separately.")
            return JSONResponse(result)

        if intent in ("communication_rewiring", "pathway_program"):
            result["answer"] = "Open the Explore or Differential tab for this job's communication or GO Terms results."
            result["status"] = "use_existing_view"
            return JSONResponse(result)

        intent = _CHAT_PROTOCOL_ALIAS.get(intent, intent)

        if intent == "markers":
            if state and state not in states:
                result["answer"] = (f"{state} is not a cell state in this dataset. "
                                    f"It holds {len(states)}: {', '.join(states[:8])}"
                                    + (" ..." if len(states) > 8 else "."))
                result["status"] = "not_found"
                return JSONResponse(result)
            rows = _chat_markers_for_states(app, meta, cache, [state or states[0]], 25)
            sources = sorted({str(row["source"]) for row in rows})
            result["answer"] = (f"{len(rows)} top marker genes of {state or states[0]}, "
                                f"from the {' and '.join(sources) or 'marker'} analysis.")
            result["table"] = {
                "columns": ["gene", "cluster", "fold", "p", "source"],
                "rows": [[r["gene"], r["cluster"], round(float(r["fold"]), 4), r["p"], r["source"]]
                         for r in rows]}
            result["plot"] = {"kind": "dotplot", "genes": [r["gene"] for r in rows[:12]]}
            return JSONResponse(result)

        if intent == "compare":
            pair = [s for s in (state, state2) if s]
            if len(pair) < 2:
                result["answer"] = "Name two cell states to compare."
                result["status"] = "clarify"
                return JSONResponse(result)
            rows = _chat_markers_for_states(app, meta, cache, pair, 30)
            result["answer"] = f"Marker genes separating {pair[0]} and {pair[1]}."
            result["table"] = {
                "columns": ["gene", "cluster", "fold", "p", "source"],
                "rows": [[r["gene"], r["cluster"], round(float(r["fold"]), 4), r["p"], r["source"]]
                         for r in rows]}
            result["plot"] = {"kind": "dotplot", "genes": [r["gene"] for r in rows[:12]]}
            return JSONResponse(result)

        if intent == "expression":
            if not genes:
                result["answer"] = "Name at least one gene."
                result["status"] = "clarify"
                return JSONResponse(result)
            stats = _gene_state_stats(cache, genes)
            found = stats["genes"]
            if not found:
                result["answer"] = f"None of {', '.join(genes)} are in this dataset."
                result["status"] = "not_found"
                return JSONResponse(result)
            table_rows = []
            for row_index, gene in enumerate(found):
                means = np.asarray(stats["mean"][row_index], dtype=float)
                for column in np.argsort(-means)[:5]:
                    table_rows.append([gene, stats["states"][int(column)],
                                       round(float(means[int(column)]), 3),
                                       round(float(stats["frac"][row_index][int(column)]), 3)])
            missing = stats["missing"]
            result["answer"] = (
                f"Highest-expressing cell states for {', '.join(found)}, by mean expression "
                f"across the {len(stats['states'])} states of this dataset."
                + (f" Not in this dataset: {', '.join(missing)}." if missing else "")
                + (" Read as an expression lookup because the question names a gene this "
                   "dataset holds." if result.get("status") == "read_from_question" else ""))
            result["table"] = {"columns": ["gene", "cell state", "mean", "fraction"],
                               "rows": table_rows}
            result["plot"] = {"kind": "dotplot", "genes": found}
            return JSONResponse(result)

        if intent == "differential":
            requested = _normalize_modality_id(reading.get("modality"), default="")
            if requested and requested != (meta.get("differential", {}).get("config", {}).get("modality") or "rna"):
                ds = UploadedGrnData(app, meta)
                selected = gnet._sibling_comparison(ds, ds.current_contrast, requested)
                if not selected:
                    result.update(status="not_run", answer=f"No completed {requested} comparison matches the selected groups. Run it in the Differential tab first.")
                    return JSONResponse(result)
                meta = dict(meta, differential=ds.runs[selected])
                contrast = _chat_contrast(meta)
            if not contrast:
                result["answer"] = (
                    "This job has not run a differential comparison yet. Open the "
                    "Differential tab, choose the two sample groups and run it; then this "
                    "question has an answer. I am not substituting a marker analysis for it.")
                result["status"] = "not_run"
                return JSONResponse(result)
            try:
                detail = _get_differential_detail_table(app, meta)
            except (FileNotFoundError, ValueError) as exc:
                result["answer"] = f"The differential result is unavailable ({exc})."
                result["status"] = "not_available"
                return JSONResponse(result)
            covered = sorted(set(detail["population"].astype(str)))
            subset = detail.loc[detail["population"].astype(str) == state] if state else detail
            if state and subset.empty:
                result["answer"] = (
                    f"{contrast['label']} is not computed for {state}. It covers "
                    f"{len(covered)} of the {len(states)} cell states in this dataset. "
                    "This is missing data, not an absence of change.")
                result["status"] = "not_covered"
                result["table"] = {"columns": ["cell state covered"], "rows": [[s] for s in covered]}
                return JSONResponse(result)
            frame = subset.copy()
            frame["fdr"] = pd.to_numeric(frame.get("fdr"), errors="coerce")
            frame["log2fc"] = pd.to_numeric(frame.get("log2fc"), errors="coerce")
            frame["pval"] = pd.to_numeric(frame.get("pval"), errors="coerce")
            frame = frame.dropna(subset=["log2fc"]).sort_values(["fdr", "pval"])
            top = frame.head(25)
            result["answer"] = (
                f"Top {_modality_definition(meta, meta.get('differential', {}).get('config', {}).get('modality', 'rna'))['feature_label']} results for {contrast['label']}" + (f" in {state}" if state else "")
                + f", from this job's cellHarmony-differential result "
                  f"({len(frame)} genes reported{' in ' + state if state else ''}).")
            result["table"] = {
                "columns": ["gene", "population", "log2fc", "fdr", "pval"],
                "rows": [[str(row.gene), str(row.population),
                          float(row.log2fc),
                          float(row.fdr) if _is_finite_number(row.fdr) else None,
                          float(row.pval) if _is_finite_number(row.pval) else None]
                         for row in top.itertuples()]}
            result["plot"] = {"kind": "volcano"}
            return JSONResponse(result)

        result["answer"] = "I did not understand that."
        return JSONResponse(result)

    @app.api_route("/api/jobs/{job_id}/marker/heatmap.tsv", methods=["GET", "HEAD", "OPTIONS"])
    def marker_heatmap_tsv(
        job_id: str,
        request: Request,
        modality: str = Query("rna"),
        # True keeps the 10-metacells-per-population run, which is the readable default.
        # False serves the all-column run when the bundle ships one.
        compact: bool = Query(True),
        cells_per_sample: Optional[int] = Query(None),
        filter1_field: Optional[str] = Query(None),
        filter1_values: List[str] = Query([]),
        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        annotations = metabolite_annotations.for_dataset(meta, modality)
        display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
        if cells_per_sample is not None:
            matrix = _sampled_marker_heatmap(app, meta, modality, display_filters, cells_per_sample)
            headers = {"Access-Control-Allow-Origin": "*", "Access-Control-Allow-Methods": "GET, HEAD, OPTIONS",
                       "Access-Control-Allow-Headers": "*", "X-Marker-Columns": str(len(matrix["col_ids"])),
                       "X-Cells-Per-Sample": str(cells_per_sample), "X-Sample-Field": matrix["sampling"]["sample_field"]}
            content = "" if request.method in {"HEAD", "OPTIONS"} else _marker_heatmap_subset_to_tsv(matrix["matrix"], matrix["row_ids"], matrix["col_ids"])
            return Response(content=metabolite_annotations.label_tsv(content, annotations), media_type="text/tab-separated-values", headers=headers)
        try:
            matrix_entry = _get_marker_heatmap_cache_entry(app, meta, modality=modality, compact=compact)
        except FileNotFoundError:
            raise HTTPException(status_code=404, detail="Marker heatmap matrix unavailable.")
        origin = str(request.headers.get("origin") or "").strip()
        referer = str(request.headers.get("referer") or "").strip()
        user_agent = str(request.headers.get("user-agent") or "").strip()
        request_method = str(request.method or "").upper()
        source_path = matrix_entry.get("source_path")
        source_size = source_path.stat().st_size if isinstance(source_path, Path) and source_path.exists() else "-"
        store.append_log(
            job_id,
            f"[marker-heatmap] heatmap.tsv requested method={request_method} "
            f"source={matrix_entry.get('source')} path={source_path or '-'} size={source_size} "
            f"filters={display_filters or '-'} origin={origin or '-'} referer={referer or '-'} ua={user_agent or '-'}",
        )
        common_headers = {
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Methods": "GET, HEAD, OPTIONS",
            "Access-Control-Allow-Headers": "*",
        }
        if request_method == "OPTIONS":
            return JSONResponse({"ok": True}, headers=common_headers)
        if display_filters:
            filtered_tsv, kept_columns = _filter_marker_heatmap_matrix(app, meta, modality, display_filters)
            store.append_log(
                job_id,
                f"[marker-heatmap] filtered matrix generated columns={kept_columns} filters={display_filters}",
            )
            if request_method == "HEAD":
                return Response(
                    status_code=200,
                    headers={**common_headers, "X-Marker-Columns": str(kept_columns)},
                    media_type="text/tab-separated-values",
                )
            return Response(
                content=metabolite_annotations.label_tsv(filtered_tsv, annotations),
                media_type="text/tab-separated-values",
                headers={**common_headers, "X-Marker-Columns": str(kept_columns)},
            )
        col_count = int(len(matrix_entry["col_ids"]))
        if request_method == "HEAD":
            return Response(
                status_code=200,
                headers={**common_headers, "X-Marker-Columns": str(col_count)},
                media_type="text/tab-separated-values",
            )
        tsv_path = matrix_entry.get("tsv_path")
        if isinstance(tsv_path, Path) and tsv_path.exists() and not annotations:
            return FileResponse(
                tsv_path,
                filename=tsv_path.name,
                media_type="text/tab-separated-values",
                headers=common_headers,
            )
        payload = _marker_heatmap_subset_to_tsv(
            matrix_entry["matrix"],
            matrix_entry["row_ids"],
            matrix_entry["col_ids"],
        )
        return Response(
            content=metabolite_annotations.label_tsv(payload, annotations),
            media_type="text/tab-separated-values",
            headers={**common_headers, "X-Marker-Columns": str(col_count)},
        )

    @app.get("/jobs/{job_id}/marker/heatmap/viewer", response_class=HTMLResponse)
    def marker_heatmap_viewer(request: Request, job_id: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        modality = _normalize_modality_id(request.query_params.get("modality"), default="rna")
        compact = str(request.query_params.get("compact", "true")).strip().lower() \
            not in {"0", "false", "no", "off"}
        try:
            matrix_entry = _get_marker_heatmap_cache_entry(app, meta, modality=modality, compact=compact)
        except FileNotFoundError:
            raise HTTPException(status_code=404, detail="Marker heatmap matrix unavailable.")
        source_path = matrix_entry.get("source_path")
        source_size = source_path.stat().st_size if isinstance(source_path, Path) and source_path.exists() else "-"
        store.append_log(
            job_id,
            f"[marker-viewer] viewer route requested source={matrix_entry.get('source')} path={source_path or '-'} size={source_size}",
        )
        try:
            content = _build_marker_heatmap_viewer_html(
                job_id=job_id,
                root_path=_normalize_root_path(request.scope.get("root_path") or ""),
            )
            store.append_log(job_id, "[marker-viewer] viewer HTML generated successfully.")
            return HTMLResponse(content=content)
        except Exception as exc:
            store.append_log(job_id, f"[marker-viewer] viewer HTML generation failed: {exc!r}")
            raise

    @app.get("/api/jobs/{job_id}/marker/heatmap.pdf")
    def marker_heatmap_pdf(
        job_id: str,
        modality: str = Query("rna"),
        # True keeps the 10-metacells-per-population run, which is the readable default.
        # False serves the all-column run when the bundle ships one.
        compact: bool = Query(True),
        cells_per_sample: Optional[int] = Query(None),
        filter1_field: Optional[str] = Query(None),
        filter1_values: List[str] = Query([]),
        filter2_field: Optional[str] = Query(None),
        filter2_values: List[str] = Query([]),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
        if cells_per_sample is not None:
            matrix = _sampled_marker_heatmap(app, meta, modality, display_filters, cells_per_sample)
            frame = pd.DataFrame(matrix["matrix"], index=matrix["row_ids"], columns=matrix["col_ids"])
        else:
            try:
                matrix_entry = _get_marker_heatmap_cache_entry(app, meta, modality=modality, compact=compact)
            except FileNotFoundError:
                raise HTTPException(status_code=404, detail="Marker heatmap matrix unavailable.")
            if display_filters:
                filtered_tsv, _ = _filter_marker_heatmap_matrix(app, meta, modality, display_filters)
                frame = pd.read_csv(io.StringIO(filtered_tsv), sep="\t", index_col=0)
            else:
                frame = pd.DataFrame(
                    matrix_entry["matrix"],
                    index=matrix_entry["row_ids"],
                    columns=matrix_entry["col_ids"],
                )
        annotations = metabolite_annotations.for_dataset(meta, modality)
        frame = frame.rename(index=lambda name: metabolite_annotations.display_name(name, annotations))
        pdf = _render_marker_heatmap_pdf(frame, title="MarkerHeatmap")
        return StreamingResponse(
            pdf,
            media_type="application/pdf",
            headers={"Content-Disposition": f'attachment; filename="{job_id}_marker_heatmap.pdf"'},
        )

    @app.get("/api/jobs/{job_id}/marker/network/pdf")
    def marker_network_pdf(job_id: str, population: str = Query(...), modality: str = Query("rna")):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        payload = _build_marker_network_payload(meta, population, modality=modality)
        pdf = _render_network_pdf(payload, f"{population} marker network")
        safe_population = re.sub(r"[^A-Za-z0-9_.-]+", "_", population).strip("._") or "population"
        return StreamingResponse(
            pdf,
            media_type="application/pdf",
            headers={"Content-Disposition": f'attachment; filename="{job_id}_{safe_population}_marker_network.pdf"'},
        )

    @app.get("/api/jobs/{job_id}/download/{artifact}")
    async def download(job_id: str, artifact: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        artifact_map = meta.get("artifacts", {})
        if artifact not in artifact_map:
            raise HTTPException(status_code=404, detail=f"Artifact '{artifact}' unavailable.")
        path = Path(artifact_map[artifact])
        if not path.exists():
            raise HTTPException(status_code=404, detail="Artifact missing on disk.")
        return FileResponse(path, filename=path.name, headers={"Access-Control-Allow-Origin": "*"})

    @app.get("/api/jobs/{job_id}/log")
    async def download_log(job_id: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        path = store.logs_dir(job_id) / "pipeline.log"
        if not path.exists():
            raise HTTPException(status_code=404, detail="Log file unavailable.")
        return FileResponse(path, filename=path.name, media_type="text/plain")

    @app.get("/api/jobs/{job_id}/differential/archive")
    async def download_differential_archive(job_id: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        path = _get_differential_artifact(meta, "archive")
        return FileResponse(path, filename=path.name, media_type="application/zip")

    @app.get("/api/jobs/{job_id}/differential/artifact/{artifact_key}")
    async def download_differential_artifact(job_id: str, artifact_key: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        path = _get_differential_artifact(meta, artifact_key)
        suffix = path.suffix.lower()
        media_type = {
            ".pdf": "application/pdf",
            ".png": "image/png",
            ".svg": "image/svg+xml",
            ".tsv": "text/tab-separated-values",
            ".h5ad": "application/octet-stream",
        }.get(suffix, "application/octet-stream")
        return FileResponse(path, filename=path.name, media_type=media_type)

    @app.get("/api/jobs/{job_id}/differential/heatmap")
    async def download_differential_heatmap(job_id: str, format: str = Query("svg")):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        if format not in {"svg", "pdf", "png"}:
            raise HTTPException(status_code=400, detail="Heatmap format must be svg, pdf, or png.")
        artifact_key = {
            "pdf": "heatmap_pdf",
            "png": "heatmap_png",
            "svg": "heatmap_svg",
        }[format]
        path = _get_differential_artifact(meta, artifact_key)
        media_type = {
            "pdf": "application/pdf",
            "png": "image/png",
            "svg": "image/svg+xml",
        }[format]
        return FileResponse(path, filename=path.name, media_type=media_type)

    @app.get("/api/jobs/{job_id}/differential/network/{network_id}")
    async def download_differential_network(job_id: str, network_id: str, format: str = Query("png")):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        differential = meta.get("differential", {})
        selected = next((entry for entry in differential.get("networks", []) if entry.get("id") == network_id), None)
        if selected is None:
            raise HTTPException(status_code=404, detail="Differential network unavailable.")
        if format not in {"png", "pdf", "tsv"}:
            raise HTTPException(status_code=400, detail="Network format must be png, pdf, or tsv.")
        raw_path = str(selected.get(format, "")).strip()
        if not raw_path:
            raise HTTPException(status_code=404, detail="Differential network artifact unavailable.")
        path = Path(raw_path)
        if not path.exists():
            raise HTTPException(status_code=404, detail="Differential network artifact missing on disk.")
        media_type = {
            "png": "image/png",
            "pdf": "application/pdf",
            "tsv": "text/tab-separated-values",
        }[format]
        return FileResponse(path, filename=path.name, media_type=media_type)

    @app.get("/api/jobs/{job_id}/differential/interactive/heatmap")
    def differential_heatmap_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_heatmap_payload(app, meta, population))

    @app.get("/api/jobs/{job_id}/differential/interactive/summary")
    def differential_summary_data(job_id: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_summary_payload(app, meta))

    @app.get("/api/jobs/{job_id}/differential/interactive/volcano")
    def differential_volcano_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_volcano_payload(app, meta, population))

    @app.get("/api/jobs/{job_id}/differential/interactive/go")
    def differential_go_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_go_payload(app, meta, population))

    @app.get("/api/jobs/{job_id}/differential/interactive/network")
    def differential_network_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_network_payload(app, meta, population, root_path=app.state.root_path))

    @app.get("/api/jobs/{job_id}/differential/interactive/table")
    def differential_table_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        modality = _normalize_modality_id((meta.get("differential", {}) or {}).get("config", {}).get("modality"), default="rna")
        if modality != "cell_communication":
            raise HTTPException(status_code=404, detail="Differential interaction table is only available for Cell communication.")
        return JSONResponse(_build_differential_cell_communication_table_payload(app, meta, population))

    @app.get("/api/jobs/{job_id}/differential/interactive/gene")
    def differential_gene_data(
        job_id: str,
        population: str = Query(...),
        gene: str = Query(...),
        feature: Optional[str] = Query(None),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_gene_detail_payload(app, meta, population, gene, feature=feature))

    @app.get("/api/jobs/{job_id}/differential/interactive/gene/pdf")
    def differential_gene_pdf(
        job_id: str,
        population: str = Query(...),
        gene: str = Query(...),
        feature: Optional[str] = Query(None),
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        payload = _build_differential_gene_detail_payload(app, meta, population, gene, feature=feature)
        pdf = _render_differential_gene_pdf(metabolite_annotations.pdf_payload(payload, meta))
        safe_gene = re.sub(r"[^A-Za-z0-9_.-]+", "_", gene).strip("._") or "gene"
        safe_population = re.sub(r"[^A-Za-z0-9_.-]+", "_", population).strip("._") or "population"
        filename = f"differential_gene_{safe_gene}_{safe_population}.pdf"
        return StreamingResponse(
            pdf,
            media_type="application/pdf",
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )

    @app.get("/api/jobs/{job_id}/differential/interactive/pdf")
    def differential_rendered_pdf(job_id: str, mode: str = Query(...), population: str = Query("")):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        mode_key = str(mode or "").strip().lower()
        if mode_key == "summary":
            payload = _build_differential_summary_payload(app, meta)
            pdf = _render_differential_summary_pdf(payload)
        elif mode_key == "heatmap":
            payload = _build_differential_heatmap_payload(app, meta, population)
            pdf = _render_differential_heatmap_pdf(metabolite_annotations.pdf_payload(payload, meta))
        elif mode_key == "volcano":
            payload = _build_differential_volcano_payload(app, meta, population)
            pdf = _render_differential_volcano_pdf(payload)
        elif mode_key == "go":
            payload = _build_differential_go_payload(app, meta, population)
            pdf = _render_differential_go_pdf(payload)
        elif mode_key == "network":
            payload = _build_differential_network_payload(app, meta, population, root_path=app.state.root_path)
            pdf = _render_network_pdf(payload, f"{population} network")
        else:
            raise HTTPException(status_code=400, detail="Differential mode must be summary, heatmap, volcano, network, or go.")
        safe_population = re.sub(r"[^A-Za-z0-9_.-]+", "_", population).strip("._") or "population"
        filename = f"differential_{mode_key}_{safe_population}.pdf"
        return StreamingResponse(
            pdf,
            media_type="application/pdf",
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )

    @app.post("/api/tools/approximate-umap")
    async def run_approximate_umap(payload: ApproximateUMAPRequest):
        if not payload.query and not payload.query_clusters_tsv:
            raise HTTPException(
                status_code=422,
                detail="Provide either 'query' (h5ad) or 'query_clusters_tsv'.",
            )
        if not payload.reference and not (
            payload.reference_coords_tsv and payload.reference_clusters_tsv
        ):
            raise HTTPException(
                status_code=422,
                detail=(
                    "Provide either 'reference' (h5ad) or both "
                    "'reference_coords_tsv' and 'reference_clusters_tsv'."
                ),
            )

        # Every path below arrives in the request body, so each one is pinned to
        # the job storage or the reference tree before it is opened. Without this
        # the route reads any file the process can reach and returns its contents.
        input_roots = _tool_input_roots(app)
        custom_colors_tsv = (
            str(_require_within_any(payload.custom_colors_tsv, input_roots, "custom_colors_tsv"))
            if payload.custom_colors_tsv
            else None
        )

        if payload.reference:
            reference_source = ad.read_h5ad(
                str(_require_within_any(payload.reference, input_roots, "reference"))
            )
        else:
            reference_source = approx_mod._load_reference_from_tsv(
                str(_require_within_any(payload.reference_coords_tsv, input_roots, "reference_coords_tsv")),
                str(_require_within_any(payload.reference_clusters_tsv, input_roots, "reference_clusters_tsv")),
                umap_key=payload.umap_key,
                cluster_key=payload.reference_cluster_key or payload.query_cluster_key,
            )

        if payload.query:
            query_path = _require_within_any(payload.query, input_roots, "query")
            query_source = str(query_path)
        else:
            query_path = _require_within_any(
                payload.query_clusters_tsv, input_roots, "query_clusters_tsv"
            )
            query_source = approx_mod._load_query_from_tsv(
                str(query_path),
                cluster_key=payload.query_cluster_key,
            )

        if payload.verbose:
            logging.getLogger().setLevel(logging.INFO)

        result = approx_mod.approximate_umap(
            query=query_source,
            reference=reference_source,
            query_cluster_key=payload.query_cluster_key,
            reference_cluster_key=payload.reference_cluster_key,
            umap_key=payload.umap_key,
            jitter=payload.jitter,
            num_reference_cells=payload.num_reference_cells,
            random_state=payload.random_state,
            custom_color_tsv=custom_colors_tsv,
            restrict_obs_field=payload.restrict_obs_field,
            restrict_obs_value=payload.restrict_obs_value,
            copy_query=False,
        )

        # Outputs stay under job storage. outdir, output_prefix, output_h5ad and
        # output_pdf are all caller-supplied, and _resolve_output_path returns an
        # absolute candidate untouched, so each result is re-checked below.
        job_storage = Path(app.state.config["JOB_STORAGE"])
        job_storage.mkdir(parents=True, exist_ok=True)
        out_dir = _require_within(Path(payload.outdir), job_storage, "outdir")
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix_name = _secure_filename(
            payload.output_prefix or f"{query_path.stem}-approximate-umap"
        )
        prefix = _require_within(out_dir / prefix_name, job_storage, "output_prefix")

        export_approx_pdfs = bool(app.state.config.get("EXPORT_APPROX_PDFS", False)) or bool(payload.output_pdf)
        output_pdf: Optional[Path] = None
        if export_approx_pdfs:
            output_pdf = _require_within(
                _resolve_output_path(
                    out_dir,
                    payload.output_pdf or f"{prefix_name}-comparison.pdf",
                ),
                job_storage,
                "output_pdf",
            )
            output_pdf.parent.mkdir(parents=True, exist_ok=True)

        save_h5ad = bool(payload.save_updated_h5ad)
        output_h5ad: Optional[Path] = None
        if payload.output_h5ad:
            output_h5ad = _require_within(
                _resolve_output_path(out_dir, payload.output_h5ad), job_storage, "output_h5ad"
            )
            save_h5ad = True
        elif save_h5ad:
            output_h5ad = out_dir / f"{prefix_name}.h5ad"

        if save_h5ad and output_h5ad is not None:
            output_h5ad.parent.mkdir(parents=True, exist_ok=True)
            compression = _normalize_h5ad_compression(app.state.config.get("H5AD_COMPRESSION", "lzf"))
            approx_mod.ensure_h5ad_compat_for_write(result.query_adata)
            result.query_adata.write(output_h5ad, compression=compression)

        text_outputs = result.write_text_outputs(prefix)
        annotated_pdf = None
        plain_pdf = None
        if export_approx_pdfs and output_pdf is not None:
            reference_for_plot = (
                reference_source if isinstance(reference_source, ad.AnnData) else ad.read_h5ad(str(reference_source))
            )
            annotated_pdf, plain_pdf = result.write_comparison_pdf(
                reference_adata=reference_for_plot,
                umap_key=payload.umap_key,
                reference_cluster_key=payload.reference_cluster_key or payload.query_cluster_key,
                query_cluster_key=payload.query_cluster_key,
                output_path=output_pdf,
                custom_color_map=result.plot_color_map,
                restrict_obs_field=payload.restrict_obs_field,
                restrict_obs_value=payload.restrict_obs_value,
            )

        outputs = {
            "coordinates": text_outputs["coordinates"],
            "augmented": text_outputs["augmented"],
            "placeholder_expression": text_outputs["placeholder_expression"],
        }
        if annotated_pdf and plain_pdf:
            outputs["comparison_pdf"] = annotated_pdf
            outputs["comparison_pdf_plain"] = plain_pdf
        if output_h5ad is not None:
            outputs["updated_h5ad"] = str(output_h5ad)

        return JSONResponse(
            {
                "status": "completed",
                "query_cluster_key": payload.query_cluster_key,
                "reference_cluster_key": payload.reference_cluster_key or payload.query_cluster_key,
                "umap_key": payload.umap_key,
                "cluster_order": result.cluster_order,
                "n_query_cells": int(result.query_adata.n_obs),
                "outputs": outputs,
            }
        )

    from .integration_routes import install as install_integrated_routes
    install_integrated_routes(app)
    return app


app = create_app()
