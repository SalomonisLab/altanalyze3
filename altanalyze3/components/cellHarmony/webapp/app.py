from __future__ import annotations

import io
import json
import logging
import re
import shutil
import threading
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

from .config import BASE_DIR, load_config


class QCSettings(BaseModel):
    min_genes: int = 500
    min_counts: int = 1000
    min_cells: int = 0
    mit_percent: int = 15
    align_cutoff: float = 0.4
    ambient_correction: str = Field(default="no", pattern="^(no|yes)$")


class JobConfigSettings(BaseModel):
    species: str
    reference: str
    soupx_option: Optional[str] = None


class DifferentialSettings(BaseModel):
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


class ClientLogRequest(BaseModel):
    message: str


def _secure_filename(filename: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", filename).strip("._")
    return cleaned or "upload"


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
    return {
        "enabled": enabled,
        "sample_names": sample_names,
        "sample_fields": sample_fields,
        "sample_values": sample_values,
        "default_sample_field": default_sample_field,
        "population_columns": population_columns,
        "default_population_col": default_population_col,
        "comparison_types": comparison_types,
        "upload_profile": upload_profile,
    }


def _build_differential_payload(app: FastAPI, job_id: str, meta: Dict, root_path: str = "") -> Dict:
    options = _differential_options(meta)
    differential = dict(meta.get("differential") or {})
    artifacts = differential.get("artifacts") or {}
    result_populations: List[str] = []
    visualization_populations: Dict[str, List[str]] = {"heatmap": [], "volcano": [], "network": [], "go": []}
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
        visualization_populations["network"] = list(
            dict.fromkeys(
                str(entry.get("population", "")).strip()
                for entry in differential.get("networks") or []
                if str(entry.get("population", "")).strip()
            )
        )
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
    return {
        **options,
        "status": status,
        "progress": int(differential.get("progress", 0) or 0),
        "message": differential.get("message") or (
            "Select a cell-state field and numerator/denominator sample groups to run cellHarmony-differential."
            if options["enabled"]
            else "Differential gene analyses between biological groups (i.e., disease versus controls) are only enabled when two or more samples (multiple h5 files or a single h5ad) are uploaded for the job."
        ),
        "config": differential.get("config") or {},
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


def _differential_result_populations(app: FastAPI, meta: Dict) -> List[str]:
    detailed = _get_differential_detail_table(app, meta)
    ordered = list(dict.fromkeys(detailed["population"].astype(str).tolist()))
    return [value for value in ordered if value]


def _differential_heatmap_populations(app: FastAPI, meta: Dict) -> List[str]:
    try:
        detailed = _get_differential_detail_table(app, meta)
    except Exception:
        return []
    ordered = list(dict.fromkeys(detailed["population"].astype(str).tolist()))
    return [value for value in ordered if value]


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
    raw_combined = str(meta.get("artifacts", {}).get("combined_h5ad", "")).strip()
    if raw_combined and Path(raw_combined).exists():
        return Path(raw_combined)
    raise FileNotFoundError("Aligned AnnData output is unavailable for differential gene detail.")


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
        cache.pop(str(job_id), None)
    locks = getattr(app.state, "expression_cache_locks", None)
    if isinstance(locks, dict):
        locks.pop(str(job_id), None)


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
        cache.pop(str(job_id), None)


def _get_marker_heatmap_cache_entry(app: FastAPI, meta: Dict) -> Dict[str, Any]:
    job_id = str(meta.get("job_id") or "").strip()
    if not job_id:
        raise ValueError("Job metadata is missing a job_id.")

    marker_analysis = meta.get("marker_analysis") or {}
    cache_path_text = str(marker_analysis.get("heatmap_cache", "")).strip()
    heatmap_tsv_text = str(marker_analysis.get("heatmap_tsv", "")).strip()
    expression_tsv_text = str(marker_analysis.get("expression_tsv", "")).strip()

    cache_path = Path(cache_path_text) if cache_path_text else None
    heatmap_tsv_path = Path(heatmap_tsv_text) if heatmap_tsv_text else None
    expression_tsv_path = Path(expression_tsv_text) if expression_tsv_text else None

    cache_signature = {
        "heatmap_cache": str(cache_path or ""),
        "heatmap_tsv": str(heatmap_tsv_path or ""),
        "expression_tsv": str(expression_tsv_path or ""),
    }
    cache = app.state.marker_heatmap_cache
    existing = cache.get(job_id)
    if existing and existing.get("signature") == cache_signature:
        return existing

    tsv_path: Optional[Path] = None
    if heatmap_tsv_path and heatmap_tsv_path.exists():
        tsv_path = heatmap_tsv_path
    elif expression_tsv_path and expression_tsv_path.exists():
        tsv_path = expression_tsv_path

    if cache_path and cache_path.exists():
        with np.load(cache_path, allow_pickle=False) as npz_data:
            matrix = np.asarray(npz_data["matrix"], dtype=np.float32)
            row_ids = np.asarray(npz_data["row_ids"], dtype=str)
            col_ids = np.asarray(npz_data["col_ids"], dtype=str)
            if "col_barcodes" in npz_data:
                col_barcodes = np.asarray(npz_data["col_barcodes"], dtype=str)
            else:
                col_barcodes = np.asarray(
                    [value.split(":", 1)[1] if ":" in value else value for value in col_ids],
                    dtype=str,
                )
        entry = {
            "signature": cache_signature,
            "source": "cache",
            "source_path": cache_path,
            "matrix": matrix,
            "row_ids": row_ids,
            "col_ids": col_ids,
            "col_barcodes": col_barcodes,
            "tsv_path": tsv_path,
        }
        cache[job_id] = entry
        return entry

    if tsv_path and tsv_path.exists():
        frame = pd.read_csv(tsv_path, sep="\t", index_col=0)
        col_ids = frame.columns.astype(str).to_numpy()
        entry = {
            "signature": cache_signature,
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
        cache[job_id] = entry
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


def _get_expression_cache(app: FastAPI, meta: Dict) -> Dict[str, Any]:
    job_id = str(meta.get("job_id") or "").strip()
    if not job_id:
        raise ValueError("Job metadata is missing a job_id.")

    artifacts = meta.get("artifacts", {})
    h5ad_path = artifacts.get("combined_h5ad")
    if not h5ad_path or not Path(h5ad_path).exists():
        raise FileNotFoundError("Combined AnnData output unavailable for this job.")

    cluster_key = meta.get("cluster_key")
    if not cluster_key:
        raise ValueError("Cluster assignments missing from AnnData output.")

    umap_path = artifacts.get("umap_coordinates")
    cache = app.state.expression_cache
    cache_entry = cache.get(job_id)
    if (
        cache_entry
        and cache_entry.get("h5ad_path") == str(h5ad_path)
        and cache_entry.get("umap_path") == str(umap_path or "")
        and cache_entry.get("cluster_key") == str(cluster_key)
    ):
        return cache_entry

    lock = _get_cache_lock(app, "expression_cache_locks", job_id)
    with lock:
        cache_entry = cache.get(job_id)
        if (
            cache_entry
            and cache_entry.get("h5ad_path") == str(h5ad_path)
            and cache_entry.get("umap_path") == str(umap_path or "")
            and cache_entry.get("cluster_key") == str(cluster_key)
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
            "h5ad_path": str(h5ad_path),
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
            "obs_filter_values": obs_filter_values,
            "display_filters_meta": {
                "fields": filter_fields,
                "values": filter_values,
                "default_primary_field": default_primary_field,
                "default_secondary_field": default_secondary_field,
                "show_secondary": bool(upload_profile.get("show_secondary_display_filter", True)),
            },
        }
        cache[job_id] = cache_entry
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
        "go_table": None,
        "network_tables": {},
        "primary_adata": None,
        "fallback_adata": None,
        "primary_var_names": None,
        "fallback_var_names": None,
    }
    cache[job_id] = entry
    return entry


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
    heatmap_exact_rows: Dict[tuple[str, str], List[Optional[float]]] = {}
    use_heatmap_columns = False

    try:
        heatmap_frame = _get_differential_heatmap_table(app, meta)
        heatmap_columns = [str(label).split(":", 1)[0] for label in heatmap_frame.columns]
        column_labels = list(dict.fromkeys(heatmap_columns))
        if len(column_labels) != len(heatmap_frame.columns):
            column_labels = [str(label) for label in heatmap_frame.columns]
        heatmap_frame = heatmap_frame.copy()
        heatmap_frame.columns = column_labels
        for raw_key, values in heatmap_frame.iterrows():
            parsed = _parse_heatmap_row_key(raw_key)
            gene = parsed["gene"]
            row_values = [float(value) if _is_finite_number(value) else 0.0 for value in values.tolist()]
            row_key = (parsed["population"], gene)
            heatmap_exact_rows.setdefault(row_key, row_values)
        use_heatmap_columns = any(parsed_population == population for parsed_population, _ in heatmap_exact_rows.keys())
    except Exception:
        use_heatmap_columns = False

    if not use_heatmap_columns:
        column_labels = list(dict.fromkeys(detailed["population"].astype(str).tolist()))

    detailed_matrix = detailed_matrix.reindex(columns=column_labels).fillna(0.0)

    rows = []
    for row in subset.itertuples():
        gene = str(row.gene)
        values = (
            heatmap_exact_rows.get((population, gene))
            or (detailed_matrix.loc[gene].tolist() if gene in detailed_matrix.index else [0.0] * len(column_labels))
        )
        row_values = [float(value) if _is_finite_number(value) else 0.0 for value in values]
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
    subset["fdr_clamped"] = subset["fdr"].fillna(subset["pval"]).clip(lower=1e-300)
    subset["score"] = -np.log10(subset["fdr_clamped"])
    subset["direction"] = np.where(subset["log2fc"] >= 0, "up", "down")
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
        "points": points,
        "default_gene": points[0]["gene"] if points else None,
    }


def _build_differential_go_payload(app: FastAPI, meta: Dict, population: str) -> Dict:
    try:
        frame = _get_differential_go_table(app, meta)
    except HTTPException:
        return {"population": population, "terms": [], "default_gene": None}
    except FileNotFoundError:
        return {"population": population, "terms": [], "default_gene": None}

    if frame.empty or "population" not in frame.columns:
        return {"population": population, "terms": [], "default_gene": None}

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
    subset["score"] = -np.log10(subset["fdr"].fillna(subset["p_value"]).clip(lower=1e-300))
    subset = subset.sort_values(["fdr", "p_value", "term_name"], ascending=[True, True, True]).head(15)
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
                "score": float(row.score) if _is_finite_number(row.score) else None,
                "overlap_genes": overlap_genes,
                "selected_gene": overlap_genes[0] if overlap_genes else None,
            }
        )
    default_gene = next((term["selected_gene"] for term in terms if term["selected_gene"]), None)
    return {"population": population, "terms": terms, "default_gene": default_gene}


def _build_differential_network_payload(app: FastAPI, meta: Dict, population: str, root_path: str = "") -> Dict:
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
                    "log2fc": float(source_fc) if _is_finite_number(source_fc) else 0.0,
                }
            }
        for target in targets:
            if target not in nodes:
                target_fc = fc_map.get(target)
                nodes[target] = {
                    "data": {
                        "id": target,
                        "label": target,
                        "log2fc": float(target_fc) if _is_finite_number(target_fc) else 0.0,
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


def _build_marker_network_payload(meta: Dict, population: str) -> Dict:
    marker_analysis = meta.get("marker_analysis") or {}
    networks = marker_analysis.get("networks") or []
    entry = next((item for item in networks if str(item.get("population", "")).strip() == str(population).strip()), None)
    if not entry:
        return {"population": population, "elements": [], "default_gene": None}

    tsv_path = Path(str(entry.get("tsv", "")).strip())
    if not tsv_path.exists():
        return {"population": population, "elements": [], "default_gene": None}

    marker_stats_path = Path(str(marker_analysis.get("markers_tsv", "")).strip())
    fc_map: Dict[str, float] = {}
    if marker_stats_path.exists():
        stats_df = pd.read_csv(marker_stats_path, sep="\t")
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
            source_fc = fc_map.get(source, 0.0)
            nodes[source] = {"data": {"id": source, "label": source, "log2fc": float(source_fc)}}
        for target in targets:
            if target not in nodes:
                target_fc = fc_map.get(target, 0.0)
                nodes[target] = {"data": {"id": target, "label": target, "log2fc": float(target_fc)}}
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


def _build_differential_gene_detail_payload(app: FastAPI, meta: Dict, population: str, gene: str) -> Dict:
    differential = meta.get("differential", {}) or {}
    config = differential.get("config", {}) or {}
    group1_samples = [str(value).strip() for value in config.get("group1_samples", []) if str(value).strip()]
    group2_samples = [str(value).strip() for value in config.get("group2_samples", []) if str(value).strip()]
    sample_field = str(config.get("sample_field", "")).strip()
    if not group1_samples or not group2_samples:
        raise ValueError("Differential comparison groups are not configured.")
    population_col = _differential_population_col(meta)
    case_label, control_label = _differential_group_labels(meta)

    expression_cache = _get_expression_cache(app, meta)
    adata = expression_cache["adata"]
    resolved_gene = _resolve_gene_name(expression_cache["var_names"], gene)
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
    mask = (population_values == population) & sample_values.isin(resolved_group1 + resolved_group2)
    if int(np.asarray(mask).sum()) == 0:
        raise HTTPException(status_code=404, detail=f"No cells were found for '{population}' in the selected comparison groups.")

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
        "groups": groups,
        "stats": stats_payload,
    }


def _validate_differential_request(meta: Dict, payload: DifferentialSettings) -> Dict:
    options = _differential_options(meta)
    if not options["enabled"]:
        raise HTTPException(status_code=400, detail="Differential analysis requires two or more uploaded samples.")
    if meta.get("status") != "completed":
        raise HTTPException(status_code=400, detail="Run the cellHarmony analysis before starting differential analysis.")

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


def _build_umap_payload(
    app: FastAPI,
    meta: Dict,
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
) -> Dict[str, List[Dict]]:
    cache_entry = _get_expression_cache(app, meta)
    obs_names = cache_entry["obs_names"]
    populations = cache_entry["populations"]
    sample_field = str(cache_entry.get("sample_field") or "").strip()
    sample_labels = cache_entry.get("sample_labels")
    umap_x = cache_entry["umap_x"]
    umap_y = cache_entry["umap_y"]
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

    reference_points = []
    ref_adata = _load_reference_adata(app, meta)
    if ref_adata is not None:
        ref_cluster_key = meta.get("reference_cluster_key") or cluster_key
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
    return {"reference": reference_points, "query": query_points, "sample_field": sample_field}


def _reference_entry_for_meta(meta: Dict) -> Dict:
    registry_path = Path(app.state.config["REFERENCE_REGISTRY"]) if "app" in globals() else None
    if registry_path is None:
        raise ValueError("Reference registry is unavailable.")
    reference_entry = pipeline_mod._lookup_reference(meta["species"], meta["reference"], registry_path)
    pipeline_mod._ensure_reference_fields(reference_entry)
    return reference_entry


def _load_reference_states_table(meta: Dict) -> pd.DataFrame:
    reference_entry = _reference_entry_for_meta(meta)
    states_path = Path(reference_entry["states_tsv"])
    if not states_path.exists():
        raise FileNotFoundError("Reference states TSV is unavailable for this job.")
    return pd.read_csv(states_path, sep="\t", index_col=0)


def _build_reference_expression_payload(meta: Dict, gene: str) -> Dict:
    reference_df = _load_reference_states_table(meta)
    resolved_gene = _resolve_gene_name(reference_df.index.astype(str), gene)
    if not resolved_gene:
        combined_path = str(meta.get("artifacts", {}).get("combined_h5ad", "")).strip()
        if combined_path and Path(combined_path).exists():
            fallback_adata = ad.read_h5ad(combined_path, backed="r")
            try:
                if len(fallback_adata.var_names):
                    fallback_gene = str(fallback_adata.var_names[0])
                    return _build_expression_payload(app, meta, fallback_gene)
            finally:
                _close_backed_adata(fallback_adata)
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
    }


def _build_expression_payload(
    app: FastAPI,
    meta: Dict,
    gene: str,
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
) -> Dict:
    cache_entry = _get_expression_cache(app, meta)
    adata = cache_entry["adata"]
    populations = cache_entry["populations"]
    obs_names = cache_entry["obs_names"]
    umap_x = cache_entry["umap_x"]
    umap_y = cache_entry["umap_y"]
    display_mask = _apply_display_filter_mask(cache_entry, display_filters)
    resolved_gene = _resolve_gene_name(cache_entry["var_names"], gene)
    if not resolved_gene:
        if len(cache_entry["var_names"]):
            fallback_gene = str(cache_entry["var_names"][0])
            values = _flatten_expr(adata[:, fallback_gene].X)
            scatter_data = [
                {"population": pop, "value": float(val)}
                for pop, val in zip(populations, values)
                if _is_finite_number(val)
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
                violin_data.append(
                    {
                        "population": pop,
                        "values": [float(v) for v in finite_values],
                        "mean": float(np.mean(finite_values)) if len(finite_values) else 0.0,
                    }
                )
            violin_data = sorted(violin_data, key=lambda x: x["mean"], reverse=True)[:10]

            return {
                "gene": fallback_gene,
                "requested_gene": gene,
                "resolved_gene": fallback_gene,
                "source": "query",
                "message": None if int(np.asarray(display_mask).sum()) else "No cells match the current Display only filters.",
                "scatter": scatter_data,
                "violin": violin_data,
                "umap": umap_points,
            }
        return _build_reference_expression_payload(meta, gene)

    values = _flatten_expr(adata[:, resolved_gene].X)
    scatter_data = [
        {"population": pop, "value": float(val)}
        for pop, val in zip(populations, values)
        if _is_finite_number(val)
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
        violin_data.append(
            {
                "population": pop,
                "values": [float(v) for v in finite_values],
                "mean": float(np.mean(finite_values)) if len(finite_values) else 0.0,
            }
        )
    violin_data = sorted(violin_data, key=lambda x: x["mean"], reverse=True)[:10]

    return {
        "gene": resolved_gene,
        "requested_gene": gene,
        "resolved_gene": resolved_gene,
        "source": "query",
        "message": None if int(np.asarray(display_mask).sum()) else "No cells match the current Display only filters.",
        "scatter": scatter_data,
        "violin": violin_data,
        "umap": umap_points,
    }


def _build_gene_suggestions_payload(app: FastAPI, meta: Dict) -> Dict:
    cache_entry = _get_expression_cache(app, meta)
    genes = [str(gene) for gene in cache_entry["var_names"].tolist()]
    return {"genes": genes}


def _build_display_filter_payload(app: FastAPI, meta: Dict) -> Dict[str, Any]:
    cache_entry = _get_expression_cache(app, meta)
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


def _filter_marker_heatmap_matrix(
    app: FastAPI,
    meta: Dict,
    display_filters: Optional[List[tuple[str, List[str]]]],
) -> tuple[str, int]:
    cache_entry = _get_expression_cache(app, meta)
    marker_matrix = _get_marker_heatmap_cache_entry(app, meta)
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
    if mode == "violin":
        violin_data = payload["violin"]
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
        ax.set_title(f"{gene} expression (top 10 states by mean)")
        ax.set_ylabel("Expression")
    else:
        umap_points = payload["umap"]
        expression_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            "expression_grey_red",
            [
                (0.0, "#e5e7eb"),
                (0.15, "#f3f4f6"),
                (0.35, "#fecaca"),
                (0.6, "#f87171"),
                (1.0, "#b91c1c"),
            ],
        )
        sc = ax.scatter(
            [p["x"] for p in umap_points],
            [p["y"] for p in umap_points],
            s=4,
            c=[p["value"] for p in umap_points],
            cmap=expression_cmap,
            linewidths=0,
        )
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label(gene)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title(f"{gene} expression")
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
    ax.set_ylabel("Normalized expression")
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
    figure_height = max(6.0, min(28.0, len(rows) * 0.18 + 2.0))
    figure_width = max(7.0, min(18.0, len(columns) * 0.65 + 2.5))
    fig, ax = plt.subplots(figsize=(figure_width, figure_height))
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-color_extent, vmax=color_extent)
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
    ax.set_ylabel("-log10(FDR)")
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

    ordered = list(reversed(terms))
    scores = [float(term.get("score", 0.0) or 0.0) for term in ordered]
    labels = [str(term.get("term_name", "")) for term in ordered]
    colors = ["#dc2626" if str(term.get("direction", "")).lower() == "up" else "#2563eb" for term in ordered]

    fig_height = max(5.5, min(18.0, len(ordered) * 0.42 + 1.6))
    fig, ax = plt.subplots(figsize=(9.0, fig_height))
    positions = np.arange(len(ordered), dtype=float)
    ax.barh(positions, scores, color=colors, edgecolor="none", height=0.78)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("-log10(FDR)")
    ax.set_title(f"GO terms: {payload.get('population', '')}")
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
    positions = _network_positions(ordered_nodes)
    fig, ax = plt.subplots(figsize=(8.2, 7.4))
    for edge in edges:
        source = str(edge.get("source", "")).strip()
        target = str(edge.get("target", "")).strip()
        if source not in positions or target not in positions:
            continue
        interaction_type = str(edge.get("interaction_type", "")).lower()
        color = "#9ca3af"
        if "transcription" in interaction_type:
            color = "#ef4444"
        elif "tbar" in interaction_type:
            color = "#60a5fa"
        start = positions[source]
        end = positions[target]
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops={
                "arrowstyle": "-|>",
                "color": color,
                "linewidth": 1.0,
                "shrinkA": 14,
                "shrinkB": 14,
                "alpha": 0.8,
            },
        )

    for node_id in ordered_nodes:
        data = node_map[node_id]
        x, y = positions[node_id]
        log2fc = float(data.get("log2fc", 0.0) or 0.0)
        color = "#fca5a5" if log2fc >= 0 else "#7dd3fc"
        ax.scatter([x], [y], s=320, c=[color], edgecolors="white", linewidths=1.0, zorder=3)
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
        js_version = int((static_dir / "app.js").stat().st_mtime) if (static_dir / "app.js").exists() else 0
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
    async def reference_preview(species: str = Query(...), reference: str = Query(...)):
        return JSONResponse(_build_reference_preview_payload(app, species, reference))

    @app.post("/api/jobs")
    async def create_job(
        species: str = Form(...),
        reference: str = Form(...),
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

        metadata = store.create_job(species, reference, soupx_option, files=[])
        job_id = metadata["job_id"]
        _invalidate_expression_cache(app, job_id)
        _invalidate_differential_cache(app, job_id)
        _invalidate_marker_heatmap_cache(app, job_id)
        uploads_dir = store.uploads_dir(job_id)
        records = []
        for sample, upload in zip(normalized_samples, files):
            if not upload.filename:
                raise HTTPException(status_code=400, detail="One uploaded file is missing a filename.")
            if not _allowed_file(app, upload.filename):
                raise HTTPException(status_code=400, detail=f"Unsupported file extension for {upload.filename}.")
            dest_name = f"{sample}_{_secure_filename(upload.filename)}"
            dest_path = uploads_dir / dest_name
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
        requested_soupx = config.soupx_option if config.soupx_option is not None else meta.get("soupx_option")
        changed = (
            str(meta.get("species") or "") != config.species
            or str(meta.get("reference") or "") != config.reference
            or str(meta.get("soupx_option") or "") != str(requested_soupx or "")
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
                soupx_option=requested_soupx,
                status="uploaded",
                progress=0,
                message="Reference updated. Configure QC and rerun alignment.",
                artifacts={},
                marker_analysis={},
                differential={},
            )
        return JSONResponse(
            {
                "job_id": job_id,
                "status": meta.get("status"),
                "species": meta.get("species"),
                "reference": meta.get("reference"),
                "soupx_option": meta.get("soupx_option"),
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
        store.update_job(job_id, marker_analysis={})
        runner.submit(job_id)
        store.append_log(job_id, "Job queued by user request.")
        store.update_job(job_id, message="Job submitted to worker.", status="queued", progress=15)
        return JSONResponse({"job_id": job_id, "status": "queued"})

    @app.get("/api/jobs/{job_id}/status")
    async def status(job_id: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
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
        meta["differential_ui"] = _build_differential_payload(app, job_id, meta, root_path=app.state.root_path)
        return JSONResponse(meta, headers={"Cache-Control": "no-store"})

    @app.get("/api/jobs/{job_id}/differential/status")
    async def differential_status(job_id: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(
            _build_differential_payload(app, job_id, meta, root_path=app.state.root_path),
            headers={"Cache-Control": "no-store"},
        )

    @app.post("/api/jobs/{job_id}/differential")
    async def run_differential(job_id: str, payload: DifferentialSettings):
        store, runner = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        config = _validate_differential_request(meta, payload)
        options = _differential_options(meta)
        _invalidate_differential_cache(app, job_id)
        store.update_job(
            job_id,
            differential={
                "status": "queued",
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
    async def umap(
        job_id: str,
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
            return JSONResponse(_build_umap_payload(app, meta, display_filters))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/expression")
    async def expression(
        job_id: str,
        gene: str = Query(...),
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
            return JSONResponse(_build_expression_payload(app, meta, gene, display_filters))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/marker/network")
    async def marker_network(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        try:
            return JSONResponse(_build_marker_network_payload(meta, population))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/genes")
    async def job_genes(job_id: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        try:
            return JSONResponse(_build_gene_suggestions_payload(app, meta))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/display-filters")
    async def job_display_filters(job_id: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        try:
            return JSONResponse(_build_display_filter_payload(app, meta))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/umap/pdf")
    async def umap_pdf(
        job_id: str,
        mode: str = Query("relative"),
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
            payload = _build_umap_payload(app, meta, display_filters)
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
    async def expression_pdf(
        job_id: str,
        gene: str = Query(...),
        mode: str = Query("umap"),
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
            payload = _build_expression_payload(app, meta, gene, display_filters)
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))
        pdf = _render_expression_pdf(payload, mode)
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

    @app.api_route("/api/jobs/{job_id}/marker/heatmap.tsv", methods=["GET", "HEAD", "OPTIONS"])
    async def marker_heatmap_tsv(
        job_id: str,
        request: Request,
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
            matrix_entry = _get_marker_heatmap_cache_entry(app, meta)
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
            filtered_tsv, kept_columns = _filter_marker_heatmap_matrix(app, meta, display_filters)
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
                content=filtered_tsv,
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
        if isinstance(tsv_path, Path) and tsv_path.exists():
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
            content=payload,
            media_type="text/tab-separated-values",
            headers={**common_headers, "X-Marker-Columns": str(col_count)},
        )

    @app.get("/jobs/{job_id}/marker/heatmap/viewer", response_class=HTMLResponse)
    async def marker_heatmap_viewer(request: Request, job_id: str):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        try:
            matrix_entry = _get_marker_heatmap_cache_entry(app, meta)
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
    async def marker_heatmap_pdf(
        job_id: str,
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
            matrix_entry = _get_marker_heatmap_cache_entry(app, meta)
        except FileNotFoundError:
            raise HTTPException(status_code=404, detail="Marker heatmap matrix unavailable.")
        if display_filters:
            filtered_tsv, _ = _filter_marker_heatmap_matrix(app, meta, display_filters)
            frame = pd.read_csv(io.StringIO(filtered_tsv), sep="\t", index_col=0)
        else:
            frame = pd.DataFrame(
                matrix_entry["matrix"],
                index=matrix_entry["row_ids"],
                columns=matrix_entry["col_ids"],
            )
        pdf = _render_marker_heatmap_pdf(frame, title="MarkerHeatmap")
        return StreamingResponse(
            pdf,
            media_type="application/pdf",
            headers={"Content-Disposition": f'attachment; filename="{job_id}_marker_heatmap.pdf"'},
        )

    @app.get("/api/jobs/{job_id}/marker/network/pdf")
    async def marker_network_pdf(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        payload = _build_marker_network_payload(meta, population)
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
    async def differential_heatmap_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_heatmap_payload(app, meta, population))

    @app.get("/api/jobs/{job_id}/differential/interactive/volcano")
    async def differential_volcano_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_volcano_payload(app, meta, population))

    @app.get("/api/jobs/{job_id}/differential/interactive/go")
    async def differential_go_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_go_payload(app, meta, population))

    @app.get("/api/jobs/{job_id}/differential/interactive/network")
    async def differential_network_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_network_payload(app, meta, population, root_path=app.state.root_path))

    @app.get("/api/jobs/{job_id}/differential/interactive/gene")
    async def differential_gene_data(job_id: str, population: str = Query(...), gene: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        return JSONResponse(_build_differential_gene_detail_payload(app, meta, population, gene))

    @app.get("/api/jobs/{job_id}/differential/interactive/gene/pdf")
    async def differential_gene_pdf(job_id: str, population: str = Query(...), gene: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        payload = _build_differential_gene_detail_payload(app, meta, population, gene)
        pdf = _render_differential_gene_pdf(payload)
        safe_gene = re.sub(r"[^A-Za-z0-9_.-]+", "_", gene).strip("._") or "gene"
        safe_population = re.sub(r"[^A-Za-z0-9_.-]+", "_", population).strip("._") or "population"
        filename = f"differential_gene_{safe_gene}_{safe_population}.pdf"
        return StreamingResponse(
            pdf,
            media_type="application/pdf",
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )

    @app.get("/api/jobs/{job_id}/differential/interactive/pdf")
    async def differential_rendered_pdf(job_id: str, mode: str = Query(...), population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        mode_key = str(mode or "").strip().lower()
        if mode_key == "heatmap":
            payload = _build_differential_heatmap_payload(app, meta, population)
            pdf = _render_differential_heatmap_pdf(payload)
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
            raise HTTPException(status_code=400, detail="Differential mode must be heatmap, volcano, network, or go.")
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

        if payload.reference:
            reference_source = ad.read_h5ad(payload.reference)
        else:
            reference_source = approx_mod._load_reference_from_tsv(
                payload.reference_coords_tsv,
                payload.reference_clusters_tsv,
                umap_key=payload.umap_key,
                cluster_key=payload.reference_cluster_key or payload.query_cluster_key,
            )

        if payload.query:
            query_source = payload.query
            query_path = Path(payload.query)
        else:
            query_source = approx_mod._load_query_from_tsv(
                payload.query_clusters_tsv,
                cluster_key=payload.query_cluster_key,
            )
            query_path = Path(payload.query_clusters_tsv)

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
            custom_color_tsv=payload.custom_colors_tsv,
            restrict_obs_field=payload.restrict_obs_field,
            restrict_obs_value=payload.restrict_obs_value,
            copy_query=False,
        )

        out_dir = Path(payload.outdir)
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix_name = payload.output_prefix or f"{query_path.stem}-approximate-umap"
        prefix = out_dir / prefix_name

        export_approx_pdfs = bool(app.state.config.get("EXPORT_APPROX_PDFS", False)) or bool(payload.output_pdf)
        output_pdf: Optional[Path] = None
        if export_approx_pdfs:
            output_pdf = _resolve_output_path(
                out_dir,
                payload.output_pdf or f"{prefix_name}-comparison.pdf",
            )
            output_pdf.parent.mkdir(parents=True, exist_ok=True)

        save_h5ad = bool(payload.save_updated_h5ad)
        output_h5ad: Optional[Path] = None
        if payload.output_h5ad:
            output_h5ad = _resolve_output_path(out_dir, payload.output_h5ad)
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

    return app


app = create_app()
