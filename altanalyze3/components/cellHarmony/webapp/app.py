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
    impute_modality: Optional[str] = "none"


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


_DEFAULT_MODALITY_DEFINITIONS: Dict[str, Dict[str, object]] = {
    "rna": {
        "id": "rna",
        "label": "RNA",
        "feature_label": "gene",
        "supports_marker_heatmap": True,
        "supports_marker_network": True,
        "supports_differential_network": True,
        "supports_differential_go": True,
    },
    "lipids": {
        "id": "lipids",
        "label": "Lipids",
        "feature_label": "lipid",
        "supports_marker_heatmap": True,
        "supports_marker_network": False,
        "supports_differential_network": False,
        "supports_differential_go": False,
    },
    "adt": {
        "id": "adt",
        "label": "ADT (CITE-seq)",
        "feature_label": "ADT",
        "supports_marker_heatmap": True,
        "supports_marker_network": False,
        "supports_differential_network": False,
        "supports_differential_go": False,
    },
    "cell_communication": {
        "id": "cell_communication",
        "label": "Cell communication",
        "feature_label": "ligand-receptor interaction",
        "supports_marker_heatmap": False,
        "supports_marker_network": False,
        "supports_differential_network": False,
        "supports_differential_go": False,
    },
}


def _normalize_modality_id(value: object, *, default: str = "rna") -> str:
    raw = str(value or "").strip().lower()
    if not raw or raw in {"none", "null", "false", "off"}:
        return default
    if raw in {"lipid", "lipids"}:
        return "lipids"
    if raw in {"adt", "adts", "cite", "cite-seq", "citeseq"}:
        return "adt"
    if raw in {"cell_communication", "cell communication", "communication", "fastcomm", "fastcomm_network"}:
        return "cell_communication"
    if raw == "rna":
        return "rna"
    return raw


def _modalities_state(meta: Dict) -> Dict[str, object]:
    stored = dict(meta.get("modalities") or {})
    available_raw = stored.get("available") or []
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
    return dict(_DEFAULT_MODALITY_DEFINITIONS["rna"])


def _modality_marker_analysis(meta: Dict, modality: object) -> Dict[str, object]:
    normalized = _normalize_modality_id(modality)
    by_modality = meta.get("marker_analysis_by_modality") or {}
    if isinstance(by_modality, dict):
        entry = by_modality.get(normalized)
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
        raw_path = str(((meta.get("modality_artifacts") or {}).get(normalized) or {}).get("h5ad", "")).strip()
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
    modalities_available = list(modalities_state["available"])
    fastcomm_analysis = meta.get("fastcomm_analysis") or {}
    if isinstance(fastcomm_analysis, dict) and fastcomm_analysis.get("enabled"):
        if not any(_normalize_modality_id(entry.get("id"), default="") == "cell_communication" for entry in modalities_available):
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
    visualization_populations: Dict[str, List[str]] = {"heatmap": [], "volcano": [], "network": [], "go": [], "table": []}
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
            {"value": "volcano", "label": "Score delta"},
            {"value": "network", "label": "Cell-state network"},
            {"value": "table", "label": "Top interaction table"},
        ]
    else:
        visualization_modes = [
            {"value": "heatmap", "label": "Heatmap"},
            {"value": "volcano", "label": "Volcano"},
        ]
        if bool(modality_info.get("supports_differential_network")):
            visualization_modes.append({"value": "network", "label": "Network"})
        if bool(modality_info.get("supports_differential_go")):
            visualization_modes.append({"value": "go", "label": "GO Terms"})
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
    modality = _normalize_modality_id((meta.get("differential", {}).get("config", {}) or {}).get("modality"), default="rna")
    try:
        return _modality_h5ad_path(meta, modality)
    except FileNotFoundError:
        pass
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


def _get_marker_heatmap_cache_entry(app: FastAPI, meta: Dict, modality: str = "rna") -> Dict[str, Any]:
    job_id = str(meta.get("job_id") or "").strip()
    if not job_id:
        raise ValueError("Job metadata is missing a job_id.")

    normalized_modality = _normalize_modality_id(modality)
    marker_analysis = _modality_marker_analysis(meta, normalized_modality) or {}
    cache_path_text = str(marker_analysis.get("heatmap_cache", "")).strip()
    heatmap_tsv_text = str(marker_analysis.get("heatmap_tsv", "")).strip()
    expression_tsv_text = str(marker_analysis.get("expression_tsv", "")).strip()

    cache_path = Path(cache_path_text) if cache_path_text else None
    heatmap_tsv_path = Path(heatmap_tsv_text) if heatmap_tsv_text else None
    expression_tsv_path = Path(expression_tsv_text) if expression_tsv_text else None

    cache_signature = {
        "modality": normalized_modality,
        "heatmap_cache": str(cache_path or ""),
        "heatmap_tsv": str(heatmap_tsv_path or ""),
        "expression_tsv": str(expression_tsv_path or ""),
    }
    cache = app.state.marker_heatmap_cache
    cache_key = f"{job_id}:{normalized_modality}"
    existing = cache.get(cache_key)
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
    artifacts = meta.get("artifacts", {})
    h5ad_path = _modality_h5ad_path(meta, normalized_modality)

    cluster_key = meta.get("cluster_key")
    if not cluster_key:
        raise ValueError("Cluster assignments missing from AnnData output.")

    umap_path = artifacts.get("umap_coordinates")
    cache = app.state.expression_cache
    cache_key = f"{job_id}:{normalized_modality}"
    cache_entry = cache.get(cache_key)
    if (
        cache_entry
        and cache_entry.get("h5ad_path") == str(h5ad_path)
        and cache_entry.get("umap_path") == str(umap_path or "")
        and cache_entry.get("cluster_key") == str(cluster_key)
        and cache_entry.get("modality") == normalized_modality
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
    sample_key = str((_fastcomm_analysis(meta).get("sample_key") or "")).strip()
    if not sample_key:
        return None
    cache_entry = _get_expression_cache(app, meta, modality="rna")
    display_mask = _apply_display_filter_mask(cache_entry, display_filters)
    if not np.any(display_mask):
        raise HTTPException(status_code=404, detail="No cells match the current display filters for cell communication.")
    sample_values = (cache_entry.get("obs_filter_values") or {}).get(sample_key)
    if sample_values is None:
        adata = cache_entry.get("adata")
        if adata is not None and sample_key in adata.obs.columns:
            sample_values = (
                adata.obs[sample_key]
                .astype(str)
                .str.strip()
                .replace({"nan": "", "None": ""})
                .to_numpy(dtype=str)
            )
    if sample_values is None:
        return None
    selected = [
        value
        for value in pd.Index(np.asarray(sample_values, dtype=str)[display_mask]).unique().tolist()
        if str(value).strip()
    ]
    return sorted(dict.fromkeys(str(value).strip() for value in selected if str(value).strip()))


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
    selected_splits = _fastcomm_selected_splits(app, meta, display_filters)
    if selected_splits:
        split_scores = _prepare_fastcomm_scores(_fastcomm_split_scores_table(app, meta))
        if "split" not in split_scores.columns:
            raise ValueError("Per-sample fastComm scores missing required column: split")
        split_scores["split"] = split_scores["split"].astype(str)
        filtered_split_scores = split_scores.loc[split_scores["split"].isin(selected_splits)].copy()
        if filtered_split_scores.empty:
            raise HTTPException(status_code=404, detail="No per-sample cell-communication scores match the current display filters.")
        scores = _aggregate_fastcomm_scores(filtered_split_scores)
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
        split_scores = _prepare_fastcomm_scores(_fastcomm_split_scores_table(app, meta))
        if "split" not in split_scores.columns:
            raise ValueError("Per-sample fastComm scores missing required column: split")
        split_scores["split"] = split_scores["split"].astype(str)
        if selected_splits:
            split_scores = split_scores.loc[split_scores["split"].isin(selected_splits)].copy()
            if split_scores.empty:
                raise HTTPException(status_code=404, detail="No per-sample cell-communication scores match the current display filters.")
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
        raise ValueError("Differential comparison groups are not configured.")
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

    expression_cache = _get_expression_cache(app, meta, modality=modality)
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
    if str(population).strip() == _POOLED_OVERALL_LABEL:
        mask = sample_values.isin(resolved_group1 + resolved_group2)
    else:
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
        "modality": modality,
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


def _build_umap_payload(
    app: FastAPI,
    meta: Dict,
    modality: str = "rna",
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
) -> Dict[str, List[Dict]]:
    cache_entry = _get_expression_cache(app, meta, modality=modality)
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
    return {"reference": reference_points, "query": query_points, "sample_field": sample_field}


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


def _build_expression_payload(
    app: FastAPI,
    meta: Dict,
    gene: str,
    modality: str = "rna",
    display_filters: Optional[List[tuple[str, List[str]]]] = None,
) -> Dict:
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
    display_mask = _apply_display_filter_mask(cache_entry, display_filters)
    resolved_gene = _resolve_gene_name(cache_entry["var_names"], gene)
    if not resolved_gene:
        if len(cache_entry["var_names"]):
            fallback_gene = str(cache_entry["var_names"][0])
            values = _flatten_expr(adata[:, fallback_gene].X)
            global_min, global_max = _expression_global_range(values)
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
                "modality": normalized_modality,
                "message": None if int(np.asarray(display_mask).sum()) else "No cells match the current Display only filters.",
                "scatter": scatter_data,
                "violin": violin_data,
                "umap": umap_points,
                "global_min": global_min,
                "global_max": global_max,
            }
        if normalized_modality == "rna":
            return _build_reference_expression_payload(app, meta, gene)
        raise KeyError(f"Feature '{gene}' not found in the selected modality output.")

    values = _flatten_expr(adata[:, resolved_gene].X)
    global_min, global_max = _expression_global_range(values)
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
    global_min = float(payload.get("global_min", 0.0) or 0.0)
    global_max = float(payload.get("global_max", 0.0) or 0.0)
    if global_max <= global_min:
        global_max = global_min + 1e-9
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
        pad = max((global_max - global_min) * 0.04, 0.05)
        ax.set_ylim(global_min - pad, global_max + pad)
    else:
        umap_points = payload["umap"]
        if _normalize_modality_id(payload.get("modality"), default="rna") in {"lipids", "adt"}:
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
            vmin=global_min,
            vmax=global_max,
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
        color = str(data.get("color", "")).strip() or ("#fca5a5" if log2fc >= 0 else "#7dd3fc")
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
        meta = store.get_job(job_id)
        registry_path = Path(app.state.config["REFERENCE_REGISTRY"])
        reference_entry = pipeline_mod._lookup_reference(meta["species"], meta["reference"], registry_path)
        supported_modalities = {
            _normalize_modality_id(value, default="")
            for value in (reference_entry.get("impute_modalities") or [])
        }
        requested_modality = _normalize_modality_id(qc.impute_modality, default="")
        if requested_modality and requested_modality not in supported_modalities:
            raise HTTPException(status_code=400, detail="Selected impute modality is not supported for this reference.")
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
            return JSONResponse(_build_umap_payload(app, meta, modality=modality, display_filters=display_filters))
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
    ):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        display_filters = _display_filter_specs(filter1_field, filter1_values, filter2_field, filter2_values)
        try:
            return JSONResponse(_build_expression_payload(app, meta, gene, modality=modality, display_filters=display_filters))
        except FileNotFoundError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc))
        except ValueError as exc:
            raise HTTPException(status_code=500, detail=str(exc))

    @app.get("/api/jobs/{job_id}/marker/network")
    async def marker_network(job_id: str, population: str = Query(...), modality: str = Query("rna")):
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

    @app.get("/api/jobs/{job_id}/fastcomm/network")
    async def fastcomm_network(
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
    async def fastcomm_plot(
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

    @app.get("/api/jobs/{job_id}/genes")
    async def job_genes(job_id: str, modality: str = Query("rna")):
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

    @app.get("/api/jobs/{job_id}/display-filters")
    async def job_display_filters(job_id: str, modality: str = Query("rna")):
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
    async def umap_pdf(
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
    async def expression_pdf(
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
            matrix_entry = _get_marker_heatmap_cache_entry(app, meta, modality=modality)
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
        modality = _normalize_modality_id(request.query_params.get("modality"), default="rna")
        try:
            matrix_entry = _get_marker_heatmap_cache_entry(app, meta, modality=modality)
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
            matrix_entry = _get_marker_heatmap_cache_entry(app, meta, modality=modality)
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
        pdf = _render_marker_heatmap_pdf(frame, title="MarkerHeatmap")
        return StreamingResponse(
            pdf,
            media_type="application/pdf",
            headers={"Content-Disposition": f'attachment; filename="{job_id}_marker_heatmap.pdf"'},
        )

    @app.get("/api/jobs/{job_id}/marker/network/pdf")
    async def marker_network_pdf(job_id: str, population: str = Query(...), modality: str = Query("rna")):
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

    @app.get("/api/jobs/{job_id}/differential/interactive/table")
    async def differential_table_data(job_id: str, population: str = Query(...)):
        store, _ = _job_resources(app)
        if not store.job_exists(job_id):
            raise HTTPException(status_code=404, detail="Job not found.")
        meta = store.get_job(job_id)
        modality = _normalize_modality_id((meta.get("differential", {}) or {}).get("config", {}).get("modality"), default="rna")
        if modality != "cell_communication":
            raise HTTPException(status_code=404, detail="Differential interaction table is only available for Cell communication.")
        return JSONResponse(_build_differential_cell_communication_table_payload(app, meta, population))

    @app.get("/api/jobs/{job_id}/differential/interactive/gene")
    async def differential_gene_data(
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
    async def differential_gene_pdf(
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
