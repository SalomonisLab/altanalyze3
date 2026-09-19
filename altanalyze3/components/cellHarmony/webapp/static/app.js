"use strict";

const MAX_SAMPLES = 15;
const APP_ROOT_PATH = normalizeRootPath(window.__APP_ROOT_PATH__ || "");
let registry = window.__REFERENCE_REGISTRY__ || { species: [] };
let sampleCount = 0;
let pollTimer = null;
const VISUALIZATION_PANELS = ["viz1", "viz2"];
const VISUALIZATION_DEFAULT_MODE = {
  viz1: "cluster",
  viz2: "expression_umap",
};
let panelPlotData = {
  viz1: null,
  viz2: null,
};
let loadedResultsJobId = null;
let loadedGeneSuggestionsSignature = "";
let loadedDisplayFiltersJobId = null;
let exploreWarmupJobId = null;
let exploreWarmupPromise = null;
let exploreResultsReadyJobId = null;
let exploreResultsPendingJobId = null;
let exploreResultsReadyPromise = null;
let exploreAutoOpenPendingJobId = null;
let currentDisplayFiltersMeta = null;
let currentJobStatus = "";
let previousJobStatus = "";
let currentJobSpecies = "";
let currentJobReference = "";
let referenceRerunPending = false;
let restoringSavedSession = Boolean(new URLSearchParams(window.location.search).get("job_id"));
const markerHeatmapRenderTokens = {};
let currentMarkerAnalysis = null;
let currentMarkerAnalysisByModality = { rna: null };
let currentFastCommAnalysis = null;
let currentModalitiesState = { default: "rna", available: [{ id: "rna", label: "RNA", feature_label: "gene", example_feature: "MPO" }] };
let currentDifferentialState = null;
let differentialVisualizationRequest = 0;
let differentialDetailRequest = 0;
let currentDifferentialGene = "";
let currentDifferentialPopulation = "";
let currentDifferentialInteraction = null;
let currentDifferentialFeatureRole = "ligand";
let differentialCy = null;
// Differential gene filter. One typed gene restricts every view of the Differential
// Explorer to that gene and to the genes it interacts with in the selected cell state.
// `differentialNetworkAdjacencyCache` holds one adjacency map per (job, cell state), so
// switching between Heatmap, Volcano, Network and GO Terms costs one network fetch.
let currentDifferentialGeneFilter = "";
let differentialNetworkAdjacencyCache = {};
let expressionCyByPanel = {
  viz1: null,
  viz2: null,
};
let lastDownloadArtifactSignature = "";
let svg2PdfLoaderPromise = null;
let cytoscapeSvgLoaderPromise = null;
const BASE_VISUALIZATION_MODES = [
  { value: "cluster", label: "UMAP cell types" },
  { value: "relative", label: "UMAP broad" },
  { value: "frequency", label: "Cell frequency" },
  { value: "expression_umap", label: "UMAP" },
  { value: "violin", label: "Violin" },
  { value: "dotplot", label: "DotPlot" },
  { value: "combplot", label: "CombPlot" },
];

// DotPlot and CombPlot take a list of genes rather than one symbol, so they show
// the gene-set box instead of the single-molecule input.
const GENE_SET_MODES = new Set(["dotplot", "combplot"]);

// The backend scopes this catalog to the source study. Original feature IDs stay
// in plot coordinates, requests and click handlers; only presentation is expanded.
const featureAnnotationCache = new Map();
const featureAnnotationRequests = new Map();
async function ensureFeatureAnnotations(jobId) {
  if (!jobId || featureAnnotationCache.has(jobId)) return;
  if (!featureAnnotationRequests.has(jobId)) {
    featureAnnotationRequests.set(jobId, (async () => {
      try {
        const response = await fetch(apiPath(`/jobs/${jobId}/feature-annotations?modality=metabolite`));
        if (!response.ok) throw new Error(`Feature annotations: ${response.status}`);
        featureAnnotationCache.set(jobId, await response.json());
      } catch (error) { console.warn("Feature source annotations unavailable", error); }
      finally { featureAnnotationRequests.delete(jobId); }
    })());
  }
  await featureAnnotationRequests.get(jobId);
}
function featureAnnotation(feature, modality) {
  if (!String(modality || "").toLowerCase().startsWith("metabolite")) return null;
  return featureAnnotationCache.get(getResultsJobId())?.[String(feature)] || null;
}
function featureDisplayName(feature, modality) {
  return featureAnnotation(feature, modality)?.label || String(feature);
}
function featurePlotTitle(feature, modality) {
  const info = featureAnnotation(feature, modality);
  return info ? `${feature}<br>m/z ${info.mz} · ${info.assay}<br>RT ${info.rt_min} min` : String(feature);
}
function featureHover(feature, modality) {
  const escape = value => String(value).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const info = featureAnnotation(feature, modality);
  if (!info) return escape(feature);
  return `${escape(info.label)}<br>${escape(info.status)} · ${escape(info.source)}<br>${escape(info.source_cell)} · DOI ${escape(info.doi)}`
    + (info.formula ? `<br>Source formula: ${escape(info.formula)}` : "");
}
function differentialFeatureModality() {
  return currentDifferentialState?.config?.modality || currentDifferentialState?.modality || "rna";
}

// The Coordinates entry that means "not an embedding: use two obs columns".
const OBS_AXES_KEY = "__obs__";
// Panels that draw one point per cell on an embedding. Each takes the X and Y
// coordinate pair. "frequency" and "violin" plot no embedding, so they are out.
const UMAP_COORD_MODES = new Set(["cluster", "relative", "expression_umap"]);
const PAIRED_COLOR_STOPS = [
  [0.0, [0.6509804129600525, 0.8078431487083435, 0.8901960849761963]],
  [0.09090909090909091, [0.12156862765550613, 0.47058823704719543, 0.7058823704719543]],
  [0.18181818181818182, [0.6980392336845398, 0.8745098114013672, 0.5411764979362488]],
  [0.2727272727272727, [0.20000000298023224, 0.6274510025978088, 0.1725490242242813]],
  [0.36363636363636365, [0.9843137264251709, 0.6039215922355652, 0.6000000238418579]],
  [0.45454545454545453, [0.8901960849761963, 0.10196078568696976, 0.10980392247438431]],
  [0.5454545454545454, [0.9921568632125854, 0.7490196228027344, 0.43529412150382996]],
  [0.6363636363636364, [1.0, 0.49803921580314636, 0.0]],
  [0.7272727272727273, [0.7921568751335144, 0.6980392336845398, 0.8392156958580017]],
  [0.8181818181818182, [0.4156862795352936, 0.239215686917305, 0.6039215922355652]],
  [0.9090909090909091, [1.0, 1.0, 0.6000000238418579]],
  [1.0, [0.6941176652908325, 0.3490196168422699, 0.1568627506494522]],
];

function normalizeRootPath(rootPath) {
  const value = String(rootPath || "").trim();
  if (!value || value === "/") {
    return "";
  }
  return `/${value.replace(/^\/+|\/+$/g, "")}`;
}

function withRootPath(path) {
  if (!path) {
    return APP_ROOT_PATH || "/";
  }
  if (/^https?:\/\//i.test(path)) {
    return path;
  }
  const normalizedPath = path.startsWith("/") ? path : `/${path}`;
  return `${APP_ROOT_PATH}${normalizedPath}`;
}

function getResultsJobId() {
  return String(document.getElementById("results-job-id")?.value || "").trim();
}

function areExploreResultsReady(jobId = getResultsJobId()) {
  const normalizedJobId = String(jobId || "").trim();
  return Boolean(normalizedJobId) && exploreResultsReadyJobId === normalizedJobId;
}

function resetExploreResultsReadiness(jobId = null) {
  const normalizedJobId = String(jobId || "").trim();
  if (!normalizedJobId || exploreResultsReadyJobId === normalizedJobId) {
    exploreResultsReadyJobId = null;
  }
  if (!normalizedJobId || exploreResultsPendingJobId === normalizedJobId) {
    exploreResultsPendingJobId = null;
  }
  if (!normalizedJobId || exploreWarmupJobId === normalizedJobId) {
    exploreWarmupJobId = null;
    exploreWarmupPromise = null;
  }
  if (!normalizedJobId || exploreAutoOpenPendingJobId === normalizedJobId) {
    exploreAutoOpenPendingJobId = null;
  }
  exploreResultsReadyPromise = null;
}

async function ensureExploreResultsReady(jobId, statusData = null) {
  const normalizedJobId = String(jobId || "").trim();
  if (!normalizedJobId) {
    return false;
  }
  if (areExploreResultsReady(normalizedJobId)) {
    return true;
  }
  if (exploreResultsReadyPromise && exploreResultsPendingJobId === normalizedJobId) {
    return exploreResultsReadyPromise;
  }

  exploreResultsPendingJobId = normalizedJobId;
  exploreResultsReadyJobId = null;
  updateWorkflowPanels(referenceRerunPending ? "uploaded" : currentJobStatus);

  const readinessPromise = (async () => {
    await populateDownloadLinks(normalizedJobId, statusData);
    await loadGeneSuggestions(normalizedJobId);
    await warmExploreResults(normalizedJobId);
    if (getResultsJobId() !== normalizedJobId) {
      return false;
    }
    exploreResultsReadyJobId = normalizedJobId;
    exploreResultsPendingJobId = null;
    void loadChatExamples(normalizedJobId);
    updateWorkflowPanels(referenceRerunPending ? "uploaded" : currentJobStatus);
    if (!referenceRerunPending && statusData) {
      document.getElementById("qc-cell-status").textContent = buildQcCellSummary(statusData);
    }
    setResultMode(
      currentDifferentialState && currentDifferentialState.status === "completed"
        ? "differential"
        : "baseline"
    );
    if (
      exploreAutoOpenPendingJobId === normalizedJobId &&
      currentJobStatus === "completed" &&
      !referenceRerunPending
    ) {
      setExplorerTab("explore");
    } else {
      syncExplorerWorkspace(activeExplorerTab);
    }
    if (exploreAutoOpenPendingJobId === normalizedJobId) {
      exploreAutoOpenPendingJobId = null;
    }
    return true;
  })().catch((error) => {
    if (exploreResultsPendingJobId === normalizedJobId) {
      exploreResultsPendingJobId = null;
    }
    updateWorkflowPanels(referenceRerunPending ? "uploaded" : currentJobStatus);
    document.getElementById("qc-cell-status").textContent =
      `Unable to load expression results: ${error.message || error}. Reload to retry.`;
    const tab = document.querySelector('.workspace-tab-btn[data-tab="explore"]');
    if (tab) {
      tab.textContent = "Explore (unavailable)";
      tab.setAttribute("aria-busy", "false");
    }
    throw error;
  }).finally(() => {
    if (exploreResultsPendingJobId !== normalizedJobId) {
      exploreResultsReadyPromise = null;
    }
  });

  exploreResultsReadyPromise = readinessPromise;
  return readinessPromise;
}

function apiPath(path) {
  if (!path) {
    return withRootPath("/api");
  }
  if (path === "/api" || path.startsWith("/api/")) {
    return withRootPath(path);
  }
  const normalizedPath = path.startsWith("/") ? path : `/${path}`;
  return withRootPath(`/api${normalizedPath}`);
}

let jsPdfLoaderPromise = null;

function loadExternalScript(src) {
  return new Promise((resolve, reject) => {
    const existing = Array.from(document.scripts).find((script) => script.src === src);
    if (existing) {
      if (existing.dataset.loaded === "true") {
        resolve();
        return;
      }
      existing.addEventListener("load", () => resolve(), { once: true });
      existing.addEventListener("error", () => reject(new Error(`Failed to load ${src}`)), { once: true });
      return;
    }
    const script = document.createElement("script");
    script.src = src;
    script.async = true;
    script.crossOrigin = "anonymous";
    script.addEventListener("load", () => {
      script.dataset.loaded = "true";
      resolve();
    }, { once: true });
    script.addEventListener("error", () => reject(new Error(`Failed to load ${src}`)), { once: true });
    document.head.appendChild(script);
  });
}

async function ensureJsPdfLoaded() {
  if (window.jspdf && window.jspdf.jsPDF) {
    return window.jspdf.jsPDF;
  }
  if (!jsPdfLoaderPromise) {
    jsPdfLoaderPromise = (async () => {
      const sources = [
        "https://cdn.jsdelivr.net/npm/jspdf@2.5.1/dist/jspdf.umd.min.js",
        "https://unpkg.com/jspdf@2.5.1/dist/jspdf.umd.min.js",
      ];
      let lastError = null;
      for (const source of sources) {
        try {
          await loadExternalScript(source);
          if (window.jspdf && window.jspdf.jsPDF) {
            return window.jspdf.jsPDF;
          }
        } catch (err) {
          lastError = err;
        }
      }
      throw lastError || new Error("Unable to load jsPDF.");
    })();
  }
  return jsPdfLoaderPromise;
}

async function ensureSvg2PdfLoaded() {
  const JsPdf = await ensureJsPdfLoaded();
  if (JsPdf && JsPdf.API && typeof JsPdf.API.svg === "function") {
    return;
  }
  if (!svg2PdfLoaderPromise) {
    svg2PdfLoaderPromise = (async () => {
      const sources = [
        "https://cdn.jsdelivr.net/npm/svg2pdf.js@2.5.0/dist/svg2pdf.umd.min.js",
        "https://unpkg.com/svg2pdf.js@2.5.0/dist/svg2pdf.umd.min.js",
      ];
      let lastError = null;
      for (const source of sources) {
        try {
          await loadExternalScript(source);
          if (JsPdf && JsPdf.API && typeof JsPdf.API.svg === "function") {
            return;
          }
        } catch (err) {
          lastError = err;
        }
      }
      throw lastError || new Error("Unable to load svg2pdf.js.");
    })();
  }
  await svg2PdfLoaderPromise;
}

async function ensureCytoscapeSvgLoaded() {
  const ready = () => window.cytoscape && typeof window.cytoscape("core", "svg") === "function";
  if (ready()) {
    return;
  }
  if (!cytoscapeSvgLoaderPromise) {
    cytoscapeSvgLoaderPromise = (async () => {
      const sources = [
        "https://cdn.jsdelivr.net/npm/cytoscape-svg@0.4.0/cytoscape-svg.js",
        "https://unpkg.com/cytoscape-svg@0.4.0/cytoscape-svg.js",
      ];
      let lastError = null;
      for (const source of sources) {
        try {
          await loadExternalScript(source);
          if (window.cytoscapeSvg && window.cytoscape) {
            window.cytoscapeSvg(window.cytoscape);
          }
          if (ready()) {
            return;
          }
        } catch (err) {
          lastError = err;
        }
      }
      throw lastError || new Error("Unable to load the Cytoscape SVG exporter.");
    })();
  }
  await cytoscapeSvgLoaderPromise;
}

function slugifyFilenamePart(value, fallback = "plot") {
  const text = String(value || "")
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "_")
    .replace(/^_+|_+$/g, "");
  return text || fallback;
}

function buildPdfFilename(parts, fallback = "plot") {
  const tokens = (parts || [])
    .map((part) => slugifyFilenamePart(part, ""))
    .filter(Boolean);
  return `${tokens.length ? tokens.join("_") : fallback}.pdf`;
}

function getImageDimensions(dataUrl) {
  return new Promise((resolve, reject) => {
    const image = new Image();
    image.onload = () => {
      resolve({
        width: image.naturalWidth || image.width || 1,
        height: image.naturalHeight || image.height || 1,
      });
    };
    image.onerror = () => reject(new Error("Unable to read the exported image."));
    image.src = dataUrl;
  });
}

async function saveImageDataUrlAsPdf(dataUrl, filename) {
  const JsPdf = await ensureJsPdfLoaded();
  const { width, height } = await getImageDimensions(dataUrl);
  const orientation = width >= height ? "landscape" : "portrait";
  const pdf = new JsPdf({
    orientation,
    unit: "pt",
    format: "a4",
    compress: true,
  });
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
  pdf.save(filename);
}

function decodeSvgPayload(svgPayload) {
  const raw = String(svgPayload || "");
  if (raw.startsWith("data:image/svg+xml;base64,")) {
    return window.atob(raw.split(",")[1] || "");
  }
  if (raw.startsWith("data:image/svg+xml")) {
    return decodeURIComponent(raw.slice(raw.indexOf(",") + 1));
  }
  return raw;
}

function getSvgIntrinsicSize(svgElement) {
  const parseLength = (value) => {
    const text = String(value || "").trim();
    if (!text) {
      return null;
    }
    const numeric = Number.parseFloat(text.replace(/px$/i, ""));
    return Number.isFinite(numeric) ? numeric : null;
  };
  const width = parseLength(svgElement.getAttribute("width"));
  const height = parseLength(svgElement.getAttribute("height"));
  const viewBox = String(svgElement.getAttribute("viewBox") || "")
    .trim()
    .split(/[\s,]+/)
    .map((token) => Number.parseFloat(token));
  if (viewBox.length === 4 && Number.isFinite(viewBox[2]) && Number.isFinite(viewBox[3]) && viewBox[2] > 0 && viewBox[3] > 0) {
    return {
      width: width || viewBox[2],
      height: height || viewBox[3],
    };
  }
  return {
    width: width || 960,
    height: height || 640,
  };
}

async function saveSvgMarkupAsPdf(svgMarkup, filename, caption = "") {
  const JsPdf = await ensureJsPdfLoaded();
  await ensureSvg2PdfLoaded();
  const parser = new DOMParser();
  const documentSvg = parser.parseFromString(svgMarkup, "image/svg+xml").documentElement;
  if (!documentSvg || String(documentSvg.nodeName || "").toLowerCase() !== "svg") {
    throw new Error("The current plot could not be converted to SVG.");
  }
  documentSvg.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  documentSvg.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
  const { width, height } = getSvgIntrinsicSize(documentSvg);
  const orientation = width >= height ? "landscape" : "portrait";
  const pdf = new JsPdf({
    orientation,
    unit: "pt",
    format: "a4",
    compress: true,
  });
  const pageWidth = pdf.internal.pageSize.getWidth();
  const pageHeight = pdf.internal.pageSize.getHeight();
  const margin = 18;
  const availableWidth = pageWidth - margin * 2;
  pdf.setFontSize(9);
  const captionLines = caption ? pdf.splitTextToSize(caption.replaceAll("→", "->"), availableWidth) : [];
  const captionHeight = captionLines.length ? captionLines.length * 11 + 12 : 0;
  if (captionLines.length) pdf.text(captionLines, margin, margin + 9);
  const availableHeight = pageHeight - margin * 2 - captionHeight;
  const scale = Math.min(availableWidth / width, availableHeight / height, 1);
  const renderWidth = width * scale;
  const renderHeight = height * scale;
  const x = (pageWidth - renderWidth) / 2;
  const y = margin + captionHeight + (availableHeight - renderHeight) / 2;
  // svg2pdf reads the SVG viewport, not the width/height options below. Keep
  // the original coordinate system while fitting the whole graph on the page.
  if (!documentSvg.hasAttribute("viewBox")) documentSvg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  documentSvg.setAttribute("width", String(renderWidth));
  documentSvg.setAttribute("height", String(renderHeight));
  await pdf.svg(documentSvg, {
    x,
    y,
    width: renderWidth,
    height: renderHeight,
  });
  pdf.save(filename);
}

async function exportPlotlyElementToPdf(elementOrId, filename) {
  const element = typeof elementOrId === "string" ? document.getElementById(elementOrId) : elementOrId;
  if (!element || !element.data || !element.layout) {
    throw new Error("The current Plotly view is not available.");
  }
  const width = Math.max(Math.round(element.clientWidth || 960), 720);
  const height = Math.max(Math.round(element.clientHeight || 640), 480);
  const dataUrl = await Plotly.toImage(element, {
    format: "png",
    width,
    height,
    scale: 2,
  });
  await saveImageDataUrlAsPdf(dataUrl, filename);
}

async function exportCytoscapeToPdf(cy, filename) {
  if (!cy) {
    throw new Error("The current network view is not available.");
  }
  const dataUrl = cy.png({
    full: true,
    scale: 2,
    bg: "#ffffff",
  });
  await saveImageDataUrlAsPdf(dataUrl, filename);
}

async function exportPlotlyElementToVectorPdf(elementOrId, filename) {
  const element = typeof elementOrId === "string" ? document.getElementById(elementOrId) : elementOrId;
  if (!element || !element.data || !element.layout) {
    throw new Error("The current Plotly view is not available.");
  }
  const width = Math.max(
    Math.round(element.clientWidth || element.offsetWidth || element.layout?.width || 960),
    720,
  );
  const height = Math.max(
    Math.round(element.clientHeight || element.offsetHeight || element.layout?.height || 640) + 24,
    520,
  );
  const svgPayload = await Plotly.toImage(element, {
    format: "svg",
    width,
    height,
    scale: 1,
  });
  const svgMarkup = decodeSvgPayload(svgPayload);
  if (!svgMarkup) {
    throw new Error("The current Plotly view could not be exported.");
  }
  await saveSvgMarkupAsPdf(svgMarkup, filename);
}

async function exportCytoscapeToVectorPdf(cy, filename) {
  if (!cy) {
    throw new Error("The current network view is not available.");
  }
  await ensureCytoscapeSvgLoaded();
  if (typeof cy.svg !== "function") {
    throw new Error("The Cytoscape SVG exporter is not available.");
  }
  const svgMarkup = cy.svg({
    full: false,
    scale: 1,
    bg: "#ffffff",
  });
  await saveSvgMarkupAsPdf(svgMarkup, filename);
}

function showDownloadError(error) {
  const message = error && error.message ? error.message : "Unable to export the current plot.";
  window.alert(message);
}

function interpolatePairedColor(t) {
  const value = Math.max(0, Math.min(1, Number(t) || 0));
  for (let index = 1; index < PAIRED_COLOR_STOPS.length; index += 1) {
    const [rightPos, rightRgb] = PAIRED_COLOR_STOPS[index];
    const [leftPos, leftRgb] = PAIRED_COLOR_STOPS[index - 1];
    if (value <= rightPos || index === PAIRED_COLOR_STOPS.length - 1) {
      const span = Math.max(rightPos - leftPos, 1e-9);
      const ratio = Math.max(0, Math.min(1, (value - leftPos) / span));
      const rgb = leftRgb.map((channel, channelIndex) => {
        return channel + (rightRgb[channelIndex] - channel) * ratio;
      });
      return `rgb(${rgb.map((channel) => Math.round(channel * 255)).join(",")})`;
    }
  }
  const fallback = PAIRED_COLOR_STOPS[PAIRED_COLOR_STOPS.length - 1][1];
  return `rgb(${fallback.map((channel) => Math.round(channel * 255)).join(",")})`;
}

function customShuffleIndices(indices) {
  const shuffled = [];
  indices.forEach((value, index) => {
    if (!shuffled.includes(value)) {
      shuffled.push(value);
    }
    const fromEnd = indices[indices.length - 1 - index] ?? indices[indices.length - 1];
    if (!shuffled.includes(fromEnd)) {
      shuffled.push(fromEnd);
    }
    const fromMiddle = indices[Math.floor((index + indices.length) / 2)] ?? indices[indices.length - 1];
    if (!shuffled.includes(fromMiddle)) {
      shuffled.push(fromMiddle);
    }
  });
  return shuffled;
}

function seededShuffle(items, seed = 0) {
  let state = (seed >>> 0) || 1;
  const nextRand = () => {
    state = (1664525 * state + 1013904223) >>> 0;
    return state / 4294967296;
  };
  const values = [...items];
  for (let index = values.length - 1; index > 0; index -= 1) {
    const swapIndex = Math.floor(nextRand() * (index + 1));
    [values[index], values[swapIndex]] = [values[swapIndex], values[index]];
  }
  return values;
}

function buildReferencePreviewColorMap(populations) {
  const ordered = [...populations];
  if (ordered.length <= 4) {
    const base = ["#ff0000", "#0000ff", "#ffff00", "#00aa00", "#ffffff", "#000000", "#ff00ff"];
    return new Map(ordered.map((population, index) => [population, base[index % base.length]]));
  }

  const indices = seededShuffle(customShuffleIndices(ordered.map((_, index) => index)), 0);
  const colors = indices.map((index) => {
    const denominator = Math.max(ordered.length - 1, 1);
    return interpolatePairedColor(index / denominator);
  });
  return new Map(ordered.map((population, index) => [population, colors[index]]));
}

function buildStableUmapPopulationOrder(umapData) {
  const ordered = [];
  const seen = new Set();
  const appendPopulation = (population) => {
    const label = String(population || "").trim();
    if (!label || seen.has(label)) {
      return;
    }
    seen.add(label);
    ordered.push(label);
  };

  (umapData?.reference || []).forEach((point) => appendPopulation(point.population));
  (umapData?.query || []).forEach((point) => appendPopulation(point.population));
  return ordered;
}

function networkEdgeColor(edge) {
  const interactionType = String(edge?.data("interaction_type") || "").toLowerCase();
  if (interactionType.includes("transcription")) {
    return "#ef4444";
  }
  if (interactionType.includes("tbar")) {
    return "#60a5fa";
  }
  return "#9ca3af";
}

function networkEdgeArrowShape(edge) {
  const interactionType = String(edge?.data("interaction_type") || "").toLowerCase();
  if (interactionType.includes("tbar")) {
    return "tee";
  }
  return "triangle";
}

function panelElementId(panelKey, suffix) {
  return `${panelKey}-${suffix}`;
}

function panelPlotId(panelKey) {
  return panelElementId(panelKey, "plot");
}

function getPanelSelectValue(panelKey, suffix) {
  return String(document.getElementById(panelElementId(panelKey, suffix))?.value || "").trim();
}

function setPanelSummary(panelKey, text) {
  const summary = document.getElementById(panelElementId(panelKey, "filter-summary"));
  if (summary) {
    summary.textContent = String(text || "");
  }
}

function resetVisualizationSurface(panelKey) {
  markerHeatmapRenderTokens[panelKey] = null;
  const plot = document.getElementById(panelPlotId(panelKey));
  if (!plot) {
    return;
  }
  const cy = expressionCyByPanel[panelKey];
  if (cy) {
    cy.destroy();
    expressionCyByPanel[panelKey] = null;
  }
  plot._integratedDispose?.();
  plot.classList.remove("integrated-view");
  try {
    Plotly.purge(plot);
  } catch (_) {
    // Ignore Plotly cleanup errors when the plot is not initialized.
  }
  plot.innerHTML = "";
  plot.style.height = "";
  plot.style.minHeight = "";
  plot.style.overflowX = "";
}

function selectedReferenceConfig() {
  const species = document.getElementById("species-select")?.value;
  const reference = document.getElementById("reference-select")?.value;
  const speciesEntry = (registry.species || []).find((entry) => String(entry.id || "") === String(species || ""));
  if (!speciesEntry) {
    return null;
  }
  return (speciesEntry.references || []).find((entry) => String(entry.id || "") === String(reference || "")) || null;
}

function normalizeModalityId(value, fallback = "rna") {
  const raw = String(value || "").trim().toLowerCase();
  if (!raw || raw === "none") {
    return fallback;
  }
  if (raw === "lipids") {
    return "lipids";
  }
  if (raw === "lipid" || raw === "lipid_aml" || raw === "aml_lipid" || raw === "rna2lipid_aml") {
    return "lipid";
  }
  if (raw === "metabolite" || raw === "metabolites" || raw === "rna2metabolite") {
    return "metabolite";
  }
  if (raw === "grn" || raw === "grns" || raw === "regulon" || raw === "regulons" || raw === "gene_regulatory_network" || raw === "rna2grn" || raw === "grn_edges" || raw === "tf_edges") {
    return "grn";
  }
  if (["grn_tf", "tf_activity", "regulator_activity"].includes(raw)) return "grn_tf";
  if (raw === "adts" || raw === "cite" || raw === "cite-seq" || raw === "citeseq") {
    return "adt";
  }
  if (raw === "cell communication" || raw === "communication" || raw === "fastcomm" || raw === "fastcomm_network") {
    return "cell_communication";
  }
  return raw;
}

function availableModalities() {
  const raw = ((currentModalitiesState || {}).available || []).filter((entry) => entry && entry.id);
  if (!raw.length) {
    return [{ id: "rna", label: "RNA", feature_label: "gene", example_feature: "MPO" }];
  }
  return raw;
}

function modalityDefinition(modalityId) {
  const normalized = normalizeModalityId(modalityId);
  let defaultFeatureLabel = "gene";
  let defaultExample = "MPO";
  if (normalized === "lipids") {
    defaultFeatureLabel = "lipid";
    defaultExample = "PE(O-16:0/22:4)";
  } else if (normalized === "adt") {
    defaultFeatureLabel = "ADT";
    defaultExample = "Hu.CD4";
  } else if (normalized === "metabolite") {
    defaultFeatureLabel = "metabolite";
    defaultExample = "Lactate";
  } else if (normalized === "lipid") {
    defaultFeatureLabel = "lipid";
    defaultExample = "Hex2Cer 18:1;2O/16:0";
  } else if (normalized === "grn_tf") {
    defaultFeatureLabel = "factor";
    defaultExample = "GATA1";
  }
  return availableModalities().find((entry) => normalizeModalityId(entry.id) === normalized)
    || {
      id: normalized,
      label: normalized.toUpperCase(),
      feature_label: defaultFeatureLabel,
      example_feature: defaultExample,
    };
}

function modalityFeatureLabel(modalityId) {
  return String(modalityDefinition(modalityId).feature_label || "gene").trim().toLowerCase();
}

function modalityExampleFeature(modalityId) {
  const normalized = normalizeModalityId(modalityId);
  let fallback = "MPO";
  if (normalized === "lipids") fallback = "PE(O-16:0/22:4)";
  else if (normalized === "adt") fallback = "CD4";
  else if (normalized === "metabolite") fallback = "Lactate";
  else if (normalized === "lipid") fallback = "Hex2Cer 18:1;2O/16:0";
  else if (normalized === "grn") fallback = "GATA1";
  return String(modalityDefinition(modalityId).example_feature || fallback).trim() || fallback;
}

function preferredFeatureForModality(modalityId, features = []) {
  const values = (features || [])
    .map((value) => String(value || "").trim())
    .filter(Boolean);
  if (!values.length) {
    return modalityExampleFeature(modalityId);
  }
  const example = modalityExampleFeature(modalityId);
  if (values.includes(example)) {
    return example;
  }
  return values[0];
}

function updatePanelFeatureInput(panelKey, options = {}) {
  const input = document.getElementById(panelElementId(panelKey, "gene-query"));
  const field = document.getElementById(panelElementId(panelKey, "gene-field"));
  const label = field ? field.querySelector("span") : null;
  const modality = panelModality(panelKey);
  const featureLabel = modalityFeatureLabel(modality);
  const exampleFeature = modalityExampleFeature(modality);
  if (label) {
    label.textContent = `Select ${featureLabel}`;
  }
  if (input) {
    input.placeholder = `e.g. ${exampleFeature}`;
    if (Object.prototype.hasOwnProperty.call(options, "value")) {
      input.value = String(options.value || "");
    }
  }
}

function panelAvailableModalities() {
  // EXPLORE DROPS A MODALITY THAT DECLARES ITSELF UNBROWSABLE.
  //
  // A precomputed bundle sets supports_explore=false for a store held at cell-state
  // granularity, because every cell of a state then carries one identical value and the
  // embedding shows flat patches rather than a gradient. bundle_meta._modalities_block
  // records why. Absent the flag the modality is shown, so nothing that never sets it
  // changes.
  //
  // ONLY THIS LIST IS FILTERED. availableModalities() also backs modalityDefinition(),
  // which resolves a label and feature noun BY ID for panels including the differential
  // one, so filtering there would strip the label off a modality that is still selectable
  // elsewhere.
  return availableModalities().filter((entry) => entry.supports_explore !== false);
}

function panelModality(panelKey) {
  const mode = getPanelSelectValue(panelKey, "mode");
  if (["fastcomm_network", "cluster", "relative", "marker_network"].includes(mode)) {
    return "rna";
  }
  const select = document.getElementById(panelElementId(panelKey, "modality"));
  if (!select || select.classList.contains("hidden")) {
    return "rna";
  }
  // While mode === "fastcomm_network" this select is reused as the fastComm plot-type
  // list. Leaving that mode reads the select before it is refilled, so a stale
  // plot-type id would be treated as a modality and would drop MarkerHeatmap,
  // MarkerNetwork and Cell communication from the plot-type options.
  if (FASTCOMM_PLOT_OPTIONS.some((entry) => entry.id === String(select.value || "").trim().toLowerCase())) {
    return "rna";
  }
  return normalizeModalityId(select.value || "rna");
}

function panelCommunicationDirection(panelKey) {
  const value = String(document.getElementById(panelElementId(panelKey, "modality"))?.value || "incoming").trim().toLowerCase();
  return value === "outgoing" ? "outgoing" : "incoming";
}

const FASTCOMM_PLOT_OPTIONS = [
  { id: "focused_incoming", label: "Focused incoming" },
  { id: "focused_outgoing", label: "Focused outgoing" },
  { id: "cell_state_network", label: "Cell-state network" },
  { id: "lr_dotplot", label: "Ligand-receptor dot plot" },
  { id: "state_heatmap", label: "Cell-state heatmap" },
  { id: "top_table", label: "Top interactions table" },
  { id: "per_sample", label: "Per-sample comparison" },
];

function panelCommunicationPlotType(panelKey) {
  const value = String(document.getElementById(panelElementId(panelKey, "modality"))?.value || "focused_incoming").trim().toLowerCase();
  return FASTCOMM_PLOT_OPTIONS.some((entry) => entry.id === value) ? value : "focused_incoming";
}

function fastCommPlotNeedsPopulation(plotType) {
  return !["cell_state_network", "state_heatmap"].includes(String(plotType || "").trim().toLowerCase());
}

function panelMarkerAnalysis(panelKey) {
  const modality = panelModality(panelKey);
  const byModality = currentMarkerAnalysisByModality || {};
  return byModality[modality] || (modality === "rna" ? currentMarkerAnalysis : null) || null;
}

function markerNetworkPopulations(panelKey) {
  return (((panelMarkerAnalysis(panelKey) || {}).networks) || [])
    .map((entry) => String(entry.population || "").trim())
    .filter((value, index, values) => value && values.indexOf(value) === index);
}

function fastCommAvailable() {
  return Boolean(currentFastCommAnalysis && currentFastCommAnalysis.enabled && currentFastCommAnalysis.status === "completed");
}

function fastCommPopulations() {
  const summary = currentFastCommAnalysis?.summary || {};
  const fromSummary = currentFastCommAnalysis?.populations || summary.populations || [];
  if (Array.isArray(fromSummary) && fromSummary.length) {
    const seen = new Set();
    return fromSummary
      .map((value) => String(value || "").trim())
      .filter((value) => {
        if (!value || seen.has(value)) {
          return false;
        }
        seen.add(value);
        return true;
      });
  }
  return [];
}

function markerHeatmapAvailable(panelKey) {
  const markerAnalysis = panelMarkerAnalysis(panelKey);
  if (!markerAnalysis || !markerAnalysis.enabled) {
    return false;
  }
  return Boolean(markerAnalysis.heatmap_tsv || markerAnalysis.heatmap_cache);
}

function panelCellsPerSample(panelKey) {
  const value = Number(document.getElementById(panelElementId(panelKey, "marker-density"))?.value ?? 10);
  return [0, 5, 10, 20, 50].includes(value) ? value : 10;
}

function availableVisualizationModes(panelKey) {
  const modality = panelModality(panelKey);
  const modalityInfo = modalityDefinition(modality);
  const modes = [...BASE_VISUALIZATION_MODES];
  if (markerHeatmapAvailable(panelKey)) {
    modes.push({ value: "marker_heatmap", label: "MarkerHeatmap" });
  }
  const networkPopulations = markerNetworkPopulations(panelKey);
  if ((currentMarkerAnalysisByModality.rna || currentMarkerAnalysis)?.networks?.length) {
    modes.push({ value: "marker_network", label: "MarkerNetwork" });
  }
  if (fastCommAvailable()) {
    modes.push({ value: "fastcomm_network", label: "Cell communication" });
  }
  if (availableModalities().some(m => m.id === "grn")) {
    modes.push({ value: "grn_network", label: "GRN edges" });
  }
  modes.push({value:"integrated_network",label:"Regulatory network"});
  if (crossPathwayContexts.marker?.jobId === getResultsJobId()) modes.push({value:"integrated_cross_pathway",label:"Pathway (cross-modality)"});
  if (["lipid", "lipids", "metabolite"].includes(modality)) modes.push({value:"integrated_pathway",label:"Pathway"});
  return modes;
}

function updateExpressionModeOptions() {
  VISUALIZATION_PANELS.forEach((panelKey) => {
    const modalities = panelAvailableModalities();
    const modeSelect = document.getElementById(panelElementId(panelKey, "mode"));
    const modalityField = document.getElementById(panelElementId(panelKey, "modality-field"));
    const modalitySelect = document.getElementById(panelElementId(panelKey, "modality"));
    const markerPopulationField = document.getElementById(panelElementId(panelKey, "marker-population-field"));
    const markerPopulationSelect = document.getElementById(panelElementId(panelKey, "marker-population"));
    const geneField = document.getElementById(panelElementId(panelKey, "gene-field"));
    if (!modeSelect || !markerPopulationField || !markerPopulationSelect || !geneField || !modalityField || !modalitySelect) {
      return;
    }
    const previousModalityValue = modalitySelect.value;
    const networkPopulations = markerNetworkPopulations(panelKey);
    const fastcommPopulations = fastCommPopulations();
    const modes = availableVisualizationModes(panelKey);
    const currentMode = modeSelect.value;
    modeSelect.innerHTML = "";
    modes.forEach((mode) => {
      const option = document.createElement("option");
      option.value = mode.value;
      option.textContent = mode.label;
      modeSelect.appendChild(option);
    });
    const defaultMode = VISUALIZATION_DEFAULT_MODE[panelKey] || "cluster";
    modeSelect.value = modes.some((mode) => mode.value === currentMode) ? currentMode : defaultMode;
    const mode = modeSelect.value;

    // Gene-set plot types swap the single-molecule input for the paste box, and
    // only the CombPlot carries a per-donor minimum.
    const geneSetField = document.getElementById(panelElementId(panelKey, "geneset-field"));
    const combMinField = document.getElementById(panelElementId(panelKey, "combmin-field"));
    const wantsGeneSet = GENE_SET_MODES.has(mode);
    if (geneSetField) geneSetField.classList.toggle("hidden", !wantsGeneSet);
    const groupByField = document.getElementById(panelElementId(panelKey, "groupby-field"));
    const groupsField = document.getElementById(panelElementId(panelKey, "groups-field"));
    if (groupByField) groupByField.classList.toggle("hidden", !wantsGeneSet);
    if (groupsField) groupsField.classList.toggle("hidden", !wantsGeneSet);
    if (wantsGeneSet) refreshGroupControls(panelKey);
    const combUnitField = document.getElementById(panelElementId(panelKey, "combunit-field"));
    if (combUnitField) combUnitField.classList.toggle("hidden", mode !== "combplot");
    if (combMinField) combMinField.classList.toggle("hidden", mode !== "combplot" || panelCombUnit(panelKey) !== "donor");
    if (wantsGeneSet) geneField.classList.add("hidden");

    // Nathan, 2026-09-01: "all UMAP plots should have the option to change umap
    // coordinates". Every panel that draws cells on an embedding therefore gets
    // the X and Y lists: the cell-type view, "UMAP broad" and the expression
    // UMAP. Colour-by stays with the cell-type view, which is the only one whose
    // colour is a label rather than a measured value.
    const wantsUmapOptions = UMAP_COORD_MODES.has(mode);
    const colorByField = document.getElementById(panelElementId(panelKey, "colorby-field"));
    const coordsField = document.getElementById(panelElementId(panelKey, "coords-field"));
    if (colorByField) colorByField.classList.toggle("hidden", mode !== "cluster");
    if (coordsField) coordsField.classList.toggle("hidden", !wantsUmapOptions);
    if (wantsUmapOptions) refreshUmapOptions(panelKey);
    syncUmapAxisFields(panelKey);

    const modalityLabel = modalityField.querySelector("span");
    modalitySelect.innerHTML = "";
    if (mode === "fastcomm_network") {
      FASTCOMM_PLOT_OPTIONS.forEach((entry) => {
        const option = document.createElement("option");
        option.value = entry.id;
        option.textContent = entry.label;
        if (entry.id === previousModalityValue || (!previousModalityValue && entry.id === "focused_incoming")) {
          option.selected = true;
        }
        modalitySelect.appendChild(option);
      });
      if (!FASTCOMM_PLOT_OPTIONS.some((entry) => entry.id === modalitySelect.value)) {
        modalitySelect.value = "focused_incoming";
      }
      if (modalityLabel) {
        modalityLabel.textContent = "Plot type";
      }
      modalityField.classList.remove("hidden");
    } else {
      const currentModality = normalizeModalityId(previousModalityValue || ((currentModalitiesState || {}).default || "rna"));
      modalities.forEach((entry) => {
        const option = document.createElement("option");
        option.value = entry.id;
        option.textContent = entry.label;
        if (normalizeModalityId(entry.id) === currentModality) {
          option.selected = true;
        }
        modalitySelect.appendChild(option);
      });
      if (!modalitySelect.value && modalities.length) {
        modalitySelect.value = modalities[0].id;
      }
      if (modalityLabel) {
        modalityLabel.textContent = "Modality";
      }
      modalityField.classList.toggle("hidden", modalities.length <= 1 || ["cluster", "relative", "marker_network", "grn_network", "integrated_network", "integrated_cross_pathway"].includes(mode));
    }

    const currentPopulation = markerPopulationSelect.value;
    markerPopulationSelect.innerHTML = "";
    const fastcommPlotType = panelCommunicationPlotType(panelKey);
    const integratedPopulations = [...new Set([
      ...(currentMarkerAnalysis?.populations || []),
      ...(currentMarkerAnalysis?.networks || []).map(entry => entry.population),
      ...(currentDisplayFiltersMeta?.values?.[currentDisplayFiltersMeta?.default_secondary_field] || []),
    ].filter(Boolean))];
    const dropdownPopulations = modeSelect.value.startsWith("integrated_") ? integratedPopulations : modeSelect.value === "fastcomm_network" ? fastcommPopulations : networkPopulations;
    dropdownPopulations.forEach((population) => {
      const option = document.createElement("option");
      option.value = population;
      option.textContent = population;
      if (population === currentPopulation) {
        option.selected = true;
      }
      markerPopulationSelect.appendChild(option);
    });
    if (!markerPopulationSelect.value && dropdownPopulations.length) {
      markerPopulationSelect.value = dropdownPopulations[0];
    }

    const showMarkerPopulation =
      ((mode === "marker_network" || mode.startsWith("integrated_")) && dropdownPopulations.length > 0)
      || (mode === "fastcomm_network" && fastcommPopulations.length > 0 && fastCommPlotNeedsPopulation(fastcommPlotType));
    const showGene = mode === "expression_umap" || mode === "violin";
    updatePanelFeatureInput(panelKey);
    const populationLabel = markerPopulationField.querySelector("span");
    if (populationLabel) {
      populationLabel.textContent = mode === "fastcomm_network" ? "Marker cell state" : "Marker cell state";
    }
    markerPopulationField.classList.toggle("hidden", !showMarkerPopulation);
    geneField.classList.toggle("hidden", !showGene);
    const showGrn = mode === "grn_network";
    if (showGrn) {
      const grnGenes = document.getElementById(panelElementId(panelKey, "grn-genes"));
      const feature = document.getElementById(panelElementId(panelKey, "gene-query"))?.value || "";
      if (grnGenes && !grnGenes.value.trim() && panelModality(panelKey) === "grn_tf") grnGenes.value = feature;
    }
    ["grn-genes-field", "grn-sample-field", "grn-cellstate-field", "grn-threshold-field", "grn-limit-field", "grn-description"].forEach((suffix) => {
      const el = document.getElementById(panelElementId(panelKey, suffix));
      if (el) {
        el.classList.toggle("hidden", !showGrn);
      }
    });
    // GRN edges use their own Sample / Cell type dropdowns — the generic display filters
    // don't apply to the network, so hide them to avoid confusion.
    const filterStack = document.getElementById(panelElementId(panelKey, "filter-stack"));
    if (filterStack) {
      filterStack.classList.toggle("hidden", showGrn || mode === "marker_network" || mode.startsWith("integrated_"));
    }
    const densityRow = document.getElementById(panelElementId(panelKey, "marker-density-row"));
    if (densityRow) {
      densityRow.hidden = !(mode === "marker_heatmap" || (mode === "combplot" && panelCombUnit(panelKey) === "cells"));
    }
  });
}

function populateGrnDropdowns(panelKey, data) {
  const fill = (suffix, values, selected, allowAny = true) => {
    const sel = document.getElementById(panelElementId(panelKey, suffix));
    if (!sel) return;
    const desired = [...(allowAny ? [""] : []), ...(values || [])];
    const existing = Array.from(sel.options).map((o) => o.value);
    if (existing.length === desired.length && existing.every((v, i) => v === desired[i])) {
      if (desired.includes(selected)) sel.value = selected;
      return;
    }
    const current = selected ?? sel.value;
    sel.innerHTML = "";
    desired.forEach((v) => {
      const o = document.createElement("option");
      o.value = v;
      o.textContent = v === "" ? "All" : v;
      sel.appendChild(o);
    });
    if (desired.includes(current)) {
      sel.value = current;
    }
  };
  fill("grn-sample", data && data.available_samples, data?.sample);
  const sampleLabel = document.querySelector(`#${panelKey}-grn-sample-field > span`);
  if (sampleLabel) sampleLabel.textContent = data?.sample_field || "Sample";
  fill("grn-cellstate", data && data.available_cell_states, data?.cell_state, false);
  const input = document.getElementById(panelElementId(panelKey, "grn-genes"));
  if (input && !input.value.trim()) input.value = (data?.genes || []).join(", ");
  const suggestions = document.getElementById(`${panelKey}-grn-suggestions`);
  if (suggestions) {
    suggestions.replaceChildren(...(data?.available_factors || []).map(tf => {
      const option = document.createElement("option"); option.value = tf; return option;
    }));
  }
}

function setPanelGeneValue(panelKey, value) {
  const input = document.getElementById(panelElementId(panelKey, "gene-query"));
  if (input) {
    input.value = String(value || "");
  }
}

function buildMarkerHeatmapViewerUrl(jobId, modality) {
  const params = new URLSearchParams();
  params.set("modality", normalizeModalityId(modality));
  return `${withRootPath(`/jobs/${jobId}/marker/heatmap/viewer`)}?${params.toString()}`;
}

async function logClientEvent(jobId, message) {
  if (!jobId || !message) {
    return;
  }
  try {
    await fetch(apiPath(`/jobs/${jobId}/client-log`), {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ message }),
    });
  } catch (err) {
    console.debug("Client log failed", err);
  }
}

async function renderMarkerHeatmapViewer(jobId, panelKey) {
  resetVisualizationSurface(panelKey);
  const token = {};
  markerHeatmapRenderTokens[panelKey] = token;
  const isCurrent = () => markerHeatmapRenderTokens[panelKey] === token;
  const plot = document.getElementById(panelPlotId(panelKey));
  const params = getDisplayFilterParams(panelKey);
  params.set("modality", panelModality(panelKey));
  params.set("cells_per_sample", String(panelCellsPerSample(panelKey)));
  const suffix = params.toString() ? `?${params.toString()}` : "";
  const datasetPath = apiPath(`/jobs/${jobId}/marker/heatmap.tsv${suffix}`);
  const datasetUrl = `${window.location.origin}${datasetPath}`;
  const viewerPath = `${withRootPath(`/jobs/${jobId}/marker/heatmap/viewer`)}${suffix}`;
  const viewerUrl = `${window.location.origin}${viewerPath}`;
  await logClientEvent(jobId, `MarkerHeatmap requested. dataset_url=${datasetUrl}`);
  await logClientEvent(jobId, `MarkerHeatmap viewer_url=${viewerUrl}`);
  try {
    const resp = await fetch(datasetPath, {
      method: "HEAD",
      cache: "no-store",
    });
    await logClientEvent(
      jobId,
      `MarkerHeatmap preflight status=${resp.status} ok=${resp.ok} content_type=${resp.headers.get("content-type") || "-"} content_length=${resp.headers.get("content-length") || "-"}`
    );
    if (resp.ok) {
      const count = resp.headers.get("X-Marker-Columns");
      const field = resp.headers.get("X-Sample-Field");
      const limit = panelCellsPerSample(panelKey);
      setPanelSummary(panelKey, `${count || ""} individual cells; ${limit ? `up to ${limit} per sample per cell type` : "all cells"}. ${field ? `Sample annotation: ${field}.` : "Dataset treated as one sample."} Colours show per-gene standardized expression.`);
    }
    if (!resp.ok) {
      throw new Error(`Marker heatmap TSV returned ${resp.status}.`);
    }
  } catch (err) {
    if (!isCurrent()) return;
    await logClientEvent(jobId, `MarkerHeatmap preflight failed: ${err.message || err}`);
    renderVisualizationMessage(
      panelKey,
      `Marker heatmap TSV preflight failed.\n${err.message || err}`,
      "MarkerHeatmap"
    );
    return;
  }
  try {
    const viewerResp = await fetch(viewerPath, {
      method: "GET",
      cache: "no-store",
    });
    await logClientEvent(
      jobId,
      `MarkerHeatmap viewer preflight status=${viewerResp.status} ok=${viewerResp.ok} content_type=${viewerResp.headers.get("content-type") || "-"}`
    );
    if (!viewerResp.ok) {
      throw new Error(`Marker heatmap viewer returned ${viewerResp.status}.`);
    }
  } catch (err) {
    if (!isCurrent()) return;
    await logClientEvent(jobId, `MarkerHeatmap viewer preflight failed: ${err.message || err}`);
    renderVisualizationMessage(
      panelKey,
      `Marker heatmap viewer preflight failed.\n${err.message || err}`,
      "MarkerHeatmap"
    );
    return;
  }
  if (!isCurrent()) return;
  const iframe = document.createElement("iframe");
  iframe.className = "morpheus-frame";
  iframe.loading = "lazy";
  iframe.referrerPolicy = "no-referrer";
  iframe.addEventListener("load", () => {
    logClientEvent(jobId, `MarkerHeatmap iframe loaded src=${iframe.src}`);
  });
  iframe.addEventListener("error", () => {
    logClientEvent(jobId, `MarkerHeatmap iframe error event fired src=${iframe.src}`);
  });
  iframe.src = viewerPath;
  plot.appendChild(iframe);
}

// A missing or zero fold is neutral; it is not evidence of upregulation.
function networkFoldColor(value) {
  if (value === null || value === undefined || value === "" || !Number.isFinite(Number(value)) || Number(value) === 0) return "#cbd5e1";
  return Number(value) > 0 ? "#fca5a5" : "#7dd3fc";
}

function measurementLabel(modality) {
  const id = normalizeModalityId(modality);
  if (id === "grn_tf") return "Imputed TF activity";
  if (id === "grn") return "Predicted edge score";
  return id === "rna" ? "Expression" : "Abundance";
}

function grnNetworkDescription(payload) {
  const genes = (payload.genes || []).join(", ");
  const scope = `${payload.cell_state || "no cell type"}; ${payload.sample_field || "samples"}: ${payload.sample || "all"}`;
  return `Predicted TF → target connections involving ${genes || "the selected genes"} (${scope}). `
    + `Showing ${payload.n_edges || 0} of ${payload.n_matching_edges || 0} matching edges, ranked by |mean score|. `
    + `Scores average ${payload.n_aggregates || 0} sample × cell-type aggregates. `
    + "This view explores predicted connections; use Regulatory network to restrict targets to positive markers.";
}

function renderExpressionNetwork(panelKey, payload) {
  resetVisualizationSurface(panelKey);
  const plot = document.getElementById(panelPlotId(panelKey));
  if (plot) {
    // resetVisualizationSurface clears the container height; cytoscape needs a sized container,
    // otherwise it mounts blank (e.g. after a Plotly "no edges" message at a high threshold).
    plot.style.height = "560px";
    plot.style.minHeight = "560px";
  }
  const elements = (payload.elements || []).map((element) => {
    if (!element.data || !element.data.id || element.data.source) {
      return element;
    }
    return {
      data: {
        ...element.data,
        color: payload.node_encoding === "role"
          ? (element.data.role === "tf" ? "#facc15" : "#ef4444")
          : networkFoldColor(element.data.log2fc),
      },
    };
  });
  const maxEdgeScore = Math.max(...elements.filter(e => e.data?.source).map(e => Math.abs(Number(e.data.score) || 0)), 1e-12);
  if (!elements.length) {
    renderVisualizationMessage(panelKey, `No marker network was available for ${payload.population}.`, `${payload.population} marker network`);
    return;
  }
  expressionCyByPanel[panelKey] = cytoscape({
    container: plot,
    elements,
    style: [
      {
        selector: "node",
        style: {
          "background-color": "data(color)",
          label: "data(label)",
          color: "#0f172a",
          "font-size": (node) => payload.node_encoding === "role" ? (node.data("role") === "tf" ? 20 : 14) : 12,
          "font-weight": (node) => payload.node_encoding === "role" && node.data("role") === "tf" ? "bold" : "normal",
          "text-valign": "center",
          "text-halign": "center",
          width: (node) => payload.node_encoding === "role" && node.data("role") === "tf" ? 52 : 26,
          height: (node) => payload.node_encoding === "role" && node.data("role") === "tf" ? 52 : 26,
          "border-width": (node) => payload.node_encoding === "role" && node.data("queried") ? 3 : 0,
          "border-color": "#0f172a",
          shape: (node) => payload.node_encoding === "role" && node.data("role") === "tf" ? "diamond" : "ellipse",
        },
      },
      {
        selector: "edge",
        style: {
          // GRN widths are relative to this graph; raw scores remain available on hover.
          width: (edge) => {
            const s = Math.abs(Number(edge.data("score")) || 0);
            if (payload.node_encoding === "role") return 1 + 4 * s / maxEdgeScore;
            return s > 0 ? Math.max(1, Math.min(9, 1 + s * 24)) : 1.8;
          },
          "line-color": payload.node_encoding === "role" ? "#94a3b8" : networkEdgeColor,
          "target-arrow-color": payload.node_encoding === "role" ? "#94a3b8" : networkEdgeColor,
          "target-arrow-shape": networkEdgeArrowShape,
          "curve-style": "bezier",
          opacity: 0.85,
        },
      },
      {
        selector: "node:selected",
        style: {
          "border-width": 3,
          "border-color": "#0f172a",
        },
      },
    ],
    layout: payload.node_encoding === "role" ? {
      name: "concentric", animate: false, fit: true, padding: 45,
      concentric: (node) => node.data("queried") ? 2 : 1,
      levelWidth: () => 1, minNodeSpacing: 10, avoidOverlap: true,
      nodeDimensionsIncludeLabels: true,
    } : {
      name: "cose",
      animate: false,
      fit: true,
      padding: 36,
      randomize: true,
      idealEdgeLength: 80,
      nodeOverlap: 8,
      componentSpacing: 90,
    },
  });
  expressionCyByPanel[panelKey].on("tap", "node", (event) => {
    const gene = event?.target?.data("id");
    if (gene) {
      setPanelGeneValue(panelKey, gene);
    }
  });
  if (payload.node_encoding === "role") {
    expressionCyByPanel[panelKey].on("mouseover", "node, edge", (event) => {
      const item = event.target;
      const text = item.isNode() ? `${item.id()} — ${item.data("role") === "tf" ? "transcription factor" : "target"}`
        : `${item.data("source")} → ${item.data("target")}; predicted edge score ${item.data("score")}`;
      setNetworkHoverTooltip(panelKey, text, event.renderedPosition);
    });
    expressionCyByPanel[panelKey].on("mouseout", "node, edge", () => setNetworkHoverTooltip(panelKey, ""));
  }
}

function setNetworkHoverTooltip(panelKey, text, renderedPosition = null) {
  const plot = document.getElementById(panelPlotId(panelKey));
  if (!plot) {
    return;
  }
  let tooltip = plot.querySelector(".network-hover-tooltip");
  if (!tooltip) {
    tooltip = document.createElement("div");
    tooltip.className = "network-hover-tooltip";
    plot.appendChild(tooltip);
  }
  const message = String(text || "").trim();
  if (!message) {
    tooltip.classList.add("hidden");
    tooltip.textContent = "";
    return;
  }
  tooltip.textContent = message;
  tooltip.classList.remove("hidden");
  if (renderedPosition) {
    tooltip.style.left = `${Math.min(plot.clientWidth - 260, Math.max(12, renderedPosition.x + 12))}px`;
    tooltip.style.top = `${Math.min(plot.clientHeight - 120, Math.max(12, renderedPosition.y + 12))}px`;
  }
}

function renderFastCommNetwork(panelKey, payload) {
  resetVisualizationSurface(panelKey);
  const plot = document.getElementById(panelPlotId(panelKey));
  const elements = payload.elements || [];
  if (!elements.length) {
    renderVisualizationMessage(panelKey, payload.message || `No fastComm interactions were available for ${payload.population}.`, "fastComm");
    return;
  }
  expressionCyByPanel[panelKey] = cytoscape({
    container: plot,
    elements,
    style: [
      {
        selector: "node",
        style: {
          "background-color": "data(color)",
          label: "data(label)",
          color: "#0f172a",
          "font-size": 12,
          "text-valign": "center",
          "text-halign": "center",
          width: 34,
          height: 34,
          "border-width": 1,
          "border-color": "#94a3b8",
        },
      },
      {
        selector: "node[node_type = 'focus']",
        style: {
          width: 52,
          height: 52,
          color: "#ffffff",
          "font-weight": 700,
          "border-width": 3,
          "border-color": "#134e4a",
        },
      },
      {
        selector: "edge",
        style: {
          width: "data(weight)",
          "line-color": "#0f766e",
          "target-arrow-color": "#0f766e",
          "target-arrow-shape": "triangle",
          "curve-style": "bezier",
          opacity: "data(edge_opacity)",
        },
      },
      {
        selector: "edge.hovered",
        style: {
          opacity: 1,
          "line-color": "#f97316",
          "target-arrow-color": "#f97316",
          label: "data(label)",
          color: "#0f172a",
          "font-size": 10,
          "text-background-color": "#ffffff",
          "text-background-opacity": 0.85,
          "text-background-padding": 3,
        },
      },
      {
        selector: "node:selected",
        style: {
          "border-width": 4,
          "border-color": "#f97316",
        },
      },
    ],
    layout: {
      name: "circle",
      animate: true,
      fit: true,
      padding: 56,
    },
  });
  expressionCyByPanel[panelKey].on("mouseover", "edge", (event) => {
    event.target.addClass("hovered");
    setNetworkHoverTooltip(panelKey, event.target.data("tooltip"), event.renderedPosition);
  });
  expressionCyByPanel[panelKey].on("mouseout", "edge", (event) => {
    event.target.removeClass("hovered");
    setNetworkHoverTooltip(panelKey, "");
  });
  expressionCyByPanel[panelKey].on("tap", "edge", (event) => {
    event.target.addClass("hovered");
    setNetworkHoverTooltip(panelKey, event.target.data("tooltip"), event.renderedPosition);
  });
  expressionCyByPanel[panelKey].on("tap", "node", (event) => {
    const label = event?.target?.data("label");
    if (label) {
      const directionLabel = payload.direction === "global" ? "Global" : (payload.direction === "outgoing" ? "Outgoing" : "Incoming");
      setPanelSummary(panelKey, `${directionLabel} cell communication view${payload.population ? ` for ${payload.population}` : ""}; selected ${label}.`);
    }
  });
}

function renderFastCommHeatmap(panelKey, payload) {
  resetVisualizationSurface(panelKey);
  const senders = payload.senders || [];
  const receivers = payload.receivers || [];
  if (!senders.length || !receivers.length) {
    renderVisualizationMessage(panelKey, payload.message || "No cell-state communication matrix is available.", "Cell communication heatmap");
    return;
  }
  Plotly.newPlot(panelPlotId(panelKey), [{
    type: "heatmap",
    x: receivers,
    y: senders,
    z: payload.z || [],
    text: payload.text || [],
    hovertemplate: "%{text}<extra></extra>",
    colorscale: [
      [0, "#eff6ff"],
      [0.4, "#38bdf8"],
      [0.75, "#0f766e"],
      [1, "#7c2d12"],
    ],
    colorbar: { title: "Summed score" },
  }], {
    title: "Cell-state communication strength",
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.94)",
    height: Math.max(520, senders.length * 18 + 130),
    margin: { t: 48, l: 150, r: 28, b: 150 },
    xaxis: { title: "Receiver state", automargin: true },
    yaxis: { title: "Sender state", automargin: true },
  }, { responsive: true });
}

function renderFastCommDotPlot(panelKey, payload) {
  resetVisualizationSurface(panelKey);
  const points = payload.points || [];
  if (!points.length) {
    renderVisualizationMessage(panelKey, payload.message || "No ligand-receptor points are available.", "Ligand-receptor dot plot");
    return;
  }
  const maxScore = Math.max(...points.map((point) => Number(point.score) || 0), 1e-9);
  Plotly.newPlot(panelPlotId(panelKey), [{
    type: "scatter",
    mode: "markers",
    x: points.map((point) => point.x),
    y: points.map((point) => point.y),
    text: points.map((point) => (
      `${point.sender_state} -> ${point.receiver_state}<br>` +
      `${point.y}<br>` +
      `score=${Number(point.score || 0).toFixed(3)}<br>` +
      `LR=${Number(point.lr_expression_score || 0).toFixed(3)}; response=${Number(point.receiver_response_score || 0).toFixed(3)}`
    )),
    marker: {
      size: points.map((point) => 8 + 26 * Math.sqrt((Number(point.score) || 0) / maxScore)),
      color: points.map((point) => Number(point.receiver_response_score) || 0),
      colorscale: [
        [0, "#dbeafe"],
        [0.45, "#14b8a6"],
        [1, "#f97316"],
      ],
      showscale: true,
      colorbar: { title: "Response score" },
      line: { color: "#0f172a", width: 0.5 },
      opacity: 0.82,
    },
    hovertemplate: "%{text}<extra></extra>",
  }], {
    title: `${payload.direction === "outgoing" ? "Outgoing" : "Incoming"} ligand-receptor evidence for ${payload.population}`,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.94)",
    height: Math.max(520, Math.min(980, points.length * 18 + 140)),
    margin: { t: 48, l: 170, r: 40, b: 130 },
    xaxis: { title: payload.direction === "outgoing" ? "Receiver state" : "Sender state", automargin: true },
    yaxis: { title: "Ligand -> receptor", automargin: true },
  }, { responsive: true });
}

function renderFastCommTable(panelKey, payload) {
  resetVisualizationSurface(panelKey);
  const rows = payload.rows || [];
  const columns = payload.columns || [];
  if (!rows.length || !columns.length) {
    renderVisualizationMessage(panelKey, payload.message || "No significant interactions are available.", "Top interactions");
    return;
  }
  const labels = columns.map((column) => column
    .replace("fastcomm_score", "score")
    .replace("lr_expression_score_scaled", "LR expression")
    .replace("receiver_response_score", "response")
    .replaceAll("_", " "));
  const values = columns.map((column) => rows.map((row) => {
    const value = row[column];
    return Number.isFinite(Number(value)) && String(value).trim() !== "" ? Number(value).toFixed(3) : String(value ?? "");
  }));
  Plotly.newPlot(panelPlotId(panelKey), [{
    type: "table",
    header: {
      values: labels,
      align: "left",
      fill: { color: "#0f766e" },
      font: { color: "white", size: 11 },
    },
    cells: {
      values,
      align: "left",
      fill: { color: rows.map((_, index) => (index % 2 ? "#f8fafc" : "#ffffff")) },
      font: { color: "#0f172a", size: 10 },
      height: 24,
    },
  }], {
    title: payload.direction === "global" ? "Top significant cell-communication interactions" : `Top significant interactions for ${payload.population}`,
    paper_bgcolor: "rgba(0,0,0,0)",
    height: Math.max(520, Math.min(900, rows.length * 26 + 110)),
    margin: { t: 46, l: 12, r: 12, b: 12 },
  }, { responsive: true });
}

function renderFastCommPerSample(panelKey, payload) {
  resetVisualizationSurface(panelKey);
  const rows = payload.rows || [];
  if (!rows.length) {
    renderVisualizationMessage(panelKey, payload.message || "Per-sample communication scores are unavailable.", "Per-sample communication");
    return;
  }
  const states = [...new Set(rows.map((row) => row.state))];
  const traces = states.map((state) => {
    const stateRows = rows.filter((row) => row.state === state);
    return {
      type: "bar",
      name: state,
      x: stateRows.map((row) => row.sample),
      y: stateRows.map((row) => row.total_score),
      text: stateRows.map((row) => `${state}<br>${row.sample}<br>total=${Number(row.total_score || 0).toFixed(3)}<br>interactions=${row.n_interactions}`),
      hovertemplate: "%{text}<extra></extra>",
    };
  });
  Plotly.newPlot(panelPlotId(panelKey), traces, {
    title: `${payload.direction === "outgoing" ? "Outgoing" : "Incoming"} per-sample communication for ${payload.population}`,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.94)",
    height: 560,
    margin: { t: 48, l: 64, r: 24, b: 120 },
    barmode: "group",
    xaxis: { title: payload.sample_key || "sample", automargin: true },
    yaxis: { title: "Summed communication score" },
    legend: { orientation: "h", y: -0.26 },
  }, { responsive: true });
}

function renderFastCommPlot(panelKey, payload) {
  const plotType = payload.plot_type || panelCommunicationPlotType(panelKey);
  if (plotType === "state_heatmap") {
    renderFastCommHeatmap(panelKey, payload);
    return;
  }
  if (plotType === "lr_dotplot") {
    renderFastCommDotPlot(panelKey, payload);
    return;
  }
  if (plotType === "top_table") {
    renderFastCommTable(panelKey, payload);
    return;
  }
  if (plotType === "per_sample") {
    renderFastCommPerSample(panelKey, payload);
    return;
  }
  renderFastCommNetwork(panelKey, payload);
}

function finiteExtent(values, fallbackMin = 0, fallbackMax = 1) {
  let min = Number.POSITIVE_INFINITY;
  let max = Number.NEGATIVE_INFINITY;
  (values || []).forEach((value) => {
    const numeric = Number(value);
    if (!Number.isFinite(numeric)) {
      return;
    }
    if (numeric < min) {
      min = numeric;
    }
    if (numeric > max) {
      max = numeric;
    }
  });
  if (!Number.isFinite(min) || !Number.isFinite(max)) {
    return [fallbackMin, fallbackMax];
  }
  return [min, max];
}

function relaxReferencePreviewLabels(labels, points, plotElement) {
  if (!Array.isArray(labels) || labels.length <= 1) {
    return labels || [];
  }

  const xValues = (points || []).map((point) => Number(point.x)).filter(Number.isFinite);
  const yValues = (points || []).map((point) => Number(point.y)).filter(Number.isFinite);
  labels.forEach((label) => {
    if (Number.isFinite(Number(label.x))) {
      xValues.push(Number(label.x));
    }
    if (Number.isFinite(Number(label.y))) {
      yValues.push(Number(label.y));
    }
  });

  const [xMin, xMax] = finiteExtent(xValues, 0, 1);
  const [yMin, yMax] = finiteExtent(yValues, 0, 1);
  const xRange = Math.max(xMax - xMin, 1);
  const yRange = Math.max(yMax - yMin, 1);
  const plotWidth = Math.max((plotElement?.clientWidth || 900) - 80, 480);
  const plotHeight = 520;
  const padX = (8 / plotWidth) * xRange;
  const padY = (5 / plotHeight) * yRange;
  const maxDx = xRange * 0.03;
  const maxDy = yRange * 0.03;

  const relaxed = labels.map((label) => {
    const text = String(label.population || "");
    const widthPx = Math.max(42, Math.min(170, text.length * 6.2));
    const heightPx = 18;
    return {
      ...label,
      x: Number(label.x),
      y: Number(label.y),
      anchorX: Number(label.x),
      anchorY: Number(label.y),
      halfWidth: (widthPx / plotWidth) * xRange * 0.5,
      halfHeight: (heightPx / plotHeight) * yRange * 0.5,
    };
  });

  for (let iteration = 0; iteration < 110; iteration += 1) {
    for (let index = 0; index < relaxed.length; index += 1) {
      const current = relaxed[index];
      for (let compareIndex = index + 1; compareIndex < relaxed.length; compareIndex += 1) {
        const other = relaxed[compareIndex];
        const dx = other.x - current.x;
        const dy = other.y - current.y;
        const overlapX = current.halfWidth + other.halfWidth + padX - Math.abs(dx);
        const overlapY = current.halfHeight + other.halfHeight + padY - Math.abs(dy);
        if (overlapX <= 0 || overlapY <= 0) {
          continue;
        }

        let pushX = overlapX * 0.35;
        let pushY = overlapY * 0.35;
        if (Math.abs(dx) < 1e-6) {
          pushX *= index % 2 === 0 ? -1 : 1;
        } else {
          pushX *= dx > 0 ? -1 : 1;
        }
        if (Math.abs(dy) < 1e-6) {
          pushY *= index % 2 === 0 ? -1 : 1;
        } else {
          pushY *= dy > 0 ? -1 : 1;
        }

        current.x += pushX;
        other.x -= pushX;
        current.y += pushY;
        other.y -= pushY;
      }
    }

    relaxed.forEach((label) => {
      label.x += (label.anchorX - label.x) * 0.14;
      label.y += (label.anchorY - label.y) * 0.14;
      label.x = Math.max(label.anchorX - maxDx, Math.min(label.anchorX + maxDx, label.x));
      label.y = Math.max(label.anchorY - maxDy, Math.min(label.anchorY + maxDy, label.y));
    });
  }

  return relaxed;
}

function buildPopulationCentroids(points) {
  const grouped = new Map();
  (points || []).forEach((point) => {
    const population = String(point.population || "").trim();
    const x = Number(point.x);
    const y = Number(point.y);
    if (!population || !Number.isFinite(x) || !Number.isFinite(y)) {
      return;
    }
    if (!grouped.has(population)) {
      grouped.set(population, { population, xs: [], ys: [] });
    }
    const entry = grouped.get(population);
    entry.xs.push(x);
    entry.ys.push(y);
  });
  return Array.from(grouped.values()).map((entry) => ({
    population: entry.population,
    x: median(entry.xs),
    y: median(entry.ys),
  }));
}

function buildSquareUmapAxes(points, paddingFraction = 0.08) {
  const finitePoints = (points || []).filter(
    (point) => Number.isFinite(Number(point?.x)) && Number.isFinite(Number(point?.y))
  );
  if (!finitePoints.length) {
    return {};
  }
  const xs = finitePoints.map((point) => Number(point.x));
  const ys = finitePoints.map((point) => Number(point.y));
  const [minX, maxX] = finiteExtent(xs, 0, 1);
  const [minY, maxY] = finiteExtent(ys, 0, 1);
  const xSpan = Math.max(maxX - minX, 1);
  const ySpan = Math.max(maxY - minY, 1);
  const xPad = xSpan * paddingFraction;
  const yPad = ySpan * Math.max(0.02, paddingFraction * 0.45);
  return {
    xaxis: {
      range: [minX - xPad, maxX + xPad],
      showgrid: false,
      zeroline: false,
      showticklabels: false,
      ticks: "",
    },
    yaxis: {
      range: [minY - yPad, maxY + yPad],
      showgrid: false,
      zeroline: false,
      showticklabels: false,
      ticks: "",
    },
  };
}

function median(values) {
  if (!Array.isArray(values) || !values.length) {
    return 0;
  }
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  if (sorted.length % 2 === 0) {
    return (sorted[middle - 1] + sorted[middle]) / 2;
  }
  return sorted[middle];
}

document.addEventListener("DOMContentLoaded", () => {
  initExplorerSandbox();
  initSpeciesSelect();
  initSampleRows();
  hookForms();
  updateWorkflowPanels(null);
  updateDifferentialUi(null);
  updateResetDataButton();
  loadReferencePreview();
  initWindowCount();
  initGeneSetBoxes();
  initChatTab();
  restoreJobFromUrl();
});

async function restoreJobFromUrl() {
  const jobId = new URLSearchParams(window.location.search).get("job_id");
  if (!jobId) return;
  try {
    const response = await fetch(apiPath(`/jobs/${encodeURIComponent(jobId)}/status`), {cache:"no-store"});
    const data = await parseApiResponse(response);
    if (!response.ok) throw new Error(data.detail || "Unable to reopen this job.");
    document.getElementById("species-select").value = data.species;
    document.getElementById("species-select").dispatchEvent(new Event("change"));
    document.getElementById("reference-select").value = data.reference;
    document.getElementById("reference-select").dispatchEvent(new Event("change"));
    for (const [name, value] of Object.entries(data.qc || {})) {
      const field = document.querySelector(`#qc-form [name="${name}"]`);
      if (field && !Array.isArray(value)) field.value = value == null ? "" : String(value);
    }
    const imputed = data.qc?.impute_modalities || [];
    if (imputed.length) document.getElementById("qc-impute-modality-select").value = imputed.includes("all") ? "all" : imputed[0];
    applyJobStatus(jobId, data);
    if (["queued", "processing"].includes(data.status) ||
        ["queued", "processing"].includes(data.differential_ui?.status)) {
      startStatusPolling(jobId);
    }
    if (data.status === "completed") {
      await ensureExploreResultsReady(jobId, data);
      setExplorerTab("explore");
    }
  } catch (error) {
    alert(error.message || "Unable to reopen this job.");
  } finally {
    restoringSavedSession = false;
    if (!currentJobStatus) loadReferencePreview();
  }
}

function initSpeciesSelect() {
  const speciesSelect = document.getElementById("species-select");
  const referenceSelect = document.getElementById("reference-select");

  registry.species.forEach((sp) => {
    const option = document.createElement("option");
    option.value = sp.id;
    option.textContent = sp.label;
    speciesSelect.appendChild(option);
  });

  speciesSelect.addEventListener("change", () => {
    const selected = registry.species.find((sp) => sp.id === speciesSelect.value);
    referenceSelect.innerHTML = "";
    if (!selected) {
      updateImputeModalityField();
      return;
    }
    selected.references.forEach((ref) => {
      const option = document.createElement("option");
      option.value = ref.id;
      option.textContent = ref.label;
      referenceSelect.appendChild(option);
    });
    updateReferenceChangeState();
    referenceSelect.dispatchEvent(new Event("change"));
  });

  referenceSelect.addEventListener("change", () => {
    updateImputeModalityField();
    updateReferenceChangeState();
    loadReferencePreview();
  });

  if (registry.species.length) {
    speciesSelect.value = registry.species[0].id;
    speciesSelect.dispatchEvent(new Event("change"));
  }
}

function updateImputeModalityField() {
  const field = document.getElementById("qc-impute-modality-field");
  const select = document.getElementById("qc-impute-modality-select");
  if (!field || !select) {
    return;
  }
  const reference = selectedReferenceConfig();
  const supported = Array.from(new Set(((reference && reference.impute_modalities) || []).map((value) => normalizeModalityId(value, "")))).filter(Boolean);
  field.classList.toggle("hidden", supported.length === 0);
  if (!supported.length) {
    select.value = "none";
    return;
  }
  Array.from(select.options).forEach((option) => {
    if (option.value === "none" || option.value === "all") {
      option.hidden = false;
      option.disabled = false;
      return;
    }
    const enabled = supported.includes(normalizeModalityId(option.value, ""));
    option.hidden = !enabled;
    option.disabled = !enabled;
  });
  const current = select.value;
  if (current !== "none" && current !== "all" && !supported.includes(normalizeModalityId(current, ""))) {
    select.value = "none";
  }
}

function selectedReferenceDiffersFromLoadedJob() {
  const species = document.getElementById("species-select").value;
  const reference = document.getElementById("reference-select").value;
  if (!currentJobSpecies || !currentJobReference) {
    return false;
  }
  return species !== currentJobSpecies || reference !== currentJobReference;
}

function updateReferenceChangeState() {
  referenceRerunPending = Boolean(
    document.getElementById("qc-job-id").value.trim()
      && currentJobStatus === "completed"
      && selectedReferenceDiffersFromLoadedJob()
  );

  if (referenceRerunPending) {
    updateWorkflowPanels("uploaded");
    document.getElementById("differential-panel").classList.add("hidden");
    document.getElementById("baseline-results-view").classList.add("hidden");
    document.getElementById("differential-results-view").classList.add("hidden");
    clearGeneSuggestions();
    clearDisplayFilters();
    resetDifferentialResults();
    currentMarkerAnalysisByModality = { rna: null };
    currentFastCommAnalysis = null;
    document.getElementById("qc-cell-status").textContent =
      "Reference changed. Click Save QC and run to realign the uploaded data against the newly selected reference.";
    return;
  }

  updateWorkflowPanels(currentJobStatus);
  if (currentJobStatus === "completed" && currentDifferentialState) {
    updateDifferentialUi(currentDifferentialState);
  }
}

function initSampleRows() {
  document.getElementById("add-sample-btn").addEventListener("click", () => {
    if (sampleCount >= MAX_SAMPLES) {
      alert(`You can only upload ${MAX_SAMPLES} files per job.`);
      return;
    }
    addSampleRow();
  });
  addSampleRow();
}

function addSampleRow() {
  const container = document.getElementById("sample-container");
  const row = document.createElement("div");
  row.className = "sample-row";
  row.innerHTML = `
    <label class="field sample-inline-field">
      <span class="sample-inline-label">Sample name</span>
      <input type="text" class="sample-name" placeholder="e.g. Sample_${sampleCount + 1}" required>
    </label>
    <label class="field sample-inline-field">
      <span class="sample-inline-label">H5/H5AD file</span>
      <input type="file" class="sample-file" accept=".h5,.h5ad" required>
    </label>
    <div class="sample-actions">
      <button type="button" class="ghost-btn remove-sample">Remove</button>
    </div>
  `;
  container.appendChild(row);
  row.querySelector(".remove-sample").addEventListener("click", () => {
    row.remove();
    sampleCount = Math.max(0, sampleCount - 1);
  });
  sampleCount += 1;
}

function hookForms() {
  document.getElementById("job-form").addEventListener("submit", handleJobSubmit);
  document.getElementById("qc-form").addEventListener("submit", handleQcSubmit);
  document.getElementById("differential-form").addEventListener("submit", handleDifferentialSubmit);
  document.getElementById("results-form").addEventListener("submit", handleResultsSubmit);
  document.getElementById("reset-data-btn").addEventListener("click", resetWorkspaceData);
  document.getElementById("results-job-id").addEventListener("change", handleResultsJobChange);
  VISUALIZATION_PANELS.forEach((panelKey) => {
    document.getElementById(panelElementId(panelKey, "mode")).addEventListener("change", () => {
      updateExpressionModeOptions();
      loadVisualizationPanel(panelKey);
    });
    document.getElementById(panelElementId(panelKey, "modality")).addEventListener("change", () => {
      const datalist = document.getElementById(`${panelKey}-feature-suggestions`);
      if (datalist) {
        datalist.innerHTML = "";
      }
      updatePanelFeatureInput(panelKey, { value: "" });
      loadedGeneSuggestionsSignature = "";
      updateExpressionModeOptions();
      void (async () => {
        await loadGeneSuggestions(getResultsJobId());
        loadVisualizationPanel(panelKey);
      })();
    });
    document.getElementById(panelElementId(panelKey, "gene-query")).addEventListener("change", () => {
      loadVisualizationPanel(panelKey);
    });
    document.getElementById(panelElementId(panelKey, "marker-population")).addEventListener("change", () => {
      loadVisualizationPanel(panelKey);
    });
    const densitySelect = document.getElementById(panelElementId(panelKey, "marker-density"));
    if (densitySelect) {
      densitySelect.addEventListener("change", () => loadVisualizationPanel(panelKey));
    }
    ["grn-genes", "grn-sample", "grn-cellstate", "grn-limit"].forEach((suffix) => {
      const el = document.getElementById(panelElementId(panelKey, suffix));
      if (el) {
        el.addEventListener("change", () => loadVisualizationPanel(panelKey));
      }
    });
    const grnThreshold = document.getElementById(panelElementId(panelKey, "grn-threshold"));
    if (grnThreshold) {
      grnThreshold.addEventListener("input", () => {
        const lbl = document.getElementById(panelElementId(panelKey, "grn-threshold-value"));
        if (lbl) {
          lbl.textContent = Number(grnThreshold.value).toFixed(3);
        }
      });
      grnThreshold.addEventListener("change", () => loadVisualizationPanel(panelKey));
    }
    document.getElementById(panelElementId(panelKey, "filter1-field")).addEventListener("change", () => {
      syncDisplayFilterValueOptions(panelKey, 1);
      loadVisualizationPanel(panelKey);
    });
    document.getElementById(panelElementId(panelKey, "filter2-field")).addEventListener("change", () => {
      syncDisplayFilterValueOptions(panelKey, 2);
      loadVisualizationPanel(panelKey);
    });
    document.getElementById(panelElementId(panelKey, "filter1-values")).addEventListener("change", () => {
      loadVisualizationPanel(panelKey);
    });
    document.getElementById(panelElementId(panelKey, "filter2-values")).addEventListener("change", () => {
      loadVisualizationPanel(panelKey);
    });
    document.getElementById(`download-${panelKey}-image-btn`).addEventListener("click", (event) => {
      event.preventDefault();
      downloadVisualizationImage(panelKey);
    });
  });
  document.getElementById("plot-dot-scale").addEventListener("change", () => {
    VISUALIZATION_PANELS.forEach((panelKey) => {
      renderVisualizationPanel(panelKey);
    });
  });
  document.getElementById("download-differential-left-btn").addEventListener("click", (event) => {
    event.preventDefault();
    downloadDifferentialLeftPdf();
  });
  document.getElementById("download-differential-gene-btn").addEventListener("click", (event) => {
    event.preventDefault();
    downloadDifferentialGenePdf();
  });
  document.getElementById("differential-viz-mode").addEventListener("change", () => {
    syncDifferentialPopulationSelect(currentDifferentialState);
    updateDifferentialDownloadButton();
    loadDifferentialVisualization();
  });
  document.getElementById("differential-population").addEventListener("change", () => {
    if (!currentDifferentialState) {
      return;
    }
    const nextState = {
      ...currentDifferentialState,
      config: {
        ...(currentDifferentialState.config || {}),
        population_col: document.getElementById("differential-population").value,
      },
    };
    currentDifferentialState = markDifferentialConfigDirty(nextState);
    updateDifferentialUi(currentDifferentialState);
  });
  document.getElementById("differential-modality").addEventListener("change", async () => {
    if (!currentDifferentialState) return;
    clearDifferentialGeneFilter();
    const modality = normalizeModalityId(document.getElementById("differential-modality").value || "rna");
    const config = {...currentDifferentialState.config, modality};
    const matching = (currentDifferentialState.completed_comparisons || []).filter(entry => {
      if (entry.modality !== modality) return false;
      let identity;
      try { identity = JSON.parse(entry.contrast); } catch (_) { return false; }
      return ["population_col", "sample_field", "comparison_type"].every(k => (identity[k] || "") === (config[k] || "")) &&
        ["group1_samples", "group2_samples"].every(k => JSON.stringify([...(identity[k] || [])].sort()) === JSON.stringify([...(config[k] || [])].sort()));
    });
    if (matching.length) {
      await selectCompletedDifferential(matching[matching.length - 1].id);
    } else {
      ++savedComparisonSelectionRequest;
      currentDifferentialState = markDifferentialConfigDirty({...currentDifferentialState, config});
      updateDifferentialUi(currentDifferentialState);
    }
  });
  document.getElementById("differential-sample-field").addEventListener("change", () => {
    if (!currentDifferentialState) {
      return;
    }
    currentDifferentialState = markDifferentialConfigDirty({
      ...currentDifferentialState,
      config: {
        ...(currentDifferentialState.config || {}),
        sample_field: document.getElementById("differential-sample-field").value,
        group1_samples: [],
        group2_samples: [],
      },
    });
    updateDifferentialUi(currentDifferentialState);
  });
  document.getElementById("differential-result-population").addEventListener("change", () => {
    currentDifferentialGene = "";
    currentDifferentialInteraction = null;
    loadDifferentialVisualization();
  });
  // The gene filter. `change` fires when a datalist suggestion is picked or the field is
  // committed; `search` fires when the native clear button of a type=search empties it.
  const geneFilterInput = differentialGeneFilterInput();
  if (geneFilterInput) {
    const applyGeneFilter = () => {
      const next = String(geneFilterInput.value || "").trim();
      if (next === currentDifferentialGeneFilter) {
        return;
      }
      currentDifferentialGeneFilter = next;
      loadDifferentialVisualization();
    };
    geneFilterInput.addEventListener("change", applyGeneFilter);
    geneFilterInput.addEventListener("search", applyGeneFilter);
    geneFilterInput.addEventListener("keydown", (event) => {
      if (event.key === "Enter") {
        event.preventDefault();
        applyGeneFilter();
      }
    });
  }
  initToggleMultiSelect(document.getElementById("differential-group1"));
  initToggleMultiSelect(document.getElementById("differential-group2"));
  attachDifferentialLrToggleHandlers();
}

let savedComparisonSelectionRequest = 0;
async function selectCompletedDifferential(contrast) {
  const request = ++savedComparisonSelectionRequest;
  ++differentialVisualizationRequest;
  ++differentialDetailRequest;
  try {
    const response = await fetch(apiPath(`/jobs/${getResultsJobId()}/differential/select?contrast=${encodeURIComponent(contrast)}`), {method: "POST"});
    const body = await response.json();
    if (request !== savedComparisonSelectionRequest) return;
    if (!response.ok) throw new Error(body.detail || "Comparison selection failed");
    updateDifferentialUi(body);
    setResultMode("differential");
  } catch (error) {
    if (request === savedComparisonSelectionRequest) document.getElementById("differential-message").textContent = error.message;
  }
}

function markDifferentialConfigDirty(state) {
  if (!state || state.status !== "completed") {
    return state;
  }
  return {
    ...state,
    status: "idle",
    message: "Differential settings changed. Run cellHarmony-differential to analyze the updated configuration.",
  };
}

function updateResetDataButton() {
  const button = document.getElementById("reset-data-btn");
  if (!button) {
    return;
  }
  const hasDataset = Boolean(document.getElementById("upload-job-id").value.trim() || currentJobStatus);
  button.classList.toggle("hidden", !hasDataset);
}

function resetWorkspaceData() {
  setSessionUrl(null);
  if (pollTimer) {
    clearInterval(pollTimer);
    pollTimer = null;
  }

  panelPlotData = {
    viz1: null,
    viz2: null,
  };
  loadedResultsJobId = null;
  loadedGeneSuggestionsSignature = "";
  loadedDisplayFiltersJobId = null;
  resetExploreResultsReadiness();
  currentDisplayFiltersMeta = null;
  previousJobStatus = "";
  currentJobStatus = "";
  currentJobSpecies = "";
  currentJobReference = "";
  referenceRerunPending = false;
  currentMarkerAnalysis = null;
  currentMarkerAnalysisByModality = { rna: null };
  currentFastCommAnalysis = null;
  currentModalitiesState = { default: "rna", available: [{ id: "rna", label: "RNA", feature_label: "gene", example_feature: "MPO" }] };
  currentDifferentialState = null;
  currentDifferentialGene = "";
  currentDifferentialPopulation = "";
  lastDownloadArtifactSignature = "";

  if (differentialCy) {
    differentialCy.destroy();
    differentialCy = null;
  }
  VISUALIZATION_PANELS.forEach((panelKey) => {
    const cy = expressionCyByPanel[panelKey];
    if (cy) {
      cy.destroy();
      expressionCyByPanel[panelKey] = null;
    }
  });

  document.getElementById("upload-job-id").value = "";
  document.getElementById("qc-job-id").value = "";
  document.getElementById("results-job-id").value = "";
  const uploadJobIdDisplay = document.getElementById("upload-job-id-display");
  uploadJobIdDisplay.textContent = "";
  uploadJobIdDisplay.classList.add("hidden");

  document.getElementById("download-links").innerHTML = "";
  document.getElementById("job-log").textContent = "";
  document.getElementById("qc-cell-status").textContent = "QC cell counts will appear here while the analysis runs.";
  document.getElementById("differential-message").textContent =
    "Differential analysis is enabled when the job contains two or more samples.";
  document.getElementById("differential-archive-link").classList.add("hidden");

  setUploadProgress(0, false);
  document.getElementById("job-progress").style.width = "0%";
  document.getElementById("job-progress-label").textContent = "0%";
  document.getElementById("differential-progress").style.width = "0%";
  document.getElementById("differential-progress-label").textContent = "0%";

  clearGeneSuggestions();
  clearDisplayFilters();
  updateExpressionModeOptions();
  resetDifferentialResults();

  VISUALIZATION_PANELS.forEach((panelKey) => {
    resetVisualizationSurface(panelKey);
    setPanelSummary(panelKey, "");
    setPanelGeneValue(panelKey, "");
  });

  const previewPlot = document.getElementById("reference-preview-plot");
  try {
    Plotly.purge(previewPlot);
  } catch (_) {
    // Ignore Plotly cleanup errors when the plot is not initialized.
  }
  previewPlot.innerHTML = "";

  const qcLivePlot = document.getElementById("qc-live-plot");
  try {
    Plotly.purge(qcLivePlot);
  } catch (_) {
    // Ignore Plotly cleanup errors when the plot is not initialized.
  }
  qcLivePlot.innerHTML = "";
  const qcLiveCaption = document.getElementById("qc-live-caption");
  if (qcLiveCaption) {
    qcLiveCaption.textContent = "Live counts parsed from alignment log.";
  }
  const ambientLivePlot = document.getElementById("ambient-live-plot");
  if (ambientLivePlot) {
    try {
      Plotly.purge(ambientLivePlot);
    } catch (_) {
      // Ignore Plotly cleanup errors when the plot is not initialized.
    }
    ambientLivePlot.innerHTML = "";
  }
  const ambientLiveCaption = document.getElementById("ambient-live-caption");
  if (ambientLiveCaption) {
    ambientLiveCaption.textContent = "Ambient RNA correction percentages will appear here when correction is performed.";
  }

  const sampleContainer = document.getElementById("sample-container");
  sampleContainer.innerHTML = "";
  sampleCount = 0;
  addSampleRow();
  updateImputeModalityField();

  setExplorerTab("run");
  updateWorkflowPanels(null);
  updateDifferentialUi(null);
  updateResetDataButton();
  loadReferencePreview();
}

async function parseApiResponse(resp) {
  const rawText = await resp.text();
  if (!rawText) {
    return {};
  }
  try {
    return JSON.parse(rawText);
  } catch (_) {
    return { detail: rawText.trim() || `HTTP ${resp.status}` };
  }
}

function setUploadProgress(percent, visible = true) {
  const shell = document.getElementById("upload-progress-shell");
  const fill = document.getElementById("upload-progress");
  const label = document.getElementById("upload-progress-label");
  const normalized = Math.max(0, Math.min(100, Math.round(Number(percent) || 0)));
  shell.classList.toggle("hidden", !visible);
  fill.style.width = `${normalized}%`;
  label.textContent = `${normalized}%`;
}

function uploadFormDataWithProgress(url, formData, onProgress) {
  return new Promise((resolve, reject) => {
    const xhr = new XMLHttpRequest();
    xhr.open("POST", url);
    xhr.responseType = "text";

    xhr.upload.onprogress = (event) => {
      if (!event.lengthComputable || !onProgress) {
        return;
      }
      onProgress((event.loaded / event.total) * 100);
    };

    xhr.onload = () => {
      const rawText = xhr.responseText || "";
      let data = {};
      if (rawText) {
        try {
          data = JSON.parse(rawText);
        } catch (_) {
          data = { detail: rawText.trim() || `HTTP ${xhr.status}` };
        }
      }
      resolve({ ok: xhr.status >= 200 && xhr.status < 300, status: xhr.status, data });
    };

    xhr.onerror = () => reject(new Error("Upload failed."));
    xhr.onabort = () => reject(new Error("Upload cancelled."));
    xhr.send(formData);
  });
}

async function handleJobSubmit(evt) {
  evt.preventDefault();
  const species = document.getElementById("species-select").value;
  const reference = document.getElementById("reference-select").value;
  const formData = new FormData();
  formData.append("species", species);
  formData.append("reference", reference);

  const rows = document.querySelectorAll("#sample-container .sample-row");
  if (!rows.length) {
    alert("Add at least one sample.");
    return;
  }

  for (const row of rows) {
    const nameInput = row.querySelector(".sample-name");
    const fileInput = row.querySelector(".sample-file");
    if (!nameInput.value || !fileInput.files.length) {
      alert("Each row needs both a sample name and a file.");
      return;
    }
    formData.append("sample_names", nameInput.value.trim());
    formData.append("files", fileInput.files[0]);
  }

  const submitBtn = document.getElementById("upload-submit-btn");
  try {
    submitBtn.disabled = true;
    setUploadProgress(0, true);
    const resp = await uploadFormDataWithProgress(apiPath("/jobs"), formData, (percent) => {
      setUploadProgress(percent, true);
    });
    if (!resp.ok) {
      throw new Error(resp.data.detail || "Failed to create job.");
    }
    const data = resp.data;
    setSessionUrl(data.job_id);
    setUploadProgress(100, true);
    document.getElementById("upload-job-id").value = data.job_id;
    const uploadJobIdDisplay = document.getElementById("upload-job-id-display");
    uploadJobIdDisplay.textContent = data.job_id;
    uploadJobIdDisplay.classList.remove("hidden");
    document.getElementById("qc-job-id").value = data.job_id;
    document.getElementById("results-job-id").value = data.job_id;
    loadedResultsJobId = null;
    resetExploreResultsReadiness();
    setResultMode("baseline");
    updateResetDataButton();
    await loadJobState(data.job_id);
  } catch (err) {
    setUploadProgress(0, false);
    alert(err.message);
  } finally {
    submitBtn.disabled = false;
  }
}

async function handleQcSubmit(evt) {
  evt.preventDefault();
  const jobId = document.getElementById("qc-job-id").value.trim();
  if (!jobId) {
    alert("Enter a job id.");
    return;
  }
  const payload = {
    min_genes: evt.target.min_genes.value,
    min_counts: evt.target.min_counts.value,
    min_cells: evt.target.min_cells.value,
    mit_percent: evt.target.mit_percent.value,
    align_cutoff: evt.target.align_cutoff.value,
    ambient_correction: evt.target.ambient_correction.value,
    impute_modalities: (() => {
      const value = evt.target.impute_modality ? evt.target.impute_modality.value : "none";
      return value && value !== "none" ? [value] : [];   // "all" expands on the backend
    })(),
  };
  try {
    if (selectedReferenceDiffersFromLoadedJob()) {
      const configureResp = await fetch(apiPath(`/jobs/${jobId}/configure`), {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          species: document.getElementById("species-select").value,
          reference: document.getElementById("reference-select").value,
        }),
      });
      const configureData = await parseApiResponse(configureResp);
      if (!configureResp.ok) {
        throw new Error(configureData.detail || "Failed to update the job reference.");
      }
      currentJobSpecies = configureData.species || currentJobSpecies;
      currentJobReference = configureData.reference || currentJobReference;
      referenceRerunPending = false;
      currentDifferentialState = null;
      lastDownloadArtifactSignature = "";
      panelPlotData = {
        viz1: null,
        viz2: null,
      };
      currentMarkerAnalysis = null;
      currentFastCommAnalysis = null;
      loadedResultsJobId = null;
      resetExploreResultsReadiness();
      clearGeneSuggestions();
      clearDisplayFilters();
      document.getElementById("download-links").innerHTML = "";
      resetDifferentialResults();
      document.getElementById("differential-panel").classList.add("hidden");
      updateExpressionModeOptions();
      updateWorkflowPanels(configureData.status || "uploaded");
    }

    let resp = await fetch(apiPath(`/jobs/${jobId}/qc`), {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    let data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "Failed to save QC.");
    }
    resp = await fetch(apiPath(`/jobs/${jobId}/run`), { method: "POST" });
    data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "Failed to queue job.");
    }
    loadedResultsJobId = null;
    resetExploreResultsReadiness(jobId);
    setResultMode("baseline");
    startStatusPolling(jobId);
  } catch (err) {
    alert(err.message);
  }
}

async function handleDifferentialSubmit(evt) {
  evt.preventDefault();
  const jobId = document.getElementById("results-job-id").value.trim();
  if (!jobId) {
    alert("Enter a completed job id in Section 3 first.");
    return;
  }
  const payload = {
    modality: document.getElementById("differential-modality").value || "rna",
    population_col: document.getElementById("differential-population").value,
    sample_field: document.getElementById("differential-sample-field").value,
    group1_samples: getMultiSelectValues(document.getElementById("differential-group1")),
    group2_samples: getMultiSelectValues(document.getElementById("differential-group2")),
    comparison_type: document.getElementById("differential-comparison-type").value || "cells",
  };

  try {
    const resp = await fetch(apiPath(`/jobs/${jobId}/differential`), {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const data = await parseApiResponse(resp);
    if (!resp.ok) {
      throw new Error(data.detail || "Failed to start differential analysis.");
    }
    updateDifferentialUi(data);
    setExplorerTab("differential");
    startStatusPolling(jobId);
  } catch (err) {
    alert(err.message);
  }
}

async function handleResultsJobChange() {
  const jobId = document.getElementById("results-job-id").value.trim();
  if (!jobId) {
    return;
  }
  loadedResultsJobId = null;
  resetExploreResultsReadiness();
  await loadJobState(jobId);
}

function handleResultsSubmit(evt) {
  evt.preventDefault();
  refreshResults();
}

async function loadJobState(jobId) {
  try {
    const resp = await fetch(apiPath(`/jobs/${jobId}/status?t=${Date.now()}`), { cache: "no-store" });
    const data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "Status request failed.");
    }
    applyJobStatus(jobId, data);
  } catch (err) {
    alert(err.message);
  }
}

function startStatusPolling(jobId) {
  if (pollTimer) {
    clearInterval(pollTimer);
  }
  pollTimer = setInterval(() => pollStatus(jobId), 2000);
  pollStatus(jobId);
}

async function pollStatus(jobId) {
  try {
    const resp = await fetch(apiPath(`/jobs/${jobId}/status?t=${Date.now()}`), { cache: "no-store" });
    const data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "Status request failed.");
    }
    applyJobStatus(jobId, data);
  } catch (err) {
    console.warn(err);
  }
}

function setSessionUrl(jobId) {
  if (window.__SCALABLE_VIEWER__) return;
  const url = new URL(window.location.href);
  if (jobId) url.searchParams.set("job_id", jobId);
  else url.searchParams.delete("job_id");
  if (url.href !== window.location.href) window.history.replaceState(null, "", url);
}

function applyJobStatus(jobId, data) {
  setSessionUrl(jobId);
  previousJobStatus = currentJobStatus;
  currentJobStatus = String(data.status || "").trim().toLowerCase();
  currentJobSpecies = String(data.species || currentJobSpecies || "");
  currentJobReference = String(data.reference || currentJobReference || "");
  currentMarkerAnalysis = data.marker_analysis && data.marker_analysis.enabled ? data.marker_analysis : null;
  currentMarkerAnalysisByModality = data.marker_analysis_by_modality || { rna: currentMarkerAnalysis };
  currentFastCommAnalysis = data.fastcomm_analysis && data.fastcomm_analysis.enabled ? data.fastcomm_analysis : null;
  currentModalitiesState = data.modalities || { default: "rna", available: [{ id: "rna", label: "RNA", feature_label: "gene", example_feature: "MPO" }] };
  updateExpressionModeOptions();
  updateWorkflowPanels(referenceRerunPending ? "uploaded" : (data.status || null));
  document.getElementById("job-progress").style.width = `${data.progress || 0}%`;
  document.getElementById("job-progress-label").textContent = `${data.progress || 0}%`;
  document.getElementById("job-log").textContent = formatPanelLogTail(data.log_head || [], data.log_tail || []);
  renderQcLiveProgress(data);
  if (!referenceRerunPending) {
    document.getElementById("qc-cell-status").textContent = buildQcCellSummary(data);
  }

  if (!document.getElementById("upload-job-id").value) {
    document.getElementById("upload-job-id").value = jobId;
    const uploadJobIdDisplay = document.getElementById("upload-job-id-display");
    uploadJobIdDisplay.textContent = jobId;
    uploadJobIdDisplay.classList.remove("hidden");
  }
  if (!document.getElementById("qc-job-id").value) {
    document.getElementById("qc-job-id").value = jobId;
  }
  if (!document.getElementById("results-job-id").value) {
    document.getElementById("results-job-id").value = jobId;
  }
  updateResetDataButton();

  updateDifferentialUi(data.differential_ui || null);
  updateReferenceChangeState();

  const alignmentJustCompleted =
    previousJobStatus !== "completed" &&
    currentJobStatus === "completed" &&
    !referenceRerunPending;

  if (alignmentJustCompleted) {
    exploreAutoOpenPendingJobId = jobId;
  }

  if (data.status === "completed" && !referenceRerunPending) {
    const artifactSignature = JSON.stringify(Object.keys(data.artifacts || {}).sort());
    if (artifactSignature !== lastDownloadArtifactSignature) {
      lastDownloadArtifactSignature = artifactSignature;
    }
    if (data.default_gene) {
      VISUALIZATION_PANELS.forEach((panelKey) => {
        const mode = getPanelSelectValue(panelKey, "mode");
        if (!requiresGeneMode(mode)) {
          return;
        }
        const input = document.getElementById(panelElementId(panelKey, "gene-query"));
        if (input && !String(input.value || "").trim()) {
          const modality = panelModality(panelKey);
          input.value = normalizeModalityId(modality) === "rna"
            ? data.default_gene
            : preferredFeatureForModality(modality);
        }
      });
    }
    void ensureExploreResultsReady(jobId, data).catch((err) => {
      console.warn(err);
    });
  } else {
    resetExploreResultsReadiness(jobId);
    currentMarkerAnalysis = null;
    currentMarkerAnalysisByModality = { rna: null };
    currentFastCommAnalysis = null;
    currentModalitiesState = data.modalities || { default: "rna", available: [{ id: "rna", label: "RNA", feature_label: "gene", example_feature: "MPO" }] };
    updateExpressionModeOptions();
    clearDisplayFilters();
  }

  const pipelineActive = data.status === "queued" || data.status === "processing";
  const differentialStatus = (data.differential_ui || {}).status;
  const differentialActive = differentialStatus === "queued" || differentialStatus === "processing";
  if (!pipelineActive && !differentialActive && pollTimer) {
    clearInterval(pollTimer);
    pollTimer = null;
  }
}

function updateWorkflowPanels(status) {
  const normalizedStatus = String(status || "").trim().toLowerCase();
  const hasUploadedJob = Boolean(normalizedStatus);
  const hasStartedQcRun = hasUploadedJob && normalizedStatus !== "uploaded";
  const hasCompletedAlignment = normalizedStatus === "completed";
  const hasExploreReady = hasCompletedAlignment && areExploreResultsReady();
  const runGrid = document.getElementById("sandbox-run-grid");
  const qcPanel = document.getElementById("qc-panel");
  const qcLivePanel = document.getElementById("qc-live-panel");
  const differentialPanel = document.getElementById("differential-panel");
  const exploreControlsPanel = document.getElementById("explore-controls-panel");
  const previewPanel = document.getElementById("reference-preview-panel");
  const baseline = document.getElementById("baseline-results-view");
  const differential = document.getElementById("differential-results-view");
  const exploreTabButton = document.querySelector('.workspace-tab-btn[data-tab="explore"]');
  const differentialTabButton = document.querySelector('.workspace-tab-btn[data-tab="differential"]');

  qcPanel.classList.toggle("hidden", !hasUploadedJob);
  qcLivePanel.classList.toggle("hidden", !hasStartedQcRun);
  differentialPanel.classList.toggle("hidden", !hasCompletedAlignment);
  exploreControlsPanel.classList.toggle("hidden", !hasExploreReady);
  previewPanel.classList.toggle("hidden", hasUploadedJob);
  if (runGrid) {
    runGrid.classList.toggle("preupload-mode", !hasUploadedJob);
    runGrid.classList.toggle("workflow-mode", hasUploadedJob);
  }
  if (exploreTabButton) {
    // Large saved sessions need time to load their expression matrices. Keep
    // the tab visible so loading does not look like missing saved results.
    exploreTabButton.classList.toggle("hidden", !hasCompletedAlignment);
    exploreTabButton.disabled = !hasExploreReady;
    exploreTabButton.textContent = hasCompletedAlignment && !hasExploreReady
      ? "Explore (loading…)" : "Explore";
    exploreTabButton.setAttribute("aria-busy", String(hasCompletedAlignment && !hasExploreReady));
  }
  if (differentialTabButton) {
    differentialTabButton.classList.toggle("hidden", !hasCompletedAlignment);
    differentialTabButton.disabled = !hasCompletedAlignment;
  }

  if (!hasCompletedAlignment) {
    baseline.classList.add("hidden");
    differential.classList.add("hidden");
  } else if (!hasExploreReady) {
    baseline.classList.add("hidden");
  }
  syncExplorerWorkspace(activeExplorerTab);
}

async function loadReferencePreview() {
  if (restoringSavedSession) return;
  const previewPanel = document.getElementById("reference-preview-panel");
  if (currentJobStatus) {
    previewPanel.classList.add("hidden");
    return;
  }

  const species = document.getElementById("species-select").value;
  const reference = document.getElementById("reference-select").value;
  const plot = document.getElementById("reference-preview-plot");
  if (!species || !reference) {
    Plotly.purge(plot);
    plot.innerHTML = "";
    return;
  }

  try {
    const resp = await fetch(
      apiPath(`/meta/reference-preview?species=${encodeURIComponent(species)}&reference=${encodeURIComponent(reference)}`)
    );
    const payload = await parseApiResponse(resp);
    if (!resp.ok) {
      throw new Error(payload.detail || "Reference preview request failed.");
    }
    renderReferencePreview(payload);
  } catch (err) {
    Plotly.purge(plot);
    plot.innerHTML = `<div class="empty-state">${err.message || "Unable to load reference preview."}</div>`;
  }
}

function renderReferencePreview(payload) {
  const plot = document.getElementById("reference-preview-plot");
  const populations = [...new Set((payload.points || []).map((point) => point.population))];
  const colorMap = buildReferencePreviewColorMap(populations);
  const relaxedLabels = relaxReferencePreviewLabels(payload.labels || [], payload.points || [], plot);

  const pointTrace = {
    type: "scattergl",
    mode: "markers",
    x: (payload.points || []).map((point) => point.x),
    y: (payload.points || []).map((point) => point.y),
    text: (payload.points || []).map((point) => point.population),
    hovertemplate: "%{text}<extra></extra>",
    marker: {
      size: 2,
      opacity: 0.5,
      color: (payload.points || []).map((point) => colorMap.get(point.population)),
    },
    showlegend: false,
  };

  Plotly.newPlot(
    plot,
    [pointTrace],
    {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      margin: { t: 8, l: 8, r: 8, b: 8 },
      height: 640,
      annotations: relaxedLabels.map((label) => ({
        x: label.x,
        y: label.y,
        text: label.population,
        showarrow: false,
        xref: "x",
        yref: "y",
        xanchor: "center",
        yanchor: "middle",
        font: {
          size: 11,
          color: "#0f172a",
        },
        bgcolor: "rgba(255,255,255,0)",
        opacity: 1,
      })),
      xaxis: { visible: false, showgrid: false, zeroline: false },
      yaxis: { visible: false, showgrid: false, zeroline: false, scaleanchor: "x", scaleratio: 1 },
      hovermode: "closest",
    },
    { responsive: true, displayModeBar: false }
  );
}

function parseProgressPercent(value) {
  if (typeof value === "number") {
    return Number.isFinite(value) ? value : 0;
  }
  if (typeof value === "string") {
    const parsed = Number.parseFloat(value.replace("%", "").trim());
    return Number.isFinite(parsed) ? parsed : 0;
  }
  return 0;
}

function getStatusLogLines(data) {
  const head = Array.isArray(data && data.log_head) ? data.log_head : [];
  const tail = Array.isArray(data && data.log_tail) ? data.log_tail : [];
  return head.concat(tail).map((line) => String(line || ""));
}

function formatPanelLogTail(headLines, tailLines) {
  const head = Array.isArray(headLines) ? headLines : [];
  const tail = Array.isArray(tailLines) ? tailLines : [];
  if (!head.length && !tail.length) {
    return "";
  }
  const noisyPatterns = [
    /\[mem\]\[timing\]/i,
    /Loading input files:\s+\d+%/i,
    /Downsampling clusters:\s+\d+%/i,
    /Computing marker statistics:\s+\d+%/i,
  ];
  const cleanHead = head.filter((line) => {
    const value = String(line || "");
    return !noisyPatterns.some((pattern) => pattern.test(value));
  });
  const filteredTail = tail.filter((line) => {
    const value = String(line || "");
    return !noisyPatterns.some((pattern) => pattern.test(value));
  });
  const selectedTail = (filteredTail.length ? filteredTail : tail).slice(-120);
  if (!cleanHead.length) {
    return selectedTail.join("");
  }
  const seen = new Set(cleanHead.map((line) => String(line)));
  const dedupedTail = selectedTail.filter((line) => !seen.has(String(line)));
  return cleanHead.concat(dedupedTail).join("");
}

function extractQcThresholdState(lines) {
  const state = {
    total: null,
    minGenesThreshold: null,
    afterMinGenes: null,
    minCountsThreshold: null,
    afterMinCounts: null,
    afterMito: null,
  };
  const values = Array.isArray(lines) ? lines : [];
  for (const rawLine of values) {
    const line = String(rawLine || "");
    let match = line.match(/(?:reimported\s+)?adata shape:\s*\((\d+),/i);
    if (match) {
      state.total = Number(match[1]);
      continue;
    }
    match = line.match(/Cells remaining after min_genes\s+([0-9.]+)\s+filtering:\s*(\d+)/i);
    if (match) {
      state.minGenesThreshold = Number(match[1]);
      state.afterMinGenes = Number(match[2]);
      continue;
    }
    match = line.match(/Cells remaining after min_counts\s+([0-9.]+)\s+filtering:\s*(\d+)/i);
    if (match) {
      state.minCountsThreshold = Number(match[1]);
      state.afterMinCounts = Number(match[2]);
      continue;
    }
    match = line.match(/Cells remaining after mito-percent filtering:\s*(\d+)/i);
    if (match) {
      state.afterMito = Number(match[1]);
    }
  }
  if (!Number.isFinite(state.total)) {
    const fallbackTotal = [state.afterMito, state.afterMinCounts, state.afterMinGenes].find((value) => Number.isFinite(value));
    if (Number.isFinite(fallbackTotal)) {
      state.total = fallbackTotal;
    }
  }
  return state;
}

function getQcThresholdInput(name) {
  const input = document.querySelector(`#qc-form [name="${name}"]`);
  const parsed = Number.parseFloat(String((input && input.value) || "").trim());
  return Number.isFinite(parsed) ? parsed : null;
}

function extractAmbientCorrectionState(lines) {
  const rows = [];
  const seen = new Map();
  for (const line of lines) {
    const match = String(line || "").match(/Auto-selected rho for library '(.+?)':\s*([0-9.]+)/i);
    if (!match) {
      continue;
    }
    const library = String(match[1] || "").trim();
    const rho = Number(match[2]);
    if (!library || !Number.isFinite(rho)) {
      continue;
    }
    seen.set(library, Math.max(0, Math.min(1, rho)));
  }
  for (const [library, rho] of seen.entries()) {
    rows.push({ library, rho });
  }
  rows.sort((a, b) => b.rho - a.rho || a.library.localeCompare(b.library));
  return rows;
}

function clearPlotEmptyState(plot) {
  if (!plot) {
    return;
  }
  if (plot.querySelector(".empty-state")) {
    plot.innerHTML = "";
  }
}

function renderQcLiveProgress(data) {
  const plot = document.getElementById("qc-live-plot");
  const caption = document.getElementById("qc-live-caption");
  const ambientPlot = document.getElementById("ambient-live-plot");
  const ambientCaption = document.getElementById("ambient-live-caption");
  if (!plot || !caption) {
    return;
  }
  const lines = getStatusLogLines(data);
  const state = extractQcThresholdState(lines);
  const ambientRows = extractAmbientCorrectionState(lines);
  const total = Number.isFinite(state.total) ? state.total : null;

  if (!Number.isFinite(total) || total <= 0) {
    try {
      Plotly.purge(plot);
    } catch (_) {
      // Ignore cleanup errors when no Plotly instance exists yet.
    }
    plot.innerHTML = "<div class=\"empty-state\">QC metrics will appear here as filters are applied.</div>";
    caption.textContent = "Live counts parsed from alignment log.";
    renderAmbientCorrectionProgress(ambientPlot, ambientCaption, ambientRows, data);
    return;
  }

  const minGenesThreshold = Number.isFinite(state.minGenesThreshold)
    ? state.minGenesThreshold
    : getQcThresholdInput("min_genes");
  const minCountsThreshold = Number.isFinite(state.minCountsThreshold)
    ? state.minCountsThreshold
    : getQcThresholdInput("min_counts");
  const mitoThreshold = getQcThresholdInput("mit_percent");

  const rows = [];
  if (Number.isFinite(state.afterMinGenes)) {
    const label = Number.isFinite(minGenesThreshold) ? `min genes \u2265 ${minGenesThreshold}` : "min genes filter";
    rows.push({ label, kept: state.afterMinGenes });
  }
  if (Number.isFinite(state.afterMinCounts)) {
    const label = Number.isFinite(minCountsThreshold) ? `min counts \u2265 ${minCountsThreshold}` : "min counts filter";
    rows.push({ label, kept: state.afterMinCounts });
  }
  if (Number.isFinite(state.afterMito)) {
    const label = Number.isFinite(mitoThreshold) ? `mito % \u2264 ${mitoThreshold}` : "mito % filter";
    rows.push({ label, kept: state.afterMito });
  }

  if (!rows.length) {
    try {
      Plotly.purge(plot);
    } catch (_) {
      // Ignore cleanup errors when no Plotly instance exists yet.
    }
    plot.innerHTML = "<div class=\"empty-state\">QC metrics will appear here as filters are applied.</div>";
    caption.textContent = `Total detected cells: ${total.toLocaleString()}.`;
    renderAmbientCorrectionProgress(ambientPlot, ambientCaption, ambientRows, data);
    return;
  }

  const labels = rows.map((row) => row.label);
  const kept = rows.map((row) => Math.max(0, Math.min(total, Number(row.kept) || 0)));
  const filtered = kept.map((value) => Math.max(0, total - value));
  const maxLabelLength = labels.reduce((max, value) => Math.max(max, String(value).length), 0);
  const marginLeft = Math.min(210, Math.max(120, 18 + maxLabelLength * 6));

  const traces = [
    {
      type: "bar",
      orientation: "h",
      y: labels,
      x: kept,
      name: "retained",
      marker: { color: "rgba(5, 150, 105, 0.85)" },
      hovertemplate: "%{y}<br>Retained: %{x:,}<extra></extra>",
    },
    {
      type: "bar",
      orientation: "h",
      y: labels,
      x: filtered,
      name: "filtered out",
      marker: { color: "rgba(59, 130, 246, 0.35)" },
      hovertemplate: "%{y}<br>Filtered out: %{x:,}<extra></extra>",
    },
  ];

  clearPlotEmptyState(plot);
  Plotly.react(
    plot,
    traces,
    {
      barmode: "stack",
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.88)",
      margin: { t: 10, l: marginLeft, r: 16, b: 48 },
      height: 320,
      xaxis: {
        title: { text: "Cells (relative to total detected)" },
        range: [0, total],
        tickformat: ",d",
        showgrid: true,
        gridcolor: "rgba(148,163,184,0.2)",
      },
      yaxis: {
        autorange: "reversed",
      },
      legend: {
        orientation: "h",
        x: 0,
        y: 1.15,
      },
    },
    { responsive: true, displayModeBar: false }
  );

  caption.textContent = `Total detected cells: ${total.toLocaleString()}.`;
  renderAmbientCorrectionProgress(ambientPlot, ambientCaption, ambientRows, data);
}

function renderAmbientCorrectionProgress(plot, caption, ambientRows, data) {
  if (!plot || !caption) {
    return;
  }
  const ambientSetting = String((((data || {}).qc || {}).ambient_correction || (((currentJobStatus || {}).qc || {}).ambient_correction || "no"))).trim().toLowerCase();
  const ambientEnabled = ambientSetting === "yes";
  plot.style.display = ambientEnabled ? "" : "none";
  caption.style.display = ambientEnabled ? "" : "none";
  if (!ambientEnabled) {
    try {
      Plotly.purge(plot);
    } catch (_) {
      // Ignore cleanup errors when no Plotly instance exists yet.
    }
    plot.innerHTML = "";
    caption.textContent = "";
    return;
  }
  if (!ambientRows.length) {
    try {
      Plotly.purge(plot);
    } catch (_) {
      // Ignore cleanup errors when no Plotly instance exists yet.
    }
    plot.innerHTML = "<div class=\"empty-state\">Ambient RNA correction percentages will appear here when correction is performed.</div>";
    caption.textContent = "No ambient correction percentages reported yet.";
    return;
  }

  const labels = ambientRows.map((row) => row.library);
  const correctedPct = ambientRows.map((row) => Number((row.rho * 100).toFixed(2)));
  const remainingPct = correctedPct.map((value) => Math.max(0, 100 - value));
  const maxLabelLength = labels.reduce((max, value) => Math.max(max, String(value).length), 0);
  const marginLeft = Math.min(210, Math.max(120, 18 + maxLabelLength * 6));

  clearPlotEmptyState(plot);
  Plotly.react(
    plot,
    [
      {
        type: "bar",
        orientation: "h",
        y: labels,
        x: correctedPct,
        name: "ambient correction",
        marker: { color: "rgba(239, 68, 68, 0.8)" },
        hovertemplate: "%{y}<br>Ambient correction: %{x:.1f}%<extra></extra>",
      },
      {
        type: "bar",
        orientation: "h",
        y: labels,
        x: remainingPct,
        name: "uncorrected",
        marker: { color: "rgba(148, 163, 184, 0.22)" },
        hovertemplate: "%{y}<br>Uncorrected remainder: %{x:.1f}%<extra></extra>",
      },
    ],
    {
      barmode: "stack",
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.88)",
      margin: { t: 2, l: marginLeft, r: 16, b: 42 },
      height: Math.max(170, 40 + labels.length * 34),
      xaxis: {
        title: { text: "Ambient RNA correction (%)" },
        range: [0, 100],
        ticksuffix: "%",
        showgrid: true,
        gridcolor: "rgba(148,163,184,0.2)",
      },
      yaxis: {
        autorange: "reversed",
      },
      legend: {
        orientation: "h",
        x: 0.5,
        xanchor: "center",
        y: 1.12,
        yanchor: "bottom",
        traceorder: "normal",
      },
    },
    { responsive: true, displayModeBar: false }
  );
  caption.textContent = `Per-sample ambient RNA correction selected automatically from the alignment log: ${ambientRows.map((row) => `${row.library} ${Math.round(row.rho * 100)}%`).join(" | ")}.`;
}

function buildQcCellSummary(data) {
  const lines = getStatusLogLines(data);
  const qcState = extractQcThresholdState(lines);
  const status = String(data.status || "").trim().toLowerCase();
  const message = String(data.message || "").trim();
  const progress = parseProgressPercent(data.progress);
  let alignmentExcluded = null;
  for (const line of lines) {
    const match = line.match(/Applied min_alignment_score=.*?Excluded\s+(\d+)\s+cells,\s+kept\s+(\d+)/i);
    if (match) {
      alignmentExcluded = Number(match[1]);
    }
  }

  const stageMarkers = [
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
  ];
  if (status !== "completed" && status !== "failed") {
    for (let index = lines.length - 1; index >= 0; index -= 1) {
      const line = String(lines[index] || "");
      const matchedStage = stageMarkers.find((marker) => line.includes(marker));
      if (matchedStage) {
        return matchedStage;
      }
      if (line.includes("Aligning cells to reference")) {
        return "Aligning cells to reference...";
      }
      if (line.includes("Normalization steps")) {
        return "Normalizing expression values...";
      }
    }
    if (message) {
      const genericMessages = new Set([
        "Preparing inputs…",
        "Preparing inputs...",
        "Job submitted to worker.",
      ]);
      if (!genericMessages.has(message)) {
        return message;
      }
    }
    if (progress >= 82) {
      return "Running approximate UMAP placement.";
    }
    if (progress >= 78) {
      return "Exporting NetPerspective marker networks.";
    }
    if (progress >= 72) {
      return "Identifying the top 50 unique markers for 100 cells per cell state.";
    }
  }

  if (status === "failed" && message) {
    return `Analysis failed. ${message}`;
  }

  if (status === "completed") {
    if (!referenceRerunPending && !areExploreResultsReady()) {
      return "Alignment completed. Finalizing Explore results and loading the combined h5ad for interactive viewing...";
    }
    const createdAt = Date.parse(data.created_at || "");
    const updatedAt = Date.parse(data.updated_at || "");
    const exclusionSuffix =
      alignmentExcluded !== null
        ? ` ${alignmentExcluded.toLocaleString()} cells excluded due to poor alignment.`
        : "";
    if (Number.isFinite(createdAt) && Number.isFinite(updatedAt) && updatedAt >= createdAt) {
      const elapsedSeconds = Math.max(0, Math.round((updatedAt - createdAt) / 1000));
      return `Analysis completed and results saved in ${elapsedSeconds} seconds.${exclusionSuffix}`;
    }
    return `Analysis completed and results saved.${exclusionSuffix}`.trim();
  }

  const segments = [];
  if (Number.isFinite(qcState.total)) {
    segments.push(`Cells before QC: ${qcState.total.toLocaleString()}`);
  }
  if (Number.isFinite(qcState.afterMinGenes)) {
    segments.push(`After min genes: ${qcState.afterMinGenes.toLocaleString()}`);
  }
  if (Number.isFinite(qcState.afterMinCounts)) {
    segments.push(`After min counts: ${qcState.afterMinCounts.toLocaleString()}`);
  }
  if (Number.isFinite(qcState.afterMito)) {
    segments.push(`After mito filter: ${qcState.afterMito.toLocaleString()}`);
  }

  if (segments.length) {
    return segments.join(" | ");
  }
  if (data.status === "queued" || data.status === "processing") {
    return "QC cell counts will appear here while the analysis runs.";
  }
  return "QC cell counts will appear here while the analysis runs.";
}

async function populateDownloadLinks(jobId, statusData = null) {
  const container = document.getElementById("download-links");
  container.innerHTML = "";
  let data = statusData;
  if (!data) {
    const resp = await fetch(apiPath(`/jobs/${jobId}/status`));
    data = await resp.json();
    if (!resp.ok) {
      return;
    }
  }
  const artifacts = data.artifacts || {};
  const labelMap = {
    assignments: "Download assignments",
    combined_h5ad: "Download combined_h5ad",
    marker_genes_zip: "Download marker genes ZIP",
    imputed_lipids_h5ad: "Download imputed_lipids_h5ad",
    lipid_marker_genes_zip: "Download lipid marker ZIP",
    imputed_lipids_results_zip: "Download Lipids results ZIP",
    imputed_adt_results_zip: "Download ADT results ZIP",
    imputed_metabolite_results_zip: "Download Metabolite results ZIP",
    imputed_lipid_results_zip: "Download Lipid (AML) results ZIP",
    imputed_grn_results_zip: "Download GRN edge results ZIP",
    imputed_grn_tf_results_zip: "Download TF activity results ZIP",
    fastcomm_archive: "Download cell communication ZIP",
  };
  Object.keys(artifacts).forEach((key) => {
    if (
      key === "umap_coordinates" ||
      key === "umap_placeholder_expression" ||
      key === "umap_pdf" ||
      key === "umap_pdf_plain" ||
      key === "imputed_lipids_summary_json" ||
      (key.startsWith("fastcomm_") && key !== "fastcomm_archive")
    ) {
      return;
    }
    const btn = document.createElement("a");
    btn.className = "download-btn";
    btn.href = apiPath(`/jobs/${jobId}/download/${key}`);
    btn.textContent = labelMap[key] || `Download ${key}`;
    container.appendChild(btn);
  });
  if (jobId && data && data.status === "completed") {
    const logBtn = document.createElement("a");
    logBtn.className = "download-btn";
    logBtn.href = apiPath(`/jobs/${jobId}/log`);
    logBtn.textContent = "Download log";
    container.appendChild(logBtn);
  }
}

function updateDifferentialUi(state) {
  if (currentDifferentialState?.run_id !== state?.run_id ||
      currentDifferentialState?.config?.modality !== state?.config?.modality) {
    ++differentialVisualizationRequest;
    ++differentialDetailRequest;
    currentDifferentialGene = "";
    clearDifferentialGeneFilter();
  }
  currentDifferentialState = state;
  const panel = document.getElementById("differential-panel");
  const emptyState = document.getElementById("differential-tab-empty");
  const resultsView = document.getElementById("differential-results-view");
  const intro = document.getElementById("differential-intro");
  const populationSelect = document.getElementById("differential-population");
  const modalityField = document.getElementById("differential-modality-field");
  const modalitySelect = document.getElementById("differential-modality");
  const sampleFieldSelect = document.getElementById("differential-sample-field");
  const group1Select = document.getElementById("differential-group1");
  const group2Select = document.getElementById("differential-group2");
  const comparisonField = document.getElementById("differential-comparison-type-field");
  const comparisonSelect = document.getElementById("differential-comparison-type");
  const runBtn = document.getElementById("differential-run-btn");
  const progress = state ? state.progress || 0 : 0;
  const message = document.getElementById("differential-message");
  const archiveLink = document.getElementById("differential-archive-link");
  const logLink = document.getElementById("differential-log-link");
  const vizModeSelect = document.getElementById("differential-viz-mode");
  const detailTitle = document.getElementById("differential-detail-title");

  const enabled = Boolean(state && state.enabled);
  const config = (state && state.config) || {};
  const populationOptions = (state && state.population_columns) || [];
  let completedField = document.getElementById("differential-completed-field");
  if (!completedField) {
    completedField = document.createElement("label");
    completedField.id = "differential-completed-field";
    completedField.className = "field";
    completedField.innerHTML = '<span>Completed comparison</span><select id="differential-completed-run"></select>';
    panel.prepend(completedField);
    document.getElementById("differential-completed-run").addEventListener("change", async (event) => {
      if (!event.target.value) return;
      event.target.disabled = true;
      try { await selectCompletedDifferential(event.target.value); }
      finally { event.target.disabled = false; }
    });
  }
  const completed = (state && state.completed_comparisons) || [];
  completedField.classList.toggle("hidden", !completed.length || Boolean(window.__SCALABLE_VIEWER__));
  populateSingleSelect(document.getElementById("differential-completed-run"),
    completed.map(c => ({value: c.id, label: `${modalityDefinition(c.modality).label}: ${c.comparison}`})), state && state.run_id);
  document.getElementById("differential-completed-run").disabled = ["queued", "processing"].includes(state && state.status);
  const modalityOptions = (state && state.modalities) || [{ id: "rna", label: "RNA", feature_label: "gene" }];
  const sampleFieldOptions = (state && state.sample_fields) || [];
  const sampleValuesMap = (state && state.sample_values) || {};
  const optionHasValue = (options, value) => {
    const target = String(value || "").trim();
    return Boolean(target) && (options || []).some((optionData) => String(optionData.value || "").trim() === target);
  };
  const differentialRunning = state && (state.status === "queued" || state.status === "processing");
  const currentPopulationValue = populationSelect.value;
  const currentModalityValue = normalizeModalityId(modalitySelect.value || "rna");
  const currentSampleFieldValue = sampleFieldSelect.value;
  const selectedModality = (
    (differentialRunning && modalityOptions.some((entry) => normalizeModalityId(entry.id) === currentModalityValue) && currentModalityValue)
      || normalizeModalityId(config.modality || (state && state.default_modality) || "rna")
  );
  const featureLabel = String(((modalityOptions.find((entry) => normalizeModalityId(entry.id) === selectedModality) || {}).feature_label) || "gene");
  const selectedPopulation = (
    (differentialRunning && optionHasValue(populationOptions, currentPopulationValue) && currentPopulationValue)
      || config.population_col
      || (state && state.default_population_col)
      || ""
  );
  const selectedSampleField = (
    (differentialRunning && optionHasValue(sampleFieldOptions, currentSampleFieldValue) && currentSampleFieldValue)
      || config.sample_field
      || (state && state.default_sample_field)
      || ""
  );
  populateSingleSelect(modalitySelect, modalityOptions.map((entry) => ({ value: entry.id, label: entry.label })), selectedModality);
  modalityField.classList.toggle("hidden", modalityOptions.length <= 1);
  populateSingleSelect(sampleFieldSelect, sampleFieldOptions, selectedSampleField);
  const sampleValues = sampleValuesMap[sampleFieldSelect.value] || [];
  const comparisonTypes = (state && state.comparison_types) || ["cells"];
  const showComparisonType = comparisonTypes.includes("pseudobulk");
  const selectedComparisonType = showComparisonType && comparisonTypes.includes(config.comparison_type)
    ? (config.comparison_type || "cells")
    : "cells";
  const showPanel = Boolean(state && state.enabled && populationOptions.length && sampleFieldOptions.length);

  panel.classList.toggle("hidden", !showPanel);

  populateSingleSelect(populationSelect, populationOptions, selectedPopulation);
  populateMultiSelect(group1Select, sampleValues, config.group1_samples || []);
  populateMultiSelect(group2Select, sampleValues, config.group2_samples || []);
  comparisonField.classList.toggle("hidden", !showComparisonType);
  Array.from(comparisonSelect.options).forEach((option) => {
    option.hidden = !comparisonTypes.includes(option.value);
    option.disabled = !comparisonTypes.includes(option.value);
  });
  comparisonSelect.value = selectedComparisonType;

  const disableInputs = !enabled || !populationOptions.length || !sampleFieldOptions.length;
  populationSelect.disabled = disableInputs;
  modalitySelect.disabled = disableInputs || modalityOptions.length <= 1;
  sampleFieldSelect.disabled = disableInputs;
  group1Select.disabled = disableInputs;
  group2Select.disabled = disableInputs;
  comparisonSelect.disabled = disableInputs || !showComparisonType;
  runBtn.disabled = disableInputs;

  document.getElementById("differential-progress").style.width = `${progress}%`;
  document.getElementById("differential-progress-label").textContent = `${progress}%`;

  if (!state) {
    modalityField.classList.add("hidden");
    if (emptyState) {
      emptyState.textContent = "Available after alignment completes for jobs with two or more samples.";
      emptyState.classList.remove("hidden");
    }
    if (resultsView) {
      resultsView.classList.add("hidden");
    }
    intro.textContent = "Available after alignment completes for jobs with two or more samples.";
    message.textContent = "Differential analysis is enabled when the job contains two or more samples.";
    archiveLink.classList.add("hidden");
    logLink.classList.add("hidden");
    resetDifferentialResults();
    setResultMode("baseline");
    return;
  }

  if (!showPanel) {
    modalityField.classList.add("hidden");
    if (emptyState) {
      emptyState.textContent = "Differential gene analyses between biological groups (i.e., disease versus controls) are only enabled when two or more samples (multiple h5 files or a single h5ad) are uploaded for the job.";
      emptyState.classList.remove("hidden");
    }
    if (resultsView) {
      resultsView.classList.add("hidden");
    }
    archiveLink.classList.add("hidden");
    logLink.classList.add("hidden");
    resetDifferentialResults();
    setResultMode("baseline");
    return;
  }

  if (emptyState) {
    emptyState.classList.add("hidden");
  }

  intro.textContent = enabled
    ? `Perform cell-state differential ${featureLabel} analysis.`
    : "Differential analysis is only enabled when two or more samples were uploaded for the job.";
  if (detailTitle) {
    detailTitle.textContent = selectedModality === "grn"
      ? "GRN Detail"
      : (featureLabel === "lipid" ? "Lipid Detail" : "Gene Detail");
  }

  let statusMessage = state.message || "";
  if (state.status === "completed" && !state.go_terms_included && (state.visualization_modes || []).some((entry) => entry.value === "go")) {
    statusMessage += " GO terms were not available for this run.";
  }
  message.textContent = statusMessage;

  const currentVizMode = vizModeSelect.value;
  const visualizationModes = [...(state.visualization_modes || [])];
  if (crossPathwayContexts.differential?.jobId === getResultsJobId() && crossPathwayContexts.differential?.contrast === state.run_id)
    visualizationModes.push({value:"integrated_cross_pathway",label:"Pathway (cross-modality)"});
  vizModeSelect.innerHTML = "";
  visualizationModes.forEach((entry) => {
    const option = document.createElement("option");
    option.value = entry.value;
    option.textContent = entry.label;
    if (entry.value === currentVizMode) {
      option.selected = true;
    }
    vizModeSelect.appendChild(option);
  });
  if (!vizModeSelect.value && visualizationModes.length) {
    vizModeSelect.value = visualizationModes[0].value;
  }

  if (state.archive_url) {
    archiveLink.href = state.archive_url;
    archiveLink.classList.remove("hidden");
  } else {
    archiveLink.classList.add("hidden");
  }
  if (loadedResultsJobId && state && state.status === "completed") {
    logLink.href = apiPath(`/jobs/${loadedResultsJobId}/log`);
    logLink.classList.remove("hidden");
  } else {
    logLink.classList.add("hidden");
  }

  if (state.status === "completed") {
    setResultMode("differential");
    if ((state.result_populations || []).length) {
      renderDifferentialResults(state);
    } else {
      resetDifferentialResults();
      renderDifferentialEmpty(
        `Differential analysis completed, but no significant DE ${featureLabel}s were available for the selected comparison.`
      );
    }
  } else {
    resetDifferentialResults();
    setResultMode("baseline");
  }
}

function populateSingleSelect(element, options, selectedValue) {
  const currentValue = selectedValue || element.value;
  element.innerHTML = "";
  options.forEach((optionData) => {
    const option = document.createElement("option");
    option.value = optionData.value;
    option.textContent = optionData.n_categories
      ? `${optionData.label} (${optionData.n_categories})`
      : optionData.label;
    if (option.value === currentValue) {
      option.selected = true;
    }
    element.appendChild(option);
  });
}

function populateMultiSelect(element, values, selectedValues) {
  const wanted = new Set(selectedValues || []);
  const existingSelections = wanted.size ? wanted : new Set(getMultiSelectValues(element));
  element.innerHTML = "";
  values.forEach((value) => {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = value;
    option.selected = existingSelections.has(value);
    element.appendChild(option);
  });
}

function initToggleMultiSelect(element) {
  if (!element || element.dataset.toggleMultiSelectBound === "true") {
    return;
  }
  element.dataset.toggleMultiSelectBound = "true";
  element.addEventListener("mousedown", (event) => {
    const target = event.target;
    if (!(target instanceof HTMLOptionElement)) {
      return;
    }
    event.preventDefault();
    if (element.disabled) {
      return;
    }
    target.selected = !target.selected;
    element.focus();
    element.dispatchEvent(new Event("change", { bubbles: true }));
  });
}

function getMultiSelectValues(element) {
  return Array.from(element.selectedOptions).map((option) => option.value);
}

function renderDifferentialResults(state) {
  syncDifferentialPopulationSelect(state);
  updateDifferentialDownloadButton();
  loadDifferentialVisualization();
}

function differentialPopulationsForMode(state, mode) {
  const mapping = (state && state.visualization_populations) || {};
  const values = mapping[mode] || [];
  return values.filter((value, index) => value && values.indexOf(value) === index);
}

function syncDifferentialPopulationSelect(state) {
  syncDifferentialGeneFilterLabel();
  const populationSelect = document.getElementById("differential-result-population");
  const mode = document.getElementById("differential-viz-mode").value;
  const populations = differentialPopulationsForMode(state, mode);
  // The counts chart draws every cell state, so neither the state select nor the
  // one-gene filter applies to it.
  populationSelect.classList.toggle("hidden", mode === "summary");
  const filterInput = differentialGeneFilterInput();
  if (filterInput) {
    filterInput.classList.toggle("hidden", mode === "summary" || mode.startsWith("integrated_"));
  }
  const fallbackPopulation = state.default_result_population || "";
  const wantedPopulation = currentDifferentialPopulation || populationSelect.value || fallbackPopulation;
  populationSelect.innerHTML = "";
  // The volcano covers every cell state the differential tested, read from the sample
  // pseudobulks. The detail violin draws this atlas's released replicate unit, which
  // covers fewer states: the COPD atlas releases 50 metacell states of the 81 tested.
  // Marking the difference stops a reader picking a state whose distribution cannot be
  // drawn, which is why "Deuterosomal" looked broken.
  populations.forEach((population) => {
    const option = document.createElement("option");
    option.value = population;
    option.textContent = (replicateStates && !replicateStates.has(population))
      ? `${population} (no replicates)`
      : population;
    if (population === wantedPopulation) {
      option.selected = true;
    }
    populationSelect.appendChild(option);
  });
  populationSelect.disabled = populations.length === 0;
  if (!populationSelect.value && populations.length) {
    populationSelect.value = populations[0];
  }
  if (!populations.includes(populationSelect.value)) {
    populationSelect.value = populations[0] || "";
  }
  currentDifferentialPopulation = populationSelect.value || "";
}

// ---------------------------------------------------------------------------------
// Differential gene filter
//
// The box beside the cell-state select takes ONE gene. That gene resolves to a gene
// SET: the gene itself plus every gene directly connected to it in the interaction
// network of the selected cell state, the network `/differential/interactive/network`
// already serves. Every view is then restricted to that set - Heatmap rows, Volcano
// points, GO terms that overlap it, and the Network subgraph. One set, four views.
//
// The set is built from the network payload, never from a new statistic. When the gene
// has no edge in that cell state, or the cell state has no network, the set is the gene
// alone and the note under the controls says so.
// ---------------------------------------------------------------------------------

function differentialGeneFilterInput() {
  return document.getElementById("differential-gene-filter");
}

// scALABLE runs the Differential tab on five modalities. RNA names a gene, ADT names a
// protein, GRN names a TF. The filter box and its note follow the modality's own noun.
function differentialFeatureNoun() {
  return String((currentDifferentialState && currentDifferentialState.feature_label) || "gene");
}

function syncDifferentialGeneFilterLabel() {
  const input = differentialGeneFilterInput();
  if (!input) {
    return;
  }
  const noun = differentialFeatureNoun();
  input.placeholder = `Filter by ${noun}`;
  input.setAttribute(
    "aria-label",
    `Filter the differential view to one ${noun} and the ${noun}s it interacts with`,
  );
}

function setDifferentialGeneFilterNote(message) {
  const note = document.getElementById("differential-gene-filter-note");
  if (!note) {
    return;
  }
  const text = String(message || "").trim();
  note.textContent = text;
  note.classList.toggle("hidden", !text);
}

function setDifferentialGeneFilterOptions(genes) {
  const datalist = document.getElementById("differential-gene-filter-options");
  if (!datalist) {
    return;
  }
  const unique = Array.from(new Set((genes || []).map((gene) => String(gene || "").trim()).filter(Boolean)));
  unique.sort((left, right) => left.localeCompare(right));
  datalist.innerHTML = "";
  unique.forEach((gene) => {
    const option = document.createElement("option");
    option.value = gene;
    datalist.appendChild(option);
  });
}

// The autocomplete universe is the genes of the view now loaded, so any gene offered
// always yields a non-empty result.
function differentialPayloadGenes(mode, payload) {
  if (!payload) {
    return [];
  }
  if (mode === "heatmap") {
    return (payload.rows || []).map((row) => row.gene);
  }
  if (mode === "volcano") {
    return (payload.points || []).map((point) => point.gene);
  }
  if (mode === "network") {
    return (payload.elements || [])
      .filter((element) => element && element.data && element.data.id && !element.data.source)
      .map((element) => element.data.id);
  }
  if (mode === "go") {
    const genes = [];
    (payload.terms || []).forEach((term) => {
      (term.overlap_genes || []).forEach((gene) => genes.push(gene));
    });
    return genes;
  }
  if (mode === "table") {
    // Cell communication: the features are the ligands and the receptors.
    const genes = [];
    (payload.rows || []).forEach((row) => {
      if (row.ligand) genes.push(row.ligand);
      if (row.receptor) genes.push(row.receptor);
    });
    return genes;
  }
  return [];
}

// The ligands and receptors one cell-communication edge carries. The payload states them
// on the edge itself and repeats them inside `top_interactions`, so both are read.
function communicationEdgeFeatures(data) {
  const features = [];
  if (data.ligand) features.push(String(data.ligand));
  if (data.receptor) features.push(String(data.receptor));
  (data.top_interactions || []).forEach((entry) => {
    if (entry && entry.ligand) features.push(String(entry.ligand));
    if (entry && entry.receptor) features.push(String(entry.receptor));
  });
  return features;
}

function communicationRowMatches(row, geneFilter) {
  if (!geneFilter) {
    return true;
  }
  return geneFilter.genes.has(String(row.ligand || ""))
    || geneFilter.genes.has(String(row.receptor || ""));
}

// A cell-communication view has no gene-gene network behind it, so the filter never
// expands to interacting genes there. The note must say what it did match.
function describeCommunicationGeneFilter(geneFilter, population, kept, total, noun) {
  if (!geneFilter) {
    return "";
  }
  return `Filtered to interactions using ${geneFilter.gene} as ligand or receptor: `
    + `${kept} of ${total} ${noun} shown.`;
}

async function differentialNetworkAdjacency(jobId, population) {
  const key = `${jobId}::${population}`;
  if (differentialNetworkAdjacencyCache[key]) {
    return differentialNetworkAdjacencyCache[key];
  }
  let payload = null;
  try {
    payload = await fetchDifferentialJson(
      apiPath(`/jobs/${jobId}/differential/interactive/network?population=${encodeURIComponent(population)}`),
    );
  } catch (error) {
    payload = null;
  }
  const adjacency = new Map();
  const geneNetwork = Boolean(payload) && payload.network_type !== "cell_communication_diff";
  if (geneNetwork) {
    (payload.elements || []).forEach((element) => {
      const data = element && element.data ? element.data : null;
      if (!data) {
        return;
      }
      if (data.source && data.target) {
        const source = String(data.source);
        const target = String(data.target);
        if (!adjacency.has(source)) adjacency.set(source, new Set());
        if (!adjacency.has(target)) adjacency.set(target, new Set());
        adjacency.get(source).add(target);
        adjacency.get(target).add(source);
      } else if (data.id) {
        const node = String(data.id);
        if (!adjacency.has(node)) adjacency.set(node, new Set());
      }
    });
  }
  const entry = { adjacency, geneNetwork, available: Boolean(payload) && adjacency.size > 0 };
  differentialNetworkAdjacencyCache[key] = entry;
  return entry;
}

// Returns null when no gene is typed, so every render path stays unfiltered by default.
async function resolveDifferentialGeneFilter(jobId, population) {
  const gene = String(currentDifferentialGeneFilter || "").trim();
  if (!gene) {
    return null;
  }
  const entry = await differentialNetworkAdjacency(jobId, population);
  const genes = new Set([gene]);
  let neighbours = 0;
  if (entry.geneNetwork && entry.adjacency.has(gene)) {
    entry.adjacency.get(gene).forEach((partner) => {
      genes.add(partner);
      neighbours += 1;
    });
  }
  return {
    gene,
    genes,
    neighbours,
    inNetwork: entry.geneNetwork && entry.adjacency.has(gene),
    geneNetwork: entry.geneNetwork,
    networkAvailable: entry.available,
  };
}

function describeDifferentialGeneFilter(filter, population, kept, total, noun) {
  if (!filter) {
    return "";
  }
  const feature = differentialFeatureNoun();
  const scope = filter.inNetwork
    ? `${filter.gene} + ${filter.neighbours} interacting ${feature}${filter.neighbours === 1 ? "" : "s"}`
    : (filter.networkAvailable
      ? `${filter.gene} alone (no interaction edge in ${population})`
      : `${filter.gene} alone (no interaction network for ${population})`);
  return `Filtered to ${scope}: ${kept} of ${total} ${noun} shown.`;
}

function clearDifferentialGeneFilter() {
  currentDifferentialGeneFilter = "";
  differentialNetworkAdjacencyCache = {};
  const input = differentialGeneFilterInput();
  if (input) {
    input.value = "";
    input.classList.remove("filter-active", "filter-unmatched");
  }
  setDifferentialGeneFilterNote("");
  setDifferentialGeneFilterOptions([]);
}

function markDifferentialGeneFilterMatched(matched) {
  const input = differentialGeneFilterInput();
  if (!input) {
    return;
  }
  const active = Boolean(String(currentDifferentialGeneFilter || "").trim());
  input.classList.toggle("filter-active", active && matched);
  input.classList.toggle("filter-unmatched", active && !matched);
}


const crossPathwayContexts = {};
async function openCrossPathway(row, spec) {
  const jobId = getResultsJobId();
  const mode = "integrated_cross_pathway";
  const context = {...spec, id:row.id, jobId, crossModal:true, kind:"integrated_pathway"};
  const setOption = (id, value, label) => {
    const select = document.getElementById(id);
    if (![...select.options].some(o=>o.value === value)) select.add(new Option(label || value,value));
    select.value=value;
  };
  if (spec.source === "differential") {
    await selectCompletedDifferential(spec.contrast);
    if (currentDifferentialState?.run_id !== spec.contrast) return;
    crossPathwayContexts.differential=context;
    setOption("differential-viz-mode", mode, "Pathway (cross-modality)");
    setOption("differential-result-population", spec.cell_state);
    document.querySelector('[data-tab="differential"]')?.click();
    await loadDifferentialVisualization();
  } else {
    crossPathwayContexts.marker=context;
    updateExpressionModeOptions();
    setOption("viz1-mode",mode,"Pathway (cross-modality)");
    updateExpressionModeOptions();
    setOption("viz1-marker-population",spec.cell_state);
    document.querySelector('[data-tab="explore"]')?.click();
    await loadVisualizationPanel("viz1");
  }
}

function integratedOptions(jobId, kind, cellState, modality, source="differential") {
  return {jobId,kind,cell_state:cellState,modality,source,externalPdf:true,availableModalities:availableModalities().map(m=>m.id),
    contrast: window.__SCALABLE_VIEWER__?.contrast || currentDifferentialState?.run_id || "",apiPath,
    ...(kind === "integrated_cross_pathway" ? {...crossPathwayContexts[source], cell_state:cellState} : {}),
    onFeature: (gene, kind) => {
      const panel="viz1",mod=document.getElementById("viz1-modality");
      if(mod && [...mod.options].some(o=>o.value===kind))mod.value=kind;
      updateExpressionModeOptions();document.getElementById("viz1-mode").value="expression_umap";
      setPanelGeneValue(panel,gene);
      document.querySelector('[data-tab="explore"]')?.click();loadVisualizationPanel(panel);
    }};
}

async function loadDifferentialVisualization() {
  const request = ++differentialVisualizationRequest;
  ++differentialDetailRequest;
  const state = currentDifferentialState;
  const plotEmpty = document.getElementById("differential-plot-empty");
  const jobId = document.getElementById("results-job-id").value.trim();
  await ensureFeatureAnnotations(jobId);
  if (!state || state.status !== "completed") {
    resetDifferentialResults();
    return;
  }

  const population = document.getElementById("differential-result-population").value;
  const mode = document.getElementById("differential-viz-mode").value;
  const isCurrent = () => request === differentialVisualizationRequest &&
    state.run_id === currentDifferentialState?.run_id &&
    jobId === document.getElementById("results-job-id").value.trim();
  const integratedMode=mode.startsWith("integrated_");
  document.getElementById("differential-results-view").classList.toggle("show-integrated",integratedMode);
  if(!integratedMode)document.getElementById("differential-plot-area")._integratedDispose?.();
  updateDifferentialDownloadButton();

  // The counts chart draws every cell state at once, so it needs no selected state.
  if (mode === "summary") {
    updateDifferentialDownloadButton();
    try {
      const summary = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/summary`));
      if (!isCurrent()) return;
      renderDifferentialSummary(summary);
      setDifferentialGeneFilterOptions([]);
      resetDifferentialGeneDetail();
    } catch (err) {
      if (!isCurrent()) return;
      destroyDifferentialNetwork();
      renderDifferentialEmpty(err.message || "Unable to load the differential counts.");
      plotEmpty.classList.remove("hidden");
    }
    return;
  }
  if (!population) {
    renderDifferentialEmpty(`No ${mode} data are available for this differential run.`);
    resetDifferentialGeneDetail();
    return;
  }

  currentDifferentialPopulation = population;
  if(mode.startsWith("integrated_")) {
    destroyDifferentialNetwork();resetDifferentialGeneDetail();
    const host=document.getElementById("differential-plot-area");
    try {Plotly.purge(host);} catch(e) {}
    host.classList.remove("hidden");plotEmpty.classList.add("hidden");
    await ScalableIntegrated.mount(host,integratedOptions(jobId,mode,population,state.config?.modality || state.modality || "rna"));
    return;
  }
  updateDifferentialDownloadButton();
  try {
    let payload = null;
    const geneFilter = await resolveDifferentialGeneFilter(jobId, population);
    if (!isCurrent()) return;
    if (mode === "heatmap") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/heatmap?population=${encodeURIComponent(population)}`));
      if (!isCurrent()) return;
      renderDifferentialHeatmap(payload, geneFilter);
    } else if (mode === "volcano") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/volcano?population=${encodeURIComponent(population)}`));
      if (!isCurrent()) return;
      renderDifferentialVolcano(payload, geneFilter);
    } else if (mode === "network") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/network?population=${encodeURIComponent(population)}`));
      if (!isCurrent()) return;
      renderDifferentialNetwork(payload, geneFilter);
    } else if (mode === "table") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/table?population=${encodeURIComponent(population)}`));
      if (!isCurrent()) return;
      renderDifferentialCommunicationTable(payload, geneFilter);
    } else if (mode === "go") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/go?population=${encodeURIComponent(population)}`));
      if (!isCurrent()) return;
      renderDifferentialGo(payload, geneFilter);
    }
    const availableGenes = differentialPayloadGenes(mode, payload);
    setDifferentialGeneFilterOptions(availableGenes);
    const payloadDefaultGene = (payload && payload.default_gene) || "";
    // A feature carried over from another modality does not exist in this one, so the
    // detail fetch failed and the panel fell back to "Select a <feature> from the
    // differential view". Switching RNA -> ADT kept an RNA gene and emptied the panel.
    // Treat a carried-over feature as stale so the payload's own default is used, which
    // the server sets to the top result by p-value.
    if (currentDifferentialGene && Array.isArray(availableGenes) && availableGenes.length
        && !availableGenes.some((g) => String(g) === String(currentDifferentialGene))) {
      currentDifferentialGene = "";
    }
    const isCellCommunicationMode = mode === "network" || mode === "table";
    let staleGene = "";
    if (isCellCommunicationMode && currentDifferentialGene) {
      const receiverPart = currentDifferentialGene.split(":")[0].split("->").slice(-1)[0] || "";
      if (receiverPart && receiverPart !== population) {
        staleGene = currentDifferentialGene;
      }
    }
    const nextGene = (staleGene ? "" : currentDifferentialGene) || payloadDefaultGene || "";
    if (nextGene) {
      currentDifferentialGene = nextGene;
      let defaultFeature = "";
      if (isCellCommunicationMode) {
        if (mode === "table" && Array.isArray(payload && payload.rows) && payload.rows.length) {
          const targetRow = payload.rows.find((row) => String(row.gene || "") === nextGene) || payload.rows[0];
          const ligand = targetRow && targetRow.ligand ? String(targetRow.ligand) : "";
          const receptor = targetRow && targetRow.receptor ? String(targetRow.receptor) : "";
          const role = currentDifferentialFeatureRole === "receptor" && receptor ? "receptor" : "ligand";
          defaultFeature = role === "ligand" ? (ligand || receptor) : (receptor || ligand);
          currentDifferentialFeatureRole = role;
        } else if (mode === "network" && Array.isArray(payload && payload.elements)) {
          for (const element of payload.elements) {
            const data = element && element.data ? element.data : null;
            if (!data || !data.source) continue;
            const tops = data.top_interactions || [];
            const match = tops.find((entry) => String(entry.gene || "") === nextGene);
            if (match) {
              const ligand = String(match.ligand || "");
              const receptor = String(match.receptor || "");
              const role = currentDifferentialFeatureRole === "receptor" && receptor ? "receptor" : "ligand";
              defaultFeature = role === "ligand" ? (ligand || receptor) : (receptor || ligand);
              currentDifferentialFeatureRole = role;
              break;
            }
          }
        }
      }
      await loadDifferentialGeneDetail(nextGene, population, defaultFeature ? { feature: defaultFeature } : undefined);
    } else {
      resetDifferentialGeneDetail();
    }
  } catch (err) {
    if (!isCurrent()) return;
    destroyDifferentialNetwork();
    renderDifferentialEmpty(err.message || "Unable to load the differential visualization.");
    plotEmpty.classList.remove("hidden");
  }
}

async function fetchDifferentialJson(url) {
  const resp = await fetch(url);
  const data = await parseApiResponse(resp);
  if (!resp.ok) {
    throw new Error(data.detail || "Differential data request failed.");
  }
  return data;
}

// The two bar colours, read from the reference figure Nathan supplied.
const DEG_COUNT_UP_COLOR = "#C75252";
const DEG_COUNT_DOWN_COLOR = "#7CC7E9";

// The first view of the Differential workspace: how many features the run called up
// and called down in each cell state. Down and up share one row, down to the left of
// zero and up to the right. Clicking either bar opens that cell state in the Volcano.
function renderDifferentialSummary(payload) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  const rows = payload.rows || [];
  setDifferentialGeneFilterNote("");
  markDifferentialGeneFilterMatched(true);
  if (!rows.length) {
    renderDifferentialEmpty("No differential counts were found for this comparison.");
    return;
  }
  // Plotly draws the first category at the bottom of a horizontal bar axis, so the
  // array is reversed to put the first cell state of the lineage order on top.
  const ordered = rows.slice().reverse();
  const labels = ordered.map((row) => String(row.population));
  const upValues = ordered.map((row) => Number(row.up) || 0);
  const downValues = ordered.map((row) => -(Number(row.down) || 0));
  const featureLabel = String(payload.feature_label || "gene");
  const limit = Math.max(1, Number(payload.max_count) || 1);
  const height = Math.max(360, Math.min(2400, rows.length * 26 + 150));
  // Absolute tick labels: the sign of a bar only says which side of zero it is on.
  // About four ticks per side, on a 1/2/5 x 10^k step.
  const tickTarget = Math.max(1, limit / 4);
  const tickExponent = Math.floor(Math.log10(tickTarget));
  let tickStep = Math.pow(10, tickExponent);
  for (const multiple of [1, 2, 5, 10]) {
    tickStep = multiple * Math.pow(10, tickExponent);
    if (tickStep >= tickTarget) break;
  }
  const tickCount = Math.ceil(limit / tickStep);
  const tickVals = [];
  for (let index = -tickCount; index <= tickCount; index += 1) {
    tickVals.push(index * tickStep);
  }
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  Plotly.newPlot(
    plot,
    [
      {
        type: "bar",
        orientation: "h",
        name: "Upregulated",
        x: upValues,
        y: labels,
        marker: { color: DEG_COUNT_UP_COLOR },
        customdata: ordered.map((row) => [row.up, row.down]),
        hovertemplate: `%{y}<br>%{customdata[0]} up ${featureLabel}s<extra></extra>`,
      },
      {
        type: "bar",
        orientation: "h",
        name: "Downregulated",
        x: downValues,
        y: labels,
        marker: { color: DEG_COUNT_DOWN_COLOR },
        customdata: ordered.map((row) => [row.up, row.down]),
        hovertemplate: `%{y}<br>%{customdata[1]} down ${featureLabel}s<extra></extra>`,
      },
    ],
    {
      barmode: "relative",
      bargap: 0.35,
      title: `${payload.case_label || "Group 1"} versus ${payload.control_label || "Group 2"}`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      height,
      margin: { t: 56, l: 170, r: 30, b: 60 },
      legend: { orientation: "h", x: 1, xanchor: "right", y: 1.06 },
      xaxis: {
        title: `Number of differential ${featureLabel}s`,
        range: [-limit * 1.08, limit * 1.08],
        zeroline: true,
        zerolinecolor: "#000000",
        zerolinewidth: 1,
        tickvals: tickVals,
        ticktext: tickVals.map((value) => String(Math.abs(value))),
      },
      yaxis: { automargin: true, tickfont: { size: 11 }, type: "category" },
    },
    { responsive: true }
  );
  plot.on("plotly_click", (event) => {
    const state = event?.points?.[0]?.y;
    if (!state) {
      return;
    }
    jumpToDifferentialVolcano(String(state));
  });
}

// Selecting a cell state in the counts chart opens that state in the Volcano view.
function jumpToDifferentialVolcano(population) {
  const populationSelect = document.getElementById("differential-result-population");
  const modeSelect = document.getElementById("differential-viz-mode");
  currentDifferentialPopulation = population;
  currentDifferentialGene = "";
  modeSelect.value = "volcano";
  syncDifferentialPopulationSelect(currentDifferentialState);
  if (populationSelect.value !== population) {
    const known = Array.from(populationSelect.options).some((option) => option.value === population);
    if (!known) {
      renderDifferentialEmpty(`The volcano view holds no result for ${population}.`);
      return;
    }
    populationSelect.value = population;
    currentDifferentialPopulation = population;
  }
  updateDifferentialDownloadButton();
  loadDifferentialVisualization();
}

function renderDifferentialHeatmap(payload, geneFilter = null) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  const allRows = payload.rows || [];
  const rows = geneFilter ? allRows.filter((row) => geneFilter.genes.has(String(row.gene))) : allRows;
  setDifferentialGeneFilterNote(
    describeDifferentialGeneFilter(geneFilter, payload.population, rows.length, allRows.length, "rows"),
  );
  markDifferentialGeneFilterMatched(rows.length > 0);
  if (!rows.length) {
    renderDifferentialEmpty(
      geneFilter
        ? `No heatmap row for ${geneFilter.gene} or its interacting genes in ${payload.population}.`
        : `No heatmap rows were found for ${payload.population}.`,
    );
    return;
  }
  const z = rows.map((row) => row.values);
  // Keep missing folds missing. A separate black layer supplies their colour and
  // tooltip without replacing them with a measured zero on the fold-change scale.
  const missing = z.map(values => values.map(value =>
    value === null || value === undefined || !Number.isFinite(Number(value)) ? 1 : null));
  const finiteValues = z
    .flat()
    .filter((value) => value !== null && Number.isFinite(Number(value)))
    .map((value) => Number(value));
  const [vmin, vmax] = finiteValues.length ? finiteExtent(finiteValues, -1, 1) : [-1, 1];
  const contrastFactor = 3.0;
  const maxAbs = Math.max(Math.abs(vmin), Math.abs(vmax));
  let colorExtent = maxAbs / contrastFactor;
  if (!(colorExtent > 0)) {
    colorExtent = 1;
  }
  const colorMin = -colorExtent;
  const colorMax = colorExtent;
  const modality = differentialFeatureModality();
  const hover = rows.map(row => payload.columns.map(() => featureHover(row.gene, modality)));
  const y = rows.map((row) => row.gene);
  const directions = rows.map((row) => payload.columns.map(() => row.direction));
  const figureHeight = Math.max(560, Math.min(2800, rows.length * 18 + 180));
  const tickStep = Math.max(1, Math.ceil(rows.length / ((figureHeight - 180) / 16)));
  const tickGenes = y.filter((_, i) => i % tickStep === 0);
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  Plotly.newPlot(
    plot,
    [
      {
        z: missing,
        text: hover,
        x: payload.columns,
        y,
        type: "heatmap",
        colorscale: [[0, "#000000"], [1, "#000000"]],
        zmin: 0,
        zmax: 1,
        showscale: false,
        hoverongaps: false,
        hovertemplate: "%{text}<br>%{x}<br>log2FC: not retained in the results<extra></extra>",
      },
      {
        z,
        text: hover,
        x: payload.columns,
        y,
        type: "heatmap",
        colorscale: [
          [0.0, "#00f0ff"],
          [0.5, "#000000"],
          [1.0, "#ffff00"],
        ],
        zmin: colorMin,
        zmax: colorMax,
        zmid: 0,
        customdata: directions,
        hoverongaps: false,
        hovertemplate: "%{text}<br>%{x}<br>log2FC=%{z:.3f}<br>%{customdata}<extra></extra>",
        colorbar: {
          title: "log2FC",
          tickmode: "array",
          tickvals: [colorMin, 0, colorMax],
          ticktext: [colorMin.toFixed(2), "0", colorMax.toFixed(2)],
        },
      },
    ],
    {
      title: `Heatmap: ${payload.population}`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      height: figureHeight,
      margin: { t: 56, l: 170, r: 30, b: 110 },
      xaxis: { tickangle: -40, automargin: true },
      yaxis: { automargin: true, tickfont: { size: 10 }, tickvals: tickGenes, ticktext: tickGenes.map(g => featureDisplayName(g, modality)) },
    },
    { responsive: true }
  );
  plot.on("plotly_click", (event) => {
    const gene = event?.points?.[0]?.y;
    if (gene) {
      currentDifferentialGene = gene;
      loadDifferentialGeneDetail(gene, payload.population);
    }
  });
}

function renderDifferentialVolcano(payload, geneFilter = null) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  const modality = differentialFeatureModality();
  const statisticLabel = payload.statistic_label || "FDR";
  const allPoints = payload.points || [];
  const points = geneFilter
    ? allPoints.filter((point) => geneFilter.genes.has(String(point.gene)))
    : allPoints;
  setDifferentialGeneFilterNote(
    describeDifferentialGeneFilter(geneFilter, payload.population, points.length, allPoints.length, "genes"),
  );
  markDifferentialGeneFilterMatched(points.length > 0);
  if (!points.length) {
    renderDifferentialEmpty(
      geneFilter
        ? `No volcano point for ${geneFilter.gene} or its interacting genes in ${payload.population}.`
        : `No volcano data were found for ${payload.population}.`,
    );
    return;
  }
  const up = points.filter((point) => point.direction === "up");
  const down = points.filter((point) => point.direction === "down");
  // A filtered volcano holds few enough points to name every one of them on the plot.
  const pointMode = geneFilter && points.length <= 60 ? "markers+text" : "markers";
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  Plotly.newPlot(
    plot,
    [
      {
        x: down.map((point) => point.log2fc),
        y: down.map((point) => point.score),
        text: down.map((point) => featureDisplayName(point.gene, modality)),
        hovertext: down.map(point => featureHover(point.gene, modality)),
        customdata: down.map((point) => [point.gene, point.fdr, point.pval]),
        type: "scattergl",
        mode: pointMode,
        textposition: "top center",
        textfont: { size: 10, color: "#1e293b" },
        name: "Down",
        marker: { color: "#2563eb", size: 16, opacity: 0.72 },
        hovertemplate: `%{hovertext}<br>log2FC=%{x:.3f}<br>-log10(${statisticLabel})=%{y:.3f}<extra></extra>`,
      },
      {
        x: up.map((point) => point.log2fc),
        y: up.map((point) => point.score),
        text: up.map((point) => featureDisplayName(point.gene, modality)),
        hovertext: up.map(point => featureHover(point.gene, modality)),
        customdata: up.map((point) => [point.gene, point.fdr, point.pval]),
        type: "scattergl",
        mode: pointMode,
        textposition: "top center",
        textfont: { size: 10, color: "#1e293b" },
        name: "Up",
        marker: { color: "#dc2626", size: 16, opacity: 0.72 },
        hovertemplate: `%{hovertext}<br>log2FC=%{x:.3f}<br>-log10(${statisticLabel})=%{y:.3f}<extra></extra>`,
      },
    ],
    {
      title: `Volcano: ${payload.population}`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      margin: { t: 56, l: 60, r: 20, b: 56 },
      height: 640,
      xaxis: { title: "log2 fold change", zeroline: true, zerolinecolor: "rgba(100,116,139,0.45)" },
      yaxis: { title: `-log10(${statisticLabel})` },
      hovermode: "closest",
    },
    { responsive: true }
  );
  plot.on("plotly_click", (event) => {
    const gene = event?.points?.[0]?.customdata?.[0];
    if (gene) {
      currentDifferentialGene = gene;
      loadDifferentialGeneDetail(gene, payload.population);
    }
  });
}

function renderDifferentialGo(payload, geneFilter = null) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  const allTerms = payload.terms || [];
  // A GO term survives the filter when its own overlapping genes meet the filter set.
  const terms = geneFilter
    ? allTerms.filter((term) => (term.overlap_genes || []).some((gene) => geneFilter.genes.has(String(gene))))
    : allTerms;
  setDifferentialGeneFilterNote(
    describeDifferentialGeneFilter(geneFilter, payload.population, terms.length, allTerms.length, "GO terms"),
  );
  markDifferentialGeneFilterMatched(terms.length > 0);
  if (!terms.length) {
    renderDifferentialEmpty(
      geneFilter
        ? `No GO term for ${geneFilter.gene} or its interacting genes in ${payload.population}.`
        : `No GO term enrichment was available for ${payload.population}.`,
    );
    return;
  }
  const ordered = [...terms].map((term) => {
    const overlapGenes = term.overlap_genes || [];
    // With a filter on, the typed gene names the term wherever the term contains it.
    const filteredPick = geneFilter && overlapGenes.includes(geneFilter.gene) ? geneFilter.gene : null;
    const selectedGene = filteredPick
      || (overlapGenes.includes(currentDifferentialGene)
        ? currentDifferentialGene
        : (term.selected_gene || overlapGenes[0] || null));
    return {
      ...term,
      selected_gene: selectedGene,
      z_score: Number(term.z_score),
      fdr_plot: Number(term.fdr_plot || term.fdr || term.p_value),
      p_value: Number(term.p_value),
      fdr: Number(term.fdr),
    };
  }).filter((term) => Number.isFinite(term.z_score) && Number.isFinite(term.fdr_plot) && term.fdr_plot > 0);
  if (!ordered.length) {
    renderDifferentialEmpty(`No GO term enrichment was available for ${payload.population}.`);
    return;
  }
  const significant = ordered.filter((term) => Boolean(term.is_selected_positive_sig));
  const background = ordered.filter((term) => !term.is_selected_positive_sig);
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  Plotly.newPlot(
    plot,
    [
      {
        x: background.map((term) => term.z_score),
        y: background.map((term) => term.fdr_plot),
        type: "scattergl",
        mode: "markers",
        name: "Other terms",
        customdata: background.map((term) => [
          term.selected_gene,
          term.direction,
          (term.overlap_genes || []).join(", "),
          term.overlap_genes || [],
          term.term_name,
          term.z_score,
          term.p_value,
          term.fdr,
        ]),
        marker: {
          color: "#d1d5db",
          size: 11,
          opacity: 0.95,
        },
        hovertemplate: (
          "%{customdata[4]}<br>Z-score=%{x:.3f}<br>FDR=%{y:.3e}<br>Fisher p=%{customdata[6]:.3e}"
          + "<br>Selected gene: %{customdata[0]}"
          + "<br>GO genes: %{customdata[2]}<extra></extra>"
        ),
      },
      {
        x: significant.map((term) => term.z_score),
        y: significant.map((term) => term.fdr_plot),
        type: "scattergl",
        mode: "markers",
        name: "GO-Elite selected terms",
        customdata: significant.map((term) => [
          term.selected_gene,
          term.direction,
          (term.overlap_genes || []).join(", "),
          term.overlap_genes || [],
          term.term_name,
          term.z_score,
          term.p_value,
          term.fdr,
        ]),
        marker: {
          color: "#1f19c7",
          size: 11,
          opacity: 0.98,
        },
        hovertemplate: (
          "%{customdata[4]}<br>Z-score=%{x:.3f}<br>FDR=%{y:.3e}<br>Fisher p=%{customdata[6]:.3e}"
          + "<br>Selected gene: %{customdata[0]}"
          + "<br>GO genes: %{customdata[2]}<extra></extra>"
        ),
      },
    ],
    {
      title: `GO terms: ${payload.population}`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      margin: { t: 56, l: 86, r: 72, b: 72 },
      height: 720,
      xaxis: {
        title: "Z-Score",
        zeroline: false,
        range: [
          Math.min(-10, Math.floor(Math.min(...ordered.map((term) => term.z_score)) - 0.5)),
          Math.max(20, Math.ceil(Math.max(...ordered.map((term) => term.z_score)) + 2.5)),
        ],
      },
      yaxis: {
        title: "Fishers FDR p",
        type: "log",
        autorange: true,
      },
      showlegend: false,
    },
    { responsive: true }
  );
  plot.on("plotly_click", (event) => {
    const selectedGene = event?.points?.[0]?.customdata?.[0];
    const overlapGenes = event?.points?.[0]?.customdata?.[3] || [];
    let gene = selectedGene;
    if (Array.isArray(overlapGenes) && overlapGenes.length > 1) {
      const response = window.prompt(
        `Select a gene from this GO term:\n${overlapGenes.join(", ")}`,
        selectedGene || overlapGenes[0] || ""
      );
      if (response === null) {
        return;
      }
      const trimmed = response.trim();
      if (trimmed && overlapGenes.includes(trimmed)) {
        gene = trimmed;
      }
    }
    if (gene) {
      currentDifferentialGene = gene;
      loadDifferentialGeneDetail(gene, payload.population);
    }
  });
}

function renderDifferentialNetwork(payload, geneFilter = null) {
  if (payload && payload.network_type === "cell_communication_diff") {
    renderDifferentialCommunicationNetwork(payload, geneFilter);
    return;
  }
  const plot = document.getElementById("differential-plot-area");
  Plotly.purge(plot);
  destroyDifferentialNetwork();
  // With a filter on, keep the typed gene and the nodes it connects to, and keep only
  // the edges whose two ends both survive.
  const allElements = payload.elements || [];
  const keptElements = geneFilter
    ? allElements.filter((element) => {
      const data = element && element.data ? element.data : null;
      if (!data) {
        return false;
      }
      if (data.source && data.target) {
        return geneFilter.genes.has(String(data.source)) && geneFilter.genes.has(String(data.target));
      }
      return data.id ? geneFilter.genes.has(String(data.id)) : false;
    })
    : allElements;
  const elements = keptElements.map((element) => {
    if (!element.data || !element.data.id || element.data.source) {
      return element;
    }
    return {
      data: {
        ...element.data,
        color: networkFoldColor(element.data.log2fc),
      },
    };
  });
  const countNodes = (list) => list.filter((element) => element && element.data
    && element.data.id && !element.data.source).length;
  setDifferentialGeneFilterNote(
    describeDifferentialGeneFilter(
      geneFilter, payload.population, countNodes(keptElements), countNodes(allElements), "nodes"),
  );
  markDifferentialGeneFilterMatched(countNodes(keptElements) > 0);
  if (!elements.length) {
    renderDifferentialEmpty(
      geneFilter
        ? `${geneFilter.gene} has no interaction network node in ${payload.population}.`
        : `No interaction network was available for ${payload.population}.`,
    );
    return;
  }
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  differentialCy = cytoscape({
    container: plot,
    elements,
    style: [
      {
        selector: "node",
        style: {
          "background-color": "data(color)",
          label: "data(label)",
          color: "#0f172a",
          "font-size": 12,
          "text-valign": "center",
          "text-halign": "center",
          width: 26,
          height: 26,
        },
      },
      {
        selector: "edge",
        style: {
          // GRN edges carry a score (regulatory activity) -> encode as width so cell-state
          // differences are visible (e.g. HLF edges thick in HSC, absent in monocytes).
          width: (edge) => {
            const s = Math.abs(Number(edge.data("score")) || 0);
            return s > 0 ? Math.max(1, Math.min(9, 1 + s * 24)) : 1.8;
          },
          "line-color": networkEdgeColor,
          "target-arrow-color": networkEdgeColor,
          "target-arrow-shape": networkEdgeArrowShape,
          "curve-style": "bezier",
          opacity: 0.85,
        },
      },
      {
        selector: "node:selected",
        style: {
          "border-width": 3,
          "border-color": "#0f172a",
        },
      },
    ],
    layout: {
      name: "cose",
      animate: true,
      fit: true,
      padding: 36,
      randomize: true,
      idealEdgeLength: 80,
      nodeOverlap: 8,
      componentSpacing: 90,
    },
  });
  differentialCy.on("tap", "node", (event) => {
    const gene = event?.target?.data("id");
    if (gene) {
      currentDifferentialGene = gene;
      loadDifferentialGeneDetail(gene, payload.population);
    }
  });
}

function setDifferentialNetworkHoverTooltip(text, renderedPosition = null) {
  const plot = document.getElementById("differential-plot-area");
  if (!plot) {
    return;
  }
  let tooltip = plot.querySelector(".network-hover-tooltip");
  if (!tooltip) {
    tooltip = document.createElement("div");
    tooltip.className = "network-hover-tooltip";
    plot.appendChild(tooltip);
  }
  const message = String(text || "").trim();
  if (!message) {
    tooltip.classList.add("hidden");
    tooltip.textContent = "";
    return;
  }
  tooltip.textContent = message;
  tooltip.classList.remove("hidden");
  if (renderedPosition) {
    tooltip.style.left = `${Math.min(plot.clientWidth - 260, Math.max(12, renderedPosition.x + 12))}px`;
    tooltip.style.top = `${Math.min(plot.clientHeight - 120, Math.max(12, renderedPosition.y + 12))}px`;
  }
}

function renderDifferentialCommunicationNetwork(payload, geneFilter = null) {
  const plot = document.getElementById("differential-plot-area");
  Plotly.purge(plot);
  destroyDifferentialNetwork();
  const allElements = payload.elements || [];
  // The nodes are cell states, so the filter works on the edges: keep an edge that uses
  // the typed gene as its ligand or its receptor, then keep the states those edges join
  // plus the focus state, so the surviving edges still have both ends.
  let elements = allElements;
  if (geneFilter) {
    const keptEdges = allElements.filter((element) => {
      const data = element && element.data ? element.data : null;
      return Boolean(data && data.source && data.target)
        && communicationEdgeFeatures(data).some((value) => geneFilter.genes.has(value));
    });
    const keepNodes = new Set();
    keptEdges.forEach((edge) => {
      keepNodes.add(String(edge.data.source));
      keepNodes.add(String(edge.data.target));
    });
    allElements.forEach((element) => {
      const data = element && element.data ? element.data : null;
      if (data && !data.source && data.node_type === "focus" && data.id) {
        keepNodes.add(String(data.id));
      }
    });
    elements = allElements.filter((element) => {
      const data = element && element.data ? element.data : null;
      if (!data) return false;
      if (data.source && data.target) return keptEdges.indexOf(element) >= 0;
      return data.id ? keepNodes.has(String(data.id)) : false;
    });
    const totalEdges = allElements.filter((e) => e && e.data && e.data.source).length;
    setDifferentialGeneFilterNote(describeCommunicationGeneFilter(
      geneFilter, payload.population, keptEdges.length, totalEdges, "interactions"));
    markDifferentialGeneFilterMatched(keptEdges.length > 0);
    // No edge survived. The focus state alone is not a network, so say so instead of
    // drawing one lonely node.
    if (!keptEdges.length) {
      elements = [];
    }
  } else {
    setDifferentialGeneFilterNote("");
    markDifferentialGeneFilterMatched(true);
  }
  if (!elements.length) {
    renderDifferentialEmpty(
      geneFilter
        ? `No cell communication interaction uses ${geneFilter.gene} in ${payload.population}.`
        : `No differential cell-state network was available for ${payload.population}.`,
    );
    return;
  }
  const usePresetLayout = elements.some((element) => element && !element.data?.source && element.position);
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  differentialCy = cytoscape({
    container: plot,
    elements,
    style: [
      {
        selector: "node",
        style: {
          "background-color": "data(color)",
          label: "data(label)",
          color: "#0f172a",
          "font-size": 12,
          "text-valign": "center",
          "text-halign": "center",
          width: 34,
          height: 34,
          "border-width": 1,
          "border-color": "#94a3b8",
        },
      },
      {
        selector: "node[node_type = 'focus']",
        style: {
          width: 54,
          height: 54,
          color: "#ffffff",
          "font-weight": 700,
          "border-width": 3,
          "border-color": "#134e4a",
        },
      },
      {
        selector: "edge",
        style: {
          width: "data(weight)",
          "line-color": "data(edge_color)",
          "target-arrow-color": "data(edge_color)",
          "target-arrow-shape": "triangle",
          "curve-style": "bezier",
          opacity: "data(edge_opacity)",
        },
      },
      {
        selector: "edge.hovered",
        style: {
          opacity: 1,
          label: "data(label)",
          color: "#0f172a",
          "font-size": 10,
          "text-background-color": "#ffffff",
          "text-background-opacity": 0.88,
          "text-background-padding": 3,
        },
      },
      {
        selector: "node:selected",
        style: {
          "border-width": 4,
          "border-color": "#f97316",
        },
      },
    ],
    layout: usePresetLayout
      ? {
        name: "preset",
        fit: true,
        padding: 56,
      }
      : {
        name: "circle",
        animate: true,
        fit: true,
        padding: 56,
      },
  });
  differentialCy.on("mouseover", "edge", (event) => {
    event.target.addClass("hovered");
    setDifferentialNetworkHoverTooltip(event.target.data("tooltip"), event.renderedPosition);
  });
  differentialCy.on("mouseout", "edge", (event) => {
    event.target.removeClass("hovered");
    setDifferentialNetworkHoverTooltip("");
  });
  differentialCy.on("tap", "edge", (event) => {
    const gene = event?.target?.data("gene");
    event.target.addClass("hovered");
    setDifferentialNetworkHoverTooltip(event.target.data("tooltip"), event.renderedPosition);
    if (!gene) {
      return;
    }
    currentDifferentialGene = gene;
    const topInteractions = event.target.data("top_interactions") || [];
    const top = topInteractions.find((entry) => String(entry.gene || "") === String(gene)) || topInteractions[0] || {};
    const ligand = String(top.ligand || "").trim();
    const receptor = String(top.receptor || "").trim();
    const role = currentDifferentialFeatureRole === "receptor" && receptor ? "receptor" : "ligand";
    const featureSymbol = role === "ligand" ? (ligand || receptor) : (receptor || ligand);
    currentDifferentialFeatureRole = role;
    loadDifferentialGeneDetail(gene, payload.population, featureSymbol ? { feature: featureSymbol } : undefined);
  });
}

function renderDifferentialCommunicationTable(payload, geneFilter = null) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  const allRows = payload.rows || [];
  // One row is one ligand-receptor interaction, so the filter keeps a row whose ligand
  // or whose receptor is the typed gene.
  const rows = geneFilter
    ? allRows.filter((row) => communicationRowMatches(row, geneFilter))
    : allRows;
  const columns = payload.columns || [];
  setDifferentialGeneFilterNote(describeCommunicationGeneFilter(
    geneFilter, payload.population, rows.length, allRows.length, "interactions"));
  markDifferentialGeneFilterMatched(rows.length > 0);
  if (!rows.length || !columns.length) {
    renderDifferentialEmpty(
      geneFilter && allRows.length
        ? `No cell communication interaction uses ${geneFilter.gene} in ${payload.population}.`
        : `No differential interaction table was available for ${payload.population}.`,
    );
    return;
  }
  const labels = columns.map((column) => column
    .replace("delta_score", "delta score")
    .replace("case_mean_score", "case mean")
    .replace("control_mean_score", "control mean")
    .replace("lr_expression_score", "LR expression")
    .replace("receiver_response_score", "response")
    .replaceAll("_", " "));
  const values = columns.map((column) => rows.map((row) => {
    const value = row[column];
    return Number.isFinite(Number(value)) && String(value).trim() !== "" ? Number(value).toFixed(3) : String(value ?? "");
  }));
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  Plotly.newPlot(plot, [{
    type: "table",
    header: {
      values: labels,
      align: "left",
      fill: { color: "#0f766e" },
      font: { color: "white", size: 11 },
    },
    cells: {
      values,
      align: "left",
      fill: { color: rows.map((_, index) => (index % 2 ? "#f8fafc" : "#ffffff")) },
      font: { color: "#0f172a", size: 10 },
      height: 24,
    },
  }], {
    title: `Top differential cell-communication interactions for ${payload.population}`,
    paper_bgcolor: "rgba(0,0,0,0)",
    height: Math.max(520, Math.min(920, rows.length * 26 + 110)),
    margin: { t: 46, l: 12, r: 12, b: 12 },
  }, { responsive: true });
  plot.on("plotly_click", (event) => {
    const pointNumber = event?.points?.[0]?.pointNumber;
    const rowIndex = Array.isArray(pointNumber) ? Number(pointNumber[0]) : Number(pointNumber);
    const row = Number.isInteger(rowIndex) ? rows[rowIndex] : null;
    const gene = row && row.gene ? String(row.gene) : "";
    if (!gene) {
      return;
    }
    currentDifferentialGene = gene;
    const ligand = row && row.ligand ? String(row.ligand) : "";
    const receptor = row && row.receptor ? String(row.receptor) : "";
    const role = currentDifferentialFeatureRole === "receptor" && receptor ? "receptor" : "ligand";
    const featureSymbol = role === "ligand" ? (ligand || receptor) : (receptor || ligand);
    currentDifferentialFeatureRole = role;
    loadDifferentialGeneDetail(gene, payload.population, featureSymbol ? { feature: featureSymbol } : undefined);
  });
}

async function loadDifferentialGeneDetail(gene, population, options) {
  const request = ++differentialDetailRequest;
  const runId = currentDifferentialState?.run_id;
  const jobId = document.getElementById("results-job-id").value.trim();
  const isCurrent = () => request === differentialDetailRequest &&
    runId === currentDifferentialState?.run_id &&
    jobId === document.getElementById("results-job-id").value.trim();
  if (!gene || !population || !jobId) {
    resetDifferentialGeneDetail();
    return;
  }
  const opts = options || {};
  let url = apiPath(`/jobs/${jobId}/differential/interactive/gene?population=${encodeURIComponent(population)}&gene=${encodeURIComponent(gene)}`);
  if (opts.feature) {
    url += `&feature=${encodeURIComponent(opts.feature)}`;
  }
  try {
    const payload = await fetchDifferentialJson(url);
    if (!isCurrent()) return;
    currentDifferentialPopulation = payload.population || population;
    if (payload.interaction_key || payload.ligand || payload.receptor) {
      currentDifferentialInteraction = {
        population: payload.population || population,
        interaction_key: payload.interaction_key || gene,
        interaction: payload.interaction || payload.gene,
        ligand: payload.ligand || "",
        receptor: payload.receptor || "",
        sender_state: payload.sender_state || "",
        receiver_state: payload.receiver_state || payload.population || population,
      };
      if (payload.feature_role === "ligand" || payload.feature_role === "receptor") {
        currentDifferentialFeatureRole = payload.feature_role;
      }
    } else {
      currentDifferentialInteraction = null;
    }
    renderDifferentialGeneDetail(payload);
  } catch (err) {
    if (!isCurrent()) return;
    resetDifferentialGeneDetail();
    document.getElementById("differential-selected-gene").textContent = gene;
    document.getElementById("differential-gene-empty").textContent = err.message || "Unable to load gene detail.";
    document.getElementById("differential-gene-empty").classList.remove("hidden");
  }
}

function renderDifferentialGeneDetail(payload) {
  const modality = differentialFeatureModality();
  const displayGene = featureDisplayName(payload.gene, modality);
  const plot = document.getElementById("differential-gene-plot");
  const empty = document.getElementById("differential-gene-empty");
  const stats = document.getElementById("differential-gene-stats");
  const downloadButton = document.getElementById("download-differential-gene-btn");
  const isCommunication = isCurrentDifferentialCommunicationAnalysis(payload);
  const isFeatureExpression = payload.view_kind === "feature_expression";
  const yTitle = payload.value_label || (isFeatureExpression
    ? "Normalized expression"
    : (isCommunication ? "fastComm communication score" : "Normalized expression"));
  const headerLabel = isFeatureExpression
    ? `${displayGene} (${payload.feature_role || "feature"}) in ${payload.feature_state || payload.population}`
    : `${displayGene} in ${payload.population}`;
  document.getElementById("differential-selected-gene").textContent = headerLabel;
  updateDifferentialLrToggle(payload);
  plot.classList.remove("hidden");
  empty.classList.add("hidden");
  downloadButton.href = "#";
  downloadButton.classList.remove("hidden");

  const traces = (payload.groups || []).map((group, index) => ({
    type: "violin",
    name: group.label,
    y: group.values,
    points: "all",
    pointpos: 0,
    jitter: 0.28,
    marker: {
      size: 2,
      opacity: 0.35,
      color: index === 0 ? "#dc2626" : "#2563eb",
    },
    line: {
      color: index === 0 ? "#dc2626" : "#2563eb",
    },
    box: { visible: true },
    meanline: { visible: true },
    hovertemplate: `${featureHover(payload.gene, modality)}<br>${group.label}<br>${isFeatureExpression ? "expr" : (isCommunication ? "score" : "expr")}=%{y:.3f}<extra></extra>`,
  }));
  const titleSubject = isFeatureExpression
    ? `${displayGene} in ${payload.feature_state || payload.population}`
    : displayGene;
  Plotly.newPlot(
    plot,
    traces,
    {
      // The y axis already names the measure, so the title carries only the subject.
      // Nathan, 2026-09-01: "Just make it the gene name. The prefix is unnecessary."
      title: {text: featureAnnotation(payload.gene, modality) ? featurePlotTitle(payload.gene, modality) : titleSubject, font:{size:14}},
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      margin: { t: featureAnnotation(payload.gene, modality) ? 90 : 52, l: 50, r: 18, b: 48 },
      height: 420,
      yaxis: { title: yTitle },
      xaxis: { automargin: true },
      showlegend: false,
    },
    { responsive: true }
  );
  stats.classList.remove("hidden");
  const pValue = formatDifferentialStat(payload.stats?.p_value);
  const fdr = formatDifferentialStat(payload.stats?.fdr);
  const log2fc = formatDifferentialStat(payload.stats?.log2fc);
  stats.innerHTML = `
    <div class="gene-stat-grid">
      <div class="gene-stat">
        <div class="gene-stat-label">p-value</div>
        <div class="gene-stat-value">${pValue}</div>
      </div>
      <div class="gene-stat">
        <div class="gene-stat-label">FDR</div>
        <div class="gene-stat-value">${fdr}</div>
      </div>
      <div class="gene-stat">
        <div class="gene-stat-label">log2FC</div>
        <div class="gene-stat-value">${log2fc}</div>
      </div>
    </div>
  `;
}

function hideDifferentialLrToggle() {
  const toggle = document.getElementById("differential-lr-toggle");
  if (!toggle) {
    return;
  }
  toggle.classList.add("hidden");
  toggle.hidden = true;
  toggle.setAttribute("aria-hidden", "true");
  toggle.querySelectorAll(".lr-toggle-btn").forEach((btn) => {
    btn.classList.add("hidden");
    btn.disabled = true;
    btn.classList.remove("active");
    btn.textContent = btn.getAttribute("data-lr-role") === "receptor" ? "Receptor" : "Ligand";
  });
  currentDifferentialInteraction = null;
}

function isCurrentDifferentialCommunicationAnalysis(payload) {
  const payloadModality = normalizeModalityId(payload && payload.modality, "");
  const configuredModality = normalizeModalityId(
    (currentDifferentialState && currentDifferentialState.config && currentDifferentialState.config.modality)
      || document.getElementById("differential-modality")?.value
      || "",
    ""
  );
  return payloadModality === "cell_communication" && configuredModality === "cell_communication";
}

function updateDifferentialLrToggle(payload) {
  const toggle = document.getElementById("differential-lr-toggle");
  if (!toggle) {
    return;
  }
  if (!isCurrentDifferentialCommunicationAnalysis(payload)) {
    hideDifferentialLrToggle();
    return;
  }
  const ligand = String((payload && payload.ligand) || "").trim();
  const receptor = String((payload && payload.receptor) || "").trim();
  if (!ligand && !receptor) {
    hideDifferentialLrToggle();
    return;
  }
  toggle.hidden = false;
  toggle.removeAttribute("aria-hidden");
  toggle.classList.remove("hidden");
  const buttons = toggle.querySelectorAll(".lr-toggle-btn");
  const activeRole = payload && payload.feature_role
    ? payload.feature_role
    : currentDifferentialFeatureRole;
  buttons.forEach((btn) => {
    const role = btn.getAttribute("data-lr-role");
    const symbol = role === "ligand" ? ligand : receptor;
    if (!symbol) {
      btn.classList.add("hidden");
      btn.disabled = true;
      btn.textContent = role === "ligand" ? "Ligand" : "Receptor";
    } else {
      btn.classList.remove("hidden");
      btn.disabled = false;
      btn.textContent = `${role === "ligand" ? "Ligand" : "Receptor"}: ${symbol}`;
    }
    btn.classList.toggle("active", role === activeRole);
  });
}

function attachDifferentialLrToggleHandlers() {
  const toggle = document.getElementById("differential-lr-toggle");
  if (!toggle || toggle.dataset.bound === "true") {
    return;
  }
  toggle.dataset.bound = "true";
  toggle.addEventListener("click", (event) => {
    const btn = event.target.closest(".lr-toggle-btn");
    if (!btn || btn.disabled) {
      return;
    }
    const role = btn.getAttribute("data-lr-role");
    if (!role || !currentDifferentialInteraction) {
      return;
    }
    const ctx = currentDifferentialInteraction;
    const featureSymbol = role === "ligand" ? ctx.ligand : ctx.receptor;
    if (!featureSymbol) {
      return;
    }
    currentDifferentialFeatureRole = role;
    loadDifferentialGeneDetail(ctx.interaction_key, ctx.population, { feature: featureSymbol });
  });
}

function formatDifferentialStat(value) {
  if (value === null || value === undefined || Number.isNaN(Number(value))) {
    return "NA";
  }
  const numeric = Number(value);
  if (Math.abs(numeric) > 0 && Math.abs(numeric) < 1e-3) {
    return numeric.toExponential(2);
  }
  return numeric.toFixed(3);
}

function updateDifferentialDownloadButton() {
  const state = currentDifferentialState;
  const button = document.getElementById("download-differential-left-btn");
  if (!state || state.status !== "completed") {
    button.classList.add("hidden");
    button.href = "#";
    return;
  }
  const mode = document.getElementById("differential-viz-mode").value;

  const population = document.getElementById("differential-result-population").value;
  if (!population && mode !== "summary") {
    button.classList.add("hidden");
    button.href = "#";
    return;
  }
  button.href = "#";
  button.classList.remove("hidden");
}

function destroyDifferentialNetwork() {
  if (differentialCy) {
    differentialCy.destroy();
    differentialCy = null;
  }
  setDifferentialNetworkHoverTooltip("");
}

function renderDifferentialEmpty(message) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  Plotly.purge(plot);
  plot.classList.add("hidden");
  const empty = document.getElementById("differential-plot-empty");
  empty.textContent = message;
  empty.classList.remove("hidden");
}

function resetDifferentialResults() {
  destroyDifferentialNetwork();
  const vizModes = ((currentDifferentialState && currentDifferentialState.visualization_modes) || []).map((entry) => entry.label);
  const vizSummary = vizModes.length ? vizModes.join(", ") : "Heatmap, Volcano, Network, or GO terms";
  const plot = document.getElementById("differential-plot-area");
  Plotly.purge(plot);
  plot.classList.add("hidden");
  const empty = document.getElementById("differential-plot-empty");
  empty.textContent = `Run cellHarmony-differential to explore ${vizSummary}.`;
  empty.classList.remove("hidden");
  document.getElementById("download-differential-left-btn").classList.add("hidden");
  const populationSelect = document.getElementById("differential-result-population");
  populationSelect.innerHTML = "";
  currentDifferentialPopulation = "";
  currentDifferentialGene = "";
  // A new differential run invalidates the gene filter and every cached adjacency.
  clearDifferentialGeneFilter();
  syncDifferentialGeneFilterLabel();
  resetDifferentialGeneDetail();
}

function resetDifferentialGeneDetail() {
  ++differentialDetailRequest;
  const featureLabel = String((currentDifferentialState && currentDifferentialState.feature_label) || "gene");
  const isGrn = featureLabel === "TF";   // GRN is the only TF-feature modality
  const detailNoun = isGrn ? "edge" : featureLabel;
  const capitalizedFeature = detailNoun.charAt(0).toUpperCase() + detailNoun.slice(1);
  const plot = document.getElementById("differential-gene-plot");
  Plotly.purge(plot);
  plot.classList.add("hidden");
  document.getElementById("differential-selected-gene").textContent = `Select a ${detailNoun} from the left panel.`;
  document.getElementById("download-differential-gene-btn").classList.add("hidden");
  document.getElementById("download-differential-gene-btn").removeAttribute("href");
  const empty = document.getElementById("differential-gene-empty");
  empty.textContent = `Select a ${detailNoun} from the differential view to compare ${detailNoun} values between groups.`;
  empty.classList.remove("hidden");
  document.getElementById("differential-gene-stats").classList.add("hidden");
  const detailTitle = document.getElementById("differential-detail-title");
  if (detailTitle) {
    detailTitle.textContent = isGrn ? "GRN Detail" : `${capitalizedFeature} Detail`;
  }
  hideDifferentialLrToggle();
}

function clearGeneSuggestions() {
  VISUALIZATION_PANELS.forEach((panelKey) => {
    const datalist = document.getElementById(`${panelKey}-feature-suggestions`);
    if (datalist) {
      datalist.innerHTML = "";
    }
  });
  loadedGeneSuggestionsSignature = "";
}

function clearDisplayFilters() {
  currentDisplayFiltersMeta = null;
  loadedDisplayFiltersJobId = null;
  exploreWarmupJobId = null;
  exploreWarmupPromise = null;
  VISUALIZATION_PANELS.forEach((panelKey) => {
    const filter1Field = document.getElementById(panelElementId(panelKey, "filter1-field"));
    const filter2Field = document.getElementById(panelElementId(panelKey, "filter2-field"));
    const filter1Values = document.getElementById(panelElementId(panelKey, "filter1-values"));
    const filter2Values = document.getElementById(panelElementId(panelKey, "filter2-values"));
    const secondaryRow = document.getElementById(panelElementId(panelKey, "filter-secondary-row"));
    if (filter1Field) {
      filter1Field.innerHTML = "";
    }
    if (filter2Field) {
      filter2Field.innerHTML = "";
    }
    if (filter1Values) {
      filter1Values.innerHTML = "";
    }
    if (filter2Values) {
      filter2Values.innerHTML = "";
    }
    if (secondaryRow) {
      secondaryRow.classList.remove("hidden");
    }
  });
}

async function loadGeneSuggestions(jobId) {
  await ensureFeatureAnnotations(jobId);
  const signature = [
    String(jobId || "").trim(),
    ...VISUALIZATION_PANELS.map((panelKey) => `${panelKey}:${panelModality(panelKey)}`),
  ].join("|");
  if (!jobId || loadedGeneSuggestionsSignature === signature) {
    return;
  }
  try {
    for (const panelKey of VISUALIZATION_PANELS) {
      const modality = panelModality(panelKey);
      const datalist = document.getElementById(`${panelKey}-feature-suggestions`);
      if (!datalist) {
        continue;
      }
      datalist.innerHTML = "";
      const resp = await fetch(apiPath(`/jobs/${jobId}/genes?modality=${encodeURIComponent(modality)}`));
      const data = await parseApiResponse(resp);
      if (!resp.ok) {
        throw new Error(data.detail || "Unable to load gene suggestions.");
      }
      const seen = new Set();
      // `accepted` is what may stay in the box; `seen` is what the box offers.
      const accepted = new Set();
      const features = [];
      (data.genes || []).forEach((gene) => {
        const value = String(gene || "").trim();
        if (!value || seen.has(value)) {
          return;
        }
        seen.add(value);
        accepted.add(value);
        features.push(value);
        const option = document.createElement("option");
        option.value = value;
        option.label = featureDisplayName(value, modality);
        datalist.appendChild(option);
      });
      // A NAME THE SERVER CAN RESOLVE IS NEVER OVERWRITTEN.
      //
      // Nathan, 2026-09-08: "every name I select defaults to CE(18:2)". The reset below
      // fires whenever the typed value is not in the suggestion list, and it replaces it
      // with the FIRST suggestion. Once lipids began showing "18:2 Cholesterol ester",
      // the list and the value could legitimately disagree: a bundle resolves both the
      // display name and the abbreviation, so `CE(20:3)` is valid and was being wiped
      // because only display names are offered.
      //
      // So the accepted set is the suggestions PLUS the keys the payload names, while the
      // datalist still shows only the reader's name. A payload without `keys` -- every
      // modality that has no display column -- leaves this exactly as it was.
      (data.keys || []).forEach((key) => {
        const value = String(key || "").trim();
        if (value) accepted.add(value);
      });
      const input = document.getElementById(panelElementId(panelKey, "gene-query"));
      if (input) {
        input.setAttribute("list", `${panelKey}-feature-suggestions`);
        const currentValue = String(input.value || "").trim();
        if (!currentValue || !accepted.has(currentValue)) {
          input.value = preferredFeatureForModality(modality, features);
        }
      }
      updatePanelFeatureInput(panelKey);
    }
    loadedGeneSuggestionsSignature = signature;
  } catch (err) {
    console.warn(err);
  }
}

function populateSelectOptions(selectEl, options, selectedValue, includeBlank = false) {
  selectEl.innerHTML = "";
  if (includeBlank) {
    const blank = document.createElement("option");
    blank.value = "";
    blank.textContent = "None";
    selectEl.appendChild(blank);
  }
  (options || []).forEach((optionData) => {
    const option = document.createElement("option");
    option.value = optionData.value;
    option.textContent = optionData.label;
    if (selectedValue && optionData.value === selectedValue) {
      option.selected = true;
    }
    selectEl.appendChild(option);
  });
}

function syncDisplayFilterValueOptions(panelKey, index) {
  const fieldSelect = document.getElementById(panelElementId(panelKey, `filter${index}-field`));
  const valueSelect = document.getElementById(panelElementId(panelKey, `filter${index}-values`));
  const field = fieldSelect.value;
  valueSelect.innerHTML = "";
  const values = (currentDisplayFiltersMeta?.values && currentDisplayFiltersMeta.values[field]) || [];
  const allOption = document.createElement("option");
  allOption.value = "";
  allOption.textContent = "All";
  valueSelect.appendChild(allOption);
  values.forEach((value) => {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = value;
    valueSelect.appendChild(option);
  });
  valueSelect.value = "";
  valueSelect.disabled = !field;
}

function populateDisplayFilterControlsForPanel(panelKey, meta) {
  const fields = (meta && meta.fields) || [];
  const primaryField = (meta && meta.default_primary_field) || "";
  const secondaryField = (meta && meta.default_secondary_field) || "";
  const showSecondary = !meta || meta.show_secondary !== false;
  populateSelectOptions(document.getElementById(panelElementId(panelKey, "filter1-field")), fields, primaryField, false);
  populateSelectOptions(document.getElementById(panelElementId(panelKey, "filter2-field")), fields, secondaryField, true);
  syncDisplayFilterValueOptions(panelKey, 1);
  syncDisplayFilterValueOptions(panelKey, 2);
  document.getElementById(panelElementId(panelKey, "filter-secondary-row")).classList.toggle("hidden", !showSecondary);
  if (!showSecondary) {
    document.getElementById(panelElementId(panelKey, "filter2-field")).value = "";
    document.getElementById(panelElementId(panelKey, "filter2-values")).innerHTML = "";
    document.getElementById(panelElementId(panelKey, "filter2-values")).disabled = true;
  }
}

function populateDisplayFilterControls(meta) {
  currentDisplayFiltersMeta = meta || null;
  VISUALIZATION_PANELS.forEach((panelKey) => {
    populateDisplayFilterControlsForPanel(panelKey, meta);
  });
  updateExpressionModeOptions();
}

async function loadDisplayFilters(jobId) {
  if (!jobId || loadedDisplayFiltersJobId === jobId) {
    return;
  }
  try {
    const resp = await fetch(apiPath(`/jobs/${jobId}/display-filters`));
    const data = await parseApiResponse(resp);
    if (!resp.ok) {
      throw new Error(data.detail || "Unable to load display filters.");
    }
    populateDisplayFilterControls(data);
    loadedDisplayFiltersJobId = jobId;
  } catch (err) {
    console.warn(err);
  }
}

async function warmExploreResults(jobId) {
  if (!jobId) {
    return;
  }
  if (exploreWarmupPromise && exploreWarmupJobId === jobId) {
    return exploreWarmupPromise;
  }
  const warmupPromise = (async () => {
    await loadDisplayFilters(jobId);
    if (loadedResultsJobId !== jobId) {
      loadedResultsJobId = jobId;
    }
    await refreshResults();
  })();
  exploreWarmupJobId = jobId;
  exploreWarmupPromise = warmupPromise.finally(() => {
    if (exploreWarmupJobId === jobId) {
      exploreWarmupJobId = null;
      exploreWarmupPromise = null;
    }
  });
  return exploreWarmupPromise;
}

function getDisplayFilterParams(panelKey) {
  const params = new URLSearchParams();
  const appendFilter = (index) => {
    const field = document.getElementById(panelElementId(panelKey, `filter${index}-field`)).value;
    const value = document.getElementById(panelElementId(panelKey, `filter${index}-values`)).value;
    if (!field || !value) {
      return;
    }
    params.append(`filter${index}_field`, field);
    params.append(`filter${index}_values`, value);
  };
  appendFilter(1);
  appendFilter(2);
  return params;
}

function getDisplayFilterSummary(panelKey) {
  const parts = [];
  const appendFilter = (index) => {
    const fieldSelect = document.getElementById(panelElementId(panelKey, `filter${index}-field`));
    const valueSelect = document.getElementById(panelElementId(panelKey, `filter${index}-values`));
    if (!fieldSelect || !valueSelect) {
      return;
    }
    const field = fieldSelect.value;
    const value = valueSelect.value;
    if (!field || !value) {
      return;
    }
    parts.push(`${field}: ${value}`);
  };
  appendFilter(1);
  appendFilter(2);
  return parts.join(" | ");
}

function updateBaselineFilterSummaries(panelKey) {
  const summary = getDisplayFilterSummary(panelKey);
  const text = summary ? `Display only: ${summary}` : "";
  setPanelSummary(panelKey, text);
}

function getPlotDotScale() {
  const value = Number(document.getElementById("plot-dot-scale")?.value || "0.5");
  return Number.isFinite(value) && value > 0 ? value * 4 : 1;
}

function setResultMode(mode) {
  const baseline = document.getElementById("baseline-results-view");
  const differential = document.getElementById("differential-results-view");
  if (currentJobStatus !== "completed" || !areExploreResultsReady()) {
    baseline.classList.add("hidden");
    differential.classList.add("hidden");
    setResultsControlsDisabled(true);
    return;
  }
  baseline.classList.remove("hidden");
  differential.classList.toggle("hidden", mode !== "differential");
  setResultsControlsDisabled(false);
}

function setResultsControlsDisabled(disabled) {
  VISUALIZATION_PANELS.forEach((panelKey) => {
    const controls = [
      document.getElementById(panelElementId(panelKey, "mode")),
      document.getElementById(panelElementId(panelKey, "modality")),
      document.getElementById(panelElementId(panelKey, "gene-query")),
      document.getElementById(panelElementId(panelKey, "marker-population")),
      document.getElementById(panelElementId(panelKey, "filter1-field")),
      document.getElementById(panelElementId(panelKey, "filter1-values")),
      document.getElementById(panelElementId(panelKey, "filter2-field")),
      document.getElementById(panelElementId(panelKey, "filter2-values")),
    ];
    controls.forEach((control) => {
      if (!control) {
        return;
      }
      control.disabled = disabled;
    });
  });
}

function isUmapMode(mode) {
  return mode === "relative" || mode === "cluster" || mode === "frequency";
}

function requiresGeneMode(mode) {
  return mode === "expression_umap" || mode === "violin";
}

function panelModeLabel(mode) {
  const labels = {
    relative: "umap_broad",
    cluster: "umap_cell_types",
    frequency: "cell_frequency",
    expression_umap: "expression_umap",
    violin: "violin",
    marker_heatmap: "marker_heatmap",
    marker_network: "marker_network",
    fastcomm_network: "fastcomm_network",
    grn_network: "grn_network",
  };
  return labels[mode] || "plot";
}

async function loadVisualizationPanel(panelKey) {
  const jobId = document.getElementById("results-job-id").value.trim();
  const mode = getPanelSelectValue(panelKey, "mode");
  const modality = panelModality(panelKey);
  if (!jobId || !mode) {
    return;
  }

  if(mode.startsWith("integrated_")) {
    const host=document.getElementById(panelElementId(panelKey,"plot"));
    resetVisualizationSurface(panelKey);
    setPanelSummary(panelKey, "");
    host.classList.remove("hidden");
    await ScalableIntegrated.mount(host,integratedOptions(jobId,mode,getPanelSelectValue(panelKey,"marker-population"),modality,"marker"));
    return;
  }

  // Gene-set plot types fetch their own payload and draw it here, so they do not
  // travel through the single-molecule expression path below.
  if (GENE_SET_MODES.has(mode)) {
    const genes = panelGeneSet(panelKey);
    const params = new URLSearchParams();
    // The DotPlot and CombPlot read the same modality the panel is set to. Without
    // this they always read RNA, so choosing ADT or Lipids left them unchanged.
    params.set("modality", modality);
    if (genes.length) {
      params.set("genes", genes.join(","));
    }
    if (mode === "combplot") {
      params.set("unit", panelCombUnit(panelKey));
      params.set("cells_per_sample", String(panelCellsPerSample(panelKey)));
      if (panelCombUnit(panelKey) === "donor") params.set("min_cells", String(panelCombMinCells(panelKey)));
    }
    // "Filter data to display" restricts these plots too, not only the UMAP and
    // violin. Without this the DotPlot and CombPlot ignored both annotation rows.
    appendGeneSetSubsetParams(params, panelKey);
    // Group by, and which levels to show. This call was lost in an edit, so
    // `group_by` stopped being sent and the plot silently stayed on cell state
    // however the control was set.
    appendGeneSetGroupParams(params, panelKey);
    const query = params.toString() ? `?${params.toString()}` : "";
    try {
      const response = await fetch(apiPath(`/api/jobs/${jobId}/${mode}${query}`));
      const payload = await response.json();
      if (!response.ok) throw new Error(payload.detail || `${mode} failed`);
      panelPlotData[panelKey] = { source: mode, payload };
    } catch (err) {
      panelPlotData[panelKey] = { source: "error", payload: { message: err.message } };
    }
    renderVisualizationPanel(panelKey);
    return;
  }
  try {
    if (isUmapMode(mode)) {
      const params = getDisplayFilterParams(panelKey);
      params.set("modality", modality);
      // Only the cell-type view offers these. "UMAP broad" draws reference
      // against query, and the reference is hidden whenever either choice is
      // off-default, which would leave that view with nothing to contrast.
      if (UMAP_COORD_MODES.has(mode)) {
        const colorBy = mode === "cluster" ? getPanelSelectValue(panelKey, "colorby") : "";
        const coordsKey = getPanelSelectValue(panelKey, "coords");
        if (colorBy) params.set("color_by", colorBy);
        // The X and Y fields are the coordinate control, so they always decide.
        const xField = getPanelSelectValue(panelKey, "xfield");
        const yField = getPanelSelectValue(panelKey, "yfield");
        if (xField && yField) {
          params.set("x_field", xField);
          params.set("y_field", yField);
        } else if (coordsKey) {
          params.set("coords", coordsKey);
        }
      }
      const suffix = params.toString() ? `?${params.toString()}` : "";
      const resp = await fetch(apiPath(`/jobs/${jobId}/umap${suffix}`));
      const data = await parseApiResponse(resp);
      if (!resp.ok) {
        throw new Error(data.detail || "UMAP not ready.");
      }
      panelPlotData[panelKey] = { source: "umap", payload: data };
      renderVisualizationPanel(panelKey);
      return;
    }

    if (mode === "marker_heatmap") {
      panelPlotData[panelKey] = { source: "marker_heatmap", payload: {} };
      renderVisualizationPanel(panelKey);
      return;
    }

    if (mode === "marker_network") {
      const population = getPanelSelectValue(panelKey, "marker-population");
      if (!population) {
        panelPlotData[panelKey] = {
          source: "marker_network",
          payload: { population: "", elements: [], message: "No marker network populations are available." },
        };
        renderVisualizationPanel(panelKey);
        return;
      }
      const resp = await fetch(apiPath(`/jobs/${jobId}/marker/network?population=${encodeURIComponent(population)}&modality=${encodeURIComponent(modality)}`));
      const data = await parseApiResponse(resp);
      if (!resp.ok) {
        throw new Error(data.detail || "Marker network unavailable.");
      }
      panelPlotData[panelKey] = { source: "marker_network", payload: data };
      renderVisualizationPanel(panelKey);
      return;
    }

    if (mode === "fastcomm_network") {
      const plotType = panelCommunicationPlotType(panelKey);
      const population = fastCommPlotNeedsPopulation(plotType) ? getPanelSelectValue(panelKey, "marker-population") : "";
      if (fastCommPlotNeedsPopulation(plotType) && !population) {
        panelPlotData[panelKey] = {
          source: "fastcomm_network",
          payload: { population: "", elements: [], message: "No fastComm cell states are available." },
        };
        renderVisualizationPanel(panelKey);
        return;
      }
      const params = getDisplayFilterParams(panelKey);
      params.set("population", population);
      params.set("plot_type", plotType);
      params.set("limit", "60");
      const resp = await fetch(apiPath(`/jobs/${jobId}/fastcomm/plot?${params.toString()}`));
      const data = await parseApiResponse(resp);
      if (!resp.ok) {
        throw new Error(data.detail || "fastComm network unavailable.");
      }
      panelPlotData[panelKey] = { source: "fastcomm_network", payload: data };
      renderVisualizationPanel(panelKey);
      return;
    }

    if (mode === "grn_network") {
      const genesRaw = String(document.getElementById(panelElementId(panelKey, "grn-genes"))?.value || "");
      const genes = genesRaw.split(/[\s,]+/).map((s) => s.trim()).filter(Boolean);
      const params = new URLSearchParams();
      genes.forEach((g) => params.append("genes", g));
      params.set("sample", getPanelSelectValue(panelKey, "grn-sample"));
      params.set("cell_state", getPanelSelectValue(panelKey, "grn-cellstate"));
      params.set("threshold", String(document.getElementById(panelElementId(panelKey, "grn-threshold"))?.value || "0"));
      params.set("max_edges", getPanelSelectValue(panelKey, "grn-limit") || "25");
      const resp = await fetch(apiPath(`/jobs/${jobId}/grn/network?${params.toString()}`));
      const data = await parseApiResponse(resp);
      if (!resp.ok) {
        throw new Error(data.detail || "GRN network unavailable.");
      }
      populateGrnDropdowns(panelKey, data);
      panelPlotData[panelKey] = { source: "grn_network", payload: data };
      renderVisualizationPanel(panelKey);
      return;
    }

    if (requiresGeneMode(mode)) {
      const geneInput = document.getElementById(panelElementId(panelKey, "gene-query"));
      const gene = String(geneInput?.value || "").trim();
      if (!gene) {
        panelPlotData[panelKey] = {
          source: "missing_gene",
          payload: { message: "Select a molecule to render this plot.", umap: [], violin: [], gene: "", requested_gene: "" },
        };
        renderVisualizationPanel(panelKey);
        return;
      }
      const params = getDisplayFilterParams(panelKey);
      params.set("gene", gene);
      params.set("modality", modality);
        // One window is twice as wide, so the violin plot draws more cell states
        // rather than leaving the extra space empty.
        params.set("violin_limit", singleWindowActive() ? "30" : "10");
        // The expression UMAP draws on the same coordinate pair as the other UMAP
        // panels. Without this the control appeared but changed nothing.
        if (UMAP_COORD_MODES.has(mode)) {
          const xField = getPanelSelectValue(panelKey, "xfield");
          const yField = getPanelSelectValue(panelKey, "yfield");
          if (xField && yField) {
            params.set("x_field", xField);
            params.set("y_field", yField);
          }
        }
      const resp = await fetch(apiPath(`/jobs/${jobId}/expression?${params.toString()}`));
      const data = await parseApiResponse(resp);
      if (!resp.ok) {
        throw new Error(data.detail || "Expression unavailable.");
      }
      if (data?.gene && data.gene !== gene && geneInput) {
        geneInput.value = data.gene;
      }
      panelPlotData[panelKey] = { source: "expression", payload: data };
      renderVisualizationPanel(panelKey);
      return;
    }
  } catch (err) {
    panelPlotData[panelKey] = {
      source: "error",
      payload: { message: err.message || "Plot unavailable." },
    };
    renderVisualizationPanel(panelKey);
  }
}

function renderVisualizationMessage(panelKey, message, title = "Visualization unavailable") {
  resetVisualizationSurface(panelKey);
  setPanelSummary(panelKey, "");
  Plotly.newPlot(panelPlotId(panelKey), [], {
    title,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.9)",
    height: 560,
    margin: { t: 48, l: 32, r: 32, b: 32 },
    xaxis: { visible: false },
    yaxis: { visible: false },
    annotations: [
      {
        text: message,
        x: 0.5,
        y: 0.5,
        yshift: -3,
        xref: "paper",
        yref: "paper",
        showarrow: false,
        font: { size: 15, color: "#64748b" },
        align: "center",
      },
    ],
  }, { displayModeBar: false });
}

/* The caption under a UMAP panel. It names the colour column and the embedding
 * whenever they are not the cellHarmony defaults, and says why the reference
 * atlas is missing, which otherwise reads as a drawing fault. */
function setUmapPanelSummary(panelKey, umapData) {
  const parts = [];
  const filters = getDisplayFilterSummary(panelKey);
  if (filters) parts.push(`Display only: ${filters}`);
  if (umapData && umapData.reference_hidden) {
    parts.push(`Colored by ${umapData.color_label}`);
    parts.push(umapData.axes_source === "obs"
      ? `Axes: ${umapData.x_label} against ${umapData.y_label}`
      : `Coordinates: ${umapData.coords_label}`);
    parts.push("Reference atlas hidden: it carries neither this annotation nor these coordinates");
  }
  // Cells the panel could not place. Silence here would read as an absence of
  // cells rather than an absence of a recorded value.
  if (umapData && Number(umapData.n_dropped_no_coordinate) > 0) {
    parts.push(`${umapData.n_dropped_no_coordinate} of ${umapData.n_cells_selected} cells not drawn: no value on both axes`);
  }
  if (parts.length) setPanelSummary(panelKey, parts.join(" | "));
}

function renderPanelUmap(panelKey, umapData, mode, dotScale) {
  updateBaselineFilterSummaries(panelKey);
  const traces = [];
  const layout = {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.9)",
    height: 560,
    margin: { t: 20, l: 40, r: 20, b: 40 },
    hovermode: "closest",
    xaxis: { showgrid: false, zeroline: false },
    yaxis: { showgrid: false, zeroline: false },
  };
  if (mode === "frequency") {
    const sampleField = String(umapData.sample_field || "sample").trim() || "sample";
    const queryPoints = (umapData.query || []).filter((point) => String(point.sample || "").trim());
    if (!queryPoints.length) {
      Plotly.newPlot(panelPlotId(panelKey), [], {
        ...layout,
        title: "Cell frequency",
        xaxis: { visible: false },
        yaxis: { visible: false },
        annotations: [
          {
            text: "Sample labels were not available for this job.",
            x: 0.5,
            y: 0.5,
            xref: "paper",
            yref: "paper",
            showarrow: false,
            font: { size: 15, color: "#64748b" },
          },
        ],
      });
      return;
    }

    const countsBySample = new Map();
    queryPoints.forEach((point) => {
      const sample = String(point.sample || "").trim();
      const population = String(point.population || "").trim();
      if (!sample || !population) {
        return;
      }
      let sampleMap = countsBySample.get(sample);
      if (!sampleMap) {
        sampleMap = new Map();
        countsBySample.set(sample, sampleMap);
      }
      sampleMap.set(population, (sampleMap.get(population) || 0) + 1);
    });

    const populationStats = new Map();
    countsBySample.forEach((sampleMap, sample) => {
      let total = 0;
      sampleMap.forEach((count) => {
        total += count;
      });
      if (!(total > 0)) {
        return;
      }
      sampleMap.forEach((count, population) => {
        const fraction = count / total;
        let entry = populationStats.get(population);
        if (!entry) {
          entry = { totalFraction: 0, sampleFractions: [] };
          populationStats.set(population, entry);
        }
        entry.totalFraction += fraction;
        entry.sampleFractions.push({ sample, fraction, count, total });
      });
    });

    const ranked = [...populationStats.entries()]
      .map(([population, entry]) => ({
        population,
        meanFraction: entry.sampleFractions.length ? entry.totalFraction / entry.sampleFractions.length : 0,
        sampleFractions: entry.sampleFractions,
      }))
      .sort((a, b) => b.meanFraction - a.meanFraction || a.population.localeCompare(b.population));
    const maxLabelLength = ranked.reduce((maxLength, entry) => Math.max(maxLength, String(entry.population || "").length), 0);
    const leftMargin = Math.max(70, Math.min(120, 36 + maxLabelLength * 5));

    Plotly.newPlot(panelPlotId(panelKey), [
      {
        type: "bar",
        orientation: "h",
        y: ranked.map((entry) => entry.population),
        x: ranked.map((entry) => entry.meanFraction),
        customdata: ranked.map((entry) => [
          entry.sampleFractions
            .sort((a, b) => b.fraction - a.fraction || a.sample.localeCompare(b.sample))
            .map((item) => `${item.sample}: ${(item.fraction * 100).toFixed(1)}% (${item.count}/${item.total})`)
            .join("<br>"),
        ]),
        marker: {
          color: ranked.map((entry) => interpolatePairedColor(Math.min(1, Math.max(0, entry.meanFraction)))),
        },
        hovertemplate: "%{y}<br>Mean normalized frequency=%{x:.3f}<br>%{customdata[0]}<extra></extra>",
      },
    ], {
      ...layout,
      title: "Cell frequency",
      height: Math.max(320, ranked.length * 17 + 90),
      margin: { t: 36, l: leftMargin, r: 20, b: 40 },
      xaxis: {
        title: `Mean fraction of filtered cells per ${sampleField}`,
        range: [0, 1],
        tickformat: ".0%",
        showgrid: false,
        zeroline: false,
      },
      yaxis: {
        tickmode: "array",
        tickvals: ranked.map((entry) => entry.population),
        ticktext: ranked.map((entry) => entry.population),
        tickfont: { size: 9 },
        automargin: false,
        autorange: "reversed",
        showgrid: false,
        zeroline: false,
      },
    });
    return;
  }
  if (mode === "relative") {
    traces.push({
      x: umapData.reference.map((p) => p.x),
      y: umapData.reference.map((p) => p.y),
      text: umapData.reference.map((p) => `${p.barcode}<br>${p.population}`),
      mode: "markers",
      type: "scattergl",
      marker: { color: "#94a3b8", size: 2 * dotScale },
      name: "Reference",
    });
    traces.push({
      x: umapData.query.map((p) => p.x),
      y: umapData.query.map((p) => p.y),
      text: umapData.query.map((p) => `${p.barcode}<br>${p.population}`),
      mode: "markers",
      type: "scattergl",
      marker: { color: "#f97316", size: 2 * dotScale },
      name: "Query",
    });
    layout.legend = { orientation: "h" };
    Plotly.newPlot(panelPlotId(panelKey), traces, layout);
    return;
  }
  setUmapPanelSummary(panelKey, umapData);
  const populations = buildStableUmapPopulationOrder(umapData);
  const colorMap = buildReferencePreviewColorMap(populations);
  const labelPoints = relaxReferencePreviewLabels(
    buildPopulationCentroids(umapData.query),
    umapData.query,
    document.getElementById(panelPlotId(panelKey))
  );
  traces.push({
    x: umapData.reference.map((p) => p.x),
    y: umapData.reference.map((p) => p.y),
    text: umapData.reference.map((p) => `${p.barcode}<br>${p.population}`),
    mode: "markers",
    type: "scattergl",
    marker: { color: "#e5e7eb", size: Math.max(0.5, 1 * dotScale), opacity: 0.3 },
    name: "Reference",
    showlegend: false,
  });
  traces.push({
    x: umapData.query.map((p) => p.x),
    y: umapData.query.map((p) => p.y),
    text: umapData.query.map((p) => `${p.barcode}<br>${p.population}`),
    mode: "markers",
    type: "scattergl",
    marker: {
      size: 2 * dotScale,
      opacity: 0.5,
      color: umapData.query.map((p) => colorMap.get(p.population)),
    },
    showlegend: false,
    name: "Query",
  });
  layout.annotations = labelPoints.map((label) => ({
    x: label.x,
    y: label.y,
    text: label.population,
    showarrow: false,
    xref: "x",
    yref: "y",
    xanchor: "center",
    yanchor: "middle",
    font: { size: 11, color: "#0f172a" },
    bgcolor: "rgba(255,255,255,0)",
    opacity: 1,
  }));
  if (umapData.axes_source === "obs") {
    // A UMAP hides its tick labels because the numbers mean nothing. A metadata
    // axis is the opposite: the value is the point of the plot, so it is
    // labelled and its ticks are shown.
    layout.xaxis = { title: { text: umapData.x_label, font: { size: 12 } }, showgrid: false, zeroline: false };
    layout.yaxis = { title: { text: umapData.y_label, font: { size: 12 } }, showgrid: false, zeroline: false };
    layout.margin = { t: 20, l: 60, r: 20, b: 55 };
  } else {
    Object.assign(layout, buildSquareUmapAxes(umapData.query, 0.06));
  }
  Plotly.newPlot(panelPlotId(panelKey), traces, layout);
}

function renderPanelExpression(panelKey, expressionData, mode, dotScale) {
  updateBaselineFilterSummaries(panelKey);
  const displayGene = featureDisplayName(expressionData?.gene, expressionData?.modality);
  const valueLabel = measurementLabel(expressionData?.modality);
  if (!expressionData?.umap?.length && !expressionData?.violin?.length) {
    renderVisualizationMessage(
      panelKey,
      expressionData?.message || `Gene '${expressionData?.requested_gene || expressionData?.gene || ""}' was not found.`,
      `${expressionData?.requested_gene || expressionData?.gene || "gene"} expression`,
    );
    return;
  }
  if (mode === "violin") {
    resetVisualizationSurface(panelKey);
    if (expressionData.source === "reference") {
      const trace = {
        type: "bar",
        x: expressionData.violin.map((entry) => entry.population),
        y: expressionData.violin.map((entry) => entry.mean),
        marker: { color: "#475569", opacity: 0.85 },
      };
      Plotly.newPlot(panelPlotId(panelKey), [trace], {
        title: `${displayGene} (reference centroid expression, top 10 states)`,
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(255,255,255,0.9)",
        height: 560,
        margin: { t: 36, l: 48, r: 20, b: 120 },
      });
      return;
    }
    const globalMin = Number(expressionData.global_min ?? 0);
    let globalMax = Number(expressionData.global_max ?? 0);
    if (!(Number.isFinite(globalMax) && Number.isFinite(globalMin))) {
      globalMax = 0;
    }
    if (!(globalMax > globalMin)) {
      globalMax = globalMin + 1e-9;
    }
    const yPad = Math.max((globalMax - globalMin) * 0.04, 0.05);
    const violinTraces = expressionData.violin.map((entry) => ({
      type: "violin",
      hovertemplate: `${featureHover(expressionData.gene, expressionData.modality)}<br>%{y:.3f}<extra>%{fullData.name}</extra>`,
      name: entry.population,
      y: entry.values,
      box: { visible: false },
      meanline: { visible: true },
      points: "all",
      jitter: 0.3,
      pointpos: 0,
      marker: { size: 3 * dotScale, opacity: 0.55 },
    }));
    Plotly.newPlot(panelPlotId(panelKey), violinTraces, {
      title: `${displayGene} (top 10 states by mean)`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.9)",
      height: 560,
      margin: { t: 36, l: 48, r: 20, b: 120 },
      yaxis: { title: valueLabel, range: [globalMin - yPad, globalMax + yPad] },
    });
    return;
  }
  const zeroPoints = expressionData.umap.filter((p) => Number(p.value) === 0);
  const expressedPoints = expressionData.umap
    .filter((p) => Number.isFinite(Number(p.value)) && Number(p.value) !== 0)
    .sort((a, b) => Number(a.value) - Number(b.value));
  const traces = [];
  if (zeroPoints.length) {
    traces.push({
      x: zeroPoints.map((p) => p.x),
      y: zeroPoints.map((p) => p.y),
      text: zeroPoints.map((p) => `${p.barcode}<br>${p.population}<br>${featureHover(expressionData.gene, expressionData.modality)}: 0.000`),
      mode: "markers",
      type: "scattergl",
      marker: {
        size: 2 * dotScale,
        color: "#e5e7eb",
        opacity: 0.9,
      },
      hoverinfo: "text",
      name: `${displayGene} = 0`,
      showlegend: false,
    });
  }
  if (expressedPoints.length) {
    let minValue = Number(expressionData.global_min ?? expressedPoints[0].value);
    let maxValue = Number(expressionData.global_max ?? expressedPoints[expressedPoints.length - 1].value);
    if (!(Number.isFinite(minValue))) {
      minValue = Number(expressedPoints[0].value);
    }
    if (!(Number.isFinite(maxValue))) {
      maxValue = Number(expressedPoints[expressedPoints.length - 1].value);
    }
    if (!(maxValue > minValue)) {
      maxValue = minValue + 1e-9;
    }
    const normalizedModality = normalizeModalityId(expressionData?.modality);
    const useImputedPalette = normalizedModality === "lipids" || normalizedModality === "adt" || normalizedModality === "metabolite" || normalizedModality === "lipid" || normalizedModality === "grn" || normalizedModality === "grn_tf";
    const colorscale = useImputedPalette
      ? [
          [0.0, "#2563eb"],
          [0.5, "#fde047"],
          [1.0, "#dc2626"],
        ]
      : [
          [0.0, "#f3f4f6"],
          [0.15, "#fecaca"],
          [0.35, "#fca5a5"],
          [0.6, "#ef4444"],
          [1.0, "#b91c1c"],
        ];
    const binCount = Math.min(12, Math.max(4, Math.ceil(Math.sqrt(expressedPoints.length / 2500))));
    const step = Math.max(1, Math.ceil(expressedPoints.length / binCount));
    for (let start = 0; start < expressedPoints.length; start += step) {
      const bin = expressedPoints.slice(start, Math.min(start + step, expressedPoints.length));
      traces.push({
        x: bin.map((p) => p.x),
        y: bin.map((p) => p.y),
        text: bin.map((p) => `${p.barcode}<br>${p.population}<br>${featureHover(expressionData.gene, expressionData.modality)}: ${p.value.toFixed(3)}`),
        mode: "markers",
        type: "scattergl",
        marker: {
          size: 2 * dotScale,
          color: bin.map((p) => p.value),
          colorscale,
          cmin: minValue,
          cmax: maxValue,
          showscale: start + step >= expressedPoints.length,
          colorbar: start + step >= expressedPoints.length ? { title: displayGene } : undefined,
        },
        name: displayGene,
        showlegend: false,
      });
    }
  }
  const expressionLayout = {
    title: expressionData.source === "reference"
      ? `${displayGene} reference expression`
      : `${displayGene} ${valueLabel.charAt(0).toLowerCase() + valueLabel.slice(1)}`,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.9)",
    height: 560,
    margin: { t: 36, l: 48, r: 20, b: 40 },
    xaxis: { showgrid: false, zeroline: false },
    yaxis: { showgrid: false, zeroline: false },
    annotations: [],
  };
  if (expressionData.message) {
    expressionLayout.annotations.push({
      text: expressionData.message,
      x: 0.5,
      y: 1.08,
      xref: "paper",
      yref: "paper",
      showarrow: false,
      font: { size: 12, color: "#64748b" },
    });
  }
  Object.assign(expressionLayout, buildSquareUmapAxes(expressionData.umap, 0.06));
  resetVisualizationSurface(panelKey);
  Plotly.newPlot(panelPlotId(panelKey), traces, expressionLayout);
}

function renderVisualizationPanel(panelKey) {
  const mode = getPanelSelectValue(panelKey, "mode");
  // Integrated views own their live SVG/Cytoscape state. Repainting an old
  // expression payload here replaces a pathway when the window count changes.
  if (mode.startsWith("integrated_")) {
    document.getElementById(panelPlotId(panelKey))?._integratedResize?.();
    return;
  }
  const data = panelPlotData[panelKey];

  // Gene-set figures are drawn from their own payload, not the expression one.
  if (GENE_SET_MODES.has(mode)) {
    if (!data || data.source === "error") {
      renderVisualizationMessage(panelKey,
        (data && data.payload && data.payload.message) || "Enter a gene set to draw this plot.");
      return;
    }
    resetVisualizationSurface(panelKey);
    setPanelSummary(panelKey, "");
    if (mode === "combplot") {
      renderCombPlotFigure(panelPlotId(panelKey), data.payload || {});
      const sampling = data.payload?.sampling;
      setPanelSummary(panelKey, sampling ? sampling.description : "");
    }
    else renderDotPlotFigure(panelPlotId(panelKey), data.payload || {});
    return;
  }
  const dotScale = getPlotDotScale();
  if (!data) {
    renderVisualizationMessage(panelKey, "Complete an alignment to visualize results.");
    return;
  }
  if (data.source === "error") {
    renderVisualizationMessage(panelKey, data.payload?.message || "Plot unavailable.");
    return;
  }
  if (mode === "marker_heatmap") {
    const jobId = document.getElementById("results-job-id").value.trim();
    const featureLabel = String(modalityDefinition(panelModality(panelKey)).feature_label || "gene");
    setPanelSummary(panelKey, `Marker heatmap view uses the exported ${featureLabel} marker analysis matrix.`);
    if (!jobId || !markerHeatmapAvailable(panelKey)) {
      renderVisualizationMessage(panelKey, "Marker heatmap output is unavailable for this job.", "MarkerHeatmap");
      return;
    }
    renderMarkerHeatmapViewer(jobId, panelKey);
    return;
  }
  if (mode === "marker_network") {
    const population = String(data.payload?.population || getPanelSelectValue(panelKey, "marker-population") || "").trim();
    setPanelSummary(panelKey, population ? `Marker network view for ${population}.` : "Marker network view.");
    if (!(data.payload?.elements || []).length) {
      renderVisualizationMessage(panelKey, data.payload?.message || "Marker network unavailable.", population ? `${population} marker network` : "MarkerNetwork");
      return;
    }
    renderExpressionNetwork(panelKey, data.payload);
    return;
  }
  if (mode === "grn_network") {
    const p = data.payload || {};
    const description = document.getElementById(panelElementId(panelKey, "grn-description"));
    if (description) description.textContent = grnNetworkDescription(p);
    setPanelSummary(panelKey, "Yellow diamond = TF; red circle = target; dark outline = your selected gene. Colors identify node roles. Arrow = TF → target. Width scales with |mean score| within this graph; hover for scores. These predictions do not indicate activation/repression or marker up/down regulation.");
    if (!(p.elements || []).length) {
      renderVisualizationMessage(panelKey, p.message || "No GRN edges matched.", "GRN edges");
      return;
    }
    renderExpressionNetwork(panelKey, p);
    return;
  }
  if (mode === "fastcomm_network") {
    const population = String(data.payload?.population || getPanelSelectValue(panelKey, "marker-population") || "").trim();
    const plotType = data.payload?.plot_type || panelCommunicationPlotType(panelKey);
    const plotLabel = (FASTCOMM_PLOT_OPTIONS.find((entry) => entry.id === plotType) || {}).label || "Cell communication";
    setPanelSummary(panelKey, population ? `${plotLabel} for ${population}.` : plotLabel);
    const hasNetwork = (data.payload?.elements || []).length > 0;
    const hasRows = (data.payload?.rows || []).length > 0;
    const hasPoints = (data.payload?.points || []).length > 0;
    const hasHeatmap = (data.payload?.z || []).length > 0;
    if (!hasNetwork && !hasRows && !hasPoints && !hasHeatmap) {
      renderVisualizationMessage(panelKey, data.payload?.message || "Cell communication plot unavailable.", population ? `${population} cell communication` : "Cell communication");
      return;
    }
    renderFastCommPlot(panelKey, data.payload);
    return;
  }
  if (isUmapMode(mode)) {
    renderPanelUmap(panelKey, data.payload || {}, mode, dotScale);
    return;
  }
  renderPanelExpression(panelKey, data.payload || {}, mode, dotScale);
}

async function refreshResults() {
  await Promise.all(VISUALIZATION_PANELS.map((panelKey) => loadVisualizationPanel(panelKey)));
}

async function downloadVisualizationImage(panelKey) {
  const jobId = document.getElementById("results-job-id").value.trim();
  const mode = getPanelSelectValue(panelKey, "mode");
  const modality = panelModality(panelKey);
  if (jobId && mode?.startsWith("integrated_")) {
    try {await document.getElementById(panelPlotId(panelKey))._integratedPdf?.();}
    catch(error){showDownloadError(error);}return;
  }
  if (!jobId || !mode || !panelPlotData[panelKey]) {
    return;
  }
  if (mode === "marker_heatmap") {
    const params = getDisplayFilterParams(panelKey);
    params.set("modality", modality);
    params.set("cells_per_sample", String(panelCellsPerSample(panelKey)));
    window.open(apiPath(`/jobs/${jobId}/marker/heatmap.pdf?${params.toString()}`), "_blank");
    return;
  }
  if (mode === "marker_network") {
    const population = getPanelSelectValue(panelKey, "marker-population");
    if (!population) {
      return;
    }
    window.open(apiPath(`/jobs/${jobId}/marker/network/pdf?population=${encodeURIComponent(population)}&modality=${encodeURIComponent(modality)}`), "_blank");
    return;
  }
  if (mode === "fastcomm_network") {
    const cy = expressionCyByPanel[panelKey];
    const plotType = panelCommunicationPlotType(panelKey);
    if (cy) {
      try {
        await ensureCytoscapeSvgLoaded();
        await saveSvgMarkupAsPdf(cy.svg({ scale: 1, full: true }),
          buildPdfFilename([jobId, panelKey, plotType], plotType));
      } catch (err) {
        showDownloadError(err);
      }
      return;
    }
    await exportPlotlyElementToVectorPdf(
      panelPlotId(panelKey),
      buildPdfFilename([jobId, panelKey, plotType], plotType),
    );
    return;
  }
  if (mode === "grn_network") {
    const cy = expressionCyByPanel[panelKey];
    if (!cy) return;
    try {
      await ensureCytoscapeSvgLoaded();
      await saveSvgMarkupAsPdf(cy.svg({ scale: 1, full: true }),
        buildPdfFilename([jobId, panelKey, getPanelSelectValue(panelKey, "grn-cellstate"), "grn_edges"], "grn_edges"),
        grnNetworkDescription(panelPlotData[panelKey]?.payload || {}) + " "
          + document.getElementById(panelElementId(panelKey, "filter-summary")).textContent);
    } catch (error) {
      showDownloadError(error);
    }
    return;
  }
  try {
    const gene = getPanelSelectValue(panelKey, "gene-query");
    await exportPlotlyElementToVectorPdf(
      panelPlotId(panelKey),
      buildPdfFilename([jobId, panelKey, panelModeLabel(mode), gene || "gene"], panelModeLabel(mode)),
    );
  } catch (error) {
    console.error(error);
    showDownloadError(error);
  }
}

async function downloadDifferentialLeftPdf() {
  const jobId = document.getElementById("results-job-id").value.trim();
  const state = currentDifferentialState;
  if (!jobId || !state || state.status !== "completed") {
    return;
  }
  const mode = document.getElementById("differential-viz-mode").value;
  const population = document.getElementById("differential-result-population").value;
  if (mode === "summary") {
    window.open(apiPath(`/jobs/${jobId}/differential/interactive/pdf?mode=summary`), "_blank");
    return;
  }
  if (!population) {
    return;
  }
  if (mode.startsWith("integrated_")) {
    try {await document.getElementById("differential-plot-area")._integratedPdf?.();}
    catch(error){showDownloadError(error);}return;
  }
  if (mode === "table") {
    const plot = document.getElementById("differential-plot-area");
    if (plot && plot.data && plot.layout) {
      try {
        await exportPlotlyElementToVectorPdf(
          plot,
          buildPdfFilename([jobId, population, "differential_interaction_table"], "differential_interaction_table"),
        );
        return;
      } catch (error) {
        console.error(error);
        showDownloadError(error);
        return;
      }
    }
  }
  // The server PDF renderer reads the whole differential result and knows nothing about
  // the gene filter. With a filter on, export what the screen shows instead, so the PDF
  // and the view always carry the same genes.
  if (String(currentDifferentialGeneFilter || "").trim()) {
    const filterTag = `filtered_${String(currentDifferentialGeneFilter).trim()}`;
    try {
      if (mode === "network") {
        await exportCytoscapeToVectorPdf(
          differentialCy,
          buildPdfFilename([jobId, population, mode, filterTag], mode),
        );
      } else {
        await exportPlotlyElementToVectorPdf(
          document.getElementById("differential-plot-area"),
          buildPdfFilename([jobId, population, mode, filterTag], mode),
        );
      }
      return;
    } catch (error) {
      console.error(error);
      showDownloadError(error);
      return;
    }
  }
  window.open(
    apiPath(`/jobs/${jobId}/differential/interactive/pdf?mode=${encodeURIComponent(mode)}&population=${encodeURIComponent(population)}`),
    "_blank",
  );
}

async function downloadDifferentialGenePdf() {
  const jobId = document.getElementById("results-job-id").value.trim();
  if (!jobId || !currentDifferentialGene) {
    return;
  }
  const plot = document.getElementById("differential-gene-plot");
  if (plot && plot.data && plot.layout) {
    try {
      await exportPlotlyElementToVectorPdf(
        plot,
        buildPdfFilename([jobId, currentDifferentialPopulation || "population", currentDifferentialGene, "gene_detail"], "gene_detail"),
      );
      return;
    } catch (error) {
      console.error(error);
      showDownloadError(error);
      return;
    }
  }
  let pdfUrl = `/jobs/${jobId}/differential/interactive/gene/pdf?population=${encodeURIComponent(currentDifferentialPopulation)}&gene=${encodeURIComponent(currentDifferentialGene)}`;
  if (currentDifferentialInteraction) {
    const role = currentDifferentialFeatureRole;
    const featureSymbol = role === "receptor"
      ? currentDifferentialInteraction.receptor
      : currentDifferentialInteraction.ligand;
    if (featureSymbol) {
      pdfUrl += `&feature=${encodeURIComponent(featureSymbol)}`;
    }
  }
  window.open(apiPath(pdfUrl), "_blank");
}

let activeExplorerTab = "run";

function initExplorerSandbox() {
  const buttons = Array.from(document.querySelectorAll(".workspace-tab-btn"));
  buttons.forEach((button) => {
    button.addEventListener("click", () => {
      const tab = String(button.dataset.tab || "run");
      setExplorerTab(tab);
    });
  });
  setExplorerTab("run");
}

function setExplorerTab(tab) {
  if (tab === "explore" && !areExploreResultsReady()) {
    tab = "run";
  }
  if (tab === "differential" && currentJobStatus !== "completed") {
    tab = "run";
  }
  const tabChanged = activeExplorerTab !== tab;
  activeExplorerTab = tab;
  document.querySelectorAll(".workspace-tab-btn").forEach((button) => {
    button.classList.toggle("active", button.dataset.tab === tab);
  });
  document.querySelectorAll(".workspace-panel").forEach((panel) => {
    panel.classList.toggle("active", panel.dataset.tabPanel === tab);
  });
  if (tabChanged && tab === "explore") {
    // Initial plots are prepared while the workspace is hidden. Resize only
    // after it is visible, so restored plots fit the actual panel dimensions.
    requestAnimationFrame(() => {
      VISUALIZATION_PANELS.forEach((panelKey) => {
        const plot = document.getElementById(panelPlotId(panelKey));
        if (plot?.data?.length && plot.offsetWidth && plot.offsetHeight) {
          Plotly.Plots.resize(plot);
        }
      });
    });
  }
}

function syncExplorerWorkspace(preferredTab = null) {
  const baseline = document.getElementById("baseline-results-view");
  const differential = document.getElementById("differential-results-view");
  const hasBaseline = Boolean(baseline && !baseline.classList.contains("hidden"));
  const hasDifferential = Boolean(differential && !differential.classList.contains("hidden"));

  if (preferredTab === "run") {
    setExplorerTab("run");
    return;
  }

  if (preferredTab === "differential") {
    setExplorerTab("differential");
    return;
  }
  if (preferredTab === "explore" && hasBaseline) {
    setExplorerTab("explore");
    return;
  }
  if (!hasBaseline && !hasDifferential && activeExplorerTab === "explore") {
      setExplorerTab("run");
  }
}

/* ---------------------------------------------------------------------------
 * Window count, gene-set plot types and the Chat tab.
 *
 * These are scALABLE features, not viewer overrides: the upload workflow and any
 * precomputed bundle get the same controls from this file.
 * ------------------------------------------------------------------------ */

/* Windows = 2 | 1. The results grid is two equal columns, so a plot that runs
 * across every cell state has half the page. Choosing 1 hides panel 2 and lets
 * panel 1 span the full width. Panel 2 is hidden, never reset, so its plot and
 * its settings are unchanged when the count goes back to 2. */
/* Whether the results grid is currently showing a single window. */
function singleWindowActive() {
  const view = document.getElementById("baseline-results-view");
  return Boolean(view && view.classList.contains("single-window"));
}

function applyWindowCount(count) {
  const view = document.getElementById("baseline-results-view");
  if (!view) return;
  view.classList.toggle("single-window", count === 1);
  // A Plotly figure is drawn to the width the panel had, so both panels are
  // redrawn at the new one.
  if (typeof renderVisualizationPanel === "function") {
    renderVisualizationPanel("viz1");
    if (count === 2) renderVisualizationPanel("viz2");
  }
}

function initWindowCount() {
  const select = document.getElementById("viz-window-count");
  if (!select || select.dataset.wired === "1") return;
  select.dataset.wired = "1";
  select.addEventListener("change", () => {
    applyWindowCount(Number(select.value) === 1 ? 1 : 2);
  });
}

/* Gene symbols out of whatever the user pasted: an Excel column arrives newline
 * separated, a row tab separated, and people also type commas and spaces.
 * Duplicates are dropped and order is kept, so the plot reads in the order the
 * genes were listed. */
function parseGeneSet(text) {
  const seen = new Set();
  const out = [];
  String(text || "").split(/[\s,;|]+/).forEach((part) => {
    const gene = part.trim().replace(/^["']|["']$/g, "");
    if (gene && !seen.has(gene)) { seen.add(gene); out.push(gene); }
  });
  return out;
}

function panelGeneSet(panelKey) {
  const box = document.getElementById(panelElementId(panelKey, "geneset"));
  return parseGeneSet(box ? box.value : "");
}

function panelCombUnit(panelKey) {
  return document.getElementById(panelElementId(panelKey, "combunit"))?.value === "donor" ? "donor" : "cells";
}

function panelCombMinCells(panelKey) {
  const select = document.getElementById(panelElementId(panelKey, "combmin"));
  const value = Number(select && select.value);
  return Number.isFinite(value) && value > 0 ? value : 5;
}

function initGeneSetBoxes() {
  VISUALIZATION_PANELS.forEach((panelKey) => {
    const box = document.getElementById(panelElementId(panelKey, "geneset"));
    const min = document.getElementById(panelElementId(panelKey, "combmin"));
    const unit = document.getElementById(panelElementId(panelKey, "combunit"));
    if (unit && unit.dataset.wired !== "1") {
      unit.dataset.wired = "1";
      unit.addEventListener("change", () => {
        updateExpressionModeOptions();
        loadVisualizationPanel(panelKey);
      });
    }
    if (box && box.dataset.wired !== "1") {
      box.dataset.wired = "1";
      let timer = null;
      box.addEventListener("input", () => {
        window.clearTimeout(timer);
        timer = window.setTimeout(() => loadVisualizationPanel(panelKey), 600);
      });
    }
    if (min && min.dataset.wired !== "1") {
      min.dataset.wired = "1";
      min.addEventListener("change", () => loadVisualizationPanel(panelKey));
    }
  });
}

/* The Chat tab.
 *
 * The question goes to /api/jobs/{id}/chat. That route asks the LungMAP
 * assistant to read the sentence into one of a few supported queries, then runs
 * the query here and returns real numbers. The model never sees the data, so
 * nothing on this panel is a generated statistic.
 */
let chatLastResult = null;

function initChatTab() {
  const send = document.getElementById("chat-send");
  if (!send || send.dataset.wired === "1") return;
  send.dataset.wired = "1";
  send.addEventListener("click", askChat);

  const box = document.getElementById("chat-question");
  if (box) {
    box.addEventListener("keydown", (event) => {
      if (event.key === "Enter" && !event.shiftKey) { event.preventDefault(); askChat(); }
    });
  }
  document.querySelectorAll(".chat-example").forEach((example) => {
    example.addEventListener("click", () => {
      document.getElementById("chat-question").value = example.textContent;
      askChat();
    });
  });
  const views = document.getElementById("chat-views");
  if (views) {
    views.querySelectorAll("button").forEach((view) => {
      view.addEventListener("click", () => {
        views.querySelectorAll("button").forEach((b) => b.classList.remove("active"));
        view.classList.add("active");
        const wantsPlot = view.dataset.view === "plot";
        document.getElementById("chat-plot").classList.toggle("hidden", !wantsPlot);
        document.getElementById("chat-table").classList.toggle("hidden", wantsPlot);
        if (wantsPlot) drawChatPlot();
      });
    });
  }
}

/* The Chat examples, built by the server from this job's reference and its own
 * cell states. The tab shipped eight fixed lung sentences, so a bone-marrow job
 * offered AT2 and COPD questions that its data cannot answer. */
async function loadChatExamples(jobId) {
  const host = document.getElementById("chat-examples");
  const box = document.getElementById("chat-question");
  if (!host || !jobId) return;
  try {
    const response = await fetch(apiPath(`/jobs/${jobId}/chat-examples`));
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail || "examples unavailable");
    host.innerHTML = "";
    const label = document.createElement("span");
    label.textContent = data.reference ? `Try (${data.reference}):` : "Try:";
    host.appendChild(label);
    (data.examples || []).forEach((question) => {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "ghost-btn chat-example";
      button.textContent = question;
      button.addEventListener("click", () => {
        box.value = question;
        askChat();
      });
      host.appendChild(button);
    });
    if (box && data.placeholder) box.placeholder = data.placeholder;
  } catch (err) {
    console.warn("chat examples unavailable:", err);
  }
}

/* Empty every panel the previous answer wrote.
 *
 * The old code cleared the answer text and the table but left the figure and
 * the follow-up buttons standing. A new question therefore showed the previous
 * question's plot for as long as the server took to reply, and the two read as
 * one answer. Plotly holds its own state on the node, so purging it is what
 * actually removes the figure; emptying `innerHTML` alone leaves the chart
 * registered and the next `newPlot` inherits its layout. */
function clearChatOutput() {
  const correlationPlot = document.getElementById("chat-correlation-plot");
  if (correlationPlot) Plotly.purge(correlationPlot);
  chatResultState = null;
  document.getElementById("chat-result-controls")?.classList.add("hidden");
  chatLastResult = null;
  try { Plotly.purge("chat-plot"); } catch (err) { /* nothing drawn yet */ }
  ["chat-answer", "chat-table", "chat-plot", "chat-followups"].forEach((id) => {
    const host = document.getElementById(id);
    if (host) host.innerHTML = "";
  });
  const plot = document.getElementById("chat-plot");
  if (plot) plot.classList.add("hidden");
  const table = document.getElementById("chat-table");
  if (table) table.classList.add("hidden");
  const views = document.getElementById("chat-views");
  if (views) views.classList.add("hidden");
}

async function askChat() {
  const question = String(document.getElementById("chat-question").value || "").trim();
  const jobIdField = document.getElementById("results-job-id");
  const jobId = jobIdField ? jobIdField.value.trim() : "";
  const status = document.getElementById("chat-status");
  if (!question) return;
  // Clear before the request, not after it, so the old figure goes at the click
  // rather than when the answer arrives.
  clearChatOutput();
  if (!jobId) { status.textContent = "Load a dataset first."; return; }
  await ensureFeatureAnnotations(jobId);
  status.textContent = "Working...";
  try {
    const response = await fetch(apiPath(`/api/jobs/${jobId}/chat`), {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ question }),
    });
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail || `chat failed (${response.status})`);
    chatLastResult = data;
    renderChatAnswer(data);
    status.textContent = "";
  } catch (err) {
    status.textContent = err.message;
  }
}

function renderChatAnswer(data) {
  const reading = data.reading || {};
  const bits = [];
  if (reading.cell_state) bits.push(`cell state <b>${reading.cell_state}</b>`);
  if (reading.cell_state_2) bits.push(`and <b>${reading.cell_state_2}</b>`);
  if (reading.contrast) bits.push(`comparison <b>${reading.contrast}</b>`);
  if ((reading.genes || []).length) bits.push(`genes <b>${reading.genes.join(", ")}</b>`);
  document.getElementById("chat-answer").innerHTML =
    `<p>${data.answer || ""}</p><p class="panel-copy">Read as <b>${data.intent}</b>`
    + (bits.length ? `, ${bits.join(", ")}.` : ".") + "</p>";

  // Follow-up questions: each is answerable by this dataset, so a click never
  // lands on a protocol with no recipe.
  const followHost = document.getElementById("chat-followups");
  if (followHost) {
    followHost.innerHTML = "";
    (data.follow_ups || []).forEach((question, index) => {
      if (index === 0) {
        const label = document.createElement("span");
        label.textContent = "Next:";
        followHost.appendChild(label);
      }
      const button = document.createElement("button");
      button.type = "button";
      button.className = "ghost-btn chat-example";
      button.textContent = question;
      button.addEventListener("click", () => {
        document.getElementById("chat-question").value = question;
        askChat();
      });
      followHost.appendChild(button);
    });
  }

  const table = data.table;
  const views = document.getElementById("chat-views");
  if (!table || !(table.rows || []).length) { views.classList.add("hidden"); return; }
  views.classList.remove("hidden");
  initChatResultControls(data);
  document.getElementById("chat-table").classList.remove("hidden");
  document.getElementById("chat-plot").classList.add("hidden");
  views.querySelectorAll("button").forEach((b, i) => b.classList.toggle("active", i === 0));
  if (["cross_modal", "modality_markers"].includes(data.plot?.kind)) {
    document.getElementById("chat-table").classList.add("hidden");
    document.getElementById("chat-plot").classList.remove("hidden");
    views.querySelectorAll("button").forEach(b => b.classList.toggle("active", b.dataset.view === "plot"));
    drawChatPlot();
  }
}

// One result set drives the table and correlation plots; filters never discard it.
let chatResultState = null;
function chatResultRows() {
  const columns = chatLastResult?.table?.columns || [];
  return (chatLastResult?.table?.rows || []).map(row => Array.isArray(row)
    ? Object.fromEntries(columns.map((key, i) => [key, row[i]])) : row);
}
function filteredChatRows() {
  const state = chatResultState;
  if (!state) return chatResultRows();
  const query = state.search.trim().toLocaleLowerCase();
  const rows = chatResultRows().filter(row => {
    if (query && ![...Object.values(row), featureDisplayName(row.feature || row.gene || "", row.modality)].some(v => String(v ?? "").toLocaleLowerCase().includes(query))) return false;
    if ((state.requiredModalities || []).some(m => !(row[m] > 0))) return false;
    if (!state.correlation) return true;
    const rho = row.rho;
    if (!Number.isFinite(rho)) return false;
    if (state.relationship === "positive" && rho <= 0) return false;
    if (state.relationship === "negative" && rho >= 0) return false;
    if (state.relationship === "near_zero" && Math.abs(rho) > state.threshold) return false;
    return (state.min === "" || rho >= Number(state.min)) && (state.max === "" || rho <= Number(state.max));
  });
  const direction = state.direction === "desc" ? -1 : 1;
  if (state.sort) rows.sort((a, b) => {
    const x = a[state.sort], y = b[state.sort];
    if (x == null) return y == null ? 0 : 1;
    if (y == null) return -1;
    return direction * (typeof x === "number" && typeof y === "number"
      ? x - y : String(x).localeCompare(String(y), undefined, {numeric: true}));
  });
  return rows;
}
function filteredChatPairs() {
  const pairs = chatLastResult?.plot?.pairs || [];
  const lookup = new Map(pairs.map(pair => [JSON.stringify([pair.feature, pair.gene]), pair]));
  return filteredChatRows().map(row => lookup.get(JSON.stringify([row.feature, row.gene]))).filter(Boolean);
}
function initChatResultControls(data) {
  const host = document.getElementById("chat-result-controls");
  const defaults = data.result_controls || {};
  const columns = data.table.columns || Object.keys(data.table.rows[0] || {});
  chatResultState = {search:"", sort:defaults.sort_by || "", direction:defaults.sort_direction || "asc",
    limit:data.plot?.kind === "cross_pathways" ? 25 : 50, page:0, requiredModalities:[], relationship:defaults.relationship || "all", threshold:defaults.near_zero_threshold ?? .2,
    min:"", max:"", correlation:data.plot?.kind === "cross_modal", selectedPair:null};
  host.replaceChildren();
  host.classList.remove("hidden");
  const field = (text, control, id) => {
    const label = document.createElement("label");
    label.className = "chat-result-field";
    const caption = document.createElement("span"); caption.textContent = text;
    control.id = id; label.append(caption, control); host.appendChild(label); return control;
  };
  const input = (label, id, key, type="text") => {
    const el = document.createElement("input"); el.type = type; el.value = chatResultState[key];
    if (type === "number") { el.min = key === "threshold" ? "0" : "-1"; el.max = "1"; el.step = "any"; }
    field(label, el, id);
    el.addEventListener("input", () => {
      if (!el.checkValidity() || (key === "threshold" && el.value === "")) return;
      chatResultState[key] = key === "threshold" ? Number(el.value) : el.value;
      chatResultState.page = 0; refreshChatResults();
    });
    return el;
  };
  const select = (label, id, key, choices) => {
    const el = document.createElement("select");
    choices.forEach(([value, text]) => { const o=document.createElement("option"); o.value=value; o.textContent=text; el.appendChild(o); });
    el.value = String(chatResultState[key]); field(label, el, id);
    el.addEventListener("change", () => {
      chatResultState[key] = key === "limit" ? Number(el.value) : el.value;
      if (key === "relationship" && el.value === "near_zero") {
        chatResultState.sort="abs_rho"; chatResultState.direction="asc";
      }
      chatResultState.page=0; refreshChatResults();
    });
    return el;
  };
  input("Search results", "chat-results-search", "search", "search");
  if (data.plot?.kind === "cross_pathways") {
    const group=document.createElement("fieldset");group.className="chat-modality-filters";
    const legend=document.createElement("legend");legend.textContent="Require a hit in every checked modality";group.appendChild(legend);
    data.plot.modalities.forEach(mod=>{
      const label=document.createElement("label"),box=document.createElement("input");box.type="checkbox";box.value=mod.id;
      box.dataset.modality=mod.id; label.append(box, document.createTextNode(mod.label));group.appendChild(label);
      box.addEventListener("change",()=>{chatResultState.requiredModalities=[...group.querySelectorAll("input:checked")].map(e=>e.value);chatResultState.page=0;refreshChatResults();});
    });host.appendChild(group);
    const pdf=document.createElement("button");pdf.type="button";pdf.className="ghost-btn";pdf.id="chat-pathways-pdf";pdf.textContent="Download PDF";
    pdf.addEventListener("click",async()=>{
      try {
        document.getElementById("chat-table").classList.add("hidden");document.getElementById("chat-plot").classList.remove("hidden");
        document.querySelectorAll("#chat-views button").forEach(b=>b.classList.toggle("active",b.dataset.view === "plot"));
        await drawChatPlot();
        await exportPlotlyElementToVectorPdf("chat-plot",buildPdfFilename([getResultsJobId(),data.plot.cell_state,"cross_modality_pathways"]));
      } catch(error) {showDownloadError(error);}
    });host.appendChild(pdf);
  }
  if (chatResultState.correlation) {
    select("Relationship", "chat-results-relationship", "relationship", [["all","All pairs"],["positive","Positive ρ"],["negative","Negative ρ"],["near_zero","Near zero"]]);
    input("Near zero: |ρ| ≤", "chat-results-threshold", "threshold", "number").disabled = chatResultState.relationship !== "near_zero";
    input("Minimum ρ", "chat-results-min", "min", "number");
    input("Maximum ρ", "chat-results-max", "max", "number");
  }
  select("Sort by", "chat-results-sort", "sort", [["","Returned order"], ...columns.map(c=>[c, c === "abs_rho" ? "|ρ| (correlation strength)" : (data.column_labels?.[c] || c)])]);
  select("Order", "chat-results-direction", "direction", [["asc","Ascending"],["desc","Descending"]]);
  select("Rows per page", "chat-results-limit", "limit", [25,50,100,200,500,1000,2000,0].map(n=>[String(n),n ? String(n) : "All"]));
  const pagination=document.createElement("div"); pagination.className="chat-result-pages";
  for (const [id,label,step] of [["prev","Previous",-1],["next","Next",1]]) {
    const button=document.createElement("button"); button.type="button"; button.className="ghost-btn";
    button.id=`chat-results-${id}`; button.textContent=label;
    button.addEventListener("click",()=>{chatResultState.page+=step; refreshChatResults();});
    pagination.appendChild(button);
  }
  const count=document.createElement("span"); count.id="chat-results-count"; count.setAttribute("aria-live","polite");
  pagination.appendChild(count);host.appendChild(pagination);
  renderChatResultTable();
}
function refreshChatResults() {
  document.getElementById("chat-results-sort").value=chatResultState.sort;
  document.getElementById("chat-results-direction").value=chatResultState.direction;
  const threshold=document.getElementById("chat-results-threshold");
  if (threshold) threshold.disabled=chatResultState.relationship !== "near_zero";
  renderChatResultTable();
  if ((chatLastResult?.plot?.filterable || chatResultState.correlation || ["modality_markers","cross_pathways"].includes(chatLastResult?.plot?.kind)) && !document.getElementById("chat-plot").classList.contains("hidden")) drawChatPlot();
}
function renderChatResultTable() {
  const host=document.getElementById("chat-table"), state=chatResultState;
  const rows=filteredChatRows(), columns=chatLastResult.table.columns || Object.keys(rows[0] || {});
  const limit=state.limit || Math.max(rows.length,1);
  state.page=Math.max(0,Math.min(state.page,Math.ceil(rows.length/limit)-1));
  const start=state.page*limit, end=Math.min(start+limit,rows.length);
  document.getElementById("chat-results-prev").disabled=state.page === 0;
  document.getElementById("chat-results-next").disabled=end >= rows.length;
  document.getElementById("chat-results-count").textContent=`${rows.length ? start+1 : 0}–${end} of ${rows.length} matching / ${chatResultRows().length} returned results`;
  host.replaceChildren();
  const table=document.createElement("table"), head=document.createElement("thead"), tr=document.createElement("tr"), body=document.createElement("tbody");
  columns.forEach(key=>{
    const th=document.createElement("th"), button=document.createElement("button");
    const active=state.sort === key;
    th.setAttribute("aria-sort",active ? (state.direction === "asc" ? "ascending" : "descending") : "none");
    button.type="button"; button.className="chat-sort-button";
    button.textContent=(key === "abs_rho" ? "|ρ|" : (chatLastResult.column_labels?.[key] || key))+(active ? (state.direction === "asc" ? " ↑" : " ↓") : " ↕");
    button.addEventListener("click",()=>{state.direction=active && state.direction === "asc" ? "desc" : "asc";state.sort=key;state.page=0;refreshChatResults();});
    th.appendChild(button);tr.appendChild(th);
  });
  head.appendChild(tr);
  rows.slice(start,end).forEach(row=>{
    const tr=document.createElement("tr");
    columns.forEach(key=>{const td=document.createElement("td");td.textContent=["rho","abs_rho"].includes(key) && Number.isFinite(row[key]) ? row[key].toFixed(5) : formatChatCell(row[key]);td.title=String(row[key] ?? "");if (["feature","gene"].includes(key) && featureAnnotation(row[key], row.modality)) { td.textContent=featureDisplayName(row[key], row.modality); const info=featureAnnotation(row[key],row.modality); td.title=`${info.status} · ${info.source} · ${info.source_cell}`; }if (key === "pathway" && chatLastResult.plot?.kind === "cross_pathways") {
      const button=document.createElement("button");button.type="button";button.className="ghost-btn";button.textContent=row.pathway;
      button.title=`Open ${row.id} in ${chatLastResult.plot.source === "marker" ? "Explore" : "Differential"}`;
      button.addEventListener("click",()=>openCrossPathway(row,chatLastResult.plot));td.replaceChildren(button);
    } else if (row.hits?.[key]) td.title=row.hits[key].join("\n");
    tr.appendChild(td);});
    if (state.correlation) {
      const td=document.createElement("td"), button=document.createElement("button");button.type="button";button.className="ghost-btn";button.textContent="Plot";
      button.setAttribute("aria-label",`Plot ${row.feature} / ${row.gene}`);
      button.addEventListener("click",()=>{state.selectedPair=JSON.stringify([row.feature,row.gene]);document.querySelector('#chat-views [data-view="plot"]').click();});
      td.appendChild(button);tr.appendChild(td);
    }
    body.appendChild(tr);
  });
  if (state.correlation) {const th=document.createElement("th");th.textContent="View";tr.appendChild(th);}
  table.append(head,body);host.appendChild(table);
  if (!rows.length) {const message=document.createElement("p");message.textContent="No results match these filters. Change the search or selected filters.";host.appendChild(message);}
}

function formatChatCell(value) {
  if (value === null || value === undefined) return "";
  if (typeof value === "number") {
    if (value !== 0 && Math.abs(value) < 0.001) return value.toExponential(2);
    return String(Math.round(value * 1000) / 1000);
  }
  return String(value);
}

/* The plot view reuses the DotPlot the Explore tab draws, so a chat answer can
 * never introduce a figure the rest of the tool does not produce. */

/* DotPlot: one dot per (gene, cell state). Colour is the mean of the expression
 * layer, size is the fraction of cells in which the gene is detected. Cell
 * states run in the dataset's canonical order, the centroid ordering, so this
 * figure reads the same way as every other plot in the tool. */
function renderDotPlotFigure(hostId, payload) {
  const genes = payload.genes || [];
  const states = payload.states || [];
  const mean = payload.mean || [];
  const frac = payload.frac || [];
  const host = document.getElementById(hostId);
  if (!host) return;
  if (!genes.length || !states.length || (payload.state_n && !payload.state_n.some(n => n > 0))) {
    try { Plotly.purge(hostId); } catch (err) { /* no previous plot */ }
    host.innerHTML = "<span class=\"warn\">No observations match the selected genes and filters.</span>";
    return;
  }
  // scALABLE writes plain HTML into this div for its status messages. Plotly
  // still believes it owns the container after that, and the redraw comes back
  // blank, so the container is released first.
  try { Plotly.purge(hostId); } catch (err) { /* nothing drawn there yet */ }
  let maxMean = 0, minMean = 0;
  mean.forEach((row) => row.forEach((v) => { if (v > maxMean) maxMean = v; if (v < minMean) minMean = v; }));
  const x = [], y = [], size = [], color = [], text = [];
  genes.forEach((gene, gi) => {
    states.forEach((state, si) => {
      if (payload.state_n && !payload.state_n[si]) return;
      x.push(si);
      y.push(gi);
      const f = (frac[gi] || [])[si] || 0;
      const m = (mean[gi] || [])[si] || 0;
      size.push(4 + f * 18);
      color.push(m);
      text.push(`${featureHover(gene, payload.modality)}<br>${state}<br>mean ${m.toFixed(3)}<br>values > 0: ${(f * 100).toFixed(0)}%`);
    });
  });
  host.style.overflowX = "auto";
  const plotWidth = Math.max(host.clientWidth || 500, 440 + 28 * states.length);
  Plotly.react(hostId, [{
    type: "scatter", mode: "markers", x, y, text, hoverinfo: "text",
    marker: {
      size, color,
      // White at no expression through to red at the highest mean, so an
      // unexpressed gene reads as absent rather than as a pale colour.
      colorscale: [[0, "#FFFFFF"], [0.5, "#F4A582"], [1, "#B2182B"]],
      cmin: minMean, cmax: maxMean > minMean ? maxMean : minMean + 1,
      colorbar: { title: { text: "mean", side: "right" }, thickness: 10 },
      line: { width: 0 },
    },
  }], {
    width: plotWidth,
    autosize: false,
    margin: { l: genes.some(g => featureAnnotation(g, payload.modality)) ? 370 : 180, r: 50, t: 10, b: 180 },
    xaxis: {
      tickvals: states.map((_, i) => i), ticktext: states,
      tickangle: -90, tickfont: { size: 9 }, showgrid: false,
      range: [-0.6, states.length - 0.4],
    },
    yaxis: {
      tickvals: genes.map((_, i) => i), ticktext: genes.map(g => featureDisplayName(g, payload.modality)),
      tickfont: { size: 10 }, showgrid: false, range: [-0.6, genes.length - 0.4],
    },
    height: Math.max(240, 22 * genes.length + 170),
  }, { responsive: true, displaylogo: false });
}

/* CombPlot shows individual observations by default; donor means are optional.
 * State names label the top colour strip. Column details belong only in hover.
 * Explicit subplot domains keep each gene and annotation band on its own row. */
function renderCombPlotFigure(hostId, payload) {
  const genes = payload.genes || [];
  const columns = payload.columns || [];
  const values = payload.values || [];
  const colors = payload.colors || [];
  const columnUnit = payload.unit === "cells" ? (payload.observation_unit || "cells") : "donor groups";
  const groupLabel = payload.group_label || "cell state";
  const host = document.getElementById(hostId);
  if (!host) return;
  if (!genes.length || !columns.length) {
    host.innerHTML = "<span class=\"warn\">No observations match the selected filters.</span>";
    return;
  }

  try { Plotly.purge(hostId); } catch (err) { /* nothing drawn there yet */ }
  const nCols = columns.length;
  const x = columns.map((_, i) => i);
  const hover = columns.map((c) =>
    payload.unit === "cells"
      ? `${c.cell}<br>${c.group}${c.donor ? `<br>${c.donor}` : ""}`
      : `${c.donor}<br>${c.group}<br>${c.n_cells} cells`);

  // Where each group's block of donors starts and ends, for the strip and ticks.
  const firstAt = new Map(), lastAt = new Map();
  columns.forEach((c, i) => {
    if (!firstAt.has(c.group)) firstAt.set(c.group, i);
    lastAt.set(c.group, i);
  });
  const blocks = [];
  firstAt.forEach((from, group) => {
    const to = lastAt.get(group);
    blocks.push({ group, center: (from + to) / 2, span: to - from + 1 });
  });
  blocks.sort((a, b) => a.center - b.center);

  /* The annotation bands under the gene rows: one band per covariate LEVEL,
   * with a mark on every column that carries that level. The server sends one
   * value per column per covariate, so a level with no column draws no band. */
  const trackTable = payload.tracks || {};
  const trackLevels = payload.track_levels || {};
  const trackNames = payload.track_names || Object.keys(trackTable);
  const trackRows = [];
  trackNames.forEach((name) => {
    const perColumn = trackTable[name] || [];
    if (perColumn.length !== nCols) return;
    const levels = trackLevels[name] || Array.from(new Set(perColumn));
    levels.forEach((level) => {
      const at = [];
      for (let i = 0; i < nCols; i += 1) if (perColumn[i] === level) at.push(i);
      if (at.length) {
        trackRows.push({ name, level, at, label: `${name}: ${level === "" ? "unrecorded" : level}` });
      }
    });
  });

  /* Left margin has to hold the longest band label, which is longer than any
   * gene name. Measured from the character count rather than guessed, so a
   * bundle with longer level names still fits. */
  const labelChars = trackRows.reduce((m, r) => Math.max(m, r.label.length), 0);
  const marginL = Math.min(190, Math.max(84, Math.round(labelChars * 4.7) + 40));
  const marginR = 14;

  /* Which state names get a label. With 39 states over one panel the -45 degree
   * labels collide wherever a state holds few donors, so labels are placed
   * widest-block-first and a name is skipped when its block centre falls within
   * MIN_LABEL_PX of a name already placed. Wide blocks therefore keep their
   * label and narrow ones lose it; every column still names its state on hover,
   * both on the colour strip and on every bar. */
  const plotW = Math.max(240, (host.clientWidth || 1000) - marginL - marginR);
  const pxPerColumn = plotW / Math.max(1, nCols);
  const MIN_LABEL_PX = 16;
  const centrePx = blocks.map((b) => (b.center + 0.5) * pxPerColumn);
  const byWidth = blocks.map((_, i) => i)
    .sort((a, b) => (blocks[b].span - blocks[a].span) || (a - b));
  const placed = [];
  byWidth.forEach((i) => {
    if (placed.every((j) => Math.abs(centrePx[i] - centrePx[j]) >= MIN_LABEL_PX)) placed.push(i);
  });
  placed.sort((a, b) => a - b);
  const tickvals = placed.map((i) => blocks[i].center);
  const ticktext = placed.map((i) => blocks[i].group);

  /* Vertical budget in pixels, converted to axis domains, so a row keeps the
   * same height whatever the gene count. Top to bottom: the colour strip, one
   * row per gene, then the annotation bands. */
  const STRIP_PX = 16, GENE_PX = 88, GAP_PX = 10, TRACK_PX = 11;
  const bodyPx = STRIP_PX + genes.length * GENE_PX + trackRows.length * TRACK_PX
    + GAP_PX * (genes.length + (trackRows.length ? 1 : 0));
  /* The top margin is exactly the height the rotated labels need. The figure
   * used a fixed 96 px, which left the names floating over empty space. */
  const longest = ticktext.reduce((m, t) => Math.max(m, t.length), 0);
  const marginT = Math.min(150, Math.max(16, Math.round(longest * 5.4 * 0.7071) + 8));
  const marginB = 42;

  const frac = (px) => px / bodyPx;
  const traces = [];
  const layout = {
    showlegend: false, bargap: 0, hovermode: "closest",
    margin: { l: marginL, r: marginR, t: marginT, b: marginB },
    height: Math.round(bodyPx + marginT + marginB),
  };

  /* The colour strip: one cell per column, coloured by its group. It sits on
   * yaxis, a name Plotly accepts. The earlier code called it "ystrip", which is
   * not a legal axis id, so the strip fell back onto the default axis and the
   * label axis anchored to nothing. */
  let top = 1;
  // The strip is one rectangle per group block, not per column: a run of the
  // same colour is a single mark. 1324 bars become 24 on this dataset.
  const sx = [], sw = [], scolor = [], stext = [];
  let blockFrom = 0;
  columns.forEach((c, i) => {
    if (i === columns.length - 1 || columns[i + 1].group !== c.group) {
      const span = i - blockFrom + 1;
      sx.push((blockFrom + i) / 2);
      sw.push(span);
      scolor.push(colors[blockFrom]);
      stext.push(`${c.group}<br>${span} ${columnUnit}`);
      blockFrom = i + 1;
    }
  });
  traces.push({
    type: "bar", x: sx, y: sx.map(() => 1), width: sw,
    marker: { color: scolor, line: { width: 0 } },
    hovertext: stext, textposition: "none", hoverinfo: "text",
    xaxis: "x", yaxis: "y",
  });
  layout.yaxis = {
    domain: [top - frac(STRIP_PX), top], range: [0, 1], autorange: false,
    showticklabels: false, showgrid: false, zeroline: false, fixedrange: true,
  };
  top -= frac(STRIP_PX);

  // The state names, on their own axis anchored to the top of the strip, so
  // they sit against the colours they annotate.
  layout.xaxis2 = {
    overlaying: "x", side: "top", anchor: "y",
    tickvals, ticktext, tickangle: -45, tickfont: { size: 9 },
    showgrid: false, zeroline: false, showline: false, ticks: "",
    range: [-0.5, nCols - 0.5], fixedrange: true,
  };
  traces.push({
    type: "scatter", mode: "markers", x: tickvals, y: tickvals.map(() => 1),
    marker: { opacity: 0 }, hoverinfo: "skip",
    xaxis: "x2", yaxis: "y", showlegend: false,
  });

  genes.forEach((gene, gi) => {
    top -= frac(GAP_PX);
    const axis = `y${gi + 2}`;
    traces.push({
      type: "bar", x, y: values[gi], width: 1,
      marker: { color: colors, line: { width: 0 } },
      hovertext: hover.map(value => `${featureHover(gene, payload.modality)}<br>${value}`), textposition: "none", hoverinfo: "text+y",
      xaxis: "x", yaxis: axis,
    });
    layout[`yaxis${gi + 2}`] = {
      domain: [Math.max(0, top - frac(GENE_PX)), top],
      title: { text: gene, font: { size: 11 } },
      rangemode: "tozero", zeroline: true, showgrid: false,
      tickfont: { size: 9 },
    };
    top -= frac(GENE_PX);
  });

  /* One band per level, drawn as bars rather than a heatmap: a Plotly heatmap
   * exports as a raster image, and these figures have to stay vector. Every
   * mark is a rectangle in the SVG. */
  const trackAxis = `y${genes.length + 2}`;
  if (trackRows.length) {
    top -= frac(GAP_PX);
    // One bar per contiguous run, not one per column. On this dataset the
    // covariate marks fall from 5,296 to 611, 88.5% fewer, because a level
    // almost always occupies a stretch of neighbouring donors. The figure is
    // identical; only the number of SVG rectangles changes.
    const tx = [], ty = [], tbase = [], ttext = [], twidth = [];
    trackRows.forEach((row, k) => {
      const floorAt = trackRows.length - 1 - k;   // first covariate on top
      let runStart = null, previous = null;
      const flush = () => {
        if (runStart === null) return;
        const span = previous - runStart + 1;
        tx.push((runStart + previous) / 2);
        ty.push(0.82);
        tbase.push(floorAt + 0.09);
        twidth.push(span);
        ttext.push(`${row.label}<br>${span} ${columnUnit}`);
        runStart = null;
      };
      row.at.forEach((i) => {
        if (previous !== null && i === previous + 1) { previous = i; return; }
        flush();
        runStart = i;
        previous = i;
      });
      flush();
    });
    traces.push({
      type: "bar", x: tx, y: ty, base: tbase, width: twidth,
      marker: { color: "#2B2B2B", line: { width: 0 } },
      hovertext: ttext, textposition: "none", hoverinfo: "text",
      xaxis: "x", yaxis: trackAxis,
    });
    layout[`yaxis${genes.length + 2}`] = {
      domain: [Math.max(0, top - frac(trackRows.length * TRACK_PX)), top],
      range: [0, trackRows.length], autorange: false,
      tickvals: trackRows.map((_, k) => trackRows.length - 1 - k + 0.5),
      ticktext: trackRows.map((r) => r.label),
      tickfont: { size: 8 }, ticks: "", showgrid: false, zeroline: false,
      fixedrange: true,
    };
  }

  const dropped = blocks.length - placed.length;
  const note = dropped > 0
    ? `; ${dropped} of ${blocks.length} state names hidden where blocks are too narrow to label`
    : "";
  layout.xaxis = {
    domain: [0, 1], anchor: trackRows.length ? trackAxis : `y${genes.length + 1}`,
    showticklabels: false, showgrid: false,
    zeroline: false, range: [-0.5, nCols - 0.5], fixedrange: true,
    title: { text: `${nCols} ${columnUnit}, ordered by ${groupLabel}${note}`,
             font: { size: 10 } },
  };
  Plotly.react(hostId, traces, layout, { responsive: true, displaylogo: false });
}

/* Group-by and level filtering for the DotPlot and CombPlot.
 *
 * The variables offered come from /plot-variables, which lists only the
 * categorical columns with a workable number of levels. The whole obs table
 * would include per-cell numerics such as n_counts, and grouping by one of
 * those would draw a column per cell.
 */
let plotVariablesCache = null;

// Cell states that carry replicate profiles, i.e. the states the released object holds.
// loadPlotVariables reports them as the cluster variable's values.
let replicateStates = null;

async function loadPlotVariables(jobId) {
  if (plotVariablesCache && plotVariablesCache.jobId === jobId) return plotVariablesCache;
  // apiPath honours CELLHARMONY_ROOT_PATH; a bare path 404s wherever the app
  // is mounted under a prefix, which is how the Group by list stayed empty.
  const response = await fetch(apiPath(`/api/jobs/${jobId}/plot-variables`));
  const data = await response.json();
  if (!response.ok) throw new Error(data.detail || "plot variables unavailable");
  plotVariablesCache = { jobId, ...data };
  const clusterVariable = (data.variables || []).find(
    (variable) => String(variable.field) === String(data.cluster_key));
  if (clusterVariable && Array.isArray(clusterVariable.values)) {
    replicateStates = new Set(clusterVariable.values.map(String));
  }
  return plotVariablesCache;
}

async function refreshGroupControls(panelKey) {
  const jobField = document.getElementById("results-job-id");
  const jobId = jobField ? jobField.value.trim() : "";
  const groupBy = document.getElementById(panelElementId(panelKey, "groupby"));
  const groups = document.getElementById(panelElementId(panelKey, "groups"));
  if (!jobId || !groupBy || !groups) return;
  let info;
  try {
    info = await loadPlotVariables(jobId);
  } catch (err) {
    // Swallowing this left an empty Group by select with nothing said, so the
    // control looked broken rather than unavailable.
    console.warn("plot variables unavailable:", err);
    groupBy.innerHTML = '<option value="">unavailable</option>';
    return;
  }

  if (!groupBy.options.length) {
    info.variables.forEach((variable) => {
      const option = document.createElement("option");
      option.value = variable.field;
      option.textContent = `${variable.field} (${variable.n})`;
      groupBy.appendChild(option);
    });
    groupBy.value = info.cluster_key;
    groupBy.addEventListener("change", () => {
      fillGroupLevels(panelKey, info);
      loadVisualizationPanel(panelKey);
    });
    groups.addEventListener("change", () => loadVisualizationPanel(panelKey));
  }
  fillGroupLevels(panelKey, info);
}

/* The levels of the chosen variable. Nothing is selected to start with, which
 * the server reads as "every level", so the plot opens complete. */
function fillGroupLevels(panelKey, info) {
  const groupBy = document.getElementById(panelElementId(panelKey, "groupby"));
  const groups = document.getElementById(panelElementId(panelKey, "groups"));
  if (!groupBy || !groups) return;
  const variable = info.variables.find((v) => v.field === groupBy.value);
  groups.innerHTML = "";
  (variable ? variable.values : []).forEach((level) => {
    const option = document.createElement("option");
    option.value = level;
    option.textContent = level;
    groups.appendChild(option);
  });
}

/* The Color by and Coordinates lists for the UMAP cell-type view.
 *
 * Both come from /plot-variables, which offers the categorical obs columns with
 * a workable number of levels and every obsm entry with two or more columns. The
 * empty value is the default in each list: the cellHarmony cell-state assignment
 * and the cellHarmony projection, which is what the panel drew before.
 */
async function refreshUmapOptions(panelKey) {
  const jobId = getResultsJobId();
  const colorBy = document.getElementById(panelElementId(panelKey, "colorby"));
  const coords = document.getElementById(panelElementId(panelKey, "coords"));
  if (!jobId || !colorBy || !coords) return;
  if (colorBy.dataset.jobId === jobId && coords.dataset.jobId === jobId) return;
  let info;
  try {
    info = await loadPlotVariables(jobId);
  } catch (err) {
    // Say it rather than leaving two empty lists that look broken.
    console.warn("UMAP options unavailable:", err);
    colorBy.innerHTML = '<option value="">unavailable</option>';
    coords.innerHTML = '<option value="">unavailable</option>';
    return;
  }

  const previousColor = colorBy.value;
  colorBy.innerHTML = "";
  const clusterOption = document.createElement("option");
  clusterOption.value = "";
  clusterOption.textContent = `${info.cluster_key} (cellHarmony)`;
  colorBy.appendChild(clusterOption);
  (info.color_variables || info.variables || []).forEach((variable) => {
    if (String(variable.field) === String(info.cluster_key)) return;
    const option = document.createElement("option");
    option.value = variable.field;
    option.textContent = `${variable.field} (${variable.n})`;
    colorBy.appendChild(option);
  });
  colorBy.value = Array.from(colorBy.options).some((o) => o.value === previousColor) ? previousColor : "";
  colorBy.dataset.jobId = jobId;

  const previousCoords = coords.value;
  coords.innerHTML = "";
  const coordOptions = (info.coords || []).length ? info.coords : [{ key: "", label: "cellHarmony UMAP" }];
  const numeric = info.numeric_variables || [];
  // Two numeric annotations are the minimum for a pair of axes.
  // Nathan, 2026-09-01: "Coordinates has to be 2 fields not 1". The X and Y
  // selects below ARE the coordinate control, so no meta-option is offered and
  // nothing hides them.
  const coordEntries = coordOptions;
  coordEntries.forEach((entry) => {
    const option = document.createElement("option");
    option.value = entry.key || "";
    option.textContent = entry.label || entry.key || "cellHarmony UMAP";
    coords.appendChild(option);
  });
  coords.value = Array.from(coords.options).some((o) => o.value === previousCoords) ? previousCoords : "";
  coords.dataset.jobId = jobId;

  // The axis lists. A field shows its range and how many cells carry a value,
  // because a metadata column recorded for part of the dataset draws a panel
  // with fewer cells than the UMAP has.
  const axisSelects = [
    [document.getElementById(panelElementId(panelKey, "xfield")), 0],
    [document.getElementById(panelElementId(panelKey, "yfield")), 1],
  ];
  axisSelects.forEach(([select, defaultIndex]) => {
    if (!select) return;
    const previous = select.value;
    select.innerHTML = "";
    numeric.forEach((entry) => {
      const option = document.createElement("option");
      option.value = entry.field;
      const missing = Number(entry.n_missing) > 0 ? `, ${entry.n_missing} without a value` : "";
      option.textContent = `${entry.label || entry.field} (${formatAxisNumber(entry.min)} to ${formatAxisNumber(entry.max)}${missing})`;
      select.appendChild(option);
    });
    if (Array.from(select.options).some((o) => o.value === previous)) {
      select.value = previous;
    } else if (numeric[defaultIndex]) {
      select.value = numeric[defaultIndex].field;
    }
    select.dataset.jobId = jobId;
    if (select.dataset.wired !== "1") {
      select.dataset.wired = "1";
      select.addEventListener("change", () => loadVisualizationPanel(panelKey));
    }
  });

  if (colorBy.dataset.wired !== "1") {
    colorBy.dataset.wired = "1";
    colorBy.addEventListener("change", () => loadVisualizationPanel(panelKey));
  }
  if (coords.dataset.wired !== "1") {
    coords.dataset.wired = "1";
    coords.addEventListener("change", () => {
      syncUmapAxisFields(panelKey);
      loadVisualizationPanel(panelKey);
    });
  }
  syncUmapAxisFields(panelKey);
}

/* The X and Y lists belong to the "obs columns" coordinate choice alone, so they
 * stay hidden while the panel draws an embedding. */
function syncUmapAxisFields(panelKey) {
  const mode = getPanelSelectValue(panelKey, "mode");
  // The X and Y lists ARE the coordinate control, so they show whenever the panel
  // draws cells. Hiding them behind a meta-option is what left Nathan looking at a
  // single "Coordinates" dropdown three times over.
  const showAxes = UMAP_COORD_MODES.has(mode);
  // The single "Coordinates" dropdown is retired: it offered one field where a
  // coordinate set needs two. It stays in the DOM because other code reads its
  // value, but the reader never sees it.
  const legacy = document.getElementById(panelElementId(panelKey, "coords-field"));
  if (legacy) legacy.classList.add("hidden");
  ["xfield-field", "yfield-field"].forEach((suffix) => {
    const field = document.getElementById(panelElementId(panelKey, suffix));
    if (field) field.classList.toggle("hidden", !showAxes);
  });
}

/* Axis ranges read as numbers a person can compare, not as 17 decimal places. */
function formatAxisNumber(value) {
  const number = Number(value);
  if (!Number.isFinite(number)) return "?";
  if (Math.abs(number) >= 1000 || (number !== 0 && Math.abs(number) < 0.01)) {
    return number.toExponential(1);
  }
  return String(Math.round(number * 100) / 100);
}

function appendGeneSetGroupParams(params, panelKey) {
  const groupBy = document.getElementById(panelElementId(panelKey, "groupby"));
  const groups = document.getElementById(panelElementId(panelKey, "groups"));
  if (groupBy && groupBy.value) params.set("group_by", groupBy.value);
  if (groups) {
    Array.from(groups.selectedOptions).forEach((option) => {
      params.append("groups", option.value);
    });
  }
}


function panelGroupParams(panelKey) {
  const groupBy = document.getElementById(panelElementId(panelKey, "groupby"));
  const groups = document.getElementById(panelElementId(panelKey, "groups"));
  const parts = [];
  if (groupBy && groupBy.value) parts.push(`group_by=${encodeURIComponent(groupBy.value)}`);
  if (groups) {
    Array.from(groups.selectedOptions).forEach((option) => {
      parts.push(`groups=${encodeURIComponent(option.value)}`);
    });
  }
  return parts.join("&");
}

/* Draw the figure for the current chat answer.
 *
 * Two defects this replaces. The container was cleared with innerHTML = "",
 * which removes the nodes but leaves Plotly believing the div is still its own,
 * so the next Plotly.react drew nothing and the panel came back blank on the
 * second visit. Plotly.purge is the supported way to release a container.
 *
 * And a differential answer asked for a volcano, which nothing here drew, so
 * that example reported "no plot" however many times it was asked.
 */
async function drawChatPlot() {
  const result = chatLastResult || {};
  const spec = result.plot;
  const host = document.getElementById("chat-plot");
  if (!host) return;

  // Release the container before reusing it, or the second draw is blank.
  try { Plotly.purge("chat-plot"); } catch (err) { /* never drawn yet */ }
  const previousCorrelation = document.getElementById("chat-correlation-plot");
  if (previousCorrelation) Plotly.purge(previousCorrelation);
  host.innerHTML = "";

  if (!spec) {
    host.innerHTML = "<span class=\"warn\">This answer has no figure; the table holds the result.</span>";
    return;
  }
  const jobId = document.getElementById("results-job-id").value.trim();

  if (spec.kind === "cross_pathways") {
    const all=filteredChatRows(), limit=chatResultState.limit || Math.max(1,all.length);
    const rows=all.slice(chatResultState.page*limit,(chatResultState.page+1)*limit);
    if (!rows.length) {host.textContent="No pathways match the selected modality requirements.";return;}
    const traces=spec.modalities.map(mod=>({type:"bar",orientation:"h",name:mod.label,
      x:rows.map(r=>r[mod.id]), y:rows.map(r=>r.id), marker:{color:mod.color},
      customdata:rows.map(r=>[r.id,r.pathway,(r.hits[mod.id] || []).join(", ")]),
      hovertemplate:"%{customdata[1]}<br>Unique features: %{x}<br>%{customdata[2]}<extra>"+mod.label+"</extra>"}));
    await Plotly.newPlot(host,traces,{barmode:"stack",height:Math.max(380,rows.length*32+160),
      margin:{l:310,r:25,t:65,b:85},title:{text:`Cross-modality pathways: ${spec.cell_state}`,font:{size:16}},
      xaxis:{title:{text:"Unique features per pathway, counted within each modality"},rangemode:"tozero",dtick:Math.max(1,Math.ceil(Math.max(...rows.map(r=>r.total_hits))/10))},
      yaxis:{tickvals:rows.map(r=>r.id),ticktext:rows.map(r=>r.pathway),categoryorder:"array",categoryarray:rows.map(r=>r.id).reverse(),automargin:true},
      legend:{orientation:"h",y:1.12}}, {responsive:true,displaylogo:false});
    host.on("plotly_click",event=>{const id=event.points?.[0]?.customdata?.[0];const row=rows.find(r=>r.id===id);if(row)openCrossPathway(row,spec);});
    return;
  }

  if (spec.kind === "modality_markers") {
    const best = new Map();
    filteredChatRows().forEach(row => {
      if (!best.has(row.modality) || best.get(row.modality).marker_r < row.marker_r) best.set(row.modality, row);
    });
    const rows = [...best.values()].sort((a,b) => b.marker_r - a.marker_r);
    if (!rows.length) { host.textContent = "No markers match these filters."; return; }
    const labels = rows.map(r => `${r.modality}: ${featureDisplayName(r.feature, r.modality)}`);
    Plotly.newPlot(host, [{type:"bar", orientation:"h", x:rows.map(r=>r.marker_r), y:labels,
      marker:{color:rows.map(r=>r.annotation === "Unidentified" ? "#64748b" : "#0f807c")},
      hovertemplate:"%{y}<br>MarkerFinder Pearson r: %{x:.4f}<extra></extra>"}],
      {height:Math.max(300,rows.length*55+110), margin:{l:270,r:25,t:45,b:60},
       title:{text:`Best retained marker per modality: ${spec.cell_state}`,font:{size:15}},
       xaxis:{title:{text:"Pearson r with cell-state membership"},rangemode:"tozero"},
       yaxis:{automargin:true,categoryorder:"array",categoryarray:labels.slice().reverse()}},
      {responsive:true,displaylogo:false});
    return;
  }
  if (spec.kind === "cross_modal") {
    const pairs = filteredChatPairs();
    if (!pairs.length) { host.textContent = "No pairs match these filters. Change the relationship, search, or correlation range."; return; }
    const select = document.createElement("select");
    select.setAttribute("aria-label", "Expression correlation pair");
    pairs.forEach((pair, i) => {
      const option = document.createElement("option");
      option.value = String(i);
      option.textContent = `${pair.feature} / ${pair.gene} (Spearman ρ = ${pair.rho.toFixed(3)})`;
      select.appendChild(option);
    });
    const selected = pairs.findIndex(pair => JSON.stringify([pair.feature, pair.gene]) === chatResultState?.selectedPair);
    if (selected >= 0) select.value = String(selected);
    const chart = document.createElement("div");
    chart.id = "chat-correlation-plot";
    host.append(select, chart);
    const draw = () => {
      const pair = pairs[Number(select.value)];
      chatResultState.selectedPair = JSON.stringify([pair.feature, pair.gene]);
      Plotly.react(chart, [{type:"scatter", mode:"markers", x:pair.x, y:pair.y,
        text:pair.labels, customdata:pair.n_cells, marker:{size:9, color:"#0f807c"},
        hovertemplate:"%{text}<br>RNA: %{x:.3f}<br>Measurement: %{y:.3f}<br>Matched observations: %{customdata}<extra></extra>"}],
        {height:440, margin:{l:85,r:25,t:45,b:70},
         title:{text:`${pair.feature} / ${pair.gene}: ${spec.unit}`,font:{size:14}},
         xaxis:{title:{text:`${pair.gene}: mean RNA expression`}},
         yaxis:{title:{text:`${pair.feature}: mean ${spec.y_label}`}}},
        {responsive:true,displaylogo:false});
    };
    select.addEventListener("change", draw);
    draw();
    return;
  }
  if (spec.kind.startsWith("integrated_")) {
    await ScalableIntegrated.mount(host,{...integratedOptions(jobId,spec.kind,spec.cell_state,spec.modality,spec.source),...spec,externalPdf:false});
    return;
  }
  if (spec.kind === "volcano") {
    drawChatVolcano(host, result);
    return;
  }
  // A CombPlot answer: severity_gradient, coexpression_module and dose_response
  // all return one. Without this the panel said "no figure" for three of the
  // protocols that do the most work.
  if (spec.kind === "combplot" && (spec.genes || []).length) {
    try {
      const bits = [`genes=${encodeURIComponent(spec.genes.join(","))}`, "min_cells=5", "unit=donor"];
      if (spec.group_by) bits.push(`group_by=${encodeURIComponent(spec.group_by)}`);
      const response = await fetch(apiPath(`/api/jobs/${jobId}/combplot?${bits.join("&")}`));
      const data = await response.json();
      if (!response.ok) throw new Error(data.detail || "combplot failed");
      renderCombPlotFigure("chat-plot", data);
    } catch (err) {
      host.innerHTML = `<span class="warn">${err.message}</span>`;
    }
    return;
  }
  // A bar chart or a heatmap is drawn from the table the answer already
  // carries, so neither needs a second request.
  if (spec.kind === "gradient") { drawChatGradient(spec); return; }
  if (spec.kind === "frequency") { drawChatFrequency(spec); return; }
  if (spec.kind === "signature") { drawChatSignature(spec); return; }
  if (spec.kind === "barchart") { drawChatBars(result, spec); return; }
  if (spec.kind === "heatmap") { drawChatHeatmap(result); return; }
  if (spec.kind === "dotplot" && (spec.genes || []).length) {
    try {
      const response = await fetch(
        apiPath(`/api/jobs/${jobId}/dotplot?genes=${encodeURIComponent(spec.genes.join(","))}`));
      const data = await response.json();
      if (!response.ok) throw new Error(data.detail || "dotplot failed");
      renderDotPlotFigure("chat-plot", data);
    } catch (err) {
      host.innerHTML = `<span class="warn">${err.message}</span>`;
    }
    return;
  }
  // A regulatory network, drawn from the nodes and edges the answer already carries.
  if (spec.kind === "network") { drawChatNetwork(result, spec); return; }
  host.innerHTML = "<span class=\"warn\">This answer has no figure; the table holds the result.</span>";
}

var chatNetworkCy = null;

/* The regulatory network a differential implies: the features that moved, and the factors
 * the cell state's regulatory model puts above them.
 *
 * THE RULES ARE THE SITE'S, NOT MINE. scalable_viewer/grn_network.py applies them and
 * records where each came from; this only draws what that returns. Nathan's own
 * instruction, quoted there: the targets are the features on screen, and "The TFs do not
 * need to be regulated themselves since their activity but not expression might be
 * changed."
 *
 * SO A NODE HAS THREE COLOUR STATES, NOT TWO. `colour_by` says which number paints it:
 * its differential EXPRESSION, or its differential ACTIVITY when only the activity moved,
 * or neither. A factor that acts without changing its own transcript is the whole reason
 * the activity layer exists, so painting it grey would hide the finding, and painting it
 * as a measured zero would invent one. Activity-coloured nodes carry a dashed rim, which
 * is what the site's legend promises.
 */
function drawChatNetwork(result, spec) {
  const host = document.getElementById("chat-plot");
  const nodes = result.nodes || [];
  const edges = result.edges || [];
  if (!nodes.length) {
    host.innerHTML = `<span class="warn">${result.note || "No regulator reaches these features."}</span>`;
    return;
  }
  if (typeof cytoscape !== "function") {
    host.innerHTML = "<span class=\"warn\">Cytoscape did not load, so the network cannot be drawn.</span>";
    return;
  }
  if (chatNetworkCy) {
    try { chatNetworkCy.destroy(); } catch (err) { /* already gone */ }
    chatNetworkCy = null;
  }

  // A DIVERGING RAMP ON THE FOLD CHANGE, GREY WHERE THERE IS NO ROW. Explicit hex, so the
  // ramp is the same on every machine, and no rainbow.
  const ramp = ["#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b"];
  const fade = "#cbd5e1";
  const colourFor = (node) => {
    const by = String(node.colour_by || "none");
    if (by === "none") return fade;
    const raw = by === "activity" ? node.activity_log2fc : node.expression_log2fc;
    if (raw === null || raw === undefined || raw === "") return fade;
    const fc = Number(raw);
    if (!Number.isFinite(fc)) return fade;
    // Saturate at |2| so one extreme feature does not flatten every other colour.
    const t = Math.max(0, Math.min(1, (fc + 2) / 4));
    return ramp[Math.min(ramp.length - 1, Math.floor(t * ramp.length))];
  };
  const four = (v) => (v !== null && v !== undefined && v !== "" && Number.isFinite(Number(v))
    ? Number(v).toPrecision(3) : "not reported");

  const elements = [];
  nodes.forEach((n) => {
    elements.push({
      data: {
        id: n.id,
        label: n.label || n.id,
        role: n.role,
        colour: colourFor(n),
        shape: n.role === "factor" ? "diamond" : "ellipse",
        // The rim states WHICH number the colour came from.
        rim: n.colour_by === "activity" ? "dashed" : "solid",
        size: n.role === "factor" ? 30 : 22,
        // The hover carries the number the colour did NOT use, so both are readable.
        tip: n.role === "factor"
          ? `${n.id}: activity log2FC ${four(n.activity_log2fc)} (FDR ${four(n.activity_fdr)}), `
            + `expression log2FC ${four(n.expression_log2fc)}, `
            + `${n.n_targets_here} target(s) here`
          : `${n.id}: expression log2FC ${four(n.expression_log2fc)} (FDR ${four(n.expression_fdr)})`,
      },
    });
  });
  edges.forEach((e, i) => {
    elements.push({ data: { id: `e${i}`, source: e.source, target: e.target, score: e.score } });
  });

  // The edge score sets the line width, scaled to THIS network's own range, so a network
  // of weak edges is still readable and a strong one is not all maximum width.
  const range = (spec && spec.edge_score_range) || result.edge_score_range || null;
  const lo = range ? Number(range[0]) : 0;
  const hi = range ? Number(range[1]) : 1;
  const span = hi > lo ? hi - lo : 1;

  chatNetworkCy = cytoscape({
    container: host,
    elements,
    style: [
      {
        selector: "node",
        style: {
          "background-color": "data(colour)",
          shape: "data(shape)",
          label: "data(label)",
          color: "#0f172a",
          "font-size": 11,
          "text-valign": "center",
          "text-halign": "center",
          width: "data(size)",
          height: "data(size)",
          "border-width": 2,
          "border-style": "data(rim)",
          "border-color": "#0f172a",
        },
      },
      {
        selector: "edge",
        style: {
          width: (edge) => {
            const s = Number(edge.data("score"));
            if (!Number.isFinite(s)) return 1.2;
            return Math.max(1, Math.min(9, 1 + ((s - lo) / span) * 8));
          },
          "line-color": "#94a3b8",
          "target-arrow-color": "#94a3b8",
          "target-arrow-shape": "triangle",
          "curve-style": "bezier",
          opacity: 0.8,
        },
      },
      { selector: "node:selected", style: { "border-width": 4, "border-color": "#0f172a" } },
    ],
    layout: {
      name: "cose", animate: true, fit: true, padding: 36, randomize: true,
      idealEdgeLength: 80, nodeOverlap: 8, componentSpacing: 90,
    },
  });
  chatNetworkCy.on("mouseover", "node", (event) => {
    setDifferentialNetworkHoverTooltip(event.target.data("tip"), event.renderedPosition);
  });
  chatNetworkCy.on("mouseout", "node", () => setDifferentialNetworkHoverTooltip(""));
}

/* A volcano drawn from the differential table the answer already carries, so it
 * shows exactly the rows the table shows and needs no second request. */
function drawChatVolcano(host, result) {
  const rows = ((result.table || {}).rows) || [];
  if (!rows.length) {
    host.innerHTML = "<span class=\"warn\">No differential rows to plot.</span>";
    return;
  }
  const get = (row, name) => (Array.isArray(row)
    ? row[((result.table.columns || []).indexOf(name))] : row[name]);
  const x = [], y = [], text = [], color = [];
  rows.forEach((row) => {
    const rawFold = get(row, "log2fc"), rawFdr = get(row, "fdr");
    if (rawFold === null || rawFold === undefined || rawFold === "" || rawFdr === null || rawFdr === undefined || rawFdr === "") return;
    const fc = Number(rawFold), fdr = Number(rawFdr);
    if (!Number.isFinite(fc) || !Number.isFinite(fdr) || fdr < 0 || fdr > 1) return;
    // A reported FDR of 0 is below the smallest float the table stores, so it is
    // drawn at the top of the axis rather than dropped.
    const safe = fdr > 0 ? fdr : 1e-300;
    x.push(fc);
    y.push(-Math.log10(safe));
    text.push(`${get(row, "gene")}<br>log2FC ${fc.toFixed(3)}<br>FDR ${fdr.toExponential(2)}`);
    color.push(fc >= 0 ? "#B2182B" : "#2166AC");
  });
  Plotly.newPlot("chat-plot", [{
    type: "scatter", mode: "markers", x, y, text, hoverinfo: "text",
    marker: { size: 16, color, opacity: 0.85, line: { width: 1, color: "#FFFFFF" } },
  }], {
    margin: { l: 60, r: 20, t: 20, b: 50 },
    xaxis: { title: { text: "log2 fold change", font: { size: 11 } }, zeroline: true },
    yaxis: { title: { text: "-log10 FDR", font: { size: 11 } } },
    height: 360, showlegend: false,
  }, { responsive: true, displaylogo: false });
}

/* The "Filter data to display" rows, translated for the DotPlot and CombPlot.
 *
 * Annotation 1 and Annotation 2 each name one obs column and one level of it.
 * The two gene-set routes read them as subset_by/subset_values and
 * subset2_by/subset2_values, and keep only the cells that satisfy both.
 *
 * A level goes across whole. preGOLD2 carries the level "I,II", so splitting a
 * value on a comma here would ask the server for two levels that do not exist
 * and the plot would come back unfiltered. URLSearchParams encodes the comma.
 *
 * "All" is the empty value, and an empty field or value adds no parameter, so a
 * panel nobody has filtered sends the request it sent before.
 */
function appendGeneSetSubsetParams(params, panelKey) {
  [[1, "subset"], [2, "subset2"]].forEach(([index, prefix]) => {
    const fieldSelect = document.getElementById(panelElementId(panelKey, `filter${index}-field`));
    const valueSelect = document.getElementById(panelElementId(panelKey, `filter${index}-values`));
    if (!fieldSelect || !valueSelect) {
      return;
    }
    const field = String(fieldSelect.value || "").trim();
    const value = String(valueSelect.value || "");
    if (!field || !value) {
      return;
    }
    params.set(`${prefix}_by`, field);
    params.append(`${prefix}_values`, value);
  });
  return params;
}

/* A bar chart from the answer's own table: the last numeric column is the
 * value, the first column is the label. Used by composition_shift, where the
 * value is a log2 ratio and the sign matters. */

/* A heatmap from the answer's table, drawn as rectangles rather than a Plotly
 * heatmap trace: a heatmap trace exports as a raster image, and these figures
 * have to stay vector. Used by annotation_concordance, where the value is the
 * fraction of a cell state carrying each label of the second annotation. */
function drawChatHeatmap(result) {
  const table = result.table || {};
  const rows = table.rows || [];
  const columns = table.columns || [];
  if (!rows.length) {
    document.getElementById("chat-plot").innerHTML =
      "<span class=\"warn\">No rows to plot.</span>";
    return;
  }
  const get = (row, i) => (Array.isArray(row) ? row[i] : row[columns[i]]);
  const yLabels = [...new Set(rows.map((r) => String(get(r, 0))))];
  const xLabels = [...new Set(rows.map((r) => String(get(r, 1))))];
  const valueAt = columns.length - 1;
  let maximum = 0;
  rows.forEach((r) => { maximum = Math.max(maximum, Number(get(r, valueAt)) || 0); });

  const shapes = [], hx = [], hy = [], htext = [];
  rows.forEach((r) => {
    const yi = yLabels.indexOf(String(get(r, 0)));
    const xi = xLabels.indexOf(String(get(r, 1)));
    const value = Number(get(r, valueAt)) || 0;
    const strength = maximum > 0 ? value / maximum : 0;
    // White at zero through to red at the maximum, the same ramp the DotPlot uses.
    const channel = Math.round(255 - strength * (255 - 178));
    shapes.push({
      type: "rect", xref: "x", yref: "y",
      x0: xi - 0.5, x1: xi + 0.5, y0: yi - 0.5, y1: yi + 0.5,
      fillcolor: `rgb(${channel > 178 ? 255 : 178}, ${channel}, ${channel})`,
      line: { width: 0.4, color: "#E6E6E6" }, layer: "below",
    });
    hx.push(xi); hy.push(yi);
    htext.push(`${get(r, 0)}<br>${get(r, 1)}<br>${columns[valueAt]}: ${value}`);
  });
  Plotly.newPlot("chat-plot", [{
    type: "scatter", mode: "markers", x: hx, y: hy,
    marker: { size: 18, opacity: 0 }, text: htext, hoverinfo: "text",
  }], {
    shapes,
    margin: { l: 200, r: 20, t: 20, b: 150 },
    xaxis: { tickvals: xLabels.map((_, i) => i), ticktext: xLabels,
             tickangle: -45, tickfont: { size: 9 }, range: [-0.5, xLabels.length - 0.5],
             showgrid: false, zeroline: false },
    yaxis: { tickvals: yLabels.map((_, i) => i), ticktext: yLabels,
             tickfont: { size: 9 }, range: [-0.5, yLabels.length - 0.5],
             showgrid: false, zeroline: false },
    height: Math.max(280, 22 * yLabels.length + 180), showlegend: false,
  }, { responsive: true, displaylogo: false });
}

/* The figure for an association that came from a statistical test.
 *
 * One panel per gene. Each point is one donor: the clinical variable across,
 * that donor's pseudobulk level of the gene up, and the fitted line through
 * them. The slope and R-squared are printed on the panel, because a rank
 * correlation alone can be high while the gene barely moves.
 *
 * This replaces a CombPlot, which drew all 39 cell states for a correlation
 * computed inside one of them. Only the donors contributing to that one cell
 * state appear here, which is what the question asked about.
 */
function drawChatGradient(spec) {
  const host = document.getElementById("chat-plot");
  const points = spec.points || {};
  const genes = Object.keys(points);
  if (!host) return;
  if (!genes.length) {
    host.innerHTML = "<span class=\"warn\">No per-donor points to plot.</span>";
    return;
  }
  const columns = Math.min(3, genes.length);
  const rows = Math.ceil(genes.length / columns);
  const traces = [];
  const layout = {
    showlegend: false,
    margin: { l: 56, r: 16, t: 44, b: 52 },
    height: Math.max(280, rows * 210 + 60),
    annotations: [],
    grid: { rows, columns, pattern: "independent",
            xgap: 0.16, ygap: 0.3 },
  };
  genes.forEach((gene, index) => {
    const entry = points[gene];
    const axis = index === 0 ? "" : String(index + 1);
    traces.push({
      type: "scatter", mode: "markers",
      x: entry.x, y: entry.y,
      marker: { size: 9, color: "#2166AC", opacity: 0.7,
                line: { width: 0.5, color: "#FFFFFF" } },
      text: (spec.donors || []).map((d, i) =>
        `${d}<br>${spec.covariate}: ${entry.x[i]}<br>${gene}: ${entry.y[i]}`),
      hoverinfo: "text",
      xaxis: `x${axis}`, yaxis: `y${axis}`,
    });
    // The fitted line, as two points: a straight line is one 2-point path.
    const lo = Math.min(...entry.x), hi = Math.max(...entry.x);
    traces.push({
      type: "scatter", mode: "lines",
      x: [lo, hi],
      y: [entry.slope * lo + entry.intercept, entry.slope * hi + entry.intercept],
      line: { color: "#B2182B", width: 2 }, hoverinfo: "skip",
      xaxis: `x${axis}`, yaxis: `y${axis}`,
    });
    layout[`xaxis${axis}`] = {
      title: { text: spec.covariate, font: { size: 10 } }, tickfont: { size: 9 },
      zeroline: false,
    };
    if(spec.levels){layout[`xaxis${axis}`].tickvals=spec.levels.map((_,i)=>i);layout[`xaxis${axis}`].ticktext=spec.levels;}
    layout[`yaxis${axis}`] = {
      title: { text: gene, font: { size: 11 } }, tickfont: { size: 9 },
      zeroline: false, rangemode: "tozero",
    };
    layout.annotations.push({
      text: `rho ${entry.rho}  slope ${entry.slope}  R² ${entry.r2}`,
      xref: `x${axis} domain`, yref: `y${axis} domain`,
      x: 0, y: 1.14, showarrow: false,
      font: { size: 9, color: "#555555" }, xanchor: "left",
    });
  });
  Plotly.newPlot("chat-plot", traces, layout,
                 { responsive: true, displaylogo: false });
}

/* Cell frequency, two groups side by side.
 *
 * One pair of bars per cell state: each group's mean share of a donor's own
 * cells. A log ratio on its own cannot say whether a state is common in both
 * groups or rare in both, and that is usually the first thing worth knowing.
 * A state whose rank test clears 0.05 is marked.
 */

/* Who carries a signature: genes down, donors across, ordered by score.
 *
 * Each cell is a gene's standardised level in one donor, drawn as a rectangle
 * so the figure stays vector; a Plotly heatmap trace exports as a raster. Above
 * the matrix sits each donor's aggregate score and a band naming their group,
 * so a signature carried by a subset shows as colour collecting at one end
 * rather than spreading across the ranking.
 */
function drawChatSignature(spec) {
  const host = document.getElementById("chat-plot");
  const genes = spec.genes || [];
  const donors = spec.donors || [];
  const z = spec.z || [];
  if (!host) return;
  if (!genes.length || !donors.length) {
    host.innerHTML = "<span class=\"warn\">No signature to plot.</span>";
    return;
  }
  const groups = spec.groups || [];
  const levels = [...new Set(groups.filter(Boolean))];
  const groupColour = ["#B2182B", "#2166AC", "#4D9221", "#9970AB"];

  // Colour scale: blue low, white at the mean, red high, clipped at +/- 2 SD.
  const shapes = [];
  const clip = 2;
  z.forEach((row, gi) => {
    row.forEach((value, di) => {
      const t = Math.max(-1, Math.min(1, value / clip));
      const colour = t >= 0
        ? `rgb(${Math.round(255 - t * 77)}, ${Math.round(255 - t * 231)}, ${Math.round(255 - t * 212)})`
        : `rgb(${Math.round(255 + t * 222)}, ${Math.round(255 + t * 152)}, ${Math.round(255 + t * 60)})`;
      shapes.push({
        type: "rect", xref: "x", yref: "y",
        x0: di - 0.5, x1: di + 0.5, y0: gi - 0.5, y1: gi + 0.5,
        fillcolor: colour, line: { width: 0 }, layer: "below",
      });
    });
  });
  // The group band, one rectangle per donor, on its own axis above the matrix.
  groups.forEach((group, di) => {
    const index = Math.max(0, levels.indexOf(group));
    shapes.push({
      type: "rect", xref: "x", yref: "paper",
      x0: di - 0.5, x1: di + 0.5, y0: 1.0, y1: 1.03,
      fillcolor: groupColour[index % groupColour.length], line: { width: 0 },
    });
  });

  const traces = [{
    // The aggregate score, on its own axis above the matrix.
    type: "scatter", mode: "markers", x: donors.map((_, i) => i), y: spec.score || [],
    marker: { size: 5, color: "#333333" },
    text: donors.map((d, i) => `${d}<br>${groups[i] || ""}<br>score ${(spec.score || [])[i]}`),
    hoverinfo: "text", xaxis: "x", yaxis: "y2",
  }, {
    // Invisible points carrying the per-cell hover for the matrix.
    type: "scatter", mode: "markers",
    x: [].concat(...z.map((row, gi) => row.map((_, di) => di))),
    y: [].concat(...z.map((row, gi) => row.map(() => gi))),
    marker: { size: 6, opacity: 0 },
    text: [].concat(...z.map((row, gi) => row.map((v, di) =>
      `${genes[gi]}<br>${donors[di]}<br>z ${v}`))),
    hoverinfo: "text", xaxis: "x", yaxis: "y",
  }];

  Plotly.newPlot("chat-plot", traces, {
    shapes, showlegend: false,
    margin: { l: 110, r: 16, t: 58, b: 30 },
    height: Math.max(360, 13 * genes.length + 190),
    xaxis: { domain: [0, 1], showticklabels: false, showgrid: false, zeroline: false,
             range: [-0.5, donors.length - 0.5],
             title: { text: `${donors.length} donors, ordered by signature score`,
                      font: { size: 10 } } },
    yaxis: { domain: [0, 0.76], tickvals: genes.map((_, i) => i), ticktext: genes,
             tickfont: { size: 8 }, showgrid: false, zeroline: false,
             range: [-0.5, genes.length - 0.5], autorange: "reversed" },
    yaxis2: { domain: [0.82, 1.0], title: { text: "score", font: { size: 10 } },
              tickfont: { size: 8 }, showgrid: false, zeroline: true, anchor: "x" },
    annotations: levels.map((level, i) => ({
      text: level, xref: "paper", yref: "paper",
      x: i * 0.16, y: 1.075, showarrow: false, xanchor: "left",
      font: { size: 9, color: groupColour[i % groupColour.length] },
    })).concat([{
      text: `top ${spec.n_up} up, then the down genes`,
      xref: "paper", yref: "paper", x: 1, y: 1.075, showarrow: false,
      xanchor: "right", font: { size: 9, color: "#555555" },
    }]),
  }, { responsive: true, displaylogo: false });
}

/* A bar chart from the answer's own table.
 *
 * The column to plot is named by the answer (`value_column`), because guessing
 * the rightmost column drew nothing for GO terms, whose last column is a gene
 * list. Where no column is named, the rightmost one that actually parses as a
 * number is used, and if none does the panel says so rather than drawing empty
 * axes.
 */
function drawChatBars(result, spec) {
  const host = document.getElementById("chat-plot");
  const table = result.table || {};
  const allRows = spec.filterable ? filteredChatRows() : (table.rows || []);
  const limit = chatResultState?.limit || Math.max(1, allRows.length);
  const rows = spec.filterable ? allRows.slice(chatResultState.page*limit,(chatResultState.page+1)*limit) : allRows;
  const columns = table.columns || [];
  if (!rows.length) {
    host.innerHTML = "<span class=\"warn\">No rows to plot.</span>";
    return;
  }
  const at = (row, index) => (Array.isArray(row) ? row[index] : row[columns[index]]);
  const numeric = (index) => rows.every((r) => Number.isFinite(Number(at(r, index))));

  let valueAt = columns.indexOf((spec || {}).value_column);
  if (valueAt < 0 || !numeric(valueAt)) {
    valueAt = -1;
    for (let i = columns.length - 1; i >= 1; i -= 1) {
      if (numeric(i)) { valueAt = i; break; }
    }
  }
  if (valueAt < 0) {
    host.innerHTML = "<span class=\"warn\">No numeric column to plot; the table holds the result.</span>";
    return;
  }
  let labelAt = columns.indexOf((spec || {}).label_column);
  if (labelAt < 0) labelAt = 0;
  const signAt = columns.indexOf((spec || {}).sign_column);

  const labels = rows.map((r) => String(at(r, labelAt)));
  const values = rows.map((r) => Number(at(r, valueAt)));
  // A direction column decides the colour when the value itself is unsigned:
  // a GO Z score is always positive, but the term may be up or down.
  const colours = rows.map((r, i) => {
    if (signAt >= 0) {
      const sign = at(r, signAt);
      if (sign === null || sign === undefined || sign === "") return "#cbd5e1";
      if (Number.isFinite(Number(sign))) {
        return Number(sign) < 0 ? "#2166AC" : (Number(sign) > 0 ? "#B2182B" : "#cbd5e1");
      }
      return String(sign).toLowerCase().startsWith("d") ? "#2166AC" : "#B2182B";
    }
    return values[i] >= 0 ? "#B2182B" : "#2166AC";
  });

  Plotly.newPlot("chat-plot", [{
    type: "bar", orientation: "h", x: values, y: labels,
    marker: { color: colours, line: { width: 0 } },
    hovertemplate: `%{y}<br>${columns[valueAt]}: %{x}<extra></extra>`,
  }], {
    margin: { l: 240, r: 20, t: 20, b: 48 },
    xaxis: { title: { text: columns[valueAt], font: { size: 11 } }, zeroline: true },
    yaxis: { automargin: true, tickfont: { size: 9 },
             categoryorder: "array", categoryarray: labels.slice().reverse() },
    height: Math.max(280, 22 * labels.length + 90), showlegend: false,
  }, { responsive: true, displaylogo: false });
}

/* Cell frequency, two groups per cell type.
 *
 * A pair of bars per cell state, sky blue for the reference group and light red
 * for the case group, so the same cell type is read side by side. Each bar
 * carries the standard error of its mean and every donor as a jittered point,
 * because a mean over 141 donors and a mean over 4 are indistinguishable
 * otherwise, and the spread is usually the thing worth seeing.
 *
 * Points are drawn as inline vector circles; nothing here rasterises.
 */
function drawChatFrequency(spec) {
  const host = document.getElementById("chat-plot");
  const states = spec.states || [];
  if (!host) return;
  if (!states.length) {
    host.innerHTML = "<span class=\"warn\">No cell states to plot.</span>";
    return;
  }
  const groups = spec.groups || ["group 1", "group 2"];
  const counts = spec.n_donors || [0, 0];
  const pvals = spec.p || [];

  // The server orders the groups reference first, taking the control side from
  // the contrast's own control_label, so the first series is always the
  // reference. An earlier version re-guessed it here with a regex for
  // control/non/healthy, which matched nothing in "GOLD I, II" against
  // "GOLD IV" and painted the milder stage red.
  const SKY = "#87CEEB";        // reference group, drawn first
  const LIGHT_RED = "#F08080";  // case group
  const colourA = SKY;
  const colourB = LIGHT_RED;

  // Widest first, so the comparison worth making is at the left.
  const LIMIT = 14;
  const shown = states.slice(0, LIMIT);
  const hidden = states.length - shown.length;
  const labels = shown.map((s, i) => (pvals[i] < 0.05 ? `${s} *` : s));

  const traces = [
    { type: "bar", name: `${groups[0]} (${counts[0]})`, x: labels,
      y: (spec.a || []).slice(0, shown.length),
      error_y: { type: "data", array: (spec.a_sem || []).slice(0, shown.length),
                 visible: true, color: "#555555", thickness: 1, width: 3 },
      marker: { color: colourA, line: { color: "#4F4F4F", width: 0.6 } },
      offsetgroup: "a",
      hovertemplate: `%{x}<br>${groups[0]}: %{y:.5f}<extra></extra>` },
    { type: "bar", name: `${groups[1]} (${counts[1]})`, x: labels,
      y: (spec.b || []).slice(0, shown.length),
      error_y: { type: "data", array: (spec.b_sem || []).slice(0, shown.length),
                 visible: true, color: "#555555", thickness: 1, width: 3 },
      marker: { color: colourB, line: { color: "#4F4F4F", width: 0.6 } },
      offsetgroup: "b",
      hovertemplate: `%{x}<br>${groups[1]}: %{y:.5f}<extra></extra>` },
  ];

  // One point per donor, jittered inside its own bar. A deterministic offset,
  // not a random one, so the figure is identical every time it is drawn.
  const jitter = (n, index) => (n <= 1 ? 0 : ((index / (n - 1)) - 0.5) * 0.26);
  [["a_points", -0.2, groups[0]], ["b_points", 0.2, groups[1]]].forEach(
    ([key, shift, name]) => {
      const px = [], py = [], ptext = [];
      (spec[key] || []).slice(0, shown.length).forEach((values, si) => {
        (values || []).forEach((value, di) => {
          px.push(si + shift + jitter(values.length, di));
          py.push(value);
          ptext.push(`${shown[si]}<br>${name}<br>${value}`);
        });
      });
      traces.push({
        type: "scatter", mode: "markers", x: px, y: py,
        marker: { size: 4, color: "#2F2F2F", opacity: 0.45,
                  line: { width: 0 } },
        text: ptext, hoverinfo: "text", showlegend: false,
        xaxis: "x2", yaxis: "y",
      });
    });

  Plotly.newPlot("chat-plot", traces, {
    barmode: "group", bargap: 0.3, bargroupgap: 0.08,
    margin: { l: 74, r: 20, t: 54, b: 168 },
    xaxis: { tickvals: labels.map((_, i) => i), ticktext: labels,
             tickangle: -45, tickfont: { size: 9 }, showgrid: false,
             range: [-0.6, shown.length - 0.4] },
    // A second x-axis on the same range carries the donor points, so their
    // numeric positions are not snapped to the category centres.
    xaxis2: { overlaying: "x", range: [-0.6, shown.length - 0.4],
              showticklabels: false, showgrid: false, zeroline: false },
    yaxis: { title: { text: "share of a donor's cells", font: { size: 11 } },
             tickfont: { size: 9 }, rangemode: "tozero" },
    height: 460,
    legend: { orientation: "h", y: 1.1, x: 0 },
    annotations: [{
      text: `each point is one donor  ·  bars are mean ± SEM  ·  * rank test below 0.05`
            + (hidden > 0 ? `  ·  ${shown.length} of ${states.length} states shown` : ""),
      xref: "paper", yref: "paper", x: 1, y: 1.1, showarrow: false,
      xanchor: "right", font: { size: 9, color: "#555555" },
    }],
  }, { responsive: true, displaylogo: false });
}
