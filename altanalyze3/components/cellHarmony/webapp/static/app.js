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
let currentMarkerAnalysis = null;
let currentMarkerAnalysisByModality = { rna: null };
let currentFastCommAnalysis = null;
let currentModalitiesState = { default: "rna", available: [{ id: "rna", label: "RNA", feature_label: "gene", example_feature: "MPO" }] };
let currentDifferentialState = null;
let currentDifferentialGene = "";
let currentDifferentialPopulation = "";
let currentDifferentialInteraction = null;
let currentDifferentialFeatureRole = "ligand";
let differentialCy = null;
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
];
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
  if (window.cytoscape && window.cytoscape.prototype && typeof window.cytoscape.prototype.svg === "function") {
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
          if (window.cytoscape && window.cytoscape.prototype && typeof window.cytoscape.prototype.svg === "function") {
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

async function saveSvgMarkupAsPdf(svgMarkup, filename) {
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
  const availableHeight = pageHeight - margin * 2;
  const scale = Math.min(availableWidth / width, availableHeight / height, 1);
  const renderWidth = width * scale;
  const renderHeight = height * scale;
  const x = (pageWidth - renderWidth) / 2;
  const y = (pageHeight - renderHeight) / 2;
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
  const plot = document.getElementById(panelPlotId(panelKey));
  if (!plot) {
    return;
  }
  const cy = expressionCyByPanel[panelKey];
  if (cy) {
    cy.destroy();
    expressionCyByPanel[panelKey] = null;
  }
  try {
    Plotly.purge(plot);
  } catch (_) {
    // Ignore Plotly cleanup errors when the plot is not initialized.
  }
  plot.innerHTML = "";
  plot.style.height = "";
  plot.style.minHeight = "";
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
  if (raw === "lipid") {
    return "lipids";
  }
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
  return availableModalities();
}

function panelModality(panelKey) {
  const mode = getPanelSelectValue(panelKey, "mode");
  if (mode === "fastcomm_network") {
    return "rna";
  }
  const select = document.getElementById(panelElementId(panelKey, "modality"));
  if (!select || select.classList.contains("hidden")) {
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

function availableVisualizationModes(panelKey) {
  const modality = panelModality(panelKey);
  const modalityInfo = modalityDefinition(modality);
  const modes = [...BASE_VISUALIZATION_MODES];
  if (markerHeatmapAvailable(panelKey)) {
    modes.push({ value: "marker_heatmap", label: "MarkerHeatmap" });
  }
  const networkPopulations = markerNetworkPopulations(panelKey);
  if (modalityInfo.supports_marker_network !== false && panelMarkerAnalysis(panelKey) && networkPopulations.length) {
    modes.push({ value: "marker_network", label: "MarkerNetwork" });
  }
  if (modality === "rna" && fastCommAvailable()) {
    modes.push({ value: "fastcomm_network", label: "Cell communication" });
  }
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
      const currentModality = normalizeModalityId(modalitySelect.value || ((currentModalitiesState || {}).default || "rna"));
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
      modalityField.classList.toggle("hidden", modalities.length <= 1);
    }

    const currentPopulation = markerPopulationSelect.value;
    markerPopulationSelect.innerHTML = "";
    const fastcommPlotType = panelCommunicationPlotType(panelKey);
    const dropdownPopulations = modeSelect.value === "fastcomm_network" ? fastcommPopulations : networkPopulations;
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
      (mode === "marker_network" && networkPopulations.length > 0)
      || (mode === "fastcomm_network" && fastcommPopulations.length > 0 && fastCommPlotNeedsPopulation(fastcommPlotType));
    const showGene = mode === "expression_umap" || mode === "violin";
    updatePanelFeatureInput(panelKey);
    const populationLabel = markerPopulationField.querySelector("span");
    if (populationLabel) {
      populationLabel.textContent = mode === "fastcomm_network" ? "Marker cell state" : "Marker cell state";
    }
    markerPopulationField.classList.toggle("hidden", !showMarkerPopulation);
    geneField.classList.toggle("hidden", !showGene);
  });
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
  const plot = document.getElementById(panelPlotId(panelKey));
  const params = getDisplayFilterParams(panelKey);
  params.set("modality", panelModality(panelKey));
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
    if (!resp.ok) {
      throw new Error(`Marker heatmap TSV returned ${resp.status}.`);
    }
  } catch (err) {
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
    await logClientEvent(jobId, `MarkerHeatmap viewer preflight failed: ${err.message || err}`);
    renderVisualizationMessage(
      panelKey,
      `Marker heatmap viewer preflight failed.\n${err.message || err}`,
      "MarkerHeatmap"
    );
    return;
  }
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

function renderExpressionNetwork(panelKey, payload) {
  resetVisualizationSurface(panelKey);
  const plot = document.getElementById(panelPlotId(panelKey));
  const elements = (payload.elements || []).map((element) => {
    if (!element.data || !element.data.id || element.data.source) {
      return element;
    }
    const log2fc = Number(element.data.log2fc || 0);
    return {
      data: {
        ...element.data,
        color: log2fc >= 0 ? "#fca5a5" : "#7dd3fc",
      },
    };
  });
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
          width: 1.8,
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
  expressionCyByPanel[panelKey].on("tap", "node", (event) => {
    const gene = event?.target?.data("id");
    if (gene) {
      setPanelGeneValue(panelKey, gene);
    }
  });
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
});

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
    if (option.value === "none") {
      option.hidden = false;
      option.disabled = false;
      return;
    }
    const enabled = supported.includes(normalizeModalityId(option.value, ""));
    option.hidden = !enabled;
    option.disabled = !enabled;
  });
  if (!supported.includes(normalizeModalityId(select.value, ""))) {
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
  document.getElementById("differential-modality").addEventListener("change", () => {
    if (!currentDifferentialState) {
      return;
    }
    currentDifferentialState = {
      ...markDifferentialConfigDirty({
        ...currentDifferentialState,
        config: {
          ...(currentDifferentialState.config || {}),
          modality: normalizeModalityId(document.getElementById("differential-modality").value || "rna"),
        },
      }),
    };
    updateDifferentialUi(currentDifferentialState);
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
  initToggleMultiSelect(document.getElementById("differential-group1"));
  initToggleMultiSelect(document.getElementById("differential-group2"));
  attachDifferentialLrToggleHandlers();
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
    impute_modality: evt.target.impute_modality ? evt.target.impute_modality.value : "none",
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

function applyJobStatus(jobId, data) {
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
    exploreTabButton.classList.toggle("hidden", !hasExploreReady);
    exploreTabButton.disabled = !hasExploreReady;
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
    imputed_adt_h5ad: "Download imputed_adt_h5ad",
    adt_marker_genes_zip: "Download ADT marker ZIP",
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
    detailTitle.textContent = featureLabel === "lipid" ? "Lipid Detail" : "Gene Detail";
  }

  let statusMessage = state.message || "";
  if (state.status === "completed" && !state.go_terms_included && (state.visualization_modes || []).some((entry) => entry.value === "go")) {
    statusMessage += " GO terms were not available for this run.";
  }
  message.textContent = statusMessage;

  const currentVizMode = vizModeSelect.value;
  const visualizationModes = (state.visualization_modes || []);
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
  const populationSelect = document.getElementById("differential-result-population");
  const mode = document.getElementById("differential-viz-mode").value;
  const populations = differentialPopulationsForMode(state, mode);
  const fallbackPopulation = state.default_result_population || "";
  const wantedPopulation = currentDifferentialPopulation || populationSelect.value || fallbackPopulation;
  populationSelect.innerHTML = "";
  populations.forEach((population) => {
    const option = document.createElement("option");
    option.value = population;
    option.textContent = population;
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

async function loadDifferentialVisualization() {
  const state = currentDifferentialState;
  const plotEmpty = document.getElementById("differential-plot-empty");
  const jobId = document.getElementById("results-job-id").value.trim();
  if (!state || state.status !== "completed") {
    resetDifferentialResults();
    return;
  }

  const population = document.getElementById("differential-result-population").value;
  const mode = document.getElementById("differential-viz-mode").value;
  if (!population) {
    renderDifferentialEmpty(`No ${mode} data are available for this differential run.`);
    resetDifferentialGeneDetail();
    return;
  }

  currentDifferentialPopulation = population;
  updateDifferentialDownloadButton();
  try {
    let payload = null;
    if (mode === "heatmap") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/heatmap?population=${encodeURIComponent(population)}`));
      renderDifferentialHeatmap(payload);
    } else if (mode === "volcano") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/volcano?population=${encodeURIComponent(population)}`));
      renderDifferentialVolcano(payload);
    } else if (mode === "network") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/network?population=${encodeURIComponent(population)}`));
      renderDifferentialNetwork(payload);
    } else if (mode === "table") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/table?population=${encodeURIComponent(population)}`));
      renderDifferentialCommunicationTable(payload);
    } else if (mode === "go") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/go?population=${encodeURIComponent(population)}`));
      renderDifferentialGo(payload);
    }
    const payloadDefaultGene = (payload && payload.default_gene) || "";
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

function renderDifferentialHeatmap(payload) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  const rows = payload.rows || [];
  if (!rows.length) {
    renderDifferentialEmpty(`No heatmap rows were found for ${payload.population}.`);
    return;
  }
  const z = rows.map((row) => row.values);
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
  const y = rows.map((row) => row.gene);
  const directions = rows.map((row) => payload.columns.map(() => row.direction));
  const figureHeight = Math.max(560, Math.min(2800, rows.length * 18 + 180));
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  Plotly.newPlot(
    plot,
    [
      {
        z,
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
        hovertemplate: "%{y}<br>%{x}<br>log2FC=%{z:.3f}<br>%{customdata}<extra></extra>",
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
      yaxis: { automargin: true, tickfont: { size: 10 } },
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

function renderDifferentialVolcano(payload) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  const points = payload.points || [];
  if (!points.length) {
    renderDifferentialEmpty(`No volcano data were found for ${payload.population}.`);
    return;
  }
  const up = points.filter((point) => point.direction === "up");
  const down = points.filter((point) => point.direction === "down");
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  Plotly.newPlot(
    plot,
    [
      {
        x: down.map((point) => point.log2fc),
        y: down.map((point) => point.score),
        text: down.map((point) => point.gene),
        customdata: down.map((point) => [point.gene, point.fdr, point.pval]),
        type: "scattergl",
        mode: "markers",
        name: "Down",
        marker: { color: "#2563eb", size: 16, opacity: 0.72 },
        hovertemplate: "%{text}<br>log2FC=%{x:.3f}<br>-log10(FDR)=%{y:.3f}<extra></extra>",
      },
      {
        x: up.map((point) => point.log2fc),
        y: up.map((point) => point.score),
        text: up.map((point) => point.gene),
        customdata: up.map((point) => [point.gene, point.fdr, point.pval]),
        type: "scattergl",
        mode: "markers",
        name: "Up",
        marker: { color: "#dc2626", size: 16, opacity: 0.72 },
        hovertemplate: "%{text}<br>log2FC=%{x:.3f}<br>-log10(FDR)=%{y:.3f}<extra></extra>",
      },
    ],
    {
      title: `Volcano: ${payload.population}`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      margin: { t: 56, l: 60, r: 20, b: 56 },
      height: 640,
      xaxis: { title: "log2 fold change", zeroline: true, zerolinecolor: "rgba(100,116,139,0.45)" },
      yaxis: { title: "-log10(FDR)" },
      hovermode: "closest",
    },
    { responsive: true }
  );
  plot.on("plotly_click", (event) => {
    const gene = event?.points?.[0]?.text;
    if (gene) {
      currentDifferentialGene = gene;
      loadDifferentialGeneDetail(gene, payload.population);
    }
  });
}

function renderDifferentialGo(payload) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  const terms = payload.terms || [];
  if (!terms.length) {
    renderDifferentialEmpty(`No GO term enrichment was available for ${payload.population}.`);
    return;
  }
  const ordered = [...terms]
    .map((term) => {
      const overlapGenes = term.overlap_genes || [];
      const selectedGene = overlapGenes.includes(currentDifferentialGene)
        ? currentDifferentialGene
        : (term.selected_gene || overlapGenes[0] || null);
      return {
        ...term,
        selected_gene: selectedGene,
      };
    })
    .reverse();
  const yKeys = ordered.map((term, index) => `${term.term_name}__${term.direction || "unknown"}__${index}`);
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  Plotly.newPlot(
    plot,
    [
      {
        x: ordered.map((term) => term.score),
        y: yKeys,
        type: "bar",
        orientation: "h",
        customdata: ordered.map((term) => [term.selected_gene, term.direction, (term.overlap_genes || []).join(", "), term.overlap_genes || []]),
        marker: {
          color: ordered.map((term) => (term.direction === "up" ? "#dc2626" : "#2563eb")),
        },
        hovertemplate: (
          "%{y}<br>-log10(FDR)=%{x:.3f}"
          + "<br>Selected gene: %{customdata[0]}"
          + "<br>GO genes: %{customdata[2]}<extra></extra>"
        ),
      },
    ],
    {
      title: `GO terms: ${payload.population}`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      margin: { t: 56, l: 260, r: 20, b: 56 },
      height: Math.max(540, ordered.length * 34 + 120),
      xaxis: { title: "-log10(FDR)" },
      yaxis: {
        automargin: true,
        tickmode: "array",
        tickvals: yKeys,
        ticktext: ordered.map((term) => term.term_name),
      },
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

function renderDifferentialNetwork(payload) {
  if (payload && payload.network_type === "cell_communication_diff") {
    renderDifferentialCommunicationNetwork(payload);
    return;
  }
  const plot = document.getElementById("differential-plot-area");
  Plotly.purge(plot);
  destroyDifferentialNetwork();
  const elements = (payload.elements || []).map((element) => {
    if (!element.data || !element.data.id || element.data.source) {
      return element;
    }
    const log2fc = Number(element.data.log2fc || 0);
    return {
      data: {
        ...element.data,
        color: log2fc >= 0 ? "#fca5a5" : "#7dd3fc",
      },
    };
  });
  if (!elements.length) {
    renderDifferentialEmpty(`No interaction network was available for ${payload.population}.`);
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
          width: 1.8,
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

function renderDifferentialCommunicationNetwork(payload) {
  const plot = document.getElementById("differential-plot-area");
  Plotly.purge(plot);
  destroyDifferentialNetwork();
  const elements = payload.elements || [];
  if (!elements.length) {
    renderDifferentialEmpty(`No differential cell-state network was available for ${payload.population}.`);
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

function renderDifferentialCommunicationTable(payload) {
  destroyDifferentialNetwork();
  const plot = document.getElementById("differential-plot-area");
  const rows = payload.rows || [];
  const columns = payload.columns || [];
  if (!rows.length || !columns.length) {
    renderDifferentialEmpty(`No differential interaction table was available for ${payload.population}.`);
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
  const jobId = document.getElementById("results-job-id").value.trim();
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
    resetDifferentialGeneDetail();
    document.getElementById("differential-selected-gene").textContent = gene;
    document.getElementById("differential-gene-empty").textContent = err.message || "Unable to load gene detail.";
    document.getElementById("differential-gene-empty").classList.remove("hidden");
  }
}

function renderDifferentialGeneDetail(payload) {
  const plot = document.getElementById("differential-gene-plot");
  const empty = document.getElementById("differential-gene-empty");
  const stats = document.getElementById("differential-gene-stats");
  const downloadButton = document.getElementById("download-differential-gene-btn");
  const isCommunication = normalizeModalityId(payload.modality) === "cell_communication";
  const isFeatureExpression = payload.view_kind === "feature_expression";
  const yTitle = isFeatureExpression
    ? "Normalized expression"
    : (isCommunication ? "fastComm communication score" : "Normalized expression");
  const headerLabel = isFeatureExpression
    ? `${payload.gene} (${payload.feature_role || "feature"}) in ${payload.feature_state || payload.population}`
    : `${payload.gene} in ${payload.population}`;
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
    hovertemplate: `${group.label}<br>${isFeatureExpression ? "expr" : (isCommunication ? "score" : "expr")}=%{y:.3f}<extra></extra>`,
  }));
  const titleSubject = isFeatureExpression
    ? `${payload.gene} in ${payload.feature_state || payload.population}`
    : payload.gene;
  Plotly.newPlot(
    plot,
    traces,
    {
      title: `${yTitle}: ${titleSubject}`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      margin: { t: 52, l: 50, r: 18, b: 48 },
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

function updateDifferentialLrToggle(payload) {
  const toggle = document.getElementById("differential-lr-toggle");
  if (!toggle) {
    return;
  }
  const ligand = String((payload && payload.ligand) || "").trim();
  const receptor = String((payload && payload.receptor) || "").trim();
  if (!ligand && !receptor) {
    toggle.classList.add("hidden");
    return;
  }
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
  const population = document.getElementById("differential-result-population").value;
  if (!population) {
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
  resetDifferentialGeneDetail();
}

function resetDifferentialGeneDetail() {
  const featureLabel = String((currentDifferentialState && currentDifferentialState.feature_label) || "gene");
  const capitalizedFeature = featureLabel.charAt(0).toUpperCase() + featureLabel.slice(1);
  const plot = document.getElementById("differential-gene-plot");
  Plotly.purge(plot);
  plot.classList.add("hidden");
  document.getElementById("differential-selected-gene").textContent = `Select a ${featureLabel} from the left panel.`;
  document.getElementById("download-differential-gene-btn").classList.add("hidden");
  document.getElementById("download-differential-gene-btn").removeAttribute("href");
  const empty = document.getElementById("differential-gene-empty");
  empty.textContent = `Select a ${featureLabel} from the differential view to compare ${featureLabel} values between groups.`;
  empty.classList.remove("hidden");
  document.getElementById("differential-gene-stats").classList.add("hidden");
  const detailTitle = document.getElementById("differential-detail-title");
  if (detailTitle) {
    detailTitle.textContent = `${capitalizedFeature} Detail`;
  }
  const lrToggle = document.getElementById("differential-lr-toggle");
  if (lrToggle) {
    lrToggle.classList.add("hidden");
  }
  currentDifferentialInteraction = null;
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
      const features = [];
      (data.genes || []).forEach((gene) => {
        const value = String(gene || "").trim();
        if (!value || seen.has(value)) {
          return;
        }
        seen.add(value);
        features.push(value);
        const option = document.createElement("option");
        option.value = value;
        datalist.appendChild(option);
      });
      const input = document.getElementById(panelElementId(panelKey, "gene-query"));
      if (input) {
        input.setAttribute("list", `${panelKey}-feature-suggestions`);
        const currentValue = String(input.value || "").trim();
        if (!currentValue || !seen.has(currentValue)) {
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
  try {
    if (isUmapMode(mode)) {
      const params = getDisplayFilterParams(panelKey);
      params.set("modality", modality);
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
  Object.assign(layout, buildSquareUmapAxes(umapData.query, 0.06));
  Plotly.newPlot(panelPlotId(panelKey), traces, layout);
}

function renderPanelExpression(panelKey, expressionData, mode, dotScale) {
  updateBaselineFilterSummaries(panelKey);
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
        title: `${expressionData.gene} (reference centroid expression, top 10 states)`,
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
      title: `${expressionData.gene} (top 10 states by mean)`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.9)",
      height: 560,
      margin: { t: 36, l: 48, r: 20, b: 120 },
      yaxis: { title: "Expression", range: [globalMin - yPad, globalMax + yPad] },
    });
    return;
  }
  const zeroPoints = expressionData.umap.filter((p) => Number(p.value) <= 0);
  const expressedPoints = expressionData.umap
    .filter((p) => Number(p.value) > 0)
    .sort((a, b) => Number(a.value) - Number(b.value));
  const traces = [];
  if (zeroPoints.length) {
    traces.push({
      x: zeroPoints.map((p) => p.x),
      y: zeroPoints.map((p) => p.y),
      text: zeroPoints.map((p) => `${p.barcode}<br>${p.population}<br>${expressionData.gene}: 0.000`),
      mode: "markers",
      type: "scattergl",
      marker: {
        size: 2 * dotScale,
        color: "#e5e7eb",
        opacity: 0.9,
      },
      hoverinfo: "text",
      name: `${expressionData.gene} = 0`,
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
    const useImputedPalette = normalizedModality === "lipids" || normalizedModality === "adt";
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
        text: bin.map((p) => `${p.barcode}<br>${p.population}<br>${expressionData.gene}: ${p.value.toFixed(3)}`),
        mode: "markers",
        type: "scattergl",
        marker: {
          size: 2 * dotScale,
          color: bin.map((p) => p.value),
          colorscale,
          cmin: minValue,
          cmax: maxValue,
          showscale: start + step >= expressedPoints.length,
          colorbar: start + step >= expressedPoints.length ? { title: expressionData.gene } : undefined,
        },
        name: expressionData.gene,
        showlegend: false,
      });
    }
  }
  const expressionLayout = {
    title: expressionData.source === "reference"
      ? `${expressionData.gene} reference expression`
      : `${expressionData.gene} expression`,
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
  const data = panelPlotData[panelKey];
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
  if (!jobId || !mode || !panelPlotData[panelKey]) {
    return;
  }
  if (mode === "marker_heatmap") {
    const params = getDisplayFilterParams(panelKey);
    params.set("modality", modality);
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
        const svg = cy.svg({ scale: 1, full: true });
        const blob = new Blob([svg], { type: "image/svg+xml;charset=utf-8" });
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = buildPdfFilename([jobId, panelKey, plotType], "svg").replace(/\.pdf$/i, ".svg");
        link.click();
        URL.revokeObjectURL(url);
      } catch (err) {
        console.warn(err);
      }
      return;
    }
    await exportPlotlyElementToVectorPdf(
      panelPlotId(panelKey),
      buildPdfFilename([jobId, panelKey, plotType], plotType),
    );
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
  if (!population) {
    return;
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
  if (mode === "network") {
    window.open(
      apiPath(`/jobs/${jobId}/differential/interactive/pdf?mode=${encodeURIComponent(mode)}&population=${encodeURIComponent(population)}`),
      "_blank",
    );
    return;
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
  activeExplorerTab = tab;
  document.querySelectorAll(".workspace-tab-btn").forEach((button) => {
    button.classList.toggle("active", button.dataset.tab === tab);
  });
  document.querySelectorAll(".workspace-panel").forEach((panel) => {
    panel.classList.toggle("active", panel.dataset.tabPanel === tab);
  });
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
