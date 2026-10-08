"use strict";
// scALABLE-discover front end. Loaded after webapp/static/app.js, which stays unedited.
// app.js is a classic script, so its top-level functions are globals that a later script
// can wrap or replace, and its top-level const arrays can be edited in place.

const DISCOVER_LAYER_COOKIE = "discover_layer";
const GOELITE_MODE = "goelite_biomarkers";
let discoverStatus = null;
let discoverGoeliteStates = { key: "", states: [] };
let discoverGoeliteRequest = null;

(function clearDiscoverResults() {
  const original = explorePayloadCache.onClear;
  explorePayloadCache.onClear = () => {
    discoverStatus = null;
    discoverGoeliteStates = { key: "", states: [] };
    discoverGoeliteRequest = null;
    document.getElementById("discover-layer-field")?.remove();
    original();
  };
})();

window.scalableQcOptions = function discoverQcOptions(form) {
  const maxK = form.elements.max_k.value.trim();
  return {
    umap_fit_mode: form.elements.umap_fit_mode.value,
    max_k: maxK === "" ? null : Number(maxK),
  };
};

// Saved jobs retain their original category IDs for requests and filtering. Only
// the label shown on their UMAP is abbreviated; new jobs write c18 directly.
function discoverCellStateLabel(value) {
  return String(value ?? "").replace(/^UNK[-_]+c?(\d+)$/i, "c$1");
}

(function cleanSavedUmapLabels() {
  const original = renderPanelUmap;
  renderPanelUmap = function discoverRenderPanelUmap(panelKey, data, mode, dotScale) {
    const needsCleanup = ["query", "reference"].some((key) =>
      (data[key] || []).some((point) => discoverCellStateLabel(point.population) !== point.population));
    if (!needsCleanup) return original.apply(this, arguments);
    const shown = { ...data };
    for (const key of ["query", "reference"]) {
      shown[key] = (data[key] || []).map((point) => {
        const population = discoverCellStateLabel(point.population);
        return population === point.population ? point : { ...point, population };
      });
    }
    for (const key of ["population_order", "populations", "lineage_order"]) {
      if (Array.isArray(data[key])) shown[key] = data[key].map(discoverCellStateLabel);
    }
    return original.call(this, panelKey, shown, mode, dotScale);
  };
})();

// Keep each option's original value so filtering saved jobs still addresses the
// original category. Abbreviate only the text that visitors see.
(function cleanSavedMenuLabels() {
  const original = populateSelectOptions;
  populateSelectOptions = function discoverPopulateSelectOptions(selectEl) {
    const result = original.apply(this, arguments);
    for (const option of selectEl.options) option.textContent = discoverCellStateLabel(option.textContent);
    return result;
  };
  const originalFilter = syncDisplayFilterValueOptions;
  syncDisplayFilterValueOptions = function discoverSyncDisplayFilterValueOptions(panelKey, index) {
    const result = originalFilter.apply(this, arguments);
    const select = document.getElementById(panelElementId(panelKey, `filter${index}-values`));
    for (const option of select.options) option.textContent = discoverCellStateLabel(option.textContent);
    return result;
  };
})();

// No reference atlas: the "UMAP broad" view draws the query over reference cells.
(function adjustPlotTypes() {
  const broad = BASE_VISUALIZATION_MODES.findIndex((mode) => mode.value === "relative");
  if (broad >= 0) BASE_VISUALIZATION_MODES.splice(broad, 1);
  const states = BASE_VISUALIZATION_MODES.find((mode) => mode.value === "cluster");
  if (states) states.label = "UMAP cell states";
})();

// "Regulatory network" draws TF-to-target edges that only GRN imputation produces, which this
// version does not run. This app has one modality, so its pathway view is plain "Pathway".
// The GO-Elite BioMarkers plot joins the list when the job holds ICGS3's BioMarkers table.
(function adjustModeList() {
  const original = availableVisualizationModes;
  availableVisualizationModes = function discoverVisualizationModes(panelKey) {
    const modes = original.apply(this, arguments)
      .filter((mode) => mode.value !== "integrated_network")
      .map((mode) => (mode.value === "integrated_cross_pathway" ? { ...mode, label: "Pathway" } : mode));
    if (discoverStatus?.icgs3_analysis?.goelite?.available) {
      const after = modes.findIndex((mode) => mode.value === "marker_network");
      modes.splice(after >= 0 ? after + 1 : modes.length, 0, { value: GOELITE_MODE, label: "GO-Elite BioMarkers" });
    }
    return modes;
  };
})();

// ---------------------------------------------------------------- cell-state layers

function discoverCookiePath() {
  return APP_ROOT_PATH || "/";
}

function renderLayerControl(data) {
  const layers = data?.cell_state_layers?.layers || [];
  const form = document.getElementById("results-form");
  if (!form || layers.length < 2) return;
  let field = document.getElementById("discover-layer-field");
  if (!field) {
    field = document.createElement("label");
    field.className = "field";
    field.id = "discover-layer-field";
    field.innerHTML = '<span>Cell-state layer</span><select id="discover-layer-select"></select>';
    form.insertBefore(field, form.querySelector("label.field"));
    field.querySelector("select").addEventListener("change", (event) => {
      // The server reads the layer from this cookie on every request, so a reload redraws
      // every panel, menu and Chat example from that layer.
      document.cookie = `${DISCOVER_LAYER_COOKIE}=${encodeURIComponent(event.target.value)}; path=${discoverCookiePath()}; SameSite=Lax`;
      window.location.reload();
    });
  }
  const select = field.querySelector("select");
  select.innerHTML = "";
  layers.forEach((layer) => {
    const option = document.createElement("option");
    option.value = layer.key;
    option.textContent = layer.label || layer.key;
    select.appendChild(option);
  });
  select.value = data.active_cell_state_layer || data.cell_state_layers.default || layers[0].key;
}

(function captureStatus() {
  const original = applyJobStatus;
  applyJobStatus = function discoverApplyJobStatus(jobId, data) {
    const result = original.apply(this, arguments);
    discoverStatus = data;
    renderLayerControl(data);
    updateExpressionModeOptions();
    void loadGoeliteStates(jobId, data);
    return result;
  };
})();

// ------------------------------------------------------------- GO-Elite BioMarkers plot

async function loadGoeliteStates(jobId, data) {
  if (!jobId || !data?.icgs3_analysis?.goelite?.available || data.status !== "completed") return;
  const generation = explorePayloadCache.generation;
  const key = `${generation}:${jobId}:${data.active_cell_state_layer || ""}`;
  if (discoverGoeliteStates.key === key) return;
  if (discoverGoeliteRequest?.key === key) return discoverGoeliteRequest.promise;
  const isCurrent = () => generation === explorePayloadCache.generation && jobId === getResultsJobId();
  discoverGoeliteStates = { key: "", states: [] };
  const request = { key, promise: null };
  request.promise = (async () => {
    try {
      const payload = await exploreMetadataCache.fetch(apiPath(`/jobs/${encodeURIComponent(jobId)}/biomarkers/goelite/states`));
      if (!isCurrent()) return;
      discoverGoeliteStates = { key, states: payload.states || [] };
      updateExpressionModeOptions();
    } catch (error) {
      if (isCurrent()) console.warn("GO-Elite BioMarkers states unavailable:", error);
    } finally {
      if (discoverGoeliteRequest === request) discoverGoeliteRequest = null;
    }
  })();
  discoverGoeliteRequest = request;
  return request.promise;
}

(function goeliteMenus() {
  const original = updateExpressionModeOptions;
  updateExpressionModeOptions = function discoverUpdateExpressionModeOptions() {
    const result = original.apply(this, arguments);
    VISUALIZATION_PANELS.forEach((panelKey) => {
      if (getPanelSelectValue(panelKey, "mode") !== GOELITE_MODE) return;
      const field = document.getElementById(panelElementId(panelKey, "marker-population-field"));
      const select = document.getElementById(panelElementId(panelKey, "marker-population"));
      if (!field || !select) return;
      const previous = select.value;
      select.innerHTML = "";
      discoverGoeliteStates.states.forEach((state) => {
        const option = document.createElement("option");
        option.value = state;
        option.textContent = state;
        select.appendChild(option);
      });
      if (discoverGoeliteStates.states.includes(previous)) select.value = previous;
      field.classList.toggle("hidden", !discoverGoeliteStates.states.length);
      const label = field.querySelector("span");
      if (label) label.textContent = "Cell state";
      ["gene-field", "modality-field", "filter-stack"].forEach((suffix) => {
        document.getElementById(panelElementId(panelKey, suffix))?.classList.add("hidden");
      });
    });
    return result;
  };
})();


(function goelitePanel() {
  const originalLoad = loadVisualizationPanel;
  loadVisualizationPanel = async function discoverLoadVisualizationPanel(panelKey) {
    if (getPanelSelectValue(panelKey, "mode") !== GOELITE_MODE) return originalLoad.apply(this, arguments);
    const jobId = getResultsJobId();
    const state = getPanelSelectValue(panelKey, "marker-population");
    const requestId = ++panelVisualizationRequest[panelKey];
    const generation = explorePayloadCache.generation;
    const isCurrent = () => requestId === panelVisualizationRequest[panelKey]
      && generation === explorePayloadCache.generation && jobId === getResultsJobId()
      && getPanelSelectValue(panelKey, "mode") === GOELITE_MODE
      && state === getPanelSelectValue(panelKey, "marker-population");
    resetVisualizationSurface(panelKey);
    if (!jobId || !state) {
      renderVisualizationMessage(panelKey, "Choose a cell state.", "GO-Elite BioMarkers");
      return;
    }
    try {
      const payload = await explorePayloadCache.fetch(apiPath(`/jobs/${encodeURIComponent(jobId)}/biomarkers/goelite?population=${encodeURIComponent(state)}`));
      if (!isCurrent()) return;
      if (!(payload.terms || []).length) {
        renderVisualizationMessage(panelKey, payload.message || payload.detail || "No BioMarkers terms.", "GO-Elite BioMarkers");
        return;
      }
      panelPlotData[panelKey] = { source: GOELITE_MODE, payload };
      renderGoeliteBiomarkers(panelKey, payload);
      const passing = payload.terms.filter((term) => term.is_selected_positive_sig).length;
      setPanelSummary(panelKey, `${payload.terms.length} BioMarkers terms overlap the markers of ${payload.population}; `
        + `${passing} pass FDR <= 0.05 and z > 2. The cell-state label comes from "${payload.terms[0].term_name}".`);
    } catch (error) {
      if (isCurrent()) renderVisualizationMessage(panelKey, error.message || "GO-Elite BioMarkers plot unavailable.", "GO-Elite BioMarkers");
    }
  };
  const originalDownload = downloadVisualizationImage;
  downloadVisualizationImage = async function discoverDownloadVisualizationImage(panelKey) {
    if (getPanelSelectValue(panelKey, "mode") !== GOELITE_MODE) return originalDownload.apply(this, arguments);
    const jobId = getResultsJobId();
    const state = getPanelSelectValue(panelKey, "marker-population");
    if (jobId && state) {
      window.open(apiPath(`/jobs/${encodeURIComponent(jobId)}/biomarkers/goelite.pdf?population=${encodeURIComponent(state)}`), "_blank");
    }
  };
})();

// Changing Windows from 2 to 1 (and any other resize) redraws each panel from its cached data
// through renderVisualizationPanel, which has no branch for this plot type and fell through to
// the gene-expression view ("Gene '' was not found."). Redraw the cached BioMarkers payload.
(function goeliteRedraw() {
  const original = renderVisualizationPanel;
  renderVisualizationPanel = function discoverRenderVisualizationPanel(panelKey) {
    if (getPanelSelectValue(panelKey, "mode") !== GOELITE_MODE) return original.apply(this, arguments);
    const data = panelPlotData[panelKey];
    if (!data || data.source !== GOELITE_MODE) {
      void loadVisualizationPanel(panelKey);
      return;
    }
    resetVisualizationSurface(panelKey);
    renderGoeliteBiomarkers(panelKey, data.payload);
  };
})();

// ------------------------------------------------------------------- labels and panels

// The default embedding and the default colouring are ICGS3's; scALABLE-web names them after
// cellHarmony. Names here carry no tool prefix.
(function relabelUmapOptions() {
  const original = refreshUmapOptions;
  refreshUmapOptions = async function discoverRefreshUmapOptions(panelKey) {
    await original.apply(this, arguments);
    document.querySelectorAll(`#${panelKey}-colorby option, #${panelKey}-coords option`).forEach((option) => {
      option.textContent = option.textContent.replace(" (cellHarmony)", "").replace("cellHarmony UMAP", "UMAP");
    });
  };
})();

// The Run tab's preview panel described the chosen reference atlas. ICGS3 has none, so the
// panel shows the ICGS3 workflow figure instead. app.js calls this function by name.
function loadReferencePreview() {
  const plot = document.getElementById("reference-preview-plot");
  if (!plot || plot.dataset.discover === "1") return;
  plot.dataset.discover = "1";
  if (window.Plotly) Plotly.purge(plot);
  plot.innerHTML = `<img src="${APP_ROOT_PATH}/discover-static/icgs3_overview.png"
    alt="ICGS3: H5 files, downsampling and iterative discovery, ordered cell states"
    style="display:block;width:100%;height:100%;object-fit:contain;">`;
}

// One line under the QC counts: how many QC-retained cells ICGS3 placed in a cluster.
// ICGS3 drops a cell whose winning SVM score fails its threshold, so this denominator matters.
function discoverClusteringSummary(data) {
  const icgs3 = data && data.icgs3_analysis;
  if (!icgs3 || icgs3.status !== "completed") return "";
  const qc = Number(icgs3.qc_cells);
  const kept = Number(icgs3.clustered_cells);
  if (!Number.isFinite(qc) || !Number.isFinite(kept) || qc <= 0) return "";
  const percent = (100 * kept / qc).toFixed(1);
  return `ICGS3 assigned ${kept.toLocaleString()} of ${qc.toLocaleString()} QC-retained cells (${percent}%) ` +
    `to ${icgs3.n_clusters} cell states; ${(qc - kept).toLocaleString()} cells received none and are not shown.`;
}

// Every path that writes the status line goes through buildQcCellSummary. While the job runs,
// the line shows the runner's current stage (QC, then each ICGS3 step) with the QC count;
// scALABLE-web's own log scan would report its last alignment-era marker instead.
(function wrapQcCellSummary() {
  const original = buildQcCellSummary;
  buildQcCellSummary = function discoverQcCellSummary(data) {
    const status = String(data && data.status || "").toLowerCase();
    if (status === "processing" && data.message) {
      const qc = extractQcThresholdState(getStatusLogLines(data));
      const parts = [String(data.message)];
      if (Number.isFinite(qc.afterMito)) parts.push(`QC-retained cells: ${qc.afterMito.toLocaleString()}`);
      return parts.join(" | ");
    }
    const text = String(original.apply(this, arguments) || "");
    const summary = status === "completed" ? discoverClusteringSummary(data) : "";
    return summary ? `${text.trim()} ${summary}`.trim() : text;
  };
})();

const DISCOVER_DOWNLOAD_LABELS = {
  assignments: "Download cluster assignments",
  marker_genes_zip: "Download MarkerFinder ZIP",
  icgs3_results_zip: "Download ICGS3 results ZIP",
  parameters_json: "Download run parameters",
};

(function wrapDownloadLinks() {
  const original = populateDownloadLinks;
  populateDownloadLinks = async function discoverPopulateDownloadLinks(jobId, statusData) {
    await original.apply(this, arguments);
    document.querySelectorAll("#download-links a.download-btn").forEach((link) => {
      const key = decodeURIComponent(String(link.getAttribute("href") || "").split("/download/")[1] || "");
      if (DISCOVER_DOWNLOAD_LABELS[key]) link.textContent = DISCOVER_DOWNLOAD_LABELS[key];
    });
  };
})();
