"use strict";

const MAX_SAMPLES = 7;
let registry = window.__REFERENCE_REGISTRY__ || { species: [] };
let sampleCount = 0;
let pollTimer = null;
let currentUmapData = null;
let currentExpressionData = null;
let loadedResultsJobId = null;
let currentDifferentialState = null;

document.addEventListener("DOMContentLoaded", () => {
  initSpeciesSelect();
  initSampleRows();
  hookForms();
  updateDifferentialUi(null);
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
      return;
    }
    selected.references.forEach((ref) => {
      const option = document.createElement("option");
      option.value = ref.id;
      option.textContent = ref.label;
      referenceSelect.appendChild(option);
    });
  });

  if (registry.species.length) {
    speciesSelect.value = registry.species[0].id;
    speciesSelect.dispatchEvent(new Event("change"));
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
    <label class="field">
      <span>Sample name</span>
      <input type="text" class="sample-name" placeholder="e.g. Sample_${sampleCount + 1}" required>
    </label>
    <label class="field">
      <span>H5/H5AD file</span>
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
  document.getElementById("results-job-id").addEventListener("change", handleResultsJobChange);
  document.getElementById("gene-query").addEventListener("change", () => refreshResults());
  document.getElementById("expression-mode").addEventListener("change", renderCurrentExpression);
  document.getElementById("umap-mode").addEventListener("change", renderCurrentUmap);
  document.getElementById("download-umap-image-btn").addEventListener("click", downloadUmapImage);
  document.getElementById("download-expression-image-btn").addEventListener("click", downloadExpressionImage);
  document.getElementById("differential-network-select").addEventListener("change", renderSelectedNetwork);
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

  try {
    const resp = await fetch("/api/jobs", { method: "POST", body: formData });
    const data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "Failed to create job.");
    }
    const notice = document.getElementById("job-created");
    notice.classList.remove("hidden");
    notice.textContent = `Job created: ${data.job_id}`;
    document.getElementById("qc-job-id").value = data.job_id;
    document.getElementById("results-job-id").value = data.job_id;
    loadedResultsJobId = null;
    setResultMode("baseline");
    await loadJobState(data.job_id);
  } catch (err) {
    alert(err.message);
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
  };
  try {
    let resp = await fetch(`/api/jobs/${jobId}/qc`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    let data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "Failed to save QC.");
    }
    resp = await fetch(`/api/jobs/${jobId}/run`, { method: "POST" });
    data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "Failed to queue job.");
    }
    loadedResultsJobId = null;
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
    population_col: document.getElementById("differential-population").value,
    group1_samples: getMultiSelectValues(document.getElementById("differential-group1")),
    group2_samples: getMultiSelectValues(document.getElementById("differential-group2")),
  };

  try {
    const resp = await fetch(`/api/jobs/${jobId}/differential`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const data = await parseApiResponse(resp);
    if (!resp.ok) {
      throw new Error(data.detail || "Failed to start differential analysis.");
    }
    updateDifferentialUi(data);
    setResultMode("baseline");
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
  await loadJobState(jobId);
}

async function loadJobState(jobId) {
  try {
    const resp = await fetch(`/api/jobs/${jobId}/status`);
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
    const resp = await fetch(`/api/jobs/${jobId}/status`);
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
  document.getElementById("job-progress").style.width = `${data.progress || 0}%`;
  document.getElementById("job-progress-label").textContent = `${data.progress || 0}%`;
  document.getElementById("job-log").textContent = (data.qc_log_tail || []).join("");

  if (!document.getElementById("qc-job-id").value) {
    document.getElementById("qc-job-id").value = jobId;
  }
  if (!document.getElementById("results-job-id").value) {
    document.getElementById("results-job-id").value = jobId;
  }

  updateDifferentialUi(data.differential_ui || null);

  if (data.status === "completed") {
    populateDownloadLinks(jobId);
    if (!document.getElementById("gene-query").value && data.default_gene) {
      document.getElementById("gene-query").value = data.default_gene;
    }
    if (loadedResultsJobId !== jobId) {
      loadedResultsJobId = jobId;
      refreshResults();
    }
  }

  const pipelineActive = data.status === "queued" || data.status === "processing";
  const differentialStatus = (data.differential_ui || {}).status;
  const differentialActive = differentialStatus === "queued" || differentialStatus === "processing";
  if (!pipelineActive && !differentialActive && pollTimer) {
    clearInterval(pollTimer);
    pollTimer = null;
  }
}

async function populateDownloadLinks(jobId) {
  const container = document.getElementById("download-links");
  container.innerHTML = "";
  const resp = await fetch(`/api/jobs/${jobId}/status`);
  const data = await resp.json();
  if (!resp.ok) {
    return;
  }
  const artifacts = data.artifacts || {};
  Object.keys(artifacts).forEach((key) => {
    if (
      key === "umap_coordinates" ||
      key === "umap_placeholder_expression" ||
      key === "umap_pdf" ||
      key === "umap_pdf_plain"
    ) {
      return;
    }
    const btn = document.createElement("a");
    btn.className = "download-btn";
    btn.href = `/api/jobs/${jobId}/download/${key}`;
    btn.textContent = `Download ${key}`;
    container.appendChild(btn);
  });
}

function updateDifferentialUi(state) {
  currentDifferentialState = state;
  const panel = document.getElementById("differential-panel");
  const intro = document.getElementById("differential-intro");
  const populationSelect = document.getElementById("differential-population");
  const group1Select = document.getElementById("differential-group1");
  const group2Select = document.getElementById("differential-group2");
  const runBtn = document.getElementById("differential-run-btn");
  const progress = state ? state.progress || 0 : 0;
  const message = document.getElementById("differential-message");
  const archiveLink = document.getElementById("differential-archive-link");

  const enabled = Boolean(state && state.enabled);
  const config = (state && state.config) || {};
  const populationOptions = (state && state.population_columns) || [];
  const sampleNames = (state && state.sample_names) || [];
  const selectedPopulation = config.population_col || (state && state.default_population_col) || "";
  const showPanel = Boolean(state && state.enabled && populationOptions.length);

  panel.classList.toggle("hidden", !showPanel);

  populateSingleSelect(populationSelect, populationOptions, selectedPopulation);
  populateMultiSelect(group1Select, sampleNames, config.group1_samples || []);
  populateMultiSelect(group2Select, sampleNames, config.group2_samples || []);

  const disableInputs = !enabled || !populationOptions.length;
  populationSelect.disabled = disableInputs;
  group1Select.disabled = disableInputs;
  group2Select.disabled = disableInputs;
  runBtn.disabled = disableInputs;

  document.getElementById("differential-progress").style.width = `${progress}%`;
  document.getElementById("differential-progress-label").textContent = `${progress}%`;

  if (!state) {
    intro.textContent = "Available after alignment completes for jobs with two or more samples.";
    message.textContent = "Differential analysis is enabled when the job contains two or more samples.";
    archiveLink.classList.add("hidden");
    resetDifferentialResults();
    setResultMode("baseline");
    return;
  }

  if (!showPanel) {
    archiveLink.classList.add("hidden");
    resetDifferentialResults();
    setResultMode("baseline");
    return;
  }

  intro.textContent = enabled
    ? "Use the aligned cellHarmony AnnData to compare two sample groups across a selected cell-state field."
    : "Differential analysis is only enabled when two or more samples were uploaded for the job.";

  let statusMessage = state.message || "";
  if (state.status === "completed" && !state.go_terms_included) {
    statusMessage += " GO terms were not available for this run.";
  }
  message.textContent = statusMessage;

  if (state.archive_url) {
    archiveLink.href = state.archive_url;
    archiveLink.classList.remove("hidden");
  } else {
    archiveLink.classList.add("hidden");
  }

  if (state.status === "completed" && state.heatmap_svg_url) {
    renderDifferentialResults(state);
    setResultMode("differential");
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

function getMultiSelectValues(element) {
  return Array.from(element.selectedOptions).map((option) => option.value);
}

function renderDifferentialResults(state) {
  const heatmapImage = document.getElementById("differential-heatmap-image");
  const heatmapEmpty = document.getElementById("differential-heatmap-empty");
  const heatmapDownload = document.getElementById("download-differential-heatmap-btn");
  heatmapImage.src = `${state.heatmap_svg_url}&run=${encodeURIComponent(state.run_id || "latest")}`;
  heatmapImage.classList.remove("hidden");
  heatmapEmpty.classList.add("hidden");
  if (state.heatmap_pdf_url) {
    heatmapDownload.href = state.heatmap_pdf_url;
    heatmapDownload.classList.remove("hidden");
  } else {
    heatmapDownload.classList.add("hidden");
  }

  const networkSelect = document.getElementById("differential-network-select");
  const networkDownload = document.getElementById("download-differential-network-btn");
  const desiredNetwork = networkSelect.value;
  networkSelect.innerHTML = "";
  (state.networks || []).forEach((network) => {
    const option = document.createElement("option");
    option.value = network.id;
    option.textContent = network.population;
    if (network.id === desiredNetwork) {
      option.selected = true;
    }
    networkSelect.appendChild(option);
  });

  if ((state.networks || []).length) {
    networkSelect.classList.remove("hidden");
    networkDownload.classList.remove("hidden");
    if (!networkSelect.value) {
      networkSelect.selectedIndex = 0;
    }
    renderSelectedNetwork();
  } else {
    networkSelect.classList.add("hidden");
    networkDownload.classList.add("hidden");
    const image = document.getElementById("differential-network-image");
    const empty = document.getElementById("differential-network-empty");
    if (state.cell_frequency_grouped_png_url) {
      image.src = `${state.cell_frequency_grouped_png_url}?run=${encodeURIComponent(state.run_id || "latest")}`;
      image.classList.remove("hidden");
      empty.classList.remove("hidden");
      empty.textContent = "Interaction network preview was unavailable for this run. Showing the cell-frequency comparison plot instead.";
    } else {
      image.classList.add("hidden");
      empty.classList.remove("hidden");
      empty.textContent = "No interaction networks were generated for this comparison.";
    }
  }
}

function renderSelectedNetwork() {
  const state = currentDifferentialState;
  if (!state || !state.networks || !state.networks.length) {
    return;
  }
  const select = document.getElementById("differential-network-select");
  const networkDownload = document.getElementById("download-differential-network-btn");
  const selected = state.networks.find((network) => network.id === select.value) || state.networks[0];
  const image = document.getElementById("differential-network-image");
  const empty = document.getElementById("differential-network-empty");
  image.src = `${selected.png_url}&run=${encodeURIComponent(state.run_id || "latest")}`;
  image.classList.remove("hidden");
  empty.classList.add("hidden");
  networkDownload.href = selected.pdf_url;
  networkDownload.classList.remove("hidden");
}

function resetDifferentialResults() {
  const heatmapImage = document.getElementById("differential-heatmap-image");
  const networkImage = document.getElementById("differential-network-image");
  const heatmapDownload = document.getElementById("download-differential-heatmap-btn");
  const networkDownload = document.getElementById("download-differential-network-btn");
  heatmapImage.removeAttribute("src");
  networkImage.removeAttribute("src");
  heatmapImage.classList.add("hidden");
  networkImage.classList.add("hidden");
  heatmapDownload.classList.add("hidden");
  networkDownload.classList.add("hidden");
  document.getElementById("differential-heatmap-empty").classList.remove("hidden");
  document.getElementById("differential-network-empty").classList.remove("hidden");
  document.getElementById("differential-network-empty").textContent = "Interaction networks will appear here when the differential run finishes.";
  document.getElementById("differential-network-select").classList.add("hidden");
}

function setResultMode(mode) {
  const baseline = document.getElementById("baseline-results-view");
  const differential = document.getElementById("differential-results-view");
  if (mode === "differential") {
    baseline.classList.add("hidden");
    differential.classList.remove("hidden");
    return;
  }
  baseline.classList.remove("hidden");
  differential.classList.add("hidden");
}

async function loadUmap() {
  const jobId = document.getElementById("results-job-id").value.trim();
  if (!jobId) {
    return;
  }
  try {
    const resp = await fetch(`/api/jobs/${jobId}/umap`);
    const data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "UMAP not ready.");
    }
    currentUmapData = data;
    renderCurrentUmap();
  } catch (err) {
    alert(err.message);
  }
}

function renderCurrentUmap() {
  if (!currentUmapData) {
    return;
  }
  const mode = document.getElementById("umap-mode").value;
  const traces = [];
  if (mode === "relative") {
    traces.push({
      x: currentUmapData.reference.map((p) => p.x),
      y: currentUmapData.reference.map((p) => p.y),
      text: currentUmapData.reference.map((p) => `${p.barcode}<br>${p.population}`),
      mode: "markers",
      type: "scattergl",
      marker: { color: "#94a3b8", size: 2 },
      name: "Reference",
    });
    traces.push({
      x: currentUmapData.query.map((p) => p.x),
      y: currentUmapData.query.map((p) => p.y),
      text: currentUmapData.query.map((p) => `${p.barcode}<br>${p.population}`),
      mode: "markers",
      type: "scattergl",
      marker: { color: "#f97316", size: 2 },
      name: "Query",
    });
  } else {
    traces.push({
      x: currentUmapData.reference.map((p) => p.x),
      y: currentUmapData.reference.map((p) => p.y),
      text: currentUmapData.reference.map((p) => `${p.barcode}<br>${p.population}`),
      mode: "markers",
      type: "scattergl",
      marker: { color: "#d1d5db", size: 1, opacity: 0.45 },
      name: "Reference",
      showlegend: false,
    });
    const populations = [...new Set(currentUmapData.query.map((p) => p.population))];
    populations.forEach((population) => {
      const subset = currentUmapData.query.filter((p) => p.population === population);
      traces.push({
        x: subset.map((p) => p.x),
        y: subset.map((p) => p.y),
        text: subset.map((p) => `${p.barcode}<br>${p.population}`),
        mode: "markers",
        type: "scattergl",
        marker: { size: 2 },
        name: population,
      });
    });
  }
  Plotly.newPlot("umap-plot", traces, {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.9)",
    margin: { t: 20, l: 40, r: 20, b: 40 },
    legend: { orientation: "h" },
    hovermode: "closest",
  });
}

async function loadExpression() {
  const jobId = document.getElementById("results-job-id").value.trim();
  const gene = document.getElementById("gene-query").value.trim();
  if (!jobId || !gene) {
    return;
  }
  try {
    const resp = await fetch(`/api/jobs/${jobId}/expression?gene=${encodeURIComponent(gene)}`);
    const data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "Expression unavailable.");
    }
    currentExpressionData = data;
    renderCurrentExpression();
  } catch (err) {
    alert(err.message);
  }
}

function renderCurrentExpression() {
  if (!currentExpressionData) {
    return;
  }
  const mode = document.getElementById("expression-mode").value;
  if (mode === "violin") {
    const violinTraces = currentExpressionData.violin.map((entry) => ({
      type: "violin",
      name: entry.population,
      y: entry.values,
      box: { visible: false },
      meanline: { visible: true },
      points: "all",
      jitter: 0.18,
      pointpos: 0,
      marker: { size: 3, opacity: 0.45 },
    }));
    Plotly.newPlot("expression-plot", violinTraces, {
      title: `${currentExpressionData.gene} (top 10 states by mean)`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.9)",
      margin: { t: 36, l: 48, r: 20, b: 120 },
    });
    return;
  }

  const umapTrace = {
    x: currentExpressionData.umap.map((p) => p.x),
    y: currentExpressionData.umap.map((p) => p.y),
    text: currentExpressionData.umap.map((p) => `${p.barcode}<br>${p.population}<br>${currentExpressionData.gene}: ${p.value.toFixed(3)}`),
    mode: "markers",
    type: "scattergl",
    marker: {
      size: 2,
      color: currentExpressionData.umap.map((p) => p.value),
      colorscale: "Viridis",
      showscale: true,
      colorbar: { title: currentExpressionData.gene },
    },
    name: currentExpressionData.gene,
  };
  Plotly.newPlot("expression-plot", [umapTrace], {
    title: `${currentExpressionData.gene} expression`,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.9)",
    margin: { t: 36, l: 48, r: 20, b: 40 },
  });
}

async function refreshResults() {
  await loadUmap();
  await loadExpression();
}

function downloadUmapImage() {
  const jobId = document.getElementById("results-job-id").value.trim();
  if (!jobId || !currentUmapData) {
    return;
  }
  const mode = document.getElementById("umap-mode").value;
  window.open(`/api/jobs/${jobId}/umap/pdf?mode=${encodeURIComponent(mode)}`, "_blank");
}

function downloadExpressionImage() {
  const jobId = document.getElementById("results-job-id").value.trim();
  const gene = document.getElementById("gene-query").value.trim();
  if (!jobId || !gene || !currentExpressionData) {
    return;
  }
  const mode = document.getElementById("expression-mode").value;
  window.open(
    `/api/jobs/${jobId}/expression/pdf?gene=${encodeURIComponent(gene)}&mode=${encodeURIComponent(mode)}`,
    "_blank",
  );
}
