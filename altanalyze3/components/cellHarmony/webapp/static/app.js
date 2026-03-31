"use strict";

const MAX_SAMPLES = 7;
let registry = window.__REFERENCE_REGISTRY__ || { species: [] };
let sampleCount = 0;
let pollTimer = null;
let currentUmapData = null;
let currentExpressionData = null;
let loadedResultsJobId = null;
let currentDifferentialState = null;
let currentDifferentialGene = "";
let currentDifferentialPopulation = "";
let differentialCy = null;
let lastDownloadArtifactSignature = "";

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
  document.getElementById("differential-viz-mode").addEventListener("change", () => {
    syncDifferentialPopulationSelect(currentDifferentialState);
    updateDifferentialDownloadButton();
    loadDifferentialVisualization();
  });
  document.getElementById("differential-result-population").addEventListener("change", () => {
    currentDifferentialGene = "";
    loadDifferentialVisualization();
  });
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
    document.getElementById("upload-job-id").value = data.job_id;
    const uploadJobIdDisplay = document.getElementById("upload-job-id-display");
    uploadJobIdDisplay.textContent = data.job_id;
    uploadJobIdDisplay.classList.remove("hidden");
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
  document.getElementById("qc-cell-status").textContent = buildQcCellSummary(data);

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

  updateDifferentialUi(data.differential_ui || null);

  if (data.status === "completed") {
    const artifactSignature = JSON.stringify(Object.keys(data.artifacts || {}).sort());
    if (artifactSignature !== lastDownloadArtifactSignature) {
      lastDownloadArtifactSignature = artifactSignature;
      populateDownloadLinks(jobId, data);
    }
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

function buildQcCellSummary(data) {
  const lines = data.log_tail || [];
  let beforeQc = null;
  let afterMinGenes = null;
  let afterMinCounts = null;
  let afterMito = null;

  for (const line of lines) {
    let match = line.match(/(?:reimported\s+)?adata shape:\s*\((\d+),/i);
    if (match) {
      beforeQc = Number(match[1]);
      continue;
    }
    match = line.match(/Cells remaining after min_genes .* filtering:\s*(\d+)/i);
    if (match) {
      afterMinGenes = Number(match[1]);
      continue;
    }
    match = line.match(/Cells remaining after min_counts .* filtering:\s*(\d+)/i);
    if (match) {
      afterMinCounts = Number(match[1]);
      continue;
    }
    match = line.match(/Cells remaining after mito-percent filtering:\s*(\d+)/i);
    if (match) {
      afterMito = Number(match[1]);
    }
  }

  const segments = [];
  if (beforeQc !== null) {
    segments.push(`Cells before QC: ${beforeQc.toLocaleString()}`);
  }
  if (afterMinGenes !== null) {
    segments.push(`After min genes: ${afterMinGenes.toLocaleString()}`);
  }
  if (afterMinCounts !== null) {
    segments.push(`After min counts: ${afterMinCounts.toLocaleString()}`);
  }
  if (afterMito !== null) {
    segments.push(`After mito filter: ${afterMito.toLocaleString()}`);
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
    const resp = await fetch(`/api/jobs/${jobId}/status`);
    data = await resp.json();
    if (!resp.ok) {
      return;
    }
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

  if (state.status === "completed" && (state.result_populations || []).length) {
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
      payload = await fetchDifferentialJson(`/api/jobs/${jobId}/differential/interactive/heatmap?population=${encodeURIComponent(population)}`);
      renderDifferentialHeatmap(payload);
    } else if (mode === "volcano") {
      payload = await fetchDifferentialJson(`/api/jobs/${jobId}/differential/interactive/volcano?population=${encodeURIComponent(population)}`);
      renderDifferentialVolcano(payload);
    } else if (mode === "network") {
      payload = await fetchDifferentialJson(`/api/jobs/${jobId}/differential/interactive/network?population=${encodeURIComponent(population)}`);
      renderDifferentialNetwork(payload);
    } else if (mode === "go") {
      payload = await fetchDifferentialJson(`/api/jobs/${jobId}/differential/interactive/go?population=${encodeURIComponent(population)}`);
      renderDifferentialGo(payload);
    }
    const nextGene = currentDifferentialGene || (payload && payload.default_gene) || "";
    if (nextGene) {
      currentDifferentialGene = nextGene;
      await loadDifferentialGeneDetail(nextGene, population);
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
          [0, "#06b6d4"],
          [0.5, "#111827"],
          [1, "#eab308"],
        ],
        zmid: 0,
        customdata: directions,
        hovertemplate: "%{y}<br>%{x}<br>log2FC=%{z:.3f}<br>%{customdata}<extra></extra>",
        colorbar: { title: "log2FC" },
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
        marker: { color: "#2563eb", size: 8, opacity: 0.72 },
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
        marker: { color: "#dc2626", size: 8, opacity: 0.72 },
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
  plot.classList.remove("hidden");
  document.getElementById("differential-plot-empty").classList.add("hidden");
  Plotly.newPlot(
    plot,
    [
      {
        x: ordered.map((term) => term.score),
        y: ordered.map((term) => term.term_name),
        type: "bar",
        orientation: "h",
        customdata: ordered.map((term) => [term.selected_gene, term.direction, (term.overlap_genes || []).join(", ")]),
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
      yaxis: { automargin: true },
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

function renderDifferentialNetwork(payload) {
  const plot = document.getElementById("differential-plot-area");
  Plotly.purge(plot);
  if (differentialCy) {
    differentialCy.destroy();
    differentialCy = null;
  }
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
          "line-color": (edge) => (edge.data("direction") === "up" ? "#fca5a5" : "#7dd3fc"),
          "target-arrow-color": (edge) => (edge.data("direction") === "up" ? "#fca5a5" : "#7dd3fc"),
          "target-arrow-shape": "triangle",
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

async function loadDifferentialGeneDetail(gene, population) {
  const jobId = document.getElementById("results-job-id").value.trim();
  if (!gene || !population || !jobId) {
    resetDifferentialGeneDetail();
    return;
  }
  try {
    const payload = await fetchDifferentialJson(
      `/api/jobs/${jobId}/differential/interactive/gene?population=${encodeURIComponent(population)}&gene=${encodeURIComponent(gene)}`
    );
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
  document.getElementById("differential-selected-gene").textContent = `${payload.gene} in ${payload.population}`;
  plot.classList.remove("hidden");
  empty.classList.add("hidden");

  const traces = (payload.groups || []).map((group, index) => ({
    type: "violin",
    name: group.label,
    y: group.values,
    points: "all",
    pointpos: 0,
    jitter: 0.28,
    marker: {
      size: 4,
      opacity: 0.35,
      color: index === 0 ? "#dc2626" : "#2563eb",
    },
    line: {
      color: index === 0 ? "#dc2626" : "#2563eb",
    },
    box: { visible: true },
    meanline: { visible: true },
    hovertemplate: `${group.label}<br>expr=%{y:.3f}<extra></extra>`,
  }));
  Plotly.newPlot(
    plot,
    traces,
    {
      title: `Normalized expression: ${payload.gene}`,
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.94)",
      margin: { t: 52, l: 50, r: 18, b: 48 },
      height: 420,
      yaxis: { title: "Normalized expression" },
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
    return;
  }
  const mode = document.getElementById("differential-viz-mode").value;
  const population = document.getElementById("differential-result-population").value;
  if (mode === "heatmap" && state.heatmap_pdf_url) {
    button.href = state.heatmap_pdf_url;
    button.classList.remove("hidden");
    return;
  }
  if (mode === "network" && population) {
    const selected = (state.networks || []).find((entry) => entry.population === population);
    if (selected?.pdf_url) {
      button.href = selected.pdf_url;
      button.classList.remove("hidden");
      return;
    }
  }
  button.classList.add("hidden");
}

function destroyDifferentialNetwork() {
  if (differentialCy) {
    differentialCy.destroy();
    differentialCy = null;
  }
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
  const plot = document.getElementById("differential-plot-area");
  Plotly.purge(plot);
  plot.classList.add("hidden");
  const empty = document.getElementById("differential-plot-empty");
  empty.textContent = "Run cellHarmony-differential to explore the integrated heatmap, volcano plot, network, or GO terms.";
  empty.classList.remove("hidden");
  document.getElementById("download-differential-left-btn").classList.add("hidden");
  const populationSelect = document.getElementById("differential-result-population");
  populationSelect.innerHTML = "";
  currentDifferentialPopulation = "";
  currentDifferentialGene = "";
  resetDifferentialGeneDetail();
}

function resetDifferentialGeneDetail() {
  const plot = document.getElementById("differential-gene-plot");
  Plotly.purge(plot);
  plot.classList.add("hidden");
  document.getElementById("differential-selected-gene").textContent = "Select a gene from the left panel.";
  const empty = document.getElementById("differential-gene-empty");
  empty.textContent = "Select a gene from the heatmap, volcano plot, network, or GO terms to compare expression between groups.";
  empty.classList.remove("hidden");
  document.getElementById("differential-gene-stats").classList.add("hidden");
}

function setResultMode(mode) {
  const baseline = document.getElementById("baseline-results-view");
  const differential = document.getElementById("differential-results-view");
  if (mode === "differential") {
    baseline.classList.add("hidden");
    differential.classList.remove("hidden");
    setResultsControlsDisabled(true);
    return;
  }
  baseline.classList.remove("hidden");
  differential.classList.add("hidden");
  setResultsControlsDisabled(false);
}

function setResultsControlsDisabled(disabled) {
  const controls = [
    document.getElementById("gene-query"),
    document.getElementById("umap-mode"),
    document.getElementById("expression-mode"),
  ];
  controls.forEach((control) => {
    if (!control) {
      return;
    }
    control.disabled = disabled;
  });
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
      marker: { color: "#e5e7eb", size: 1, opacity: 0.3 },
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
