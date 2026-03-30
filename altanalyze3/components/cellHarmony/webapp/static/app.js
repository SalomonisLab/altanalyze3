"use strict";

const MAX_SAMPLES = 7;
let registry = window.__REFERENCE_REGISTRY__ || { species: [] };
let sampleCount = 0;
let pollTimer = null;
let currentUmapData = null;
let currentExpressionData = null;

document.addEventListener("DOMContentLoaded", () => {
  initSpeciesSelect();
  initSampleRows();
  hookForms();
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
  document.getElementById("gene-query").addEventListener("change", () => refreshResults());
  document.getElementById("expression-mode").addEventListener("change", renderCurrentExpression);
  document.getElementById("umap-mode").addEventListener("change", renderCurrentUmap);
  document.getElementById("download-umap-image-btn").addEventListener("click", downloadUmapImage);
  document.getElementById("download-expression-image-btn").addEventListener("click", downloadExpressionImage);
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
    startStatusPolling(jobId);
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
    const progress = data.progress || 0;
    document.getElementById("job-progress").style.width = `${progress}%`;
    document.getElementById("job-progress-label").textContent = `${progress}%`;
    document.getElementById("job-log").textContent = (data.qc_log_tail || []).join("");
    if (data.status === "completed" && pollTimer) {
      clearInterval(pollTimer);
      populateDownloadLinks(jobId);
      document.getElementById("results-job-id").value = jobId;
      if (!document.getElementById("gene-query").value && data.default_gene) {
        document.getElementById("gene-query").value = data.default_gene;
      }
      refreshResults();
    }
  } catch (err) {
    console.warn(err);
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
    populations.forEach((population, index) => {
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
