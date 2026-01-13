"use strict";

const MAX_SAMPLES = 7;
let registry = window.__REFERENCE_REGISTRY__ || { species: [] };
let sampleCount = 0;
let pollTimer = null;

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
  row.dataset.index = sampleCount;
  row.innerHTML = `
    <div class="mb-2">
      <label>Sample name</label>
      <input type="text" class="form-control sample-name" placeholder="e.g. Sample_${sampleCount + 1}" required>
    </div>
    <div class="mb-2">
      <label>H5/H5AD file</label>
      <input type="file" class="form-control sample-file" accept=".h5,.h5ad" required>
    </div>
    <button type="button" class="btn btn-sm btn-outline-danger remove-sample">Remove</button>
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
  document.getElementById("load-umap-btn").addEventListener("click", () => loadUmap());
  document.getElementById("load-expression-btn").addEventListener("click", () => loadExpression());
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

  rows.forEach((row) => {
    const nameInput = row.querySelector(".sample-name");
    const fileInput = row.querySelector(".sample-file");
    if (!nameInput.value || !fileInput.files.length) {
      return;
    }
    formData.append("sample_names[]", nameInput.value.trim());
    formData.append("files[]", fileInput.files[0]);
  });

  try {
    const resp = await fetch("/api/jobs", { method: "POST", body: formData });
    const data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.error || "Failed to create job.");
    }
    document.getElementById("job-created").classList.remove("d-none");
    document.getElementById("job-created").textContent = `Job created: ${data.job_id}`;
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
  };
  try {
    let resp = await fetch(`/api/jobs/${jobId}/qc`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    let data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.error || "Failed to save QC.");
    }
    resp = await fetch(`/api/jobs/${jobId}/run`, { method: "POST" });
    data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.error || "Failed to queue job.");
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
      throw new Error(data.error || "Status request failed.");
    }
    document.getElementById("job-progress").style.width = `${data.progress || 0}%`;
    document.getElementById("job-progress").textContent = `${data.progress || 0}%`;
    document.getElementById("job-log").textContent = (data.log_tail || []).join("");
    if (data.status === "completed" && pollTimer) {
      clearInterval(pollTimer);
      populateDownloadLinks(jobId);
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
    const btn = document.createElement("a");
    btn.className = "btn btn-outline-success";
    btn.href = `/api/jobs/${jobId}/download/${key}`;
    btn.textContent = `Download ${key}`;
    container.appendChild(btn);
  });
}

async function loadUmap() {
  const jobId = document.getElementById("results-job-id").value.trim();
  if (!jobId) {
    alert("Enter a job id.");
    return;
  }
  try {
    const resp = await fetch(`/api/jobs/${jobId}/umap`);
    const data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.error || "UMAP not ready.");
    }
    const refTrace = {
      x: data.reference.map((p) => p.x),
      y: data.reference.map((p) => p.y),
      text: data.reference.map((p) => `${p.barcode}<br>${p.population}`),
      mode: "markers",
      type: "scattergl",
      marker: { color: "#94a3b8", size: 4 },
      name: "Reference",
    };
    const qryTrace = {
      x: data.query.map((p) => p.x),
      y: data.query.map((p) => p.y),
      text: data.query.map((p) => `${p.barcode}<br>${p.population}`),
      mode: "markers",
      type: "scattergl",
      marker: { color: "#e879f9", size: 4 },
      name: "Query",
    };
    Plotly.newPlot("umap-plot", [refTrace, qryTrace], {
      margin: { t: 20 },
      legend: { orientation: "h" },
      hovermode: "closest",
    });
  } catch (err) {
    alert(err.message);
  }
}

async function loadExpression() {
  const jobId = document.getElementById("results-job-id").value.trim();
  const gene = document.getElementById("gene-query").value.trim();
  if (!jobId || !gene) {
    alert("Enter a job id and gene symbol.");
    return;
  }
  try {
    const resp = await fetch(`/api/jobs/${jobId}/expression?gene=${encodeURIComponent(gene)}`);
    const data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.error || "Expression unavailable.");
    }

    const scatter = {
      x: data.scatter.map((d) => d.population),
      y: data.scatter.map((d) => d.value),
      mode: "markers",
      marker: { color: "#3b82f6", size: 8 },
      type: "scatter",
      name: "Cells",
    };

    const violinTraces = data.violin.map((group) => ({
      type: "violin",
      name: group.population,
      y: group.values,
      points: "none",
      box: { visible: true },
      meanline: { visible: true },
      spanmode: "hard",
    }));

    Plotly.newPlot("expression-plot", [...violinTraces, scatter], {
      title: `${data.gene} expression`,
      margin: { t: 40 },
      xaxis: { title: "Population" },
      yaxis: { title: "Expression (a.u.)" },
    });
  } catch (err) {
    alert(err.message);
  }
}
