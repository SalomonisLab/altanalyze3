"use strict";

const MAX_SAMPLES = 7;
const APP_ROOT_PATH = normalizeRootPath(window.__APP_ROOT_PATH__ || "");
let registry = window.__REFERENCE_REGISTRY__ || { species: [] };
let sampleCount = 0;
let pollTimer = null;
let currentUmapData = null;
let currentExpressionData = null;
let loadedResultsJobId = null;
let loadedGeneSuggestionsJobId = null;
let currentJobStatus = "";
let currentJobSpecies = "";
let currentJobReference = "";
let referenceRerunPending = false;
let currentDifferentialState = null;
let currentDifferentialGene = "";
let currentDifferentialPopulation = "";
let differentialCy = null;
let lastDownloadArtifactSignature = "";
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
  initSpeciesSelect();
  initSampleRows();
  hookForms();
  updateWorkflowPanels(null);
  updateDifferentialUi(null);
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
    updateReferenceChangeState();
    loadReferencePreview();
  });

  if (registry.species.length) {
    speciesSelect.value = registry.species[0].id;
    speciesSelect.dispatchEvent(new Event("change"));
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
    clearGeneSuggestions();
    resetDifferentialResults();
    setResultMode("baseline");
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
  document.getElementById("results-form").addEventListener("submit", handleResultsSubmit);
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
    setResultMode("baseline");
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
      lastDownloadArtifactSignature = "";
      currentUmapData = null;
      currentExpressionData = null;
      loadedResultsJobId = null;
      clearGeneSuggestions();
      document.getElementById("download-links").innerHTML = "";
      resetDifferentialResults();
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

function handleResultsSubmit(evt) {
  evt.preventDefault();
  refreshResults();
}

async function loadJobState(jobId) {
  try {
    const resp = await fetch(apiPath(`/jobs/${jobId}/status`));
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
    const resp = await fetch(apiPath(`/jobs/${jobId}/status`));
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
  currentJobStatus = String(data.status || "").trim().toLowerCase();
  currentJobSpecies = String(data.species || currentJobSpecies || "");
  currentJobReference = String(data.reference || currentJobReference || "");
  updateWorkflowPanels(referenceRerunPending ? "uploaded" : (data.status || null));
  document.getElementById("job-progress").style.width = `${data.progress || 0}%`;
  document.getElementById("job-progress-label").textContent = `${data.progress || 0}%`;
  document.getElementById("job-log").textContent = (data.qc_log_tail || []).join("");
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

  updateDifferentialUi(data.differential_ui || null);
  updateReferenceChangeState();

  if (data.status === "completed" && !referenceRerunPending) {
    loadGeneSuggestions(jobId);
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

function updateWorkflowPanels(status) {
  const normalizedStatus = String(status || "").trim().toLowerCase();
  const hasUploadedJob = Boolean(normalizedStatus);
  const hasCompletedAlignment = normalizedStatus === "completed";
  const qcPanel = document.getElementById("qc-panel");
  const resultsPanel = document.getElementById("results-panel");
  const previewPanel = document.getElementById("reference-preview-panel");
  const baseline = document.getElementById("baseline-results-view");
  const differential = document.getElementById("differential-results-view");

  qcPanel.classList.toggle("hidden", !hasUploadedJob);
  resultsPanel.classList.toggle("hidden", !hasCompletedAlignment);
  previewPanel.classList.toggle("hidden", hasUploadedJob);

  if (!hasCompletedAlignment) {
    baseline.classList.add("hidden");
    differential.classList.add("hidden");
  }
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
    const resp = await fetch(apiPath(`/jobs/${jobId}/status`));
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
    btn.href = apiPath(`/jobs/${jobId}/download/${key}`);
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
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/heatmap?population=${encodeURIComponent(population)}`));
      renderDifferentialHeatmap(payload);
    } else if (mode === "volcano") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/volcano?population=${encodeURIComponent(population)}`));
      renderDifferentialVolcano(payload);
    } else if (mode === "network") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/network?population=${encodeURIComponent(population)}`));
      renderDifferentialNetwork(payload);
    } else if (mode === "go") {
      payload = await fetchDifferentialJson(apiPath(`/jobs/${jobId}/differential/interactive/go?population=${encodeURIComponent(population)}`));
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
  const finiteValues = z
    .flat()
    .filter((value) => value !== null && Number.isFinite(Number(value)))
    .map((value) => Number(value));
  const [vmin, vmax] = finiteValues.length ? finiteExtent(finiteValues, -1, 1) : [-1, 1];
  const contrastFactor = 2.0;
  let colorMin = vmin / contrastFactor;
  let colorMax = vmax / contrastFactor;
  if (!(colorMin < colorMax)) {
    colorMin = -1;
    colorMax = 1;
  }
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
          [0, "#87ceeb"],
          [0.5, "#111827"],
          [1, "#ffff00"],
        ],
        zmin: colorMin,
        zmax: colorMax,
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

async function loadDifferentialGeneDetail(gene, population) {
  const jobId = document.getElementById("results-job-id").value.trim();
  if (!gene || !population || !jobId) {
    resetDifferentialGeneDetail();
    return;
  }
  try {
    const payload = await fetchDifferentialJson(
      apiPath(`/jobs/${jobId}/differential/interactive/gene?population=${encodeURIComponent(population)}&gene=${encodeURIComponent(gene)}`)
    );
    currentDifferentialPopulation = payload.population || population;
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
  document.getElementById("differential-selected-gene").textContent = `${payload.gene} in ${payload.population}`;
  plot.classList.remove("hidden");
  empty.classList.add("hidden");
  downloadButton.href = apiPath(
    `/jobs/${document.getElementById("results-job-id").value.trim()}/differential/interactive/gene/pdf?population=${encodeURIComponent(payload.population)}&gene=${encodeURIComponent(payload.gene)}`
  );
  downloadButton.classList.remove("hidden");

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
  document.getElementById("download-differential-gene-btn").classList.add("hidden");
  document.getElementById("download-differential-gene-btn").removeAttribute("href");
  const empty = document.getElementById("differential-gene-empty");
  empty.textContent = "Select a gene from the heatmap, volcano plot, network, or GO terms to compare expression between groups.";
  empty.classList.remove("hidden");
  document.getElementById("differential-gene-stats").classList.add("hidden");
}

function clearGeneSuggestions() {
  document.getElementById("gene-suggestions").innerHTML = "";
  loadedGeneSuggestionsJobId = null;
}

async function loadGeneSuggestions(jobId) {
  if (!jobId || loadedGeneSuggestionsJobId === jobId) {
    return;
  }
  try {
    const resp = await fetch(apiPath(`/jobs/${jobId}/genes`));
    const data = await parseApiResponse(resp);
    if (!resp.ok) {
      throw new Error(data.detail || "Unable to load gene suggestions.");
    }
    const datalist = document.getElementById("gene-suggestions");
    datalist.innerHTML = "";
    (data.genes || []).forEach((gene) => {
      const option = document.createElement("option");
      option.value = gene;
      datalist.appendChild(option);
    });
    loadedGeneSuggestionsJobId = jobId;
  } catch (err) {
    console.warn(err);
  }
}

function setResultMode(mode) {
  const baseline = document.getElementById("baseline-results-view");
  const differential = document.getElementById("differential-results-view");
  if (currentJobStatus !== "completed") {
    baseline.classList.add("hidden");
    differential.classList.add("hidden");
    setResultsControlsDisabled(false);
    return;
  }
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
    const resp = await fetch(apiPath(`/jobs/${jobId}/umap`));
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
  const layout = {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.9)",
    margin: { t: 20, l: 40, r: 20, b: 40 },
    hovermode: "closest",
    xaxis: { showgrid: false, zeroline: false },
    yaxis: { showgrid: false, zeroline: false },
  };
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
    const populations = [...new Set(currentUmapData.query.map((p) => p.population))];
    const colorMap = buildReferencePreviewColorMap(populations);
    const labelPoints = relaxReferencePreviewLabels(
      buildPopulationCentroids(currentUmapData.query),
      currentUmapData.query,
      document.getElementById("umap-plot")
    );
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
    traces.push({
      x: currentUmapData.query.map((p) => p.x),
      y: currentUmapData.query.map((p) => p.y),
      text: currentUmapData.query.map((p) => `${p.barcode}<br>${p.population}`),
      mode: "markers",
      type: "scattergl",
      marker: {
        size: 2,
        opacity: 0.5,
        color: currentUmapData.query.map((p) => colorMap.get(p.population)),
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
      font: {
        size: 11,
        color: "#0f172a",
      },
      bgcolor: "rgba(255,255,255,0)",
      opacity: 1,
    }));
    Object.assign(layout, buildSquareUmapAxes(currentUmapData.query, 0.06));
  }
  if (mode === "relative") {
    layout.legend = { orientation: "h" };
  }
  Plotly.newPlot("umap-plot", traces, layout);
}

async function loadExpression() {
  const jobId = document.getElementById("results-job-id").value.trim();
  const geneInput = document.getElementById("gene-query");
  const gene = geneInput.value.trim();
  if (!jobId || !gene) {
    return;
  }
  try {
    const resp = await fetch(apiPath(`/jobs/${jobId}/expression?gene=${encodeURIComponent(gene)}`));
    const data = await resp.json();
    if (!resp.ok) {
      throw new Error(data.detail || "Expression unavailable.");
    }
    if (data?.gene && data.gene !== gene) {
      geneInput.value = data.gene;
    }
    currentExpressionData = data;
    renderCurrentExpression();
  } catch (err) {
    currentExpressionData = {
      gene,
      requested_gene: gene,
      resolved_gene: null,
      source: "missing",
      message: err.message || "Expression unavailable.",
      scatter: [],
      violin: [],
      umap: [],
    };
    renderCurrentExpression();
  }
}

function renderExpressionMessage(message, title = "Expression unavailable") {
  Plotly.newPlot("expression-plot", [], {
    title,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.9)",
    margin: { t: 48, l: 32, r: 32, b: 32 },
    xaxis: { visible: false },
    yaxis: { visible: false },
    annotations: [
      {
        text: message,
        x: 0.5,
        y: 0.5,
        xref: "paper",
        yref: "paper",
        showarrow: false,
        font: { size: 15, color: "#64748b" },
        align: "center",
      },
    ],
  }, { displayModeBar: false });
}

function renderCurrentExpression() {
  if (!currentExpressionData) {
    return;
  }
  if (!currentExpressionData.umap?.length && !currentExpressionData.violin?.length) {
    renderExpressionMessage(
      currentExpressionData.message || `Gene '${currentExpressionData.requested_gene || currentExpressionData.gene}' was not found.`,
      `${currentExpressionData.requested_gene || currentExpressionData.gene} expression`,
    );
    return;
  }
  const mode = document.getElementById("expression-mode").value;
  if (mode === "violin") {
    if (currentExpressionData.source === "reference") {
      const trace = {
        type: "bar",
        x: currentExpressionData.violin.map((entry) => entry.population),
        y: currentExpressionData.violin.map((entry) => entry.mean),
        marker: { color: "#475569", opacity: 0.85 },
      };
      Plotly.newPlot("expression-plot", [trace], {
        title: `${currentExpressionData.gene} (reference centroid expression, top 10 states)`,
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(255,255,255,0.9)",
        margin: { t: 36, l: 48, r: 20, b: 120 },
      });
      return;
    }
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
      colorscale: [
        [0.0, "#e5e7eb"],
        [0.15, "#f3f4f6"],
        [0.35, "#fecaca"],
        [0.6, "#f87171"],
        [1.0, "#b91c1c"],
      ],
      showscale: true,
      colorbar: { title: currentExpressionData.gene },
    },
    name: currentExpressionData.gene,
  };
  const expressionLayout = {
    title: currentExpressionData.source === "reference"
      ? `${currentExpressionData.gene} reference expression`
      : `${currentExpressionData.gene} expression`,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(255,255,255,0.9)",
    margin: { t: 36, l: 48, r: 20, b: 40 },
    xaxis: { showgrid: false, zeroline: false },
    yaxis: { showgrid: false, zeroline: false },
    annotations: currentExpressionData.message
      ? [
          {
            text: currentExpressionData.message,
            x: 0.5,
            y: 1.08,
            xref: "paper",
            yref: "paper",
            showarrow: false,
            font: { size: 12, color: "#64748b" },
          },
        ]
      : [],
  };
  Object.assign(expressionLayout, buildSquareUmapAxes(currentExpressionData.umap, 0.06));
  Plotly.newPlot("expression-plot", [umapTrace], expressionLayout);
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
  window.open(apiPath(`/jobs/${jobId}/umap/pdf?mode=${encodeURIComponent(mode)}`), "_blank");
}

function downloadExpressionImage() {
  const jobId = document.getElementById("results-job-id").value.trim();
  const gene = document.getElementById("gene-query").value.trim();
  if (!jobId || !gene || !currentExpressionData) {
    return;
  }
  const mode = document.getElementById("expression-mode").value;
  window.open(
    apiPath(`/jobs/${jobId}/expression/pdf?gene=${encodeURIComponent(gene)}&mode=${encodeURIComponent(mode)}`),
    "_blank",
  );
}
