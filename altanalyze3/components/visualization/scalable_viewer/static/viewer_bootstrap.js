"use strict";
/* scalable_viewer bootstrap.
 *
 * scALABLE's own app.js is loaded before this file and provides every plot. This
 * script only:
 *   1. removes the upload workflow a precomputed bundle cannot use,
 *   2. auto-loads the first catalog entry (there is no Load button),
 *   3. adds a dataset selector and a differential-contrast selector,
 *   4. adds a DotPlot mode, which scALABLE does not have, drawn with the same
 *      Plotly layout conventions and exported through scALABLE's own vector-PDF
 *      helper (app.js:460 exportPlotlyElementToVectorPdf).
 * It draws no UMAP, violin, heatmap, volcano, network or GO plot of its own.
 */

(function () {
  const SV = { datasets: [], current: null, contrasts: [], contrast: null };
  window.__SCALABLE_VIEWER__ = SV;

  // ---------------------------------------------------------------- utilities

  function el(id) { return document.getElementById(id); }

  function hide(node) { if (node) node.classList.add("hidden"); }

  async function getJson(url) {
    const resp = await fetch(url, { cache: "no-store" });
    const data = await resp.json();
    if (!resp.ok) throw new Error(data.detail || `${url} failed (${resp.status})`);
    return data;
  }

  // ------------------------------------------------------- chrome adjustments

  function stripRunWorkflow() {
    // Hidden, never removed: scALABLE's own code queries these nodes on every status
    // update (app.js updateWorkflowPanels), so the DOM has to stay intact.
    const runTab = document.querySelector('.workspace-tab-btn[data-tab="run"]');
    if (runTab) runTab.style.display = "none";
    const runPanel = document.querySelector('.workspace-panel[data-tab-panel="run"]');
    if (runPanel) runPanel.style.display = "none";
    hide(el("reset-data-btn"));
    // scALABLE blocks Explore when the Run tab's reference select no longer matches the
    // loaded job (app.js updateReferenceChangeState -> selectedReferenceDiffersFromLoadedJob).
    // A bundle has no Run tab and cannot be realigned, so that check never applies here.
    if (typeof selectedReferenceDiffersFromLoadedJob === "function") {
      window.selectedReferenceDiffersFromLoadedJob = () => false;
    }
    // A bundle's differentials are precomputed; nothing can be queued or reset.
    if (typeof updateResetDataButton === "function") {
      window.updateResetDataButton = () => hide(el("reset-data-btn"));
    }
    const originalDifferentialUi = window.updateDifferentialUi;
    if (typeof originalDifferentialUi === "function") {
      window.updateDifferentialUi = function (state) {
        originalDifferentialUi(state);
        stripDifferentialRunConfig();
        hide(el("differential-run-btn"));
        ["differential-sample-field", "differential-group1", "differential-group2",
         "differential-comparison-type", "differential-population-col"].forEach((id) => {
          const node = el(id);
          if (node) node.disabled = true;
        });
      };
    }
    const lede = document.querySelector(".hero-copy .lede");
    if (lede) {
      lede.textContent =
        "Precomputed single-cell atlases: explore aligned cell states, marker programs, " +
        "cell communication and group comparisons.";
    }
  }

  /* CHANGE 2: take the run-configuration controls out of the "4. Differential"
   * panel. A bundle's contrasts are precomputed, so "Group values from",
   * "Comparison Type", the two group list boxes, the run progress bar and the
   * "run cellHarmony-differential" notice configure nothing.
   *
   * The nodes are marked, not deleted. scALABLE's updateDifferentialUi writes into
   * #differential-group1, #differential-group2, #differential-sample-field,
   * #differential-progress and #differential-message on every status poll, so the
   * DOM has to stay intact. `.sv-removed` is `display: none !important`, which takes
   * them out of the layout and beats the `hidden` class scALABLE toggles.
   *
   * Kept: "Cell-state aligned to" (#differential-population) and the precomputed
   * contrast selector (#sv-contrast). */
  function stripDifferentialRunConfig() {
    const panel = el("differential-panel");
    if (!panel) return;
    const kill = (node) => { if (node) node.classList.add("sv-removed"); };
    const wrapperOf = (id) => { const node = el(id); return node ? node.closest(".two-col") : null; };
    kill(wrapperOf("differential-sample-field"));   // "Group values from"
    kill(el("differential-comparison-type-field")); // "Comparison Type"
    kill(wrapperOf("differential-group1"));         // Group 1 + Group 2 share one .two-col
    kill(panel.querySelector(".progress-shell"));   // progress bar + its % label
    kill(el("differential-message"));               // "Differential settings changed..." notice
    kill(el("differential-run-btn"));
    // A hidden control that is still `required` blocks HTML form validation.
    ["differential-sample-field", "differential-group1", "differential-group2"]
      .forEach((id) => { const node = el(id); if (node) node.required = false; });
  }

  /* Panel headings replaced by an inline Dot size control.
   *
   * scALABLE has ONE global #plot-dot-scale, which lived in the deleted "3. Results"
   * panel and drove both plots (app.js:4525 getPlotDotScale, app.js:5099). Each
   * visualization panel now carries its own select, placed on the heading's own row
   * so the header does not get taller. getPlotDotScale is re-pointed at the select
   * of the panel currently rendering; it keeps scALABLE's own `value * 4` scaling. */
  function installPanelDotSize() {
    ["viz1", "viz2"].forEach((panelKey) => {
      const topline = document.querySelector(
        `#baseline-results-view .panel:nth-of-type(${panelKey === "viz1" ? 1 : 2}) .panel-head-topline`);
      if (!topline || el(`sv-dotsize-${panelKey}`)) return;
      const heading = topline.querySelector("h2");
      if (heading) heading.classList.add("sv-removed");
      const wrap = document.createElement("label");
      wrap.className = "sv-dotsize";
      wrap.innerHTML =
        `<span>Dot size</span><select id="sv-dotsize-${panelKey}" class="compact-select">`
        + '<option value="0.25">0.25</option><option value="0.5" selected>0.5</option>'
        + '<option value="0.75">0.75</option><option value="1">1</option></select>';
      topline.insertBefore(wrap, topline.firstChild);
      el(`sv-dotsize-${panelKey}`).addEventListener("change", () => {
        renderVisualizationPanel(panelKey);
      });
    });
    if (SV.dotSizePatchInstalled) return;
    const originalDotScale = window.getPlotDotScale;
    if (typeof originalDotScale !== "function") return;
    window.getPlotDotScale = function () {
      const select = el(`sv-dotsize-${SV.activePanel || "viz1"}`);
      const value = Number(select && select.value);
      if (Number.isFinite(value) && value > 0) return value * 4;
      return originalDotScale();
    };
    SV.dotSizePatchInstalled = true;
  }

  function buildDatasetSelector() {
    const host = document.querySelector(".hero-links");
    if (!host || el("sv-dataset")) return;
    const wrap = document.createElement("label");
    wrap.className = "sv-chip";
    wrap.innerHTML = '<span>Dataset</span><select id="sv-dataset"></select>';
    host.prepend(wrap);
    el("sv-dataset").addEventListener("change", (event) => {
      loadDataset(event.target.value).catch((err) => console.warn(err));
    });
  }

  function buildContrastSelector() {
    const host = el("differential-panel");
    if (!host || el("sv-contrast")) return;
    const wrap = document.createElement("label");
    wrap.className = "field sv-contrast-field";
    wrap.innerHTML =
      '<span>Precomputed comparison</span><select id="sv-contrast" class="compact-select"></select>';
    host.insertBefore(wrap, host.firstChild.nextSibling);
    el("sv-contrast").addEventListener("change", async (event) => {
      const jobId = el("results-job-id").value.trim();
      if (!jobId) return;
      await fetch(`/api/jobs/${jobId}/differential/select?contrast=${encodeURIComponent(event.target.value)}`,
        { method: "POST" });
      SV.contrast = event.target.value;
      await pollStatus(jobId);
    });
  }

  function fillContrastSelector(contrasts, selected) {
    const select = el("sv-contrast");
    if (!select) return;
    select.innerHTML = "";
    contrasts.forEach((entry) => {
      const option = document.createElement("option");
      option.value = entry.id;
      option.textContent = `${entry.comparison} (${entry.n_rows} rows)`;
      if (entry.id === selected) option.selected = true;
      select.appendChild(option);
    });
    select.parentElement.classList.toggle("hidden", contrasts.length === 0);
  }

  // ------------------------------------------------------------- state colours

  /* scALABLE colours cell types with a generated paired ramp
   * (app.js buildReferencePreviewColorMap). A bundle ships explicit colours in its
   * metadata JSON, and the server applies the same map when it renders the PDF.
   * This wrapper uses the bundle map when every label is a known cell state, and
   * otherwise defers to scALABLE's own function. No expression, heatmap or volcano
   * gradient is touched. */
  function installBundleStateColors() {
    if (SV.colorPatchInstalled) return;
    const original = window.buildReferencePreviewColorMap;
    if (typeof original !== "function") return;
    window.buildReferencePreviewColorMap = function (populations) {
      const map = SV.stateColors;
      const labels = (populations || []).map((p) => String(p));
      if (map && labels.length && labels.every((label) => map[label])) {
        return new Map(labels.map((label) => [label, map[label]]));
      }
      return original(populations);
    };
    SV.colorPatchInstalled = true;
  }

  // ------------------------------------------------------------------ dot plot

  function installDotPlotMode() {
    if (typeof BASE_VISUALIZATION_MODES === "undefined") return;
    if (BASE_VISUALIZATION_MODES.some((mode) => mode.value === "dotplot")) return;
    BASE_VISUALIZATION_MODES.push({ value: "dotplot", label: "DotPlot" });

    // The DotPlot takes a free-text gene list. Blank means "top marker of every cell
    // state", so the box is cleared whenever the panel switches into DotPlot and the
    // gene field is shown, which scALABLE hides for non-gene modes.
    const originalModeOptions = window.updateExpressionModeOptions;
    window.updateExpressionModeOptions = function () {
      originalModeOptions();
      ["viz1", "viz2"].forEach((panelKey) => {
        const modeSelect = el(panelElementId(panelKey, "mode"));
        const geneField = el(panelElementId(panelKey, "gene-field"));
        const geneInput = el(panelElementId(panelKey, "gene-query"));
        if (!modeSelect || !geneField || !geneInput) return;
        if (modeSelect.value !== "dotplot") return;
        geneField.classList.remove("hidden");
        geneInput.placeholder = "blank = top marker of every cell state";
        const label = geneField.querySelector("span");
        if (label) label.textContent = "Genes (comma separated)";
        if (geneInput.dataset.svDotplot !== "1") {
          geneInput.value = "";
          geneInput.dataset.svDotplot = "1";
        }
      });
      ["viz1", "viz2"].forEach((panelKey) => {
        const modeSelect = el(panelElementId(panelKey, "mode"));
        const geneInput = el(panelElementId(panelKey, "gene-query"));
        if (modeSelect && geneInput && modeSelect.value !== "dotplot") {
          delete geneInput.dataset.svDotplot;
        }
      });
    };

    const originalLoad = window.loadVisualizationPanel;
    window.loadVisualizationPanel = async function (panelKey) {
      const mode = getPanelSelectValue(panelKey, "mode");
      if (mode !== "dotplot") return originalLoad(panelKey);
      const jobId = el("results-job-id").value.trim();
      if (!jobId) return;
      try {
        const genes = String(el(panelElementId(panelKey, "gene-query"))?.value || "").trim();
        const query = genes ? `?genes=${encodeURIComponent(genes)}` : "";
        const data = await getJson(`/api/jobs/${jobId}/dotplot${query}`);
        panelPlotData[panelKey] = { source: "dotplot", payload: data };
        renderVisualizationPanel(panelKey);
      } catch (err) {
        panelPlotData[panelKey] = { source: "error", payload: { message: err.message } };
        renderVisualizationPanel(panelKey);
      }
    };

    const originalRender = window.renderVisualizationPanel;
    window.renderVisualizationPanel = function (panelKey) {
      // getPlotDotScale takes no panel argument (app.js:4525), so record which panel
      // is rendering before delegating; the per-panel Dot size select reads this.
      SV.activePanel = panelKey;
      const mode = getPanelSelectValue(panelKey, "mode");
      if (mode !== "dotplot") return originalRender(panelKey);
      const data = panelPlotData[panelKey];
      if (!data || data.source === "error") return originalRender(panelKey);
      renderDotPlot(panelKey, data.payload || {});
    };
  }

  function renderDotPlot(panelKey, payload) {
    const genes = payload.genes || [];
    const states = payload.states || [];
    const mean = payload.mean || [];
    const frac = payload.frac || [];
    const x = [], y = [], size = [], color = [], text = [];
    let maxMean = 0;
    mean.forEach((row) => row.forEach((v) => { if (v > maxMean) maxMean = v; }));
    genes.forEach((gene, gi) => {
      states.forEach((state, si) => {
        const m = (mean[gi] || [])[si] || 0;
        const f = (frac[gi] || [])[si] || 0;
        x.push(state);
        y.push(gene);
        size.push(4 + 18 * f);
        color.push(m);
        text.push(`${gene}<br>${state}<br>mean=${m.toFixed(3)}<br>detected=${(100 * f).toFixed(1)}%`
          + `<br>n=${(payload.state_n || [])[si] || 0} cells`);
      });
    });
    const height = Math.max(420, 22 * genes.length + 220);
    Plotly.newPlot(panelPlotId(panelKey), [{
      type: "scattergl",
      mode: "markers",
      x, y, text,
      hovertemplate: "%{text}<extra></extra>",
      marker: {
        size, color,
        // scALABLE's expression ramp (app.py:3869 "expression_grey_red").
        colorscale: [[0, "#e5e7eb"], [0.15, "#f3f4f6"], [0.35, "#fecaca"],
                     [0.6, "#f87171"], [1, "#b91c1c"]],
        cmin: 0, cmax: maxMean || 1,
        line: { width: 0.4, color: "#475569" },
        colorbar: { title: { text: "mean", side: "right" }, thickness: 10 },
      },
    }], {
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(255,255,255,0.9)",
      height,
      margin: { t: 24, l: 150, r: 30, b: 170 },
      xaxis: { tickangle: -55, automargin: true, showgrid: true, gridcolor: "#eef2f7" },
      yaxis: { automargin: true, showgrid: true, gridcolor: "#eef2f7", autorange: "reversed" },
      hovermode: "closest",
    });
    const label = payload.is_default
      ? `DotPlot: top marker of every cell state (${genes.length} genes x ${states.length} states). `
        + "Dot size = fraction of cells detected; colour = mean expression."
      : `DotPlot: ${genes.length} gene(s) x ${states.length} states.`;
    setPanelSummary(panelKey, label);
  }

  // -------------------------------------------------------------- auto loading

  async function loadDataset(datasetId) {
    const entry = SV.datasets.find((d) => d.id === datasetId) || SV.datasets[0];
    if (!entry) return;
    SV.current = entry;
    resetExploreResultsReadiness();
    el("results-job-id").value = entry.id;
    const upload = el("upload-job-id"); if (upload) upload.value = entry.id;
    const qc = el("qc-job-id"); if (qc) qc.value = entry.id;
    SV.contrasts = entry.contrasts || [];
    fillContrastSelector(SV.contrasts, SV.contrast);
    try {
      SV.stateColors = (await getJson(`/api/jobs/${entry.id}/state-colors`)).colors || null;
      installBundleStateColors();
    } catch (err) { console.warn("state colours unavailable", err); }
    await pollStatus(entry.id);
    // pollStatus starts scALABLE's own async warm-up (ensureExploreResultsReady).
    // Open Explore as soon as that reports ready; there is no Load button.
    for (let attempt = 0; attempt < 600; attempt += 1) {
      if (areExploreResultsReady(entry.id)) break;
      await new Promise((resolve) => setTimeout(resolve, 250));
    }
    setExplorerTab("explore");
  }

  async function boot() {
    stripRunWorkflow();
    stripDifferentialRunConfig();
    installPanelDotSize();
    buildDatasetSelector();
    buildContrastSelector();
    installDotPlotMode();
    const catalog = await getJson("/api/catalog");
    SV.datasets = catalog.datasets || [];
    const select = el("sv-dataset");
    SV.datasets.forEach((entry) => {
      const option = document.createElement("option");
      option.value = entry.id;
      option.textContent = `${entry.label} (${Number(entry.n_cells).toLocaleString()} cells)`;
      select.appendChild(option);
    });
    if (!SV.datasets.length) return;
    select.parentElement.classList.toggle("hidden", SV.datasets.length < 2);
    await loadDataset(SV.datasets[0].id);
  }

  document.addEventListener("DOMContentLoaded", () => {
    // scALABLE's own DOMContentLoaded handler runs first; it initializes the shell.
    setTimeout(() => { boot().catch((err) => console.error("[scalable_viewer]", err)); }, 0);
  });
}());
