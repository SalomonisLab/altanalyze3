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
 *      helper (app.js:460 exportPlotlyElementToVectorPdf),
 *   5. redraws the differential GO scatter with three colour tiers and a legend, on
 *      scALABLE's own axes and click behaviour (installGoSignificanceTiers).
 * It draws no UMAP, violin, heatmap, volcano or network plot of its own.
 */

(function () {
  const SV = { datasets: [], current: null, contrasts: [], contrast: null };
  window.__SCALABLE_VIEWER__ = SV;

  // ---------------------------------------------------------------- utilities

  // scALABLE is the analysis tool; this deployment serves precomputed bundles, so the
  // page is named for what it is. The server already rewrites <title> and <h1> on the
  // way out (scalable_app.py _install_index_override); this repeats it in the DOM so a
  // page served from a browser cache still carries the new name.
  const VIEWER_NAME = "scALABLE-viewer";

  function applyViewerName() {
    if (document.title.trim() === "scALABLE") document.title = VIEWER_NAME;
    document.querySelectorAll("h1").forEach((node) => {
      if (node.textContent.trim() === "scALABLE") node.textContent = VIEWER_NAME;
    });
  }

  /* The header's "How to Use" and "ReadMe" links. index.html:30-31 sends both to
   * scALABLE's own documents, which describe the upload-and-run workflow this deployment
   * removes; the viewer ships its own pair. The server already rewrites the hrefs on the
   * way out (scalable_app.py DOC_LINKS + _install_index_override); this repeats it in the
   * DOM so a page served from a browser cache also lands on the right document. The
   * shared template is not edited. */
  const GITHUB_BLOB =
    "https://github.com/SalomonisLab/altanalyze3/blob/master/altanalyze3/components";
  const DOC_LINKS = {
    [`${GITHUB_BLOB}/cellHarmony/webapp/HOW_TO_USE.md`]:
      `${GITHUB_BLOB}/visualization/scalable_viewer/HOW_TO_USE.md`,
    [`${GITHUB_BLOB}/cellHarmony/webapp/README.md`]:
      `${GITHUB_BLOB}/visualization/scalable_viewer/README.md`,
  };

  function applyDocLinks() {
    document.querySelectorAll("a.docs-link[href]").forEach((node) => {
      const next = DOC_LINKS[node.getAttribute("href")];
      if (next) node.setAttribute("href", next);
    });
  }

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

  /* GO scatter: three colour tiers and a legend.
   *
   * scALABLE drew two traces, blue for `is_selected_positive_sig` and grey for the
   * rest, with `showlegend: false` (app.js:3606-3607, 3687). Blue is GO-Elite's
   * `selected` flag, which is DAG pruning (goelite/prio.py:34-54), not the FDR call, so
   * grey held both "significant but pruned" and "not significant" with nothing on the
   * page to say so. The server now stamps `color_tier` on every term
   * (scalable_app.py _install_goelite_significance_tiers); this draws it.
   *
   * Axes, ranges, click behaviour and hover fields are scALABLE's; the extra hover
   * field is the overlapping-gene count, because a one-gene term can reach Z=38. */
  function installGoSignificanceTiers() {
    const original = window.renderDifferentialGo;
    if (typeof original !== "function") return;
    const FALLBACK_TIERS = [
      { key: "representative", color: "#1f19c7", label: "GO-Elite representative" },
      { key: "significant", color: "#60a5fa", label: "Significant, not representative" },
      { key: "other", color: "#d1d5db", label: "Not significant" },
    ];

    // `geneFilter` is the second argument scALABLE's own renderer now takes
    // (webapp/static/app.js renderDifferentialGo). It must be forwarded, or the
    // Differential gene filter would silently do nothing on the GO view of this viewer.
    window.renderDifferentialGo = function (payload, geneFilter = null) {
      const significance = payload && payload.significance;
      if (!significance || !Array.isArray(significance.tiers)) {
        original(payload, geneFilter);  // no tier stamp: scALABLE's own plot, unchanged
        return;
      }
      destroyDifferentialNetwork();
      const plot = el("differential-plot-area");
      const allTerms = payload.terms || [];
      // Same rule as scALABLE's renderer: a term survives when its overlapping genes
      // meet the filter set (the typed gene plus the genes it interacts with).
      const terms = geneFilter
        ? allTerms.filter((term) => (term.overlap_genes || [])
          .some((gene) => geneFilter.genes.has(String(gene))))
        : allTerms;
      setDifferentialGeneFilterNote(describeDifferentialGeneFilter(
        geneFilter, payload.population, terms.length, allTerms.length, "GO terms"));
      markDifferentialGeneFilterMatched(terms.length > 0);
      if (!terms.length) {
        renderDifferentialEmpty(geneFilter
          ? `No GO term for ${geneFilter.gene} or its interacting genes in ${payload.population}.`
          : `No GO term enrichment was available for ${payload.population}.`);
        return;
      }
      const ordered = terms.map((term) => {
        const overlapGenes = term.overlap_genes || [];
        const filteredPick = geneFilter && overlapGenes.includes(geneFilter.gene)
          ? geneFilter.gene
          : null;
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
      }).filter((term) => Number.isFinite(term.z_score)
        && Number.isFinite(term.fdr_plot) && term.fdr_plot > 0);
      if (!ordered.length) {
        renderDifferentialEmpty(`No GO term enrichment was available for ${payload.population}.`);
        return;
      }
      const tiers = significance.tiers.length ? significance.tiers : FALLBACK_TIERS;
      // The legend reads "<tier count> of <total>". With the gene filter on, the total
      // must be the filtered term count, or the two numbers carry different denominators.
      const total = geneFilter ? terms.length : (Number(significance.n_terms) || terms.length);
      const hover = (
        "%{customdata[4]}<br>Z-score=%{x:.3f}<br>FDR=%{y:.3e}<br>Fisher p=%{customdata[6]:.3e}"
        + "<br>Overlapping genes: %{customdata[8]}"
        + "<br>Selected gene: %{customdata[0]}"
        + "<br>GO genes: %{customdata[2]}<extra></extra>"
      );
      // Weakest evidence first so the representative terms sit on top.
      const traces = [...tiers].reverse().map((tier) => {
        const rows = ordered.filter((term) => (term.color_tier || "other") === tier.key);
        return {
          x: rows.map((term) => term.z_score),
          y: rows.map((term) => term.fdr_plot),
          type: "scattergl",
          mode: "markers",
          name: `${tier.label} (${rows.length} of ${total})`,
          legendrank: 1000 + tiers.findIndex((entry) => entry.key === tier.key),
          customdata: rows.map((term) => [
            term.selected_gene,
            term.direction,
            (term.overlap_genes || []).join(", "),
            term.overlap_genes || [],
            term.term_name,
            term.z_score,
            term.p_value,
            term.fdr,
            (term.overlap_genes || []).length,
          ]),
          marker: { color: tier.color, size: 11, opacity: tier.key === "other" ? 0.95 : 0.98 },
          hovertemplate: hover,
        };
      });
      plot.classList.remove("hidden");
      hide(el("differential-plot-empty"));
      Plotly.newPlot(
        plot,
        traces,
        {
          title: `GO terms: ${payload.population}`,
          paper_bgcolor: "rgba(0,0,0,0)",
          plot_bgcolor: "rgba(255,255,255,0.94)",
          margin: { t: 56, l: 86, r: 72, b: 200 },
          height: 720,
          xaxis: {
            title: "Z-Score",
            zeroline: false,
            range: [
              Math.min(-10, Math.floor(Math.min(...ordered.map((term) => term.z_score)) - 0.5)),
              Math.max(20, Math.ceil(Math.max(...ordered.map((term) => term.z_score)) + 2.5)),
            ],
          },
          yaxis: { title: "Fishers FDR p", type: "log", autorange: true },
          showlegend: true,
          // Legend under the x-axis label, note under the legend: the same stacking the
          // PDF uses, so the on-screen plot and the download read identically.
          legend: {
            orientation: "h",
            x: 0, xanchor: "left",
            y: -0.11, yanchor: "top",
            traceorder: "normal",
            font: { size: 11 },
            bgcolor: "rgba(0,0,0,0)",
          },
          annotations: [{
            text: String(significance.note || "").replace(/\n/g, "<br>"),
            xref: "paper", yref: "paper",
            x: 0, xanchor: "left",
            y: -0.25, yanchor: "top",
            showarrow: false, align: "left",
            font: { size: 10, color: "#4b5563" },
          }],
          hovermode: "closest",
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
          if (response === null) return;
          const trimmed = response.trim();
          if (trimmed && overlapGenes.includes(trimmed)) gene = trimmed;
        }
        if (gene) {
          currentDifferentialGene = gene;
          loadDifferentialGeneDetail(gene, payload.population);
        }
      });
    };
  }

  /* Panel headings replaced by an inline Dot size control.
   *
   * scALABLE has ONE global #plot-dot-scale, which lived in the deleted "3. Results"
   * panel and drove both plots (app.js:4525 getPlotDotScale, app.js:5099). Each
   * visualization panel now carries its own select, placed on the heading's own row
   * so the header does not get taller. getPlotDotScale is re-pointed at the select
   * of the panel currently rendering; it keeps scALABLE's own `value * 4` scaling. */
  const LARGE_ATLAS_CELLS = 100000;

  function defaultDotSize() {
    const entry = SV.current
      || (SV.datasets && SV.datasets.length ? SV.datasets[0] : null);
    const n = Number(entry && entry.n_cells);
    return Number.isFinite(n) && n > LARGE_ATLAS_CELLS ? "0.25" : "0.5";
  }

  function installPanelDotSize() {
    const defaultDot = defaultDotSize();
    ["viz1", "viz2"].forEach((panelKey) => {
      const topline = document.querySelector(
        `#baseline-results-view .panel:nth-of-type(${panelKey === "viz1" ? 1 : 2}) .panel-head-topline`);
      if (!topline || el(`sv-dotsize-${panelKey}`)) return;
      const heading = topline.querySelector("h2");
      if (heading) heading.classList.add("sv-removed");
      const wrap = document.createElement("label");
      wrap.className = "sv-dotsize";
      wrap.innerHTML =
        // A large atlas overplots at 0.5: the dots merge and structure disappears.
        // Above 100k cells the default drops to 0.25. The user can still pick any value.
        `<span>Dot size</span><select id="sv-dotsize-${panelKey}" class="compact-select">`
        + `<option value="0.25"${defaultDot === "0.25" ? " selected" : ""}>0.25</option>`
        + `<option value="0.5"${defaultDot === "0.5" ? " selected" : ""}>0.5</option>`
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

  // ------------------------------------------------- violin grouping covariate
  //
  /* scALABLE's Violin groups cells by one column only, the cluster key: the server walks
   * `cache_entry["populations"]` (app.py:3552), which is `adata.obs[cluster_key]`. This
   * adds a "Covariate" select immediately to the right of "Select molecule" that names
   * the grouping column instead.
   *
   * Options are the SAME list "Annotation 1" offers - `fields` from
   * GET /api/jobs/{id}/display-filters (app.js:4430 populates Annotation 1 from it), so
   * the two controls can never disagree about what this dataset carries.
   *
   * The default is the bundle's cluster key, so the plot on first paint is the plot
   * scALABLE drew before. The server proves it rather than assuming it: a request that
   * names the cluster key still takes the regrouping path, and must return the payload a
   * request with no `covariate` parameter returns.
   *
   * One select per panel, like Dot size, so panel 1 and panel 2 stay independent.
   *
   * ONLY the violin request carries the parameter. A UMAP, heatmap or network request is
   * sent exactly as scALABLE sent it before, because the covariate regroups nothing there.
   *
   * High cardinality: the cap is scALABLE's own - order the groups by mean expression
   * descending and draw at most 10 (app.py:3564). It is not raised, because keeping it is
   * what makes the default grouping identical. It is not silent either: the server returns
   * the group total and the plot title and caption both name it, so "top 10 of 18 pool
   * groups by mean" is on screen whenever groups are left out. */

  const VIZ_PANELS = ["viz1", "viz2"];

  function covariateValue(panelKey) {
    const select = el(`sv-covariate-${panelKey}`);
    return select ? String(select.value || "").trim() : "";
  }

  function syncCovariateVisibility() {
    VIZ_PANELS.forEach((panelKey) => {
      const field = el(`sv-covariate-field-${panelKey}`);
      if (!field) return;
      field.classList.toggle("hidden", getPanelSelectValue(panelKey, "mode") !== "violin");
    });
  }

  function buildCovariateControls() {
    VIZ_PANELS.forEach((panelKey) => {
      if (el(`sv-covariate-field-${panelKey}`)) return;
      const geneField = el(panelElementId(panelKey, "gene-field"));
      if (!geneField || !geneField.parentElement) return;
      const wrap = document.createElement("label");
      // Same classes the gene field carries, so it takes the next cell of the same
      // .viz-controls-grid row (styles.css:925) instead of starting a row of its own.
      wrap.className = "field inline-field compact-inline-field hidden";
      wrap.id = `sv-covariate-field-${panelKey}`;
      wrap.innerHTML = '<span>Covariate</span>'
        + `<select id="sv-covariate-${panelKey}" class="compact-select"></select>`;
      geneField.parentElement.insertBefore(wrap, geneField.nextSibling);
      el(`sv-covariate-${panelKey}`).addEventListener("change", () => {
        loadVisualizationPanel(panelKey);
      });
    });
    syncCovariateVisibility();
  }

  function fillCovariateSelectors(fields, clusterKey) {
    const options = fields || [];
    VIZ_PANELS.forEach((panelKey) => {
      const select = el(`sv-covariate-${panelKey}`);
      if (!select) return;
      const previous = select.value;
      select.innerHTML = "";
      options.forEach((entry) => {
        const option = document.createElement("option");
        option.value = entry.value;
        option.textContent = entry.value === clusterKey
          ? `${entry.label} (cell type)` : entry.label;
        select.appendChild(option);
      });
      const known = (value) => options.some((entry) => entry.value === value);
      // Default = cell type. A panel the user already set keeps its own choice.
      if (known(previous)) select.value = previous;
      else if (known(clusterKey)) select.value = clusterKey;
      else if (options.length) select.value = options[0].value;
    });
    syncCovariateVisibility();
  }

  function installViolinCovariate() {
    buildCovariateControls();
    if (SV.covariatePatchInstalled) return;

    const originalParams = window.getDisplayFilterParams;
    if (typeof originalParams !== "function") return;
    window.getDisplayFilterParams = function (panelKey) {
      const params = originalParams(panelKey);
      if (getPanelSelectValue(panelKey, "mode") === "violin") {
        const value = covariateValue(panelKey);
        if (value) params.set("covariate", value);
      }
      return params;
    };

    // The control belongs to the Violin, so it is shown for the Violin. scALABLE hides a
    // mode's own controls the same way (app.js:970 geneField.classList.toggle).
    const originalModeOptions = window.updateExpressionModeOptions;
    if (typeof originalModeOptions === "function") {
      window.updateExpressionModeOptions = function () {
        originalModeOptions();
        syncCovariateVisibility();
      };
    }

    /* Title and caption. scALABLE hard-codes "(top 10 states by mean)" (app.js:4993),
     * which is right for the cluster key and wrong for anything else. The plot is drawn
     * by scALABLE, then relabelled - no violin trace, colour or axis is rebuilt here. The
     * default covariate is left completely alone, wording included. */
    const originalRender = window.renderPanelExpression;
    if (typeof originalRender === "function") {
      window.renderPanelExpression = function (panelKey, data, mode, dotScale) {
        originalRender(panelKey, data, mode, dotScale);
        const payload = data || {};
        if (mode !== "violin" || !payload.violin_covariate) return;
        if (payload.violin_covariate_error) {
          setPanelSummary(panelKey, payload.violin_covariate_error);
          return;
        }
        if (payload.violin_covariate_is_default) return;
        const covariate = payload.violin_covariate;
        const title = `${payload.gene} by ${covariate} (${payload.violin_note})`;
        const relayout = () => Plotly.relayout(panelPlotId(panelKey), { title });
        try { relayout(); } catch (err) { setTimeout(relayout, 0); }
        const unlabeled = Number(payload.violin_n_cells_unlabeled || 0);
        const parts = [`Grouped by ${covariate}: ${payload.violin_note}.`];
        if (Number(payload.violin_n_groups_total) > Number(payload.violin_n_groups_shown)) {
          parts.push(`${Number(payload.violin_n_groups_total)
            - Number(payload.violin_n_groups_shown)} lower-mean groups are not drawn.`);
        }
        if (unlabeled > 0) {
          parts.push(`${unlabeled.toLocaleString()} cells have no ${covariate} value and are excluded.`);
        }
        const existing = document.getElementById(panelElementId(panelKey, "filter-summary"));
        const prefix = existing && existing.textContent ? `${existing.textContent} ` : "";
        setPanelSummary(panelKey, prefix + parts.join(" "));
      };
    }
    SV.covariatePatchInstalled = true;
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
      const stateColors = await getJson(`/api/jobs/${entry.id}/state-colors`);
      SV.stateColors = stateColors.colors || null;
      SV.clusterKey = stateColors.cluster_key || "";
      installBundleStateColors();
    } catch (err) { console.warn("state colours unavailable", err); }
    try {
      // The same list the "Annotation 1" selector is built from (app.js:4430).
      const filters = await getJson(`/api/jobs/${entry.id}/display-filters`);
      SV.covariateFields = filters.fields || [];
      fillCovariateSelectors(SV.covariateFields, SV.clusterKey);
    } catch (err) { console.warn("covariate list unavailable", err); }
    await pollStatus(entry.id);
    // pollStatus starts scALABLE's own async warm-up (ensureExploreResultsReady).
    // Open Explore as soon as that reports ready; there is no Load button.
    for (let attempt = 0; attempt < 600; attempt += 1) {
      if (areExploreResultsReady(entry.id)) break;
      await new Promise((resolve) => setTimeout(resolve, 250));
    }
    setExplorerTab("explore");
  }

  // ------------------------------------------------------------------ Study tab
  // First tab, but NOT the landing tab: loadDataset() ends on setExplorerTab("explore"),
  // so Explore still opens by default. Every value comes from GET /api/study, which
  // reads the LungMAP site database (breath.sqlite) read-only; this file holds no study
  // prose, no accession and no sample name of its own. Layout mirrors a LungMAP study
  // page: title, description, metadata grid, then Tools / Downloads / Samples.

  const SAMPLE_PAGE_SIZE = 25;
  const SAMPLE_SKIP_KEYS = ["record_id", "target_id", "display_order", "applies_to_dataset"];

  function esc(value) {
    return String(value === null || value === undefined ? "" : value)
      .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;");
  }

  function studyRow(k, v) {
    return `<div class="sv-study-row"><span class="sv-study-key">${esc(k)}</span>`
         + `<span class="sv-study-val">${v}</span></div>`;
  }

  // A row that spans every column of the study grid.
  //
  // The grid is `repeat(auto-fit, minmax(300px, 1fr))`, so an ordinary row is
  // about 300px wide. Two fields do not fit that: a paper title runs to six
  // wrapped lines and a 30-author list to twenty-five, turning one cell into a
  // narrow column of text beside a column of white space. Both are given the
  // full width instead, and placed after the compact fields so the grid above
  // them stays a tidy block.
  function studyWideRow(k, v) {
    return `<div class="sv-study-row sv-study-row--wide">`
         + `<span class="sv-study-key">${esc(k)}</span>`
         + `<span class="sv-study-val">${v}</span></div>`;
  }

  function studyLink(url, label) {
    const text = esc(label || url || "");
    if (!url) return text || "—";
    return `<a href="${esc(url)}" target="_blank" rel="noopener">${text}</a>`;
  }

  function prettyKey(key) {
    return String(key).replace(/_/g, " ").replace(/\./g, " ")
      .replace(/\b[a-z]/g, (c) => c.toUpperCase());
  }

  function asNumber(value) {
    const n = Number(String(value).replace(/,/g, ""));
    return Number.isFinite(n) && String(value).trim() !== "" ? n : null;
  }

  function fmtCount(value) {
    const n = asNumber(value);
    return n === null ? (value || "—") : n.toLocaleString();
  }

  /* A scrollable, sortable, paged table with a filter box. Used for the sample list,
   * which is 112 rows from the site database or 178 meta-samples from the bundle. */
  function renderSampleTable(host, columns, rows) {
    // key stays null until a header is clicked, so the first draw keeps the order the
    // endpoint sent: lmdb:has_display_order for site samples, category order for the
    // bundle's meta-samples.
    const state = { key: null, dir: 1, page: 0, filter: "" };

    host.innerHTML =
      '<div class="sv-sample-tools">'
      + '<input type="search" class="sv-sample-filter" placeholder="Filter rows">'
      + '<span class="sv-sample-count"></span>'
      + '<span class="sv-sample-nav">'
      + '<button type="button" class="sv-nav-btn" data-nav="first">&laquo;</button>'
      + '<button type="button" class="sv-nav-btn" data-nav="prev">Prev</button>'
      + '<span class="sv-sample-page"></span>'
      + '<button type="button" class="sv-nav-btn" data-nav="next">Next</button>'
      + '<button type="button" class="sv-nav-btn" data-nav="last">&raquo;</button>'
      + '</span></div>'
      + '<div class="sv-table-scroll"><table class="sv-table">'
      + '<thead><tr>'
      + columns.map((c) => `<th data-col="${esc(c)}"><span>${esc(prettyKey(c))}</span>`
                           + '<i class="sv-sort"></i></th>').join("")
      + '</tr></thead><tbody></tbody></table></div>';

    const body = host.querySelector("tbody");
    const countNode = host.querySelector(".sv-sample-count");
    const pageNode = host.querySelector(".sv-sample-page");
    const filterNode = host.querySelector(".sv-sample-filter");

    function filtered() {
      if (!state.filter) return rows;
      const needle = state.filter.toLowerCase();
      return rows.filter((r) => columns.some(
        (c) => String(r[c] === undefined ? "" : r[c]).toLowerCase().indexOf(needle) >= 0));
    }

    function sorted(list) {
      const key = state.key;
      if (!key) return list;
      const copy = list.slice();
      copy.sort((a, b) => {
        const av = a[key], bv = b[key];
        const an = asNumber(av), bn = asNumber(bv);
        if (an !== null && bn !== null) return (an - bn) * state.dir;
        return String(av === undefined ? "" : av)
          .localeCompare(String(bv === undefined ? "" : bv)) * state.dir;
      });
      return copy;
    }

    function draw() {
      const list = sorted(filtered());
      const pages = Math.max(1, Math.ceil(list.length / SAMPLE_PAGE_SIZE));
      if (state.page >= pages) state.page = pages - 1;
      if (state.page < 0) state.page = 0;
      const start = state.page * SAMPLE_PAGE_SIZE;
      const slice = list.slice(start, start + SAMPLE_PAGE_SIZE);
      body.innerHTML = slice.map((r) => "<tr>" + columns.map(
        (c) => `<td>${esc(r[c] === undefined || r[c] === null ? "" : r[c])}</td>`).join("")
        + "</tr>").join("");
      countNode.textContent = list.length === rows.length
        ? `${rows.length} rows`
        : `${list.length} of ${rows.length} rows`;
      pageNode.textContent = list.length
        ? `${start + 1}–${Math.min(start + SAMPLE_PAGE_SIZE, list.length)} / ${list.length}`
        : "0 / 0";
      host.querySelectorAll("th").forEach((th) => {
        const mark = th.getAttribute("data-col") === state.key
          ? (state.dir > 0 ? "▲" : "▼") : "";
        th.querySelector(".sv-sort").textContent = mark;
      });
    }

    host.querySelectorAll("th").forEach((th) => {
      th.addEventListener("click", () => {
        const key = th.getAttribute("data-col");
        if (state.key === key) state.dir = -state.dir;
        else { state.key = key; state.dir = 1; }
        state.page = 0;
        draw();
      });
    });
    host.querySelectorAll(".sv-nav-btn").forEach((b) => {
      b.addEventListener("click", () => {
        const list = filtered();
        const pages = Math.max(1, Math.ceil(list.length / SAMPLE_PAGE_SIZE));
        const nav = b.getAttribute("data-nav");
        if (nav === "first") state.page = 0;
        else if (nav === "prev") state.page -= 1;
        else if (nav === "next") state.page += 1;
        else state.page = pages - 1;
        draw();
      });
    });
    filterNode.addEventListener("input", () => {
      state.filter = filterNode.value.trim();
      state.page = 0;
      draw();
    });
    draw();
  }

  function sampleColumns(rows) {
    const seen = [];
    rows.forEach((row) => Object.keys(row).forEach((key) => {
      if (SAMPLE_SKIP_KEYS.indexOf(key) >= 0) return;
      if (seen.indexOf(key) < 0) seen.push(key);
    }));
    const rest = seen.filter((k) => k !== "name");
    return (seen.indexOf("name") >= 0 ? ["name"] : []).concat(rest);
  }

  function renderStudyPanel(panel, study) {
    const viewer = study.viewer_stats || {};
    const grid =
      studyRow("Dataset ID", esc(study.dataset_id || ""))
      + studyRow("Assay Type", esc(study.assay || "—"))
      + studyRow("Organism", esc(study.organism || "—"))
      + studyRow("Technology", esc((study.technologies || []).join(", ") || "—"))
      + studyRow("Cell Count", esc(fmtCount(study.cell_count)))
      // Every declared accession, not only the first: a study can carry a GEO series
      // and a dbGaP study at once, and hiding the second would understate the release.
      + studyRow("Raw data", (study.accessions || []).length
          ? (study.accessions || []).map((a) => studyLink(a.url, a.name)).join(", ")
          : (study.raw_data && study.raw_data.url
              ? studyLink(study.raw_data.url, study.raw_data.name) : "—"))
      + studyRow("Sample Type", esc(study.sample_type || "—"))
      + studyRow("Age Range", esc((study.age_ranges || []).join(", ") || "—"))
      // The analysis protocol is a field of its own. It used to appear only as one
      // row among the downloads, where a reader looking for the method never saw it.
      + studyRow("Protocol", study.protocol && study.protocol.url
          ? studyLink(study.protocol.url, study.protocol.name) : "—")
      + studyRow("Served in this viewer",
          viewer.n_cells
            ? esc(`${Number(viewer.n_cells).toLocaleString()} metacells x `
                  + `${Number(viewer.n_genes).toLocaleString()} genes, `
                  + `${viewer.n_states} cell states`)
            : "—")
      // Reference and Researchers span the full grid width. Both are long
      // free text and neither is scannable in a 300px column.
      + studyWideRow("Reference", study.reference && study.reference.url
          ? studyLink(study.reference.url, study.reference.name) : "—")
      // The author list. It shipped empty for a full day.
      + studyWideRow("Researchers", (study.researchers || []).length
          ? esc((study.researchers || []).map((r) => r.name).join(", "))
          : "—");

    // Two different notices. is_fallback would mean a DIFFERENT study is on screen.
    // source === "source_tables" means this study's OWN record, read from the TSVs the
    // site database is built from, because the database has not been rebuilt yet.
    const notice = study.is_fallback
      ? '<p class="sv-study-notice">Record ' + esc(study.requested_id)
        + ' is not in the site database yet. Showing ' + esc(study.study_id)
        + ' instead.</p>'
      : (study.source === "source_tables" && study.notice
          ? `<p class="sv-study-notice">${esc(study.notice)}</p>`
          : "");

    const tools = (study.tools || []);
    const toolsHtml = tools.length
      ? '<ul class="sv-study-list">' + tools.map((t) => "<li>"
          + studyLink(t.url || t.path, t.name)
          + (t.description ? ` <span class="sv-dim">${esc(t.description)}</span>` : "")
          + "</li>").join("") + "</ul>"
      : '<p class="sv-dim">The site record lists no tools.</p>';

    const files = (study.files || []);
    const filesHtml = files.length
      ? '<table class="sv-table sv-file-table"><thead><tr><th>File</th><th>Size</th>'
        + "<th>Description</th></tr></thead><tbody>"
        + files.map((f) => "<tr><td>" + studyLink(f.url, f.name) + "</td>"
            + `<td>${esc(f.size || "")}</td><td>${esc(f.description || "")}</td></tr>`).join("")
        + "</tbody></table>"
      : '<p class="sv-dim">The site record lists no files.</p>';

    panel.innerHTML =
      `<h2 class="sv-study-title">${esc(study.title || study.study_id || "")}</h2>`
      + notice
      + `<div class="sv-study-card"><p class="sv-study-desc">${esc(study.description || "")}</p>`
      + `<div class="sv-study-grid">${grid}</div></div>`
      + '<div class="sv-study-tabs">'
      + '<button type="button" class="sv-study-subtab active" data-sub="tools">Tools</button>'
      + '<button type="button" class="sv-study-subtab" data-sub="downloads">Downloads</button>'
      + '<button type="button" class="sv-study-subtab" data-sub="samples">Samples</button></div>'
      + `<div class="sv-study-sub" data-sub-panel="tools">${toolsHtml}</div>`
      + `<div class="sv-study-sub hidden" data-sub-panel="downloads">${filesHtml}</div>`
      + '<div class="sv-study-sub hidden" data-sub-panel="samples">'
      + '<p class="sv-sample-source"></p><div class="sv-sample-table"></div></div>';
    // No provenance footer. A reader does not need this machine's database path, and
    // printing an absolute local path on a study page is noise at best.

    panel.querySelectorAll(".sv-study-subtab").forEach((b) => {
      b.addEventListener("click", () => {
        panel.querySelectorAll(".sv-study-subtab").forEach((x) => x.classList.remove("active"));
        b.classList.add("active");
        panel.querySelectorAll(".sv-study-sub").forEach((p) => {
          p.classList.toggle("hidden",
            p.getAttribute("data-sub-panel") !== b.getAttribute("data-sub"));
        });
      });
    });

    // Samples: the site record when it is the requested study and carries sample rows;
    // otherwise the bundle's own meta-samples. A fallback record's samples belong to a
    // different study, so they are never shown as this viewer's samples. The caption
    // always says which list is on screen and where it came from.
    const host = panel.querySelector(".sv-sample-table");
    const caption = panel.querySelector(".sv-sample-source");
    const dbSamples = study.samples || [];
    const viewerSamples = (study.viewer_samples || {});
    const viewerRows = viewerSamples.rows || [];
    if (dbSamples.length && !study.is_fallback) {
      caption.textContent = `${dbSamples.length} samples from the LungMAP record `
        + `${study.study_id}. Click a column to sort; type to filter.`;
      renderSampleTable(host, sampleColumns(dbSamples), dbSamples);
    } else if (viewerRows.length) {
      caption.textContent = (dbSamples.length
          ? `The ${dbSamples.length} samples of stand-in record ${study.study_id} belong to a `
            + "different study and are not shown. "
          : "The site record carries no sample rows. ")
        + `Showing the ${viewerRows.length} meta-samples of the loaded bundle `
        + `(obs column "${viewerSamples.column}"). Click a column to sort; type to filter.`;
      renderSampleTable(host, viewerSamples.columns || sampleColumns(viewerRows), viewerRows);
    } else {
      caption.textContent = "No samples in the site record and none in the bundle.";
      host.innerHTML = `<p class="sv-dim">${esc(viewerSamples.error || "")}</p>`;
    }
    return panel;
  }

  async function installStudyTab() {
    const bar = document.querySelector(".workspace-tab-btn")?.parentElement;
    const anyPanel = document.querySelector(".workspace-panel[data-tab-panel]");
    if (!bar || !anyPanel) return;

    const btn = document.createElement("button");
    btn.className = "workspace-tab-btn";
    btn.setAttribute("data-tab", "study");
    btn.type = "button";
    btn.textContent = "Study";
    bar.insertBefore(btn, bar.firstChild);          // FIRST tab

    const panel = document.createElement("section");
    panel.className = "workspace-panel hidden";
    panel.setAttribute("data-tab-panel", "study");
    panel.innerHTML = '<p class="sv-dim">Loading the study record…</p>';
    anyPanel.parentElement.insertBefore(panel, anyPanel);

    // scALABLE hides a panel by removing .active, not by adding .hidden, and it
    // re-asserts that on every tab click. Mirror BOTH so the Study panel actually
    // swaps with Explore instead of sitting behind it.
    function showStudy(on) {
      if (on) {
        document.querySelectorAll(".workspace-tab-btn").forEach((b) => b.classList.remove("active"));
        document.querySelectorAll(".workspace-panel[data-tab-panel]").forEach((p) => {
          p.classList.add("hidden"); p.classList.remove("active");
        });
        btn.classList.add("active");
        panel.classList.remove("hidden"); panel.classList.add("active");
      } else {
        btn.classList.remove("active");
        panel.classList.add("hidden"); panel.classList.remove("active");
      }
    }
    SV.showStudy = showStudy;
    btn.addEventListener("click", () => showStudy(true));
    document.querySelectorAll('.workspace-tab-btn:not([data-tab="study"])').forEach((b) => {
      b.addEventListener("click", () => showStudy(false));
    });
    showStudy(false);   // Explore is the landing tab; Study starts closed

    let study = null;
    try {
      study = await getJson("/api/study");
    } catch (err) {
      panel.innerHTML = `<p class="sv-study-error">The study record could not be read: `
        + `${esc(err.message)}</p>`;
      console.warn("[scalable_viewer] /api/study failed", err);
      return;
    }
    SV.study = study;
    if (!study.ok) {
      panel.innerHTML = `<p class="sv-study-error">${esc(study.error || "no study record")}</p>`
        + `<p class="sv-dim">Database: ${esc(study.db_path || "")}<br>Tried: `
        + `${esc((study.candidates || []).join(", "))}</p>`;
      return;
    }
    renderStudyPanel(panel, study);
  }

  async function boot() {
    applyViewerName();
    applyDocLinks();
    stripRunWorkflow();
    stripDifferentialRunConfig();
    installGoSignificanceTiers();
    await installStudyTab();
    buildDatasetSelector();
    buildContrastSelector();
    installDotPlotMode();
    installViolinCovariate();
    const catalog = await getJson("/api/catalog");
    SV.datasets = catalog.datasets || [];
    // Built only now: the default dot size depends on the atlas size, which the
    // catalog carries. Called before the catalog loads it would always see 0 cells.
    installPanelDotSize();
    const select = el("sv-dataset");
    SV.datasets.forEach((entry) => {
      const option = document.createElement("option");
      option.value = entry.id;
      option.textContent = `${entry.label} (${Number(entry.n_cells).toLocaleString()} cells)`;
      select.appendChild(option);
    });
    if (!SV.datasets.length) return;
    select.parentElement.classList.toggle("hidden", SV.datasets.length < 2);

    // A dataset can be named in the URL, so another page can link straight to
    // one atlas instead of landing on whichever happens to be first. Both the
    // bundle id (`?dataset=COPD-metacells`) and the LungMAP dataset accession
    // (`?dataset_id=LMEX0000009416`) are accepted, because the site knows a
    // dataset by its accession and the catalog knows it by its bundle id.
    const wanted = new URLSearchParams(window.location.search);
    const askedBundle = (wanted.get("dataset") || "").trim();
    const askedRecord = (wanted.get("dataset_id") || "").trim().toUpperCase();
    let chosen = SV.datasets[0].id;
    if (askedBundle) {
      const hit = SV.datasets.filter((d) => d.id === askedBundle)[0];
      if (hit) chosen = hit.id;
      else console.warn("[scalable_viewer] no bundle named", askedBundle);
    } else if (askedRecord) {
      // The catalog entry carries the accession under whichever of these the
      // bundle metadata happened to set.
      const hit = SV.datasets.filter((d) =>
        [d.dataset_id, d.lungmap_id, d.accession, d.study_id]
          .filter(Boolean)
          .map((v) => String(v).toUpperCase())
          .indexOf(askedRecord) >= 0)[0];
      if (hit) chosen = hit.id;
      else console.warn("[scalable_viewer] no bundle for record", askedRecord);
    }
    select.value = chosen;
    await loadDataset(chosen);
  }

  document.addEventListener("DOMContentLoaded", () => {
    // scALABLE's own DOMContentLoaded handler runs first; it initializes the shell.
    setTimeout(() => { boot().catch((err) => console.error("[scalable_viewer]", err)); }, 0);
  });
}());
