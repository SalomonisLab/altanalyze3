/* ISV web frontend.
 *   Tab 1 "Cell type × covariate"  -> /api/isoforms : clustered isoforms, covariate-blocked heat strip
 *                                     + isoform structures (IGV-style, top axis).
 *   Tab 2 "Molecule view" (read-level) -> /api/reads : the ENGINE's own per-molecule rows
 *                                     (plot_isoform_structures_by_conditions output), drawn as a
 *                                     per-covariate read pileup with a shared gene model + axis at the
 *                                     bottom -- the interactive analogue of the ISV PDF.
 * Shared: genomic ruler + reference gene-model track, wheel-zoom / drag-pan, mouseover (exon region +
 * coords + protein length), left-click protein / mRNA / ORF FASTA export.
 */
(function () {
  "use strict";
  const $ = (s) => document.querySelector(s);
  const $$ = (s) => Array.from(document.querySelectorAll(s));
  const api = (p, o) => fetch(p, o).then((r) => { if (!r.ok) throw new Error(r.status + " " + p); return r; });
  const state = {
    catalog: null, tab: "heatmap", gene: null, combineBy: "cell_type_x_covariate", panelBy: "covariate", netPanelBy: "covariate",
    samples: new Set(), groups: new Set(), cellTypes: new Set(), cellTypeOrder: [], junctions: new Set(),
    last: null,
  };
  let view = null;
  // any control change auto-renders (debounced so toggling several chips/sliders = one render). No Render button.
  let _renderTimer = null;
  const autoRender = () => { clearTimeout(_renderTimer); _renderTimer = setTimeout(() => { if (state.tab === "network" || state.gene || geneInput.value.trim()) render(); }, 300); };

  const KNOWN = "#1d5fa8", NOVEL = "#b0306b", REF_FILL = "#3d4a5c";
  const cyanYellow = d3.interpolateRgb("#00FFFF", "#FFFF00");
  // professional qualitative palette (RGB hex, non-rainbow) keyed by cluster index, consistent across panels
  const CLUSTER_COLORS = ["#2b6cb0", "#c0531f", "#2f8a4e", "#7a4fb5", "#c2306b", "#0e8a9c", "#b8902a",
    "#4a5bbf", "#cf5aa0", "#2f9b8a", "#9a5b2e", "#56657a", "#1f7a5a", "#a23c5e", "#3d7ec2", "#8a7a1f",
    "#6d4ca0", "#1d8f74", "#bf6a2a", "#8a2f52", "#3a6ea5", "#7b9b3a", "#a64d79", "#2d7d8a"];
  const clusterColor = (i) => CLUSTER_COLORS[(((i | 0) % CLUSTER_COLORS.length) + CLUSTER_COLORS.length) % CLUSTER_COLORS.length];
  // color by the backend's stable structure-keyed color_index (consistent across panels + modes);
  // fall back to per-panel cluster_index if absent.
  const colorFor = (m) => clusterColor(m.color_index != null ? m.color_index : m.cluster_index);

  // ===================================================================== catalog / menus
  async function loadCatalog() {
    const c = await api("/api/catalog").then((r) => r.json());
    state.catalog = c;
    renderMulti("#samples", c.samples.map((s) => ({ id: s.library, label: s.library, sub: s.group })), state.samples);
    renderMulti("#groups", c.groups.map((g) => ({ id: g, label: g })), state.groups);
    state.cellTypeOrder = c.cell_types.slice();
    renderMulti("#cell-types", c.cell_types.map((t) => ({ id: t, label: t })), state.cellTypes);
    setStatus(`${c.samples.length} samples · ${c.groups.length} groups · ${c.cell_types.length} cell types · ${c.gene_count.toLocaleString()} genes`);
  }
  function renderMulti(sel, items, store) {
    const host = $(sel); host.innerHTML = "";
    items.forEach((it) => {
      const b = document.createElement("span");
      b.className = "chip"; b.dataset.id = it.id;
      b.innerHTML = it.sub ? `${it.label}<em>${it.sub}</em>` : it.label;
      if (store.has(it.id)) b.classList.add("on");
      b.onclick = () => { b.classList.toggle("on"); b.classList.contains("on") ? store.add(it.id) : store.delete(it.id); autoRender(); };
      host.appendChild(b);
    });
  }
  const selectedCellTypes = () => state.cellTypeOrder.filter((t) => state.cellTypes.has(t));

  // ===================================================================== tabs
  $$("#tabbar .tab").forEach((t) => t.addEventListener("click", () => {
    if (t.classList.contains("on")) return;
    $$("#tabbar .tab").forEach((x) => x.classList.remove("on")); t.classList.add("on");
    state.tab = t.dataset.tab; syncTabVisibility();
    if (state.tab === "network" || state.gene || geneInput.value.trim()) render();   // network is gene-free; others need a gene
  }));
  function syncTabVisibility() {
    $$(".tab-controls").forEach((el) => { el.hidden = el.dataset.for !== state.tab; });
    $$("[data-for]").forEach((el) => {
      if (el.classList.contains("tab-controls")) return;
      el.style.display = el.dataset.for === state.tab ? "" : "none";
    });
    const gf = $("#gene-field"); if (gf) gf.style.display = (state.tab === "network") ? "none" : "";   // network is gene-free
    renderLegend();
  }

  // ===================================================================== gene autocomplete
  let suggestTimer = null;
  const geneInput = $("#gene-input");
  geneInput.addEventListener("input", (e) => {
    clearTimeout(suggestTimer);
    const q = e.target.value.trim();
    if (q.length < 1) { $("#gene-suggest").innerHTML = ""; return; }
    suggestTimer = setTimeout(async () => {
      const { matches } = await api("/api/genes?q=" + encodeURIComponent(q)).then((r) => r.json());
      const box = $("#gene-suggest"); box.innerHTML = "";
      matches.forEach((m) => {
        const d = document.createElement("div");
        d.innerHTML = `<b>${m.symbol || m.gene}</b>${m.symbol ? " <span class='muted'>" + m.gene + "</span>" : ""}`;
        d.onclick = () => { geneInput.value = m.symbol || m.gene; state.gene = m.gene; box.innerHTML = ""; loadJunctions(m.gene); render(); };
        box.appendChild(d);
      });
    }, 140);
  });
  geneInput.addEventListener("keydown", (e) => { if (e.key === "Enter") { $("#gene-suggest").innerHTML = ""; render(); } });
  async function loadJunctions(gene) {
    try {
      const { junctions } = await api("/api/junctions?gene=" + encodeURIComponent(gene)).then((r) => r.json());
      state.junctions.clear();
      renderMulti("#junctions", junctions.map((j) => ({ id: j, label: j })), state.junctions);
    } catch (_) { $("#junctions").innerHTML = ""; }
  }

  // ===================================================================== controls
  $("#combine-seg").addEventListener("click", (e) => {
    const b = e.target.closest(".seg-btn"); if (!b) return;
    $$("#combine-seg .seg-btn").forEach((x) => x.classList.remove("on")); b.classList.add("on");
    state.combineBy = b.dataset.v; if (state.gene || geneInput.value.trim()) render();
  });
  $("#panel-seg").addEventListener("click", (e) => {
    const b = e.target.closest(".seg-btn"); if (!b) return;
    $$("#panel-seg .seg-btn").forEach((x) => x.classList.remove("on")); b.classList.add("on");
    state.panelBy = b.dataset.v; if (state.gene || geneInput.value.trim()) render();
  });
  $("#net-panel-seg").addEventListener("click", (e) => {
    const b = e.target.closest(".seg-btn"); if (!b) return;
    $$("#net-panel-seg .seg-btn").forEach((x) => x.classList.remove("on")); b.classList.add("on");
    state.netPanelBy = b.dataset.v; updateNetPanelNote(); if (state.tab === "network") renderNetwork();
  });
  $("#sim").addEventListener("input", (e) => { $("#sim-val").textContent = e.target.value; autoRender(); });
  $("#celltype-filter").addEventListener("input", (e) => {
    const f = e.target.value.toLowerCase();
    $$("#cell-types .chip").forEach((c) => { c.style.display = c.textContent.toLowerCase().includes(f) ? "" : "none"; });
  });
  $("#ct-all").onclick = () => { state.cellTypeOrder.forEach((t) => state.cellTypes.add(t)); refreshChips("#cell-types", state.cellTypes); autoRender(); };
  $("#ct-none").onclick = () => { state.cellTypes.clear(); refreshChips("#cell-types", state.cellTypes); autoRender(); };
  $("#grp-all").onclick = () => { state.catalog.groups.forEach((g) => state.groups.add(g)); refreshChips("#groups", state.groups); autoRender(); };
  function refreshChips(sel, store) { $$(sel + " .chip").forEach((c) => c.classList.toggle("on", store.has(c.dataset.id))); }
  // every threshold/option control auto-renders on change (replaces the removed Render button)
  ["#maxreads", "#msf", "#strategy", "#mincount", "#maxiso", "#mode", "#introns", "#rownorm"].forEach((sel) => {
    const el = $(sel); if (el) el.addEventListener("change", autoRender);
  });
  $("#zoom-in").onclick = () => view && view.svg.transition().duration(160).call(view.zoom.scaleBy, 1.6);
  $("#zoom-out").onclick = () => view && view.svg.transition().duration(160).call(view.zoom.scaleBy, 1 / 1.6);
  $("#zoom-reset").onclick = () => view && view.svg.transition().duration(220).call(view.zoom.transform, d3.zoomIdentity);
  window.addEventListener("resize", () => { if (state.last) draw(state.last); });

  // ===================================================================== query dispatch
  function render() {
    if (state.tab === "network") return renderNetwork();   // network view is gene-free (uses contexts), check first
    const gene = state.gene || geneInput.value.trim();
    if (!gene) { setStatus("enter a gene"); return; }
    return state.tab === "molecule" ? renderReads(gene) : state.tab === "lines" ? renderLines(gene) : renderHeatmap(gene);
  }

  async function renderHeatmap(gene) {
    const body = {
      gene, samples: [...state.samples], groups: [...state.groups], cell_types: selectedCellTypes(),
      combine_by: state.combineBy,
      cluster_similarity_threshold: +$("#sim").value, min_split_fraction: +$("#msf").value,
      min_count: +$("#mincount").value, max_isoforms: +$("#maxiso").value,
      cluster_strategy: $("#strategy").value, cluster_mode: $("#mode").value,
      filter_junctions: [...state.junctions], include_introns: $("#introns").checked,
    };
    if (state.combineBy !== "group" && body.cell_types.length === 0) { setStatus("select cell types for this column mode"); return; }
    await run("/api/isoforms", body, gene, "heatmap", false);
  }

  async function renderReads(gene) {
    const cell_types = selectedCellTypes();
    if (!cell_types.length) { setStatus("select cell types to build the read pileup"); return; }
    const body = {
      gene, cell_types, conditions: state.groups.size ? [...state.groups] : null,
      max_isoforms: +$("#maxreads").value || 300, panel_by: state.panelBy,
    };
    await run("/api/reads", body, gene, "reads", true);
  }

  // Isoform-usage LINE view: reuse /api/isoforms with covariate×cell-type columns (so every (group,cell-type)
  // is a column); usage = isoform reads / column total. Same default filters as the heatmap tab.
  async function renderLines(gene) {
    const body = {
      gene, samples: [...state.samples], groups: [...state.groups], cell_types: selectedCellTypes(),
      combine_by: "cell_type_x_covariate",
      cluster_similarity_threshold: +$("#sim").value, min_split_fraction: +$("#msf").value,
      min_count: +$("#mincount").value, max_isoforms: +$("#maxiso").value,
      cluster_strategy: $("#strategy").value, cluster_mode: $("#mode").value,
      filter_junctions: [...state.junctions], include_introns: $("#introns").checked,
    };
    if (body.cell_types.length === 0) { setStatus("select cell types for the line view"); return; }
    await run("/api/isoforms", body, gene, "lines", false);
  }

  // ---- Interaction network tab (PPI/PDI contrast) ----
  let _netMeta = null, _netT = null, _netWired = false;
  async function netInit() {
    if (!_netMeta || _netMeta.available === false) {   // retry a transient meta failure rather than caching it
      try { _netMeta = await api("/api/network/meta").then((r) => r.json()); } catch (e) { _netMeta = { available: false }; }
    }
    const a = $("#net-ctx-a"), b = $("#net-ctx-b");
    if (a && b && _netMeta.contexts && _netMeta.contexts.length) {
      const opts = _netMeta.contexts.map((c) => `<option value="${c}">${c.replace("||", " · ")}</option>`).join("");
      a.innerHTML = opts; b.innerHTML = opts; netDefaults();
    }
    if (!_netWired) {
      const et = $("#net-edgetype"); if (et) { et.value = "PDI"; et.onchange = () => renderNetwork(); }   // default PDI
      ["#net-genes", "#net-thresh"].forEach((s) => { const el = $(s); if (el) el.oninput = () => { clearTimeout(_netT); _netT = setTimeout(renderNetwork, 350); }; });
      updateNetPanelNote();
      _netWired = true;   // the shared cell-type/covariate chips re-render via autoRender(); no Contrast button
    }
  }
  // ensure the shared chips carry a usable default for the contrast (default cell types + both covariates)
  function netDefaults() {
    if (!selectedCellTypes().length) { DEFAULT_CELLS.forEach((c) => state.cellTypes.add(c)); refreshChips("#cell-types", state.cellTypes); }
    if (state.groups.size < 2 && state.catalog) { (state.catalog.groups || []).forEach((g) => state.groups.add(g)); refreshChips("#groups", state.groups); }
  }
  function updateNetPanelNote() {
    const el = $("#net-panel-note"); if (!el) return;
    el.textContent = state.netPanelBy === "cell_type"
      ? "Windows = the first 2 selected cell types; covariates aggregated (chips below)"
      : "Windows = the first 2 selected covariates; cell types summed (chips below)";
  }
  async function renderNetwork() {
    if (!_netMeta || _netMeta.available === false) await netInit();
    if (!_netMeta || _netMeta.available === false) { setStatus("no interaction network precomputed for this dataset (run precompute_isvweb)"); return; }
    netDefaults();
    const cts = selectedCellTypes(), gs = [...state.groups], byCell = state.netPanelBy === "cell_type";
    if (byCell) {
      if (gs.length < 1) { setStatus("select at least 1 covariate"); return; }
      if (cts.length < 2) { setStatus("select 2 cell types to contrast (one per window)"); return; }
    } else {
      if (!cts.length) { setStatus("select cell types to build the interaction networks"); return; }
      if (gs.length < 2) { setStatus("select 2 covariates to contrast (e.g. young + AML-NPM1)"); return; }
    }
    const genes = ($("#net-genes").value || "").split(",").map((s) => s.trim()).filter(Boolean);
    const body = { cell_types: cts, groups: gs, threshold: +$("#net-thresh").value || 1,
                   genes: genes.length ? genes : null, edge_type: $("#net-edgetype").value, panel_by: state.netPanelBy };
    const title = (byCell ? cts : gs).slice(0, 2).join(" vs ");
    await run("/api/network", body, title, "network", false);
  }

  async function run(url, body, gene, mode, slow) {
    const msg = slow ? "building read-level pileup for " + gene + "…" : "querying " + gene + "…";
    setStatus(slow ? "generating read-level view for " + gene + " (first time for a selection can take a few seconds)…" : "querying " + gene + "…");
    startLoading(msg);
    const t0 = performance.now();
    try {
      const res = await api(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) }).then((r) => r.json());
      res._mode = mode; state.last = res; draw(res);
      setStatus(`rendered in ${Math.round(performance.now() - t0)} ms`);
    } catch (e) { setStatus("error: " + e.message); }
    finally { stopLoading(); }
  }

  // ----- loading overlay (animated DNA helix); delayed so instant/cached responses don't flash it -----
  let _loadTimer = null;
  function startLoading(msg) {
    clearTimeout(_loadTimer);
    _loadTimer = setTimeout(() => { $("#loading-txt").textContent = msg || "rendering…"; $("#loading").classList.remove("hidden"); }, 250);
  }
  function stopLoading() { clearTimeout(_loadTimer); $("#loading").classList.add("hidden"); }
  (function buildDNA() {
    const d = $("#dna"); if (!d) return;
    for (let i = 0; i < 18; i++) { const s = document.createElement("span"); s.style.setProperty("--d", (i * 0.07) + "s"); d.appendChild(s); }
  })();

  // ===================================================================== geometry helpers
  function mergeSegments(segs) {
    const es = (segs || []).filter((s) => s.start != null && s.end != null && (s.type || "E").toUpperCase().startsWith("E"))
      .map((s) => ({ a: Math.min(s.start, s.end), b: Math.max(s.start, s.end), label: (s.exon_id || s.label || "").split(".")[0] }))
      .sort((p, q) => p.a - q.a);
    const out = [];
    es.forEach((s) => {
      const last = out[out.length - 1];
      if (last && s.a <= last.b + 2 && s.label === last.label) last.b = Math.max(last.b, s.b);
      else if (last && s.a <= last.b) { last.b = Math.max(last.b, s.b); if (s.label && !last.label) last.label = s.label; }
      else out.push({ a: s.a, b: s.b, label: s.label });
    });
    const introns = (segs || []).filter((s) => s.start != null && s.end != null && (s.type || "").toUpperCase().startsWith("I"))
      .map((s) => ({ a: Math.min(s.start, s.end), b: Math.max(s.start, s.end) }));
    return { exons: out, introns };
  }
  function extentOf(items, gm) {
    let lo = Infinity, hi = -Infinity;
    if (gm && gm.extent) { lo = Math.min(lo, gm.extent[0]); hi = Math.max(hi, gm.extent[1]); }
    (gm && gm.blocks || []).forEach((b) => { lo = Math.min(lo, b.start); hi = Math.max(hi, b.end); });
    items.forEach((it) => (it.merged.exons.concat(it.merged.introns)).forEach((s) => { lo = Math.min(lo, s.a); hi = Math.max(hi, s.b); }));
    if (!isFinite(lo)) { lo = 0; hi = 1; }
    if (hi <= lo) hi = lo + 1;
    return [lo, hi];
  }

  // shared zoom installer (horizontal only; hover + right-click preserved -- no covering overlay)
  function installZoom(svg, x0, trackX0, trackX1, height, drawX) {
    let raf = false, lastT = null;
    const onZoom = (t) => { lastT = t; if (raf) return; raf = true; requestAnimationFrame(() => { raf = false; drawX(lastT.rescaleX(x0)); }); };
    const zoom = d3.zoom().scaleExtent([1, 600])
      .extent([[trackX0, 0], [trackX1, height]]).translateExtent([[trackX0, 0], [trackX1, height]])
      .filter((ev) => {
        if (ev.ctrlKey) return false;
        if (d3.pointer(ev, svg.node())[0] < trackX0) return false;
        return ev.type === "wheel" || ev.button === 0 || ev.button == null;
      })
      .on("zoom", (ev) => onZoom(ev.transform));
    svg.call(zoom).on("dblclick.zoom", null);
    drawX(x0);
    view = { svg, zoom, x0, drawX };
  }

  // genomic ruler + reference gene-model track (exon boxes + intron spine + strand arrows + exon labels)
  function drawAxis(g, xz, gm, lo, hi, trackX0, trackX1, rulerBaseY, refTopY, refH, withLabels) {
    const ticks = xz.ticks(Math.max(4, Math.round((trackX1 - trackX0) / 120)));
    g.append("line").attr("class", "ruler-axis").attr("x1", trackX0).attr("x2", trackX1).attr("y1", rulerBaseY).attr("y2", rulerBaseY);
    ticks.forEach((tk) => {
      const px = xz(tk); if (px < trackX0 - 0.5 || px > trackX1 + 0.5) return;
      g.append("line").attr("class", "ruler-tick").attr("x1", px).attr("x2", px).attr("y1", rulerBaseY - 5).attr("y2", rulerBaseY);
      g.append("text").attr("class", "ruler-lbl").attr("x", px).attr("y", rulerBaseY - 8).text(fmt(tk));
    });
    const cy = refTopY + refH / 2;
    const refLo = Math.max(trackX0, xz(lo)), refHi = Math.min(trackX1, xz(hi));
    g.append("line").attr("class", "ref-spine").attr("x1", refLo).attr("x2", refHi).attr("y1", cy).attr("y2", cy);
    strandArrows(g, refLo, refHi, cy, gm && gm.strand);
    let lastLbl = null;
    (gm && gm.blocks || []).forEach((b) => {
      const xa = xz(b.start), xb = xz(b.end);
      if (xb < trackX0 || xa > trackX1) return;
      g.append("rect").attr("class", "ref-exon").attr("x", xa).attr("y", refTopY + 3).attr("width", Math.max(1.6, xb - xa)).attr("height", refH - 6)
        .on("mousemove", (ev) => showTip(ev, `<b>${b.exon_id || "exon"}</b><br>${gm.chrom ? gm.chrom + ":" : ""}${fmt(b.start)}–${fmt(b.end)}<br><span class="muted">reference gene model</span>`))
        .on("mouseleave", hideTip);
      if (withLabels) {
        const lbl = (b.exon_id || "").split(".")[0];
        if (lbl && lbl !== lastLbl && (xb - xa) > 7) {
          g.append("text").attr("class", "exon-lbl").attr("x", (xa + xb) / 2).attr("y", refTopY + refH + 10).text(lbl);
          lastLbl = lbl;
        }
      }
    });
  }
  function strandArrows(g, x1, x2, y, strand) {
    if (!strand || x2 - x1 < 40) return;
    const dir = strand === "-" ? -1 : 1, step = 46;
    for (let px = x1 + 24; px < x2 - 6; px += step)
      g.append("path").attr("class", "strand-arrow").attr("d", `M${px - 3 * dir},${y - 3} L${px},${y} L${px - 3 * dir},${y + 3}`);
  }

  // ===================================================================== draw dispatch
  function draw(res) {
    const plot = $("#plot"); plot.innerHTML = ""; view = null;
    if (res._mode === "network") return drawNetwork(res);
    titleBlock(res);
    return res._mode === "reads" ? drawReads(res) : res._mode === "lines" ? drawLines(res) : drawHeatmap(res);
  }

  // ============================ READ-LEVEL PILEUP (molecule view) ============================
  const R = { LEFT: 152, RIGHT: 24, READ_PX: 1.6, ROW_MIN: 0.9, ROW_CAP: 60, MAXPANEL: 340, PANEL_GAP: 14,
              PAD: 8, HDR: 22, RULER_H: 30, REF_H: 26, MAX_ROWS: 4000, RMARGIN: 16 };
  function drawReads(res) {
    const plot = $("#plot");
    const panels = (res.panels || []).filter((p) => p.n_molecules > 0);
    if (!panels.length) { plot.innerHTML = "<p class='empty'>No molecules for this selection. Pick cell types (and covariates) and Render.</p>"; $("#zoom-tools").classList.add("hidden"); return; }
    $("#zoom-tools").classList.remove("hidden");
    const gm = res.gene_model || {};

    // prep each panel: merge segments, sort by cluster then count desc, cap rows
    let truncated = 0;
    panels.forEach((p) => {
      p.mol = p.molecules.map((m) => ({ ...m, merged: mergeSegments(m.exon_segments) }))
        .sort((a, b) => (a.cluster_index - b.cluster_index) || (b.count - a.count) || (a.isoform_id < b.isoform_id ? -1 : 1));
      if (p.mol.length > R.MAX_ROWS) { truncated += p.mol.length - R.MAX_ROWS; p.mol = p.mol.slice(0, R.MAX_ROWS); }
    });
    const allMol = panels.flatMap((p) => p.mol);
    const [lo, hi] = extentOf(allMol, gm);

    const containerW = plot.clientWidth || 1100;
    const width = Math.max(containerW - 4, 900);
    const trackX0 = R.LEFT, trackX1 = width - R.RIGHT - R.RMARGIN;
    const x0 = d3.scaleLinear().domain([lo, hi]).range([trackX0, trackX1]);

    // vertical layout: stack panels, then axis at the bottom. Each panel's rows are compressed to fit
    // MAXPANEL so dense pileups stay compact (greater vertical compression) while small panels stay readable.
    panels.forEach((p) => {
      const raw = p.mol.map((m) => Math.min(R.ROW_CAP, Math.max(R.ROW_MIN, m.count * R.READ_PX)));
      const rawSum = raw.reduce((a, b) => a + b, 0) || 1;
      // target panel height: compress dense pileups to <= MAXPANEL, but EXPAND sparse ones so every distinct
      // isoform cluster gets >= ~13px of vertical space -- otherwise (low-read genes like LAT) the per-cluster
      // labels collapse into one illegible stack and reads are too thin to hover individually.
      const nClusters = new Set(p.mol.map((m) => m.cluster_index)).size;
      const targetH = Math.min(R.MAXPANEL, Math.max(rawSum, nClusters * 13));
      const scale = targetH / rawSum;
      let h = 0; p._rows = [];
      p.mol.forEach((m, i) => { const hh = Math.max(0.55, raw[i] * scale); p._rows.push({ m, dy: h, h: hh }); h += hh; });
      p._bodyH = Math.max(18, h);
      // per-structure (isoform cluster) read totals IN THIS PANEL: many molecules collapse to one
      // structure, so the tooltip reports the cluster's summed reads (not the single hovered molecule's).
      const clTot = new Map(), clMol = new Map();
      p.mol.forEach((m) => {
        clTot.set(m.cluster_index, (clTot.get(m.cluster_index) || 0) + (+m.count || 0));
        clMol.set(m.cluster_index, (clMol.get(m.cluster_index) || 0) + 1);
      });
      p._clTot = clTot; p._clMol = clMol;
    });
    let y = 10;
    panels.forEach((p) => { p._top = y; p._y = y + R.HDR; y += R.HDR + p._bodyH + R.PAD + R.PANEL_GAP; });
    const axisTop = y + 2;
    const height = axisTop + R.RULER_H + R.REF_H + 22;
    const rulerBaseY = axisTop + R.RULER_H - 4;
    const refTopY = axisTop + R.RULER_H;

    const svg = d3.select(plot).append("svg").attr("width", width).attr("height", height);
    svg.append("defs").append("clipPath").attr("id", "rclip").append("rect")
      .attr("x", trackX0).attr("y", 0).attr("width", trackX1 - trackX0).attr("height", height);
    const gStatic = svg.append("g");
    // exon/backbone layer is drawn ON TOP; make it transparent to the mouse so hovering an exon falls through
    // to the per-panel hit overlay below (so exons, labels, and gaps all trigger the isoform tooltip).
    const gx = svg.append("g").attr("clip-path", "url(#rclip)").style("pointer-events", "none");

    // static per-panel chrome: box, top-left title, sparse molecule labels, mouse overlay
    panels.forEach((p) => {
      const top = p._top, boxH = R.HDR + p._bodyH + R.PAD;
      gStatic.append("rect").attr("class", "panel-box").attr("x", 6).attr("y", top).attr("width", width - 12).attr("height", boxH).attr("rx", 7);
      gStatic.append("rect").attr("class", "panel-hdr-band").attr("x", 6).attr("y", top).attr("width", width - 12).attr("height", R.HDR).attr("rx", 7);
      gStatic.append("text").attr("class", "panel-label").attr("x", 16).attr("y", top + 15)
        .html(`${p.condition}<tspan class="panel-sub-inline"> · ${p.n_molecules.toLocaleString()} reads</tspan>`);
      // one left label per isoform cluster, VERTICALLY CENTERED in that cluster's block of reads
      const clExtent = new Map();   // cluster_index -> [yTop, yBot, representative molecule]
      p._rows.forEach((r) => {
        const e = clExtent.get(r.m.cluster_index);
        if (!e) clExtent.set(r.m.cluster_index, [r.dy, r.dy + r.h, r.m]);
        else { e[0] = Math.min(e[0], r.dy); e[1] = Math.max(e[1], r.dy + r.h); }
      });
      clExtent.forEach(([yTop, yBot, m]) => {
        gStatic.append("text").attr("class", "mol-label " + (m.known ? "known" : "novel"))
          .attr("x", R.LEFT - 8).attr("y", p._y + (yTop + yBot) / 2).attr("dominant-baseline", "central")
          .text(truncate(m.final_isoform_id || m.isoform_id, 20));   // labelled by its final collapsed isoform
      });
      // one mouse overlay per panel, spanning the LABEL column + the exon track. Rows are contiguous, so map
      // EVERY y to the row whose bottom is just past it (no gaps) and extend well below the last read. Hover
      // shows the isoform tooltip (anywhere: label, exons, gaps); LEFT-click (or right-click) opens the menu.
      const rowAt = (my) => p._rows.find((rr) => my < rr.dy + rr.h) || p._rows[p._rows.length - 1];
      const openMenu = function (ev) {
        ev.preventDefault(); ev.stopPropagation();
        const r = rowAt(d3.pointer(ev, this)[1]);
        if (r) showCtx(ev, r.m, p.mol.map((m) => m.final_isoform_id || m.isoform_id));
      };
      gStatic.append("rect").attr("class", "panel-hit").attr("x", 6).attr("y", p._y)
        .attr("width", trackX1 - 6).attr("height", p._bodyH + R.PAD + R.PANEL_GAP)
        .attr("fill", "transparent").style("cursor", "pointer")
        .on("mousemove", function (ev) {
          const r = rowAt(d3.pointer(ev, this)[1]);
          if (r) showTip(ev, molTip(r.m, p.condition, p));
        })
        .on("mouseleave", hideTip)
        .on("click", openMenu)
        .on("contextmenu", (ev) => ev.preventDefault());
    });

    function drawX(xz) {
      gx.selectAll("*").remove();
      panels.forEach((p) => {
        p._rows.forEach((r) => {
          const m = r.m, yTop = p._y + r.dy, color = colorFor(m);
          const exs = m.merged.exons;
          const eh = Math.max(0.7, r.h > 3 ? r.h - 0.8 : r.h);
          const ey = yTop + (r.h - eh) / 2;
          if (exs.length) {
            const a = xz(Math.min(...exs.map((e) => e.a))), b = xz(Math.max(...exs.map((e) => e.b)));
            gx.append("line").attr("class", "backbone").attr("x1", a).attr("x2", b).attr("y1", yTop + r.h / 2).attr("y2", yTop + r.h / 2).attr("stroke", color).attr("opacity", 0.55);
          }
          m.merged.introns.forEach((s) => {
            const xa = xz(s.a), xb = xz(s.b);
            gx.append("rect").attr("x", xa).attr("y", yTop + r.h / 2 - eh * 0.3).attr("width", Math.max(1, xb - xa)).attr("height", Math.max(1, eh * 0.6)).attr("fill", color).attr("opacity", 0.45);
          });
          exs.forEach((s) => {
            const xa = xz(s.a), xb = xz(s.b);
            gx.append("rect").attr("class", "mol-exon").attr("x", xa).attr("y", ey).attr("width", Math.max(1.4, xb - xa)).attr("height", eh).attr("fill", color);
          });
        });
      });
      drawAxis(gx, xz, gm, lo, hi, trackX0, trackX1, rulerBaseY, refTopY, R.REF_H, true);
      gx.append("text").attr("class", "axis-chrom").attr("x", trackX0).attr("y", refTopY + R.REF_H + 11).text(gm.chrom || "");
    }
    installZoom(svg, x0, trackX0, trackX1, height, drawX);
    if (truncated) setStatus(`note: ${truncated.toLocaleString()} low-rank molecules hidden (cap ${R.MAX_ROWS}/panel) — raise "Max molecules / panel"`);
    renderLegend(res);
  }
  function molTip(m, cond, p) {
    const fid = m.final_isoform_id || m.isoform_id;
    const det = m.detections || [];
    // total reads assigned to this isoform STRUCTURE in this cell-state/condition panel, from the
    // authoritative counts matrix (server: panel.structure_reads). The molecule pileup is top-N capped,
    // so summing displayed molecules undercounts -- fall back to that sum only if the total is absent.
    const sr = (p && p.structure_reads) ? p.structure_reads[fid] : undefined;
    const structTotal = (sr != null) ? sr
      : ((p && p._clTot && p._clTot.has(m.cluster_index)) ? p._clTot.get(m.cluster_index) : m.count);
    let where;
    if (det.length) {
      const lines = det.slice(0, 8).map((d) =>
        `&nbsp;&nbsp;${d.cell_type} · <span class="muted">${d.sample}</span>${d.count != null ? ` (${fmtCount(d.count)})` : ""}`);
      where = `<br>detected in:<br>${lines.join("<br>")}`
        + (det.length > 8 ? `<br><span class="muted">&nbsp;&nbsp;+${det.length - 8} more</span>` : "");
    } else {
      where = m.sample ? `<br>sample: ${m.sample}` : "";
    }
    return `<b>${fid}</b> <span class="muted">${m.known ? "known" : "novel"}</span>`
      + (fid !== m.isoform_id ? `<br><span class="muted">read: ${m.isoform_id}</span>` : "")
      + `<br>${cond} · isoform group ${m.cluster_index}`
      + where
      + `<br>reads (this structure · ${cond}): <b>${fmtCount(structTotal)}</b>`
      + `<br>protein length: ${m.protein_length != null ? m.protein_length + " aa" : "n/a"}`
      + (m.nmd_status ? `<br>NMD: ${m.nmd_status}` : "");
  }

  // ============================ HEATMAP (cell type × covariate) ============================
  const H = { ROW_H: 20, ROW_GAP: 6, LEFT: 210, CELL_W: 22, GAP: 30, RULER_H: 30, REF_H: 26, BLOCK_LBL: 20, COLHDR_H: 92, RMARGIN: 28 };
  function drawHeatmap(res) {
    const plot = $("#plot");
    if (!res.isoforms || !res.isoforms.length) { plot.innerHTML = "<p class='empty'>No isoforms for this selection.</p>"; $("#zoom-tools").classList.add("hidden"); return; }
    $("#zoom-tools").classList.remove("hidden");
    const cols = res.columns || [];
    const isos = res.isoforms.map((iso) => ({ ...iso, merged: mergeSegments(iso.exon_segments) }));
    const gm = res.gene_model || {};
    const [lo, hi] = extentOf(isos, gm);

    // covariate "block" labels only when they ADD info (cell-type × covariate: group=covariate spans the
    // cell-type columns). In the covariate-only view group == column label, so the block label just
    // duplicates the column header below -> suppress it.
    const showBlocks = cols.some((c) => c.group && c.group !== c.label);
    const annoHdrH = (showBlocks ? H.BLOCK_LBL : 0) + H.COLHDR_H;
    const HEADER = Math.max(annoHdrH, H.RULER_H + H.REF_H, 40);
    const annoW = Math.max(H.CELL_W, cols.length * H.CELL_W);
    const exprX0 = H.LEFT, trackX0 = H.LEFT + annoW + H.GAP;
    const containerW = plot.clientWidth || 1100;
    const width = Math.max(containerW - 4, trackX0 + 480);
    const trackX1 = width - H.RMARGIN;
    const bodyH = isos.length * (H.ROW_H + H.ROW_GAP);
    const height = HEADER + bodyH + 22;
    const refTop = HEADER - H.REF_H, rulerBase = refTop - 6;
    const x0 = d3.scaleLinear().domain([lo, hi]).range([trackX0, trackX1]);

    const svg = d3.select(plot).append("svg").attr("width", width).attr("height", height);
    svg.append("defs").append("clipPath").attr("id", "hclip").append("rect").attr("x", trackX0).attr("y", 0).attr("width", trackX1 - trackX0).attr("height", height);
    const gStatic = svg.append("g");
    const gx = svg.append("g").attr("clip-path", "url(#hclip)");

    isos.forEach((iso, r) => {
      const y = HEADER + r * (H.ROW_H + H.ROW_GAP);
      const row = gStatic.append("g").attr("transform", `translate(0,${y})`);
      const exportIso = (ev) => { ev.preventDefault(); ev.stopPropagation(); showCtx(ev, iso, isos.map((i) => i.isoform_id), iso.cluster_id, res.isoforms); };
      row.append("rect").attr("class", "rowhit").attr("x", 0).attr("y", 0).attr("width", width).attr("height", H.ROW_H).attr("fill", "transparent").style("cursor", "pointer")
        .on("mousemove", (ev) => { if (ev.target.classList.contains("rowhit")) showTip(ev, isoTip(iso)); })
        .on("mouseleave", hideTip)
        .on("click", exportIso)             // LEFT-click opens the sequence-export menu (like the molecule view)
        .on("contextmenu", (ev) => ev.preventDefault());
      row.append("text").attr("class", "iso-label " + (iso.known ? "known" : "novel")).attr("x", 12).attr("y", H.ROW_H * 0.7).text(truncate(iso.isoform_id, 28));
    });
    drawHeatStrip(gStatic, isos, cols, exprX0, HEADER, showBlocks, bodyH);

    function drawX(xz) {
      gx.selectAll("*").remove();
      drawAxis(gx, xz, gm, lo, hi, trackX0, trackX1, rulerBase, refTop, H.REF_H, false);
      isos.forEach((iso, r) => {
        const y = HEADER + r * (H.ROW_H + H.ROW_GAP);
        const g = gx.append("g").attr("transform", `translate(0,${y})`);
        const color = iso.known ? KNOWN : NOVEL, exs = iso.merged.exons;
        if (exs.length) {
          const a = xz(Math.min(...exs.map((e) => e.a))), b = xz(Math.max(...exs.map((e) => e.b)));
          g.append("line").attr("class", "backbone").attr("x1", a).attr("x2", b).attr("y1", H.ROW_H / 2).attr("y2", H.ROW_H / 2).attr("stroke", color).attr("opacity", 0.5);
        }
        const exportIso = (ev) => { ev.preventDefault(); ev.stopPropagation(); showCtx(ev, iso, isos.map((i) => i.isoform_id), iso.cluster_id, res.isoforms); };
        iso.merged.introns.forEach((s) => {
          const xa = xz(s.a), xb = xz(s.b);
          g.append("rect").attr("class", "seg intron").attr("x", xa).attr("y", H.ROW_H / 2 - 2.5).attr("width", Math.max(1.5, xb - xa)).attr("height", 5).attr("fill", color).attr("opacity", 0.5).style("cursor", "pointer")
            .on("mousemove", (ev) => showTip(ev, `<b>intron retention</b><br>${fmt(s.a)}–${fmt(s.b)}<br><span class="muted">${iso.isoform_id}</span>`)).on("mouseleave", hideTip)
            .on("click", exportIso).on("contextmenu", (ev) => ev.preventDefault());
        });
        // UCSC genome-browser style: coding (CDS) exon portions full-height (tall), non-coding (5'/3'-UTR)
        // portions half-height (thin). cds_min/cds_max = the isoform's coding genomic span (null => no CDS).
        const exH = H.ROW_H - 6, exY = 3;
        const utrH = Math.max(3, exH * 0.5), utrY = exY + (exH - utrH) / 2;
        const cdsLo = iso.cds_min, cdsHi = iso.cds_max;
        const exonRect = (s, ga, gb, y, h) => g.append("rect").attr("class", "seg exon")
          .attr("x", xz(ga)).attr("y", y).attr("width", Math.max(1.5, xz(gb) - xz(ga))).attr("height", h).attr("fill", color).style("cursor", "pointer")
          .on("mousemove", (ev) => showTip(ev, `<b>${s.label || "exon"}</b><br>${gm.chrom ? gm.chrom + ":" : ""}${fmt(s.a)}–${fmt(s.b)}<br>${iso.protein_length != null ? iso.protein_length + " aa · " : ""}<span class="muted">${iso.isoform_id}</span>`))
          .on("mouseleave", hideTip).on("click", exportIso).on("contextmenu", (ev) => ev.preventDefault());
        exs.forEach((s) => {
          if (cdsLo != null && cdsHi != null) {
            exonRect(s, s.a, s.b, utrY, utrH);                          // UTR (thin) across the whole exon
            const ca = Math.max(s.a, cdsLo), cb = Math.min(s.b, cdsHi);
            if (cb > ca) exonRect(s, ca, cb, exY, exH);                 // CDS (tall) over the coding overlap
          } else {
            exonRect(s, s.a, s.b, exY, exH);                           // no CDS info -> full height (unchanged)
          }
        });
      });
    }
    installZoom(svg, x0, trackX0, trackX1, height, drawX);
    renderLegend(res);
  }

  // ============================ ISOFORM USAGE LINE PLOTS (per covariate × cell type) ============================
  // Top panels: one stacked line panel per covariate (condition); x = cell types, y = % isoform usage, one
  // colored line per isoform. Bottom: the SAME structure layout as the heatmap tab (linear genomic axis +
  // reference track + per-isoform exon rows with UCSC CDS-tall/UTR-thin), color-matched to the lines, zoomable.
  const LP = { LEFT: 214, RMARGIN: 28, YAX_W: 34, PANEL_H: 132, PANEL_GAP: 26, TITLE_H: 20, XLAB_H: 50,
               RULER_H: 30, REF_H: 26, ROW_H: 20, ROW_GAP: 6, TOPN: 10 };
  function drawLines(res) {
    const plot = $("#plot");
    if (!res.isoforms || !res.isoforms.length) { plot.innerHTML = "<p class='empty'>No isoforms for this selection.</p>"; $("#zoom-tools").classList.add("hidden"); return; }
    const cols = res.columns || [];
    if (!cols.some((c) => c.group && c.group !== c.label)) {
      plot.innerHTML = "<p class='empty'>The line view needs covariate × cell-type columns — select cell types and ≥1 covariate group, then Render.</p>"; $("#zoom-tools").classList.add("hidden"); return; }
    $("#zoom-tools").classList.remove("hidden");
    const gm = res.gene_model || {};
    const isosAll = res.isoforms.map((iso) => ({ ...iso, merged: mergeSegments(iso.exon_segments) }));
    const colTotal = {}; cols.forEach((c) => { colTotal[c.key] = isosAll.reduce((s, iso) => s + (iso.expression[c.key] || 0), 0); });
    const tot = (iso) => cols.reduce((s, c) => s + (iso.expression[c.key] || 0), 0);
    const isos = isosAll.slice().sort((a, b) => tot(b) - tot(a)).slice(0, LP.TOPN);   // top isoforms by reads
    const color = new Map(isos.map((iso, i) => [iso.isoform_id, CLUSTER_COLORS[i % CLUSTER_COLORS.length]]));
    const usage = (iso, c) => (colTotal[c.key] ? 100 * (iso.expression[c.key] || 0) / colTotal[c.key] : 0);
    const conds = []; cols.forEach((c) => { const g = c.group || "all"; const last = conds[conds.length - 1]; if (!last || last.group !== g) conds.push({ group: g, cols: [c] }); else last.cols.push(c); });
    const [lo, hi] = extentOf(isos, gm);

    const containerW = plot.clientWidth || 1100;
    const width = Math.max(containerW - 4, 980);
    const plotX0 = LP.LEFT, plotX1 = width - LP.RMARGIN, trackX0 = LP.LEFT, trackX1 = plotX1;
    let y = 12;
    conds.forEach((cd) => { cd._titleY = y; cd._plotTop = y + LP.TITLE_H; cd._plotBot = cd._plotTop + LP.PANEL_H; y = cd._plotBot + LP.XLAB_H + LP.PANEL_GAP; });
    const structTop = y + 4, rulerBase = structTop + LP.RULER_H - 6, refTop = structTop + LP.RULER_H;
    const rowsTop = refTop + LP.REF_H + 10;
    const height = rowsTop + isos.length * (LP.ROW_H + LP.ROW_GAP) + 22;
    const x0 = d3.scaleLinear().domain([lo, hi]).range([trackX0, trackX1]);

    const svg = d3.select(plot).append("svg").attr("width", width).attr("height", height);
    svg.append("defs").append("clipPath").attr("id", "lclip").append("rect").attr("x", trackX0).attr("y", 0).attr("width", trackX1 - trackX0).attr("height", height);
    const gStatic = svg.append("g");
    const gx = svg.append("g").attr("clip-path", "url(#lclip)");

    // ---- one line panel per condition (covariate) ----
    conds.forEach((cd) => {
      const px0 = plotX0 + LP.YAX_W, px1 = plotX1;
      const xc = (i) => (cd.cols.length <= 1 ? (px0 + px1) / 2 : px0 + (i * (px1 - px0)) / (cd.cols.length - 1));
      const yv = (u) => cd._plotBot - (clamp(u, 0, 100) / 100) * LP.PANEL_H;
      const gp = gStatic.append("g");
      gp.append("text").attr("class", "panel-label").attr("x", plotX0).attr("y", cd._titleY + 14).text(cd.group);
      [0, 25, 50, 75, 100].forEach((u) => {
        const yy = yv(u);
        gp.append("line").attr("x1", px0).attr("x2", px1).attr("y1", yy).attr("y2", yy).attr("stroke", "#e6e9f0").attr("stroke-width", 1);
        gp.append("text").attr("x", px0 - 6).attr("y", yy + 3).attr("text-anchor", "end").attr("font-size", 9.5).attr("fill", "#8a94a4").text(u);
      });
      gp.append("text").attr("x", plotX0 + 2).attr("y", yv(50)).attr("transform", `rotate(-90 ${plotX0 + 2} ${yv(50)})`).attr("text-anchor", "middle").attr("font-size", 10).attr("fill", "#6b7686").text("% usage");
      cd.cols.forEach((c, i) => {
        const cx = xc(i);
        gp.append("text").attr("x", cx).attr("y", cd._plotBot + 12).attr("transform", `rotate(-40 ${cx} ${cd._plotBot + 12})`).attr("text-anchor", "end").attr("font-size", 10).attr("fill", "#6b7686").text(truncate(c.label, 16));
      });
      isos.forEach((iso) => {
        const col = color.get(iso.isoform_id);
        const pts = cd.cols.map((c, i) => [xc(i), yv(usage(iso, c))]);
        gp.append("path").attr("d", pts.map((p, i) => (i ? "L" : "M") + p[0].toFixed(1) + "," + p[1].toFixed(1)).join("")).attr("fill", "none").attr("stroke", col).attr("stroke-width", 2).attr("opacity", 0.95);
        pts.forEach((p, i) => gp.append("circle").attr("cx", p[0]).attr("cy", p[1]).attr("r", 2.7).attr("fill", col).style("cursor", "pointer")
          .on("mousemove", (ev) => showTip(ev, `<b>${iso.isoform_id}</b><br>${cd.group} · ${cd.cols[i].label}<br>usage: <b>${usage(iso, cd.cols[i]).toFixed(1)}%</b> · ${iso.expression[cd.cols[i].key] || 0} reads`)).on("mouseleave", hideTip));
      });
    });

    // ---- isoform structures (heatmap-tab layout; color-matched, zoomable) ----
    gStatic.append("text").attr("x", plotX0).attr("y", structTop - 3).attr("font-size", 11).attr("font-weight", 700).attr("fill", "#8a94a4").text("ISOFORM STRUCTURES");
    isos.forEach((iso, r) => {
      const yy = rowsTop + r * (LP.ROW_H + LP.ROW_GAP);
      const g = gStatic.append("g").attr("transform", `translate(0,${yy})`);
      const exportIso = (ev) => { ev.preventDefault(); ev.stopPropagation(); showCtx(ev, iso, isos.map((i) => i.isoform_id), iso.cluster_id, res.isoforms); };
      g.append("rect").attr("class", "rowhit").attr("x", 0).attr("y", 0).attr("width", width).attr("height", LP.ROW_H).attr("fill", "transparent").style("cursor", "pointer")
        .on("mousemove", (ev) => { if (ev.target.classList.contains("rowhit")) showTip(ev, isoTip(iso)); }).on("mouseleave", hideTip)
        .on("click", exportIso).on("contextmenu", (ev) => ev.preventDefault());
      g.append("rect").attr("x", 12).attr("y", LP.ROW_H / 2 - 4.5).attr("width", 9).attr("height", 9).attr("rx", 2).attr("fill", color.get(iso.isoform_id));
      g.append("text").attr("class", "iso-label " + (iso.known ? "known" : "novel")).attr("x", 26).attr("y", LP.ROW_H * 0.7).text(truncate(iso.isoform_id, 24));
    });
    function drawX(xz) {
      gx.selectAll("*").remove();
      drawAxis(gx, xz, gm, lo, hi, trackX0, trackX1, rulerBase, refTop, LP.REF_H, true);
      isos.forEach((iso, r) => {
        const yy = rowsTop + r * (LP.ROW_H + LP.ROW_GAP);
        const g = gx.append("g").attr("transform", `translate(0,${yy})`);
        const col = color.get(iso.isoform_id), exs = iso.merged.exons;
        const exportIso = (ev) => { ev.preventDefault(); ev.stopPropagation(); showCtx(ev, iso, isos.map((i) => i.isoform_id), iso.cluster_id, res.isoforms); };
        if (exs.length) {
          const a = xz(Math.min(...exs.map((e) => e.a))), b = xz(Math.max(...exs.map((e) => e.b)));
          g.append("line").attr("class", "backbone").attr("x1", a).attr("x2", b).attr("y1", LP.ROW_H / 2).attr("y2", LP.ROW_H / 2).attr("stroke", col).attr("opacity", 0.5);
        }
        iso.merged.introns.forEach((s) => {
          g.append("rect").attr("class", "seg intron").attr("x", xz(s.a)).attr("y", LP.ROW_H / 2 - 2.5).attr("width", Math.max(1.5, xz(s.b) - xz(s.a))).attr("height", 5).attr("fill", col).attr("opacity", 0.5).style("cursor", "pointer")
            .on("mousemove", (ev) => showTip(ev, `<b>intron retention</b><br>${fmt(s.a)}–${fmt(s.b)}<br><span class="muted">${iso.isoform_id}</span>`)).on("mouseleave", hideTip).on("click", exportIso).on("contextmenu", (ev) => ev.preventDefault());
        });
        const exH = LP.ROW_H - 6, exY = 3, utrH = Math.max(3, exH * 0.5), utrY = exY + (exH - utrH) / 2;
        const cdsLo = iso.cds_min, cdsHi = iso.cds_max;
        const exonRect = (s, ga, gb, ry, h) => g.append("rect").attr("class", "seg exon").attr("x", xz(ga)).attr("y", ry).attr("width", Math.max(1.5, xz(gb) - xz(ga))).attr("height", h).attr("fill", col).style("cursor", "pointer")
          .on("mousemove", (ev) => showTip(ev, `<b>${s.label || "exon"}</b><br>${gm.chrom ? gm.chrom + ":" : ""}${fmt(s.a)}–${fmt(s.b)}<br>${iso.protein_length != null ? iso.protein_length + " aa · " : ""}<span class="muted">${iso.isoform_id}</span>`)).on("mouseleave", hideTip).on("click", exportIso).on("contextmenu", (ev) => ev.preventDefault());
        exs.forEach((s) => {
          if (cdsLo != null && cdsHi != null) { exonRect(s, s.a, s.b, utrY, utrH); const ca = Math.max(s.a, cdsLo), cb = Math.min(s.b, cdsHi); if (cb > ca) exonRect(s, ca, cb, exY, exH); }
          else exonRect(s, s.a, s.b, exY, exH);
        });
      });
    }
    installZoom(svg, x0, trackX0, trackX1, height, drawX);
    renderLegend(res);
  }

  // ============================ INTERACTION NETWORK (PPI/PDI, side-by-side contrast) ============================
  // Two expression-supported subgraphs (ctxA | ctxB) on a SHARED force layout (so the same nodes sit in the
  // same place and differences pop). Source isoforms = diamonds (edge color = isoform); targets = circles;
  // node fill/size = per-context CPM; PDI = directed arrow, PPI = line. Mirrors the iso2function ego figures.
  let _fcoseReg = false;
  function ensureFcose() { if (_fcoseReg) return; _fcoseReg = true; try { if (window.cytoscapeFcose) cytoscape.use(window.cytoscapeFcose); } catch (e) {} }
  function setNetTitle(res) {
    const tb = $("#title-block"); if (!tb) return;   // replace the stale per-gene title with a network title (no BID)
    const et = (res.edge_type || "all").toUpperCase();
    tb.innerHTML = `<span class="tb-gene">Interactions</span> `
      + `<span class="tb-ensg">${res.A.ctx.replace("||", " · ")}  vs  ${res.B.ctx.replace("||", " · ")}</span> `
      + `<span class="pill">${et}</span> <span class="pill">${res.mode}</span> `
      + `<span class="pill">${res.A.edges.length} | ${res.B.edges.length} edges</span>`;
  }
  function drawNetwork(res) {
    const plot = $("#plot"); $("#zoom-tools").classList.add("hidden");
    const A = res.A, B = res.B; setNetTitle(res);
    if (!A.edges.length && !B.edges.length) {
      plot.innerHTML = "<p class='empty'>No expression-supported interactions for these contexts. Lower the CPM threshold, or search a regulator (e.g. MAX, DDIT3).</p>"; renderNetLegend(new Map()); return;
    }
    if (typeof cytoscape === "undefined") { plot.innerHTML = "<p class='empty'>Cytoscape failed to load.</p>"; renderNetLegend(new Map()); return; }
    const PAL = CLUSTER_COLORS, isoColor = new Map(); let ci = 0;
    [A, B].forEach((s) => s.edges.forEach((e) => { if (!isoColor.has(e.src_iso)) isoColor.set(e.src_iso, PAL[ci++ % PAL.length]); }));
    const nodeMap = new Map();
    [A, B].forEach((s) => s.nodes.forEach((n) => { if (!nodeMap.has(n.id)) nodeMap.set(n.id, { id: n.id, kind: n.kind, gene: n.gene, enst: n.enst, aa: n.aa, nmd: n.nmd }); }));
    const cpmA = new Map(A.nodes.map((n) => [n.id, n.cpm])), cpmB = new Map(B.nodes.map((n) => [n.id, n.cpm]));
    const vmax = Math.max(1, ...[...cpmA.values(), ...cpmB.values(), 1]);
    const _spectral = (typeof d3.interpolateSpectral === "function") ? d3.interpolateSpectral   // matplotlib "Spectral"
      : d3.interpolateRgbBasis(["#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"]);
    const col = (c) => c == null ? "#e9e9e9" : _spectral(1 - Math.log1p(c) / Math.log1p(vmax));   // reversed: low CPM = blue, high = red
    const sz = (c) => c == null ? 8 : 8 + 34 * (Math.log1p(c) / Math.log1p(vmax));   // diameter proportional to log-scaled relative expression (log1p(CPM)/log1p(max))
    const ekey = (e) => e.src_iso + "|" + e.tgt_gene + "|" + e.type;
    const inA = new Set(), inB = new Set();   // only nodes incident to an edge IN that side are shown there
    A.edges.forEach((e) => { inA.add(e.src_iso); inA.add(e.tgt_gene); });
    B.edges.forEach((e) => { inB.add(e.src_iso); inB.add(e.tgt_gene); });
    const inEither = (id) => inA.has(id) || inB.has(id);
    // ---- layout: SERVER igraph Fruchterman-Reingold (visualization/NetPerspective.py); fcose only as fallback ----
    let pos = new Map();
    if (res.layout && Object.keys(res.layout).length) {
      Object.entries(res.layout).forEach(([id, xy]) => { if (inEither(id)) pos.set(id, { x: +xy[0], y: +xy[1] }); });
    }
    if (!pos.size) {
      ensureFcose();
      const layoutEls = [];
      nodeMap.forEach((n) => { if (inEither(n.id)) layoutEls.push({ data: { id: n.id } }); });
      const seen = new Set();
      [A, B].forEach((s) => s.edges.forEach((e) => { const k = ekey(e); if (!seen.has(k)) { seen.add(k); layoutEls.push({ data: { id: k, source: e.src_iso, target: e.tgt_gene } }); } }));
      const byGene = {};
      nodeMap.forEach((n) => { if (n.kind === "isoform" && inEither(n.id)) (byGene[n.gene] = byGene[n.gene] || []).push(n.id); });
      Object.values(byGene).forEach((arr) => { for (let i = 1; i < arr.length; i++) layoutEls.push({ data: { id: "glue|" + i + "|" + arr[i], source: arr[i - 1], target: arr[i], glue: 1 } }); });
      const lay = cytoscape({ headless: true, elements: layoutEls });
      const lname = window.cytoscapeFcose ? "fcose" : "cose";
      const Nn = layoutEls.reduce((a, e) => a + (e.data.source ? 0 : 1), 0); const big = Nn > 120;
      try {
        lay.layout({ name: lname, quality: big ? "default" : "proof", animate: false, randomize: true,
          nodeRepulsion: big ? 30000 : 10000, idealEdgeLength: (e) => e.data("glue") ? 24 : (big ? 230 : 110),
          edgeElasticity: (e) => e.data("glue") ? 0.9 : (big ? 0.08 : 0.3), gravity: big ? 0.03 : 0.18, gravityRange: 5.5,
          numIter: big ? 4500 : 1800, nodeSeparation: big ? 340 : 140, packComponents: true, tile: true, nodeDimensionsIncludeLabels: false }).run();
      } catch (e) { try { lay.layout({ name: "cose", animate: false, nodeRepulsion: big ? 700000 : 200000, idealEdgeLength: big ? 120 : 70 }).run(); } catch (e2) {} }
      lay.nodes().forEach((n) => pos.set(n.id(), { x: n.position("x"), y: n.position("y") })); lay.destroy();
    }
    // normalize to a comfortable model-space density so node separation is consistent (FR coords are unit-scale)
    if (pos.size) {
      const xs = [...pos.values()].map((p) => p.x), ys = [...pos.values()].map((p) => p.y);
      const cur = Math.max(1, Math.hypot(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys)));
      const k = (150 * Math.sqrt(pos.size)) / cur;
      pos.forEach((p) => { p.x *= k; p.y *= k; });
    }
    // ---- two rendered Cytoscape panels: shared positions, each only its edge-incident nodes ----
    const H = Math.max(560, Math.min(880, 140 + Math.ceil(Math.sqrt(nodeMap.size)) * 46));
    plot.innerHTML = `<div style="display:flex;flex-direction:column;height:${H + 26}px">`
      + `<div style="display:flex;font-weight:700;font-size:12px;color:#1f2a3a;padding:4px 0">`
      + `<div style="flex:1;text-align:center">${A.ctx.replace("||", " · ")} <span style="color:#8a94a4;font-weight:400">· ${A.edges.length} edges</span></div>`
      + `<div style="flex:1;text-align:center">${B.ctx.replace("||", " · ")} <span style="color:#8a94a4;font-weight:400">· ${B.edges.length} edges</span></div></div>`
      + `<div style="display:flex;flex:1;min-height:0"><div id="net-a" style="flex:1;border-right:1px dashed #dfe4ec"></div><div id="net-b" style="flex:1"></div></div></div>`;
    const style = [
      { selector: "node", style: { "label": "data(id)", "font-size": 8, "text-valign": "top", "text-margin-y": -1, "color": "#333", "text-background-color": "#ffffff", "text-background-opacity": 0.55, "text-background-padding": 1, "border-width": 0.6, "border-color": "#555", "min-zoomed-font-size": 6 } },
      { selector: 'node[kind="isoform"]', style: { "shape": "diamond", "border-color": "#000", "border-width": 1 } },
      { selector: "edge", style: { "width": 1.5, "line-color": "data(ec)", "opacity": 0.7, "curve-style": "straight" } },
      { selector: 'edge[dir="PDI"]', style: { "target-arrow-shape": "triangle", "target-arrow-color": "data(ec)", "arrow-scale": 0.8 } },
    ];
    // shared bounding box over the WHOLE combined layout (both panels fit to THIS, not their own subset)
    const allXY = [...pos.values()];
    const bx0 = Math.min(...allXY.map((p) => p.x)), bx1 = Math.max(...allXY.map((p) => p.x));
    const by0 = Math.min(...allXY.map((p) => p.y)), by1 = Math.max(...allXY.map((p) => p.y));
    const bcx = (bx0 + bx1) / 2, bcy = (by0 + by1) / 2, bw = Math.max(1, bx1 - bx0), bh = Math.max(1, by1 - by0);
    // tooltip: constituent isoform id(s) + amino-acid length (or NMD); partner genes show gene + CPM
    function netTipHTML(id, ctxLabel, c, isos) {
      const m = nodeMap.get(id) || {};
      const cpmRow = `CPM: <b>${c != null ? c : "n/d"}</b>`;
      const isNmd = (s) => s === "NMD" || s === "Potential-NMD";   // NOT "Not-NMD" (which contains 'NMD')
      if (m.kind === "isoform") {
        const prot = isNmd(m.nmd) ? `<span style="color:#c0392b">${m.nmd}</span>`
          : (m.aa != null ? `${m.aa} aa` : "length n/d");
        return `<b>${id}</b><br><span class="muted">isoform · ${m.gene}</span><br>`
          + (m.enst ? `transcript: <b>${m.enst}</b><br>` : "")
          + `protein: <b>${prot}</b><br>${ctxLabel}<br>${cpmRow}`;
      }
      // partner / target gene: list the isoforms underlying its expression (ranked by CPM in THIS panel)
      let isoRows = "";
      if (isos && isos.length) {
        isoRows = `<br><span class="muted">isoforms (CPM):</span><br>` + isos.map((x) =>
          `&nbsp;${x.enst}${x.aa != null ? " · " + x.aa + " aa" : ""} · <b>${x.cpm}</b>`).join("<br>");
      }
      return `<b>${id}</b><br><span class="muted">partner gene</span><br>${ctxLabel}<br>${cpmRow}${isoRows}`;
    }
    function panel(sel, side, present, cpm) {
      const els = [];
      const sideIsos = new Map(side.nodes.map((n) => [n.id, n.isos]));   // per-panel gene-isoform breakdown
      nodeMap.forEach((n) => { if (!present.has(n.id) || !pos.get(n.id)) return; els.push({ data: { id: n.id, kind: n.kind, gene: n.gene }, position: { ...pos.get(n.id) } }); });
      side.edges.forEach((e) => els.push({ data: { id: ekey(e), source: e.src_iso, target: e.tgt_gene, ec: isoColor.get(e.src_iso), dir: e.type } }));
      const cy = cytoscape({ container: $(sel), elements: els, style, layout: { name: "preset", fit: false }, wheelSensitivity: 0.25, minZoom: 0.05, maxZoom: 5, boxSelectionEnabled: false, autoungrabify: true });   // nodes fixed (no dragging); pan/zoom + click still work
      cy.nodes().forEach((n) => { const c = cpm.get(n.id()); n.style({ "background-color": col(c), "width": sz(c), "height": sz(c) }); });   // size = CPM, MODEL units (constant across panels)
      cy.on("mouseover", "node", (ev) => { const oe = ev.originalEvent || {}; const id = ev.target.id(); const c = cpm.get(id); showTip({ pageX: oe.pageX || 0, pageY: oe.pageY || 0 }, netTipHTML(id, side.ctx.replace("||", " · "), c, sideIsos.get(id))); });
      cy.on("mouseout", "node", hideTip);
      return cy;
    }
    const cyA = panel("#net-a", A, inA, cpmA), cyB = panel("#net-b", B, inB, cpmB);
    // IDENTICAL viewport: fit the SHARED combined bbox -> same zoom + pan -> shared nodes at the same screen
    // position AND CPM-scaled sizes at the same apparent scale on both sides.
    function fitShared(cy) {
      const pad = 42, cw = cy.width() || 500, ch = cy.height() || 500;
      let z = Math.min((cw - 2 * pad) / bw, (ch - 2 * pad) / bh);
      z = Math.max(0.5, Math.min(z, 1.5));   // never shrink below 0.5 -> nodes stay separated (pan/zoom to explore the rest)
      cy.zoom(z); cy.pan({ x: cw / 2 - bcx * z, y: ch / 2 - bcy * z });
    }
    fitShared(cyA); fitShared(cyB);
    let sync = false;
    const linkv = (s, d) => s.on("viewport", () => { if (sync) return; sync = true; d.viewport({ zoom: s.zoom(), pan: { ...s.pan() } }); sync = false; });
    linkv(cyA, cyB); linkv(cyB, cyA);
    // click a node -> its first-degree interactions at full opacity, everything else dimmed to 0.2 (both panels,
    // synced by node id); click empty space to reset.
    function netHighlight(id) {
      [cyA, cyB].forEach((cy) => {
        cy.elements().style("opacity", 0.2);
        const n = cy.getElementById(id);
        if (!n.empty()) n.closedNeighborhood().style("opacity", 1);   // node + its edges + neighbor nodes
      });
    }
    function netClear() { [cyA, cyB].forEach((cy) => cy.elements().style("opacity", 1)); }
    [cyA, cyB].forEach((cy) => {
      cy.on("tap", "node", (ev) => netHighlight(ev.target.id()));
      cy.on("tap", (ev) => { if (ev.target === cy) netClear(); });
    });
    renderNetLegend(isoColor);
  }
  const SPECTRAL_CSS = "linear-gradient(90deg,#5e4fa2,#3288bd,#66c2a5,#abdda4,#e6f598,#ffffbf,#fee08b,#fdae61,#f46d43,#d53e4f,#9e0142)";   // low (blue) -> high (red)
  function renderNetLegend(isoColor) {
    const el = $("#legend"); if (!el) return;
    const arrow = '<svg width="30" height="10" style="vertical-align:middle"><line x1="1" y1="5" x2="20" y2="5" stroke="#7a7f8a" stroke-width="2"/><path d="M20,1.5 L29,5 L20,8.5 Z" fill="#7a7f8a"/></svg>';
    const line = '<svg width="30" height="10" style="vertical-align:middle"><line x1="1" y1="5" x2="29" y2="5" stroke="#7a7f8a" stroke-width="2"/></svg>';
    const sw = [...isoColor.entries()].slice(0, 10).map(([iso, c]) => `<div class="lg-row"><span class="sw" style="background:${c}"></span>${iso}</div>`).join("");
    el.innerHTML =
      `<div class="lg-row"><span>CPM</span><div style="flex:1;height:10px;margin-left:8px;border-radius:3px;background:${SPECTRAL_CSS}"></div></div>`
      + `<div class="lg-row"><span class="muted">low</span><span class="muted" style="margin-left:auto">high</span></div>`
      + `<div class="lg-row">${arrow}&nbsp;PDI · regulator → DNA</div>`
      + `<div class="lg-row">${line}&nbsp;PPI · physical partner</div>`
      + `<div class="lg-row"><span style="font-size:13px;line-height:1">◆</span>&nbsp;isoform&nbsp;&nbsp;<span style="font-size:13px;line-height:1">●</span>&nbsp;gene&nbsp;·&nbsp;size = CPM</div>`
      + (sw ? `<div class="lg-row muted" style="margin-top:4px">edge color = source isoform</div>${sw}` : "");
  }

  function drawHeatStrip(g, isos, cols, x0, HEADER, showBlocks, bodyH) {
    const hdrBase = HEADER - 6;
    if (showBlocks) {
      const groups = [];
      cols.forEach((c, i) => {
        const gn = c.group || "";
        if (!groups.length || groups[groups.length - 1].group !== gn) groups.push({ group: gn, i0: i, i1: i });
        else groups[groups.length - 1].i1 = i;
      });
      const gcolor = d3.scaleOrdinal().domain(groups.map((g) => g.group)).range(CLUSTER_COLORS);
      groups.forEach((gp, gi) => {
        const xa = x0 + gp.i0 * H.CELL_W, xb = x0 + (gp.i1 + 1) * H.CELL_W;
        g.append("rect").attr("class", "block-band").attr("x", xa).attr("y", HEADER - 4).attr("width", xb - xa).attr("height", bodyH + 6).attr("fill", gcolor(gp.group)).attr("opacity", 0.06);
        g.append("text").attr("class", "block-label").attr("x", (xa + xb) / 2).attr("y", H.BLOCK_LBL - 6).attr("fill", gcolor(gp.group)).text(gp.group);
        if (gi > 0) g.append("line").attr("class", "block-sep").attr("x1", xa).attr("x2", xa).attr("y1", H.BLOCK_LBL + 2).attr("y2", HEADER + bodyH);
      });
    }
    cols.forEach((c, i) => {
      const cx = x0 + i * H.CELL_W + H.CELL_W / 2;
      g.append("text").attr("class", "colhdr").attr("x", cx).attr("y", hdrBase).attr("transform", `rotate(-45 ${cx} ${hdrBase})`).text(truncate(c.label, 18));
    });
    const rownorm = $("#rownorm").checked;
    let gmax = 1;
    if (!rownorm) isos.forEach((iso) => cols.forEach((c) => { gmax = Math.max(gmax, Math.log1p(iso.expression[c.key] || 0)); }));
    isos.forEach((iso, r) => {
      const y = HEADER + r * (H.ROW_H + H.ROW_GAP);
      let med = 0, dev = 1;
      if (rownorm) { const vals = cols.map((c) => iso.expression[c.key] || 0); med = d3.median(vals) || 0; dev = d3.max(vals.map((v) => Math.abs(v - med))) || 1; }
      cols.forEach((c, i) => {
        const v = iso.expression[c.key] || 0;
        let fill = "#eef1f5";
        if (v > 0) { const t = rownorm ? clamp(0.5 + 0.5 * (v - med) / dev, 0, 1) : Math.log1p(v) / gmax; fill = cyanYellow(t); }
        g.append("rect").attr("class", "expr-cell").attr("x", x0 + i * H.CELL_W).attr("y", y + 3).attr("width", H.CELL_W - 2).attr("height", H.ROW_H - 6).attr("rx", 2).attr("fill", fill)
          .on("mousemove", (ev) => showTip(ev, `<b>${c.group ? c.group + " · " : ""}${c.label}</b><br>${iso.isoform_id}<br>count: <b>${v}</b>`))
          .on("mouseleave", hideTip);
      });
    });
  }
  function isoTip(iso) {
    return `<b>${iso.isoform_id}</b> <span class="muted">${iso.known ? "known" : "novel"}</span>`
      + `<br>protein length: <b>${iso.protein_length != null ? iso.protein_length + " aa" : "n/a"}</b>`
      + `<br>NMD: ${iso.nmd_status || "n/a"}`
      + (iso.intron_retention && iso.intron_retention !== "False" ? `<br>intron retention: ${iso.intron_retention}` : "")
      + `<br>total count: ${fmtCount(iso.total_count)} · cluster ${iso.cluster_id}`;
  }

  // ===================================================================== tooltip / context menu
  function showTip(ev, html) {
    const t = $("#tooltip"); t.innerHTML = html; t.classList.remove("hidden");
    const pad = 16, w = t.offsetWidth, h = t.offsetHeight;
    let lx = ev.pageX + pad, ly = ev.pageY + 12;
    if (lx + w > window.innerWidth) lx = ev.pageX - w - pad;
    if (ly + h > window.innerHeight) ly = ev.pageY - h - 12;
    t.style.left = lx + "px"; t.style.top = ly + "px";
  }
  function hideTip() { $("#tooltip").classList.add("hidden"); }
  function showCtx(ev, iso, visibleIds, clusterId, allIsoforms) {
    const m = $("#ctxmenu"); m.innerHTML = "";
    // FASTA link-outs key off the FINAL collapsed isoform, so every read in a group shares one RNA + one
    // protein (molecule view). Falls back to isoform_id where there's no collapse mapping (heatmap view).
    const fid = iso.final_isoform_id || iso.isoform_id;
    const head = document.createElement("div"); head.className = "ctx-head";
    head.innerHTML = `${fid}<span>${iso.protein_length != null ? iso.protein_length + " aa" : ""}</span>`;
    m.appendChild(head);
    const item = (label, fn) => { const d = document.createElement("div"); d.textContent = label; d.onclick = () => { fn(); hideCtx(); }; m.appendChild(d); };
    item("Export protein sequence (FASTA)", () => saveFasta("/api/isoform/" + encodeURIComponent(fid) + "/protein", safeName(fid) + ".fasta"));
    item("Export ORF / CDS sequence (FASTA)", () => saveFasta("/api/isoform/" + encodeURIComponent(fid) + "/orf", safeName(fid) + ".orf.fasta"));
    item("Export mRNA sequence (FASTA)", () => saveFasta("/api/isoform/" + encodeURIComponent(fid) + "/mrna", safeName(fid) + ".mrna.fasta"));
    if (clusterId != null && allIsoforms) item("Export all proteins in cluster", () =>
      saveFasta("/api/proteins", "cluster_proteins.fasta", allIsoforms.filter((i) => i.cluster_id === clusterId).map((i) => i.isoform_id)));
    if (visibleIds) item("Export all visible proteins", () => saveFasta("/api/proteins", "visible_proteins.fasta", [...new Set(visibleIds)]));
    item("Copy isoform id", () => navigator.clipboard && navigator.clipboard.writeText(fid));
    m.classList.remove("hidden");
    m.style.left = Math.min(ev.pageX, window.innerWidth - 260) + "px";
    m.style.top = Math.min(ev.pageY, window.innerHeight - m.offsetHeight - 8) + "px";
  }
  function hideCtx() { $("#ctxmenu").classList.add("hidden"); }
  document.addEventListener("click", hideCtx);
  function dl(url) { window.location = url; }
  async function dlPost(url, body, fname) {
    const r = await fetch(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
    if (!r.ok) { setStatus("nothing to export"); return; }
    const blob = await r.blob(); const u = URL.createObjectURL(blob);
    const a = document.createElement("a"); a.href = u; a.download = fname; a.click(); URL.revokeObjectURL(u);
  }
  function safeName(s) { return String(s).replace(/[^A-Za-z0-9._-]+/g, "_"); }
  // Save a FASTA via the native "Save As" picker (choose folder + filename). Open the picker FIRST while the
  // right-click->menu-click user gesture is still valid; fetching before the picker would consume the gesture
  // and make showSaveFilePicker throw silently (the old `dl`/window.location just downloaded with no dialog).
  // Safari/Firefox (no File System Access API) fall back to a normal download. `postBody` set => POST.
  async function saveFasta(url, suggestedName, postBody) {
    let handle = null;
    if (window.showSaveFilePicker) {
      try {
        handle = await window.showSaveFilePicker({
          suggestedName,
          types: [{ description: "FASTA sequence", accept: { "text/plain": [".fasta", ".fa", ".txt"] } }],
        });
      } catch (e) {
        if (e && e.name === "AbortError") { setStatus("export canceled"); return; }
        handle = null;   // picker unavailable/denied -> download fallback below
      }
    }
    setStatus("fetching sequence…");
    let blob;
    try {
      const opt = postBody !== undefined
        ? { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(postBody) }
        : undefined;
      const r = await fetch(url, opt);
      if (!r.ok) { setStatus("nothing to export (HTTP " + r.status + ")"); return; }
      blob = await r.blob();
      if (!blob.size) { setStatus("empty sequence — nothing to export"); return; }
    } catch (e) { setStatus("export error: " + ((e && e.message) || e)); return; }
    if (handle) {
      try {
        const ws = await handle.createWritable(); await ws.write(blob); await ws.close();
        setStatus("saved → " + handle.name); return;
      } catch (e) { setStatus("save error: " + ((e && e.message) || e)); return; }
    }
    const u = URL.createObjectURL(blob);   // fallback: standard download to the browser's default folder
    const a = document.createElement("a"); a.href = u; a.download = suggestedName; a.click(); URL.revokeObjectURL(u);
    setStatus("downloaded " + suggestedName);
  }

  // ===================================================================== PDF export (right-click)
  // RIGHT-click the splicing graph -> "Export graph as PDF…". The displayed D3/SVG is rendered to a
  // VECTOR pdf (svg2pdf, editable text) in the browser, then a NATIVE filesystem picker
  // (showSaveFilePicker: choose folder + filename) writes it. Safari/Firefox fall back to a download.
  const _plotEl = $("#plot");
  if (_plotEl) _plotEl.addEventListener("contextmenu", (ev) => {
    if (!_plotEl.querySelector("svg")) return;          // nothing rendered yet -> let native menu show
    ev.preventDefault();
    showGraphMenu(ev);
  });
  function showGraphMenu(ev) {
    const m = $("#ctxmenu"); m.innerHTML = "";
    const head = document.createElement("div"); head.className = "ctx-head";
    head.innerHTML = `Splicing graph<span>${state.tab || ""}</span>`;
    m.appendChild(head);
    const d = document.createElement("div"); d.textContent = "Export graph as PDF…";
    d.onclick = () => { exportGraphPdf(); hideCtx(); }; m.appendChild(d);
    m.classList.remove("hidden");
    m.style.left = Math.min(ev.pageX, window.innerWidth - 260) + "px";
    m.style.top = Math.min(ev.pageY, window.innerHeight - m.offsetHeight - 8) + "px";
  }
  // svg2pdf does NOT resolve our external stylesheet, so elements colored via CSS classes (panel-box,
  // mol-exon, labels, …) export with the SVG default fill = BLACK (the "black background"). Build an export
  // clone with computed styles inlined and a white background rect painted first, so the PDF matches screen.
  function svgForExport(svg) {
    const clone = svg.cloneNode(true);
    const PROPS = ["fill","fill-opacity","stroke","stroke-width","stroke-opacity","stroke-dasharray","opacity",
                   "font-family","font-size","font-weight","font-style","text-anchor","dominant-baseline","letter-spacing"];
    const src = [svg, ...svg.querySelectorAll("*")], dst = [clone, ...clone.querySelectorAll("*")];
    for (let i = 0; i < src.length; i++) {
      const cs = getComputedStyle(src[i]); let st = dst[i].getAttribute("style") || "";
      for (const p of PROPS) { const v = cs.getPropertyValue(p); if (v) st += p + ":" + v + ";"; }
      dst[i].setAttribute("style", st);
    }
    const NS = "http://www.w3.org/2000/svg", bg = document.createElementNS(NS, "rect");
    bg.setAttribute("x", 0); bg.setAttribute("y", 0);
    bg.setAttribute("width", svg.getAttribute("width") || svg.clientWidth || 1000);
    bg.setAttribute("height", svg.getAttribute("height") || svg.clientHeight || 700);
    bg.setAttribute("fill", "#ffffff");
    clone.insertBefore(bg, clone.firstChild);
    return clone;
  }
  async function exportGraphPdf() {
    const svg = $("#plot").querySelector("svg");
    if (!svg) { setStatus("no graph to export"); return; }
    if (!window.jspdf || !window.jspdf.jsPDF) { setStatus("PDF library not loaded"); return; }
    const gene = ((geneInput && geneInput.value) || "graph").trim() || "graph";
    const fname = (`${gene}_${state.tab || "graph"}.pdf`).replace(/[^A-Za-z0-9._-]+/g, "_");
    // Open the native picker FIRST, while the right-click user gesture is still valid. Rendering the
    // (async) PDF before this would consume the gesture and make showSaveFilePicker throw silently.
    let handle = null;
    if (window.showSaveFilePicker) {
      try {
        handle = await window.showSaveFilePicker({
          suggestedName: fname,
          types: [{ description: "PDF document", accept: { "application/pdf": [".pdf"] } }],
        });
      } catch (e) {
        if (e && e.name === "AbortError") { setStatus("PDF export canceled"); return; }
        handle = null;   // picker unavailable/failed -> download fallback below
      }
    }
    setStatus("rendering PDF…");
    let blob;
    try {
      const w = Math.ceil(+svg.getAttribute("width") || svg.clientWidth || 1000);
      const h = Math.ceil(+svg.getAttribute("height") || svg.clientHeight || 700);
      const doc = new window.jspdf.jsPDF({ unit: "pt", format: [w, h], compress: true });
      const exp = svgForExport(svg);   // inline computed styles + white bg (svg2pdf ignores external CSS -> black)
      if (typeof doc.svg === "function") await doc.svg(exp, { x: 0, y: 0, width: w, height: h });
      else if (window.svg2pdf) await window.svg2pdf(exp, doc, { x: 0, y: 0, width: w, height: h });
      else throw new Error("svg2pdf not loaded");
      blob = doc.output("blob");
    } catch (e) { setStatus("PDF render error: " + ((e && e.message) || e)); return; }
    if (handle) {
      try {
        const ws = await handle.createWritable(); await ws.write(blob); await ws.close();
        setStatus(`saved PDF → ${handle.name}`); return;
      } catch (e) { setStatus("PDF save error: " + ((e && e.message) || e)); return; }
    }
    // FALLBACK (Safari/Firefox, or picker unavailable): normal download.
    const u = URL.createObjectURL(blob);
    const a = document.createElement("a"); a.href = u; a.download = fname; a.click();
    URL.revokeObjectURL(u);
    setStatus(`downloaded ${fname}`);
  }

  // ===================================================================== header / legend
  function titleBlock(res) {
    const gm = res.gene_model || {};
    const strand = gm.strand ? ` <span class="muted">(${gm.chrom || ""}${gm.chrom ? " " : ""}${gm.strand})</span>` : "";
    let pills;
    if (res._mode === "reads") {
      pills = `<span class="pill">${(res.panels || []).length} panels</span><span class="pill">${fmtCount(res.total_molecules)} molecules</span><span class="pill">${res.n_clusters} clusters</span><span class="pill">${(res.cell_states || []).join("+")}</span>`;
    } else {
      pills = `<span class="pill">${res.isoforms.length} isoforms</span><span class="pill">${res.cluster_count} clusters</span><span class="pill">${(res.columns || []).length} columns</span>`;
    }
    $("#title-block").innerHTML = `<span class="gene-title">${res.symbol}</span> <span class="muted">${res.gene}</span>${strand} ${pills}`;
  }
  function renderLegend(res) {
    const el = $("#legend"); if (!el) return;
    if (state.tab === "molecule") {
      const n = res && res.n_clusters ? res.n_clusters : 4;
      const sw = Array.from({ length: Math.min(n, 8) }, (_, i) => `<span class="sw" style="background:${clusterColor(i)}"></span>`).join("");
      el.innerHTML = `<div class="lg-row"><span>clusters</span>${sw}</div><div class="lg-row muted">each row = one molecule · row height ∝ read count</div><div class="lg-row muted">left-click → protein / mRNA FASTA · right-click → graph PDF</div>`;
    } else if (state.tab === "lines") {
      el.innerHTML = `<div class="lg-row muted">% isoform usage by cell type · one panel per covariate</div>
        <div class="lg-row muted">top ${LP.TOPN} isoforms by reads · color shared by line + structure</div>
        <div class="lg-row muted">left-click a structure row → protein / mRNA FASTA · right-click graph → PDF</div>`;
    } else if (state.tab === "network") {
      renderNetLegend(new Map());
    } else {
      const stops = [0, .25, .5, .75, 1].map((t) => `${cyanYellow(t)} ${t * 100}%`).join(",");
      el.innerHTML = `<div class="lg-row"><span>expression</span><div class="lg-grad" style="background:linear-gradient(90deg,${stops})"></div></div>
        <div class="lg-row"><span class="muted">low</span><span class="muted" style="margin-left:auto">high</span></div>
        <div class="lg-row"><span class="sw" style="background:${KNOWN}"></span>known<span class="sw" style="background:${NOVEL};margin-left:10px"></span>novel</div>`;
    }
  }

  // ===================================================================== utils
  function truncate(s, n) { s = String(s); return s.length > n ? s.slice(0, n - 1) + "…" : s; }
  function fmt(v) { return v == null ? "?" : Math.round(v).toLocaleString(); }
  function fmtCount(v) { return v == null ? "0" : (Number.isInteger(v) ? v.toLocaleString() : (+v).toLocaleString(undefined, { maximumFractionDigits: 1 })); }
  function clamp(v, a, b) { return Math.max(a, Math.min(b, v)); }
  function setStatus(s) { $("#status").textContent = s; }

  // deep link: ?gene=BID&tab=molecule&cells=HSC-1,HSC-2,MPP-1&groups=young,AML-NPM1  (auto-renders).
  // With no params, DEFAULT to gene BID in the molecule view with HSC-1/HSC-2/MPP-1/MPP-MEP selected.
  const DEFAULT_GENE = "BID";
  const DEFAULT_CELLS = ["HSC-1", "HSC-2", "MPP-1"];
  async function applyDeepLink() {
    const q = new URLSearchParams(location.search);
    const hasGene = !!q.get("gene");
    const gene = q.get("gene") || DEFAULT_GENE;
    const tab = q.get("tab") || (hasGene ? null : "molecule");
    if (tab === "molecule" || tab === "heatmap" || tab === "lines" || tab === "network") {
      state.tab = tab; $$("#tabbar .tab").forEach((x) => x.classList.toggle("on", x.dataset.tab === tab)); syncTabVisibility();
    }
    const panel = q.get("panel");
    if (panel === "cell_type" || panel === "covariate") {
      state.panelBy = panel; $$("#panel-seg .seg-btn").forEach((x) => x.classList.toggle("on", x.dataset.v === panel));
    }
    let cells = (q.get("cells") || "").split(",").map((s) => s.trim()).filter(Boolean);
    if (!hasGene && !cells.length) cells = DEFAULT_CELLS.slice();
    cells.forEach((c) => state.cellTypes.add(c)); if (cells.length) refreshChips("#cell-types", state.cellTypes);
    const groups = (q.get("groups") || "").split(",").map((s) => s.trim()).filter(Boolean);
    groups.forEach((g) => state.groups.add(g)); if (groups.length) refreshChips("#groups", state.groups);
    if (state.tab === "network") {   // gene-free: uses the cell-type/covariate chips above; apply ?npanel/ntype/nthr/ngenes
      await netInit();
      const np = q.get("npanel");
      if (np === "cell_type" || np === "covariate") {
        state.netPanelBy = np; $$("#net-panel-seg .seg-btn").forEach((x) => x.classList.toggle("on", x.dataset.v === np)); updateNetPanelNote();
      }
      const et = $("#net-edgetype"); if (et && q.get("ntype")) et.value = q.get("ntype");
      const th = $("#net-thresh"); if (th && q.get("nthr")) th.value = q.get("nthr");
      const g2 = $("#net-genes"); if (g2 && q.get("ngenes") != null) g2.value = q.get("ngenes");
      render(); return;
    }
    geneInput.value = gene; loadJunctions(gene).finally(() => render());
  }

  syncTabVisibility();
  if (new URLSearchParams(location.search).get("spin") === "1") { $("#loading-txt").textContent = "rendering…"; $("#loading").classList.remove("hidden"); }
  loadCatalog().then(applyDeepLink).catch((e) => setStatus("catalog error: " + e.message));
})();
