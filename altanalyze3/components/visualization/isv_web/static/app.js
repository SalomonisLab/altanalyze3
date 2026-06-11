/* ISV web frontend.
 *   Tab 1 "Cell type × covariate"  -> /api/isoforms : clustered isoforms, covariate-blocked heat strip
 *                                     + isoform structures (IGV-style, top axis).
 *   Tab 2 "Molecule view" (read-level) -> /api/reads : the ENGINE's own per-molecule rows
 *                                     (plot_isoform_structures_by_conditions output), drawn as a
 *                                     per-covariate read pileup with a shared gene model + axis at the
 *                                     bottom -- the interactive analogue of the ISV PDF.
 * Shared: genomic ruler + reference gene-model track, wheel-zoom / drag-pan, mouseover (exon region +
 * coords + protein length), right-click protein / mRNA FASTA export.
 */
(function () {
  "use strict";
  const $ = (s) => document.querySelector(s);
  const $$ = (s) => Array.from(document.querySelectorAll(s));
  const api = (p, o) => fetch(p, o).then((r) => { if (!r.ok) throw new Error(r.status + " " + p); return r; });
  const state = {
    catalog: null, tab: "heatmap", gene: null, combineBy: "cell_type_x_covariate", panelBy: "covariate",
    samples: new Set(), groups: new Set(), cellTypes: new Set(), cellTypeOrder: [], junctions: new Set(),
    last: null,
  };
  let view = null;

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
      b.onclick = () => { b.classList.toggle("on"); b.classList.contains("on") ? store.add(it.id) : store.delete(it.id); };
      host.appendChild(b);
    });
  }
  const selectedCellTypes = () => state.cellTypeOrder.filter((t) => state.cellTypes.has(t));

  // ===================================================================== tabs
  $$("#tabbar .tab").forEach((t) => t.addEventListener("click", () => {
    if (t.classList.contains("on")) return;
    $$("#tabbar .tab").forEach((x) => x.classList.remove("on")); t.classList.add("on");
    state.tab = t.dataset.tab; syncTabVisibility();
    if (state.gene) render();
  }));
  function syncTabVisibility() {
    $$(".tab-controls").forEach((el) => { el.hidden = el.dataset.for !== state.tab; });
    $$("[data-for]").forEach((el) => {
      if (el.classList.contains("tab-controls")) return;
      el.style.display = el.dataset.for === state.tab ? "" : "none";
    });
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
    state.combineBy = b.dataset.v; if (state.gene) render();
  });
  $("#panel-seg").addEventListener("click", (e) => {
    const b = e.target.closest(".seg-btn"); if (!b) return;
    $$("#panel-seg .seg-btn").forEach((x) => x.classList.remove("on")); b.classList.add("on");
    state.panelBy = b.dataset.v; if (state.gene) render();
  });
  $("#sim").addEventListener("input", (e) => $("#sim-val").textContent = e.target.value);
  $("#celltype-filter").addEventListener("input", (e) => {
    const f = e.target.value.toLowerCase();
    $$("#cell-types .chip").forEach((c) => { c.style.display = c.textContent.toLowerCase().includes(f) ? "" : "none"; });
  });
  $("#ct-all").onclick = () => { state.cellTypeOrder.forEach((t) => state.cellTypes.add(t)); refreshChips("#cell-types", state.cellTypes); };
  $("#ct-none").onclick = () => { state.cellTypes.clear(); refreshChips("#cell-types", state.cellTypes); };
  $("#grp-all").onclick = () => { state.catalog.groups.forEach((g) => state.groups.add(g)); refreshChips("#groups", state.groups); };
  function refreshChips(sel, store) { $$(sel + " .chip").forEach((c) => c.classList.toggle("on", store.has(c.dataset.id))); }
  $("#render-btn").addEventListener("click", render);
  $("#zoom-in").onclick = () => view && view.svg.transition().duration(160).call(view.zoom.scaleBy, 1.6);
  $("#zoom-out").onclick = () => view && view.svg.transition().duration(160).call(view.zoom.scaleBy, 1 / 1.6);
  $("#zoom-reset").onclick = () => view && view.svg.transition().duration(220).call(view.zoom.transform, d3.zoomIdentity);
  window.addEventListener("resize", () => { if (state.last) draw(state.last); });

  // ===================================================================== query dispatch
  function render() {
    const gene = state.gene || geneInput.value.trim();
    if (!gene) { setStatus("enter a gene"); return; }
    return state.tab === "molecule" ? renderReads(gene) : renderHeatmap(gene);
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
      g.append("rect").attr("class", "ref-exon").attr("x", xa).attr("y", refTopY + 3).attr("width", Math.max(1.6, xb - xa)).attr("height", refH - 6).attr("rx", 1.5)
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
    titleBlock(res);
    return res._mode === "reads" ? drawReads(res) : drawHeatmap(res);
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
      const rawSum = raw.reduce((a, b) => a + b, 0);
      const scale = rawSum > R.MAXPANEL ? R.MAXPANEL / rawSum : 1;
      let h = 0; p._rows = [];
      p.mol.forEach((m, i) => { const hh = Math.max(0.55, raw[i] * scale); p._rows.push({ m, dy: h, h: hh }); h += hh; });
      p._bodyH = Math.max(18, h);
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
    const gx = svg.append("g").attr("clip-path", "url(#rclip)");

    // static per-panel chrome: box, top-left title, sparse molecule labels, mouse overlay
    panels.forEach((p) => {
      const top = p._top, boxH = R.HDR + p._bodyH + R.PAD;
      gStatic.append("rect").attr("class", "panel-box").attr("x", 6).attr("y", top).attr("width", width - 12).attr("height", boxH).attr("rx", 7);
      gStatic.append("rect").attr("class", "panel-hdr-band").attr("x", 6).attr("y", top).attr("width", width - 12).attr("height", R.HDR).attr("rx", 7);
      gStatic.append("text").attr("class", "panel-label").attr("x", 16).attr("y", top + 15)
        .html(`${p.condition}<tspan class="panel-sub-inline"> · ${p.n_molecules.toLocaleString()} reads</tspan>`);
      // sparse left labels: first row of each cluster, or rows tall enough
      let seenCl = new Set();
      p._rows.forEach((r) => {
        const first = !seenCl.has(r.m.cluster_index); seenCl.add(r.m.cluster_index);
        if (r.h >= 7 || first) {
          gStatic.append("text").attr("class", "mol-label " + (r.m.known ? "known" : "novel"))
            .attr("x", R.LEFT - 8).attr("y", p._y + r.dy + Math.min(r.h, 9))
            .text(truncate(r.m.isoform_id, 20));
        }
      });
      // one mouse overlay per panel -> resolve molecule by y (avoids thousands of hit rects)
      gStatic.append("rect").attr("class", "panel-hit").attr("x", trackX0).attr("y", p._y).attr("width", trackX1 - trackX0).attr("height", p._bodyH).attr("fill", "transparent")
        .on("mousemove", function (ev) {
          const my = d3.pointer(ev, this)[1];
          const r = p._rows.find((rr) => my >= rr.dy && my < rr.dy + rr.h);
          if (r) showTip(ev, molTip(r.m, p.condition));
        })
        .on("mouseleave", hideTip)
        .on("contextmenu", function (ev) {
          const my = d3.pointer(ev, this)[1];
          const r = p._rows.find((rr) => my >= rr.dy && my < rr.dy + rr.h);
          if (r) { ev.preventDefault(); showCtx(ev, r.m, p.mol.map((m) => m.isoform_id)); }
        });
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
  function molTip(m, cond) {
    return `<b>${m.isoform_id}</b> <span class="muted">${m.known ? "known" : "novel"}</span>`
      + (m.sample ? `<br>sample: ${m.sample}` : "") + `<br>${cond} · cluster ${m.cluster_index}`
      + `<br>reads: <b>${fmtCount(m.count)}</b>`
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

    const showBlocks = cols.some((c) => c.group);
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
      row.append("rect").attr("class", "rowhit").attr("x", 0).attr("y", 0).attr("width", width).attr("height", H.ROW_H).attr("fill", "transparent")
        .on("mousemove", (ev) => { if (ev.target.classList.contains("rowhit")) showTip(ev, isoTip(iso)); })
        .on("mouseleave", hideTip)
        .on("contextmenu", (ev) => { ev.preventDefault(); showCtx(ev, iso, isos.map((i) => i.isoform_id), iso.cluster_id, res.isoforms); });
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
        iso.merged.introns.forEach((s) => {
          const xa = xz(s.a), xb = xz(s.b);
          g.append("rect").attr("class", "seg intron").attr("x", xa).attr("y", H.ROW_H / 2 - 2.5).attr("width", Math.max(1.5, xb - xa)).attr("height", 5).attr("fill", color).attr("opacity", 0.5)
            .on("mousemove", (ev) => showTip(ev, `<b>intron retention</b><br>${fmt(s.a)}–${fmt(s.b)}<br><span class="muted">${iso.isoform_id}</span>`)).on("mouseleave", hideTip);
        });
        exs.forEach((s) => {
          const xa = xz(s.a), xb = xz(s.b);
          g.append("rect").attr("class", "seg exon").attr("x", xa).attr("y", 3).attr("width", Math.max(2, xb - xa)).attr("height", H.ROW_H - 6).attr("rx", 2).attr("fill", color)
            .on("mousemove", (ev) => showTip(ev, `<b>${s.label || "exon"}</b><br>${gm.chrom ? gm.chrom + ":" : ""}${fmt(s.a)}–${fmt(s.b)}<br>${iso.protein_length != null ? iso.protein_length + " aa · " : ""}<span class="muted">${iso.isoform_id}</span>`))
            .on("mouseleave", hideTip)
            .on("contextmenu", (ev) => { ev.preventDefault(); showCtx(ev, iso, isos.map((i) => i.isoform_id), iso.cluster_id, res.isoforms); });
        });
      });
    }
    installZoom(svg, x0, trackX0, trackX1, height, drawX);
    renderLegend(res);
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
    const head = document.createElement("div"); head.className = "ctx-head";
    head.innerHTML = `${iso.isoform_id}<span>${iso.protein_length != null ? iso.protein_length + " aa" : ""}</span>`;
    m.appendChild(head);
    const item = (label, fn) => { const d = document.createElement("div"); d.textContent = label; d.onclick = () => { fn(); hideCtx(); }; m.appendChild(d); };
    item("Export protein sequence (FASTA)", () => dl("/api/isoform/" + encodeURIComponent(iso.isoform_id) + "/protein"));
    item("Export ORF / CDS sequence (FASTA)", () => dl("/api/isoform/" + encodeURIComponent(iso.isoform_id) + "/orf"));
    item("Export mRNA sequence (FASTA)", () => dl("/api/isoform/" + encodeURIComponent(iso.isoform_id) + "/mrna"));
    if (clusterId != null && allIsoforms) item("Export all proteins in cluster", () =>
      dlPost("/api/proteins", allIsoforms.filter((i) => i.cluster_id === clusterId).map((i) => i.isoform_id), "cluster_proteins.fasta"));
    if (visibleIds) item("Export all visible proteins", () => dlPost("/api/proteins", visibleIds, "visible_proteins.fasta"));
    item("Copy isoform id", () => navigator.clipboard && navigator.clipboard.writeText(iso.isoform_id));
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
      el.innerHTML = `<div class="lg-row"><span>clusters</span>${sw}</div><div class="lg-row muted">each row = one molecule · row height ∝ read count</div><div class="lg-row muted">right-click → export protein / mRNA</div>`;
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

  // deep link: ?gene=HOPX&tab=molecule&cells=HSC-1,HSC-2,MPP-1&groups=young,AML-NPM1  (auto-renders)
  function applyDeepLink() {
    const q = new URLSearchParams(location.search);
    const gene = q.get("gene"); if (!gene) return;
    const tab = q.get("tab");
    if (tab === "molecule" || tab === "heatmap") {
      state.tab = tab; $$("#tabbar .tab").forEach((x) => x.classList.toggle("on", x.dataset.tab === tab)); syncTabVisibility();
    }
    const panel = q.get("panel");
    if (panel === "cell_type" || panel === "covariate") {
      state.panelBy = panel; $$("#panel-seg .seg-btn").forEach((x) => x.classList.toggle("on", x.dataset.v === panel));
    }
    const cells = (q.get("cells") || "").split(",").map((s) => s.trim()).filter(Boolean);
    cells.forEach((c) => state.cellTypes.add(c)); if (cells.length) refreshChips("#cell-types", state.cellTypes);
    const groups = (q.get("groups") || "").split(",").map((s) => s.trim()).filter(Boolean);
    groups.forEach((g) => state.groups.add(g)); if (groups.length) refreshChips("#groups", state.groups);
    geneInput.value = gene; loadJunctions(gene).finally(() => render());
  }

  syncTabVisibility();
  if (new URLSearchParams(location.search).get("spin") === "1") { $("#loading-txt").textContent = "rendering…"; $("#loading").classList.remove("hidden"); }
  loadCatalog().then(applyDeepLink).catch((e) => setStatus("catalog error: " + e.message));
})();
