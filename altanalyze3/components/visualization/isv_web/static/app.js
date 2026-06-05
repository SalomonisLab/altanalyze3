/* ISV web frontend: selection controls + D3 SVG isoform tracks with mouseover + right-click export. */
(function () {
  "use strict";
  const $ = (s) => document.querySelector(s);
  const api = (p, o) => fetch(p, o).then((r) => { if (!r.ok) throw new Error(r.status + " " + p); return r; });
  const state = { catalog: null, samples: new Set(), groups: new Set(), cellTypes: new Set(),
                  cellTypeOrder: [], junctions: new Set(), gene: null, last: null };

  // ----- catalog / menus -----
  async function loadCatalog() {
    const c = await api("/api/catalog").then((r) => r.json());
    state.catalog = c;
    renderMulti("#samples", c.samples.map((s) => ({ id: s.library, label: `${s.library}`, sub: s.group })), state.samples);
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
      b.onclick = () => { b.classList.toggle("on"); b.classList.contains("on") ? store.add(it.id) : store.delete(it.id); };
      host.appendChild(b);
    });
  }
  // preserve cell-type selection order (matters for the column layout)
  function selectedCellTypes() { return state.cellTypeOrder.filter((t) => state.cellTypes.has(t)); }

  // ----- gene autocomplete -----
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

  // ----- threshold widgets -----
  $("#sim").addEventListener("input", (e) => $("#sim-val").textContent = e.target.value);
  $("#msf").addEventListener("input", (e) => $("#msf-val").textContent = e.target.value);
  $("#celltype-filter").addEventListener("input", (e) => {
    const f = e.target.value.toLowerCase();
    $("#cell-types").querySelectorAll(".chip").forEach((c) => {
      c.style.display = c.textContent.toLowerCase().includes(f) ? "" : "none";
    });
  });
  $("#render-btn").addEventListener("click", render);

  function queryBody() {
    return {
      gene: state.gene || geneInput.value.trim(),
      samples: [...state.samples], groups: [...state.groups], cell_types: selectedCellTypes(),
      combine_by: $("#combine-by").value,
      cluster_similarity_threshold: +$("#sim").value,
      min_split_fraction: +$("#msf").value,
      min_count: +$("#mincount").value, max_isoforms: +$("#maxiso").value,
      cluster_strategy: $("#strategy").value, cluster_mode: $("#mode").value,
      filter_junctions: [...state.junctions], include_introns: $("#introns").checked,
    };
  }

  async function render() {
    const body = queryBody();
    if (!body.gene) { setStatus("enter a gene"); return; }
    if (body.combine_by !== "group" && body.cell_types.length === 0) {
      setStatus("select cell types for this column mode"); return;
    }
    setStatus("querying " + body.gene + "…");
    const t0 = performance.now();
    try {
      const res = await api("/api/isoforms", {
        method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body),
      }).then((r) => r.json());
      state.last = res;
      draw(res);
      const ms = Math.round(performance.now() - t0);
      $("#header-bar").innerHTML = `<span class="gene-title">${res.symbol}</span>
        <span class="muted">${res.gene}</span>
        <span class="pill">${res.isoforms.length} isoforms</span>
        <span class="pill">${res.cluster_count} clusters</span>
        <span class="pill">${res.columns.length} columns</span>`;
      setStatus(`rendered in ${ms} ms`);
    } catch (e) { setStatus("error: " + e.message); }
  }

  // ----- D3 SVG isoform tracks -----
  const ROW_H = 22, ROW_GAP = 7, LEFT = 250, EXPR_W = 26, PAD_L = 14;
  const HDR = 96;                 // rotated column-label band height (prevents header/grid overlap)
  const BLOCK_LBL = 22;           // covariate block-label band above the headers
  const TRACK_PAD = 44;
  const isoColor = d3.scaleOrdinal(d3.schemeTableau10);   // distinct color PER ISOFORM
  const groupColor = d3.scaleOrdinal().range(["#2b6cb0", "#c05621", "#2f855a", "#6b46c1", "#b83280"]);

  // Merge a sorted list of {start,end,type,label} into contiguous exon BLOCKS (same exon base
  // touching/overlapping regions become one box). Introns stay as thin connectors derived from gaps.
  function mergeSegments(segs) {
    const es = segs.filter((s) => s.start != null && s.end != null
                && (s.type || "E").toUpperCase().startsWith("E"))
      .map((s) => ({ a: Math.min(s.start, s.end), b: Math.max(s.start, s.end),
                     label: (s.exon_id || s.label || "").split(".")[0] }))
      .sort((p, q) => p.a - q.a);
    const out = [];
    es.forEach((s) => {
      const last = out[out.length - 1];
      if (last && s.a <= last.b + 2 && s.label === last.label) { last.b = Math.max(last.b, s.b); }
      else if (last && s.a <= last.b) { last.b = Math.max(last.b, s.b); if (s.label && !last.label) last.label = s.label; }
      else out.push({ a: s.a, b: s.b, label: s.label });
    });
    // explicit intron-retention segments (type I) kept separately
    const introns = segs.filter((s) => s.start != null && s.end != null
                && (s.type || "").toUpperCase().startsWith("I"))
      .map((s) => ({ a: Math.min(s.start, s.end), b: Math.max(s.start, s.end) }));
    return { exons: out, introns };
  }

  function draw(res) {
    const plot = $("#plot"); plot.innerHTML = "";
    if (!res.isoforms.length) { plot.innerHTML = "<p class='empty'>No isoforms for this selection.</p>"; return; }
    const cols = res.columns;                 // [{key,label,group}]

    // merge each isoform's raw exon regions into clean contiguous blocks (fixes the fragmented look)
    const isos = res.isoforms.map((iso) => ({ ...iso, merged: mergeSegments(iso.exon_segments) }));

    // genomic extent (from merged blocks + introns)
    let lo = Infinity, hi = -Infinity;
    isos.forEach((iso) => {
      iso.merged.exons.forEach((s) => { lo = Math.min(lo, s.a); hi = Math.max(hi, s.b); });
      iso.merged.introns.forEach((s) => { lo = Math.min(lo, s.a); hi = Math.max(hi, s.b); });
    });
    if (!isFinite(lo)) { lo = 0; hi = 1; }

    const showBlocks = cols.some((c) => c.group);
    const topBlock = showBlocks ? BLOCK_LBL : 0;
    const headerY = topBlock + HDR;            // baseline where the grid + tracks begin
    const exprX0 = LEFT;
    const trackX0 = exprX0 + cols.length * EXPR_W + TRACK_PAD;
    const width = Math.max(960, trackX0 + 780);
    const x = d3.scaleLinear().domain([lo, hi]).range([trackX0, width - 26]);
    const bodyH = isos.length * (ROW_H + ROW_GAP);
    const height = headerY + bodyH + 24;
    const svg = d3.select(plot).append("svg").attr("width", width).attr("height", height);

    // ---- covariate block bands + labels (only in cell_type_x_covariate mode) ----
    if (showBlocks) {
      const groups = [];
      cols.forEach((c, i) => {
        const g = c.group || "";
        if (!groups.length || groups[groups.length - 1].group !== g) groups.push({ group: g, i0: i, i1: i });
        else groups[groups.length - 1].i1 = i;
      });
      groups.forEach((g, gi) => {
        const xa = exprX0 + g.i0 * EXPR_W, xb = exprX0 + (g.i1 + 1) * EXPR_W;
        svg.append("rect").attr("class", "block-band")
          .attr("x", xa).attr("y", headerY - 4).attr("width", xb - xa).attr("height", bodyH + 6)
          .attr("fill", groupColor(g.group)).attr("opacity", 0.06);
        svg.append("text").attr("class", "block-label")
          .attr("x", (xa + xb) / 2).attr("y", BLOCK_LBL - 6).attr("fill", groupColor(g.group)).text(g.group);
        if (gi > 0) svg.append("line").attr("class", "block-sep")
          .attr("x1", xa).attr("x2", xa).attr("y1", BLOCK_LBL + 2).attr("y2", headerY + bodyH);
      });
    }

    // ---- column headers: rotate around the baseline so text rises into the header band, never over the grid ----
    cols.forEach((c, i) => {
      const cx = exprX0 + i * EXPR_W + EXPR_W / 2;
      svg.append("text").attr("class", "colhdr")
        .attr("x", cx).attr("y", headerY - 6)
        .attr("transform", `rotate(-45 ${cx} ${headerY - 6})`)
        .text(truncate(c.label, 18));
    });

    // ---- expression heat scale (per-column log-ish via global max) ----
    const maxExpr = d3.max(isos, (iso) => d3.max(cols, (c) => iso.expression[c.key] || 0)) || 1;
    const heat = d3.scaleSequential(d3.interpolateBlues).domain([0, Math.max(1, maxExpr)]);

    isos.forEach((iso, r) => {
      const y = headerY + r * (ROW_H + ROW_GAP);
      const row = svg.append("g").attr("class", "isoform-row").attr("transform", `translate(0,${y})`);
      const color = isoColor(iso.isoform_id);        // distinct color per isoform

      row.append("text").attr("class", "iso-label " + (iso.known ? "known" : "novel"))
        .attr("x", PAD_L).attr("y", ROW_H * 0.68).text(truncate(iso.isoform_id, 32));

      // expression heat cells
      cols.forEach((c, i) => {
        const v = iso.expression[c.key] || 0;
        row.append("rect").attr("class", "expr-cell")
          .attr("x", exprX0 + i * EXPR_W).attr("y", 3).attr("width", EXPR_W - 2).attr("height", ROW_H - 6)
          .attr("rx", 2).attr("fill", v > 0 ? heat(v) : "#f5f6f8")
          .on("mousemove", (ev) => showTip(ev, `<b>${c.group ? c.group + " · " : ""}${c.label}</b><br>${iso.isoform_id}<br>count: <b>${v}</b>`))
          .on("mouseleave", hideTip);
      });

      // intron backbone: thin line spanning the isoform's own exon extent
      const exs = iso.merged.exons;
      if (exs.length) {
        const x0 = x(Math.min(...exs.map((e) => e.a))), x1 = x(Math.max(...exs.map((e) => e.b)));
        row.append("line").attr("class", "backbone")
          .attr("x1", x0).attr("x2", x1).attr("y1", ROW_H / 2).attr("y2", ROW_H / 2)
          .attr("stroke", color).attr("opacity", 0.45);
      }
      // intron-retention boxes (thin, in exon color, lighter)
      iso.merged.introns.forEach((s) => {
        const xa = x(s.a), xb = x(s.b);
        row.append("rect").attr("class", "seg intron")
          .attr("x", xa).attr("y", ROW_H / 2 - 2.5).attr("width", Math.max(1.5, xb - xa)).attr("height", 5)
          .attr("fill", color).attr("opacity", 0.5)
          .on("mousemove", (ev) => showTip(ev, `<b>intron retention</b><br>${fmt(s.a)}–${fmt(s.b)}<br><span class="muted">${iso.isoform_id}</span>`))
          .on("mouseleave", hideTip);
      });
      // exon blocks (merged), colored per isoform
      exs.forEach((s) => {
        const xa = x(s.a), xb = x(s.b);
        row.append("rect").attr("class", "seg exon")
          .attr("x", xa).attr("y", 3).attr("width", Math.max(2, xb - xa)).attr("height", ROW_H - 6).attr("rx", 2)
          .attr("fill", color)
          .on("mousemove", (ev) => showTip(ev, `<b>${s.label || "exon"}</b><br>${fmt(s.a)}–${fmt(s.b)}<br><span class="muted">${iso.isoform_id}</span>`))
          .on("mouseleave", hideTip);
      });

      row.append("rect").attr("class", "rowhit")
        .attr("x", 0).attr("y", 0).attr("width", width).attr("height", ROW_H).attr("fill", "transparent")
        .on("mousemove", (ev) => { if (ev.target.classList.contains("rowhit")) showTip(ev, isoTip(iso)); })
        .on("mouseleave", hideTip)
        .on("contextmenu", (ev) => { ev.preventDefault(); showCtx(ev, iso, res); });
    });
  }

  function segTip(iso, seg) {
    return `<b>${seg.exon_id || seg.label || "?"}</b> <span class="muted">(${(seg.type||'E').toUpperCase().startsWith('E')?'exon':'intron'})</span>`
      + `<br>${seg.strand || ""} ${fmt(seg.start)}–${fmt(seg.end)}`
      + `<br><span class="muted">${iso.isoform_id}</span>`;
  }
  function isoTip(iso) {
    return `<b>${iso.isoform_id}</b> <span class="muted">${iso.known ? "known" : "novel"}</span>`
      + `<br>protein length: <b>${iso.protein_length != null ? iso.protein_length + " aa" : "n/a"}</b>`
      + `<br>NMD: ${iso.nmd_status || "n/a"}`
      + (iso.intron_retention && iso.intron_retention !== "False" ? `<br>intron retention: ${iso.intron_retention}` : "")
      + `<br>total count: ${iso.total_count} · cluster ${iso.cluster_id}`;
  }

  // ----- tooltip -----
  function showTip(ev, html) {
    const t = $("#tooltip"); t.innerHTML = html; t.classList.remove("hidden");
    const pad = 16, w = t.offsetWidth, h = t.offsetHeight;
    let lx = ev.pageX + pad, ly = ev.pageY + 12;
    if (lx + w > window.innerWidth) lx = ev.pageX - w - pad;
    if (ly + h > window.innerHeight) ly = ev.pageY - h - 12;
    t.style.left = lx + "px"; t.style.top = ly + "px";
  }
  function hideTip() { $("#tooltip").classList.add("hidden"); }

  // ----- right-click context menu (protein + mRNA export) -----
  function showCtx(ev, iso, res) {
    const m = $("#ctxmenu"); m.innerHTML = "";
    const head = document.createElement("div"); head.className = "ctx-head";
    head.innerHTML = `${iso.isoform_id}<span>${iso.protein_length != null ? iso.protein_length + " aa" : ""}</span>`;
    m.appendChild(head);
    const item = (label, fn) => { const d = document.createElement("div"); d.textContent = label; d.onclick = () => { fn(); hideCtx(); }; m.appendChild(d); };
    item("Export protein sequence (FASTA)", () => dl("/api/isoform/" + encodeURIComponent(iso.isoform_id) + "/protein"));
    item("Export mRNA sequence (FASTA)", () => dl("/api/isoform/" + encodeURIComponent(iso.isoform_id) + "/mrna"));
    item("Export all proteins in cluster", () => {
      const ids = res.isoforms.filter((i) => i.cluster_id === iso.cluster_id).map((i) => i.isoform_id);
      dlPost("/api/proteins", ids, "cluster_proteins.fasta");
    });
    item("Copy isoform id", () => navigator.clipboard && navigator.clipboard.writeText(iso.isoform_id));
    m.classList.remove("hidden");
    const mw = 250;
    m.style.left = Math.min(ev.pageX, window.innerWidth - mw) + "px";
    m.style.top = ev.pageY + "px";
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

  // ----- utils -----
  function truncate(s, n) { s = String(s); return s.length > n ? s.slice(0, n - 1) + "…" : s; }
  function fmt(v) { return v == null ? "?" : Number(v).toLocaleString(); }
  function setStatus(s) { $("#status").textContent = s; }

  loadCatalog().catch((e) => setStatus("catalog error: " + e.message));
})();
