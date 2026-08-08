/* scALABLE viewer - vanilla JS client.
 *
 * Everything the client draws comes from a precomputed bundle. Cell states always
 * keep the canonical order the server sends; the client never sorts them.
 */
"use strict";

const S = {
  ds: null, meta: null,
  emb: null, embN: 0, bounds: null,
  color: null,        // {key, kind, categories, colors, values:TypedArray}
  overlay: null,      // {gene, values: Float32Array(N), max}
  ptSize: 2,
};

const $ = (id) => document.getElementById(id);
const status = (t) => { $("status").textContent = t || ""; };
const esc = (s) => String(s == null ? "" : s).replace(/[&<>"]/g,
  (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const fmt = (v, d = 3) => (v === null || v === undefined || Number.isNaN(v)) ? "" : Number(v).toPrecision(d);

async function jget(url) {
  const r = await fetch(url, { cache: "no-store" });
  if (!r.ok) throw new Error(`${url} -> HTTP ${r.status}: ${(await r.text()).slice(0, 200)}`);
  return r.json();
}
async function bget(url) {
  const r = await fetch(url, { cache: "no-store" });
  if (!r.ok) throw new Error(`${url} -> HTTP ${r.status}`);
  return { buf: await r.arrayBuffer(), headers: r.headers };
}

/* ------------------------------------------------------------------ tabs */
function showTab(name) {
  const b = document.querySelector(`.tab[data-tab="${name}"]`);
  if (!b) return;
  document.querySelectorAll(".tab").forEach((x) => x.classList.remove("active"));
  document.querySelectorAll(".panel").forEach((x) => x.classList.remove("active"));
  b.classList.add("active");
  $("panel-" + name).classList.add("active");
  if (name === "embed") drawUmap();
}
document.querySelectorAll(".tab").forEach((b) =>
  b.addEventListener("click", () => showTab(b.dataset.tab)));

/* Deep links. Every view is reachable by URL, so a link can be shared or bookmarked:
 *   ?ds=<id>&tab=expr&gene=SFTPC
 *   ?tab=expr&genes=SFTPC,SCGB1A1
 *   ?tab=deg&deg=<comparison id>&fdr=0.05
 *   ?tab=heat&top=5
 */
async function applyDeepLink() {
  const q = new URLSearchParams(window.location.search);
  if (!q.toString()) return;
  const tab = q.get("tab");
  if (q.get("gene")) {
    $("geneInput").value = q.get("gene");
    $("violinGene").value = q.get("gene");
    await setOverlay(q.get("gene"));
    if (tab === "expr") $("violinGo").click();
  }
  if (q.get("genes")) { $("dotGenes").value = q.get("genes"); $("dotGo").click(); }
  if (q.get("colorby")) { $("colorby").value = q.get("colorby"); await loadColor(q.get("colorby")); }
  if (q.get("top")) $("heatTop").value = q.get("top");
  if (q.get("fdr")) $("degFdr").value = q.get("fdr");
  if (q.get("state")) $("degState").value = q.get("state");
  if (q.get("deg")) { $("degSelect").value = q.get("deg"); $("degGo").click(); }
  if (tab === "heat") $("heatGo").click();
  if (tab) showTab(tab);
}

/* --------------------------------------------------------------- catalog */
async function init() {
  status("loading catalog ...");
  const cat = await jget("/api/catalog");
  const sel = $("dataset");
  sel.innerHTML = "";
  cat.datasets.forEach((d) => {
    const o = document.createElement("option");
    o.value = d.id;
    o.textContent = `${d.label}  (${(d.n_cells || 0).toLocaleString()} cells)`;
    sel.appendChild(o);
  });
  if (cat.load_errors && cat.load_errors.length) {
    status(`${cat.load_errors.length} bundle(s) failed to load - see Provenance`);
  }
  sel.addEventListener("change", () => loadDataset(sel.value));
  const want = new URLSearchParams(window.location.search).get("ds");
  const first = (want && cat.datasets.some((d) => d.id === want)) ? want
    : (cat.datasets.length ? cat.datasets[0].id : null);
  if (!first) { status("no dataset in the catalog"); return; }
  sel.value = first;
  await loadDataset(first);
  await applyDeepLink();
}

async function loadDataset(id) {
  S.ds = id; S.overlay = null;
  status("loading dataset ...");
  const t0 = performance.now();
  S.meta = await jget(`/api/dataset/${id}/meta`);
  $("dsinfo").textContent =
    `${S.meta.n_cells.toLocaleString()} cells x ${S.meta.n_genes.toLocaleString()} genes, ` +
    `${S.meta.states.length} cell states, layer "${S.meta.layer}"`;

  const cb = $("colorby");
  cb.innerHTML = "";
  const optState = document.createElement("option");
  optState.value = "cell_state";
  optState.textContent = `cell state (${S.meta.cluster_key})`;
  cb.appendChild(optState);
  Object.keys(S.meta.covariates).sort().forEach((k) => {
    const o = document.createElement("option");
    o.value = k; o.textContent = `${k} (${S.meta.covariates[k].kind})`;
    cb.appendChild(o);
  });
  cb.onchange = () => loadColor(cb.value);

  const { buf, headers } = await bget(`/api/dataset/${id}/embedding.bin`);
  S.emb = new Float32Array(buf);
  S.embN = S.emb.length / 2;
  if (S.embN !== S.meta.n_cells) {
    status(`WARNING: embedding has ${S.embN} points but metadata says ${S.meta.n_cells} cells`);
  }
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (let i = 0; i < S.embN; i++) {
    const x = S.emb[2 * i], y = S.emb[2 * i + 1];
    if (x < x0) x0 = x; if (x > x1) x1 = x;
    if (y < y0) y0 = y; if (y > y1) y1 = y;
  }
  S.bounds = { x0, x1, y0, y1 };
  await loadColor("cell_state");
  await loadDegList();
  await loadCcc();
  renderAbout();
  $("heatTsv").href = `/api/dataset/${id}/heatmap.tsv?top=5&row_scale=true`;
  status(`ready in ${Math.round(performance.now() - t0)} ms`);
}

async function loadColor(key) {
  const legend = await jget(`/api/dataset/${S.ds}/colorby?key=${encodeURIComponent(key)}`);
  const { buf } = await bget(`/api/dataset/${S.ds}/colorby.bin?key=${encodeURIComponent(key)}`);
  legend.values = legend.dtype === "int16" ? new Int16Array(buf)
    : legend.dtype === "int32" ? new Int32Array(buf) : new Float32Array(buf);
  S.color = legend;
  drawLegend();
  drawUmap();
}

/* ------------------------------------------------------------------ UMAP */
function hex2rgb(h) {
  const s = (h || "#BBBBBB").replace("#", "");
  return [parseInt(s.slice(0, 2), 16), parseInt(s.slice(2, 4), 16), parseInt(s.slice(4, 6), 16)];
}
// Cyan -> yellow ramp for continuous overlays (project heatmap convention).
function ramp(t) {
  t = Math.max(0, Math.min(1, t));
  const a = [0x00, 0xFF, 0xFF], b = [0xFF, 0xFF, 0x00];
  return [Math.round(a[0] + (b[0] - a[0]) * t), Math.round(a[1] + (b[1] - a[1]) * t),
          Math.round(a[2] + (b[2] - a[2]) * t)];
}

function drawUmap() {
  const cv = $("umap");
  if (!S.emb || !cv.clientWidth) return;
  const dpr = window.devicePixelRatio || 1;
  const W = cv.clientWidth, H = cv.clientHeight;
  cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
  const g = cv.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, W, H);

  const pad = 18, b = S.bounds;
  const sx = (W - 2 * pad) / Math.max(b.x1 - b.x0, 1e-9);
  const sy = (H - 2 * pad) / Math.max(b.y1 - b.y0, 1e-9);
  const s = Math.min(sx, sy);
  const ox = pad + (W - 2 * pad - (b.x1 - b.x0) * s) / 2;
  const oy = pad + (H - 2 * pad - (b.y1 - b.y0) * s) / 2;
  const r = S.ptSize;

  const t0 = performance.now();
  // Bucket points by colour so fillStyle changes once per colour, not once per point.
  const buckets = new Map();
  const push = (k, i) => { let a = buckets.get(k); if (!a) { a = []; buckets.set(k, a); } a.push(i); };

  if (S.overlay) {
    const v = S.overlay.values, mx = S.overlay.max || 1;
    for (let i = 0; i < S.embN; i++) push(Math.min(31, Math.floor((v[i] / mx) * 31.999)), i);
  } else if (S.color.kind === "categorical") {
    const c = S.color.values;
    for (let i = 0; i < S.embN; i++) push(c[i], i);
  } else {
    const v = S.color.values, lo = S.color.min, hi = S.color.max, span = (hi - lo) || 1;
    for (let i = 0; i < S.embN; i++) {
      push(Number.isFinite(v[i]) ? Math.min(31, Math.floor(((v[i] - lo) / span) * 31.999)) : -1, i);
    }
  }
  // Zero / lowest bucket first so the signal draws on top.
  const keys = [...buckets.keys()].sort((a, z) => a - z);
  for (const k of keys) {
    let col;
    if (S.overlay) { const c = ramp(k / 31); col = `rgb(${c[0]},${c[1]},${c[2]})`; }
    else if (S.color.kind === "categorical") col = (S.color.colors && S.color.colors[k]) || "#BBBBBB";
    else if (k < 0) col = "#E0E0E0";
    else { const c = ramp(k / 31); col = `rgb(${c[0]},${c[1]},${c[2]})`; }
    g.fillStyle = col;
    for (const i of buckets.get(k)) {
      g.fillRect(ox + (S.emb[2 * i] - b.x0) * s - r / 2, oy + (b.y1 - S.emb[2 * i + 1]) * s - r / 2, r, r);
    }
  }
  const ms = Math.round(performance.now() - t0);
  $("umapNote").innerHTML =
    `${S.embN.toLocaleString()} points drawn in ${ms} ms &middot; ` +
    (S.overlay ? `overlay <b>${esc(S.overlay.gene)}</b> (max ${fmt(S.overlay.max)}, layer ${esc(S.meta.layer)})`
               : `coloured by <b>${esc(S.color.key)}</b>`) +
    ` &middot; ${esc(S.meta.embedding_method || "")}`;
}

function drawLegend() {
  const el = $("legend");
  if (!S.color) { el.innerHTML = ""; return; }
  if (S.color.kind === "categorical") {
    const rows = S.color.categories.map((c, i) =>
      `<div class="row"><span class="sw" style="background:${esc(S.color.colors[i])}"></span>` +
      `<span class="lbl" title="${esc(c)}">${esc(c)}</span>` +
      `<span class="cnt">${(S.color.counts ? S.color.counts[i] : 0).toLocaleString()}</span></div>`).join("");
    el.innerHTML = `<b>${esc(S.color.key)}</b> &mdash; ${S.color.categories.length} categories` +
      (S.color.order_source ? `<div class="note">order: ${esc(S.color.order_source)}</div>` : "") +
      rows;
  } else {
    el.innerHTML = `<b>${esc(S.color.key)}</b><div class="ramp" ` +
      `style="background:linear-gradient(90deg,#00FFFF,#FFFF00)"></div>` +
      `<div class="row"><span>${fmt(S.color.min)}</span><span class="cnt">${fmt(S.color.max)}</span></div>` +
      `<div class="note">${S.color.n_missing.toLocaleString()} of ${S.embN.toLocaleString()} cells have no value (grey)</div>`;
  }
}

/* --------------------------------------------------------- gene overlay */
async function setOverlay(gene) {
  if (!gene) return;
  status(`loading ${gene} ...`);
  const t0 = performance.now();
  let res;
  try { res = await bget(`/api/dataset/${S.ds}/expression.bin?gene=${encodeURIComponent(gene)}`); }
  catch (e) { status(`gene not found: ${gene}`); return; }
  const dv = new DataView(res.buf);
  const n = dv.getInt32(0, true);
  const idx = new Uint32Array(res.buf, 4, n);
  const val = new Float32Array(res.buf, 4 + 4 * n, n);
  const full = new Float32Array(S.embN);          // zeros are real values, not missing
  let mx = 0;
  for (let i = 0; i < n; i++) { full[idx[i]] = val[i]; if (val[i] > mx) mx = val[i]; }
  S.overlay = { gene: res.headers.get("X-Gene") || gene, values: full, max: mx };
  drawUmap();
  status(`${gene}: ${n.toLocaleString()} of ${S.embN.toLocaleString()} cells detected, ` +
         `${Math.round(performance.now() - t0)} ms`);
}

$("geneInput").addEventListener("input", async (e) => {
  const q = e.target.value.trim();
  const box = $("geneSuggest");
  if (q.length < 2) { box.innerHTML = ""; return; }
  const r = await jget(`/api/dataset/${S.ds}/genes?q=${encodeURIComponent(q)}&limit=12`);
  box.innerHTML = r.matches.map((m) => `<div data-g="${esc(m)}">${esc(m)}</div>`).join("");
  box.querySelectorAll("div").forEach((d) => d.addEventListener("click", () => {
    $("geneInput").value = d.dataset.g; box.innerHTML = ""; setOverlay(d.dataset.g);
  }));
});
$("geneInput").addEventListener("keydown", (e) => {
  if (e.key === "Enter") { $("geneSuggest").innerHTML = ""; setOverlay(e.target.value.trim()); }
});
$("clearOverlay").addEventListener("click", () => { S.overlay = null; drawUmap(); status(""); });
$("ptSize").addEventListener("input", (e) => { S.ptSize = parseFloat(e.target.value); drawUmap(); });
window.addEventListener("resize", () => drawUmap());

/* -------------------------------------------------------- violin + dots */
$("violinGo").addEventListener("click", async () => {
  const gene = $("violinGene").value.trim();
  if (!gene) return;
  status("computing violins ...");
  let v;
  try { v = await jget(`/api/dataset/${S.ds}/violin?gene=${encodeURIComponent(gene)}`); }
  catch (e) { $("exprInfo").innerHTML = `<span class="warn">${esc(e.message)}</span>`; status(""); return; }
  const st = v.states, n = st.length;
  const W = Math.max(760, n * 26), H = 470, padL = 52, padB = 170, padT = 16;
  const vmax = Math.max(v.value_max, 1e-9);
  const bw = (W - padL - 10) / n;
  const yOf = (val) => padT + (H - padT - padB) * (1 - val / vmax);
  let out = `<h3>${esc(v.gene)} per cell state &mdash; layer ${esc(v.layer)}</h3><svg width="${W}" height="${H}">`;
  // axis: one 2-point path per line
  out += `<path d="M${padL} ${padT} L${padL} ${H - padB}" stroke="#1A1A1A" fill="none"/>`;
  out += `<path d="M${padL} ${H - padB} L${W - 8} ${H - padB}" stroke="#1A1A1A" fill="none"/>`;
  for (let t = 0; t <= 4; t++) {
    const val = (vmax * t) / 4, y = yOf(val);
    out += `<path d="M${padL - 4} ${y} L${padL} ${y}" stroke="#1A1A1A" fill="none"/>`;
    out += `<text x="${padL - 7}" y="${y + 3}" font-size="9" text-anchor="end" fill="#6B6B6B">${fmt(val, 2)}</text>`;
  }
  st.forEach((s, i) => {
    const cx = padL + bw * (i + 0.5), half = bw * 0.42;
    const col = v.state_colors[i] || "#BBBBBB";
    // Each violin is scaled to its OWN largest bin, so shapes stay comparable across
    // states. A shared scale would let the zero bin of the biggest state flatten
    // every other violin to a line. Absolute counts stay in the table below.
    const mh = Math.max(1, ...(s.hist.length ? s.hist : [0]));
    let up = "", dn = "";
    s.hist.forEach((h, k) => {
      const y = yOf((v.bin_edges[k] + v.bin_edges[k + 1]) / 2);
      const w = half * (h / mh);
      up += `${up ? "L" : "M"}${(cx + w).toFixed(1)} ${y.toFixed(1)}`;
      dn = `L${(cx - w).toFixed(1)} ${y.toFixed(1)}` + dn;
    });
    if (up) out += `<path d="${up}${dn}Z" fill="${esc(col)}" fill-opacity="0.75" stroke="${esc(col)}"/>`;
    out += `<path d="M${cx} ${yOf(s.p05)} L${cx} ${yOf(s.p95)}" stroke="#1A1A1A" fill="none"/>`;
    out += `<circle cx="${cx}" cy="${yOf(s.median)}" r="1.8" fill="#1A1A1A"/>`;
    out += `<text x="${cx}" y="${H - padB + 6}" font-size="9" fill="#1A1A1A" text-anchor="end" ` +
           `transform="rotate(-60 ${cx} ${H - padB + 6})">${esc(s.state)}</text>`;
  });
  out += `</svg>`;
  out += `<table><tr><th>cell state</th><th>n</th><th>detected</th><th>frac</th><th>mean</th>` +
         `<th>median</th><th>p95</th></tr>` +
    st.map((s) => `<tr><td>${esc(s.state)}</td><td class="num">${s.n.toLocaleString()}</td>` +
      `<td class="num">${s.n_detected.toLocaleString()}</td><td class="num">${fmt(s.frac)}</td>` +
      `<td class="num">${fmt(s.mean)}</td><td class="num">${fmt(s.median)}</td>` +
      `<td class="num">${fmt(s.p95)}</td></tr>`).join("") + `</table>`;
  $("violin").innerHTML = out;
  $("exprInfo").innerHTML = `violin computed server-side in ${v.elapsed_ms} ms. ` +
    `Zero values are kept, so n = detected + not detected. Each violin is width-scaled to ` +
    `its own largest bin; compare heights and the table, not widths across states. ` +
    `Bar = p05 to p95, dot = median.`;
  status("");
});

$("dotGo").addEventListener("click", async () => {
  const genes = $("dotGenes").value.trim();
  if (!genes) return;
  status("loading dot plot ...");
  let d;
  try { d = await jget(`/api/dataset/${S.ds}/dotplot?genes=${encodeURIComponent(genes)}`); }
  catch (e) { $("dotplot").innerHTML = `<span class="warn">${esc(e.message)}</span>`; status(""); return; }
  const nS = d.states.length, nG = d.genes.length;
  const cw = 17, rh = 19, padL = 150, padT = 150;
  const W = padL + nS * cw + 20, H = padT + nG * rh + 20;
  const perGene = $("dotScale").checked;
  let mx = 0; d.mean.forEach((r) => r.forEach((v) => { if (v > mx) mx = v; }));
  // Per-gene scaling makes each gene's cell-state pattern readable. It changes only the
  // colour: the tooltip and the units below always report the absolute mean.
  const rowMax = d.mean.map((r) => Math.max(1e-12, ...r));
  let out = `<h3>Dot plot &mdash; ${nG} gene(s) x ${nS} cell states (canonical order)</h3><svg width="${W}" height="${H}">`;
  d.states.forEach((s, j) => {
    const x = padL + cw * (j + 0.5);
    out += `<text x="${x}" y="${padT - 6}" font-size="9" text-anchor="start" fill="#1A1A1A" ` +
           `transform="rotate(-90 ${x} ${padT - 6})">${esc(s)}</text>`;
  });
  d.genes.forEach((g, i) => {
    const y = padT + rh * (i + 0.5);
    out += `<text x="${padL - 6}" y="${y + 3}" font-size="10" text-anchor="end">${esc(g)}</text>`;
    d.states.forEach((s, j) => {
      const m = d.mean[i][j], f = d.frac[i][j];
      const denom = perGene ? rowMax[i] : mx;
      const c = ramp(denom > 0 ? m / denom : 0);
      const r = 1 + 7 * Math.sqrt(Math.max(f, 0));
      out += `<circle cx="${padL + cw * (j + 0.5)}" cy="${y}" r="${r.toFixed(2)}" ` +
             `fill="rgb(${c[0]},${c[1]},${c[2]})" stroke="#8C8C8C" stroke-width="0.4">` +
             `<title>${esc(g)} / ${esc(s)}\nmean ${fmt(m)}  frac ${fmt(f)}  n ${d.state_n[j]}</title></circle>`;
    });
  });
  out += `</svg><div class="note">Colour = mean of layer ${esc(d.layer)}, cyan low to yellow high, ` +
    (perGene ? `each gene scaled to its own maximum (hover for the absolute value).`
             : `one shared scale, maximum ${fmt(mx)}.`) +
    ` Radius = fraction of cells detected. ${d.n_returned} of ${d.n_requested} genes found` +
    (d.n_missing ? `; missing: ${esc(d.missing.join(", "))}` : "") +
    `. Cell states are in canonical order.</div>`;
  $("dotplot").innerHTML = out;
  status("");
});

/* --------------------------------------------------------------- heatmap */
$("heatGo").addEventListener("click", () => {
  const top = parseInt($("heatTop").value, 10) || 5;
  const scale = $("heatScale").checked;
  const scheme = $("heatScheme").value;
  const qs = `top=${top}&row_scale=${scale}&scheme=${scheme}`;
  $("morpheus").src = `/api/dataset/${S.ds}/morpheus?${qs}`;
  $("heatTsv").href = `/api/dataset/${S.ds}/heatmap.tsv?top=${top}&row_scale=${scale}`;
  fetch(`/api/dataset/${S.ds}/heatmap.tsv?top=${top}&row_scale=${scale}`, { method: "GET" })
    .then((r) => { $("heatInfo").textContent =
      `${r.headers.get("X-Rows")} marker genes x ${r.headers.get("X-Cols")} cell states; ` +
      `${r.headers.get("X-Missing")} requested genes absent. Morpheus loads from the Broad CDN.`; });
});

/* ------------------------------------------------------------------- DEG */
async function loadDegList() {
  const r = await jget(`/api/dataset/${S.ds}/deg`);
  const sel = $("degSelect");
  sel.innerHTML = r.comparisons.map((c) =>
    `<option value="${esc(c.id)}">${esc(c.comparison)} - ${esc(c.kind)} (${c.n_rows} rows)</option>`).join("");
  $("degState").innerHTML = `<option value="">(all cell states)</option>` +
    S.meta.states.map((s) => `<option value="${esc(s)}">${esc(s)}</option>`).join("");
  $("degInfo").textContent = r.n_comparisons
    ? `${r.n_comparisons} precomputed table(s) in this bundle.`
    : "no precomputed differential table in this bundle.";
}

$("degGo").addEventListener("click", async () => {
  const id = $("degSelect").value;
  if (!id) return;
  status("loading DEG table ...");
  const fdr = parseFloat($("degFdr").value);
  const st = $("degState").value;
  let q = `/api/dataset/${S.ds}/deg/table?id=${encodeURIComponent(id)}&max_rows=500`;
  if (Number.isFinite(fdr)) q += `&fdr_max=${fdr}`;
  if (st) q += `&state=${encodeURIComponent(st)}`;
  const d = await jget(q);
  drawVolcano(d);
  $("degTable").innerHTML =
    `<table><tr><th>gene</th><th>log2FC</th><th>FDR</th><th>p</th><th>cell state</th></tr>` +
    d.rows.map((r) => `<tr><td>${esc(r.gene)}</td><td class="num">${fmt(r.log2fc)}</td>` +
      `<td class="num">${fmt(r.fdr, 2)}</td><td class="num">${fmt(r.pval, 2)}</td>` +
      `<td>${esc(r.population || "")}</td></tr>`).join("") + `</table>`;
  $("degInfo").innerHTML =
    `<b>${esc(d.comparison)}</b> (${esc(d.kind)})<br>` +
    `${d.n_rows_file.toLocaleString()} rows in the file, ` +
    `${d.n_rows_after_filters.toLocaleString()} after filters, ` +
    `${d.n_rows_returned.toLocaleString()} shown.<br>` +
    `${d.n_rows_no_fdr.toLocaleString()} of ${d.n_rows_file.toLocaleString()} rows carry no FDR.<br>` +
    `source: <code>${esc(d.source)}</code>`;
  status("");
});

function drawVolcano(d) {
  const cv = $("volcano");
  const dpr = window.devicePixelRatio || 1;
  const W = cv.clientWidth || 700, H = 420;
  cv.width = W * dpr; cv.height = H * dpr;
  const g = cv.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, W, H);
  const pts = d.volcano.filter((p) => p.x !== null && p.y !== null && Number.isFinite(p.y));
  const padL = 46, padB = 34, padT = 14, padR = 12;
  if (!pts.length) {
    g.fillStyle = "#6B6B6B"; g.font = "12px Arial";
    g.fillText("no point has both a log2FC and a positive FDR", padL, H / 2);
    return;
  }
  let xa = 0, ym = 0;
  pts.forEach((p) => { xa = Math.max(xa, Math.abs(p.x)); ym = Math.max(ym, p.y); });
  xa = xa || 1; ym = ym || 1;
  const X = (v) => padL + ((v + xa) / (2 * xa)) * (W - padL - padR);
  const Y = (v) => H - padB - (v / ym) * (H - padT - padB);
  g.strokeStyle = "#1A1A1A"; g.lineWidth = 1;
  g.beginPath(); g.moveTo(padL, padT); g.lineTo(padL, H - padB); g.stroke();
  g.beginPath(); g.moveTo(padL, H - padB); g.lineTo(W - padR, H - padB); g.stroke();
  g.strokeStyle = "#DCDCDC";
  g.beginPath(); g.moveTo(X(0), padT); g.lineTo(X(0), H - padB); g.stroke();
  pts.forEach((p) => {
    g.fillStyle = p.x > 0 ? "#C44E52" : "#4C72B0";
    g.fillRect(X(p.x) - 1.3, Y(p.y) - 1.3, 2.6, 2.6);
  });
  g.fillStyle = "#1A1A1A"; g.font = "11px Arial";
  g.fillText("log2 fold change", W / 2 - 40, H - 8);
  g.save(); g.translate(13, H / 2 + 30); g.rotate(-Math.PI / 2);
  g.fillText("-log10 FDR", 0, 0); g.restore();
  g.fillStyle = "#6B6B6B"; g.font = "10px Arial";
  g.fillText((-xa).toFixed(1), padL, H - padB + 12);
  g.fillText(xa.toFixed(1), W - padR - 20, H - padB + 12);
  g.fillText(ym.toFixed(1), 16, padT + 8);
  g.fillText(`${pts.length} points`, W - 90, padT + 10);
}

/* ------------------------------------------------------------------- CCC */
async function loadCcc() {
  const d = await jget(`/api/dataset/${S.ds}/ccc`);
  if (!d.available) {
    $("cccBody").innerHTML = `<h3>Cell-cell communication</h3><div class="note">` +
      `Not available for this dataset: ${esc(d.reason || "no table in the bundle")}.<br>` +
      `Add one by rerunning precompute.py with <code>--ccc &lt;edges.tsv&gt;</code>.</div>`;
    return;
  }
  $("cccBody").innerHTML = `<h3>Cell-cell communication</h3>` +
    `<div class="note">${d.n_rows.toLocaleString()} edges; source <code>${esc(d.source)}</code></div>` +
    `<table><tr>${d.columns.map((c) => `<th>${esc(c)}</th>`).join("")}</tr>` +
    d.rows.slice(0, 500).map((r) => `<tr>${d.columns.map((c) => `<td>${esc(r[c])}</td>`).join("")}</tr>`).join("") +
    `</table>`;
}

/* ------------------------------------------------------------ provenance */
function renderAbout() {
  const m = S.meta;
  const cov = Object.entries(m.covariates).map(([k, v]) =>
    `<tr><td>${esc(k)}</td><td>${esc(v.kind)}</td><td class="num">` +
    `${v.kind === "categorical" ? v.categories.length + " categories" : fmt(v.min) + " to " + fmt(v.max)}` +
    `</td><td class="num">${(v.n_missing || 0).toLocaleString()}</td></tr>`).join("");
  $("aboutBody").innerHTML = `
    <h3>${esc(m.label)}</h3>
    <table>
      <tr><td>cells</td><td class="num">${m.n_cells.toLocaleString()}</td></tr>
      <tr><td>genes</td><td class="num">${m.n_genes.toLocaleString()}</td></tr>
      <tr><td>non-zero values stored</td><td class="num">${(m.nnz || 0).toLocaleString()}</td></tr>
      <tr><td>cell states present</td><td class="num">${m.states.length}</td></tr>
      <tr><td>cell states in the canonical order</td><td class="num">${m.n_states_in_canonical_order}</td></tr>
      <tr><td>canonical order source</td><td>${esc(m.canonical_order_source)}</td></tr>
      <tr><td>layer served</td><td>${esc(m.layer)}</td></tr>
      <tr><td>bundle built (UTC)</td><td>${esc(m.built_utc)}</td></tr>
      <tr><td>bundle directory</td><td><code>${esc(m.bundle_dir)}</code></td></tr>
      <tr><td>source h5ad</td><td><code>${esc(m.source_h5ad)}</code></td></tr>
    </table>
    <h4>Method</h4>
    <div class="note">embedding: ${esc(m.embedding_method)}<br>
      statistics: ${esc(m.stats_method)}<br>HVG: ${esc(m.hvg_method)}</div>
    <h4>Cell states in the canonical order but absent from this object
      (${m.states_in_canonical_order_not_present.length})</h4>
    <div class="note">${esc(m.states_in_canonical_order_not_present.join(", ")) || "none"}</div>
    <h4>Cell states with no colour (${m.states_missing_color.length})</h4>
    <div class="note">${esc(m.states_missing_color.join(", ")) || "none"}</div>
    <h4>Precompute warnings (${m.warnings.length})</h4>
    <div class="note ${m.warnings.length ? "warn" : ""}">${m.warnings.map(esc).join("<br>") || "none"}</div>
    <h4>Covariates (${Object.keys(m.covariates).length})</h4>
    <table><tr><th>column</th><th>kind</th><th>range</th><th>missing</th></tr>${cov}</table>`;
}

init().catch((e) => { status("startup failed: " + e.message); console.error(e); });
