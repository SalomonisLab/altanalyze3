/* Network and WikiPathways renderers ported from LungMAP Discover, 2026-09-16.
 * Same geometry, colours, layout and node/edge scaling; shared vector PDF export.
 */
(() => {
  const esc = v => String(v ?? "").replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  const fmt = (v, n) => Number(v).toFixed(n);
  const LABEL_DARK = "#1b2733", LABEL_CROSSOVER = .223;
  function foldColour(v, span) {
    if (v === null || v === undefined) return "#b7c0c9";
    var t = Math.max(-1, Math.min(1, v / (span || 1)));
    var lo = [44, 111, 181], mid = [241, 243, 245], hi = [192, 57, 43];
    var a = t < 0 ? lo : hi, k = Math.abs(t);
    return "rgb(" + [0, 1, 2].map(function (i) {
      return Math.round(mid[i] + (a[i] - mid[i]) * k);
    }).join(",") + ")";
  }
  function lipidColour(v, span) {
    if (v === null || v === undefined) return "#b7c0c9";
    var t = Math.max(-1, Math.min(1, v / (span || 1)));
    // A LIGHTER GREEN. Nathan, 2026-09-11: "the green is too dark and can't be seen with
    // the black text." rgb(27,120,55) has relative luminance 0.16, so the near-black node
    // label sat at contrast 3.0:1 on it -- below the 4.5:1 a reader needs. rgb(64,150,84)
    // reaches 0.25 and stays unmistakably green. The purple end is lightened with it, for
    // the same reason and by the same amount; leaving it alone would have made the two
    // halves of one ramp disagree about how saturated a full-scale fold looks.
    var lo = [128, 86, 174], mid = [244, 243, 246], hi = [64, 150, 84];
    var a = t < 0 ? lo : hi, k = Math.abs(t);
    return "rgb(" + [0, 1, 2].map(function (i) {
      return Math.round(mid[i] + (a[i] - mid[i]) * k);
    }).join(",") + ")";
  }
  function labelOn(fill) {
    var m = /rgb\((\d+),\s*(\d+),\s*(\d+)\)/.exec(fill);
    if (!m) { return LABEL_DARK; }              // a hex fill: #ffffff and #eef1f4 are pale
    var lum = [1, 2, 3].map(function (i) {
      var c = parseInt(m[i], 10) / 255;
      return c <= 0.03928 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
    });
    var L = 0.2126 * lum[0] + 0.7152 * lum[1] + 0.0722 * lum[2];
    return L < LABEL_CROSSOVER ? "#ffffff" : LABEL_DARK;
  }
  function isLipidNode(n) {
    return /^(lipid|metabolite)/.test(String(n.measured_as || ""));
  }

let mountSequence = 0;
function modalityColour(value, mod) {
  const span=mod.scale?.maximum || 1;
  const t=Math.min(1,Math.abs(value)/span);
  const hex=value < 0 ? "#78a0cf" : mod.color;
  const rgb=[1,3,5].map(i=>parseInt(hex.slice(i,i+2),16));
  return `rgb(${rgb.map(v=>Math.round(247+(v-247)*t)).join(",")})`;
}
function modalityScaleSvg(mod, id) {
  const scale=mod.scale || {}, min=scale.minimum || 0, max=scale.maximum || 0;
  const number=v=>Number(v).toPrecision(3).replace(/\.?0+$/,m=>m.includes(".") ? "" : m);
  const stops=scale.available ? Array.from({length:11},(_,i)=>`<stop offset="${i/10}" stop-color="${modalityColour(min+(max-min)*i/10,mod)}"/>`).join('') : '<stop offset="0" stop-color="#eef1f4"/><stop offset="1" stop-color="#eef1f4"/>';
  return `<svg xmlns="http://www.w3.org/2000/svg" class="modality-colorbar" data-modality="${esc(mod.id)}" data-min="${scale.available?min:''}" data-max="${scale.available?max:''}" width="350" height="52" viewBox="0 0 350 52" role="img" aria-label="${esc(mod.label)}: ${scale.available ? `${esc(scale.statistic)}, ${min} to ${max}` : 'no mapped values'}">
    <defs><linearGradient id="${id}">${stops}</linearGradient></defs>
    <g font-family="Helvetica" font-size="11" fill="#1b2733"><text x="0" y="12">${esc(mod.label)} · ${esc(scale.statistic || '')}</text>
    <rect x="42" y="22" width="254" height="12" fill="url(#${id})" stroke="#aab4c0" stroke-width=".5"/>
    ${scale.available ? `<text x="36" y="32" text-anchor="end">${number(min)}</text><text x="301" y="32">${number(max)}</text>${min<0 ? '<text x="169" y="47" text-anchor="middle">0</text>' : ''}` : '<text x="169" y="47" text-anchor="middle">No mapped values</text>'}</g></svg>`;
}

async function mount(host, options) {
  host._integratedDispose?.();
  const instanceId=`integrated-${++mountSequence}`;
  const S = {}; let serial = 0;
  host._integratedDispose = () => {
    ++serial; S.cy?.destroy(); S.cy=null; S.ready=false;
    host._integratedPdf=null; host._integratedResize=null;
    // Dispose tore down cytoscape but left this panel's MARKUP in the host, so the
    // "Regulatory network" heading, its Selection/Show/Gene fold/Edge fold/TF expression
    // row and its "both differential expression and grn imputed differentials are
    // required" note stayed under whatever drew next. app.js:4017 disposes when the
    // Differential Explorer leaves an integrated mode, which is why an ADT volcano showed
    // that block beneath it. mount() rebuilds the markup at integrated.js:77, so clearing
    // here is safe; app.js:719 removed the class by hand for the same reason.
    host.classList.remove("integrated-view");
    host.innerHTML = "";
  };
  host._integratedResize = () => { S.cy?.resize(); };
  host.classList.add("integrated-view");
  const pathway = options.kind === "integrated_pathway";
  const marker = !pathway && options.source === "marker";
  const modality = options.modality;
  const select = (name, label, values, value) => `<label>${label} <select data-setting="${name}">${values.map(v => `<option value="${esc(v[0])}" ${String(v[0])===String(value)?"selected":""}>${esc(v[1])}</option>`).join("")}</select></label>`;
  const vals = a => a.map(v => [v,v]);
  host.innerHTML = `<div class="integrated-controls">` +
    (pathway || marker ? "" :
      select("significance","Selection",[["reported","Reported calls"],["fdr","FDR"],["pval","Raw p"]],options.significance || "reported") +
      select("max_fdr","p",vals([.001,.01,.05,.1,1]),.05)) +
    (pathway ? '<label>Pathway <select data-setting="id"></select></label>'  :
      select("limit","Show",vals([25,50,100,200,500,1000,2000]),options.limit || 50) +
      select("gene_fold",marker ? "Marker fold" : "Gene fold",vals([1,1.2,1.5,2]),1.2) +
      (marker ? select("min_score","|Edge score| ≥",vals([0,.01,.05,.1,.25,.5,1]),0)
              : select("min_fold","Edge fold",vals([1,1.2,1.5,2]),1.2)) +
      select("min_expression","TF expression",[[0,">0"],[.5,">0.5"],[1,">1"],[2,">2"]],.5)) +
    (options.externalPdf ? "" : '<button class="ghost-btn" type="button" data-action="pdf" disabled>Download PDF</button>') +
    (pathway ? '' : '<button class="ghost-btn" type="button" data-action="tsv">Download TSV</button>') +
    '</div><h3>'+ (pathway ? (options.crossModal ? 'Cross-modality pathway' : 'Lipid / metabolite and gene pathway') : 'Regulatory network') + '</h3><p data-role="why"></p><p data-role="note" role="status"></p><div data-role="figure"></div>';
  const el = id => host.querySelector(`[data-role="${id === "dsc-view-why" ? "why" : "note"}"]`);
  const figure=host.querySelector('[data-role="figure"]');
  const get = key => host.querySelector(`[data-setting="${key}"]`)?.value;
  function params() {
    const p = new URLSearchParams({contrast: options.contrast || "", cell_state: options.cell_state || "", source: options.source || "differential"});
    if(pathway)p.set("modality",modality || "");
    host.querySelectorAll('[data-setting]').forEach(e => p.set(e.dataset.setting,e.value));
    if (options.features?.length) p.set("features", options.features.join(","));
    return p;
  }
  function url(name, p=params()) { return (options.apiPath || (v=>v))(`/api/jobs/${encodeURIComponent(options.jobId)}/integrated/${name}?${p}`); }
  async function fetchData(name,p) { const r=await fetch(url(name,p)); const d=await r.json(); if(!r.ok) throw new Error(typeof d.detail === "string" ? d.detail : "Unable to read integrated results."); return d; }
  async function draw(reloadList=true) {
    const ticket=++serial; S.cy?.destroy(); S.cy=null; S.cyExport=null;
    S.data=null;S.ready=false;
    const pdfButton=host.querySelector('[data-action="pdf"]');if(pdfButton)pdfButton.disabled=true;
    figure.replaceChildren(); el('dsc-view-why').textContent='';
    const threshold=host.querySelector('[data-setting="max_fdr"]');
    if(threshold)threshold.closest('label').hidden=get('significance')==='reported';
    el('dsc-view-note').textContent=marker ? 'Reading RNA markers and cell-state GRN scores…' : 'Reading existing differential results…';
    try {
      if(pathway && !options.crossModal && !['lipid','lipids','metabolite'].includes(modality)) {
        el('dsc-view-note').textContent='Select a lipid or metabolite modality to view its pathways.';return;
      }
      const p=params();
      if(pathway && reloadList) {
        const list=await fetchData(options.crossModal?'cross-pathways':'pathways',p); if(ticket!==serial) return;
        if(!list.available) {host.querySelector('[data-setting="id"]').replaceChildren(); el('dsc-view-note').textContent=list.note;return; }
        const picker=host.querySelector('[data-setting="id"]');const old=picker.value || options.id;
        picker.innerHTML=(list.pathways || []).map(d=>`<option value="${esc(d.id)}">${esc(d.name)} (${options.crossModal ? `${d.modality_count} modalities, ${d.total_hits} features` : `${d.n_genes_painted || 0} genes, ${d.n_lipids_painted || 0} metabolites`})</option>`).join('');
        if([...picker.options].some(o=>o.value===old))picker.value=old;
        if(!picker.value){el('dsc-view-note').textContent=list.note || 'No pathway contains a regulated lipid or metabolite feature.';return;}
        p.set('id',picker.value);
      }
      const d=await fetchData(pathway?(options.crossModal?'cross-pathway':'pathway'):'network',p); if(ticket!==serial) return;
      el('dsc-view-note').textContent=d.note || d.error || '';
      if(d.available===false || !d.nodes?.length) return;
      S.data=d;
      if(pathway)drawPathwayDiagram(figure,d);else drawCytoscape(figure,d);
      S.ready=true;if(pdfButton)pdfButton.disabled=false;
    } catch(e) { if(ticket===serial)el('dsc-view-note').textContent=e.message; }
  }
  host.querySelectorAll('[data-setting]').forEach(e=>e.addEventListener('change',()=>draw(e.dataset.setting!=='id')));
  host._integratedPdf=async () => {
    if(!S.ready)throw new Error('There is no completed figure to export.');
    const markup=S.cyExport?.() || figure.querySelector(pathway ? '.dsc-scroll > svg' : 'svg')?.outerHTML;
    if(!markup)throw new Error('There is no figure to export.');
    const svg=new DOMParser().parseFromString(markup,'image/svg+xml').documentElement;
    const size=pathway ? {width:S.data.width,height:S.data.height} : getSvgIntrinsicSize(svg);
    svg.setAttribute('width',size.width);svg.setAttribute('height',size.height);
    // Export the complete diagram independently of on-screen zoom/scroll position.
    svg.setAttribute('viewBox',`0 0 ${size.width} ${size.height}`);
    const cross=S.data.cross_modal;
    const headerHeight=cross ? 88+Math.ceil(S.data.modalities.length/2)*54 : 110;
    svg.setAttribute('x','0');svg.setAttribute('y',String(headerHeight));
    const width=Math.max(730,size.width),height=size.height+headerHeight+6;
    const scaleLegend=cross ? S.data.modalities.map((m,i)=>{
      const legend=new DOMParser().parseFromString(modalityScaleSvg(m,`${instanceId}-pdf-scale-${i}`),'image/svg+xml').documentElement;
      legend.setAttribute('x',String(12+(i%2)*360));legend.setAttribute('y',String(76+Math.floor(i/2)*54));
      return new XMLSerializer().serializeToString(legend);
    }).join('') : '';
    const title=pathway ? S.data.name : 'Regulatory network';
    const context=[S.data.cell_state || options.cell_state, options.crossModal && options.source === 'marker' ? 'Cell-state markers' : marker ? 'RNA markers versus other cell states' : S.data.comparison || options.contrast].filter(Boolean).join(' — ');
    const payload=`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">
      <rect width="100%" height="100%" fill="white"/>
      <g font-family="Helvetica" fill="#1b2733"><text x="12" y="22" font-size="16">${esc(title)}</text>
      <text x="12" y="42" font-size="11">${esc(context)}</text>
      ${(cross ? ["Each modality scaled separately; strongest absolute feature score per stripe; gray = no mapped hit."] : S.pdfLegend).map((line,i)=>`<text x="12" y="${62+16*i}" font-size="10">${esc(line)}</text>`).join('')}</g>${scaleLegend}${new XMLSerializer().serializeToString(svg)}</svg>`;
    await saveSvgMarkupAsPdf(payload,buildPdfFilename([options.jobId,options.cell_state,title],pathway?'pathway':'regulatory_network'));
  };
  const pdfButton=host.querySelector('[data-action="pdf"]');
  if(pdfButton)pdfButton.onclick=async()=>{try{await host._integratedPdf();}catch(error){showDownloadError(error);}};
  const tsv=host.querySelector('[data-action="tsv"]');if(tsv)tsv.onclick=()=>{const a=document.createElement('a');a.href=url('network.tsv');a.click();};
  function drawCytoscape(node, d) {
    node.innerHTML = "";
    var span = 0;
    d.nodes.forEach(function (n) {
      var v = (n.log2fc === null || n.log2fc === undefined) ? n.activity_log2fc : n.log2fc;
      if (v !== null && v !== undefined) span = Math.max(span, Math.abs(v));
    });
    span = span || 1;
    var scores = d.edges.map(function (e) { return Math.abs(e.score); });
    var sLo = Math.min.apply(null, scores), sHi = Math.max.apply(null, scores);
    S.pdfLegend = [
      marker
        ? `Upregulated RNA markers: positive log2 fold, light to dark red (0 to ${span.toPrecision(3)}); gray = regulator without a qualifying positive marker fold.`
        : `Node differential log2FC: blue (${-span.toPrecision(3)}) to red (${span.toPrecision(3)}); gray = no reported differential call.`,
      'Diamond = factor; circle = target; arrow = regulatory edge.' + (marker ? '' : ' Dashed rim = activity fold.'),
      `Line width = |edge score| ${sLo.toPrecision(3)} to ${sHi.toPrecision(3)}.`
    ];

    var elements = [];
    d.nodes.forEach(function (n) {
      // EXPRESSION FIRST, ACTIVITY SECOND. A factor whose expression did not move but
      // whose activity did is exactly the case Nathan named, so it is coloured by its
      // activity and marked with a dashed rim rather than dropped or greyed.
      var byExpression = (n.log2fc !== null && n.log2fc !== undefined);
      var value = byExpression ? n.log2fc : n.activity_log2fc;
      var degree = n.kind === "regulator" ? (n.n_targets || 1)
                                          : d.edges.filter(function (e) {
                                              return e.target === n.id; }).length;
      elements.push({ data: {
        id: n.id, label: n.id, kind: n.kind,
        colour: foldColour(value, span),
        border: (!byExpression && value !== null && value !== undefined)
          ? "#5b6b7b" : "#8a949e",
        borderStyle: (!byExpression && value !== null && value !== undefined)
          ? "dashed" : "solid",
        size: n.kind === "regulator" ? Math.min(46, 20 + 4 * degree)
                                     : Math.min(52, 22 + 3.2 * degree),
        tip: n.id + " — " +
          (n.kind === "regulator"
            ? "regulates " + n.n_targets + " of the features shown"
            : "regulated by " + degree + " of the factors drawn") +
          (n.expression !== null && n.expression !== undefined
            ? "; expression " + fmt(n.expression, 2) + (marker ? " (" + d.expression_scale + ")" : " log2(1+CP10k)") : "") +
          (byExpression
            ? (marker ? "; marker log2 fold versus other cell states " : "; differential expression log2 ") + fmt(n.log2fc, 2) +
              " at FDR " + (n.fdr === null || n.fdr === undefined ? "not reported" : Number(n.fdr).toExponential(1))
            : marker ? "; no qualifying positive marker fold" : "; no significant differential expression here") +
          ((n.activity_log2fc !== null && n.activity_log2fc !== undefined)
            ? "; activity log2 " + fmt(n.activity_log2fc, 2) : "")
      }});
    });
    d.edges.forEach(function (e, i) {
      elements.push({ data: {
        id: "e" + i, source: e.source, target: e.target,
        // EDGE THICKNESS IS THE EDGE SCORE. It carried a range of 0.5 to 2.1 pixels
        // before, which no reader could see. The drawn range now spans 0.7 to 5.
        width: 0.7 + 4.3 * ((Math.abs(e.score) - sLo) / ((sHi - sLo) || 1)),
        // Nathan, 2026-09-11: "Hovering over the edge should show the edge differentials."
        // The edge score is the model's STATIC strength and says nothing about this
        // comparison. The fold and the p are what the gate actually admitted the edge on,
        // so the hover reports them first and the score second. A missing fold is stated
        // in words, never left blank, so "not differential" cannot be read as "not looked
        // up".
        tip: e.source + " \u2192 " + e.target +
             (marker ? ", cell-state predicted regulatory activity" : (e.log2fc === null || e.log2fc === undefined)
               ? ", edge differential not available"
               : ", edge log2FC " + fmt(e.log2fc, 3) + ", " +
                 (e.significance || "FDR") + " " +
                 ((e.fdr === null || e.fdr === undefined)
                   ? "n/a" : Number(e.fdr).toExponential(2))) +
             ", edge score " + fmt(e.score, 3)
      }});
    });

    var frame = document.createElement("div");
    frame.className = "dsc-cy";
    frame.style.height = Math.max(420, Math.min(760, 26 * d.nodes.length)) + "px";
    node.appendChild(frame);

    var cy = cytoscape({
      container: frame,
      elements: elements,
      style: [
        { selector: "node", style: {
            "background-color": "data(colour)",
            "border-color": "data(border)",
            "border-width": 1.4,
            "border-style": "data(borderStyle)",
            "width": "data(size)", "height": "data(size)",
            "label": "data(label)", "font-size": "11px",
            "color": "#1b2733", "text-valign": "center",
            "text-halign": "center", "text-outline-color": "#ffffff",
            "text-outline-width": 2 } },
        { selector: 'node[kind = "regulator"]', style: { "shape": "diamond" } },
        { selector: 'node[kind = "target"]', style: { "shape": "ellipse" } },
        { selector: "edge", style: {
            "width": "data(width)", "line-color": "#a8b4c0",
            "curve-style": "bezier", "opacity": 0.75,
            "target-arrow-shape": "triangle", "arrow-scale": 0.6,
            "target-arrow-color": "#a8b4c0" } },
        { selector: "node:selected", style: { "border-color": "#14315c",
            "border-width": 3 } }
      ],
      layout: { name: "cose", animate: false, nodeRepulsion: 9000,
                idealEdgeLength: 90, nodeOverlap: 14, gravity: 0.5,
                numIter: 900, randomize: false, fit: true, padding: 24 },
      wheelSensitivity: 0.2
    });
    // A VECTOR EXPORT, NOT A SCREENSHOT. Cytoscape draws to a canvas, and the page's
    // PDF download keeps vector paths and editable text. The laid-out positions are readable from the graph,
    // so the same figure is written out as real paths and circles.
    S.cyExport = function () {
      var pad = 30;
      var xs = cy.nodes().map(function (n) { return n.position("x"); });
      var ys = cy.nodes().map(function (n) { return n.position("y"); });
      var x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs);
      var y0 = Math.min.apply(null, ys), y1 = Math.max.apply(null, ys);
      var w = (x1 - x0) + 2 * pad + 80, h = (y1 - y0) + 2 * pad + 40;
      var out = ['<svg xmlns="http://www.w3.org/2000/svg" width="' + Math.round(w) +
                 '" height="' + Math.round(h) + '">',
                 '<rect width="100%" height="100%" fill="#ffffff"/>'];
      function px(v) { return (v - x0 + pad + 40).toFixed(1); }
      function py(v) { return (v - y0 + pad + 20).toFixed(1); }
      cy.edges().forEach(function (e) {
        var a = e.source().position(), b = e.target().position();
        const dx=b.x-a.x,dy=b.y-a.y,len=Math.hypot(dx,dy)||1;
        const radius=Number(e.target().data("size"))/2;
        const tip={x:b.x-dx/len*radius,y:b.y-dy/len*radius};
        const back={x:tip.x-dx/len*6,y:tip.y-dy/len*6};
        out.push(`<path d="M${px(tip.x)},${py(tip.y)} L${px(back.x-dy/len*3)},${py(back.y+dx/len*3)} L${px(back.x+dy/len*3)},${py(back.y-dx/len*3)} Z" fill="#a8b4c0"/>`);
        out.push('<path d="M' + px(a.x) + "," + py(a.y) + "L" + px(tip.x) + "," +
                 py(tip.y) + '" stroke="#a8b4c0" stroke-width="' +
                 Number(e.data("width")).toFixed(2) + '" fill="none" opacity="0.75"/>');
      });
      cy.nodes().forEach(function (n) {
        var p = n.position(), r = Number(n.data("size")) / 2;
        if (n.data("kind") === "regulator") {
          out.push('<path d="M' + px(p.x) + "," + (py(p.y) - r).toFixed(1) +
                   "L" + (Number(px(p.x)) + r).toFixed(1) + "," + py(p.y) +
                   "L" + px(p.x) + "," + (Number(py(p.y)) + r).toFixed(1) +
                   "L" + (Number(px(p.x)) - r).toFixed(1) + "," + py(p.y) +
                   'Z" fill="' + n.data("colour") + '" stroke="' + n.data("border") +
                   '" stroke-width="1.4"' + (n.data("borderStyle")==="dashed" ? ' stroke-dasharray="3 2"' : "") + '/>');
        } else {
          out.push('<circle cx="' + px(p.x) + '" cy="' + py(p.y) + '" r="' +
                   r.toFixed(1) + '" fill="' + n.data("colour") + '" stroke="' +
                   n.data("border") + '" stroke-width="1.4"/>');
        }
        out.push('<text x="' + px(p.x) + '" y="' + (Number(py(p.y)) + 3.5).toFixed(1) +
                 '" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" ' +
                 'font-size="10" fill="#1b2733">' + esc(n.id()) + "</text>");
      });
      out.push("</svg>");
      return out.join("\n");
    };
    cy.nodes().forEach(function (n) { n.qtip = n.data("tip"); });
    cy.on("mouseover", "node, edge", function (ev) {
      frame.setAttribute("title", ev.target.data("tip") || "");
    });
    cy.on("tap", "node", ev => options.onFeature?.(ev.target.id(), "rna"));
    cy.on("tap", "edge", ev => options.onFeature?.(ev.target.data("source"), "rna"));
    S.cy = cy;

    var legend = document.createElement("div");
    legend.className = "dsc-legend";
    var stops = [];
    for (var q = 0; q <= 10; q++) {
      stops.push(foldColour(marker ? span * (q / 10) : -span + 2 * span * (q / 10), span) + " " + (q * 10) + "%");
    }
    legend.innerHTML =
      '<span class="dsc-legend-lab">' + (marker ? 'upregulated markers: log2 fold vs other states' : 'log2 fold change') + '</span>' +
      '<span class="dsc-legend-num">' + fmt(marker ? 0 : -span, 1) + "</span>" +
      '<span class="dsc-legend-bar" style="background:linear-gradient(to right,' +
      stops.join(",") + ')"></span>' +
      '<span class="dsc-legend-num">' + fmt(span, 1) + "</span>" +
      '<span class="dsc-legend-lab"><i style="background:#b7c0c9"></i>' +
      (marker ? "no qualifying positive marker fold</span>" : "no significant result here</span>") +
      '<span class="dsc-legend-lab">diamond = factor, circle = regulated feature</span>' +
      (marker ? '' : '<span class="dsc-legend-lab">dashed rim = coloured by activity, not expression</span>') +
      '<span class="dsc-legend-lab">line width = |edge score| ' + fmt(sLo, 2) +
      " to " + fmt(sHi, 2) + "</span>";
    node.insertBefore(legend, node.firstChild);

    el("dsc-view-why").textContent = marker
      ? `Upregulated RNA markers for ${d.cell_state} and their expressed regulators. Colour is each gene's marker fold versus other cell states. Edges show the model's state-specific scores. Click a node to explore expression; drag to rearrange.`
      :
      "The features this query returned, and the transcription factors the " +
      d.cell_state + " regulatory model puts above them. Colour is each node's own log2 " +
      "fold change in this comparison. Click a feature to explore its expression; click a factor or " +
      "an edge for that factor's expression. Drag a node to pull it out.";
    el("dsc-view-note").textContent = d.note || "";
  }
  function drawPathwayDiagram(node, d) {
        // EACH MEASUREMENT IS SCALED BY ITS OWN SPAN. A gene fold change and a lipid
        // class fold change are different quantities; one shared maximum would make
        // whichever moves less look flat.
        var geneSpan = 0, lipidSpan = 0, nGene = 0, nLipid = 0;
        d.nodes.forEach(function (n) {
          if (n.log2fc === null || n.log2fc === undefined) { return; }
          if (isLipidNode(n)) {
            lipidSpan = Math.max(lipidSpan, Math.abs(n.log2fc)); nLipid += 1;
          } else {
            geneSpan = Math.max(geneSpan, Math.abs(n.log2fc)); nGene += 1;
          }
        });
        geneSpan = geneSpan || 1;
        lipidSpan = lipidSpan || 1;
        S.pdfLegend = [
          `Gene log2FC: blue (${-geneSpan.toPrecision(3)}) to red (${geneSpan.toPrecision(3)}).`,
          `Lipid/metabolite log2FC: purple (${-lipidSpan.toPrecision(3)}) to green (${lipidSpan.toPrecision(3)}).`,
          'Rectangle = gene; rounded = metabolite; dashed = no reported call; gray = unmeasured.'
        ];
        node.innerHTML = "";
        var legend = document.createElement("div");
        legend.className = "dsc-legend";
        var ramp = function (fn, sp) {
          var out = [];
          for (var q = 0; q <= 10; q++) {
            out.push(fn(-sp + 2 * sp * (q / 10), sp) + " " + (q * 10) + "%");
          }
          return out.join(",");
        };
        legend.innerHTML =
          '<span class="dsc-legend-lab">gene, log2 fold change (' + nGene + ')</span>' +
          '<span class="dsc-legend-num">' + fmt(-geneSpan, 1) + "</span>" +
          '<span class="dsc-legend-bar" style="background:linear-gradient(to right,' +
          ramp(foldColour, geneSpan) + ')"></span>' +
          '<span class="dsc-legend-num">' + fmt(geneSpan, 1) + "</span>" +
          '<span class="dsc-legend-lab">lipid class / metabolite, log2 fold change (' + nLipid + ')</span>' +
          '<span class="dsc-legend-num">' + fmt(-lipidSpan, 1) + "</span>" +
          '<span class="dsc-legend-bar" style="background:linear-gradient(to right,' +
          ramp(lipidColour, lipidSpan) + ')"></span>' +
          '<span class="dsc-legend-num">' + fmt(lipidSpan, 1) + "</span>" +
          '<span class="dsc-legend-lab">'
            + '<i style="background:#ffffff;border:1px dashed #8a97a4"></i>' +
          "measured here, no reported differential call</span>" +
          '<span class="dsc-legend-lab"><i style="background:#eef1f4"></i>' +
          "not measured by this atlas</span>" +
          '<span class="dsc-legend-lab">rectangle = gene, rounded = metabolite</span>';
        if (d.cross_modal) {
          S.pdfLegend = ['Each modality has its own numeric scale; strongest absolute mapped score per stripe.',
            'Gray = no mapped retained hit; hover for individual feature scores.'];
          legend.classList.add('cross-modality-legend');
          legend.innerHTML=d.modalities.map((m,i)=>modalityScaleSvg(m,`${instanceId}-scale-${i}`)).join('');
        }
        node.appendChild(legend);

        // ZOOM. Nathan, 2026-09-11: "I would like to zoom into the pathway view."
        // WikiPathways lays a map out at its own coordinates, and the whole map was drawn
        // at a fixed 0.55 so it would fit the card. On a 165-metabolite overview that put
        // node labels at about 5.5 pixels, which is not readable at all.
        //
        // The viewBox is the pathway's own coordinate system and never changes; only the
        // rendered width and height do, so zooming is pure vector scaling. Nothing is
        // re-fetched and nothing is re-laid-out -- the same SVG is simply asked for more
        // pixels, and .dsc-scroll pans it. S.pwZoom persists across redraws so changing
        // the comparison or the cell type does not throw away the reader's zoom.
        var bar = document.createElement("div");
        bar.className = "dsc-zoom";
        var svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        svg.setAttribute("viewBox", "0 0 " + d.width + " " + d.height);

        var frame = document.createElement("div");
        frame.className = "dsc-scroll dsc-scroll-xy";

        var readout = document.createElement("span");
        readout.className = "dsc-zoom-pct";

        function applyZoom(z) {
          S.pwZoom = Math.max(0.2, Math.min(6, z));
          svg.setAttribute("width", Math.round(d.width * S.pwZoom));
          svg.setAttribute("height", Math.round(d.height * S.pwZoom));
          readout.textContent = Math.round(S.pwZoom * 100) + "%";
        }
        function button(text, title, fn) {
          var b = document.createElement("button");
          b.type = "button"; b.className = "dsc-zoom-b";
          b.textContent = text; b.title = title;
          b.addEventListener("click", fn);
          return b;
        }
        bar.appendChild(document.createTextNode("Zoom "));
        bar.appendChild(button("−", "Zoom out", function () {
          applyZoom(S.pwZoom / 1.3);
        }));
        bar.appendChild(readout);
        bar.appendChild(button("+", "Zoom in", function () {
          applyZoom(S.pwZoom * 1.3);
        }));
        bar.appendChild(button("Fit", "Scale the whole pathway to the width of the card",
          function () {
            applyZoom((frame.clientWidth - 4) / d.width);
          }));
        bar.appendChild(button("1:1", "Draw at the pathway's own coordinates",
          function () { applyZoom(1); }));
        node.appendChild(bar);
        node.appendChild(frame);
        frame.appendChild(svg);

        // Ctrl or Command with the wheel zooms, which is what a reader already does in
        // every map. A bare wheel keeps scrolling the page, so the pathway never traps it.
        frame.addEventListener("wheel", function (ev) {
          if (!ev.ctrlKey && !ev.metaKey) { return; }
          ev.preventDefault();
          applyZoom(S.pwZoom * (ev.deltaY < 0 ? 1.12 : 1 / 1.12));
        }, { passive: false });

        applyZoom(S.pwZoom || 0.55);

        // ONE 2-POINT PATH PER SEGMENT, never a sampled polyline.
        (d.edges || []).forEach(function (e) {
          for (var i = 0; i + 1 < e.points.length; i++) {
            var a = e.points[i], b = e.points[i + 1];
            var line = document.createElementNS("http://www.w3.org/2000/svg", "line");
            line.setAttribute("x1", a.x); line.setAttribute("y1", a.y);
            line.setAttribute("x2", b.x); line.setAttribute("y2", b.y);
            line.setAttribute("stroke", "#9aa7b3");
            line.setAttribute("stroke-width", "1.6");
            svg.appendChild(line);
          }
        });

        (d.labels || []).forEach(function (t) {
          var text = document.createElementNS("http://www.w3.org/2000/svg", "text");
          text.setAttribute("x", t.x); text.setAttribute("y", t.y);
          text.setAttribute("text-anchor", "middle");
          text.setAttribute("dominant-baseline", "central");
          text.setAttribute("font-size", Math.max(10, t.size || 10));
          text.setAttribute("fill", "#6b7885");
          text.setAttribute("font-style", "italic");
          text.textContent = t.text;
          svg.appendChild(text);
        });

        d.nodes.forEach(function (n) {
          var known = d.cross_modal ? !!n.hits?.length : (n.log2fc !== null && n.log2fc !== undefined);
          var wrap = svg;
          if (known) {
            var link = document.createElementNS("http://www.w3.org/2000/svg", "a");
            // A PATHWAY MIXES MEASUREMENTS, so each node links under ITS OWN one. A
            // metabolite opens the lipid bundle; a gene node opens gene expression even
            // while the finder ranks lipids, because the node IS a gene. `data-kind`
            // states which, so a check can test the pair rather than assume the view
            // carries a single measurement.
            link.setAttribute("data-kind",
              n.type === "Metabolite" ? "metabolite" : "gene");
            // LINK TO A FEATURE THAT EXISTS. A metabolite node is labelled with a head
            // group -- "Galactosylceramide" -- and this atlas measures species, so the
            // label itself opened a page reading "not measured in either atlas". The
            // server names the measured species the colour came from; a node with no
            // such feature is drawn without a link rather than with a dead one.
            link.setAttribute("href", "#");
            link.addEventListener("click", ev => { ev.preventDefault(); options.onFeature?.(n.link_feature || n.label, n.link_modality || (n.type === "Metabolite" ? modality : "rna")); });
            svg.appendChild(link);
            wrap = link;
          }
          var box = document.createElementNS("http://www.w3.org/2000/svg", "rect");
          box.setAttribute("x", n.x - n.w / 2);
          box.setAttribute("y", n.y - n.h / 2);
          box.setAttribute("width", n.w);
          box.setAttribute("height", n.h);
          box.setAttribute("rx", n.type === "Metabolite" ? Math.min(10, n.h / 2) : 2);
          // THREE STATES, NOT TWO. Filled means a result here; outlined means the
          // atlas measures it and this comparison did not move it; flat grey means
          // the atlas does not measure it at all.
          var fill = known
            ? (isLipidNode(n) ? lipidColour(n.log2fc, lipidSpan)
                              : foldColour(n.log2fc, geneSpan))
            : (n.in_atlas ? "#ffffff" : "#eef1f4");
          if (d.cross_modal && known) {
            const mods=d.modalities.filter(m=>n.hits.some(h=>h.modality===m.id));
            const gradient=document.createElementNS("http://www.w3.org/2000/svg","linearGradient");
            const id=instanceId+"-cross-node-"+d.id+"-"+d.nodes.indexOf(n);gradient.id=id;
            mods.forEach((m,i)=>[i/mods.length,(i+1)/mods.length].forEach(at=>{
              const stop=document.createElementNS("http://www.w3.org/2000/svg","stop");stop.setAttribute("offset",at);stop.setAttribute("stop-color",modalityColour(n.modality_values[m.id].value,m));gradient.appendChild(stop);
            }));svg.appendChild(gradient);fill=`url(#${id})`;
          }
          box.setAttribute("fill", fill);
          if (!known && n.in_atlas) {
            box.setAttribute("stroke", "#8a97a4");
            box.setAttribute("stroke-dasharray", "3 2");
          }
          box.setAttribute("stroke", known ? "#5b6b7b" : "#c6ced6");
          box.setAttribute("stroke-width", known ? "1.4" : "0.9");
          var tip = document.createElementNS("http://www.w3.org/2000/svg", "title");
          tip.textContent = n.label + " (" + (n.type || "node") + ")" +
            (known ? ": log2 " + fmt(n.log2fc, 2) + ", " + n.measured_as
                   : (n.in_atlas ? ": measured, no reported differential call" : ": this atlas does not measure it"));
          if (d.cross_modal) tip.textContent=n.label + (known ? "\n"+n.hits.map(h=>`${h.modality}: ${h.feature} — ${h.statistic} ${Number(h.value).toPrecision(4)}${h.fdr != null ? `; FDR ${h.fdr}` : ""}`).join("\n") : ": no mapped retained hit");
          box.appendChild(tip);
          wrap.appendChild(box);
          var label = document.createElementNS("http://www.w3.org/2000/svg", "text");
          label.setAttribute("x", n.x); label.setAttribute("y", n.y);
          label.setAttribute("text-anchor", "middle");
          label.setAttribute("dominant-baseline", "central");
          label.setAttribute("font-size", "10");
          // The label reads against ITS OWN node's fill, not against one colour chosen
          // for the whole figure. See labelOn().
          label.setAttribute("fill", labelOn(fill));
          label.textContent = n.label.length > 18
            ? n.label.slice(0, 17) + "\u2026" : n.label;
          label.style.pointerEvents="none";
          wrap.appendChild(label);
        });

        el("dsc-view-why").textContent =
          d.name + ", drawn at WikiPathways own coordinates. Colour is the log2 fold " +
          "change of the stored differential calls in the selected comparison and cell type. " +
          "Outlined nodes are measured without a reported call; grey nodes are unmeasured. Click a coloured node for its own page.";
        if (d.cross_modal) el("dsc-view-why").textContent=`${d.name}, at WikiPathways coordinates. ${d.cell_state} · ${d.source === "marker" ? "cell-state markers" : d.comparison}. Each stripe uses its modality’s color scale (strongest absolute retained feature score). Hover for feature IDs and all scores. Gray means no mapped retained hit.`;
        el("dsc-view-note").textContent = d.note || "";
  }
  await draw();
}
window.ScalableIntegrated={mount};
})();
