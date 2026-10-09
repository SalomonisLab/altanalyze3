const assert = require('node:assert/strict'), fs = require('node:fs'), vm = require('node:vm'), path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/cellHarmony/webapp/static/app.js'), 'utf8');
const template = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/cellHarmony/webapp/templates/index.html'), 'utf8');
const selectors = Object.fromEntries(['viz1', 'viz2'].map(key => [`${key}-association-x`, {
  value: 'unsupervised_cluster', replaceChildren(...options) { this.options = options; }
}]));
const plots = [];
const ctx = vm.createContext({ document: {
  getElementById: id => selectors[id], createElement: () => ({})
}, panelElementId: (key, suffix) => `${key}-${suffix}`, panelPlotId: key => key,
  Plotly: {newPlot: (...args) => plots.push(args)}, setPanelSummary: () => {} });
const start = source.indexOf('function renderClusterAssociations(');
vm.runInContext(source.slice(start, source.indexOf('function renderGoeliteBiomarkers(', start)), ctx);
const fields = [{value:'unsupervised_cluster', label:'Unsupervised clusters'},
  {value:'unsupervised_state', label:'Unsupervised cell states (GO-Elite)'}];
for (const key of ['viz1', 'viz2']) {
  assert.ok(template.includes(`id="${key}-association-x"`));
  const payload = {rows:[{cluster:'Macrophage_c1', state:'Unaligned', cells:2, cluster_cells:2, percent:100}],
    clusters:['Macrophage_c1'], states:['Unaligned'], cells:2,
    x_by:'unsupervised_state', x_label:fields[1].label, x_fields:fields};
  ctx.renderClusterAssociations(key, payload);
  assert.equal(selectors[`${key}-association-x`].value, 'unsupervised_state');
  assert.equal(selectors[`${key}-association-x`].options.length, 2);
  const [id, traces, layout] = plots.at(-1);
  assert.equal(id, key);
  assert.equal(traces[0].x[0], 'Macrophage_c1');
  assert.equal(traces[0].y[0], 'Unaligned');
  assert.equal(layout.xaxis.title.text, fields[1].label);
  assert.equal(layout.xaxis.tickmode, 'array');
  assert.equal(layout.xaxis.ticktext[0], 'Macrophage_c1');
  ctx.renderClusterAssociations(key, {...payload, x_by:'unsupervised_cluster', x_label:fields[0].label,
    rows:[{...payload.rows[0], cluster:'c1'}], clusters:['c1']});
  assert.equal(selectors[`${key}-association-x`].value, 'unsupervised_cluster');
  assert.equal(plots.at(-1)[1][0].x[0], 'c1');
}
const requestStart = source.indexOf('if (mode === "cluster_associations") {', source.indexOf('async function loadVisualizationPanel('));
assert.ok(source.slice(requestStart, requestStart+450).includes('params.set("x_by", getPanelSelectValue(panelKey, "association-x")'));
console.log('Cluster association axis regression passed: both panels, inferred names, cluster IDs, Unaligned counts and request selection.');
