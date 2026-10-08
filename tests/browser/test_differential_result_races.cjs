const assert = require('node:assert/strict'), fs = require('node:fs'), vm = require('node:vm'), path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/cellHarmony/webapp/static/app.js'), 'utf8');
let annotations = [], mounted = 0, rendered = 0, resets = 0;
const nodes = new Map();
const node = id => {
  if (!nodes.has(id)) nodes.set(id, {value: '', classList: {toggle() {}, add() {}, remove() {}}});
  return nodes.get(id);
};
node('results-job-id').value = 'job'; node('differential-viz-mode').value = 'integrated_network';
node('differential-result-population').value = 'c1';
const ctx = vm.createContext({console, Promise, document: {getElementById: node},
  differentialVisualizationRequest: 0, differentialDetailRequest: 0,
  currentDifferentialState: {status: 'completed', run_id: 'run', config: {modality: 'rna'}},
  currentDifferentialGene: '', currentDifferentialPopulation: '', currentDifferentialFeatureRole: 'ligand',
  getResultsJobId: () => node('results-job-id').value,
  explorePayloadCache: {generation: 1}, differentialPayloadCache: {generation: 1},
  ensureFeatureAnnotations: () => new Promise(resolve => annotations.push(resolve)),
  resetDifferentialResults: () => resets++, updateDifferentialDownloadButton() {}, destroyDifferentialNetwork() {}, resetDifferentialGeneDetail() {},
  ScalableIntegrated: {mount: async () => mounted++}, integratedOptions: () => ({}), Plotly: {purge() {}},
  apiPath: url => url, resolveDifferentialGeneFilter: async () => null,
  fetchDifferentialJson: async () => ({}), renderDifferentialVolcano: () => rendered++,
  differentialPayloadGenes: () => [], setDifferentialGeneFilterOptions() {},
});
const start = source.indexOf('async function loadDifferentialVisualization('), end = source.indexOf('async function fetchDifferentialJson(', start);
vm.runInContext(source.slice(start, end), ctx);
(async () => {
  const obsolete = ctx.loadDifferentialVisualization();
  node('differential-viz-mode').value = 'volcano';
  const current = ctx.loadDifferentialVisualization();
  annotations[1](); await current; assert.equal(rendered, 1);
  annotations[0](); await obsolete;
  assert.equal(mounted, 0, 'late annotation readiness must not mount an obsolete integrated view');
  assert.equal(rendered, 1, 'an obsolete request cannot refetch/redraw the current mode');
  // A same-job/run reset also invalidates delayed annotation completion.
  const afterReset = ctx.loadDifferentialVisualization(); ctx.explorePayloadCache.generation++;
  annotations[2](); await afterReset; assert.equal(rendered, 1);
  // An old missing state must not reset a newer completed view.
  ctx.currentDifferentialState = null; const missing = ctx.loadDifferentialVisualization();
  ctx.currentDifferentialState = {status:'completed', run_id:'next'};
  const next = ctx.loadDifferentialVisualization(); annotations[4](); await next;
  annotations[3](); await missing; assert.equal(resets, 0);
  // Exercise production differential teardown, including an integrated renderer.
  let disposed = 0, purged = [], filterClears = 0, detailClears = 0;
  node('differential-plot-area')._integratedDispose = () => disposed++;
  Object.assign(ctx, {clearDifferentialGeneFilter: () => filterClears++, syncDifferentialGeneFilterLabel() {},
    resetDifferentialGeneDetail: () => detailClears++, Plotly: {purge: plot => purged.push(plot)}});
  const resetStart = source.indexOf('function resetDifferentialResults(');
  vm.runInContext(source.slice(resetStart, source.indexOf('function resetDifferentialGeneDetail(', resetStart)), ctx);
  const previousRequest = ctx.differentialVisualizationRequest;
  ctx.resetDifferentialResults();
  assert.equal(disposed, 1); assert.equal(purged.length, 1); assert.equal(filterClears, 1); assert.equal(detailClears, 1);
  assert.equal(ctx.differentialVisualizationRequest, previousRequest + 1);
  console.log('Differential request races passed: annotation delay, mode change, result reset and integrated teardown preserve current view.');
})().catch(error => { console.error(error); process.exitCode = 1; });
