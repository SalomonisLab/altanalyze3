const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/cellHarmony/scalable_discover/static/discover.js'), 'utf8');
const replies = [];
let calls = 0, updates = 0, renders = 0, messages = 0, mode = 'goelite_biomarkers', state = 'c1';
const request = () => { calls++; return new Promise(resolve => replies.push(resolve)); };
const ctx = vm.createContext({console, Promise,
  explorePayloadCache: {generation: 1, fetch: request}, exploreMetadataCache: {fetch: request},
  getResultsJobId: () => 'job', apiPath: x => x, updateExpressionModeOptions: () => updates++,
  panelVisualizationRequest: {viz1: 0}, panelPlotData: {},
  getPanelSelectValue: (_, field) => field === 'mode' ? mode : state,
  resetVisualizationSurface() {}, renderVisualizationMessage: () => messages++,
  renderGoeliteBiomarkers: () => renders++, setPanelSummary() {},
  loadVisualizationPanel: async panel => { ctx.panelVisualizationRequest[panel]++; },
});
vm.runInContext('const GOELITE_MODE="goelite_biomarkers"; let discoverStatus=null; let discoverGoeliteStates={key:"",states:[]}; let discoverGoeliteRequest=null;', ctx);
const begin = source.indexOf('async function loadGoeliteStates(');
vm.runInContext(source.slice(begin, source.indexOf('(function goeliteMenus()', begin)), ctx);
const plotBegin = source.indexOf('  const originalLoad = loadVisualizationPanel;');
vm.runInContext(source.slice(plotBegin, source.indexOf('  const originalDownload =', plotBegin)), ctx);
const data = {status: 'completed', icgs3_analysis: {goelite: {available: true}}};
(async () => {
  const old = ctx.loadGoeliteStates('job', data);
  const duplicate = ctx.loadGoeliteStates('job', data);
  assert.equal(calls, 1, 'status polls share the pending states request');
  ctx.explorePayloadCache.generation++;
  const current = ctx.loadGoeliteStates('job', data);
  replies[1]({states: ['current']}); await current;
  replies[0]({states: ['obsolete']}); await Promise.all([old, duplicate]);
  assert.equal(vm.runInContext('discoverGoeliteStates.states[0]', ctx), 'current');
  assert.equal(updates, 1);
  await ctx.loadGoeliteStates('job', data); assert.equal(calls, 2, 'completed states are reused');
  const oldPlot = ctx.loadVisualizationPanel('viz1');
  state = 'c2'; const newPlot = ctx.loadVisualizationPanel('viz1');
  replies[3]({population: 'c2', terms: [{term_name: 'new', is_selected_positive_sig: true}]}); await newPlot;
  replies[2]({population: 'c1', terms: [{term_name: 'old'}]}); await oldPlot;
  assert.equal(ctx.panelPlotData.viz1.payload.population, 'c2'); assert.equal(renders, 1);
  const leaving = ctx.loadVisualizationPanel('viz1'); mode = 'violin'; await ctx.loadVisualizationPanel('viz1');
  replies[4]({terms: []}); await leaving;
  assert.equal(messages, 0, 'late GO-Elite replies cannot overwrite a different plot mode');
  let clears = 0, removed = 0;
  ctx.explorePayloadCache.onClear = () => clears++;
  ctx.document = {getElementById: () => ({remove: () => removed++})};
  const clearBegin = source.indexOf('(function clearDiscoverResults()');
  vm.runInContext(source.slice(clearBegin, source.indexOf('window.scalableQcOptions', clearBegin)), ctx);
  ctx.applyJobStatus = () => { ctx.explorePayloadCache.onClear(); return 'applied'; };
  ctx.renderLayerControl = () => {};
  const statusBegin = source.indexOf('(function captureStatus()');
  vm.runInContext(source.slice(statusBegin, source.indexOf('// ------------------------------------------------------------- GO-Elite', statusBegin)), ctx);
  const processing = {status: 'processing'};
  assert.equal(ctx.applyJobStatus('job', processing), 'applied');
  assert.equal(vm.runInContext('discoverStatus.status', ctx), 'processing', 'the new status survives result invalidation');
  assert.equal(vm.runInContext('discoverGoeliteStates.states.length', ctx), 0);
  assert.equal(clears, 1); assert.equal(removed, 1);
  console.log('Discover result races passed: deduplicated states, rerun ownership and latest plot/state wins.');
})().catch(error => { console.error(error); process.exitCode = 1; });
