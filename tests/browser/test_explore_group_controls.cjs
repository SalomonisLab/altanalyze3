// Exercise request ordering and selection retention in the production loaders.
const assert = require('node:assert/strict');
const fs = require('node:fs'), vm = require('node:vm'), path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/cellHarmony/webapp/static/app.js'), 'utf8');
function select() {
  return {options: [], value: '', dataset: {}, listeners: {},
    set innerHTML(value) { this.options = []; },
    appendChild(option) { this.options.push(option); },
    addEventListener(name, handler) { this.listeners[name] = handler; },
    get selectedOptions() { return this.options.filter(option => option.selected); }};
}
const groupBy = select(), groups = select();
const elements = {'results-job-id': {value: 'job'}, 'viz1-groupby': groupBy, 'viz1-groups': groups};
let resolveMetadata, requests = [];
const metadata = new Promise(resolve => {resolveMetadata = resolve;});
const ctx = vm.createContext({console, URLSearchParams,
  document: {getElementById: id => elements[id], createElement: () => ({})},
  panelElementId: (panel, suffix) => `${panel}-${suffix}`,
  loadPlotVariables: () => metadata, getResultsJobId: () => elements['results-job-id'].value, getPanelSelectValue: () => 'dotplot',
  panelModality: () => 'rna', panelVisualizationRequest: {viz1: 0},
  GENE_SET_MODES: new Set(['dotplot', 'combplot']), panelGeneSet: () => [],
  appendGeneSetSubsetParams: () => {}, panelPlotData: {}, renderVisualizationPanel: () => {},
  explorePayloadCache: {generation: 1, fetch: async url => {requests.push(url); return {}; }}, apiPath: url => url,
});
for (const [start, end] of [
  ['async function refreshGroupControls(', '/* The Color by'],
  ['function appendGeneSetGroupParams(', 'function panelGroupParams('],
  ['async function loadVisualizationPanel(', 'function renderVisualizationMessage('],
]) {
  const offset = source.indexOf(start);
  vm.runInContext(source.slice(offset, source.indexOf(end, offset)), ctx);
}
(async () => {
  const first = ctx.loadVisualizationPanel('viz1');
  assert.equal(requests.length, 0, 'first plot must wait for its default grouping');
  resolveMetadata({cluster_key: 'state', variables: [{field: 'state', n: 2, values: ['c0', 'c1']}]});
  await first;
  assert.equal(requests[0], '/api/jobs/job/dotplot?modality=rna&group_by=state');
  groups.options[1].selected = true;
  await ctx.loadVisualizationPanel('viz1');
  assert.equal(requests[1], '/api/jobs/job/dotplot?modality=rna&group_by=state&groups=c1');
  await ctx.refreshGroupControls('viz1');
  assert.equal(groups.selectedOptions.length, 1, 'switching views must preserve the group filter');
  assert.equal(groups.selectedOptions[0].value, 'c1');
  ctx.explorePayloadCache.generation++;
  ctx.loadPlotVariables=async()=>({cluster_key:'new_state',variables:[{field:'new_state',n:2,values:['new0','new1']}]});
  await ctx.refreshGroupControls('viz1');
  assert.equal(groupBy.value,'new_state','same-job rerun must replace stale grouping metadata');
  assert.equal(groups.options.map(o=>o.value).join(','),'new0,new1');
  const currentOptions=groups.options;
  let staleReply;
  ctx.explorePayloadCache.generation++;
  ctx.loadPlotVariables=()=>new Promise(resolve=>{staleReply=resolve});
  const stale=ctx.refreshGroupControls('viz1');
  ctx.explorePayloadCache.generation++;
  staleReply({cluster_key:'obsolete',variables:[{field:'obsolete',n:1,values:['wrong']}]});
  await stale;
  assert.equal(groups.options,currentOptions,'obsolete metadata cannot overwrite controls');
  // Axis/control metadata must also refresh when the same job changes layers.
  for(const suffix of ['colorby','coords','xfield','yfield'])elements[`viz1-${suffix}`]=select();
  ctx.syncUmapAxisFields=()=>{};ctx.formatAxisNumber=String;
  const uStart=source.indexOf('async function refreshUmapOptions('),uEnd=source.indexOf('function syncUmapAxisFields(',uStart);
  vm.runInContext(source.slice(uStart,uEnd),ctx);
  const umapInfo=min=>({cluster_key:'state',variables:[],coords:[{key:'',label:'cellHarmony UMAP'}],numeric_variables:[{field:'umap_1',label:'UMAP 1',min,max:min+1},{field:'umap_2',label:'UMAP 2',min,max:min+1}]});
  ctx.loadPlotVariables=async()=>umapInfo(1);
  await ctx.refreshUmapOptions('viz1');
  assert.match(elements['viz1-xfield'].options[0].textContent,/1 to 2/);
  ctx.explorePayloadCache.generation++;ctx.loadPlotVariables=async()=>umapInfo(100);
  await ctx.refreshUmapOptions('viz1');
  assert.match(elements['viz1-xfield'].options[0].textContent,/100 to 101/,'axis ranges must update on same-job rerun');
  // A temporary metadata error must not leave an unrecoverable placeholder.
  ctx.explorePayloadCache.generation++;
  ctx.loadPlotVariables=async()=>{throw new Error('temporary unavailable')};
  await ctx.refreshGroupControls('viz1');
  ctx.loadPlotVariables=async()=>({cluster_key:'recovered',variables:[{field:'recovered',n:1,values:['r0']}]});
  await ctx.refreshGroupControls('viz1');assert.equal(groupBy.value,'recovered');
  console.log('Group control regression passed: first request waits and subsequent selections are retained.');
})().catch(error => {console.error(error); process.exitCode = 1;});
