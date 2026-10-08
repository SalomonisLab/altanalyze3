const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/visualization/scalable_viewer/static/viewer_bootstrap.js'), 'utf8');
let replies = [], polls = [], opened = [], clears = 0, ready = true, fields = [], writes = [], writeReplies = [];
const nodes = new Map(['results-job-id', 'upload-job-id', 'qc-job-id'].map(id => [id, {value: ''}]));
const ctx = vm.createContext({console, Promise,
  el: id => nodes.get(id), api: x => x, explorePayloadCache: {clear: () => clears++},
  resetExploreResultsReadiness() {}, fillContrastSelector() {}, installDifferentialModalitySwitch() {},
  pruneDifferentialModalities() {}, installBundleStateColors() {},
  getJson: url => new Promise(resolve => replies.push({url, resolve})),
  fillCovariateSelectors: values => { fields = values; },
  pollStatus: async id => { polls.push(id); }, areExploreResultsReady: () => ready,
  setExplorerTab: tab => opened.push(tab), setTimeout: callback => callback(),
  ensureExploreResultsReady: async () => true, setResultMode: mode => opened.push(mode),
  fetch: url => { writes.push(url); return new Promise(resolve => writeReplies.push(resolve)); },
});
vm.runInContext('const SV={datasets:[{id:"old"},{id:"new"}],current:null}; let datasetLoadRequest=0, contrastRequest=0; let contrastWritePromise=Promise.resolve();', ctx);
const begin = source.indexOf('  async function loadDataset(');
vm.runInContext(source.slice(begin, source.indexOf('  // ------------------------------------------------------------------ Study tab', begin)), ctx);
const contrastBegin = source.indexOf('  async function selectViewerContrast(');
vm.runInContext(source.slice(contrastBegin, source.indexOf('    /* The Differential', contrastBegin)), ctx);
const tick = () => new Promise(setImmediate);
(async () => {
  const old = ctx.loadDataset('old'); const current = ctx.loadDataset('new');
  replies[1].resolve({colors: {c1: 'new'}, cluster_key: 'new_state'}); await tick();
  replies[2].resolve({fields: ['new_state']}); await current;
  replies[0].resolve({colors: {c1: 'old'}, cluster_key: 'old_state'}); await old;
  assert.equal(vm.runInContext('SV.stateColors.c1', ctx), 'new');
  assert.deepEqual(fields, ['new_state']); assert.deepEqual(polls, ['new']);
  assert.deepEqual(opened, ['explore']); assert.equal(clears, 2);
  // A readiness timeout must not open a page whose data never became ready.
  ready = false; const unavailable = ctx.loadDataset('new');
  replies[3].resolve({colors: {c1: 'new'}}); await tick(); replies[4].resolve({fields: []}); await unavailable;
  assert.equal(opened.length, 1);
  // In-flight POSTs finish before the next POST starts; superseded queued work is skipped.
  const first = ctx.selectViewerContrast('new', 'first'); await tick();
  const middle = ctx.selectViewerContrast('new', 'middle');
  const last = ctx.selectViewerContrast('new', 'last');
  assert.equal(writes.length, 1);
  writeReplies[0]({ok: true}); await tick(); assert.equal(writes.length, 2);
  assert.ok(writes[1].endsWith('contrast=last'));
  writeReplies[1]({ok: true}); await Promise.all([first, middle, last]);
  assert.deepEqual(opened, ['explore', 'differential']);
  // A dataset switch invalidates a pending comparison's UI work.
  const stale = ctx.selectViewerContrast('new', 'stale'); await tick();
  vm.runInContext('datasetLoadRequest++', ctx); nodes.get('results-job-id').value = 'old';
  writeReplies[2]({ok: true}); await stale;
  assert.equal(opened.length, 2);
  console.log('Viewer races passed: dataset metadata ownership, readiness timeout and ordered latest comparison.');
})().catch(error => { console.error(error); process.exitCode = 1; });
