// Exercise the production cache without browser/network dependencies.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/cellHarmony/webapp/static/app.js'), 'utf8');
const start = source.indexOf('class ExplorePayloadCache');
const end = source.indexOf('const explorePayloadCache', start);
let calls = [], responder, now = 0;
let lastSignal;
const ctx = vm.createContext({AbortController,performance: {now: () => now},
  fetch: async (url,options) => {calls.push(url);lastSignal=options.signal; return responder(url);},
  setTimeout: callback => {queueMicrotask(callback);}});
vm.runInContext(source.slice(start, end) + '\nthis.Cache = ExplorePayloadCache;', ctx);
const ok = (value, status = 200) => ({ok: true, status, text: async () => JSON.stringify(value)});

(async () => {
  const cache = new ctx.Cache(1000, 2, 10);
  let release;
  responder = () => new Promise(resolve => {release = resolve;});
  const one = cache.fetch('/job?gene=A'), two = cache.fetch('/job?gene=A');
  assert.equal(calls.length, 1);
  release(ok({gene: 'A'}));
  assert.equal(await one, await two);
  await cache.fetch('/job?gene=A');
  assert.equal(calls.length, 1);
  responder = url => ok({url});
  await cache.fetch('/job?gene=B');
  await cache.fetch('/job?gene=A'); // make A recent
  await cache.fetch('/job?gene=C');
  assert.equal(cache.get('/job?gene=B'), null);
  assert.notEqual(cache.get('/job?gene=A'), null);
  now = 11;
  assert.equal(cache.get('/job?gene=A'), null);
  cache.put('oversized', {}, 1001);
  assert.equal(cache.get('oversized'), null);
  cache.put('1', {}, 600); cache.put('2', {}, 600);
  assert.equal(cache.get('1'), null);
  assert.ok(cache.bytes <= cache.maxBytes);

  cache.clear();
  responder = () => new Promise(resolve => {release = resolve;});
  const old = cache.fetch('/old');
  const obsoleteSignal=lastSignal;
  cache.setIdentity('new-job-state-layer');
  assert.ok(obsoleteSignal.aborted,'result changes must cancel unfinished downloads');
  release(ok({stale: true}));
  await assert.rejects(old,/Visualization changed/);
  assert.equal(cache.get('/old'), null);

  let fail = true;
  responder = () => fail ? {ok: false, status: 404, text: async () => '{"detail":"not ready"}'} : ok({ready: true});
  await assert.rejects(cache.fetch('/retry'), /not ready/);
  assert.equal(cache.get('/retry'), null);
  fail = false;
  assert.equal((await cache.fetch('/retry')).ready, true);

  let polls = 0, messages = [];
  responder = () => ++polls < 3 ? ok({detail: 'Preparing CombPlot…'}, 202) : ok({ready: true});
  assert.equal((await cache.fetch('/async', message => messages.push(message))).ready, true);
  assert.equal(polls, 3); assert.equal(messages.length, 2);
  assert.equal(cache.get('/async').ready, true);

  const expandStart = source.indexOf('function expandPlotColumns(');
  const expandEnd = source.indexOf('async function loadVisualizationPanel(', expandStart);
  vm.runInContext(source.slice(expandStart, expandEnd), ctx);
  const encoded = {query: {encoding: 'columns-v1', length: 2,
    columns: {barcode: ['α', 'b'], population: [0, 0], x: [0, 1]}, dictionaries: {population: ['A']}}};
  const expanded = ctx.expandPlotColumns(encoded);
  assert.equal(encoded.query.encoding, 'columns-v1');
  assert.equal(expanded.query[0].barcode, 'α');
  assert.equal(expanded.query[1].population, 'A');
  const barsStart = source.indexOf('function combGeneTraces(');
  const barsEnd = source.indexOf('function renderCombPlotFigure(', barsStart);
  ctx.featureHover = name => name;
  vm.runInContext(source.slice(barsStart, barsEnd), ctx);
  const values = [0, 3, -2, -0, 0, 5];
  const colors = ['a', 'b', 'c', 'd', 'e', 'f'];
  const ids = ['cell0', 'cell1', 'cell2', 'cell3', 'cell4', 'cell5'];
  const traces = ctx.combGeneTraces('G', values, colors, ids, 'y2', 'rna');
  const restored = new Map();
  traces.forEach(trace => trace.x.forEach((x, i) => {
    assert.equal(trace.text[i], ids[x]);
    assert.equal(trace.yaxis, 'y2');
    assert.equal(trace.y[i], values[x] || 0);
    restored.set(x, trace.y[i]);
    if (trace.type === 'bar') assert.equal(trace.marker.color[i], colors[x]);
  }));
  assert.equal(restored.size, values.length);
  assert.equal(traces[0].x.length, 3);
  assert.equal(traces[1].type, 'scatter'); // SVG, not a raster export
  assert.equal(traces[1].mode, 'lines');  // no per-zero SVG marker nodes
  console.log('Explore cache checks passed: deduplication, LRU, byte budget, TTL, context changes, errors, 202 polling and immutable expansion.');
})().catch(error => {console.error(error); process.exitCode = 1;});
