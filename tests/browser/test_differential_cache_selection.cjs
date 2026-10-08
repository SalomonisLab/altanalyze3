const assert = require('node:assert/strict'), fs = require('node:fs'), vm = require('node:vm'), path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/cellHarmony/webapp/static/app.js'), 'utf8');
const cacheCode = source.slice(source.indexOf('class ExplorePayloadCache'), source.indexOf('const explorePayloadCache'));
const networkCode = source.slice(source.indexOf('async function differentialNetworkAdjacency('), source.indexOf('// Returns null when no gene is typed', source.indexOf('async function differentialNetworkAdjacency(')));
let reads = 0, now = 0, reply = null, fail = true;
const fullNetwork = {elements: Array.from({length:1000}, (_, i) => ({data:{source:'A', target:'B'+i}}))};
const response = (data, status=200) => ({ok: status<400, status, text: async()=>JSON.stringify(data)});
const ctx = vm.createContext({console, AbortController, Promise, performance: {now:()=>now},
  explorePayloadCache: {generation:1}, currentDifferentialState: {run_id:'run1'}, getResultsJobId:()=> 'job',
  differentialFeatureModality:()=> 'rna', apiPath:x=>x, setTimeout:callback=>callback(),
  fetch: async()=> {reads++; return reply ? reply() : (fail ? response({detail:'temporary failure'},503) : response(fullNetwork));}
});
vm.runInContext(cacheCode + '\nconst differentialPayloadCache = new ExplorePayloadCache(32*1024*1024,8); const differentialNetworkAdjacencyCache = new ExplorePayloadCache(8*1024*1024,8); const differentialNetworkRequests = new Map();',ctx);
const jsonStart=source.indexOf('async function fetchDifferentialJson(');
vm.runInContext(source.slice(jsonStart,source.indexOf('// The two bar colours',jsonStart))+networkCode,ctx);
const tick=()=>new Promise(setImmediate);
(async()=>{
  const failed=await ctx.differentialNetworkAdjacency('job','c1'); assert.equal(failed.available,false);
  fail=false;
  const recovered=await ctx.differentialNetworkAdjacency('job','c1');
  assert.equal(recovered.adjacency.get('A').size,1000,'every edge is retained'); assert.equal(reads,2,'failure must not poison the cache');
  assert.equal(await ctx.differentialNetworkAdjacency('job','c1'),recovered);assert.equal(reads,2);
  for(let i=0;i<12;i++)await ctx.differentialNetworkAdjacency('job','c'+(i+2));
  assert.ok(vm.runInContext('differentialNetworkAdjacencyCache.entries.size<=8 && differentialNetworkAdjacencyCache.bytes<=8*1024*1024',ctx));
  assert.ok(vm.runInContext('differentialPayloadCache.entries.size<=8 && differentialPayloadCache.bytes<=32*1024*1024',ctx));
  now=300001;assert.equal(vm.runInContext('differentialNetworkAdjacencyCache.get(JSON.stringify(["job","run1","rna","c13"]))',ctx),null);
  let release;
  reply=()=>new Promise(resolve=>{release=resolve});
  const old=ctx.differentialNetworkAdjacency('job','pending'); const duplicate=ctx.differentialNetworkAdjacency('job','pending'); await tick();
  const oldReads=reads;
  // Execute the cache-invalidation portion of the production filter reset.
  const clearStart=source.indexOf('function clearDifferentialGeneFilter('), clearEnd=source.indexOf('  const input =',clearStart);
  vm.runInContext(source.slice(clearStart,clearEnd)+'\n}',ctx);
  ctx.clearDifferentialGeneFilter();ctx.explorePayloadCache.generation++;
  release(response(fullNetwork)); await Promise.all([old,duplicate]);
  assert.equal(reads,oldReads);assert.equal(vm.runInContext('differentialNetworkAdjacencyCache.entries.size',ctx),0,'obsolete network cannot repopulate the cleared adjacency cache');
  // A new comparison invalidates successful raw payloads as well as in-flight reads.
  vm.runInContext('differentialPayloadCache.setIdentity("run2")',ctx);ctx.currentDifferentialState={run_id:'run2'};reply=null;
  await ctx.differentialNetworkAdjacency('job','c1'); assert.equal(reads,oldReads+1);
  // Exercise the production state-update prefix rather than assuming callers
  // invalidate the cache when a comparison changes.
  ctx.differentialVisualizationRequest=0;ctx.differentialDetailRequest=0;ctx.currentDifferentialGene='';
  const identityStart=source.indexOf('function updateDifferentialUi(');
  vm.runInContext(source.slice(identityStart,source.indexOf('  const panel =',identityStart))+'\n}',ctx);
  ctx.updateDifferentialUi({run_id:'run2',status:'completed',config:{modality:'rna'}});
  const cacheGeneration=vm.runInContext('differentialPayloadCache.generation',ctx);
  ctx.updateDifferentialUi({run_id:'run2',status:'completed',config:{modality:'rna'}});
  assert.equal(vm.runInContext('differentialPayloadCache.generation',ctx),cacheGeneration,'unchanged status polls reuse payloads');
  ctx.updateDifferentialUi({run_id:'run3',status:'completed',config:{modality:'rna'}});
  assert.equal(vm.runInContext('differentialPayloadCache.generation',ctx),cacheGeneration+1);
  assert.equal(vm.runInContext('differentialPayloadCache.entries.size',ctx),0);

  let postReplies=[],writes=[],applied=[],modeChanges=[];
  const job={value:'job'};
  const selection=vm.createContext({console,Promise,explorePayloadCache:{generation:1},getResultsJobId:()=>job.value,
    savedComparisonSelectionRequest:0,savedComparisonWritePromise:Promise.resolve(),differentialVisualizationRequest:0,differentialDetailRequest:0,
    document:{getElementById:()=>({textContent:''})},apiPath:x=>x,
    fetch:url=>{writes.push(url);return new Promise(resolve=>postReplies.push(resolve));},
    updateDifferentialUi:data=>applied.push(data.run_id),setResultMode:mode=>modeChanges.push(mode)
  });
  const selectionStart=source.indexOf('async function selectCompletedDifferential(');
  vm.runInContext(source.slice(selectionStart,source.indexOf('function markDifferentialConfigDirty(',selectionStart)),selection);
  const first=selection.selectCompletedDifferential('first');await tick();
  const middle=selection.selectCompletedDifferential('middle');const latest=selection.selectCompletedDifferential('latest');
  assert.equal(writes.length,1);postReplies[0]({ok:true,json:async()=>({run_id:'first'})});await tick();
  assert.equal(writes.length,2);assert.ok(writes[1].endsWith('contrast=latest'));
  postReplies[1]({ok:true,json:async()=>({run_id:'latest'})});
  assert.deepEqual(await Promise.all([first,middle,latest]),[false,false,true]);assert.deepEqual(applied,['latest']);
  const stale=selection.selectCompletedDifferential('stale');await tick();selection.explorePayloadCache.generation++;
  postReplies[2]({ok:true,json:async()=>({run_id:'stale'})});assert.equal(await stale,false);assert.deepEqual(applied,['latest']);
  const previousJob=selection.selectCompletedDifferential('oldjob');await tick();job.value='other';
  postReplies[3]({ok:false,json:async()=>({detail:'obsolete error'})});assert.equal(await previousJob,false);
  assert.deepEqual(modeChanges,['differential']);
  console.log('Differential cache/selection checks passed: full edges, retry, bounded caches, deduplication, invalidation and ordered latest selection.');
})().catch(error=>{console.error(error);process.exitCode=1});
