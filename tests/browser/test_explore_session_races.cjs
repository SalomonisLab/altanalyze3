// Result readiness and UI ownership across reruns: exercise production bodies.
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
let replies=[];
const ui={textContent:'unchanged',setAttribute(){}};
const ctx=vm.createContext({console,Promise,
  document:{getElementById:()=>ui,querySelector:()=>ui},
  getResultsJobId:()=> 'job',explorePayloadCache:{generation:1},
  exploreResultsReadyJobId:null,exploreResultsPendingJobId:null,exploreResultsReadyPromise:null,
  exploreWarmupJobId:null,exploreWarmupPromise:null,exploreAutoOpenPendingJobId:null,
  currentJobStatus:'completed',currentDifferentialState:null,referenceRerunPending:false,activeExplorerTab:'explore',
  updateWorkflowPanels(){},populateDownloadLinks:async()=>{},loadDisplayFilters:async()=>{},
  loadGeneSuggestions:()=>new Promise(resolve=>replies.push(resolve)),warmExploreResults:async()=>{},
  loadChatExamples:async()=>{},setResultMode(){},syncExplorerWorkspace(){},buildQcCellSummary:()=>'',setExplorerTab(){}
});
const start=source.indexOf('function areExploreResultsReady('),end=source.indexOf('function apiPath(',start);
vm.runInContext(source.slice(start,end),ctx);
(async()=>{
  const old=ctx.ensureExploreResultsReady('job');
  ctx.explorePayloadCache.generation++;ctx.resetExploreResultsReadiness();
  const current=ctx.ensureExploreResultsReady('job');
  const currentPromise=ctx.exploreResultsReadyPromise;
  replies[0]();
  assert.equal(await old,false,'old generation cannot mark a rerun ready');
  assert.equal(ctx.exploreResultsReadyJobId,null);
  assert.equal(ctx.exploreResultsReadyPromise,currentPromise,'old completion cannot clear the current readiness request');
  replies[1]();assert.equal(await current,true);assert.equal(ctx.exploreResultsReadyJobId,'job');
  // A late failing request must not replace a current page with an error.
  ctx.resetExploreResultsReadiness();ctx.loadGeneSuggestions=()=>Promise.reject(new Error('obsolete failure'));
  const failed=ctx.ensureExploreResultsReady('job');
  ctx.explorePayloadCache.generation++;ctx.resetExploreResultsReadiness();
  assert.equal(await failed,false);
  assert.equal(ui.textContent,'unchanged');
  console.log('Session readiness regression passed: obsolete completion/error cannot mutate a newer result.');
})().catch(e=>{console.error(e);process.exitCode=1});
