// Production DOM lifecycle and metadata loaders, with deterministic small DOM mocks.
const assert=require('node:assert/strict'), fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
const elements=new Map(), timers=new Map(); let timerId=0, requests=[], htmlRequests=[], htmlStatus=200, message;
function element(){return {listeners:{},addEventListener(name,handler){this.listeners[name]=handler},dataset:{},children:[],value:'',attributes:{},get outerHTML(){return `<base href="${this.href}">`},classList:{classes:new Set(),add(x){this.classes.add(x)},remove(x){this.classes.delete(x)},contains(x){return this.classes.has(x)}},appendChild(e){this.children.push(e)},replaceChildren(f){this.children=f.children},setAttribute(k,v){this.attributes[k]=v},after(e){this.sibling=e},remove(){this.removed=true},contentWindow:{}};}
for(const panel of ['viz1','viz2'])for(const suffix of ['plot','feature-suggestions','gene-query'])elements.set(`${panel}-${suffix}`,element());
const ctx=vm.createContext({document:{getElementById:id=>elements.get(id),createElement:()=>element(),createDocumentFragment:()=>element()},window:{location:{origin:'http://local',href:'http://local/'},addEventListener:(name,cb)=>{if(name==='message')message=cb}},
 URL,AbortController,fetch:async url=>{htmlRequests.push(url);return {ok:htmlStatus===200,status:htmlStatus,text:async()=>'<html><head></head><body><div id="morpheus-target"></div></body></html>'}},
 setTimeout:cb=>{timers.set(++timerId,cb);return timerId},clearTimeout:id=>timers.delete(id),console,
 markerHeatmapViews:new Map(),explorePayloadCache:{generation:1},panelPlotId:p=>`${p}-plot`,panelElementId:(p,s)=>`${p}-${s}`,
 getDisplayFilterParams:()=>new URLSearchParams(), panelModality:()=> 'rna',panelCellsPerSample:()=>10,withRootPath:x=>x,setPanelSummary:()=>{},
 resetVisualizationSurface:p=>{ctx.hideMarkerHeatmapView(p); elements.get(`${p}-plot`).classList.remove('hidden');},URLSearchParams,
 VISUALIZATION_PANELS:['viz1','viz2'],loadedGeneSuggestionsSignature:'',getResultsJobId:()=> 'job',ensureFeatureAnnotations:async()=>{},apiPath:x=>x,
 exploreMetadataCache:{fetch:async url=>{requests.push(url);return {genes:Array.from({length:32000},(_,i)=>`G${i}`),keys:['ID']}}},
 featureDisplayName:x=>x,preferredFeatureForModality:(_,features)=>features[0],updatePanelFeatureInput:()=>{}});
let clearStart=source.indexOf('function clearGeneSuggestions('), clearEnd=source.indexOf('function clearDisplayFilters(',clearStart);
vm.runInContext(source.slice(clearStart,clearEnd),ctx);
let a=source.indexOf('function removeMarkerHeatmapView('),b=source.indexOf('// A missing or zero fold',a);
vm.runInContext(source.slice(a,b),ctx);
a=source.indexOf('async function loadGeneSuggestions(');b=source.indexOf('function populateSelectOptions(',a);
vm.runInContext(source.slice(a,b),ctx);
(async()=>{
 await ctx.renderMarkerHeatmapViewer('job','viz1');
 const one=ctx.markerHeatmapViews.get('viz1');
 assert.equal(htmlRequests.length,1,'one HTML GET, with no iframe duplicate request');
 assert.equal(htmlRequests[0],'/jobs/job/marker/heatmap/viewer?modality=rna&cells_per_sample=10');
 assert.match(one.iframe.srcdoc,/<base href="http:\/\/local\/jobs\/job\/marker\/heatmap\/viewer\?modality=rna&cells_per_sample=10">/);
 message({origin:'http://wrong',source:one.iframe.contentWindow,data:{type:'marker-heatmap-dataset',columns:100,rows:20}});
 assert.equal(one.bytes,Infinity);
 message({origin:'http://local',source:one.iframe.contentWindow,data:{type:'marker-heatmap-dataset',columns:100,rows:20}});
 ctx.hideMarkerHeatmapView('viz1'); assert.ok(one.host.classList.contains('hidden'));
 await ctx.renderMarkerHeatmapViewer('job','viz1');
 assert.equal(ctx.markerHeatmapViews.get('viz1'),one);assert.ok(!one.host.classList.contains('hidden')); assert.ok(!one.host.removed);
 assert.equal(htmlRequests.length,1,'reusing the live frame must not refetch HTML');
 ctx.explorePayloadCache.generation++;
 await ctx.renderMarkerHeatmapViewer('job','viz1'); assert.ok(one.host.removed);
 const two=ctx.markerHeatmapViews.get('viz1');
 message({origin:'http://local',source:two.iframe.contentWindow,data:{type:'marker-heatmap-dataset',columns:100000,rows:1000}});
 ctx.hideMarkerHeatmapView('viz1');assert.equal(ctx.markerHeatmapViews.size,0);
 await ctx.renderMarkerHeatmapViewer('job','viz1');
 const three=ctx.markerHeatmapViews.get('viz1');
 message({origin:'http://local',source:three.iframe.contentWindow,data:{type:'marker-heatmap-dataset',columns:100,rows:20}});
 ctx.hideMarkerHeatmapView('viz1'); timers.get(three.timer)(); assert.equal(ctx.markerHeatmapViews.size,0);
 await ctx.renderMarkerHeatmapViewer('job','viz1');
 const failed=ctx.markerHeatmapViews.get('viz1');
 message({origin:'http://local',source:failed.iframe.contentWindow,data:{type:'marker-heatmap-error',message:'dataset request failed'}});
 await ctx.renderMarkerHeatmapViewer('job','viz1');
 assert.ok(failed.host.removed,'failed heatmaps must be rebuilt, not served from the live-view cache');
 assert.notEqual(ctx.markerHeatmapViews.get('viz1'),failed);
 ctx.disposeMarkerHeatmapViews();htmlStatus=404;
 await ctx.renderMarkerHeatmapViewer('job','viz1');
 const badHtml=ctx.markerHeatmapViews.get('viz1');
 assert.equal(badHtml.failed,true,'HTTP error documents must not remain reusable');
 htmlStatus=200;
 await ctx.renderMarkerHeatmapViewer('job','viz1');assert.ok(badHtml.host.removed);
  const recovered=ctx.markerHeatmapViews.get('viz1');
  badHtml.iframe.listeners.error();
  assert.equal(ctx.markerHeatmapViews.get('viz1'),recovered,'late events from disposed frames must not remove replacements');
  assert.ok(!recovered.failed,'late errors must not mark the recovered viewer as failed');
  // A pending HTML response must not recreate an already-disposed frame.
  ctx.disposeMarkerHeatmapViews();
  let lateReply;
  const normalFetch=ctx.fetch;
  ctx.fetch=()=>new Promise(resolve=>{lateReply=resolve});
  const oldLoad=ctx.renderMarkerHeatmapViewer('job','viz1');
  const oldView=ctx.markerHeatmapViews.get('viz1');
  ctx.disposeMarkerHeatmapViews();
  assert.ok(oldView.controller.signal.aborted,'disposing a view cancels its pending HTML request');
  ctx.fetch=normalFetch;
  await ctx.renderMarkerHeatmapViewer('job','viz1');
  const replacement=ctx.markerHeatmapViews.get('viz1');
  lateReply({ok:true,text:async()=>'<html><head></head><body>obsolete</body></html>'});
  await oldLoad;
  assert.equal(ctx.markerHeatmapViews.get('viz1'),replacement);
  assert.equal(oldView.iframe.srcdoc,undefined,'late response cannot initialize disposed iframe');
  await ctx.loadGeneSuggestions('job');
 assert.equal(requests.length,1,'one catalog for two RNA panels');
 assert.equal(elements.get('viz1-feature-suggestions').children.length,32000);
 assert.equal(elements.get('viz2-feature-suggestions').children.length,32000);
 const choices=elements.get('viz1-feature-suggestions').children;
 await ctx.loadGeneSuggestions('job');assert.equal(requests.length,1);assert.equal(elements.get('viz1-feature-suggestions').children,choices);
 ctx.clearGeneSuggestions();
 for(const panel of ['viz1','viz2'])elements.get(`${panel}-feature-suggestions`).children=[];
 await ctx.loadGeneSuggestions('job');
 assert.equal(elements.get('viz1-feature-suggestions').children.length,32000,'cleared list must repopulate in the same generation');
 console.log('Explore metadata and live heatmap lifecycle tests passed (complete 32,000-feature catalogs).');
})().catch(e=>{console.error(e);process.exitCode=1});
