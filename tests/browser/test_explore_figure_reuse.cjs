const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
const children=[],timers=new Map();let timer=0,purges=0,resizes=0,networkDisposed=0,integratedDisposed=0;
function node(id=''){return {id,clientWidth:500,className:'plot-area',svg:true,classList:{hidden:false,add(){this.hidden=true},remove(){this.hidden=false}},querySelector(){return this.svg?{}:null},after(n){children.splice(children.indexOf(this)+1,0,n)},remove(){children.splice(children.indexOf(this),1)}};}
const doc={getElementById:id=>children.find(n=>n.id===id),createElement:()=>node()};
children.push(node('viz1-plot'),node('viz2-plot'));
const ctx=vm.createContext({document:doc,retainedPanelFigures:new Map(),retainedFigureSerial:0,expressionCyByPanel:{},VISUALIZATION_PANELS:['viz1','viz2'],panelPlotId:p=>p+'-plot',
 clearTimeout:id=>timers.delete(id),setTimeout:cb=>{timers.set(++timer,cb);return timer},
 Plotly:{purge:p=>{if(p){p.svg=false;purges++}},Plots:{resize:()=>resizes++}}});
const a=source.indexOf('function removePanelFigure('),b=source.indexOf('function resetVisualizationSurface(',a);vm.runInContext(source.slice(a,b),ctx);
const violin=doc.getElementById('viz1-plot');violin._reusableFigure={key:'violin',bytes:150000*1024+4*1024*1024};
ctx.parkPanelFigure('viz1');assert.equal(ctx.retainedPanelFigures.size,1);assert.ok(violin.classList.hidden);assert.equal(purges,0);
const dot=doc.getElementById('viz1-plot');dot._reusableFigure={key:'dot',bytes:5000000};
assert.equal(ctx.restorePanelFigure('viz1','violin'),true);assert.equal(doc.getElementById('viz1-plot'),violin);assert.ok(!violin.classList.hidden);assert.ok(violin.svg);assert.equal(resizes,0);
dot.clientWidth=800;assert.equal(ctx.restorePanelFigure('viz1','dot'),true);assert.equal(doc.getElementById('viz1-plot'),dot);assert.ok(dot.svg);assert.equal(resizes,1);
ctx.parkPanelFigure('viz1');let entry=ctx.retainedPanelFigures.get('dot');timers.get(entry.timer)();assert.equal(ctx.retainedPanelFigures.has('dot'),false);
const giant=doc.getElementById('viz1-plot');giant._reusableFigure={key:'million',bytes:1024**3};ctx.parkPanelFigure('viz1');assert.equal(doc.getElementById('viz1-plot'),giant);assert.equal(giant._reusableFigure,null);
ctx.disposePanelFigures();assert.equal(ctx.retainedPanelFigures.size,0);assert.ok(doc.getElementById('viz1-plot'));
// A live network/integrated view must be disposed before a cached SVG replaces it.
const cached=doc.getElementById('viz1-plot'); cached.svg=true;
cached._reusableFigure={key:'cached-after-network',bytes:5000000};ctx.parkPanelFigure('viz1');
const network=doc.getElementById('viz1-plot');
network._integratedDispose=()=>integratedDisposed++;
ctx.expressionCyByPanel.viz1={destroy:()=>networkDisposed++};
assert.ok(ctx.restorePanelFigure('viz1','cached-after-network'));
assert.equal(networkDisposed,1);assert.equal(integratedDisposed,1);
assert.equal(ctx.expressionCyByPanel.viz1,null);
assert.equal(doc.getElementById('viz1-plot'),cached);

console.log('Rendered SVG reuse passed: exact node identity, restoration, budget, expiration and disposal.');
// Entering an integrated view must mount on the replacement host, not the SVG
// node that reset just parked. Exercise the production request branch itself.
(async()=>{
 const oldPlot=doc.getElementById('viz1-plot');oldPlot.svg=true;
 oldPlot._reusableFigure={key:'preserved-svg',bytes:5000000};
 const lookup=doc.getElementById;doc.getElementById=id=>id==='results-job-id'?{value:'job'}:lookup(id);
 Object.assign(ctx,{panelVisualizationRequest:{viz1:0},explorePayloadCache:{generation:1},getPanelSelectValue:()=> 'integrated_network',panelModality:()=> 'rna',
   panelElementId:(p,s)=>p+'-'+s,hideMarkerHeatmapView:()=>{},setPanelSummary:()=>{},integratedOptions:()=>({}),
   resetVisualizationSurface:p=>{ctx.parkPanelFigure(p);ctx.releaseVisualizationResources(p,doc.getElementById(p+'-plot'));},
   ScalableIntegrated:{mount:async host=>{host.svg=false;ctx.mountedHost=host}}});
 const begin=source.indexOf('async function loadVisualizationPanel('),end=source.indexOf('function renderVisualizationMessage(',begin);
 vm.runInContext(source.slice(begin,end),ctx);
 await ctx.loadVisualizationPanel('viz1');
 assert.equal(ctx.mountedHost,doc.getElementById('viz1-plot'),'integrated renderer must use current host after reset');
 assert.ok(oldPlot.svg,'cached SVG must not be overwritten by integrated markup');
 assert.ok(oldPlot.classList.hidden,'cached SVG must remain hidden');
 console.log('Integrated mount regression passed: replacement host used and cached SVG preserved.');
 // Execute the actual result-invalidation hook, including active resources.
 let resetCount=0,readinessReset=0;
 Object.assign(ctx,{loadedGeneSuggestionsSignature:'old',loadedDisplayFiltersJobId:'old',plotVariablesCache:{jobId:'old'},
  featureAnnotationCache:new Map([['old',{}]]),exploreMetadataCache:{clear(){}},disposeMarkerHeatmapViews(){},clearChatOutput(){},
  resetExploreResultsReadiness(){readinessReset++},resetDifferentialResults(){},savedComparisonSelectionRequest:0,
  differentialVisualizationRequest:0,differentialPayloadCache:{clear(){}},currentDifferentialInteraction:null,crossPathwayContexts:{},replicateStates:new Set(['old']),loadedResultsJobId:'old',
  panelPlotData:{viz1:{payload:{old:true}},viz2:{payload:{old:true}}},setPanelSummary(){},
  resetVisualizationSurface(panel){resetCount++;ctx.releaseVisualizationResources(panel,doc.getElementById(panel+'-plot'));}});
 const clearStart=source.indexOf('explorePayloadCache.onClear = () => {'),clearEnd=source.indexOf('const panelVisualizationRequest',clearStart);
 vm.runInContext(source.slice(clearStart,clearEnd),ctx);
 ctx.explorePayloadCache.onClear();
 assert.equal(readinessReset,1);assert.equal(resetCount,2);
 assert.equal(ctx.retainedPanelFigures.size,0);
 assert.equal(ctx.loadedResultsJobId,null);assert.equal(ctx.replicateStates,null);
 assert.equal(ctx.panelPlotData.viz1,null);assert.equal(ctx.panelPlotData.viz2,null);
 console.log('Result invalidation regression passed: retained and active resources, payloads and readiness released.');
})().catch(e=>{console.error(e);process.exitCode=1});
