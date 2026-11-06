// Chat must mount the existing renderers with the saved comparison selection.
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
const calls=[],figures=[];
const host={innerHTML:'',textContent:'',classList:{add(){},remove(){}},on(){},replaceChildren(...v){this.children=v}};
const ctx=vm.createContext({console,URL,URLSearchParams,window:{location:{href:'http://localhost/'}},
 document:{getElementById:id=>id==='chat-plot'?host:id==='results-job-id'?{value:'job'}:null,createElement:tag=>({tag,style:{},outerHTML:'<base href="/">'})},
 chatLastResult:null,chatPlotRequest:0,chatController:null,explorePayloadCache:{generation:1},
 getResultsJobId:()=> 'job',releaseChatVisualization(){},apiPath:v=>'/api'+v,withRootPath:v=>'/scalable'+v,
 panelElementId:(key,suffix)=>`${key}-${suffix}`,
 Plotly:{newPlot:(target,traces,layout)=>figures.push({target,traces,layout})},
 integratedOptions:(jobId,kind,cell_state,modality,source)=>({jobId,kind,cell_state,modality,source}),
 ScalableIntegrated:{mount:async(target,spec)=>calls.push({target,spec})},
 renderPanelUmap:(...args)=>calls.push(args),renderPanelExpression:(...args)=>calls.push(args),
 renderDotPlotFigure:(...args)=>calls.push(args),renderCombPlotFigure:(...args)=>calls.push(args),
 expandPlotColumns:v=>v,getPlotDotScale:()=>1,
 featureHover:v=>v,featureDisplayName:v=>v,finiteExtent:v=>[Math.min(...v),Math.max(...v)],
 differentialFeatureModality(){throw new Error('Chat must use its own saved modality')},
 destroyDifferentialNetwork(){throw new Error('Chat must not tear down the Differential UI')},
 fetch:async(url)=>{assert.match(url,/^\/scalable\/jobs\/job\/marker\/heatmap\/viewer\?/);
  assert.match(url,/state_layer=supervised_state/);
  assert.match(url,/cells_per_sample=10/,'Chat must share the native heatmap density default');
  return {ok:true,text:async()=>'<html><head></head><body>Native Morpheus</body></html>'}}
});
for(const [start,end] of [['function panelCellsPerSample(', 'function availableVisualizationModes('],
 ['function renderDifferentialHeatmap(', 'function renderDifferentialGo('],
 ['async function drawChatPlot(', 'var chatNetworkCy']])
 vm.runInContext(source.slice(source.indexOf(start),source.indexOf(end,source.indexOf(start))),ctx);
(async()=>{
 for(const view of ['network','pathway']){
  ctx.chatLastResult={plot:{kind:'existing_view',view,cell_state:'AT1',modality:'lipid',source:'differential',contrast:'run-selected'}};
  await ctx.drawChatPlot();const call=calls.at(-1);
  assert.equal(call.target,host);assert.equal(call.spec.kind,view==='network'?'integrated_network':'integrated_pathway');
  assert.equal(call.spec.contrast,'run-selected');assert.equal(call.spec.cell_state,'AT1');assert.equal(call.spec.externalPdf,false);
 }
 ctx.chatLastResult={plot:{kind:'existing_view',view:'volcano',modality:'rna',payload:{population:'AT1',points:[
  {gene:'G1',log2fc:2,score:3,fdr:.001,pval:.0005,direction:'up'}]}}};
 await ctx.drawChatPlot();assert.equal(figures.at(-1).target,host);
 assert.equal(figures.at(-1).traces[1].x[0],2);assert.equal(figures.at(-1).layout.yaxis.title,'-log10(FDR)');
 ctx.chatLastResult={plot:{kind:'existing_view',view:'differential_heatmap',modality:'rna',payload:{population:'AT1',columns:['AT1','AT2'],
 rows:[{gene:'G1',values:[2,null],direction:'up'}]}}};
 await ctx.drawChatPlot();const heatmap=figures.at(-1);
 assert.equal(heatmap.target,host);assert.equal(heatmap.traces[1].z[0][1],null,'missing data stay missing');
 assert.equal(heatmap.traces[0].z[0][1],1,'native missing-value overlay retained');
 ctx.chatLastResult={plot:{kind:'existing_view',view:'marker_heatmap',modality:'rna',state_layer:'supervised_state'}};
 await ctx.drawChatPlot();assert.equal(host.children[0].tag,'iframe');assert.match(host.children[0].srcdoc,/Native Morpheus/);
 // An obsolete heatmap fetch must not replace a newer question's output.
 let resolve;ctx.fetch=()=>new Promise(r=>{resolve=r});
 const old=ctx.drawChatPlot();ctx.chatPlotRequest++;
 host.children=['new answer'];resolve({ok:true,text:async()=>'<head></head>old heatmap'});await old;
 assert.equal(host.children[0],'new answer');
 console.log('Native Chat views passed: shared pathway/network, volcano/heatmap geometry, missing values, Morpheus and request ownership.');
})().catch(e=>{console.error(e);process.exitCode=1});
