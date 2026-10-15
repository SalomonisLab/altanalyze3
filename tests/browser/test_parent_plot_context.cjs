const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const code=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
let plots=[];
const ctx=vm.createContext({console,Map,Set,Math,Number,Array,Object,String,
 document:{getElementById:()=>({})},updateBaselineFilterSummaries:()=>{},panelPlotId:p=>p,
 setUmapPanelSummary:()=>{},relaxReferencePreviewLabels:()=>[],buildPopulationCentroids:()=>[],currentAnalysisMode:'supervised',
 Plotly:{newPlot:(panel,traces,layout)=>plots.push({traces,layout})},
 finiteExtent:(values)=>[Math.min(...values),Math.max(...values)],
 featureDisplayName:g=>g,measurementLabel:()=> 'Expression',featureHover:g=>g,normalizeModalityId:x=>x,
 resetVisualizationSurface:()=>{},renderVisualizationMessage:()=>{throw Error('Unexpected empty plot')},
 buildReferencePreviewColorMap:()=>{throw Error('Authoritative parent palette must be used')}});
vm.runInContext(code.slice(code.indexOf('function buildStableUmapPopulationOrder('),code.indexOf('function networkEdgeColor(')),ctx);
vm.runInContext(code.slice(code.indexOf('function buildSquareUmapAxes('),code.indexOf('function median(')),ctx);
vm.runInContext(code.slice(code.indexOf('function renderPanelUmap('),code.indexOf('function renderVisualizationPanel(')),ctx);
const full={population_order:['R','A','B','C','D','E'],population_colors:{R:'#111111',A:'#ff0000',B:'#0000ff',C:'#00ff00',D:'#ffff00',E:'#abcdef'},
 coordinate_bounds:{x:[-40,80],y:[-100,40]},reference:[{barcode:'r',population:'R',x:-40,y:40}],
 query:[{barcode:'a',population:'A',sample:'one',x:0,y:0},{barcode:'b',population:'B',sample:'one',x:1,y:2},{barcode:'c',population:'C',sample:'two',x:60,y:-80}],sample_field:'Library'};
const subset={...full,query:full.query.slice(1,2)};
const run=(payload,mode)=>{ctx.renderPanelUmap('panel',payload,mode,1);return plots.at(-1)};
const a=run(full,'cluster'),b=run(subset,'cluster');
assert.equal(a.traces[1].marker.color[1],b.traces[1].marker.color[0]);
assert.deepEqual(a.layout.xaxis.range,b.layout.xaxis.range);assert.deepEqual(a.layout.yaxis.range,b.layout.yaxis.range);
assert.equal(b.layout.xaxis.autorange,false);
assert(b.layout.xaxis.range[0]<-40 && b.layout.xaxis.range[1]>80,'parent reference and projected cells define the range');
for(const mode of ['cluster','relative']){
 const f=run({...full,axes_source:'obs',x_label:'X',y_label:'Y'},mode),s=run({...subset,axes_source:'obs',x_label:'X',y_label:'Y'},mode);
 assert.deepEqual(f.layout.xaxis.range,s.layout.xaxis.range);assert.deepEqual(f.layout.yaxis.range,s.layout.yaxis.range);
}
const frequency=run(full,'frequency'),filtered=run(subset,'frequency');
const barColors=p=>new Map(p.traces[0].y.map((name,i)=>[name,p.traces[0].marker.color[i]]));
assert.equal(barColors(frequency).get('B'),barColors(filtered).get('B'));
assert.equal(barColors(filtered).get('B'),b.traces[1].marker.color[0],'frequency and UMAP use the same population color');
for(const points of [full.query,subset.query]){
 ctx.renderPanelExpression('panel',{...full,gene:'G',modality:'rna',source:'query',umap:points.map(p=>({...p,value:1})),global_min:0,global_max:2},'expression_umap',1);
 assert.deepEqual(plots.at(-1).layout.xaxis.range,a.layout.xaxis.range);
 assert.deepEqual(plots.at(-1).layout.yaxis.range,a.layout.yaxis.range);
}
console.log('Parent plot context passed: stable categorical/frequency colors, rare-state filtering, broad/metadata/expression UMAP parent ranges.');
