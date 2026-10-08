// Native saved embeddings and categorical defaults; no model inference.
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
const nodes=new Map();
function node(id){if(!nodes.has(id))nodes.set(id,{value:'',dataset:{},options:[],children:[],appendChild(e){this.children.push(e)},addEventListener(name,fn){this[name]=fn}});return nodes.get(id);}
let loads=[];
const ctx=vm.createContext({console,currentAnalysisMode:'both',VISUALIZATION_PANELS:['viz1','viz2'],document:{getElementById:node,createElement:()=>({addEventListener(name,fn){this[name]=fn}})},panelElementId:(p,s)=>`${p}-${s}`,loadVisualizationPanel:p=>loads.push(p),getDisplayFilterSummary:()=>'',setPanelSummary:(p,text)=>{node(`${p}-filter-summary`).textContent=text;node(`${p}-filter-summary`).children=[]}});
let begin=source.indexOf('function initializeBothPanelModes('),end=source.indexOf('function renderAnalysisLayers(',begin);vm.runInContext(source.slice(begin,end),ctx);
begin=source.indexOf('function setUmapPanelSummary(');end=source.indexOf('function renderPanelUmap(',begin);vm.runInContext(source.slice(begin,end),ctx);
const numeric=[{field:'__umap_1',label:'Unsupervised UMAP 1'},{field:'__umap_2',label:'Unsupervised UMAP 2'},{field:'umap_supervised_x',label:'Supervised UMAP 1'},{field:'umap_supervised_y',label:'Supervised UMAP 2'}];
for(const p of ['viz1','viz2']){
 node(`${p}-mode`).value='expression_umap';node(`${p}-colorby`).options=[{value:''},{value:'supervised_state'}];
 for(const axis of ['xfield','yfield'])node(`${p}-${axis}`).options=numeric.map(n=>({value:n.field,textContent:n.label+' (0 to 1)'}));
}
ctx.initializeBothPanelModes('job');
for(const p of ['viz1','viz2'])ctx.applyBothUmapDefaults(p,'job',{cluster_key:'unsupervised_cluster',numeric_variables:numeric});
assert.equal(node('viz1-mode').value,'cluster');assert.equal(node('viz2-mode').value,'cluster');
assert.equal(node('viz1-colorby').value,'supervised_state');assert.equal(node('viz2-colorby').value,'');
assert.equal(node('viz1-xfield').value,'umap_supervised_x');assert.equal(node('viz2-xfield').value,'__umap_1');
// Polling must not overwrite subsequent user choices.
node('viz1-mode').value='violin';node('viz1-colorby').value='';node('viz1-xfield').value='__umap_1';
ctx.initializeBothPanelModes('job');ctx.applyBothUmapDefaults('viz1','job',{cluster_key:'unsupervised_cluster',numeric_variables:numeric});
assert.equal(node('viz1-mode').value,'violin');assert.equal(node('viz1-xfield').value,'__umap_1');
ctx.setUmapPanelSummary('viz1',{unplaced_population_counts:{Unaligned:36},color_label:'supervised_state',n_dropped_no_coordinate:36,n_cells_selected:3126});
const button=node('viz1-filter-summary').children[0];assert.match(button.textContent,/36 unaligned cells/);button.click();
assert.equal(node('viz1-xfield').value,'__umap_1');assert.equal(node('viz1-yfield').value,'__umap_2');assert.deepEqual(loads,['viz1']);
ctx.setUmapPanelSummary('viz1',{unplaced_population_counts:{},color_label:'supervised_state'});
assert.equal(node('viz1-filter-summary').children.length,0,'no missing cells means no alternate-coordinate prompt');
ctx.currentAnalysisMode='supervised';node('viz2-mode').value='expression_umap';ctx.initializeBothPanelModes('different');assert.equal(node('viz2-mode').value,'expression_umap');
console.log('Both defaults passed: supervised left, unsupervised right, polling preserves choices, retained unaligned cells use saved alternate coordinates.');
const labelStart=source.indexOf('function completedDifferentialLabel('),labelEnd=source.indexOf('function initializeBothPanelModes(',labelStart);
ctx.modalityDefinition=id=>({label:id.toUpperCase()});vm.runInContext(source.slice(labelStart,labelEnd),ctx);
const comparison={modality:'rna',comparison:'Case vs Control',contrast:JSON.stringify({population_col:'supervised_state'})};
assert.match(ctx.completedDifferentialLabel(comparison,[{value:'supervised_state',label:'Supervised cell states'}]),/Supervised cell states/);
comparison.contrast=JSON.stringify({population_col:'unsupervised_cluster'});
assert.match(ctx.completedDifferentialLabel(comparison,[{value:'unsupervised_cluster',label:'Unsupervised clusters'}]),/Unsupervised clusters/);
const dirtyStart=source.indexOf('function markDifferentialConfigDirty('),dirtyEnd=source.indexOf('function updateResetDataButton(',dirtyStart);vm.runInContext(source.slice(dirtyStart,dirtyEnd),ctx);
const saved={status:'completed',run_id:'old-run',progress:100,archive_url:'/old/archive',heatmap_svg_url:'/old.svg',result_populations:['old'],completed_comparisons:[comparison],config:{population_col:'supervised_state',modality:'adt'}};
const dirty=ctx.markDifferentialConfigDirty(saved);
assert.equal(dirty.status,'idle');assert.equal(dirty.run_id,null);assert.equal(dirty.archive_url,null);assert.equal(dirty.heatmap_svg_url,null);assert.equal(dirty.result_populations.length,0);
assert.equal(dirty.completed_comparisons,saved.completed_comparisons,'saved history remains available');assert.equal(saved.run_id,'old-run','dirty form must not mutate a saved result');
assert.equal(dirty.config.modality,'adt');assert.equal(dirty.config.population_col,'supervised_state');
console.log('Saved comparisons distinguish layers, and changed settings remove obsolete result links while retaining history.');
ctx.populateSingleSelect=(select,options,selected)=>{select.options=options;select.value=selected||select.value};
const choices={value:'previous-result',options:[]};
ctx.populateCompletedDifferentialChoices(choices,dirty,[]);
assert.equal(choices.value,'');assert.equal(choices.options[0].label,'Select a saved comparison');
assert.equal(choices.options.length,2,'older saved comparison remains selectable, without representing the changed form');
