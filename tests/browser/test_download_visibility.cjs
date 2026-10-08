const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
const container={children:[],set innerHTML(v){this.children=[];},appendChild(e){this.children.push(e);}};
const context=vm.createContext({document:{getElementById:()=>container,createElement:tag=>({tag})},
 explorePayloadCache:{generation:1},getResultsJobId:()=> 'job',apiPath:p=>p});
const start=source.indexOf('async function populateDownloadLinks('),end=source.indexOf('function updateDifferentialUi(',start);
vm.runInContext(source.slice(start,end),context);
(async()=>{
 const metadata={status:'completed',model_versions:{lipids:{model_version_id:'fixture'}},artifacts:{
  assignments:'a',model_provenance:'p',unaligned_cells_h5ad:'u',combined_h5ad:'h',
  imputed_lipids_results_zip:'l',imputed_adt_results_zip:'p',cluster_associations:'c',unsupervised_marker_genes_zip:'m'}};
 for(const mode of ['supervised','unsupervised','both']){
  await context.populateDownloadLinks('job',{...metadata,analysis_mode:mode});
  assert.ok(!container.children.some(e=>e.tag==='details'||e.href.endsWith('/model_provenance')));
  assert.equal(container.children.some(e=>e.href.endsWith('/unaligned_cells_h5ad')),mode!=='both');
  const logs=container.children.filter(e=>e.href.endsWith('/log'));
  assert.equal(logs.length,1, 'methods replace the existing log download, with no extra button');
  assert.equal(logs[0].textContent,'Download log & methods ZIP');
  assert.ok(container.children.some(e=>e.textContent==='Download Lipids results ZIP'));
  assert.ok(container.children.some(e=>e.textContent==='Download ADT results ZIP'));
  assert.ok(container.children.some(e=>e.textContent==='Download cluster associations'));
 }
 console.log('Download visibility passed: provenance in log, Both union download, all modality archives retained.');
})().catch(e=>{console.error(e);process.exitCode=1;});
