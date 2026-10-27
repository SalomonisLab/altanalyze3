const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const source = fs.readFileSync(require('node:path').join(__dirname, '../cellHarmony/webapp/static/app.js'), 'utf8');
const select = {dataset: {}, value: '', disabled:false, options:[], replaceChildren(...options) {this.options=options;this.value=options[0]?.value || ''},add(option){this.options.push(option)}};
const status = {textContent:''};
const counts = {...select, dataset:{}, options:[]};
const layout = {...select, dataset:{}, options:[]};
const mode = {value:'existing'};
const button = {disabled:false,setAttribute(){}};
const job = {value:'first'};
const replies = [];
let updates = 0;
const ctx=vm.createContext({console,Promise,Number,String,
  Option:class {constructor(text,value){this.text=text;this.value=value}},
  jobSessionGeneration:1, document:{querySelector(){return null},getElementById(id){return ({'existing-state-column':select,'existing-annotation-status':status,'qc-job-id':job,'existing-counts-source':counts,'existing-embedding-key':layout,'analysis-mode':mode})[id]}},
  apiPath:path=>path, updateAnalysisRunControl(){updates++},
  fetch:()=>new Promise(resolve=>replies.push(resolve)), parseApiResponse:async response=>response.body,
});
const start=source.indexOf('function existingInspectionText(');
vm.runInContext(source.slice(start,source.indexOf('function syncAnalysisControls()',start)),ctx);
const report={selected:'cell_type',method:'column_metadata',candidates:[{key:'cell_type',categories:2}],warnings:[],embedding_key:'X_tSNE',embedding_candidates:[{key:'X_tSNE',label:'Saved t-SNE'}],count_candidates:[{key:'raw/X',label:'raw.X'}]};
(async()=>{
  const first=ctx.loadExistingAnnotations('first');
  assert.equal(select.disabled,true);assert.equal(select.dataset.loading,'true');assert.equal(select.value,'');
  assert.match(status.textContent,/Inspecting/);assert.equal(updates,1);
  replies.shift()({ok:true,body:report});await first;
  assert.equal(select.value,'cell_type');assert.equal(select.disabled,false);assert.equal(select.dataset.loading,'false');
  const old=ctx.loadExistingAnnotations('first');
  ctx.jobSessionGeneration++;job.value='second';
  const current=ctx.loadExistingAnnotations('second');
  replies.shift()({ok:true,body:report});await old;
  assert.equal(select.disabled,true,'stale first response cannot enable a new session');
  replies.shift()({ok:true,body:{...report,selected:null,warnings:['Review required']}});await current;
  assert.equal(select.value,'','previous session selection must not leak into unresolved annotations');
  assert.match(status.textContent,/could not be identified/);assert.equal(select.disabled,false);
  assert.equal(layout.value,'X_tSNE'); assert.equal(counts.value,'auto');
  const runStart=source.indexOf('function updateAnalysisRunControl()');
  vm.runInContext(source.slice(runStart,source.indexOf('function updateWorkflowPanels(',runStart)),ctx);
  Object.assign(ctx,{qcSubmissionJobId:'',currentJobStatus:'uploaded',currentDifferentialState:null,
    stopSubmissionInFlight:false,uploadSelectionDirty:false,uploadSubmissionInFlight:false});
  ctx.document.querySelector=()=>button;
  select.value='cell_type'; select.dataset.loading='false';
  ctx.updateAnalysisRunControl(); assert.equal(button.disabled,false,'valid t-SNE allows Run');
  layout.value=''; ctx.updateAnalysisRunControl(); assert.equal(button.disabled,true,'no saved coordinates blocks Run');
  layout.value='X_tSNE'; select.value=''; ctx.updateAnalysisRunControl(); assert.equal(button.disabled,true,'missing annotation blocks Run');
  console.log('Use Existing controls passed: loading state, blocked submission state, selection and session fencing.');
})().catch(error=>{console.error(error);process.exitCode=1});
