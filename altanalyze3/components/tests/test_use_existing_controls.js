const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const source = fs.readFileSync(require('node:path').join(__dirname, '../cellHarmony/webapp/static/app.js'), 'utf8');
const select = {dataset: {}, value: '', disabled:false, options:[], replaceChildren(...options) {this.options=options;this.value=options[0]?.value || ''},add(option){this.options.push(option)}};
const status = {textContent:''};
const job = {value:'first'};
const replies = [];
let updates = 0;
const ctx=vm.createContext({console,Promise,Number,String,
  Option:class {constructor(text,value){this.text=text;this.value=value}},
  jobSessionGeneration:1, document:{getElementById(id){return ({'existing-state-column':select,'existing-annotation-status':status,'qc-job-id':job})[id]}},
  apiPath:path=>path, updateAnalysisRunControl(){updates++},
  fetch:()=>new Promise(resolve=>replies.push(resolve)), parseApiResponse:async response=>response.body,
});
const start=source.indexOf('async function loadExistingAnnotations(');
vm.runInContext(source.slice(start,source.indexOf('function syncAnalysisControls()',start)),ctx);
const report={selected:'cell_type',method:'column_metadata',candidates:[{key:'cell_type',categories:2}],warnings:[]};
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
  console.log('Use Existing controls passed: loading state, blocked submission state, selection and session fencing.');
})().catch(error=>{console.error(error);process.exitCode=1});
