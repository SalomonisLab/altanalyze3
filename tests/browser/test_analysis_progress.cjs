// Exercise the real status renderer with unfinished metadata, and delayed submissions.
const assert = require('node:assert/strict');
const fs = require('node:fs'), vm = require('node:vm'), path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/cellHarmony/webapp/static/app.js'), 'utf8');
const elements = new Map();
function element(id) {
  if (!elements.has(id)) elements.set(id, {value:'', textContent:'', style:{}, disabled:false,
    classList:{add(){},remove(){},toggle(){}}, setAttribute(name,value){this[name]=value;}});
  return elements.get(id);
}
const button = element('submit');
let alerts=[], requests=[], polls=[];
const noop=()=>{};
const ctx = vm.createContext({console, window:{}, document:{getElementById:element,
  querySelector:selector=>selector.includes('qc-form button') ? button : null},
  currentJobStatus:'uploaded', previousJobStatus:'', qcSubmissionJobId:'', currentJobSpecies:'',
  currentJobReference:'', currentAnalysisMode:'', currentAnalysisLayers:null, activeAnalysisLayer:'',
  referenceRerunPending:false, currentBiomarkerStates:[], clusterAssociationsAvailable:false,
  currentBiomarkersAvailable:false, currentMarkerAnalysis:null,currentMarkerAnalysisByModality:{},
  currentFastCommAnalysis:null,currentModalitiesState:null,lastDownloadArtifactSignature:'',
  exploreAutoOpenPendingJobId:null, loadedResultsJobId:null,pollTimer:null,VISUALIZATION_PANELS:[],
  explorePayloadCache:{setIdentity:noop}, setSessionUrl:noop, renderAnalysisLayers:noop,
  updateExpressionModeOptions:noop, updateWorkflowPanels:()=>ctx.updateAnalysisRunControl(),
  formatPanelLogTail:(head,tail)=>[...head,...tail].join('\n'), renderQcLiveProgress:noop,
  buildQcCellSummary:data=>data.message, updateResetDataButton:noop, updateDifferentialUi:noop,
  updateReferenceChangeState:noop, resetExploreResultsReadiness:noop, clearDisplayFilters:noop,
  ensureExploreResultsReady:async()=>{}, clearInterval:noop,
  selectedReferenceDiffersFromLoadedJob:()=>false, apiPath:p=>p,
  setResultMode:noop,startStatusPolling:id=>polls.push(id),alert:m=>alerts.push(m),
  fetch:(url,options)=>new Promise(resolve=>requests.push({url,options,resolve}))});
function load(start,end) {vm.runInContext(source.slice(source.indexOf(start),source.indexOf(end,source.indexOf(start))), ctx);}
load('function updateAnalysisRunControl(', 'function updateWorkflowPanels(');
load('function applyJobStatus(', 'function updateAnalysisRunControl(');
load('function initializeBothPanelModes(', 'function renderAnalysisLayers(');
load('async function handleQcSubmit(', 'async function handleDifferentialSubmit(');
element('qc-job-id').value='job'; element('results-job-id').value='job';
const response=(data,status=200)=>({ok:status<400,status,json:async()=>data});
const evt={preventDefault(){},target:Object.fromEntries(['min_genes','min_counts','min_cells','mit_percent','align_cutoff','ambient_correction'].map(k=>[k,{value:'0'}]))};
async function settle(){await new Promise(resolve=>setImmediate(resolve));}
(async()=>{
  // Previously threw at undefined layers.find, before progress/log/QC were rendered.
  for (const mode of ['supervised','unsupervised','both']) {
    ctx.applyJobStatus('job',{status:'processing',progress:44,message:'ICGS3 starting',
      analysis_mode:mode,cell_state_layers:{},icgs3_analysis:{},log_tail:['QC complete']});
    assert.equal(element('job-progress-label').textContent,'44%');
    assert.equal(element('job-log').textContent,'QC complete');
    assert.equal(element('qc-cell-status').textContent,'ICGS3 starting');
    assert.equal(button.textContent,'Running…'); assert.equal(button.disabled,true);
  }
  ctx.applyJobStatus('job',{status:'processing',progress:68,message:'UMAP fitting',
    cluster_key:'unsupervised_state',cell_state_layers:null,icgs3_analysis:{clusters:['c1']}});
  assert.equal(element('job-progress-label').textContent,'68%');
  ctx.applyJobStatus('job',{status:'completed',progress:100,message:'Completed',cell_state_layers:{}});
  assert.equal(element('job-progress-label').textContent,'100%'); assert.equal(button.disabled,false);
  ctx.currentJobStatus='uploaded';
  const first=ctx.handleQcSubmit(evt);
  assert.equal(button.disabled,true); assert.equal(button.textContent,'Starting…');
  await ctx.handleQcSubmit(evt);
  assert.equal(requests.length,1,'double click must not duplicate QC requests');
  requests[0].resolve(response({})); await settle();
  assert.equal(requests.length,2); assert.match(requests[1].url,/\/run$/);
  requests[1].resolve(response({status:'queued'})); await first;
  assert.equal(button.textContent,'Queued…'); assert.equal(button.disabled,true);
  assert.match(element('qc-cell-status').textContent,/Analysis queued/);
  await ctx.handleQcSubmit(evt); assert.equal(requests.length,2);
  assert.equal(alerts.length,0);
  // A different browser tab can start the same job between upload and QC submit.
  ctx.currentJobStatus='uploaded'; requests=[];
  const duplicate=ctx.handleQcSubmit(evt);
  requests[0].resolve(response({detail:'already running'},409)); await settle();
  assert.match(requests[1].url,/\/status$/);
  requests[1].resolve(response({status:'processing',progress:36,message:'QC running',cell_state_layers:{}}));
  await duplicate;
  assert.equal(alerts.length,0); assert.equal(element('job-progress-label').textContent,'36%');
  assert.equal(button.textContent,'Running…'); assert.equal(requests.length,2);
  // Genuine validation failures must remain visible and permit correction.
  ctx.currentJobStatus='uploaded'; requests=[];
  const invalid=ctx.handleQcSubmit(evt);
  requests[0].resolve(response({detail:'Invalid Max K'},422)); await invalid;
  assert.deepEqual(alerts,['Invalid Max K']); assert.equal(button.disabled,false);
  // Reset/new upload while saving QC must not start or modify the old job.
  requests=[]; alerts=[]; const stale=ctx.handleQcSubmit(evt);
  element('qc-job-id').value='other';
  requests[0].resolve(response({})); await stale;
  assert.equal(requests.length,1); assert.equal(button.disabled,false); assert.equal(alerts.length,0);
  console.log('Analysis progress passed: all modes, empty metadata, 44→68→100%, duplicate clicks, active 409, validation and reset.');
})().catch(err=>{console.error(err);process.exitCode=1;});
