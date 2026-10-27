'use strict';
const fs = require('fs');
const vm = require('vm');
const assert = require('assert/strict');
const source = fs.readFileSync(require('path').join(__dirname, '../cellHarmony/webapp/static/app.js'), 'utf8');
function body(name) {
  const start = source.search(new RegExp(`^(?:async )?function ${name}\\(`, 'm'));
  assert(start >= 0, name);
  const tail = source.slice(start);
  const next = tail.slice(1).search(/\n(?:async )?function /);
  return next < 0 ? tail : tail.slice(0, next + 1);
}
async function testStaleResponsesAndPollDeduplication() {
  let releases = [], applied = [], calls = 0;
  const context = vm.createContext({console, Date, Map,
    fetch: () => { calls++; return new Promise(resolve => releases.push(resolve)); },
    apiPath: p=>p, applyJobStatus: (job,data)=>applied.push(job), alert: msg=>{throw Error(msg);},
  });
  vm.runInContext(`let jobSessionGeneration=0, jobStatusRequest=0, appliedJobStatusRequest=0;
    const statusLoads=new Map(); const explorePayloadCache={generation:0};
    let current='old'; function getResultsJobId(){return current;}
    ${body('loadJobState')}\n${body('loadJobStateOnce')}`,context);
  const first = vm.runInContext('loadJobState("old")',context);
  const duplicate = vm.runInContext('loadJobState("old")',context);
  assert.equal(calls,1);
  vm.runInContext('jobSessionGeneration++; current="new";',context);
  releases.shift()({ok:true,json:async()=>({status:'processing'})});
  await Promise.all([first,duplicate]); assert.deepEqual(applied,[]);
  const second=vm.runInContext('loadJobState("new")',context);
  releases.shift()({ok:true,json:async()=>({status:'uploaded'})});
  await second; assert.deepEqual(applied,['new']);
}
async function testUploadRestoresRunControl() {
  const nodes={}; const node=id=>nodes[id]||(nodes[id]={value:'',textContent:'',disabled:false,
    classList:{add(){},remove(){},toggle(){}},setAttribute(){}});
  node('species-select').value='human';node('reference-select').value='reference';
  node('analysis-mode').value='supervised';node('qc-job-id').value='old';
  const row={querySelector: selector=> selector==='.sample-name'?{value:'new sample'}:{files:[{name:'new.h5ad'}]}};
  const context=vm.createContext({console,Set,Array,Map,Number,Date,Boolean,String,
    FormData:class{append(){}},
    document:{getElementById:node,querySelectorAll:()=>[row],querySelector:()=>node('run')},
    apiPath:p=>p, getResultsJobId:()=>node('results-job-id').value,
    stopCurrentAnalysis:async()=>{},
    resetWorkspaceData:()=>{node('qc-job-id').value='';node('results-job-id').value='';},
    setUploadProgress(){},setSessionUrl(){},resetExploreResultsReadiness(){},setResultMode(){},updateResetDataButton(){},
    alert: msg=>{throw Error(msg);},
    uploadFormDataWithProgress:async()=>({ok:true,data:{job_id:'new'}}),
    loadJobState:async()=>vm.runInContext('currentJobStatus="uploaded"; updateAnalysisRunControl()',context),
  });
  vm.runInContext(`let uploadSubmissionInFlight=false, uploadSelectionDirty=true, stopSubmissionInFlight=false;
    let jobSessionGeneration=0, qcSubmissionJobId='', currentJobStatus='', currentDifferentialState=null;
    let loadedResultsJobId=null;
    ${body('existingUploadError')}\n${body('handleJobSubmit')}\n${body('updateAnalysisRunControl')}`,context);
  await vm.runInContext('handleJobSubmit({preventDefault(){}})',context);
  assert.equal(node('qc-job-id').value,'new');
  assert.equal(node('run').disabled,false,'Run must become enabled after upload finally');
  assert.equal(node('upload-submit-btn').disabled,false);
}
async function testJobScopedBrowserLeasesAndReloadIdentity() {
  const storage=new Map([['scalable_browser_client_id','copied-from-another-tab']]);
  const requests=[],beacons=[]; let navigation='navigate';
  const context=vm.createContext({console,Date,Math,JSON,Blob,
    window:{crypto:{randomUUID:()=> 'fresh-independent-tab'}},
    performance:{getEntriesByType:()=>[{type:navigation}]},
    sessionStorage:{getItem:key=>storage.get(key),setItem:(key,value)=>storage.set(key,value)},
    navigator:{sendBeacon:(url,body)=>{beacons.push([url,body]);return true;}},
    setInterval:()=>1, clearInterval(){}, apiPath:p=>p,
    fetch:async(url,options)=>{requests.push([url,JSON.parse(options.body)]);return {ok:true};},
  });
  vm.runInContext(`let browserPresenceJobId='',browserPresenceTimer=null,browserPresenceSequence=0;
    ${body('initializeBrowserPresenceClientId')}\n${body('nextBrowserPresenceSequence')}
    const browserPresenceClientId=initializeBrowserPresenceClientId();
    ${body('syncBrowserPresence')}\n${body('releaseBrowserPresence')}`,context);
  assert.equal(storage.get('scalable_browser_client_id'),'fresh-independent-tab');
  navigation='reload';
  assert.equal(vm.runInContext('initializeBrowserPresenceClientId()',context),'fresh-independent-tab');
  await vm.runInContext('syncBrowserPresence("")',context);
  assert.equal(requests.length,0,'A new visitor without a job ID must not affect any job');
  await vm.runInContext('syncBrowserPresence("job-a")',context);
  await vm.runInContext('syncBrowserPresence("job-a")',context);
  assert.equal(requests.length,1,'Polling status must not spawn duplicate heartbeat timers');
  await vm.runInContext('syncBrowserPresence("job-b")',context);
  assert.equal(beacons[0][0],'/jobs/job-a/browser-presence');
  vm.runInContext('releaseBrowserPresence()',context);
  assert.equal(beacons[1][0],'/jobs/job-b/browser-presence');
  const leaving=JSON.parse(await beacons[1][1].text());
  assert.equal(leaving.state,'leaving');
  assert(leaving.sequence>requests[1][1].sequence);
}
(async()=>{
 await testStaleResponsesAndPollDeduplication(); await testUploadRestoresRunControl();
 await testJobScopedBrowserLeasesAndReloadIdentity();
 console.log('Session response fencing, poll deduplication and post-upload controls passed.');
})().catch(err=>{console.error(err);process.exit(1);});
