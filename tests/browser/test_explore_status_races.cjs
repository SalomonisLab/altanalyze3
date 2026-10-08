// Status requests arriving out of order or after reset must not change results.
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
let replies=[],applied=[],alerts=[];
const ctx=vm.createContext({console,Date,jobStatusRequest:0,appliedJobStatusRequest:0,
  explorePayloadCache:{generation:1},activeJob:'job',getResultsJobId:()=>ctx.activeJob,
  fetch:()=>new Promise(resolve=>replies.push(resolve)),apiPath:url=>url,
  alert:message=>alerts.push(message),applyJobStatus:(job,data)=>applied.push([job,data.status]),
  pollTimer:null,setInterval(){},clearInterval(){}
});
const start=source.indexOf('async function loadJobState('),end=source.indexOf('function setSessionUrl(',start);
vm.runInContext(source.slice(start,end),ctx);
const response=status=>({ok:true,json:async()=>({status})});
(async()=>{
  const older=ctx.pollStatus('job'),newer=ctx.loadJobState('job');
  replies[1](response('completed'));await newer;
  replies[0](response('processing'));await older;
  assert.deepEqual(applied,[['job','completed']],'older poll cannot overwrite newer completed status');
  const previous=ctx.loadJobState('job');ctx.activeJob='other';
  replies[2](response('completed'));await previous;
  assert.equal(applied.length,1,'previous job status cannot alter current result');
  ctx.activeJob='job';const reset=ctx.loadJobState('job');
  ctx.explorePayloadCache.generation++;ctx.activeJob='';
  replies[3]({ok:false,json:async()=>({detail:'obsolete error'})});await reset;
  assert.equal(alerts.length,0,'reset must suppress stale status errors');
  ctx.activeJob='job';
  const beforeFailure=ctx.pollStatus('job'),failure=ctx.loadJobState('job');
  replies[5]({ok:false,json:async()=>({detail:'job no longer available'})});await failure;
  replies[4](response('completed'));await beforeFailure;
  assert.equal(applied.length,1,'an older success cannot resurrect a job after a newer status error');
  assert.deepEqual(alerts,['job no longer available']);
  console.log('Status race regression passed: latest response wins and job/reset boundaries are respected.');
})().catch(e=>{console.error(e);process.exitCode=1});
