const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
let replies=[],calls=0;
const cache=new Map();
const request=()=>{calls++;return new Promise(resolve=>replies.push(resolve))};
const ctx=vm.createContext({console,featureAnnotationCache:cache,featureAnnotationRequests:new Map(),explorePayloadCache:{generation:1},getResultsJobId:()=> 'job',apiPath:x=>x,
 fetch:async()=>{const data=await request();return {ok:true,json:async()=>data}},exploreMetadataCache:{fetch:request}});
const a=source.indexOf('async function ensureFeatureAnnotations('),b=source.indexOf('function featureAnnotation(',a);vm.runInContext(source.slice(a,b),ctx);
(async()=>{
 const old=ctx.ensureFeatureAnnotations('job');ctx.explorePayloadCache.generation++;
 const current=ctx.ensureFeatureAnnotations('job');
 assert.equal(calls,2,'a rerun must issue its own annotation request');
 replies[1]({M1:{label:'current source'}});await current;
 replies[0]({M1:{label:'old source'}});await old;
 assert.equal(cache.get('job').M1.label,'current source','old response cannot replace annotations for the new source');
 console.log('Annotation rerun race passed: current request independent and stale response discarded.');
})().catch(e=>{console.error(e);process.exitCode=1});
