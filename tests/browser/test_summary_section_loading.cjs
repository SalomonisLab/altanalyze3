// Section headings trigger one requested interpretation and show real loading state.
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const code=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
class Node {
 constructor(tag='details'){this.tag=tag;this.children=[];this.attrs={};this.events={};this.open=false;this.isConnected=true;this.value='';}
 append(n){this.children.push(n)}
 replaceChildren(){this.children=[]}
 setAttribute(k,v){this.attrs[k]=v}
 removeAttribute(k){delete this.attrs[k]}
 addEventListener(k,f){this.events[k]=f}
 fire(k){return this.events[k]?.({target:this})}
}
let calls=[],renders=[];
const ctx=vm.createContext({Map,URLSearchParams,AbortController,Math,chatRequest:1,getResultsJobId:()=> 'job',
 document:{createElement:tag=>new Node(tag)},summaryText:(parent,tag,value)=>{const n=new Node(tag);n.textContent=value;parent.append(n);return n;},
 summaryPoll:(url,signal,progress)=>new Promise((resolve,reject)=>{calls.push({url,signal,progress,resolve,reject})}),
 appendBiologicalSectionSynthesis:(body,digest,key)=>{renders.push([digest.section_interpretation,key]);body.append(new Node('result'))}});
const start=code.indexOf('function summarySectionInference('),end=code.indexOf('function appendBiologicalSectionSynthesis(',start);
vm.runInContext(code.slice(start,end),ctx);
const tick=()=>new Promise(resolve=>setImmediate(resolve));
(async()=>{
 const host=new Node(),digest={summary_identity:'stamp',interpretation_model:'fixture',deferred_sections:{pathways:[
  {source_id:'P1',state:'coreg_global',title:'global term'},{source_id:'P2',state:'AM',title:'cell-state term'}]}};
 ctx.summarySectionInference(host,digest,'pathways');
 const block=host.children[0],select=block.children[1].children[0],body=block.children[2];
 assert.equal(calls.length,0,'rendering the unopened summary must not run a model');
 host.open=true;host.fire('toggle');
 assert.equal(calls.length,1,'opening the main heading must directly request its interpretation');
 assert.match(calls[0].url,/section=pathways.*source_id=P2/);
 assert.equal(block.attrs['aria-busy'],'true');assert.equal(body.children[0].attrs.role,'status');
 assert.equal(body.children[0].className,'summary-lazy-loading');
 calls[0].progress({detail:'Model inference running',elapsed_seconds:2});assert.match(body.children[0].textContent,/2 s elapsed/);
 calls[0].resolve({review:{status:'completed',model:'fixture',sections:{}}});await tick();
 assert.equal(renders.length,1);assert.equal(block.attrs['aria-busy'],undefined);
 host.open=false;host.fire('toggle');host.open=true;host.fire('toggle');assert.equal(calls.length,1,'reopening reuses the completed interpretation');
 select.value='P1';select.fire('change');assert.equal(calls.length,2);
 select.value='P2';select.fire('change');assert(calls[1].signal.aborted);assert.equal(block.attrs['aria-busy'],undefined);
 calls[1].resolve({review:{model:'stale'}});await tick();assert.equal(renders.at(-1)[0].model,'fixture','stale requests cannot overwrite the selected result');
 const empty=new Node();ctx.summarySectionInference(empty,{deferred_sections:{}},'grn');
 assert.match(empty.children[0].textContent,/No model request was made/);
 console.log('Summary section loading passed: main-heading trigger, spinner, one-source request, cache reuse, source switching, stale-response guard, unavailable-source message.');
})().catch(error=>{console.error(error);process.exitCode=1});
