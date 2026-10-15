const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const code=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
class Node {
 constructor(tag){this.tag=tag;this.children=[];this.attrs={};this.events={};this.open=false;this.isConnected=true;}
 append(n){this.children.push(n)} replaceChildren(){this.children=[]}
 setAttribute(k,v){this.attrs[k]=v}removeAttribute(k){delete this.attrs[k]}
 addEventListener(k,f){this.events[k]=f}fire(k){return this.events[k]?.({target:this})}
}
let calls=[],rendered=[];
const ctx=vm.createContext({URLSearchParams,AbortController,Math,chatRequest:1,getResultsJobId:()=> 'job',
 document:{createElement:tag=>new Node(tag)},summaryPoll:(url,signal,progress)=>new Promise(resolve=>calls.push({url,progress,resolve}))});
vm.runInContext(code.slice(code.indexOf('function summarySection('),code.indexOf('function summaryInterpret(')),ctx);
ctx.summaryPoll=(url,signal,progress)=>new Promise(resolve=>calls.push({url,progress,resolve}));
const digest={summary_identity:'stamp'},host=new Node('host');
for(const params of [{section:'synthesis'},{section:'synthesis',source_id:'rna'},{section:'literature'}]){
 ctx.summaryLazySection(host,params.source_id || params.section,digest,params,(body,result)=>rendered.push(result));
}
assert.equal(calls.length,0,'initial summary does not preload any interpretation');
(async()=>{
 for(let i=0;i<host.children.length;i++){
  const panel=host.children[i];panel.open=true;panel.fire('toggle');assert.equal(calls.length,i+1);
  assert.equal(panel.attrs['aria-busy'],'true');const status=panel.children[1].children[0];assert.equal(status.attrs.role,'status');
  calls[i].progress({detail:'Writing biological synthesis',elapsed_seconds:5});assert.match(status.textContent,/5 s elapsed/);
  calls[i].resolve({narrative:{model:'actual',sections:[{text:'Substantive model prose'}]}});await new Promise(resolve=>setImmediate(resolve));
  assert.equal(panel.attrs['aria-busy'],undefined);panel.open=false;panel.fire('toggle');panel.open=true;panel.fire('toggle');
  assert.equal(calls.length,i+1,'reopening uses the completed result');
 }
 assert.equal(rendered.length,3);assert.match(calls[1].url,/source_id=rna/);assert.match(calls[2].url,/section=literature/);
 assert(code.includes('renderIntegratedBiologyNarrative(body,result.narrative,digest);renderSummaryLiterature'));
 assert(!code.includes('summarySection(host,"Systems-level synthesis",true)'));
 console.log('Integrative first-access loading passed: no preload, per-section queries, spinner, elapsed progress, cached reopening, literature narrative rendering.');
})().catch(e=>{console.error(e);process.exitCode=1});
