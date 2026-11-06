// Chat output must belong to the newest question and current dataset.
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
let replies=[],signals=[],answers=[],destroyed=0;
const elements=new Map();
for(const id of ['chat-question','results-job-id','chat-status','chat-answer','chat-table','chat-plot','chat-followups','chat-views','chat-result-controls'])
 elements.set(id,{value:'',innerHTML:'',textContent:'',classList:{add(){}}});
elements.get('chat-question').value='first';elements.get('results-job-id').value='job';
const ctx=vm.createContext({console,AbortController,TextDecoder,Uint8Array,
 document:{getElementById:id=>elements.get(id)},Plotly:{purge(){}},
 explorePayloadCache:{generation:1},chatLastResult:null,chatResultState:null,chatRequest:0,chatPlotRequest:0,chatController:null,
 chatNetworkCy:{destroy(){destroyed++}},releaseVisualizationResources(){},
 chatPendingRequest:null,jobSessionGeneration:1,currentJobStatus:'completed',areExploreResultsReady:()=>true,
 ensureFeatureAnnotations:async()=>{},apiPath:url=>url,getResultsJobId:()=>elements.get('results-job-id').value,
 ensureBiologySummaryButton(){},
 fetch:(_,options)=>{signals.push(options.signal);return new Promise(resolve=>replies.push(resolve))},
 renderChatAnswer:data=>answers.push(data.answer)
});
let start=source.indexOf('function clearChatOutput('),end=source.indexOf('// Chat answers carry labels',start);
const helper=source.indexOf('function releaseChatVisualization(');
if(helper>=0)start=helper;
vm.runInContext(source.slice(start,end),ctx);
const response=answer=>({ok:true,headers:{get:()=> 'application/json'},json:async()=>({answer})});
(async()=>{
 const old=ctx.askChat();await new Promise(setImmediate);
 elements.get('chat-question').value='second';
 const current=ctx.askChat();await new Promise(setImmediate);
 assert.ok(signals[0]?.aborted,'the next question must cancel the older request');
 replies[1](response('current'));await current;
 replies[0](response('obsolete'));await old;
 assert.deepEqual(answers,['current'],'obsolete answer cannot overwrite the current answer');
 assert.equal(destroyed,1,'clearing output releases the previous network');
 const previousDataset=ctx.askChat();await new Promise(setImmediate);
 ctx.explorePayloadCache.generation++;ctx.clearChatOutput();
 replies[2](response('previous dataset'));await previousDataset;
 assert.deepEqual(answers,['current']);assert.equal(ctx.chatLastResult,null);
 // Progress and final JSON may split anywhere across UTF-8/network chunks.
 const streamed=ctx.askChat();await new Promise(setImmediate);
 const raw=new TextEncoder().encode(JSON.stringify({type:'progress',message:'Reading data…'})+'\n'+JSON.stringify({type:'result',result:{answer:'streamed answer'}})+'\n');
 const chunks=[raw.slice(0,57),raw.slice(57,58),raw.slice(58)];let position=0;
 const seen=[];
 replies[3]({ok:true,headers:{get:()=> 'application/x-ndjson'},body:{getReader(){return {async read(){
   await new Promise(setImmediate);seen.push(elements.get('chat-status').textContent);
   return position<chunks.length?{value:chunks[position++],done:false}:{done:true};
 }}}}});
 await streamed;assert.deepEqual(answers,['current','streamed answer']);
 assert.ok(seen.includes('Reading data…'),'streamed progress must be visible before the answer');
 // A same-job cache reset during initial restore must not silently drop Ask.
 let ready;ctx.areExploreResultsReady=()=>false;
 ctx.ensureExploreResultsReady=()=>new Promise(resolve=>{ready=resolve});
 const queued=ctx.askChat();await new Promise(setImmediate);
 assert.match(elements.get('chat-status').textContent,/Loading the dataset/);
 ctx.explorePayloadCache.generation++;ctx.clearChatOutput(true);
 assert.match(elements.get('chat-status').textContent,/Loading the dataset/);
 ready(true);await new Promise(setImmediate);
 replies[4](response('queued answer'));await queued;
 assert.deepEqual(answers,['current','streamed answer','queued answer']);
 // A plot fetch has its own ownership guard even within the Chat workspace.
 const plotStart=source.indexOf('async function drawChatPlot('),plotEnd=source.indexOf('var chatNetworkCy',plotStart);
 vm.runInContext(source.slice(plotStart,plotEnd),ctx);
 const plotReplies=[],drawn=[];
 ctx.explorePayloadCache.fetch=()=>new Promise(resolve=>plotReplies.push(resolve));
 ctx.renderDotPlotFigure=(_,data)=>drawn.push(data.id);ctx.escapeHtml=String;
 ctx.chatLastResult={plot:{kind:'dotplot',genes:['G0']}};
 const oldPlot=ctx.drawChatPlot();
 ctx.chatLastResult={plot:{kind:'dotplot',genes:['G1']}};
 const newPlot=ctx.drawChatPlot();
 plotReplies[1]({id:'new plot'});await newPlot;
 plotReplies[0]({id:'old plot'});await oldPlot;
 assert.deepEqual(drawn,['new plot']);
 console.log('Chat regression passed: question/dataset ownership, cancellation, streamed UTF-8 progress and final answer.');
})().catch(error=>{console.error(error);process.exitCode=1});
