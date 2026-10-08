// Failed CDN requests must not poison subsequent vector export attempts.
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const source=fs.readFileSync(path.join(__dirname,'../../altanalyze3/components/cellHarmony/webapp/static/app.js'),'utf8');
const scripts=[];let fail=true,loads=0;
const ctx=vm.createContext({window:{},console,
  jsPdfLoaderPromise:null,svg2PdfLoaderPromise:null,cytoscapeSvgLoaderPromise:null,
  document:{scripts,createElement:()=>({dataset:{},listeners:{},addEventListener(name,cb){this.listeners[name]=cb},remove(){const i=scripts.indexOf(this);if(i>=0)scripts.splice(i,1)}}),
    head:{appendChild(script){scripts.push(script);loads++;queueMicrotask(()=>{
      if(fail){script.listeners.error();return;}
      if(script.src.includes('/jspdf@'))ctx.window.jspdf={jsPDF:{API:{}}};
      if(script.src.includes('/svg2pdf.js@'))ctx.window.jspdf.jsPDF.API.svg=()=>{};
      script.listeners.load();
    });}}}
});
const start=source.indexOf('function loadExternalScript('),end=source.indexOf('function slugifyFilenamePart(',start);
vm.runInContext(source.slice(start,end),ctx);
(async()=>{
  await assert.rejects(ctx.ensureJsPdfLoaded(),/Failed to load/);
  assert.equal(scripts.length,0,'failed script elements must be removed so retry starts a new request');
  fail=false;
  await ctx.ensureJsPdfLoaded();assert.ok(ctx.window.jspdf.jsPDF);
  assert.equal(loads,3,'second export retries after both initial CDN failures');
  fail=true;
  await assert.rejects(ctx.ensureSvg2PdfLoaded(),/Failed to load/);
  fail=false;
  await ctx.ensureSvg2PdfLoaded();assert.equal(typeof ctx.window.jspdf.jsPDF.API.svg,'function');
  let installed=false;
  ctx.window.cytoscape=()=>installed?()=>{}:undefined;
  ctx.window.cytoscapeSvg=()=>{installed=true};
  fail=true;await assert.rejects(ctx.ensureCytoscapeSvgLoaded(),/Failed to load/);
  fail=false;await ctx.ensureCytoscapeSvgLoaded();assert.ok(installed);
  console.log('Export retry regression passed: failed scripts and promises released for all three exporters.');
})().catch(error=>{console.error(error);process.exitCode=1});
