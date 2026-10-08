// Render the status API's native ICGS3 steps in both shared and discover interfaces.
const assert=require('node:assert/strict'),fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const base=path.join(__dirname,'../../altanalyze3/components/cellHarmony');
const app=fs.readFileSync(path.join(base,'webapp/static/app.js'),'utf8');
const discover=fs.readFileSync(path.join(base,'scalable_discover/static/discover.js'),'utf8');
const start=app.indexOf('function buildQcCellSummary(');
const following=app.slice(start+1).search(/\nfunction \w+\(/);
const context=vm.createContext({getStatusLogLines:d=>d.log_tail||[],extractQcThresholdState:()=>({afterMito:3017}),
  parseProgressPercent:value=>Number(value),referenceRerunPending:false,areExploreResultsReady:()=>true,
  currentAnalysisMode:'unsupervised',discoverClusteringSummary:()=>''});
vm.runInContext(app.slice(start,start+1+following),context);
const stages=[
 'loading the QC-retained counts','normalizing (log1p CP10K)','protein-coding gene filter',
 'downsampling, variable genes and PageRank','variable guide-gene selection and sNMF rank',
 'sNMF clustering','MarkerFinder on the sNMF clusters','SVM reclassification of every cell',
 'GO-Elite BioMarkers cell-state predictions','fitting UMAP'];
function check(interfaceName){
 for (const [i,description] of stages.entries()) {
   const message=`Unsupervised: ICGS3 step ${i+1} of 10: ${description}`;
   const shown=context.buildQcCellSummary({status:'processing',message,progress:40+i,
     log_tail:['Normalization steps: 100%','[ICGS3] continuing current operation']});
   assert.ok(shown.startsWith(message),`${interfaceName} must show current ICGS3 step ${i+1}`);
   assert.ok(!shown.startsWith('Normalizing expression'),`${interfaceName} cannot show earlier QC`);
 }
 assert.equal(context.buildQcCellSummary({status:'queued',message:'Analysis queued.',
   log_tail:['Normalization steps: 100%']}),'Analysis queued.');
}
check('shared web');
const wrapperStart=discover.indexOf('(function wrapQcCellSummary()');
const wrapperEnd=discover.indexOf('\nconst DISCOVER_DOWNLOAD_LABELS',wrapperStart);
vm.runInContext(discover.slice(wrapperStart,wrapperEnd),context);
check('discover');
console.log('Native ICGS3 steps 1–10 render correctly in shared web and discover, despite earlier normalization logs.');
