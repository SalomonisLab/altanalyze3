// Exercise the production menu gates for successful and failed communication jobs.
const assert = require('node:assert/strict'), fs = require('node:fs'), vm = require('node:vm'), path = require('node:path');
const source = fs.readFileSync(path.join(__dirname, '../../altanalyze3/components/cellHarmony/webapp/static/app.js'), 'utf8');
const ctx = vm.createContext({
  currentFastCommAnalysis: null, currentMarkerAnalysisByModality: {}, currentMarkerAnalysis: null,
  currentAnalysisMode: 'supervised', clusterAssociationsAvailable: false, currentBiomarkersAvailable: false,
  window: {}, crossPathwayContexts: {}, currentModality: 'rna',
  panelModality: () => ctx.currentModality, modalityDefinition: () => ({}), availableModalities: () => [],
  panelMarkerAnalysis: () => null, markerNetworkPopulations: () => [],
  BASE_VISUALIZATION_MODES: [{value: 'cluster', label: 'UMAP cell types'}], getResultsJobId: () => 'fixture',
});
vm.runInContext(source.slice(source.indexOf('function fastCommAvailable('), source.indexOf('function updateExpressionModeOptions(')), ctx);
const offered = () => ctx.availableVisualizationModes('viz1').some(mode => mode.value === 'fastcomm_network');
assert.equal(offered(), false);
ctx.currentFastCommAnalysis = {enabled: true, status: 'completed', populations: ['sender', 'receiver']};
for (const modality of ['rna', 'adt', 'lipids', 'grn_tf']) {
  ctx.currentModality = modality;
  assert.equal(offered(), true, `${modality} must not hide a completed communication result`);
}
ctx.currentFastCommAnalysis.per_sample = {status: 'failed', message: 'explicit split failure'};
assert.equal(ctx.fastCommAvailable(), true, 'a split failure must not hide completed global scores');
ctx.currentFastCommAnalysis = {enabled: false, status: 'failed', message: 'explicit global failure'};
assert.equal(ctx.fastCommAvailable(), false, 'failed scores are not presented as usable');
assert.equal(offered(), true, 'failed communication remains discoverable with its diagnostic');
ctx.currentFastCommAnalysis = {enabled: false};
assert.equal(offered(), false, 'datasets that never ran communication do not claim a result');
console.log('Communication menu regressions passed across RNA, ADT, lipids, TF activity and failure states.');
