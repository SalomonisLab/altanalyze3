import json
import sys
from playwright.sync_api import sync_playwright
URL=sys.argv[1] if len(sys.argv)>1 else 'http://127.0.0.1:8000/?job_id=1b0faa928a9848fbb3cd225ca26a0e42'
with sync_playwright() as p:
 b=p.chromium.launch(channel='chrome',headless=True)
 page=b.new_page(viewport={'width':1800,'height':1150});errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
 page.goto(URL,wait_until='domcontentloaded');page.wait_for_function('areExploreResultsReady()',timeout=120000)
 page.evaluate('ensureFeatureAnnotations(getResultsJobId())')
 print('annotation',page.evaluate('featureAnnotation("Unknown 0235","metabolite")'))
 page.click('[data-tab="differential"]')
 page.wait_for_function('currentDifferentialState?.config?.modality === "metabolite"',timeout=30000)
 result=page.evaluate('''async () => {
 const base=apiPath(`/jobs/${getResultsJobId()}/differential/interactive`);
 const payload=await (await fetch(base+"/volcano?population=MEP-Eryth-2")).json();
 document.getElementById("differential-viz-mode").value="volcano";
 renderDifferentialVolcano(payload);
 const plot=document.getElementById("differential-plot-area");
 const trace=plot.data.find(t=>t.customdata.some(r=>r[0]==="Unknown 0235"));
 const at=trace.customdata.findIndex(r=>r[0]==="Unknown 0235");
 const detail=await (await fetch(base+"/gene?population=MEP-Eryth-2&gene=Unknown%200235")).json();
 renderDifferentialGeneDetail(detail);
 return {hover:trace.hovertext[at],label:trace.text[at],raw:trace.customdata[at][0],detail:document.getElementById("differential-selected-gene").textContent,points:payload.points.length};
}''');assert '234.09706' in result['hover'] and 'Table22!A464' in result['hover'];assert result['raw']=='Unknown 0235';assert 'm/z' in result['detail'];print('volcano',result)
 page.screenshot(path='/tmp/metabolite-volcano-labels.png')
 r=page.evaluate('''async()=>{ const payload=await(await fetch(apiPath(`/jobs/${getResultsJobId()}/differential/interactive/heatmap?population=MEP-Eryth-2`))).json();document.getElementById("differential-viz-mode").value="heatmap";renderDifferentialHeatmap(payload);const p=document.getElementById("differential-plot-area");return {labels:p.layout.yaxis.ticktext.filter(x=>x.includes("m/z")).slice(0,2),hover:p.data[1].text.flat().find(x=>x.includes("PDC000561")),raw:p.data[1].y.filter(x=>x.startsWith("Unknown")).slice(0,2)};}''');assert r['labels'] and r['hover'];assert 'm/z' not in r['raw'][0];print('heatmap',r)
 page.screenshot(path='/tmp/metabolite-heatmap-labels.png')
 page.click('[data-tab="explore"]')
 r=page.evaluate('''async()=>{const payload=await(await fetch(apiPath(`/jobs/${getResultsJobId()}/dotplot?modality=metabolite`))).json(); renderDotPlotFigure("viz1-plot",payload); const p=document.getElementById("viz1-plot");return {genes:payload.genes,labels:p.layout?.yaxis?.ticktext,hover:p.data?.[0]?.text?.find(x=>x.includes("PDC000561"))};}''');print('dotplot',r)
 assert r.get('hover') and any('m/z' in x for x in r['labels'])
 assert not errors,errors

 r=page.evaluate('''async()=>{
 const base=apiPath(`/jobs/${getResultsJobId()}`);
 const single=await(await fetch(base+"/dotplot?modality=metabolite&genes=Unknown%200235")).json();
 const comb=await(await fetch(base+"/combplot?modality=metabolite&genes=Unknown%200235&unit=donor")).json();
 renderCombPlotFigure("viz1-plot",comb);
 const hover=document.getElementById("viz1-plot").data.flatMap(t=>t.hovertext||[]).find(x=>x.includes("PDC000561"));
 const expr=await(await fetch(base+"/expression?modality=metabolite&gene=Unknown%200235")).json();
 renderPanelExpression("viz2",expr,"violin",1);
 return {single:single.genes,comb:comb.genes,modality:comb.modality,combHover:hover,violinHover:document.getElementById("viz2-plot").data[0].hovertemplate};
}''');assert r['single']==r['comb']==['Unknown 0235'];assert 'PDC000561' in r['combHover'] and 'PDC000561' in r['violinHover'];print('expression',r)
 assert not errors,errors
 b.close()
