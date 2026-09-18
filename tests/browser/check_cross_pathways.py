"""Exercise real Chat counts, AND filters, chart and context-preserving navigation."""
import json
import sys
from playwright.sync_api import sync_playwright

URL=sys.argv[1] if len(sys.argv)>1 else 'http://127.0.0.1:8000/?job_id=1b0faa928a9848fbb3cd225ca26a0e42'
STATE=sys.argv[2] if len(sys.argv)>2 else 'HSC-1'
with sync_playwright() as p:
    browser=p.chromium.launch(channel='chrome',headless=True)
    page=browser.new_page(viewport={'width':1700,'height':1100})
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto(URL,wait_until='domcontentloaded')
    page.wait_for_function('areExploreResultsReady()',timeout=120000)
    page.click('[data-tab="chat"]')
    page.fill('#chat-question',f'What pathways have the best cross-modality {STATE} representation?')
    page.click('#chat-send')
    page.wait_for_function('chatLastResult?.intent === "cross_modality_pathways"',timeout=60000)
    result=page.evaluate('chatLastResult')
    assert result['status']=='ok',result
    assert result['plot']['source']=='marker'
    assert page.locator('#chat-table tbody tr').count()>0
    page.check('input[data-modality="rna"]')
    if STATE=='HSC-1':page.check('input[data-modality="lipid"]')
    filtered=page.evaluate('filteredChatRows()')
    assert filtered and all(r['rna']>0 and (r['lipid']>0 if STATE=='HSC-1' else True) for r in filtered)
    assert all(r['combined_score']==next(x for x in result['table']['rows'] if x['id']==r['id'])['combined_score'] for r in filtered)
    page.click('#chat-views [data-view="plot"]')
    page.wait_for_function('document.getElementById("chat-plot").layout?.barmode === "stack"')
    chart=page.evaluate('document.getElementById("chat-plot").data.map(t=>({name:t.name,x:t.x,y:t.y}))')
    for mod,trace in zip(result['plot']['modalities'],chart):
        assert trace['x']==[r[mod['id']] for r in filtered[:25]]
    page.screenshot(path=f'/tmp/cross-pathways-chart-{STATE}.png')
    with page.expect_download() as download:
        page.click('#chat-pathways-pdf')
    download.value.save_as(f'/tmp/cross-pathways-chart-{STATE}.pdf')
    bar=page.locator('#chat-plot .barlayer .trace .point path').first
    bar.scroll_into_view_if_needed()
    bounds=bar.bounding_box()
    page.mouse.click(bounds['x']+bounds['width']/2,bounds['y']+bounds['height']/2)
    page.wait_for_function('document.querySelector("#viz1-plot [data-role=figure] svg")',timeout=60000)
    assert page.locator('#viz1-plot [data-setting=id]').input_value()==filtered[0]['id']
    page.click('[data-tab="chat"]')
    # Table selection opens the Explore diagram with the selected pathway/state.
    page.click('#chat-views [data-view="table"]')
    expected=filtered[0]['id']
    page.locator('#chat-table tbody tr').first.locator('button').first.click()
    page.wait_for_function('document.querySelector("#viz1-plot [data-role=figure] svg")',timeout=60000)
    assert page.locator('#viz1-mode').input_value()=='integrated_cross_pathway'
    assert page.locator('#viz1-marker-population').input_value()==STATE
    assert page.locator('#viz1-plot [data-setting=id]').input_value()==expected
    assert page.locator('#viz1-plot [data-role=figure] linearGradient').count()>0
    assert not page.locator('#viz1-modality').is_visible()
    page.screenshot(path=f'/tmp/cross-pathways-diagram-{STATE}.png')
    # Full vector diagram PDF follows the application's shared PDF export path.
    page.evaluate('window.__pdfCheck=null; window.__savePdfOriginal=saveSvgMarkupAsPdf; saveSvgMarkupAsPdf=async (markup,name)=>{window.__pdfCheck={markup,name};}')
    page.evaluate('document.getElementById("viz1-plot")._integratedPdf()')
    pdf=page.evaluate('window.__pdfCheck')
    assert '<linearGradient' in pdf['markup'] and 'Cell-state markers' in pdf['markup']
    page.evaluate('() => {saveSvgMarkupAsPdf=window.__savePdfOriginal;}')
    with page.expect_download() as download:
        page.evaluate('document.getElementById("viz1-plot")._integratedPdf()')
    download.value.save_as(f'/tmp/cross-pathways-diagram-{STATE}.pdf')
    page.click('[data-tab="chat"]')
    contrast_state='MPP-1' if STATE=='HSC-1' else STATE
    page.fill('#chat-question',f'What pathways have the best cross-modality {contrast_state} representation for contrasts?')
    page.click('#chat-send')
    page.wait_for_function('chatLastResult?.plot?.source === "differential"',timeout=120000)
    differential=page.evaluate('chatLastResult')
    assert differential['status']=='ok',differential
    assert differential['plot']['contrast']
    page.locator('#chat-table tbody tr').first.locator('button').first.click()
    page.wait_for_function('document.querySelector("#differential-plot-area [data-role=figure] svg")',timeout=120000)
    assert page.locator('#differential-viz-mode').input_value()=='integrated_cross_pathway'
    assert page.locator('#differential-result-population').input_value()==contrast_state
    assert page.locator('#differential-plot-area [data-setting=id]').input_value()==differential['table']['rows'][0]['id']
    assert page.evaluate('currentDifferentialState.run_id')==differential['plot']['contrast']
    assert not errors,errors
    print(json.dumps({'passed':True,'state':STATE,'markers':len(result['table']['rows']),'differential':len(differential['table']['rows']),'top':{k:v for k,v in result['table']['rows'][0].items() if k!='hits'}}))
    browser.close()
