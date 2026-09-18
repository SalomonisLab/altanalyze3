"""Regression: integrated panel survives 1↔2 windows with numeric modality scales."""
import json,sys
from playwright.sync_api import sync_playwright
url,state=sys.argv[1:3]
with sync_playwright() as p:
    browser=p.chromium.launch(channel='chrome',headless=True)
    page=browser.new_page(viewport={'width':1750,'height':1150})
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto(url,wait_until='domcontentloaded')
    page.wait_for_function('areExploreResultsReady()',timeout=180000)
    page.click('[data-tab="chat"]')
    page.fill('#chat-question',f'What pathways have the best cross-modality {state} representation?')
    page.click('#chat-send')
    page.wait_for_function('chatLastResult?.plot?.kind === "cross_pathways"',timeout=90000)
    page.locator('#chat-table tbody tr').first.locator('button').first.click()
    page.wait_for_function('document.querySelector("#viz1-plot .dsc-scroll > svg")',timeout=90000)
    options=page.locator('#viz1-plot [data-setting=id] option').evaluate_all('(opts)=>opts.map(o=>o.value)')
    if 'WP5541' in options:
        page.select_option('#viz1-plot [data-setting=id]','WP5541')
        page.wait_for_function('document.querySelector("#viz1-plot [data-role=why]").textContent.startsWith("BDNF")')
    selected=page.locator('#viz1-plot [data-setting=id]').input_value()
    scales=page.locator('#viz1-plot .modality-colorbar').evaluate_all('(nodes)=>nodes.map(n=>({modality:n.dataset.modality,min:n.dataset.min,max:n.dataset.max,label:n.getAttribute("aria-label")}))')
    assert len(scales)>=5 and any(s['max'] and float(s['max'])>0 for s in scales)
    assert all(s['min']=='0' for s in scales if s['max'])
    assert all('MarkerFinder r' in s['label'] for s in scales if s['max'])
    page.evaluate('window.__pathwaySvg=document.querySelector("#viz1-plot .dsc-scroll > svg")')
    page.select_option('#viz2-mode','integrated_network')
    page.wait_for_function('document.querySelector("#viz2-plot.integrated-view h3")?.textContent === "Regulatory network"')
    second_state=page.locator('#viz2-marker-population').input_value()
    for count in ('1','2','1','2'):
        page.select_option('#viz-window-count',count)
        assert page.evaluate('window.__pathwaySvg === document.querySelector("#viz1-plot .dsc-scroll > svg")')
        assert page.locator('#viz1-plot [data-setting=id]').input_value()==selected
        assert page.locator('#viz1-marker-population').input_value()==state
        assert page.locator('#viz2-mode').input_value()=='integrated_network'
        assert page.locator('#viz2-marker-population').input_value()==second_state
        assert 'Gene' not in page.locator('#viz1-plot').inner_text()
        assert page.locator('#viz1-plot .modality-colorbar').count()==len(scales)
    # Two copies of one pathway must not collide through globally scoped SVG IDs.
    page.select_option('#viz2-mode','integrated_cross_pathway')
    page.wait_for_function('document.querySelector("#viz2-plot .dsc-scroll > svg")',timeout=90000)
    ids=page.locator('#viz1-plot [id], #viz2-plot [id]').evaluate_all('(ns)=>ns.map(n=>n.id)')
    assert len(ids)==len(set(ids))
    page.screenshot(path=f'/tmp/pathway-window-scales-{state}.png')
    with page.expect_download() as download:
        page.evaluate('document.getElementById("viz1-plot")._integratedPdf()')
    download.value.save_as(f'/tmp/pathway-window-scales-{state}.pdf')
    assert not errors,errors
    print(json.dumps({'passed':True,'state':state,'pathway':selected,'scales':scales}))
    browser.close()
