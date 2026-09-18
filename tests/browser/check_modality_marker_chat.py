"""Verify shared multi-modality marker Chat against a live saved session."""
import json
import sys
from playwright.sync_api import sync_playwright


def check(url,state):
    with sync_playwright() as p:
        browser=p.chromium.launch(channel='chrome',headless=True)
        page=browser.new_page(viewport={'width':1700,'height':1150})
        errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(url,wait_until='domcontentloaded')
        page.wait_for_function('areExploreResultsReady()',timeout=120000)
        page.click('[data-tab="chat"]')
        page.fill('#chat-question',f'What is the best modality marker of {state}')
        page.click('#chat-send')
        page.wait_for_function('chatLastResult?.intent === "modality_markers" && document.getElementById("chat-plot")?.data?.[0]?.type === "bar"',timeout=60000)
        result=page.evaluate('({status:chatLastResult.status,best:chatLastResult.best,leaders:chatLastResult.best_by_modality,rows:chatLastResult.table.rows.length})')
        assert result['status']=='ok'
        assert page.locator('#chat-plot').is_visible()
        assert page.evaluate('document.getElementById("chat-plot").data[0].x.length')==len(result['leaders'])
        if state=='HSC-1':
            assert result['best']['feature']=='IFIT2' and result['best']['marker_r']==.5082
            assert len(result['leaders'])==5 and result['rows']>50
            page.fill('#chat-results-search','Named')
            assert page.evaluate('filteredChatRows().every(r=>r.annotation === "Named")')
            assert 'Unknown' not in page.evaluate('document.getElementById("chat-plot").data[0].y.join(" ")')
            page.fill('#chat-results-search','')
        page.click('#chat-views [data-view="table"]')
        page.select_option('#chat-results-limit','0')
        assert page.locator('#chat-table tbody tr').count()==result['rows']
        page.screenshot(path=f'/tmp/modality-markers-{state}.png')
        assert not errors,errors
        print(json.dumps({'passed':True,'state':state,**result}))
        browser.close()

if __name__=='__main__':check(sys.argv[1],sys.argv[2])
