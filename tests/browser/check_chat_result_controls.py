"""Exercise all retained Chat pairs and table/plot controls on a live session."""
import json
import sys
from playwright.sync_api import sync_playwright


def check(url):
    with sync_playwright() as p:
        browser=p.chromium.launch(channel='chrome',headless=True)
        page=browser.new_page(viewport={'width':1700,'height':1100})
        errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
        page.goto(url,wait_until='domcontentloaded')
        page.wait_for_function('areExploreResultsReady()',timeout=120000)
        page.click('[data-tab="chat"]')
        page.fill('#chat-question','Correlate TF activity with matching gene expression across cell states')
        page.click('#chat-send')
        page.wait_for_function('chatLastResult?.plot?.kind === "cross_modal" && document.getElementById("chat-correlation-plot")?.data?.length',timeout=120000)
        total=page.evaluate('chatLastResult.table.rows.length')
        assert total > 50,total
        assert page.locator('[aria-label="Expression correlation pair"] option').count()==total
        page.click('#chat-views [data-view="table"]')
        assert page.locator('#chat-table tbody tr').count()==50
        page.click('#chat-results-next')
        assert page.evaluate('chatResultState.page')==1
        page.select_option('#chat-results-limit','0')
        assert page.locator('#chat-table tbody tr').count()==total
        page.locator('#chat-table th button').filter(has_text='rho').first.click()
        assert page.evaluate('filteredChatRows().every((r,i,a)=>!i || a[i-1].rho <= r.rho)')
        page.locator('#chat-table th button').filter(has_text='rho').first.click()
        assert page.evaluate('filteredChatRows().every((r,i,a)=>!i || a[i-1].rho >= r.rho)')
        page.select_option('#chat-results-sort','gene')
        page.select_option('#chat-results-direction','asc')
        assert page.locator('#chat-table th[aria-sort="ascending"]').inner_text().startswith('gene')
        page.select_option('#chat-results-relationship','near_zero')
        page.fill('#chat-results-threshold','0.2')
        weak=page.evaluate('filteredChatRows()')
        assert weak and all(abs(r['rho']) <= .2 for r in weak)
        assert [abs(r['rho']) for r in weak] == sorted(abs(r['rho']) for r in weak)
        page.click('#chat-views [data-view="plot"]')
        assert page.locator('[aria-label="Expression correlation pair"] option').count()==len(weak)
        page.select_option('#chat-results-relationship','all')
        page.fill('#chat-results-min','-0.1')
        page.fill('#chat-results-max','0.1')
        assert page.evaluate('filteredChatRows().every(r=>r.rho >= -.1 && r.rho <= .1)')
        page.fill('#chat-results-min','');page.fill('#chat-results-max','')
        # A pair formerly outside the backend's first 50 is now searchable and plottable.
        target=page.evaluate('chatLastResult.table.rows.at(-1)')
        page.fill('#chat-results-search',target['feature'])
        page.click('#chat-views [data-view="table"]')
        assert page.locator('#chat-table tbody tr').count()>=1
        page.get_by_role('button',name=f"Plot {target['feature']} / {target['gene']}",exact=True).click()
        assert page.locator('#chat-plot').is_visible()
        assert page.evaluate('chatResultState.selectedPair')==json.dumps([target['feature'],target['gene']],separators=(',',':'))
        page.fill('#chat-results-search','NO_MATCH_THIS_FEATURE')
        assert 'No pairs match' in page.locator('#chat-plot').inner_text()
        assert page.locator('#chat-result-controls').is_visible()
        page.fill('#chat-question','Which TFs have activity uncorrelated with gene expression across cell states?')
        page.click('#chat-send')
        page.wait_for_function('chatLastResult?.result_controls?.relationship === "near_zero"',timeout=60000)
        assert page.locator('#chat-results-search').input_value()==''
        assert page.locator('#chat-results-relationship').input_value()=='near_zero'
        assert page.evaluate('chatLastResult.table.rows.length')==total
        page.click('#chat-views [data-view="table"]')
        page.screenshot(path='/tmp/chat-result-controls-'+('viewer' if '8080' in url else 'upload')+'.png')
        assert not errors,errors
        print(json.dumps({'passed':True,'url':url,'all_pairs':total,'near_zero_pairs':len(weak),'pair_beyond_first_50':target['feature'],'sortable':True,'pagination':True,'filtered_plot':True}))
        browser.close()

if __name__=='__main__':check(sys.argv[1])
