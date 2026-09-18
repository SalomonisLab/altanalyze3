"""Check that named TF activity questions show actual cell-state ranks in both apps."""
import json,sys
from playwright.sync_api import sync_playwright
url=sys.argv[1]
with sync_playwright() as p:
    browser=p.chromium.launch(channel='chrome',headless=True)
    page=browser.new_page(viewport={'width':1700,'height':1100})
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto(url,wait_until='domcontentloaded')
    page.wait_for_function('areExploreResultsReady()',timeout=180000)
    page.click('[data-tab="chat"]')
    page.fill('#chat-question','Where is SPI1 TF-activity most enriched.....')
    page.click('#chat-send')
    page.wait_for_function('chatLastResult?.by_cell_state && chatLastResult?.table?.rows?.length',timeout=120000)
    r=page.evaluate('chatLastResult')
    assert r['status']=='ok'
    rows=r['table']['rows']
    assert len(rows)>25 and all(row['factor']=='SPI1' and row['cell_state'] for row in rows)
    assert rows==sorted(rows,key=lambda row:(-row['activity'],row['cell_state'],row['factor']))
    assert rows[0]['cell_state'] in r['answer']
    page.click('#chat-views [data-view="plot"]')
    page.wait_for_function('document.getElementById("chat-plot")?.data?.[0]?.type === "bar"')
    assert page.evaluate('document.getElementById("chat-plot").data[0].y[0]')==rows[0]['cell_state']
    page.fill('#chat-results-search',rows[0]['cell_state'])
    page.wait_for_function('document.getElementById("chat-plot").data[0].y.length === filteredChatRows().length')
    assert page.evaluate('document.getElementById("chat-plot").data[0].y')==[row['cell_state'] for row in page.evaluate('filteredChatRows()')]
    page.fill('#chat-results-search','')
    page.select_option('#chat-results-limit','0')
    page.wait_for_function('document.getElementById("chat-plot").data[0].y.length === chatLastResult.table.rows.length')
    page.select_option('#chat-results-direction','asc')
    page.wait_for_function('document.getElementById("chat-plot").data[0].x[0] === filteredChatRows()[0].activity')
    assert page.evaluate('document.getElementById("chat-plot").data[0].x[0]')==min(row['activity'] for row in rows)
    assert not errors,errors
    print(json.dumps({'passed':True,'url':url,'states':len(rows),'leaders':rows[:5]}))
    browser.close()
