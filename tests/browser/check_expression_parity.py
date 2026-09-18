"""Live saved-session UI regression; reads results and selects saved comparisons only."""
import json
import sys
from playwright.sync_api import sync_playwright


def check(base, job):
    with sync_playwright() as p:
        browser = p.chromium.launch(channel='chrome', headless=True)
        page = browser.new_page(viewport={'width':1700,'height':1150})
        errors, runs = [], []
        page.on('pageerror', lambda e:errors.append(str(e)))
        page.on('request', lambda r:runs.append(r.url) if r.method=='POST' and r.url.endswith(('/run','/qc','/differential')) else None)
        page.goto(f'{base}/?job_id={job}', wait_until='domcontentloaded')
        page.wait_for_function('areExploreResultsReady() && activeExplorerTab === "explore"', timeout=90000)
        assert not page.locator('#viz1-modality-field').is_visible()
        page.select_option('#viz2-mode', 'marker_network')
        page.select_option('#viz2-marker-population', 'ASDC')
        page.wait_for_function('panelPlotData.viz2?.payload?.population === "ASDC" && panelPlotData.viz2.payload.elements?.length > 0')
        for suffix in ('modality-field','gene-field','geneset-field','filter-stack'):
            assert not page.locator('#viz2-'+suffix).is_visible(), suffix
        coverage = page.evaluate('''() => {
          const nodes = panelPlotData.viz2.payload.elements.filter(e => !e.data.source);
          return {nodes:nodes.length,colored:nodes.filter(e => e.data.log2fc !== null).length};
        }''')
        assert coverage['colored'] == coverage['nodes'] and coverage['nodes'] > 0
        page.screenshot(path='/tmp/scalable-marker-network-fixed.png')
        page.select_option('#viz2-mode','dotplot')
        page.wait_for_function('document.getElementById("viz2-groupby").value === "Hs-MarrowAtlas-L3M"')
        assert page.locator('#viz2-groups option').count() == 86
        page.wait_for_function('document.getElementById("viz2-plot").data?.[0]?.mode === "markers"', timeout=60000)
        size = page.evaluate('''() => {
          const el=document.getElementById("viz2-plot");
          return {visible:el.clientWidth,full:el.scrollWidth,plot:el.layout.width,
                  columns:el.layout.xaxis.tickvals.length,overflow:el.style.overflowX};
        }''')
        assert size['overflow'] == 'auto' and size['full'] > size['visible'], size
        assert size['plot'] >= 250+28*size['columns'], size
        page.screenshot(path='/tmp/scalable-dotplot-fixed.png')
        # A new run is not needed to revisit either completed modality.
        page.click('[data-tab="differential"]')
        initial = page.evaluate('currentDifferentialState.config.modality')
        for modality in ('rna','lipid','rna','adt'):
            page.select_option('#differential-modality', modality)
            if modality == 'lipid':
                page.wait_for_function('currentDifferentialState.status === "idle"')
            else:
                page.wait_for_function('(m) => currentDifferentialState.status === "completed" && currentDifferentialState.config.modality === m', arg=modality, timeout=60000)
                assert page.locator('#differential-results-view').is_visible()
        if initial != 'adt':
            page.select_option('#differential-modality', initial)
        page.click('[data-tab="chat"]')
        page.locator('#chat-question').fill('Which TFs have activity discordant with gene expression across cell states?')
        page.click('#chat-send')
        page.wait_for_function('document.getElementById("chat-correlation-plot")?.data?.[0]?.x?.length > 0', timeout=90000)
        assert page.locator('#chat-plot').is_visible()
        assert not page.locator('#chat-table').is_visible()
        assert page.locator('[aria-label="Expression correlation pair"]').count() == 1
        chat = page.evaluate('({status:chatLastResult.status, rows:filteredChatRows(), units:chatLastResult.provenance.unit})')
        assert chat['status'] == 'ok' and all(r['rho'] < 0 for r in chat['rows'])
        page.screenshot(path='/tmp/scalable-expression-chat.png')
        # Job creation calls this same URL helper, with no page navigation.
        page.evaluate('setSessionUrl("test-session-url")')
        assert 'job_id=test-session-url' in page.url
        page.evaluate('(job) => setSessionUrl(job)',job)
        assert f'job_id={job}' in page.url
        assert not errors, errors
        assert not runs, runs
        print(json.dumps({'passed':True,'marker_network':coverage,'dotplot':size,'chat':chat,'restored_modalities':['rna','adt'],'analysis_reruns':0}))
        browser.close()


if __name__ == '__main__':
    check(sys.argv[1].rstrip('/'),sys.argv[2])
