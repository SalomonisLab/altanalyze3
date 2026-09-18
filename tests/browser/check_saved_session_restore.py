"""Live browser regression: python check_saved_session_restore.py BASE_URL JOB_ID.

Requires Playwright and Chrome. Reads an existing completed job; never runs it.
"""
import json
import sys
from pathlib import Path
from playwright.sync_api import sync_playwright


def check_saved_session(base_url, job_id):
    with sync_playwright() as p:
        browser = p.chromium.launch(channel='chrome', headless=True)
        page = browser.new_page(viewport={'width': 1700, 'height': 1150})
        errors, requests, held = [], [], []
        page.on('pageerror', lambda e: errors.append(str(e)))
        page.on('dialog', lambda d: errors.append(d.message) or d.dismiss())
        page.on('request', lambda r: requests.append((r.method, r.url)))
        pause = True

        def delay_genes(route):
            if pause:
                held.append(route)
            else:
                route.continue_()

        page.route('**/genes?*', delay_genes)
        page.goto(f'{base_url}/?job_id={job_id}', wait_until='domcontentloaded')
        tab = page.locator('[data-tab="explore"]')
        page.wait_for_function('currentJobStatus === "completed"', timeout=60000)
        assert tab.is_visible() and tab.is_disabled()
        assert 'loading' in tab.inner_text()
        assert tab.get_attribute('aria-busy') == 'true'
        assert not [url for _, url in requests if '/reference-preview?' in url]
        pause = False
        for route in held:
            route.continue_()
        page.wait_for_function('areExploreResultsReady() && activeExplorerTab === "explore"', timeout=60000)

        def verify_plots():
            page.wait_for_function('''() => ["viz1", "viz2"].every(k => {
                const plot = document.getElementById(k + "-plot");
                return plot.offsetWidth > 0 && plot.data?.some(t => t.x?.length > 0);
            })''', timeout=60000)
            assert page.locator('[data-tab-panel="explore"]').is_visible()
            assert tab.is_enabled() and tab.inner_text() == 'Explore'
            assert page.evaluate('!referenceRerunPending')
            return page.evaluate('''() => ({
                reference: document.getElementById("reference-select").value,
                gene: document.getElementById("viz2-gene-query").value,
                traces: ["viz1", "viz2"].map(k => document.getElementById(k + "-plot").data.length)
            })''')

        initial = verify_plots()
        if page.evaluate('["queued", "processing"].includes(currentDifferentialState?.status)'):
            assert page.evaluate('pollTimer !== null')
            page.wait_for_timeout(2500)
            assert len([url for _, url in requests if '/status?t=' in url]) >= 2
        page.reload(wait_until='domcontentloaded')
        page.wait_for_function('areExploreResultsReady() && activeExplorerTab === "explore"', timeout=60000)
        restored = verify_plots()
        if page.evaluate('["queued", "processing"].includes(currentDifferentialState?.status)'):
            assert page.evaluate('pollTimer !== null')
        assert initial == restored
        assert not errors, errors
        assert not [url for method, url in requests if method == 'POST' and url.endswith(('/run', '/qc'))]
        assert not [url for _, url in requests if '/reference-preview?' in url]
        page.screenshot(path='/tmp/scalable-session-restored.png')
        print(json.dumps({'passed': True, 'job_id': job_id, 'restored': restored,
                          'visible_loading_state': True, 'reload_verified': True,
                          'unnecessary_reference_previews': 0, 'analysis_reruns': 0}))
        browser.close()


if __name__ == '__main__':
    check_saved_session(sys.argv[1].rstrip('/'), sys.argv[2])
