"""Exercise every offered Explore mode, both windows, and marker Chat in Chrome.

Requires Playwright and a local Chrome installation. Pass a saved-session URL.
Reports browser exceptions and failed API responses rather than treating HTTP 200
as proof that a plot rendered.
"""
import argparse
import json
from pathlib import Path
from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('url')
    parser.add_argument('--report', type=Path, default=Path('/tmp/scalable_browser.json'))
    parser.add_argument('--state', default='HSC-1')
    args = parser.parse_args()
    with sync_playwright() as p:
        browser = p.chromium.launch(channel='chrome', headless=True)
        page = browser.new_page(viewport={'width':1700, 'height':1150}, ignore_https_errors=True)
        errors, failures, views = [], [], []
        page.on('pageerror', lambda e: errors.append(str(e)))
        page.on('response', lambda r: failures.append([r.status, r.url])
                if r.status >= 400 and '/api/' in r.url else None)
        page.goto(args.url, wait_until='domcontentloaded')
        page.wait_for_function('areExploreResultsReady()', timeout=120000)
        page.click('[data-tab="explore"]')
        modalities = page.locator('#viz1-modality').evaluate('(s)=>Array.from(s.options).map(o=>o.value)')
        for modality in modalities:
            page.evaluate('''async(m)=>{
                document.getElementById('viz1-mode').value='expression_umap';
                updateExpressionModeOptions();
                document.getElementById('viz1-modality').value=m;
                updatePanelFeatureInput('viz1',{value:''});
                loadedGeneSuggestionsSignature='';
                updateExpressionModeOptions();
                await loadGeneSuggestions(getResultsJobId());
            }''', modality)
            modes = page.locator('#viz1-mode').evaluate('(s)=>Array.from(s.options).map(o=>o.value)')
            for mode in modes:
                actual = page.evaluate('''async({m,mode})=>{
                    document.getElementById('viz1-mode').value=mode;
                    updateExpressionModeOptions();
                    const select=document.getElementById('viz1-modality');
                    if(Array.from(select.options).some(o=>o.value===m)) select.value=m;
                    updatePanelFeatureInput('viz1',{value:''});
                    loadedGeneSuggestionsSignature='';
                    updateExpressionModeOptions();
                    await loadGeneSuggestions(getResultsJobId());
                    await loadVisualizationPanel('viz1');
                    return {modality:panelModality('viz1'), mode:getPanelSelectValue('viz1','mode'),
                            feature:document.getElementById('viz1-gene-query').value,
                            rendered:document.getElementById('viz1-plot').children.length};
                }''', {'m':modality, 'mode':mode})
                views.append(dict(requested_modality=modality, **actual))
                print(json.dumps(views[-1]), flush=True)
        page.select_option('#viz-window-count', '1')
        assert page.evaluate('singleWindowActive()')
        page.select_option('#viz-window-count', '2')
        assert not page.evaluate('singleWindowActive()')
        page.click('[data-tab="chat"]')
        page.fill('#chat-question', f'What is the best modality marker of {args.state}')
        page.click('#chat-send')
        page.wait_for_function('chatLastResult?.intent === "modality_markers" && document.getElementById("chat-plot")?.data?.[0]?.type === "bar"', timeout=60000)
        result = dict(views=views, page_errors=errors, failed_api_requests=failures,
                      chat=page.evaluate('({status:chatLastResult.status,intent:chatLastResult.intent})'))
        args.report.write_text(json.dumps(result, indent=2))
        browser.close()
        assert not errors and not failures, (errors, failures)
        assert result['chat']['status'] == 'ok'
        print(json.dumps(dict(passed=True, views=len(views), chat=result['chat'])), flush=True)


if __name__ == '__main__':
    main()
