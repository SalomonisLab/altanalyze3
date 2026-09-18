"""Verify optional output controls on a saved session without changing/rerunning it."""
import json
import sys
from playwright.sync_api import sync_playwright

url = sys.argv[1]
with sync_playwright() as p:
    browser = p.chromium.launch(channel='chrome', headless=True)
    page = browser.new_page(viewport={'width': 1700, 'height': 1100})
    errors = []
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto(url, wait_until='domcontentloaded')
    page.wait_for_function('areExploreResultsReady()', timeout=180000)
    page.click('[data-tab="run"]')
    page.get_by_text('Optional marker heatmap exports', exact=True).click()
    render = page.locator('[name="marker_render_heatmap"]')
    svg = page.locator('[name="marker_write_svg"]')
    dpi = page.locator('[name="marker_heatmap_dpi"]')
    count = page.locator('[name="marker_cells_per_cluster"]')
    # scALABLE now defaults to skipping the static render, so SVG/DPI start disabled.
    assert render.input_value() == 'false'
    assert svg.input_value() == 'true'
    assert dpi.input_value() == ''
    assert count.input_value() == '100'
    assert svg.is_disabled() and dpi.is_disabled()
    assert count.is_enabled()
    render.select_option('true')
    assert svg.is_enabled() and dpi.is_enabled()
    dpi.select_option('600')
    svg.select_option('false')
    count.select_option('25')
    render.select_option('false')
    assert svg.is_disabled() and dpi.is_disabled()
    assert count.is_enabled()
    render.select_option('true')
    assert svg.is_enabled() and dpi.is_enabled()
    assert dpi.input_value() == '600' and svg.input_value() == 'false'
    # Submission payload is inspected at the browser boundary; no live job is changed.
    captured = []
    def capture(route):
        captured.append(route.request.post_data_json)
        route.fulfill(status=400, content_type='application/json', body='{"detail":"Test submission intercepted"}')
    page.route('**/api/jobs/*/qc', capture)
    page.on('dialog', lambda dialog: dialog.dismiss())
    render.select_option('false')
    page.locator('#qc-form button[type="submit"]').click()
    page.wait_for_function('!document.querySelector("#qc-form button[type=submit]").disabled')
    assert captured and captured[0]['marker_render_heatmap'] is False
    assert captured[0]['marker_write_svg'] is False
    assert captured[0]['marker_heatmap_dpi'] == 600
    assert captured[0]['marker_cells_per_cluster'] == 25
    page.screenshot(path='/tmp/scalable-marker-output-options.png')
    assert not errors, errors
    print(json.dumps({'passed': True, 'submitted_options': captured[0]}))
    browser.close()
