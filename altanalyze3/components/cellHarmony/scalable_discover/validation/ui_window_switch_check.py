"""UI check in headless Chrome: each Explore plot type survives Windows 2 -> 1, and a pathway
opened from Chat is cleared when the panel switches to Cell frequency.

  python ui_window_switch_check.py --base-url http://127.0.0.1:8010 --job-id <id> --out <json>

A headless page is never a background tab, so Chrome does not throttle its timers.
"""
import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright

SNAP = """() => {
  const plot = document.getElementById("viz1-plot");
  const layout = plot && plot.layout;
  const title = layout ? (typeof layout.title === "string" ? layout.title : (layout.title && layout.title.text) || "") : "";
  const ticks = layout && layout.yaxis && layout.yaxis.ticktext || [];
  return {
    title: String(title).slice(0, 80),
    annotations: (layout && layout.annotations || []).map(a => a.text).join(" | ").slice(0, 160),
    traces: (plot && plot.data || []).length,
    iframe: !!(plot && plot.querySelector("iframe")),
    cytoscape: !!expressionCyByPanel.viz1,
    integrated_class: !!(plot && plot.classList.contains("integrated-view")),
    pathway_markup_left: !!(plot && plot.querySelector(".integrated-controls")),
    left_margin: layout && layout.margin ? layout.margin.l : null,
    longest_tick_chars: ticks.reduce((m, t) => Math.max(m, String(t).length), 0),
    text: (plot && plot.innerText || "").replace(/\\s+/g, " ").slice(0, 120),
  };
}"""
BAD = ("not found", "unavailable", "complete an alignment", "visualization unavailable")


def is_bad(snap):
    blob = " ".join(str(snap.get(k, "")) for k in ("title", "annotations", "text")).lower()
    return any(word in blob for word in BAD)


ap = argparse.ArgumentParser()
ap.add_argument("--base-url", required=True)
ap.add_argument("--job-id", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
out_path = Path(args.out)
result = {"job_id": args.job_id, "modes": [], "pathway_then_frequency": None, "console_errors": []}

with sync_playwright() as p:
    browser = p.chromium.launch(channel="chrome", headless=True)
    page = browser.new_page(viewport={"width": 1600, "height": 1000})
    page.on("pageerror", lambda e: result["console_errors"].append(f"pageerror: {e}"))
    page.on("console", lambda m: result["console_errors"].append(m.text) if m.type == "error" else None)
    page.goto(f"{args.base_url.rstrip('/')}/?job_id={args.job_id}")
    page.wait_for_function("() => (document.getElementById('viz1-mode')?.options.length || 0) > 3", timeout=90000)
    page.wait_for_timeout(4000)
    default_gene = page.evaluate("() => discoverStatus && discoverStatus.default_gene") or "SFTPC"
    modes = page.eval_on_selector_all("#viz1-mode option", "els => els.map(e => e.value)")
    for mode in modes:
        errors_before = len(result["console_errors"])
        page.select_option("#viz-window-count", "2")
        page.wait_for_timeout(800)
        page.select_option("#viz1-mode", mode)
        page.wait_for_timeout(1500)
        if mode in ("expression_umap", "violin"):
            page.fill("#viz1-gene-query", default_gene)
            page.dispatch_event("#viz1-gene-query", "change")
        page.wait_for_timeout(8000)
        before = page.evaluate(SNAP)
        page.select_option("#viz-window-count", "1")
        page.wait_for_timeout(4000)
        after = page.evaluate(SNAP)
        result["modes"].append({
            "mode": mode, "before": before, "after": after,
            "same_title": before["title"] == after["title"],
            "same_trace_count": before["traces"] == after["traces"],
            "error_after_switch": is_bad(after) and not is_bad(before),
            "error_before_switch": is_bad(before),
            "console_errors": result["console_errors"][errors_before:],
        })
    page.select_option("#viz-window-count", "2")
    page.wait_for_timeout(800)

    # Chat -> pathway in panel 1 -> Cell frequency, the path reported on 2026-10-05.
    page.click('[data-tab="chat"]')
    page.wait_for_timeout(1500)
    example = page.locator("#chat-examples button", has_text="pathways").first
    scenario = {"example": example.inner_text() if example.count() else None}
    if example.count():
        example.click()
        page.wait_for_selector("#chat-table td button.ghost-btn", timeout=90000)
        scenario["pathway"] = page.locator("#chat-table td button.ghost-btn").first.inner_text()
        page.locator("#chat-table td button.ghost-btn").first.click()
        page.wait_for_timeout(8000)
        scenario["pathway_open"] = page.evaluate(SNAP)
        scenario["mode_label"] = page.evaluate(
            "() => document.getElementById('viz1-mode').selectedOptions[0]?.textContent")
        page.select_option("#viz1-mode", "frequency")
        page.wait_for_timeout(6000)
        scenario["after_frequency"] = page.evaluate(SNAP)
        page.screenshot(path=str(out_path.with_name("ui_frequency_after_pathway.png")), full_page=False)
    result["pathway_then_frequency"] = scenario
    browser.close()

result["summary"] = {
    "modes_checked": len(result["modes"]),
    "modes_with_new_error_after_switch": [m["mode"] for m in result["modes"] if m["error_after_switch"]],
    "modes_title_changed": [m["mode"] for m in result["modes"] if not m["same_title"]],
    "pathway_markup_left_after_frequency": (result["pathway_then_frequency"] or {}).get("after_frequency", {}).get("pathway_markup_left"),
    "console_errors": len(result["console_errors"]),
}
out_path.write_text(json.dumps(result, indent=2))
print(json.dumps(result["summary"], indent=2))
