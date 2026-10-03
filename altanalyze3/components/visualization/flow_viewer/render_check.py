#!/usr/bin/env python3
"""Headless proof that the viewer DRAWS. An HTTP 200 is not a rendered plot.

Each check loads the page in a real browser, drives the controls a user would use, and
reads pixels back off the canvas. A mode passes only when its canvas holds the ink the
mode is supposed to produce:

  density      many distinct colours (a pseudocolour ramp), not one flat fill
  labels       several distinct colours, from the categorical palette
  gate outlines the overlay canvas is non-empty once a gate's own plot is opened
  drawn gate   a dragged rectangle returns a count from the server
  CITE space   switching space redraws with the CITE cell count

Screenshots land beside the bundle so the result can be looked at, not just asserted.
"""
import argparse, asyncio, json, os, sys


JS_INK = """(sel)=>{const c=document.querySelector(sel);const g=c.getContext('2d');
 const d=g.getImageData(0,0,c.width,c.height).data;const s=new Set();let nonwhite=0;
 for(let i=0;i<d.length;i+=4){const k=(d[i]<<16)|(d[i+1]<<8)|d[i+2];
   if(d[i+3]>0&&k!==0xffffff){nonwhite++;s.add(k);} }
 return {colours:s.size,nonwhite:nonwhite,w:c.width,h:c.height};}"""


async def run(url, out_dir, headed=False):
    from playwright.async_api import async_playwright
    os.makedirs(out_dir, exist_ok=True)
    results, errors = [], []

    def rec(name, ok, detail):
        results.append((name, ok, detail))
        print("  %-28s %s  %s" % (name, "PASS" if ok else "FAIL", detail))

    async with async_playwright() as pw:
        br = await pw.chromium.launch(
            headless=not headed,
            executable_path="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
        pg = await br.new_page(viewport={"width": 1600, "height": 950})
        pg.on("console", lambda m: errors.append(m.text) if m.type == "error" else None)
        pg.on("pageerror", lambda e: errors.append(str(e)))
        await pg.goto(url, wait_until="networkidle")
        await pg.wait_for_timeout(2500)

        ink = await pg.evaluate(JS_INK, "#cv")
        rec("density plot draws", ink["colours"] > 20 and ink["nonwhite"] > 5000,
            "%d colours, %d inked px, canvas %dx%d"
            % (ink["colours"], ink["nonwhite"], ink["w"], ink["h"]))
        await pg.screenshot(path=os.path.join(out_dir, "01_density.png"))

        # aspect: the streak bug showed as a canvas bitmap stretched by CSS
        box = await pg.evaluate("""()=>{const c=document.querySelector('#cv');
            const r=c.getBoundingClientRect();const d=window.devicePixelRatio||1;
            return {cssw:Math.round(r.width),cssh:Math.round(r.height),
                    bw:c.width,bh:c.height,dpr:d};}""")
        ok = (abs(box["bw"] - box["cssw"] * box["dpr"]) <= 2 and
              abs(box["bh"] - box["cssh"] * box["dpr"]) <= 2)
        rec("canvas not stretched", ok,
            "bitmap %dx%d vs css %dx%d @dpr %g"
            % (box["bw"], box["bh"], box["cssw"], box["cssh"], box["dpr"]))

        n_labels = await pg.evaluate(
            """()=>{const s=document.querySelector('#col');
                    return [...s.options].map(o=>o.value).filter(v=>!v.startsWith('marker:')&&v!=='density').length;}""")
        await pg.select_option("#col", index=1)
        await pg.wait_for_timeout(1800)
        ink2 = await pg.evaluate(JS_INK, "#cv")
        rec("label colouring draws", ink2["colours"] > 5 and ink2["nonwhite"] > 5000,
            "%d colours, %d inked px (%d label sets offered)"
            % (ink2["colours"], ink2["nonwhite"], n_labels))
        await pg.screenshot(path=os.path.join(out_dir, "02_labels.png"))

        n_nodes = await pg.evaluate("()=>document.querySelectorAll('.node').length")
        rec("gate tree listed", n_nodes >= 50, "%d populations in the tree" % n_nodes)

        # open a manual population: FlowJo behaviour is that its own plot comes up
        await pg.evaluate("""()=>{const n=[...document.querySelectorAll('.node')]
            .find(x=>x.textContent.includes('DN1'));n.click();}""")
        await pg.wait_for_timeout(2500)
        ov = await pg.evaluate(JS_INK, "#ov")
        axes = await pg.evaluate("()=>[document.querySelector('#xc').value,document.querySelector('#yc').value]")
        hint = await pg.evaluate("()=>document.querySelector('#axhint').textContent")
        rec("gate outline drawn", ov["nonwhite"] > 200,
            "overlay %d px, axes now %s x %s" % (ov["nonwhite"], axes[0], axes[1]))
        rec("population drill-down", "showing" in hint and "gate" in hint, hint.strip()[:70])
        await pg.screenshot(path=os.path.join(out_dir, "03_population_gate.png"))

        # drag a rectangle gate and read the server's count back
        await pg.click("#t_rect")
        b = await pg.evaluate("""()=>{const r=document.querySelector('#ov').getBoundingClientRect();
            return {x:r.x,y:r.y,w:r.width,h:r.height};}""")
        await pg.mouse.move(b["x"] + b["w"] * .35, b["y"] + b["h"] * .30)
        await pg.mouse.down()
        await pg.mouse.move(b["x"] + b["w"] * .72, b["y"] + b["h"] * .70, steps=12)
        await pg.mouse.up()
        await pg.wait_for_timeout(2200)
        txt = await pg.evaluate("()=>document.querySelector('#stats')?document.querySelector('#stats').textContent:''")
        rec("drawn gate returns stats", "events" in txt, " ".join(txt.split())[:70])
        await pg.screenshot(path=os.path.join(out_dir, "04_drawn_gate.png"))

        # a virtual gate: its own plot, its outline, and its score against the manual tree
        vg = await pg.evaluate("""()=>{const g=S.gatesets;
            for(const [k,v] of Object.entries(g))if(v.kind==='virtual')return [k,v.nodes[0].path];
            return null;}""")
        if vg:
            await pg.evaluate("""([k,path])=>{const nodes=[...document.querySelectorAll('.node')];
                const n=nodes.find(x=>x.title.startsWith(path));if(n)n.click();}""", vg)
            await pg.wait_for_timeout(2500)
            ov2 = await pg.evaluate(JS_INK, "#ov")
            panel = await pg.evaluate("()=>document.querySelector('#right').textContent")
            rec("virtual gate drawn", ov2["nonwhite"] > 200 and "F1 vs manual" in panel,
                "%s · overlay %d px · %s" % (vg[1][:28], ov2["nonwhite"],
                                             "scored vs manual" if "F1 vs manual" in panel else "NOT scored"))
            await pg.screenshot(path=os.path.join(out_dir, "06_virtual_gate.png"))
        else:
            rec("virtual gate drawn", False, "no virtual gate set in the bundle")

        # CITE space
        spaces = await pg.evaluate("()=>[...document.querySelector('#space').options].map(o=>o.value)")
        cite = [s for s in spaces if s.startswith("cite")]
        if cite:
            await pg.select_option("#space", cite[0])
            await pg.wait_for_timeout(1200)
            await pg.select_option("#mode", "embedding")
            await pg.wait_for_timeout(2500)
            ink3 = await pg.evaluate(JS_INK, "#cv")
            n = await pg.evaluate("()=>document.querySelector('#n').textContent")
            rec("CITE embedding draws", ink3["colours"] > 10 and ink3["nonwhite"] > 3000,
                "%s | %d colours, %d inked px" % (n.strip(), ink3["colours"], ink3["nonwhite"]))
            await pg.screenshot(path=os.path.join(out_dir, "05_cite_umap.png"))
        else:
            rec("CITE embedding draws", False, "no cite space in the bundle")

        await br.close()

    real = [e for e in errors if "favicon" not in e]
    print("\n  console errors: %d%s" % (len(real), (" -> " + real[0][:90]) if real else ""))
    bad = [r for r in results if not r[1]]
    print("  %d of %d checks pass; screenshots in %s" % (len(results) - len(bad), len(results), out_dir))
    return 1 if (bad or real) else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--url", default="http://127.0.0.1:8085")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sys.exit(asyncio.run(run(a.url, a.out)))
