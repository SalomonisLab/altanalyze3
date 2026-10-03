"""Exercise real RNA/ADT views, population-driven gates and the marrow query in Chrome."""
import argparse,asyncio,json,os
from pathlib import Path
from playwright.async_api import async_playwright

INK="""()=>{const c=document.querySelector('#cv'),d=c.getContext('2d').getImageData(0,0,c.width,c.height).data;let n=0;for(let i=0;i<d.length;i+=4)if(d[i]<245||d[i+1]<245||d[i+2]<245)n++;return n;}"""

async def run(url,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);results=[];errors=[]
    async with async_playwright() as pw:
        browser=await pw.chromium.launch(headless=True,executable_path='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome')
        page=await browser.new_page(viewport={'width':1800,'height':1100})
        page.on('pageerror',lambda e:errors.append(str(e)))
        page.on('console',lambda m:errors.append(m.text) if m.type=='error' else None)
        await page.goto(url,wait_until='networkidle')
        catalog=await page.evaluate('S.man')
        for space,spec in catalog['spaces'].items():
            await page.select_option('#space',space);await page.wait_for_timeout(800)
            if not spec['embeddings']:continue
            await page.select_option('#mode','embedding')
            for embedding in spec['embeddings']:
                await page.select_option('#emb',embedding);await page.wait_for_timeout(650)
                ink=await page.evaluate(INK);assert ink>3000,(space,embedding,ink)
                assert await page.locator('#targetPanel').count()==1
                results.append(dict(check='embedding',space=space,embedding=embedding,ink=ink))
            if 'RNA:Spi1' in spec['features']:
                await page.select_option('#col','marker:RNA:Spi1');await page.wait_for_timeout(600)
                assert await page.evaluate(INK)>3000
                results.append(dict(check='measured PU.1 RNA',space=space))
            if space in ['cite_grimes','cite_chinese','cite_marrow_ADT195','cite_marrow_ADT112_reference_coordinates']:
                await page.screenshot(path=str(out/(space+'.png')))
        await page.select_option('#space','cite_marrow_ADT195');await page.wait_for_timeout(800)
        await page.select_option('#targetSet','Population');await page.wait_for_timeout(600)
        await page.select_option('#targetPop',['ML-1a']);await page.click('#highlightTarget');await page.wait_for_timeout(500)
        await page.click('#compareProtein');await page.wait_for_timeout(800)
        assert 'cite_marrow_ADT195' in await page.locator('#stats').inner_text()
        await page.screenshot(path=str(out/'marrow_ML1a_ADT_PU1_IRF.png'))
        await page.select_option('#targetMethod','kde_k5_average');await page.click('#useExampleGate');await page.click('#optimizeTarget')
        await page.wait_for_function("S.optimization&&S.optimization.best.method==='constrained_sequential'",timeout=30000)
        await page.wait_for_timeout(1200)
        assert await page.locator('#space').input_value()=='flow'
        assert len(await page.evaluate('S.optimization.best.conditions'))==8
        assert await page.evaluate('S.gateHighlight.reduce((a,b)=>a+b,0)')>0
        assert 'untouched test events' in (await page.locator('#optreport').inner_text()).lower()
        results.append(dict(check='marrow → flow constrained ML1a gate',result=await page.evaluate('S.optimization.best.test')))
        await page.screenshot(path=str(out/'ML1a_required_gate.png'))
        await page.select_option('#mode','embedding');await page.select_option('#emb','cite_marrow_ADT195_projected_flow');await page.wait_for_timeout(800)
        assert await page.evaluate(INK)>3000
        results.append(dict(check='selected gate persists on flow marrow projection'))
        await page.screenshot(path=str(out/'ML1a_flow_marrow_projection.png'))
        async with page.expect_download() as dl:
            await page.get_by_role('button',name='Export gate JSON',exact=True).click()
        await (await dl.value).save_as(str(out/'browser_exported_gate.json'))
        async with page.expect_download() as dl:
            await page.get_by_role('button',name='Export event IDs',exact=True).click()
        await (await dl.value).save_as(str(out/'browser_exported_events.tsv'))
        gate=json.loads((out/'browser_exported_gate.json').read_text());events=(out/'browser_exported_events.tsv').read_text().splitlines()
        assert len(events)-1==gate['best']['all_events']['n_gated']
        results.append(dict(check='exported gate and event IDs agree',n_events=len(events)-1))
        await page.click('#viewValidation');await page.wait_for_timeout(600);assert 'kde_k5_average' in await page.locator('#stats').inner_text()
        await browser.close()
    (out/'browser_checks.json').write_text(json.dumps(dict(results=results,errors=errors),indent=2))
    print(f'{len(results)} browser checks; console/page errors: {errors}',flush=True)
    if errors:raise RuntimeError(errors)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--url',default='http://127.0.0.1:8085');p.add_argument('--out',required=True);a=p.parse_args();asyncio.run(run(a.url,a.out))
