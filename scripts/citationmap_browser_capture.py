#!/usr/bin/env python3
from __future__ import annotations
import asyncio, json, time, traceback
from pathlib import Path
from urllib.parse import urlparse
from playwright.async_api import async_playwright

OUT=Path('citationmap_browser_output'); OUT.mkdir(exist_ok=True)
PROFILE='https://citationmap.com/profile/NSpr644AAAAJ'
network=[]

async def main():
    async with async_playwright() as p:
        browser=await p.chromium.launch(headless=True)
        context=await browser.new_context(
            user_agent='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/131.0.0.0 Safari/537.36',
            locale='en-US',
            viewport={'width':1440,'height':1100},
        )
        page=await context.new_page()

        def on_request(req):
            try:
                h=req.headers
                network.append({'kind':'request','ts':time.time(),'method':req.method,'url':req.url,
                    'resource_type':req.resource_type,'post_data':req.post_data,
                    'next_action':h.get('next-action'),'content_type':h.get('content-type')})
            except Exception as e:
                network.append({'kind':'request_error','error':repr(e)})

        async def capture_response(resp):
            try:
                ct=(resp.headers.get('content-type') or '').lower()
                item={'kind':'response','ts':time.time(),'status':resp.status,'url':resp.url,'content_type':ct}
                if any(x in ct for x in ['application/json','text/x-component','text/plain']) or '/api/' in resp.url:
                    try:
                        body=await resp.text()
                        item['body']=body[:2000000]
                    except Exception as e:
                        item['body_error']=repr(e)
                network.append(item)
            except Exception as e:
                network.append({'kind':'response_error','error':repr(e)})

        page.on('request',on_request)
        page.on('response',lambda resp: asyncio.create_task(capture_response(resp)))

        console=[]
        page.on('console',lambda msg: console.append({'type':msg.type,'text':msg.text,'ts':time.time()}))
        page.on('pageerror',lambda exc: console.append({'type':'pageerror','text':str(exc),'ts':time.time()}))

        await page.goto(PROFILE,wait_until='domcontentloaded',timeout=120000)
        await page.wait_for_timeout(10000)
        initial_text=await page.locator('body').inner_text()
        (OUT/'initial_body.txt').write_text(initial_text,encoding='utf-8')
        (OUT/'initial.html').write_text(await page.content(),encoding='utf-8')
        await page.screenshot(path=str(OUT/'initial.png'),full_page=True)

        # Expand non-payment diagnostic controls and open the citations tab. These actions are
        # read-only and often trigger the profile's lazy citation fetch.
        for pattern in ['Citation coverage by paper','Citations']:
            try:
                loc=page.get_by_text(pattern,exact=False).first
                if await loc.count():
                    await loc.click(timeout=10000)
                    await page.wait_for_timeout(3000)
            except Exception as e:
                console.append({'type':'click_error','text':f'{pattern}: {e}','ts':time.time()})

        checkpoints=[]
        for i in range(24):  # four minutes total
            await page.wait_for_timeout(10000)
            body=await page.locator('body').inner_text()
            checkpoints.append({'elapsed_seconds':(i+1)*10,'length':len(body),'coverage_lines':[line for line in body.splitlines() if '/' in line and ('299' in line or 'citation' in line.lower())][:20]})
            if any(token in body for token in ['299 / 299','299/299','Crawl complete','All citations crawled']):
                break
            if i in (5,11,17):
                try:
                    await page.reload(wait_until='domcontentloaded',timeout=120000)
                    await page.wait_for_timeout(5000)
                except Exception as e:
                    console.append({'type':'reload_error','text':str(e),'ts':time.time()})

        final_html=await page.content(); final_text=await page.locator('body').inner_text()
        (OUT/'final.html').write_text(final_html,encoding='utf-8')
        (OUT/'final_body.txt').write_text(final_text,encoding='utf-8')
        await page.screenshot(path=str(OUT/'final.png'),full_page=True)
        storage=await page.evaluate('''() => ({localStorage:{...localStorage},sessionStorage:{...sessionStorage},href:location.href})''')
        (OUT/'storage.json').write_text(json.dumps(storage,ensure_ascii=False,indent=2),encoding='utf-8')
        (OUT/'network.json').write_text(json.dumps(network,ensure_ascii=False,indent=2),encoding='utf-8')
        (OUT/'console.json').write_text(json.dumps(console,ensure_ascii=False,indent=2),encoding='utf-8')
        (OUT/'checkpoints.json').write_text(json.dumps(checkpoints,ensure_ascii=False,indent=2),encoding='utf-8')
        await context.storage_state(path=str(OUT/'storage_state.json'))
        await browser.close()
        print(json.dumps({'network_events':len(network),'console_events':len(console),'checkpoints':checkpoints[-3:] if checkpoints else []},ensure_ascii=False,indent=2))

if __name__=='__main__':
    asyncio.run(main())
