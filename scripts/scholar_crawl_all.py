#!/usr/bin/env python3
from __future__ import annotations
import asyncio, json, random, re, time, traceback
from pathlib import Path
from urllib.parse import urljoin, urlparse, parse_qs, urlencode, urlunparse
from playwright.async_api import async_playwright, TimeoutError as PlaywrightTimeoutError

OUT=Path('scholar_crawl_output'); OUT.mkdir(exist_ok=True)
PROFILE='NSpr644AAAAJ'
PAPERS=[
 {'id':'P01','token':'Tyk-4Ss8FVUC','count':138,'title':'Multichannel Flexible Pulse Perception Array for Intelligent Disease Diagnosis System'},
 {'id':'P02','token':'YsMSGLbcyi4C','count':33,'title':'Influence of Surface Defects on the Thermal Conductivity of Hexagonal Boron Nitride/Poly(dimethylsiloxane) Nanocomposites: A Molecular Dynamics Simulation'},
 {'id':'P03','token':'u5HHmVD_uO8C','count':28,'title':'Effect of Shape and Size of Nanofillers on the Viscoelasticity of Polymer Nanocomposites'},
 {'id':'P04','token':'2osOgNQ5qMEC','count':28,'title':'Molecular Dynamics Simulation of Fracture Mechanism in the Double Interpenetrated Cross-Linked Polymer'},
 {'id':'P05','token':'zYLM7Y9cAGgC','count':24,'title':'Rheological Mechanism of Polymer Nanocomposites Filled with Spherical Nanoparticles: Insight from Molecular Dynamics Simulation'},
 {'id':'P06','token':'W7OEmFMy1HYC','count':11,'title':'Molecular Insights into the Topological Transition, Fracture, and Self-Healing Behavior of Vitrimer Composites with Exchangeable Interfaces'},
 {'id':'P07','token':'UeHWp8X0CEIC','count':11,'title':'Percolated Network of Mixed Nanoparticles with Different Sizes in Polymer Nanocomposites: A Coarse-Grained Molecular Dynamics Simulation'},
 {'id':'P08','token':'IjCSPb-OGe4C','count':10,'title':'Percolation of Polydisperse Nanorods in Polymer Nanocomposites: Insights from Molecular Dynamics Simulation'},
 {'id':'P09','token':'qjMakFHDy7sC','count':9,'title':'Structure and Dynamics Behavior During the Glass Transition of Polyisoprene in the Presence of Pressure: A Molecular Dynamics Simulation'},
 {'id':'P10','token':'9yKSN-GCB0IC','count':4,'title':'Manipulating the Percolated Network of Nanorods in Polymer Matrix by Adding Non-Conductive Nanospheres: A Molecular Dynamics Simulation'},
 {'id':'P11','token':'u-x6o8ySG0sC','count':3,'title':'Improving the Fracture Property of Polymer Nanocomposites by Employing Nanoparticles as Cross-Linking Points'},
]

def save(name,obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')

def norm_url(href:str|None,base='https://scholar.google.com/'):
    return urljoin(base,href) if href else None

async def accept_consent(page):
    for label in ['I agree','Accept all','Accept','Agree']:
        try:
            btn=page.get_by_role('button',name=label,exact=False).first
            if await btn.count():
                await btn.click(timeout=3000); await page.wait_for_timeout(1500); return True
        except Exception: pass
    return False

async def blocked(page):
    text=(await page.locator('body').inner_text()).lower()
    return any(x in text for x in ['unusual traffic','not a robot','automated queries','recaptcha','your computer or network may be sending automated queries'])

async def parse_results(page):
    rows=[]
    loc=page.locator('div.gs_r.gs_or.gs_scl')
    n=await loc.count()
    for i in range(n):
        node=loc.nth(i)
        title_el=node.locator('.gs_rt').first
        title=await title_el.inner_text() if await title_el.count() else ''
        title=re.sub(r'^\s*\[[^\]]+\]\s*','',title).strip()
        a=title_el.locator('a').first if await title_el.count() else None
        href=await a.get_attribute('href') if a and await a.count() else None
        authors_el=node.locator('.gs_a').first
        authors=await authors_el.inner_text() if await authors_el.count() else ''
        snippet_el=node.locator('.gs_rs').first
        snippet=await snippet_el.inner_text() if await snippet_el.count() else ''
        footer=node.locator('.gs_fl a')
        footer_links=[]
        for j in range(await footer.count()):
            e=footer.nth(j)
            footer_links.append({'text':(await e.inner_text()).strip(),'href':norm_url(await e.get_attribute('href'))})
        rows.append({'title':title,'url':href,'authors_venue':authors,'snippet':snippet,'footer_links':footer_links})
    return rows

async def main():
  state={'started_at':time.time(),'papers':{},'errors':[],'blocked':False}
  async with async_playwright() as p:
    browser=await p.chromium.launch(headless=True,args=[
      '--disable-blink-features=AutomationControlled','--disable-dev-shm-usage','--no-sandbox',
      '--lang=en-US,en','--window-size=1365,900'])
    context=await browser.new_context(
      user_agent='Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36',
      locale='en-US',timezone_id='America/Chicago',viewport={'width':1365,'height':900},
      extra_http_headers={'Accept-Language':'en-US,en;q=0.9'})
    await context.add_init_script("Object.defineProperty(navigator, 'webdriver', {get: () => undefined});")
    page=await context.new_page()
    # Establish Scholar cookies with the author profile.
    try:
      await page.goto(f'https://scholar.google.com/citations?user={PROFILE}&hl=en&pagesize=100&view_op=list_works',wait_until='domcontentloaded',timeout=90000)
      await accept_consent(page); await page.wait_for_timeout(3000)
      (OUT/'profile.html').write_text(await page.content(),encoding='utf-8')
      await page.screenshot(path=str(OUT/'profile.png'),full_page=True)
      state['profile_url']=page.url; state['profile_blocked']=await blocked(page)
    except Exception as e: state['errors'].append({'stage':'profile','error':repr(e),'trace':traceback.format_exc(limit=2)})

    for pi,paper in enumerate(PAPERS):
      pid=paper['id']; rec={'source':paper,'detail':{},'pages':[],'results':[]}
      state['papers'][pid]=rec
      detail=f'https://scholar.google.com/citations?view_op=view_citation&hl=en&user={PROFILE}&citation_for_view={PROFILE}:{paper["token"]}'
      try:
        await page.goto(detail,wait_until='domcontentloaded',timeout=90000)
        await accept_consent(page); await page.wait_for_timeout(random.randint(2200,4200))
        rec['detail']['url']=page.url; rec['detail']['blocked']=await blocked(page)
        rec['detail']['body']=(await page.locator('body').inner_text())[:50000]
        (OUT/f'{pid}_detail.html').write_text(await page.content(),encoding='utf-8')
        # Find exact Cited by link and cluster id.
        links=page.locator('a')
        candidates=[]
        for i in range(await links.count()):
          e=links.nth(i); txt=(await e.inner_text()).strip(); href=await e.get_attribute('href')
          if href and ('cites=' in href or txt.lower().startswith('cited by')):
            candidates.append({'text':txt,'href':norm_url(href,page.url)})
        rec['detail']['cited_candidates']=candidates
        chosen=next((x for x in candidates if 'cites=' in (x.get('href') or '')),None)
        if not chosen:
          raise RuntimeError('No cited-by cluster link found')
        cited_url=chosen['href']; rec['detail']['cited_url']=cited_url
      except Exception as e:
        rec['detail']['error']=repr(e); rec['detail']['trace']=traceback.format_exc(limit=3)
        save('crawl_state.json',state); continue

      # Scholar normally returns 10 results/page; follow the next button until expected count or no new records.
      start=0; seen=set(); expected=paper['count']; base_url=rec['detail']['cited_url']
      while start < expected + 20 and len(rec['results']) < expected:
        parsed=urlparse(base_url); q=parse_qs(parsed.query); q['start']=[str(start)]; q['hl']=['en'];
        page_url=urlunparse(parsed._replace(query=urlencode(q,doseq=True)))
        page_info={'start':start,'url':page_url}
        try:
          await page.goto(page_url,wait_until='domcontentloaded',timeout=90000)
          await accept_consent(page); await page.wait_for_timeout(random.randint(3500,6500))
          page_info['final_url']=page.url; page_info['blocked']=await blocked(page)
          page_info['title']=await page.title(); page_info['body_head']=(await page.locator('body').inner_text())[:4000]
          html=await page.content(); (OUT/f'{pid}_start_{start:03d}.html').write_text(html,encoding='utf-8')
          rows=await parse_results(page); page_info['parsed_count']=len(rows)
          new=0
          for row in rows:
            key=(row['title'].lower(),row.get('authors_venue','').lower())
            if key not in seen:
              seen.add(key); row['source_id']=pid; row['scholar_start']=start; rec['results'].append(row); new+=1
          page_info['new_count']=new
          if page_info['blocked']:
            state['blocked']=True; rec['pages'].append(page_info); break
          if len(rows)==0 or new==0:
            rec['pages'].append(page_info); break
          rec['pages'].append(page_info)
          start+=10
          save('crawl_state.json',state)
          await page.wait_for_timeout(random.randint(2500,5000))
        except Exception as e:
          page_info['error']=repr(e); page_info['trace']=traceback.format_exc(limit=2); rec['pages'].append(page_info); break
      rec['retrieved']=len(rec['results']); rec['expected']=expected; rec['gap']=max(0,expected-len(rec['results']))
      save(f'{pid}_results.json',rec)
      save('crawl_state.json',state)
      if state['blocked']:
        break
      await page.wait_for_timeout(random.randint(4000,8000))

    state['finished_at']=time.time(); state['duration_seconds']=state['finished_at']-state['started_at']
    state['total_retrieved']=sum(x.get('retrieved',0) for x in state['papers'].values())
    state['total_expected']=sum(x['source']['count'] for x in state['papers'].values())
    save('crawl_state.json',state)
    await context.storage_state(path=str(OUT/'storage_state.json'))
    await browser.close()
  print(json.dumps({'total_retrieved':state.get('total_retrieved'),'total_expected':state.get('total_expected'),'blocked':state.get('blocked'),'gaps':{k:v.get('gap') for k,v in state['papers'].items()}},ensure_ascii=False,indent=2))

if __name__=='__main__': asyncio.run(main())
