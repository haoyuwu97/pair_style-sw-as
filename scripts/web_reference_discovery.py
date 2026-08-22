#!/usr/bin/env python3
from __future__ import annotations
import hashlib, html, io, json, random, re, time, traceback
from pathlib import Path
from urllib.parse import quote_plus, urljoin, urlparse, parse_qs

import requests
from bs4 import BeautifulSoup
from pypdf import PdfReader

OUT=Path('web_reference_discovery_output'); OUT.mkdir(exist_ok=True)
PAPERS=[
('P01','10.1021/acsnano.2c11897','Multichannel Flexible Pulse Perception Array for Intelligent Disease Diagnosis System'),
('P02','10.1021/acs.langmuir.1c01697','Influence of Surface Defects on the Thermal Conductivity of Hexagonal Boron Nitride/Poly(dimethylsiloxane) Nanocomposites: A Molecular Dynamics Simulation'),
('P03','10.1016/j.polymer.2022.124750','Effect of Shape and Size of Nanofillers on the Viscoelasticity of Polymer Nanocomposites'),
('P04','10.1016/j.polymer.2020.122571','Molecular Dynamics Simulation of Fracture Mechanism in the Double Interpenetrated Cross-Linked Polymer'),
('P05','10.1016/j.polymer.2021.124129','Rheological Mechanism of Polymer Nanocomposites Filled with Spherical Nanoparticles: Insight from Molecular Dynamics Simulation'),
('P06','10.1021/acs.macromol.4c01541','Molecular Insights into the Topological Transition, Fracture, and Self-Healing Behavior of Vitrimer Composites with Exchangeable Interfaces'),
('P07','10.3390/ma14123301','Percolated Network of Mixed Nanoparticles with Different Sizes in Polymer Nanocomposites: A Coarse-Grained Molecular Dynamics Simulation'),
('P08','10.1016/j.compscitech.2020.108208','Percolation of Polydisperse Nanorods in Polymer Nanocomposites: Insights from Molecular Dynamics Simulation'),
('P09','10.1016/j.polymer.2021.124433','Structure and Dynamics Behavior During the Glass Transition of Polyisoprene in the Presence of Pressure: A Molecular Dynamics Simulation'),
('P10','10.1016/j.compscitech.2022.109694','Manipulating the Percolated Network of Nanorods in Polymer Matrix by Adding Non-Conductive Nanospheres: A Molecular Dynamics Simulation'),
('P11','10.1016/j.engfracmech.2020.107229','Improving the Fracture Property of Polymer Nanocomposites by Employing Nanoparticles as Cross-Linking Points'),
]
S=requests.Session(); S.headers.update({'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/131 Safari/537.36','Accept-Language':'en-US,en;q=0.9'})

def clean_url(u):
    if not u:return ''
    if u.startswith('//'):u='https:'+u
    if 'bing.com/ck/a' in u:
        # Keep redirect URL; requests will resolve it during verification.
        return u
    if 'duckduckgo.com/l/?' in u:
        q=parse_qs(urlparse(u).query); return q.get('uddg',[u])[0]
    return u

def search_bing(q,pages=4):
    rows=[]
    for pg in range(pages):
        url='https://www.bing.com/search'
        params={'q':q,'count':50,'first':pg*50+1,'setlang':'en-us'}
        r=S.get(url,params=params,timeout=60)
        (OUT/f'_bing_{hashlib.sha1((q+str(pg)).encode()).hexdigest()[:12]}.html').write_text(r.text,encoding='utf-8',errors='replace')
        soup=BeautifulSoup(r.text,'lxml')
        found=0
        for li in soup.select('li.b_algo'):
            a=li.select_one('h2 a')
            if not a: continue
            sn=li.select_one('.b_caption p')
            rows.append({'engine':'Bing','query':q,'title':a.get_text(' ',strip=True),'url':clean_url(a.get('href')),'snippet':sn.get_text(' ',strip=True) if sn else ''})
            found+=1
        if found==0:break
        time.sleep(1.5+random.random())
    return rows

def search_duck(q,pages=3):
    rows=[]
    for pg in range(pages):
        r=S.get('https://html.duckduckgo.com/html/',params={'q':q,'s':pg*30},timeout=60)
        (OUT/f'_duck_{hashlib.sha1((q+str(pg)).encode()).hexdigest()[:12]}.html').write_text(r.text,encoding='utf-8',errors='replace')
        soup=BeautifulSoup(r.text,'lxml')
        found=0
        for res in soup.select('.result'):
            a=res.select_one('.result__a')
            if not a:continue
            sn=res.select_one('.result__snippet')
            rows.append({'engine':'DuckDuckGo','query':q,'title':a.get_text(' ',strip=True),'url':clean_url(a.get('href')),'snippet':sn.get_text(' ',strip=True) if sn else ''})
            found+=1
        if found==0:break
        time.sleep(1.5+random.random())
    return rows

def extract_text(resp):
    ct=(resp.headers.get('content-type') or '').lower()
    data=resp.content[:12_000_000]
    if 'pdf' in ct or resp.url.lower().split('?')[0].endswith('.pdf') or data[:4]==b'%PDF':
        try:
            reader=PdfReader(io.BytesIO(data))
            texts=[]
            for page in reader.pages[:250]:
                try:texts.append(page.extract_text() or '')
                except Exception:pass
            return '\n'.join(texts),'pdf'
        except Exception as e:
            return '',f'pdf_error:{e!r}'
    try:
        txt=resp.text
    except Exception:
        txt=data.decode('utf-8','replace')
    soup=BeautifulSoup(txt,'lxml')
    for x in soup(['script','style','noscript']):x.decompose()
    return soup.get_text('\n',strip=True),'html'

def verify(url,doi,title):
    try:
        r=S.get(url,timeout=75,allow_redirects=True,headers={'Range':'bytes=0-11999999'})
        text,kind=extract_text(r)
        low=text.lower()
        doi_hit=doi.lower() in low or ('https://doi.org/'+doi.lower()) in low
        # Require a substantial title phrase, not just one keyword.
        core=re.sub(r'[^a-z0-9 ]+',' ',title.lower())
        core=' '.join(core.split())
        title_hit=core in re.sub(r'[^a-z0-9 ]+',' ',low)
        if not title_hit:
            words=core.split()
            phrase=' '.join(words[:min(12,len(words))])
            title_hit=phrase in re.sub(r'[^a-z0-9 ]+',' ',low)
        page_title=''
        if kind=='html':
            try:
                soup=BeautifulSoup(r.text,'lxml'); page_title=(soup.title.get_text(' ',strip=True) if soup.title else '')
            except Exception:pass
        return {'ok':bool(doi_hit or title_hit),'status':r.status_code,'final_url':r.url,'kind':kind,'doi_hit':doi_hit,'title_hit':title_hit,'length':len(r.content),'page_title':page_title,'text_head':text[:1200]}
    except Exception as e:
        return {'ok':False,'error':repr(e),'trace':traceback.format_exc(limit=1)}

summary={}; all_verified=[]; all_search=[]
for pid,doi,title in PAPERS:
    queries=[f'"{doi}" -doi.org',f'"{title}" -doi.org',f'"{doi}" references',f'"{title}" bibliography']
    raw=[]
    for q in queries:
        raw.extend(search_bing(q))
        raw.extend(search_duck(q))
    # Deduplicate URLs and remove obvious source/metadata-only landing pages.
    seen=set(); candidates=[]
    blocked_hosts={'doi.org','pubs.acs.org','sciencedirect.com','www.sciencedirect.com','api.crossref.org','search.crossref.org'}
    for x in raw:
        u=x.get('url','')
        if not u or u in seen:continue
        seen.add(u)
        host=urlparse(u).netloc.lower().removeprefix('www.')
        if host in blocked_hosts:continue
        if doi.lower() in u.lower():
            # Usually the source article or DOI resolver rather than a citing document.
            continue
        candidates.append(x)
    verified=[]
    for i,x in enumerate(candidates[:120]):
        v=verify(x['url'],doi,title)
        y={**x,'source_paper_id':pid,'source_doi':doi,'source_title':title,'verification':v}
        if v.get('ok'):
            verified.append(y); all_verified.append(y)
        all_search.append(y)
        time.sleep(.35+random.random()*.25)
    summary[pid]={'queries':queries,'raw_results':len(raw),'unique_candidates':len(candidates),'verified_mentions':len(verified)}
    (OUT/f'{pid}_verified.json').write_text(json.dumps(verified,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'all_verified.json').write_text(json.dumps(all_verified,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'all_search_results.json').write_text(json.dumps(all_search,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
