#!/usr/bin/env python3
from __future__ import annotations
import json, re, time, traceback
from pathlib import Path
from urllib.parse import quote, urljoin
import requests
from bs4 import BeautifulSoup

OUT=Path('audit_output_fast'); OUT.mkdir(exist_ok=True)
PROFILE_ID='NSpr644AAAAJ'
PAPERS=[
('P01','10.1021/acsnano.2c11897',138),('P02','10.1021/acs.langmuir.1c01697',33),
('P03','10.1016/j.polymer.2022.124750',28),('P04','10.1016/j.polymer.2020.122571',28),
('P05','10.1016/j.polymer.2021.124129',24),('P06','10.1021/acs.macromol.4c01541',11),
('P07','10.3390/ma14123301',11),('P08','10.1016/j.compscitech.2020.108208',10),
('P09','10.1016/j.polymer.2021.124433',9),('P10','10.1016/j.compscitech.2022.109694',4),
('P11','10.1016/j.engfracmech.2020.107229',3)]
UA='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/131.0.0.0 Safari/537.36'
S=requests.Session(); S.headers.update({'User-Agent':UA,'Accept-Language':'en-US,en;q=0.9'})
summary={'profile_id':PROFILE_ID,'errors':[],'openalex':{},'semantic_scholar':{},'opencitations':{},'scholar':{},'citationmap':{}}

def req(url,params=None,timeout=25):
    try:
        return S.get(url,params=params,timeout=timeout,allow_redirects=True)
    except Exception as e:
        summary['errors'].append({'url':url,'error':repr(e)})
        return None

def dump(name,obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')

def text(name,s):
    (OUT/name).write_text(s,encoding='utf-8',errors='replace')

# Scholar direct profile (single bounded request)
r=req(f'https://scholar.google.com/citations?user={PROFILE_ID}&hl=en&pagesize=100&view_op=list_works',timeout=30)
if r:
    text('scholar_profile.html',r.text)
    soup=BeautifulSoup(r.text,'html.parser')
    pubs=[]
    for row in soup.select('tr.gsc_a_tr'):
        a=row.select_one('a.gsc_a_at'); c=row.select_one('a.gsc_a_ac')
        if a: pubs.append({'title':a.get_text(' ',strip=True),'detail_href':a.get('href'),'cited_by_text':c.get_text(' ',strip=True) if c else '', 'cited_by_href':c.get('href') if c else None})
    summary['scholar']={'status':r.status_code,'length':len(r.content),'title':soup.title.get_text(' ',strip=True) if soup.title else None,'publications':pubs}

# CitationMap pages and frontend route discovery
for label,url in [('home','https://citationmap.com/'),('profile',f'https://citationmap.com/profile/{PROFILE_ID}'),('embed',f'https://citationmap.com/embed/{PROFILE_ID}')]:
    r=req(url,timeout=25)
    if r:
        text(f'citationmap_{label}.html',r.text)
        summary['citationmap'][label]={'status':r.status_code,'url':r.url,'length':len(r.content)}
home=(OUT/'citationmap_home.html').read_text(encoding='utf-8',errors='replace') if (OUT/'citationmap_home.html').exists() else ''
soup=BeautifulSoup(home,'html.parser'); routes=set(); scripts=[]
for tag in soup.find_all('script',src=True): scripts.append(urljoin('https://citationmap.com/',tag['src']))
for i,u in enumerate(scripts[:20]):
    rr=req(u,timeout=25)
    if not rr or rr.status_code!=200: continue
    text(f'citationmap_script_{i:02d}.js',rr.text)
    for pat in [r'https?://[^"\'`\\\s]+',r'/api/[A-Za-z0-9_./?=&{}:-]+',r'/[A-Za-z0-9_-]*(?:crawl|scholar|citation|profile|generate)[A-Za-z0-9_./?=&{}:-]*']:
        routes.update(x for x in re.findall(pat,rr.text,flags=re.I) if len(x)<500)
summary['citationmap']['scripts']=scripts; summary['citationmap']['candidate_routes']=sorted(routes)
text('citationmap_candidate_routes.txt','\n'.join(sorted(routes)))

# OpenAlex complete incoming-citation lists
for pid,doi,gs in PAPERS:
    r=req('https://api.openalex.org/works/'+quote('https://doi.org/'+doi,safe=':/'),params={'mailto':'hwu24@nd.edu'},timeout=30)
    if not r or r.status_code!=200:
        summary['openalex'][pid]={'status':None if not r else r.status_code}; continue
    w=r.json(); wid=w.get('id','').rsplit('/',1)[-1]; rows=[]; cursor='*'; pages=0
    while cursor and pages<20:
        rr=req('https://api.openalex.org/works',params={'filter':f'cites:{wid}','per-page':200,'cursor':cursor,'mailto':'hwu24@nd.edu','select':'id,doi,title,display_name,publication_year,publication_date,type,authorships,primary_location,open_access,ids,cited_by_count'},timeout=35)
        if not rr or rr.status_code!=200: break
        d=rr.json(); rows.extend(d.get('results',[])); nxt=d.get('meta',{}).get('next_cursor'); pages+=1
        if not nxt or nxt==cursor: break
        cursor=nxt
    dump(f'openalex_{pid}.json',rows)
    summary['openalex'][pid]={'work_id':wid,'source_reported':w.get('cited_by_count'),'retrieved':len(rows),'pages':pages,'gs_baseline':gs}

# Semantic Scholar citation lists, with bounded retries
for pid,doi,gs in PAPERS:
    url=f'https://api.semanticscholar.org/graph/v1/paper/DOI:{quote(doi,safe="")}/citations'
    r=None
    for attempt in range(3):
        r=req(url,params={'limit':1000,'fields':'title,year,venue,authors,externalIds,url,publicationTypes,publicationDate,openAccessPdf,citationCount'},timeout=30)
        if r and r.status_code==200: break
        time.sleep(2+attempt*2)
    if r and r.status_code==200:
        d=r.json(); dump(f'semanticscholar_{pid}.json',d)
        summary['semantic_scholar'][pid]={'status':200,'total':d.get('total'),'retrieved':len(d.get('data',[])),'gs_baseline':gs}
    else:
        summary['semantic_scholar'][pid]={'status':None if not r else r.status_code,'gs_baseline':gs}
        if r: text(f'semanticscholar_{pid}_error.txt',r.text)
    time.sleep(0.7)

# OpenCitations incoming DOI relations
for pid,doi,gs in PAPERS:
    r=req(f'https://api.opencitations.net/index/v2/citations/doi:{doi}',params={'format':'json'},timeout=35)
    if r and r.status_code==200:
        try:d=r.json()
        except Exception:d=[]
        dump(f'opencitations_{pid}.json',d)
        summary['opencitations'][pid]={'status':200,'retrieved':len(d) if isinstance(d,list) else None,'gs_baseline':gs}
    else:
        summary['opencitations'][pid]={'status':None if not r else r.status_code,'gs_baseline':gs}
        if r:text(f'opencitations_{pid}_error.txt',r.text)

dump('fast_summary.json',summary)
print(json.dumps(summary,ensure_ascii=False,indent=2))
