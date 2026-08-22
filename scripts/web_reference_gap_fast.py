#!/usr/bin/env python3
from __future__ import annotations
import concurrent.futures as cf, hashlib, io, json, random, re, time, traceback
from pathlib import Path
from urllib.parse import urlparse, parse_qs
import requests
from bs4 import BeautifulSoup
from pypdf import PdfReader

OUT=Path('web_gap_fast_output'); OUT.mkdir(exist_ok=True)
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
]
UA='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/131 Safari/537.36'

def session():
 s=requests.Session();s.headers.update({'User-Agent':UA,'Accept-Language':'en-US,en;q=0.9'});return s

def unwrap(u):
 if not u:return ''
 if 'duckduckgo.com/l/?' in u:
  return parse_qs(urlparse(u).query).get('uddg',[u])[0]
 return u

def search(pid,doi,title):
 s=session();out=[]
 queries=[f'"{doi}"',f'"{title}"',f'"{doi}" filetype:pdf',f'"{title}" thesis OR dissertation OR proceedings']
 for q in queries:
  try:
   r=s.get('https://www.bing.com/search',params={'q':q,'count':50,'setlang':'en-us'},timeout=45)
   soup=BeautifulSoup(r.text,'lxml')
   for li in soup.select('li.b_algo'):
    a=li.select_one('h2 a');sn=li.select_one('.b_caption p')
    if a:out.append({'pid':pid,'doi':doi,'source_title':title,'engine':'Bing','query':q,'result_title':a.get_text(' ',strip=True),'url':unwrap(a.get('href')),'snippet':sn.get_text(' ',strip=True) if sn else ''})
  except Exception:pass
  try:
   r=s.get('https://html.duckduckgo.com/html/',params={'q':q},timeout=45)
   soup=BeautifulSoup(r.text,'lxml')
   for x in soup.select('.result'):
    a=x.select_one('.result__a');sn=x.select_one('.result__snippet')
    if a:out.append({'pid':pid,'doi':doi,'source_title':title,'engine':'DuckDuckGo','query':q,'result_title':a.get_text(' ',strip=True),'url':unwrap(a.get('href')),'snippet':sn.get_text(' ',strip=True) if sn else ''})
  except Exception:pass
 return out

def text_norm(x):return re.sub(r'[^a-z0-9]+',' ',(x or '').lower()).strip()

def verify(row):
 u=row['url'];doi=row['doi'];title=row['source_title'];host=urlparse(u).netloc.lower().removeprefix('www.')
 bad={'doi.org','pubs.acs.org','sciencedirect.com','api.crossref.org','search.crossref.org','openalex.org','semanticscholar.org'}
 if not u or host in bad:return None
 try:
  s=session();r=s.get(u,timeout=55,allow_redirects=True,headers={'User-Agent':UA,'Range':'bytes=0-7999999'})
  data=r.content[:8_000_000];ct=(r.headers.get('content-type') or '').lower();kind='html';text=''
  if data[:4]==b'%PDF' or 'pdf' in ct or r.url.lower().split('?')[0].endswith('.pdf'):
   kind='pdf'
   try:
    rd=PdfReader(io.BytesIO(data));text='\n'.join((p.extract_text() or '') for p in rd.pages[:180])
   except Exception as e:return {**row,'verified':False,'status':r.status_code,'final_url':r.url,'kind':f'pdf-error:{e!r}'}
  else:
   soup=BeautifulSoup(r.text,'lxml');
   for x in soup(['script','style','noscript']):x.decompose()
   text=soup.get_text('\n',strip=True)
  low=text.lower();compact=text_norm(text);core=text_norm(title)
  words=core.split();phrase=' '.join(words[:min(12,len(words))])
  ok=(doi.lower() in low or ('https://doi.org/'+doi.lower()) in low or core in compact or (len(phrase)>30 and phrase in compact))
  return {**row,'verified':ok,'status':r.status_code,'final_url':r.url,'kind':kind,'doi_hit':doi.lower() in low,'title_hit':core in compact or phrase in compact,'text_head':text[:1500]}
 except Exception as e:return {**row,'verified':False,'error':repr(e)}

raw=[]
with cf.ThreadPoolExecutor(max_workers=6) as ex:
 for result in ex.map(lambda p:search(*p),PAPERS):raw.extend(result)
seen=set();candidates=[]
for x in raw:
 key=(x['pid'],x['url'])
 if not x['url'] or key in seen:continue
 seen.add(key);candidates.append(x)
# cap per paper, prefer PDF/edu/repository looking hits and non-source titles
selected=[]
for pid,doi,title in PAPERS:
 items=[x for x in candidates if x['pid']==pid]
 items.sort(key=lambda x:(not ('.pdf' in x['url'].lower() or '.edu' in x['url'].lower() or 'repository' in x['url'].lower() or 'thesis' in x['result_title'].lower()), x['result_title'].lower()==title.lower()))
 selected.extend(items[:80])
with cf.ThreadPoolExecutor(max_workers=12) as ex:
 verified_rows=list(ex.map(verify,selected))
verified=[x for x in verified_rows if x and x.get('verified')]
summary={pid:sum(1 for x in verified if x['pid']==pid) for pid,_,_ in PAPERS}
(OUT/'summary.json').write_text(json.dumps({'counts':summary,'raw':len(raw),'unique_candidates':len(candidates),'selected':len(selected),'verified':len(verified)},ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'verified.json').write_text(json.dumps(verified,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'all_checked.json').write_text(json.dumps([x for x in verified_rows if x],ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'counts':summary,'raw':len(raw),'unique_candidates':len(candidates),'selected':len(selected),'verified':len(verified)},ensure_ascii=False,indent=2))
