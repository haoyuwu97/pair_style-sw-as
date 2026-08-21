#!/usr/bin/env python3
from __future__ import annotations
import base64, concurrent.futures as cf, io, json, re, time
from pathlib import Path
from urllib.parse import urlparse, parse_qs, unquote
import requests
from bs4 import BeautifulSoup
from pypdf import PdfReader

OUT=Path('web_gap_fast_v2_output'); OUT.mkdir(exist_ok=True)
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

def sess():
 s=requests.Session();s.headers.update({'User-Agent':UA,'Accept-Language':'en-US,en;q=0.9'});return s

def decode_bing(u):
 try:
  q=parse_qs(urlparse(u).query)
  val=q.get('u',[''])[0]
  if val.startswith('a1'):
   raw=val[2:]; raw += '='*((4-len(raw)%4)%4)
   out=base64.urlsafe_b64decode(raw).decode('utf-8','replace')
   if out.startswith('http'): return out
 except Exception:pass
 return u

def ntext(x): return re.sub(r'[^a-z0-9]+',' ',(x or '').lower()).strip()

def search_one(p):
 pid,doi,title=p;s=sess();rows=[]
 words=ntext(title).split();short=' '.join(words[:min(10,len(words))])
 queries=[f'"{title}"',f'"{title}" references',f'"{short}" "{doi}"',f'"{doi}" polymer nanocomposite']
 for qi,q in enumerate(queries):
  for first in (1,51):
   try:
    r=s.get('https://www.bing.com/search',params={'q':q,'count':50,'first':first,'setlang':'en-us'},timeout=50)
    soup=BeautifulSoup(r.text,'lxml')
    got=0
    for li in soup.select('li.b_algo'):
     a=li.select_one('h2 a');sn=li.select_one('.b_caption p')
     if not a:continue
     u=decode_bing(a.get('href',''))
     rows.append({'pid':pid,'doi':doi,'source_title':title,'query':q,'result_title':a.get_text(' ',strip=True),'url':u,'snippet':sn.get_text(' ',strip=True) if sn else ''});got+=1
    if got==0:break
   except Exception:break
 return rows

def extract(resp):
 data=resp.content[:12_000_000];ct=(resp.headers.get('content-type') or '').lower()
 if data[:4]==b'%PDF' or 'pdf' in ct or resp.url.lower().split('?')[0].endswith('.pdf'):
  try:
   rd=PdfReader(io.BytesIO(data));return '\n'.join((p.extract_text() or '') for p in rd.pages[:220]),'pdf'
  except Exception as e:return '',f'pdf-error:{e!r}'
 soup=BeautifulSoup(resp.text,'lxml');
 for x in soup(['script','style','noscript']):x.decompose()
 return soup.get_text('\n',strip=True),'html'

def verify(row):
 u=row['url'];pid=row['pid'];doi=row['doi'];title=row['source_title'];host=urlparse(u).netloc.lower().removeprefix('www.')
 if not u or host in {'doi.org','api.crossref.org','openalex.org','semanticscholar.org'}:return None
 # Exclude obvious source article pages.
 ul=u.lower()
 if doi.lower() in ul and any(h in host for h in ['sciencedirect','pubs.acs','mdpi']):return None
 try:
  r=sess().get(u,timeout=70,allow_redirects=True,headers={'User-Agent':UA,'Range':'bytes=0-11999999'})
  text,kind=extract(r);low=text.lower();compact=ntext(text);core=ntext(title);words=core.split();phrase=' '.join(words[:min(12,len(words))])
  doi_hit=doi.lower() in low or ('https://doi.org/'+doi.lower()) in low
  title_hit=core in compact or (len(phrase)>30 and phrase in compact)
  # Avoid counting the source paper itself even after redirects.
  page_title=''
  if kind=='html':
   try:
    sp=BeautifulSoup(r.text,'lxml');page_title=sp.title.get_text(' ',strip=True) if sp.title else ''
   except Exception:pass
  is_source=(ntext(page_title)==core or doi.lower() in r.url.lower()) and (host in {'sciencedirect.com','pubs.acs.org','mdpi.com','mdpi.org'} or 'doi.org' in r.url.lower())
  return {**row,'verified':bool((doi_hit or title_hit) and not is_source),'status':r.status_code,'final_url':r.url,'kind':kind,'doi_hit':doi_hit,'title_hit':title_hit,'is_source':is_source,'page_title':page_title,'text_head':text[:1800]}
 except Exception as e:return {**row,'verified':False,'error':repr(e)}

raw=[]
with cf.ThreadPoolExecutor(max_workers=5) as ex:
 for xs in ex.map(search_one,PAPERS):raw.extend(xs)
seen=set();cand=[]
for x in raw:
 key=(x['pid'],x['url'])
 if not x['url'] or key in seen:continue
 seen.add(key);cand.append(x)
# Prefer likely document/repository pages and non-source result titles.
selected=[]
for pid,doi,title in PAPERS:
 xs=[x for x in cand if x['pid']==pid]
 xs.sort(key=lambda x:(not any(k in x['url'].lower() for k in ['.pdf','.edu','repository','handle.net','researchgate','semanticscholar','ouci','core.ac']), ntext(x['result_title'])==ntext(title)))
 selected.extend(xs[:140])
with cf.ThreadPoolExecutor(max_workers=16) as ex:checked=[x for x in ex.map(verify,selected) if x]
verified=[x for x in checked if x.get('verified')]
summary={'raw':len(raw),'unique':len(cand),'selected':len(selected),'checked':len(checked),'verified':len(verified),'counts':{pid:sum(x['pid']==pid for x in verified) for pid,_,_ in PAPERS}}
(OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'verified.json').write_text(json.dumps(verified,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'checked.json').write_text(json.dumps(checked,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
