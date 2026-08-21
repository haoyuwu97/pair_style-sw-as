#!/usr/bin/env python3
from __future__ import annotations
import json,re
from pathlib import Path
from urllib.parse import urljoin
import requests
from bs4 import BeautifulSoup
OUT=Path('citationmap_js_output');OUT.mkdir(exist_ok=True)
url='https://citationmap.com/profile/NSpr644AAAAJ'
s=requests.Session();s.headers['User-Agent']='Mozilla/5.0 Chrome/131.0.0.0 Safari/537.36'
r=s.get(url,timeout=60);r.raise_for_status();(OUT/'profile.html').write_text(r.text,encoding='utf-8')
soup=BeautifulSoup(r.text,'html.parser');scripts=[urljoin(url,t['src']) for t in soup.find_all('script',src=True)]
all_hits=[];all_routes=set()
keywords=['paper-coverage','citing-papers','refresh-scholar','active-crawl','api.citationmap.com','citation_for_view','citations_crawled','crawl','scrapingdog','scholar']
for i,u in enumerate(scripts):
    try:rr=s.get(u,timeout=60)
    except Exception as e:continue
    if rr.status_code!=200:continue
    txt=rr.text;(OUT/f'script_{i:02d}.js').write_text(txt,encoding='utf-8')
    for pat in [r'https?://[^"\'`\\\s]+',r'/api/[A-Za-z0-9_./?=&{}:$-]+',r'/scholars/[A-Za-z0-9_./?=&{}:$-]+']:
        all_routes.update(x for x in re.findall(pat,txt,flags=re.I) if len(x)<500)
    low=txt.lower()
    for kw in keywords:
        start=0
        while True:
            pos=low.find(kw.lower(),start)
            if pos<0:break
            all_hits.append({'script':u,'keyword':kw,'context':txt[max(0,pos-1000):pos+2000]})
            start=pos+len(kw)
(OUT/'routes.txt').write_text('\n'.join(sorted(all_routes)),encoding='utf-8')
(OUT/'hits.json').write_text(json.dumps(all_hits,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'summary.json').write_text(json.dumps({'scripts':scripts,'route_count':len(all_routes),'hit_count':len(all_hits)},indent=2),encoding='utf-8')
print(json.dumps({'script_count':len(scripts),'route_count':len(all_routes),'hit_count':len(all_hits)},indent=2))
