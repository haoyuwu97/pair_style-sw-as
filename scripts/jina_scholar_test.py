#!/usr/bin/env python3
from pathlib import Path
import json, urllib.parse, requests
OUT=Path('jina_scholar_test_output'); OUT.mkdir(exist_ok=True)
cluster='18159090922494394375'
targets=[
    f'https://r.jina.ai/http://scholar.google.com/scholar?cites={cluster}&hl=en&num=20',
    f'https://r.jina.ai/https://scholar.google.com/scholar?cites={cluster}&hl=en&num=20',
    f'https://r.jina.ai/http://scholar.google.com/citations?user=NSpr644AAAAJ&hl=en&pagesize=100&view_op=list_works',
    f'https://r.jina.ai/https://scholar.google.com/citations?user=NSpr644AAAAJ&hl=en&pagesize=100&view_op=list_works',
]
results=[]
for i,url in enumerate(targets):
    try:
        r=requests.get(url,timeout=120,headers={'User-Agent':'Mozilla/5.0','Accept':'text/plain,*/*'})
        (OUT/f'response_{i}.txt').write_text(r.text,encoding='utf-8',errors='replace')
        results.append({'url':url,'status':r.status_code,'final_url':r.url,'length':len(r.content),'content_type':r.headers.get('content-type'),'head':r.text[:1000]})
    except Exception as e:
        results.append({'url':url,'error':repr(e)})
(OUT/'summary.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(results,ensure_ascii=False,indent=2))
