#!/usr/bin/env python3
from __future__ import annotations
import json,time,traceback
from pathlib import Path
import requests
OUT=Path('citationmap_poll_output');OUT.mkdir(exist_ok=True)
BASE='https://api.citationmap.com/scholars/NSpr644AAAAJ'
S=requests.Session();S.headers.update({'User-Agent':'Mozilla/5.0 Chrome/131.0.0.0 Safari/537.36','Accept':'application/json'})
log=[]; final_cov=None; final_citing=None
for i in range(80):
    item={'iteration':i,'elapsed_seconds':i*15}
    try:
        r=S.get(BASE+'/paper-coverage',timeout=45); item['coverage_status']=r.status_code
        if r.status_code==200:
            cov=r.json(); final_cov=cov
            item['crawled']=cov.get('total_cites_crawled');item['gs']=cov.get('total_cites_gs')
            item['per_paper']=[{'title':p.get('title'),'gs':p.get('cited_by_count'),'crawled':p.get('citations_crawled'),'gap':p.get('gap'),'shown':p.get('citing_shown')} for p in cov.get('papers',[])]
            (OUT/f'coverage_{i:03d}.json').write_text(json.dumps(cov,ensure_ascii=False,indent=2),encoding='utf-8')
        r2=S.get(BASE+'/citing-papers?limit=1000',timeout=60);item['citing_status']=r2.status_code
        if r2.status_code==200:
            cp=r2.json();final_citing=cp;item['unique_citing_total']=cp.get('total');item['citing_rows']=len(cp.get('papers',[]))
            (OUT/f'citing_{i:03d}.json').write_text(json.dumps(cp,ensure_ascii=False,indent=2),encoding='utf-8')
    except Exception as e:
        item['error']=repr(e);item['trace']=traceback.format_exc(limit=2)
    log.append(item);print(json.dumps(item,ensure_ascii=False),flush=True)
    if item.get('crawled',0)>=item.get('gs',299) and item.get('citing_status')==200:
        break
    time.sleep(15)
if final_cov is not None:(OUT/'paper_coverage_final.json').write_text(json.dumps(final_cov,ensure_ascii=False,indent=2),encoding='utf-8')
if final_citing is not None:(OUT/'citing_papers_final.json').write_text(json.dumps(final_citing,ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'poll_log.json').write_text(json.dumps(log,ensure_ascii=False,indent=2),encoding='utf-8')
