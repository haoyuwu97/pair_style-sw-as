#!/usr/bin/env python3
from __future__ import annotations
import json,time,traceback
from pathlib import Path
import requests
OUT=Path('citationmap_force_output');OUT.mkdir(exist_ok=True)
SCHOLAR='NSpr644AAAAJ'
S=requests.Session();S.headers.update({'User-Agent':'Mozilla/5.0 Chrome/131.0.0.0 Safari/537.36','Accept':'application/json','Content-Type':'application/json'})
log=[]
try:
 r=S.post('https://citationmap.com/api/refresh-scholar?force=1',json={'scholarId':SCHOLAR},timeout=120)
 result={'status_code':r.status_code,'headers':dict(r.headers),'text':r.text}
 try:result['json']=r.json()
 except Exception:pass
except Exception as e:
 result={'error':repr(e),'trace':traceback.format_exc(limit=3)}
(OUT/'force_response.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
# Poll public progress for up to 20 minutes. Stop if complete or clearly denied/payment-gated.
for i in range(80):
 item={'iteration':i,'elapsed_seconds':i*15}
 for name,path in [('active','active-crawl'),('coverage','paper-coverage')]:
  try:
   rr=S.get(f'https://api.citationmap.com/scholars/{SCHOLAR}/{path}',timeout=60)
   item[name+'_status']=rr.status_code
   if rr.status_code==200:
    data=rr.json();item[name]=data
    (OUT/f'{name}_{i:03d}.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
  except Exception as e:item[name+'_error']=repr(e)
 log.append(item);print(json.dumps(item,ensure_ascii=False),flush=True)
 cov=item.get('coverage') or {}
 if cov.get('total_cites_crawled',0)>=cov.get('total_cites_gs',299):break
 if result.get('status_code') in (401,402,403) and not item.get('active',{}).get('task_id'):break
 time.sleep(15)
(OUT/'poll_log.json').write_text(json.dumps(log,ensure_ascii=False,indent=2),encoding='utf-8')
# Save final complete citing list/coverage snapshot.
for name,path in [('coverage_final','paper-coverage'),('citing_final','citing-papers?limit=1000'),('stats_final','stats'),('map_final','citing-map')]:
 try:
  rr=S.get(f'https://api.citationmap.com/scholars/{SCHOLAR}/{path}',timeout=90)
  data=rr.json() if rr.status_code==200 else {'status_code':rr.status_code,'text':rr.text}
 except Exception as e:data={'error':repr(e)}
 (OUT/f'{name}.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
