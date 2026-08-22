#!/usr/bin/env python3
from pathlib import Path
import json,time,requests
OUT=Path('citationmap_force_one_output');OUT.mkdir(exist_ok=True)
SID='NSpr644AAAAJ';S=requests.Session();S.headers.update({'User-Agent':'Mozilla/5.0','Accept':'application/json,text/plain,*/*','Content-Type':'application/json','Referer':f'https://citationmap.com/profile/{SID}','Origin':'https://citationmap.com'})
responses=[]
for url,body in [
 (f'https://citationmap.com/api/refresh-scholar?force=1',{'scholarId':SID,'source':'manual'}),
 (f'https://citationmap.com/api/scholars/{SID}/retry-free-crawl?force=1',None),
]:
 try:
  r=S.post(url,json=body,timeout=90) if body is not None else S.post(url,timeout=90)
  responses.append({'url':url,'status':r.status_code,'headers':dict(r.headers),'text':r.text[:20000]})
 except Exception as e:responses.append({'url':url,'error':repr(e)})
(OUT/'force_responses.json').write_text(json.dumps(responses,ensure_ascii=False,indent=2),encoding='utf-8')
log=[]
for i in range(45):
 row={'iteration':i,'elapsed':i*15}
 for name,url in [
  ('active',f'https://api.citationmap.com/scholars/{SID}/active-crawl'),
  ('coverage',f'https://api.citationmap.com/scholars/{SID}/paper-coverage'),
  ('stats',f'https://api.citationmap.com/scholars/{SID}/stats'),
 ]:
  try:
   r=S.get(url,timeout=60);row[name+'_status']=r.status_code;row[name]=r.json() if r.headers.get('content-type','').startswith('application/json') else r.text[:3000]
  except Exception as e:row[name+'_error']=repr(e)
 log.append(row)
 (OUT/f'poll_{i:02d}.json').write_text(json.dumps(row,ensure_ascii=False,indent=2),encoding='utf-8')
 cov=row.get('coverage') if isinstance(row.get('coverage'),dict) else {}
 if cov.get('total_cites_crawled',0)>=cov.get('total_cites_gs',10**9):break
 active=row.get('active') if isinstance(row.get('active'),dict) else {}
 # after at least 3 minutes, exit once no task is active and coverage has stopped changing for 4 polls
 if i>=12 and not active.get('task_id'):
  vals=[(x.get('coverage') or {}).get('total_cites_crawled') for x in log[-5:] if isinstance(x.get('coverage'),dict)]
  if len(vals)>=5 and len(set(vals))==1:break
 time.sleep(15)
(OUT/'poll_log.json').write_text(json.dumps(log,ensure_ascii=False,indent=2),encoding='utf-8')
# final full payloads
for name,url in [
 ('coverage_final',f'https://api.citationmap.com/scholars/{SID}/paper-coverage'),
 ('citing_final',f'https://api.citationmap.com/scholars/{SID}/citing-papers?limit=1000'),
 ('map_final',f'https://api.citationmap.com/scholars/{SID}/citing-map'),
 ('stats_final',f'https://api.citationmap.com/scholars/{SID}/stats'),
]:
 try:
  r=S.get(url,timeout=90);(OUT/(name+'.json')).write_text(r.text,encoding='utf-8',errors='replace')
 except Exception as e:(OUT/(name+'_error.txt')).write_text(repr(e),encoding='utf-8')
print(json.dumps({'responses':responses,'last':log[-1] if log else None},ensure_ascii=False,indent=2))
