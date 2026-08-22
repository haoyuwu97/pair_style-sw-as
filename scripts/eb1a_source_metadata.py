#!/usr/bin/env python3
from __future__ import annotations
import json, traceback
from pathlib import Path
from urllib.parse import quote
import requests

OUT=Path('source_metadata_output'); OUT.mkdir(exist_ok=True)
DOIS={
'P01':'10.1021/acsnano.2c11897','P02':'10.1021/acs.langmuir.1c01697','P03':'10.1016/j.polymer.2022.124750',
'P04':'10.1016/j.polymer.2020.122571','P05':'10.1016/j.polymer.2021.124129','P06':'10.1021/acs.macromol.4c01541',
'P07':'10.3390/ma14123301','P08':'10.1016/j.compscitech.2020.108208','P09':'10.1016/j.polymer.2021.124433',
'P10':'10.1016/j.compscitech.2022.109694','P11':'10.1016/j.engfracmech.2020.107229',
'MISSING_OC_P01':'10.1002/adrr.202500027','MISSING_OC_P06':'10.2139/ssrn.5527620'}
S=requests.Session(); S.headers.update({'User-Agent':'HaoyuWu-EB1A-citation-audit/1.0 (mailto:hwu24@nd.edu)'})
summary={}
for key,doi in DOIS.items():
    item={}
    try:
        r=S.get('https://api.openalex.org/works/'+quote('https://doi.org/'+doi,safe=':/'),params={'mailto':'hwu24@nd.edu'},timeout=45)
        item['openalex_status']=r.status_code
        if r.status_code==200:
            data=r.json(); (OUT/f'openalex_source_{key}.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
            item['openalex_id']=data.get('id'); item['authors']=[{'id':(a.get('author') or {}).get('id'),'name':(a.get('author') or {}).get('display_name'),'orcid':(a.get('author') or {}).get('orcid')} for a in data.get('authorships') or []]
    except Exception as e:item['openalex_error']=repr(e)
    try:
        r=S.get('https://api.crossref.org/works/'+quote(doi,safe=''),params={'mailto':'hwu24@nd.edu'},timeout=45)
        item['crossref_status']=r.status_code
        if r.status_code==200:
            data=r.json().get('message',{}); (OUT/f'crossref_{key}.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
            item['crossref_title']=(data.get('title') or [None])[0]; item['crossref_authors']=data.get('author')
    except Exception as e:item['crossref_error']=repr(e)
    summary[key]=item
(OUT/'source_metadata_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
