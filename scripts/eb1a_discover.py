#!/usr/bin/env python3
"""Discover the most reliable machine-readable citation routes for Haoyu Wu's 11 papers.

This script is intentionally read-only with respect to scholarly services. It saves raw
responses and compact diagnostics as a GitHub Actions artifact so the next audit pass can
choose the most complete source rather than silently inventing Scholar records.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
import traceback
from pathlib import Path
from urllib.parse import urljoin, quote

import requests
from bs4 import BeautifulSoup

OUT = Path("audit_output")
OUT.mkdir(exist_ok=True)

PROFILE_ID = "NSpr644AAAAJ"
SOURCE_PAPERS = [
    {"id":"P01","doi":"10.1021/acsnano.2c11897","title":"Multichannel Flexible Pulse Perception Array for Intelligent Disease Diagnosis System","gs":138},
    {"id":"P02","doi":"10.1021/acs.langmuir.1c01697","title":"Influence of Surface Defects on the Thermal Conductivity of Hexagonal Boron Nitride/Poly(dimethylsiloxane) Nanocomposites: A Molecular Dynamics Simulation","gs":33},
    {"id":"P03","doi":"10.1016/j.polymer.2022.124750","title":"Effect of Shape and Size of Nanofillers on the Viscoelasticity of Polymer Nanocomposites","gs":28},
    {"id":"P04","doi":"10.1016/j.polymer.2020.122571","title":"Molecular Dynamics Simulation of Fracture Mechanism in the Double Interpenetrated Cross-Linked Polymer","gs":28},
    {"id":"P05","doi":"10.1016/j.polymer.2021.124129","title":"Rheological Mechanism of Polymer Nanocomposites Filled with Spherical Nanoparticles: Insight from Molecular Dynamics Simulation","gs":24},
    {"id":"P06","doi":"10.1021/acs.macromol.4c01541","title":"Molecular Insights into the Topological Transition, Fracture, and Self-Healing Behavior of Vitrimer Composites with Exchangeable Interfaces","gs":11},
    {"id":"P07","doi":"10.3390/ma14123301","title":"Percolated Network of Mixed Nanoparticles with Different Sizes in Polymer Nanocomposites: A Coarse-Grained Molecular Dynamics Simulation","gs":11},
    {"id":"P08","doi":"10.1016/j.compscitech.2020.108208","title":"Percolation of Polydisperse Nanorods in Polymer Nanocomposites: Insights from Molecular Dynamics Simulation","gs":10},
    {"id":"P09","doi":"10.1016/j.polymer.2021.124433","title":"Structure and Dynamics Behavior During the Glass Transition of Polyisoprene in the Presence of Pressure: A Molecular Dynamics Simulation","gs":9},
    {"id":"P10","doi":"10.1016/j.compscitech.2022.109694","title":"Manipulating the Percolated Network of Nanorods in Polymer Matrix by Adding Non-Conductive Nanospheres: A Molecular Dynamics Simulation","gs":4},
    {"id":"P11","doi":"10.1016/j.engfracmech.2020.107229","title":"Improving the Fracture Property of Polymer Nanocomposites by Employing Nanoparticles as Cross-Linking Points","gs":3},
]

UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"
S = requests.Session()
S.headers.update({"User-Agent": UA, "Accept-Language":"en-US,en;q=0.9", "Accept":"text/html,application/json;q=0.9,*/*;q=0.8"})

summary: dict[str, object] = {"profile_id": PROFILE_ID, "papers": SOURCE_PAPERS, "errors": []}

def save_text(name: str, text: str) -> None:
    (OUT / name).write_text(text, encoding="utf-8", errors="replace")

def save_json(name: str, obj: object) -> None:
    (OUT / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")

def get(url: str, *, timeout: int = 45, params=None, headers=None):
    try:
        r = S.get(url, timeout=timeout, params=params, headers=headers, allow_redirects=True)
        return r
    except Exception as e:
        summary["errors"].append({"url":url,"error":repr(e),"trace":traceback.format_exc(limit=2)})
        return None

# 1. CitationMap frontend and profile discovery.
cm: dict[str, object] = {}
for label, url in [
    ("home", "https://citationmap.com/"),
    ("profile", f"https://citationmap.com/profile/{PROFILE_ID}"),
    ("embed", f"https://citationmap.com/embed/{PROFILE_ID}"),
]:
    r = get(url)
    if r is None:
        cm[label] = {"error":"request failed"}
        continue
    cm[label] = {"status":r.status_code,"url":r.url,"length":len(r.content),"content_type":r.headers.get("content-type")}
    save_text(f"citationmap_{label}.html", r.text)

home_html = (OUT / "citationmap_home.html").read_text(encoding="utf-8", errors="replace") if (OUT / "citationmap_home.html").exists() else ""
soup = BeautifulSoup(home_html, "html.parser")
script_urls = []
for tag in soup.find_all("script", src=True):
    script_urls.append(urljoin("https://citationmap.com/", tag.get("src")))
cm["script_urls"] = script_urls
api_strings = set()
for idx, url in enumerate(script_urls[:30]):
    r = get(url, timeout=60)
    if r is None or r.status_code != 200:
        continue
    text = r.text
    save_text(f"citationmap_script_{idx:02d}.js", text)
    # Extract likely endpoints, route fragments, and Scholar-related strings.
    patterns = [
        r"https?://[^\"'`\\\s]+",
        r"/api/[A-Za-z0-9_./?=&{}:-]+",
        r"/[A-Za-z0-9_-]*(?:crawl|scholar|citation|profile|generate)[A-Za-z0-9_./?=&{}:-]*",
    ]
    for pat in patterns:
        for m in re.findall(pat, text, flags=re.I):
            if len(m) < 500:
                api_strings.add(m)
cm["candidate_routes"] = sorted(api_strings)
save_text("citationmap_candidate_routes.txt", "\n".join(sorted(api_strings)))
summary["citationmap"] = cm

# 2. Direct Scholar profile request; save raw page and parse visible publications/cited-by links.
scholar_url = f"https://scholar.google.com/citations?user={PROFILE_ID}&hl=en&pagesize=100&view_op=list_works"
r = get(scholar_url, timeout=60)
scholar: dict[str, object] = {}
if r is not None:
    scholar.update({"status":r.status_code,"url":r.url,"length":len(r.content),"content_type":r.headers.get("content-type")})
    save_text("scholar_profile.html", r.text)
    sp = BeautifulSoup(r.text, "html.parser")
    scholar["page_title"] = sp.title.get_text(" ", strip=True) if sp.title else None
    pubs = []
    for row in sp.select("tr.gsc_a_tr"):
        title_el = row.select_one("a.gsc_a_at")
        cited_el = row.select_one("a.gsc_a_ac")
        if not title_el:
            continue
        pubs.append({
            "title": title_el.get_text(" ", strip=True),
            "detail_href": title_el.get("href"),
            "cited_by_text": cited_el.get_text(" ", strip=True) if cited_el else "",
            "cited_by_href": cited_el.get("href") if cited_el else None,
        })
    scholar["parsed_publications"] = pubs
summary["scholar_direct"] = scholar

# 3. Try scholarly package as a second Scholar client.
try:
    from scholarly import scholarly
    author = scholarly.search_author_id(PROFILE_ID)
    # Fill only basic/publication data to keep request volume small.
    author = scholarly.fill(author, sections=["basics", "indices", "counts", "publications"])
    slim = {
        "name":author.get("name"),
        "citedby":author.get("citedby"),
        "hindex":author.get("hindex"),
        "i10index":author.get("i10index"),
        "publications":[{
            "bib":p.get("bib"),
            "author_pub_id":p.get("author_pub_id"),
            "num_citations":p.get("num_citations")
        } for p in author.get("publications",[])],
    }
    summary["scholarly"] = {"ok":True,"data":slim}
    save_json("scholarly_author.json", slim)
except Exception as e:
    summary["scholarly"] = {"ok":False,"error":repr(e),"trace":traceback.format_exc(limit=4)}

# 4. OpenAlex: resolve each DOI and list every citing work through cursor paging.
openalex_counts = {}
for p in SOURCE_PAPERS:
    doi = p["doi"]
    resolved = get("https://api.openalex.org/works/" + quote("https://doi.org/" + doi, safe=":/"), params={"mailto":"hwu24@nd.edu"})
    if resolved is None or resolved.status_code != 200:
        openalex_counts[p["id"]] = {"resolve_status": None if resolved is None else resolved.status_code}
        continue
    work = resolved.json()
    wid = work.get("id", "").rsplit("/",1)[-1]
    rows = []
    cursor = "*"
    pages = 0
    while cursor and pages < 50:
        rr = get("https://api.openalex.org/works", params={
            "filter":f"cites:{wid}",
            "per-page":200,
            "cursor":cursor,
            "mailto":"hwu24@nd.edu",
            "select":"id,doi,title,display_name,publication_year,publication_date,type,authorships,primary_location,locations_count,open_access,ids,cited_by_count"
        }, timeout=60)
        if rr is None or rr.status_code != 200:
            break
        data = rr.json()
        rows.extend(data.get("results", []))
        new_cursor = data.get("meta",{}).get("next_cursor")
        pages += 1
        if not new_cursor or new_cursor == cursor:
            break
        cursor = new_cursor
    save_json(f"openalex_{p['id']}.json", rows)
    openalex_counts[p["id"]] = {
        "work_id":wid,
        "reported_cited_by_count":work.get("cited_by_count"),
        "retrieved":len(rows),
        "pages":pages,
    }
summary["openalex_counts"] = openalex_counts

# 5. Semantic Scholar Graph API citation lists.
s2_counts = {}
for p in SOURCE_PAPERS:
    doi = p["doi"]
    url = f"https://api.semanticscholar.org/graph/v1/paper/DOI:{quote(doi, safe='')}/citations"
    rr = get(url, params={
        "limit":1000,
        "fields":"title,year,venue,authors,externalIds,url,publicationTypes,publicationDate,openAccessPdf,citationCount"
    }, timeout=90)
    if rr is None:
        s2_counts[p["id"]] = {"status":None}
        continue
    item = {"status":rr.status_code,"length":len(rr.content)}
    if rr.status_code == 200:
        data = rr.json()
        save_json(f"semanticscholar_{p['id']}.json", data)
        item.update({"total":data.get("total"),"retrieved":len(data.get("data",[])),"offset":data.get("offset"),"next":data.get("next")})
    else:
        save_text(f"semanticscholar_{p['id']}_error.txt", rr.text)
    s2_counts[p["id"]] = item
    time.sleep(1.1)
summary["semantic_scholar_counts"] = s2_counts

# 6. OpenCitations COCI/Index v2. The response identifies citing DOI/IDs; metadata can be enriched later.
oc_counts = {}
for p in SOURCE_PAPERS:
    doi = p["doi"]
    url = f"https://api.opencitations.net/index/v2/citations/doi:{doi}"
    rr = get(url, params={"format":"json"}, timeout=90)
    if rr is None:
        oc_counts[p["id"]] = {"status":None}
        continue
    item = {"status":rr.status_code,"length":len(rr.content)}
    if rr.status_code == 200:
        try:
            data = rr.json()
        except Exception:
            data = []
        save_json(f"opencitations_{p['id']}.json", data)
        item["retrieved"] = len(data) if isinstance(data,list) else None
    else:
        save_text(f"opencitations_{p['id']}_error.txt", rr.text)
    oc_counts[p["id"]] = item
summary["opencitations_counts"] = oc_counts

save_json("discovery_summary.json", summary)
print(json.dumps(summary, ensure_ascii=False, indent=2))
