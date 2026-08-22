#!/usr/bin/env python3
"""Retrieve the exact public Google Scholar cited-by result pages for 11 source papers.

The script makes a small number of sequential requests with conservative delays. It does not
attempt to solve or bypass CAPTCHAs. If Scholar returns a block page, that fact is preserved in
the output rather than silently fabricating records.
"""
from __future__ import annotations

import csv
import hashlib
import json
import random
import re
import time
import unicodedata
from pathlib import Path
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup

OUT = Path("scholar_exact_output")
OUT.mkdir(exist_ok=True)

PAPERS = [
    {"id":"P01","target":138,"cluster":"18159090922494394375","title":"Multichannel Flexible Pulse Perception Array for Intelligent Disease Diagnosis System"},
    {"id":"P02","target":33,"cluster":"1633036088071490547","title":"Influence of Surface Defects on the Thermal Conductivity of Hexagonal Boron Nitride/Poly(dimethylsiloxane) Nanocomposites: A Molecular Dynamics Simulation"},
    {"id":"P03","target":28,"cluster":"934276711164682768","title":"Effect of Shape and Size of Nanofillers on the Viscoelasticity of Polymer Nanocomposites"},
    {"id":"P04","target":28,"cluster":"9272803533019284630","title":"Molecular Dynamics Simulation of Fracture Mechanism in the Double Interpenetrated Cross-Linked Polymer"},
    {"id":"P05","target":24,"cluster":"11455914691638211269","title":"Rheological Mechanism of Polymer Nanocomposites Filled with Spherical Nanoparticles: Insight from Molecular Dynamics Simulation"},
    {"id":"P06","target":11,"cluster":"11130286448459684035","title":"Molecular Insights into the Topological Transition, Fracture, and Self-Healing Behavior of Vitrimer Composites with Exchangeable Interfaces"},
    {"id":"P07","target":11,"cluster":"10009415780122137402","title":"Percolated Network of Mixed Nanoparticles with Different Sizes in Polymer Nanocomposites: A Coarse-Grained Molecular Dynamics Simulation"},
    {"id":"P08","target":10,"cluster":"10619014777101595467","title":"Percolation of Polydisperse Nanorods in Polymer Nanocomposites: Insights from Molecular Dynamics Simulation"},
    {"id":"P09","target":9,"cluster":"9485134047812485534","title":"Structure and Dynamics Behavior During the Glass Transition of Polyisoprene in the Presence of Pressure: A Molecular Dynamics Simulation"},
    {"id":"P10","target":4,"cluster":"8982616081193836646","title":"Manipulating the Percolated Network of Nanorods in Polymer Matrix by Adding Non-Conductive Nanospheres: A Molecular Dynamics Simulation"},
    {"id":"P11","target":3,"cluster":"3325302876559785853","title":"Improving the Fracture Property of Polymer Nanocomposites by Employing Nanoparticles as Cross-Linking Points"},
]

UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"
S = requests.Session()
S.headers.update({
    "User-Agent": UA,
    "Accept-Language": "en-US,en;q=0.9",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Referer": "https://scholar.google.com/",
})


def norm(s: str) -> str:
    s = unicodedata.normalize("NFKC", s or "").lower()
    s = re.sub(r"\[[^\]]+\]", " ", s)
    s = re.sub(r"[^\w]+", " ", s, flags=re.UNICODE)
    return re.sub(r"\s+", " ", s).strip()


def blocked(text: str) -> bool:
    low = text.lower()
    phrases = [
        "our systems have detected unusual traffic",
        "not a robot",
        "please show you're not a robot",
        "sorry...",
        "/sorry/",
        "recaptcha",
    ]
    return any(p in low for p in phrases)


def parse_result(node, paper_id: str, start: int, rank_on_page: int) -> dict:
    ri = node.select_one(".gs_ri") or node
    title_node = ri.select_one("h3.gs_rt")
    title = title_node.get_text(" ", strip=True) if title_node else ""
    title = re.sub(r"^\[[^\]]+\]\s*", "", title).strip()
    anchor = title_node.find("a") if title_node else None
    result_url = urljoin("https://scholar.google.com", anchor.get("href")) if anchor and anchor.get("href") else ""
    authors_line_node = ri.select_one(".gs_a")
    authors_line = authors_line_node.get_text(" ", strip=True) if authors_line_node else ""
    snippet_node = ri.select_one(".gs_rs")
    snippet = snippet_node.get_text(" ", strip=True) if snippet_node else ""
    footer = ri.select_one(".gs_fl")
    footer_text = footer.get_text(" ", strip=True) if footer else ""
    cited_by_url = ""
    versions_url = ""
    related_url = ""
    if footer:
        for a in footer.find_all("a", href=True):
            txt = a.get_text(" ", strip=True).lower()
            href = urljoin("https://scholar.google.com", a.get("href"))
            if txt.startswith("cited by"):
                cited_by_url = href
            elif "versions" in txt:
                versions_url = href
            elif txt == "related articles":
                related_url = href
    parent = node if "gs_r" in (node.get("class") or []) else node.find_parent(class_="gs_r")
    data_cid = (parent.get("data-cid") if parent else None) or node.get("data-cid") or ""
    data_rp = (parent.get("data-rp") if parent else None) or ""
    raw_key = "|".join([paper_id, data_cid, norm(title), norm(authors_line)])
    stable_id = hashlib.sha1(raw_key.encode("utf-8")).hexdigest()[:20]
    return {
        "stable_id": stable_id,
        "source_paper_id": paper_id,
        "scholar_cluster_result_id": data_cid,
        "scholar_data_rp": data_rp,
        "source_page_start": start,
        "rank_on_page": rank_on_page,
        "citing_title": title,
        "citing_url": result_url,
        "authors_venue_year_line": authors_line,
        "snippet": snippet,
        "footer_text": footer_text,
        "cited_by_url": cited_by_url,
        "versions_url": versions_url,
        "related_url": related_url,
    }

# Seed the public-profile cookie conservatively.
seed = S.get(
    "https://scholar.google.com/citations",
    params={"user":"NSpr644AAAAJ","hl":"en","pagesize":100,"view_op":"list_works"},
    timeout=60,
)
(OUT / "seed_profile.html").write_text(seed.text, encoding="utf-8", errors="replace")

summary = {"seed_status": seed.status_code, "seed_blocked": blocked(seed.text), "papers": {}, "total_target": sum(p["target"] for p in PAPERS)}
all_rows = []

for paper_index, paper in enumerate(PAPERS, start=1):
    pid = paper["id"]
    target = paper["target"]
    rows = []
    seen = set()
    diagnostics = []
    paper_blocked = False

    for start in range(0, target, 20):
        params = {
            "oi": "bibs",
            "hl": "en",
            "cites": paper["cluster"],
            "start": start,
            "num": 20,
            "scipsc": 1,
            "as_sdt": "2005",
            "sciodt": "0,5",
        }
        response = None
        for attempt in range(2):
            response = S.get("https://scholar.google.com/scholar", params=params, timeout=75)
            is_blocked = blocked(response.text)
            diagnostics.append({
                "start": start,
                "attempt": attempt + 1,
                "status": response.status_code,
                "url": response.url,
                "length": len(response.content),
                "blocked": is_blocked,
            })
            (OUT / f"{pid}_start_{start:03d}_attempt_{attempt+1}.html").write_text(response.text, encoding="utf-8", errors="replace")
            if response.status_code == 200 and not is_blocked:
                break
            if attempt == 0:
                time.sleep(18 + random.random() * 8)
        if response is None or response.status_code != 200 or blocked(response.text):
            paper_blocked = True
            break

        soup = BeautifulSoup(response.text, "lxml")
        result_nodes = soup.select("div.gs_r.gs_or.gs_scl")
        if not result_nodes:
            result_nodes = soup.select("div.gs_r")
        page_added = 0
        for idx, node in enumerate(result_nodes, start=1):
            row = parse_result(node, pid, start, idx)
            key = row["scholar_cluster_result_id"] or (norm(row["citing_title"]), norm(row["authors_venue_year_line"]))
            if not row["citing_title"] or key in seen:
                continue
            seen.add(key)
            rows.append(row)
            page_added += 1
        diagnostics[-1]["parsed_results"] = page_added

        # Stop if Scholar returned no rows or the requested target is reached.
        if page_added == 0 or len(rows) >= target:
            break
        time.sleep(4.5 + random.random() * 3.0)

    # Preserve everything actually returned; do not fabricate rows to match the target.
    for row in rows:
        row["source_target_count"] = target
        row["source_title"] = paper["title"]
        all_rows.append(row)
    paper_out = {
        "paper": paper,
        "retrieved": len(rows),
        "blocked": paper_blocked,
        "diagnostics": diagnostics,
        "results": rows,
    }
    (OUT / f"{pid}_scholar_results.json").write_text(json.dumps(paper_out, ensure_ascii=False, indent=2), encoding="utf-8")
    summary["papers"][pid] = {
        "target": target,
        "retrieved": len(rows),
        "blocked": paper_blocked,
        "pages_attempted": len({d["start"] for d in diagnostics}),
        "diagnostics": diagnostics,
    }
    # Longer inter-paper pause; the total request volume remains small.
    if paper_index < len(PAPERS):
        time.sleep(6.0 + random.random() * 4.0)

summary["total_retrieved"] = len(all_rows)
summary["exact_target_met"] = len(all_rows) == summary["total_target"] and all(v["retrieved"] == v["target"] for v in summary["papers"].values())
(OUT / "scholar_exact_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
(OUT / "scholar_exact_all.json").write_text(json.dumps(all_rows, ensure_ascii=False, indent=2), encoding="utf-8")

headers = [
    "stable_id","source_paper_id","source_target_count","source_title","scholar_cluster_result_id",
    "scholar_data_rp","source_page_start","rank_on_page","citing_title","citing_url",
    "authors_venue_year_line","snippet","footer_text","cited_by_url","versions_url","related_url"
]
with (OUT / "scholar_exact_all.csv").open("w", newline="", encoding="utf-8-sig") as fh:
    writer = csv.DictWriter(fh, fieldnames=headers, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(all_rows)

print(json.dumps(summary, ensure_ascii=False, indent=2))
