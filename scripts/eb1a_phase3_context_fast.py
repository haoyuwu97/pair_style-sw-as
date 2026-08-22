#!/usr/bin/env python3
from __future__ import annotations

import csv
import html
import io
import json
import re
import runpy
import time
import unicodedata
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import fitz
import requests
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "scripts" / "eb1a_phase3_metadata_fast.py"
ns = runpy.run_path(str(BASE_SCRIPT), run_name="phase3_context_base")
OUTPUT = ROOT / "audit_output_phase3_context_fast"
OUTPUT.mkdir(parents=True, exist_ok=True)
HEADERS = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/125 Safari/537.36 Haoyu-Wu-Citation-Audit/1.0",
    "Accept": "text/html,application/xhtml+xml,application/pdf,application/xml,text/plain;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.9",
}
MAX_BYTES = 35 * 1024 * 1024
TIMEOUT = 12
PRIORITIES = {"High", "Medium"}


def clean(v: Any) -> str:
    if v is None:
        return ""
    s = str(v).strip()
    return "" if s.lower() in {"nan", "none", "null"} else s


def norm_ascii(s: Any) -> str:
    x = unicodedata.normalize("NFKD", clean(s)).encode("ascii", "ignore").decode().lower()
    return " ".join(re.sub(r"[^a-z0-9]+", " ", x).split())


def norm_doi(s: Any) -> str:
    x = clean(s).lower()
    x = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", x)
    x = re.sub(r"^doi:\s*", "", x)
    return x.rstrip(" .;,)")


def compact(s: Any) -> str:
    return re.sub(r"\s+", " ", clean(s)).strip()


def trim_25(s: str) -> str:
    words = compact(s).split()
    return " ".join(words[:25]) + (" …" if len(words) > 25 else "")


def sentence_split(text: str) -> list[str]:
    text = compact(text)
    return [x.strip() for x in re.split(r"(?<=[.!?])\s+(?=[A-Z0-9])", text) if x.strip()]


def surname(source: dict[str, Any]) -> str:
    authors = source.get("authors") or []
    if not authors:
        return ""
    a = authors[0]
    name = a.get("name", "") if isinstance(a, dict) else str(a)
    parts = norm_ascii(name).split()
    return parts[-1] if parts else ""


def title_tokens(title: str) -> list[str]:
    stop = {"the","of","and","in","on","a","an","to","for","with","from","by","into","via","effect","study","simulation","molecular","polymer"}
    toks = [t for t in norm_ascii(title).split() if len(t) >= 5 and t not in stop]
    return sorted(set(toks), key=lambda x: (-len(x), x))[:10]


def jina_url(url: str) -> str:
    if url.startswith("https://"):
        return "https://r.jina.ai/https://" + url[len("https://"):]
    if url.startswith("http://"):
        return "https://r.jina.ai/http://" + url[len("http://"):]
    return ""


def candidate_urls(row: dict[str, str]) -> list[str]:
    urls = []
    for u in [row.get("Full_Text_URL"), row.get("Article_URL")]:
        u = clean(u)
        if u and u not in urls:
            urls.append(u)
    doi = norm_doi(row.get("Citing_DOI"))
    if doi:
        u = f"https://doi.org/{doi}"
        if u not in urls:
            urls.append(u)
    # Prefer original open/full-text URL, then Jina reader fallbacks.
    out = []
    for u in urls:
        if u not in out:
            out.append(u)
    for u in urls:
        ju = jina_url(u)
        if ju and ju not in out:
            out.append(ju)
    return out[:5]


@dataclass
class Retrieved:
    status: str
    url: str = ""
    kind: str = ""
    text: str = ""
    html_doc: str = ""
    pages: list[str] | None = None
    error: str = ""


def fetch_bytes(url: str) -> tuple[bytes, str, str]:
    r = requests.get(url, headers=HEADERS, timeout=TIMEOUT, allow_redirects=True, stream=True)
    r.raise_for_status()
    chunks, total = [], 0
    for chunk in r.iter_content(128 * 1024):
        if not chunk:
            continue
        total += len(chunk)
        if total > MAX_BYTES:
            raise RuntimeError("response exceeds 35 MB")
        chunks.append(chunk)
    return b"".join(chunks), clean(r.headers.get("content-type")).lower(), r.url


def html_text(data: bytes, content_type: str) -> tuple[str, str]:
    enc = "utf-8"
    m = re.search(r"charset=([^;\s]+)", content_type)
    if m:
        enc = m.group(1).strip('"\'')
    try:
        doc = data.decode(enc, errors="replace")
    except LookupError:
        doc = data.decode("utf-8", errors="replace")
    # Jina markdown is plain text; preserve it without HTML parsing.
    if not data.lstrip().startswith(b"<"):
        return doc, ""
    soup = BeautifulSoup(doc, "lxml-xml" if "xml" in content_type else "lxml")
    for tag in soup(["script", "style", "noscript", "svg", "nav", "header", "footer"]):
        tag.decompose()
    return soup.get_text("\n", strip=True), doc


def retrieve(row: dict[str, str]) -> Retrieved:
    errors = []
    for url in candidate_urls(row):
        try:
            data, ctype, final_url = fetch_bytes(url)
            if data[:5] == b"%PDF-" or "application/pdf" in ctype:
                doc = fitz.open(stream=data, filetype="pdf")
                pages = [p.get_text("text") or "" for p in doc]
                text = "\n\n".join(pages)
                if len(norm_ascii(text)) < 700:
                    raise RuntimeError("too little extractable PDF text")
                return Retrieved("retrieved", final_url, "pdf", text, pages=pages)
            if any(x in ctype for x in ("text/html", "application/xhtml", "xml", "text/plain", "text/markdown")) or data.lstrip().startswith(b"<") or "r.jina.ai" in final_url:
                text, html_doc = html_text(data, ctype)
                low = text.lower()[:2500]
                if len(norm_ascii(text)) < 700:
                    raise RuntimeError("too little extracted text")
                if any(x in low for x in ("enable javascript and cookies", "just a moment", "access denied", "verify you are human")):
                    raise RuntimeError("anti-bot interstitial")
                return Retrieved("retrieved", final_url, "html_or_text", text, html_doc=html_doc)
            errors.append(f"{url}: unsupported {ctype}")
        except Exception as e:
            errors.append(f"{url}: {type(e).__name__}: {e}")
    return Retrieved("not_retrieved", error=" | ".join(errors[:5]))


def find_reference(text: str, source: dict[str, Any]) -> tuple[int, int, str]:
    low = text.lower()
    doi = norm_doi(source.get("doi"))
    if doi:
        positions = [m.start() for m in re.finditer(re.escape(doi), low)]
        if positions:
            pos = positions[-1]
            return pos, pos + len(doi), "source DOI"
    toks = title_tokens(source.get("title", ""))
    if toks:
        # Locate a window containing at least four distinctive source-title tokens.
        for m in re.finditer(r"\n|$", text):
            end = m.start()
            start = max(0, end - 700)
            window = norm_ascii(text[start:end])
            if sum(t in window for t in toks[:7]) >= min(4, len(toks)):
                return start, end, "source title tokens"
    first = surname(source)
    year = clean(source.get("year"))
    if first and year:
        ms = list(re.finditer(rf"\b{re.escape(first)}\b.{{0,220}}\b{re.escape(year)}\b", low, flags=re.I | re.S))
        if ms:
            m = ms[-1]
            return m.start(), m.end(), "source first-author/year"
    return -1, -1, ""


def reference_marker(text: str, ref_pos: int) -> tuple[str, str]:
    if ref_pos < 0:
        return "", ""
    chunk = text[max(0, ref_pos - 700): ref_pos + 250]
    patterns = [
        r"(?:^|\n)\s*\[(\d{1,4})\]\s*[^\n]{0,500}$",
        r"(?:^|\n)\s*(\d{1,4})[.)]\s*[^\n]{0,500}$",
        r"(?:^|\n)\s*\[?(\d{1,4})\]?\s+[^\n]{0,500}$",
    ]
    for patt in patterns:
        ms = list(re.finditer(patt, chunk, flags=re.M))
        if ms:
            return ms[-1].group(1), compact(ms[-1].group(0))
    return "", trim_25(chunk[-300:])


def references_start(text: str) -> int:
    ms = list(re.finditer(r"(?:^|\n)\s*(references|bibliography|literature cited)\s*(?:\n|$)", text, flags=re.I))
    return ms[-1].start() if ms else int(len(text) * 0.78)


def structured_html_context(doc: str, source: dict[str, Any]) -> tuple[str, str, str, str]:
    if not doc:
        return "", "", "", ""
    soup = BeautifulSoup(doc, "lxml-xml" if doc.lstrip().startswith("<?xml") else "lxml")
    doi = norm_doi(source.get("doi"))
    toks = title_tokens(source.get("title", ""))
    ref = None
    for el in soup.find_all(["ref", "li", "div", "p"]):
        txt = compact(el.get_text(" ", strip=True))
        nt = norm_ascii(txt)
        if doi and doi in txt.lower():
            ref = el; break
        if toks and sum(t in nt for t in toks[:7]) >= min(4, len(toks)):
            ref = el; break
    if ref is None:
        return "", "", "", ""
    rid = clean(ref.get("id"))
    if not rid:
        child = ref.find(attrs={"id": True}) or ref.find(attrs={"name": True})
        rid = clean(child.get("id") or child.get("name")) if child else ""
    if not rid:
        return "", "", "", "structured reference found without anchor"
    links = soup.find_all(lambda tag: tag.name in {"xref", "a"} and (clean(tag.get("rid")) == rid or clean(tag.get("href")).lstrip("#") == rid))
    for link in links:
        parent = link.find_parent(["p", "div", "sec"])
        if not parent:
            continue
        txt = compact(parent.get_text(" ", strip=True))
        if len(txt) < 25:
            continue
        link_txt = compact(link.get_text(" ", strip=True))
        sentences = sentence_split(txt)
        chosen = next((s for s in sentences if link_txt and link_txt in s), sentences[0] if sentences else txt)
        return rid, trim_25(chosen), txt[:1600], "structured HTML/XML anchor"
    return rid, "", "", "structured reference found; body anchor not resolved"


def body_context(text: str, marker: str, source: dict[str, Any], pages: list[str] | None) -> tuple[str, str, str, str]:
    body_end = references_start(text)
    body = text[:body_end]
    patterns = []
    if marker:
        patterns.extend([
            re.compile(rf"\[{re.escape(marker)}(?:\s*[-–,]\s*\d+)*\]"),
            re.compile(rf"\({re.escape(marker)}(?:\s*[-–,]\s*\d+)*\)"),
        ])
    first = surname(source)
    year = clean(source.get("year"))
    if first and year:
        patterns.extend([
            re.compile(rf"\b{re.escape(first)}\b[^.\n]{{0,100}}\b{re.escape(year)}\b", re.I),
            re.compile(rf"\b{re.escape(first)}\s+et\s+al\.?", re.I),
        ])
    match = None
    method = ""
    for patt in patterns:
        ms = list(patt.finditer(body))
        if ms:
            match = ms[0]
            method = "numeric marker" if marker and patt.pattern.startswith("\\[") else "author-year marker"
            break
    if not match:
        return "", "", "", ""
    lo, hi = max(0, match.start() - 1000), min(len(body), match.end() + 1000)
    chunk = compact(body[lo:hi])
    sentences = sentence_split(chunk)
    token = compact(match.group(0))
    chosen = next((s for s in sentences if token and token in s), sentences[0] if sentences else chunk)
    location = ""
    if pages:
        cumulative = 0
        for i, page in enumerate(pages, 1):
            if match.start() <= cumulative + len(page):
                location = f"PDF extracted-text page {i}"
                break
            cumulative += len(page) + 2
    return trim_25(chosen), chunk[:1600], location, method


def classify(text: str) -> tuple[str, int, str]:
    t = norm_ascii(text)
    adoption = ["based on", "following", "adopted", "employed", "we use", "we used", "we employ", "extended", "built upon", "using the method"]
    comparison = ["consistent with", "in agreement with", "similar to", "compared with", "validated", "confirmed", "corroborat"]
    specific = ["demonstrated", "reported", "showed", "found", "revealed", "observed", "proposed", "developed", "achieved"]
    background = ["recent studies", "has been widely", "have been widely", "for example", "such as", "various studies", "among others"]
    if any(k in t for k in adoption):
        return "Method adoption/extension candidate", 3, "Use/follow/extension wording detected; manual verification required."
    if any(k in t for k in comparison):
        return "Validation/comparison candidate", 3, "Comparison/agreement wording detected; manual verification required."
    if any(k in t for k in specific):
        return "Specific result/factual-support candidate", 2, "Specific finding or method attribution wording detected."
    if any(k in t for k in background):
        return "General background/review candidate", 1, "Broad background framing detected."
    return "Context located; manual interpretation required", 1, "Body citation located without a decisive automated function label."


def paraphrase(context: str, source: dict[str, Any]) -> str:
    t = norm_ascii(context)
    first = surname(source) or "the cited work"
    if any(k in t for k in ("based on", "following", "adopted", "employed", "we use", "we used", "extended")):
        return f"The passage appears to use, follow, or extend the approach associated with {first}."
    if any(k in t for k in ("consistent with", "in agreement with", "similar to", "compared with", "validated", "confirmed")):
        return f"The passage appears to compare its result with, or state agreement with, {first}."
    if any(k in t for k in ("demonstrated", "reported", "showed", "found", "revealed", "observed", "proposed", "developed")):
        return f"The passage appears to attribute a specific finding, result, or method to {first}."
    return "A body citation was located; its evidentiary function requires manual technical review."


def process(row: dict[str, str], source: dict[str, Any]) -> dict[str, Any]:
    ret = retrieve(row)
    out = {
        "Edge_ID": row["Edge_ID"], "Source_ID": row["Source_ID"],
        "Citing_Title": row["Citing_Title"], "Context_Priority": row.get("Context_Priority", ""),
        "Citing_DOI": row.get("Citing_DOI", ""), "Article_URL": row.get("Article_URL", ""),
        "Full_Text_URL": row.get("Full_Text_URL", ""), "Retrieval_Status": ret.status,
        "URL_Used": ret.url, "Retrieved_Kind": ret.kind, "Retrieval_Error": ret.error,
        "Reference_Found": "No", "Reference_Marker": "", "Reference_Fragment_25w": "",
        "Body_Context_Found": "No", "Citation_Excerpt_25w": "", "Context_Paraphrase": "",
        "Context_Location": "", "Extraction_Method": "", "Context_Class_Preliminary": "Unlocated",
        "Context_Strength_Preliminary": 0, "Classification_Rationale": "",
        "Manual_Review_Status": "Required",
    }
    if ret.status != "retrieved":
        return out

    marker, excerpt, window, method = structured_html_context(ret.html_doc, source)
    if excerpt:
        cclass, score, rationale = classify(excerpt + " " + window)
        out.update({
            "Reference_Found": "Yes", "Reference_Marker": marker,
            "Body_Context_Found": "Yes", "Citation_Excerpt_25w": excerpt,
            "Context_Paraphrase": paraphrase(window, source), "Extraction_Method": method,
            "Context_Class_Preliminary": cclass, "Context_Strength_Preliminary": score,
            "Classification_Rationale": rationale,
        })
        return out

    ref_pos, _, ref_method = find_reference(ret.text, source)
    if ref_pos < 0:
        return out
    marker, ref_fragment = reference_marker(ret.text, ref_pos)
    out.update({
        "Reference_Found": "Yes", "Reference_Marker": marker,
        "Reference_Fragment_25w": trim_25(ref_fragment), "Extraction_Method": ref_method,
        "Context_Class_Preliminary": "Reference found; body context unlocated",
        "Classification_Rationale": "Source reference located in extracted text; body marker not yet mapped.",
    })
    body_excerpt, window, loc, body_method = body_context(ret.text, marker, source, ret.pages)
    if body_excerpt:
        cclass, score, rationale = classify(body_excerpt + " " + window)
        out.update({
            "Body_Context_Found": "Yes", "Citation_Excerpt_25w": body_excerpt,
            "Context_Paraphrase": paraphrase(window, source), "Context_Location": loc,
            "Extraction_Method": f"{ref_method}; {body_method}",
            "Context_Class_Preliminary": cclass, "Context_Strength_Preliminary": score,
            "Classification_Rationale": rationale,
        })
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys()) if rows else ["No records"]
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def main() -> None:
    rows, _ = ns["decode_existing_seed"]()
    sources = ns["extract_sources"](ns["decode_existing_script"](), rows)
    targets = [r for r in rows if clean(r.get("Context_Priority")) in PRIORITIES and clean(r.get("Full_Text_URL"))]
    targets.sort(key=lambda r: (0 if r.get("Context_Priority") == "High" else 1, r.get("Edge_ID", "")))
    print(f"Targeting {len(targets)} High/Medium records with explicit full-text URLs", flush=True)
    results = []
    with ThreadPoolExecutor(max_workers=12) as ex:
        futs = {ex.submit(process, row, sources[row["Source_ID"]]): row["Edge_ID"] for row in targets}
        for done, fut in enumerate(as_completed(futs), 1):
            try:
                results.append(fut.result())
            except Exception as e:
                results.append({"Edge_ID": futs[fut], "Retrieval_Status": "error", "Retrieval_Error": f"{type(e).__name__}: {e}", "Reference_Found": "No", "Body_Context_Found": "No", "Manual_Review_Status": "Required"})
            if done % 10 == 0:
                print(f"processed {done}/{len(targets)}", flush=True)
    results.sort(key=lambda r: r.get("Edge_ID", ""))
    write_csv(OUTPUT / "phase3_context_fast_candidates.csv", results)
    candidates = [r for r in results if r.get("Body_Context_Found") == "Yes"]
    write_csv(OUTPUT / "phase3_context_fast_body_candidates.csv", candidates)
    refs_only = [r for r in results if r.get("Reference_Found") == "Yes" and r.get("Body_Context_Found") != "Yes"]
    write_csv(OUTPUT / "phase3_context_fast_reference_only.csv", refs_only)
    summary = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "targets": len(targets),
        "retrieved_fulltext": sum(r.get("Retrieval_Status") == "retrieved" for r in results),
        "reference_found": sum(r.get("Reference_Found") == "Yes" for r in results),
        "body_context_candidates": len(candidates),
        "reference_only_candidates": len(refs_only),
        "preliminary_classes": dict(Counter(r.get("Context_Class_Preliminary", "") for r in results)),
        "notes": [
            "All excerpts are capped at 25 words.",
            "Every automated context classification is preliminary and requires manual source verification.",
            "Automated candidates do not replace the separate manually verified context layer.",
        ],
    }
    (OUTPUT / "phase3_context_fast_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
