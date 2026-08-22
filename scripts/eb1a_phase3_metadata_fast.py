#!/usr/bin/env python3
from __future__ import annotations

import ast
import base64
import csv
import io
import json
import os
import re
import runpy
import time
import unicodedata
import zlib
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any
from urllib.parse import quote

import requests
from rapidfuzz import fuzz

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "audit_inputs" / "phase3_bundle"
OUTPUT = ROOT / "audit_output_phase3_metadata_fast"
OUTPUT.mkdir(parents=True, exist_ok=True)
OPENALEX = "https://api.openalex.org/works"
CROSSREF = "https://api.crossref.org/works"
MAILTO = os.environ.get("OPENALEX_MAILTO", "m18129227366@163.com")
HEADERS = {
    "User-Agent": f"Haoyu-Wu-EB1A-citation-audit/1.0 (mailto:{MAILTO})",
    "Accept": "application/json",
}
HAOYU_NAME = "haoyu wu"
HAOYU_ORCID = "0000-0002-2805-4911"
PIDS = [f"P{i:02d}" for i in range(1, 12)]


def clean(v: Any) -> str:
    if v is None:
        return ""
    s = str(v).strip()
    return "" if s.lower() in {"nan", "none", "null"} else s


def norm_ascii(s: Any) -> str:
    x = unicodedata.normalize("NFKD", clean(s)).encode("ascii", "ignore").decode().lower()
    return " ".join(re.sub(r"[^a-z0-9]+", " ", x).split())


def norm_name(s: Any) -> str:
    return norm_ascii(s)


def norm_orcid(s: Any) -> str:
    return clean(s).lower().replace("https://orcid.org/", "")


def norm_doi(s: Any) -> str:
    x = clean(s).lower()
    x = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", x)
    x = re.sub(r"^doi:\s*", "", x)
    return x.rstrip(" .;,)")


def title_score(a: str, b: str) -> float:
    return float(fuzz.token_set_ratio(norm_ascii(a), norm_ascii(b)))


def read_parts(pattern: str) -> str:
    parts = sorted(BUNDLE.glob(pattern))
    if not parts:
        raise RuntimeError(f"No bundle parts matching {pattern}")
    return re.sub(r"\s+", "", "".join(p.read_text(encoding="utf-8") for p in parts))


def decode_existing_seed() -> tuple[list[dict[str, str]], bytes]:
    encoded = read_parts("phase3_seed_v2.part*.b64")
    raw = zlib.decompress(base64.b64decode(encoded, validate=True))
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if rows and "\ufeffEdge_ID" in rows[0]:
        for row in rows:
            row["Edge_ID"] = row.pop("\ufeffEdge_ID")
    if len(rows) != 299:
        raise RuntimeError(f"Expected 299 exact Scholar edges, got {len(rows)}")
    return rows, raw


def decode_existing_script() -> bytes:
    encoded = read_parts("phase3_script.part*.b64")
    return zlib.decompress(base64.b64decode(encoded, validate=True))


def normalize_source_dict(obj: Any) -> dict[str, Any] | None:
    if not isinstance(obj, dict) or not all(pid in obj for pid in PIDS):
        return None
    out: dict[str, Any] = {}
    for pid in PIDS:
        src = obj.get(pid)
        if not isinstance(src, dict):
            return None
        raw_authors = src.get("authors") or src.get("source_authors") or src.get("author_list")
        if not isinstance(raw_authors, list) or not raw_authors:
            return None
        authors = []
        for a in raw_authors:
            if isinstance(a, str):
                authors.append({"name": a, "id": "", "orcid": ""})
            elif isinstance(a, dict):
                name = clean(a.get("name") or a.get("display_name") or a.get("author_name"))
                if name:
                    authors.append({
                        "name": name,
                        "id": clean(a.get("id") or a.get("author_id") or a.get("openalex_id")),
                        "orcid": norm_orcid(a.get("orcid") or a.get("ORCID")),
                    })
        if not authors:
            return None
        out[pid] = {
            "title": clean(src.get("title") or src.get("source_title")),
            "doi": norm_doi(src.get("doi") or src.get("source_doi")),
            "year": src.get("year") or src.get("source_year") or "",
            "authors": authors,
        }
    return out


def extract_sources(script_bytes: bytes, rows: list[dict[str, str]]) -> dict[str, Any]:
    generated = ROOT / "scripts" / "_phase3_fast_import_target.py"
    generated.write_bytes(script_bytes)
    ns = runpy.run_path(str(generated), run_name="phase3_fast_import_target")
    candidates = []
    for name, value in ns.items():
        normalized = normalize_source_dict(value)
        if normalized:
            candidates.append((name, normalized))
    if candidates:
        candidates.sort(key=lambda x: (0 if "source" in x[0].lower() else 1, len(x[0])))
        print(f"Using source metadata object: {candidates[0][0]}", flush=True)
        return candidates[0][1]

    tree = ast.parse(script_bytes.decode("utf-8"))
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            value_node = node.value
            try:
                value = ast.literal_eval(value_node)
            except Exception:
                continue
            normalized = normalize_source_dict(value)
            if normalized:
                print("Using source metadata recovered by AST literal evaluation", flush=True)
                return normalized

    possible_fields = [
        "Source_Authors_JSON", "Source_Authors", "Source_Author_List",
        "Source_OpenAlex_Authors", "Source_Metadata_JSON",
    ]
    by_pid: dict[str, Any] = {}
    for pid in PIDS:
        row = next(r for r in rows if r.get("Source_ID") == pid)
        authors: list[dict[str, str]] = []
        for field in possible_fields:
            val = clean(row.get(field))
            if not val:
                continue
            try:
                parsed = json.loads(val)
            except Exception:
                parsed = [x.strip() for x in val.split(";") if x.strip()]
            if isinstance(parsed, list):
                for a in parsed:
                    if isinstance(a, str):
                        authors.append({"name": a, "id": "", "orcid": ""})
                    elif isinstance(a, dict):
                        name = clean(a.get("name") or a.get("display_name"))
                        if name:
                            authors.append({"name": name, "id": clean(a.get("id")), "orcid": norm_orcid(a.get("orcid"))})
            if authors:
                break
        if not authors:
            raise RuntimeError("Could not recover source-author metadata from existing Phase 3 bundle")
        by_pid[pid] = {
            "title": row.get("Source_Title", ""),
            "doi": norm_doi(row.get("Source_DOI", "")),
            "year": row.get("Source_Year", ""),
            "authors": authors,
        }
    return by_pid


def request_json(url: str, params: dict[str, Any] | None = None, tries: int = 3) -> Any:
    for attempt in range(tries):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=15)
            if r.status_code == 429:
                time.sleep(1.5 * (attempt + 1))
                continue
            if r.status_code >= 400:
                return None
            return r.json()
        except Exception:
            if attempt + 1 == tries:
                return None
            time.sleep(0.75 * (attempt + 1))
    return None


def fetch_openalex(row: dict[str, str]) -> tuple[dict[str, Any] | None, float, str]:
    doi = norm_doi(row.get("Citing_DOI"))
    title = clean(row.get("Citing_Title"))
    if doi:
        data = request_json(OPENALEX, {"filter": f"doi:{doi}", "per-page": 2, "mailto": MAILTO})
        results = data.get("results", []) if isinstance(data, dict) else []
        if results:
            rec = max(results, key=lambda x: title_score(title, x.get("display_name", "")))
            return rec, title_score(title, rec.get("display_name", "")), "OpenAlex DOI"
    data = request_json(OPENALEX, {"search": title, "per-page": 6, "mailto": MAILTO})
    results = data.get("results", []) if isinstance(data, dict) else []
    if not results:
        return None, 0.0, ""
    y = int(row["Citing_Year"]) if clean(row.get("Citing_Year")).isdigit() else None
    ranked = []
    for rec in results:
        score = title_score(title, rec.get("display_name", ""))
        ry = rec.get("publication_year")
        if y and ry:
            score -= min(abs(int(ry) - y), 5) * 2.0
        ranked.append((score, rec))
    ranked.sort(key=lambda x: x[0], reverse=True)
    if ranked[0][0] < 82:
        return None, ranked[0][0], "OpenAlex title search below threshold"
    return ranked[0][1], ranked[0][0], "OpenAlex title search"


def fetch_crossref(row: dict[str, str]) -> dict[str, Any] | None:
    doi = norm_doi(row.get("Citing_DOI"))
    if not doi:
        return None
    data = request_json(f"{CROSSREF}/{quote(doi, safe='')}")
    if not isinstance(data, dict):
        return None
    msg = data.get("message")
    return msg if isinstance(msg, dict) else None


def oa_authors(rec: dict[str, Any]) -> tuple[list[dict[str, str]], list[str], list[str]]:
    authors, insts, countries = [], set(), set()
    for au in rec.get("authorships") or []:
        a = au.get("author") or {}
        name = clean(a.get("display_name"))
        if name:
            authors.append({"name": name, "id": clean(a.get("id")), "orcid": norm_orcid(a.get("orcid"))})
        for inst in au.get("institutions") or []:
            if clean(inst.get("display_name")):
                insts.add(clean(inst.get("display_name")))
            if clean(inst.get("country_code")):
                countries.add(clean(inst.get("country_code")))
        for c in au.get("countries") or []:
            if clean(c):
                countries.add(clean(c))
    return authors, sorted(insts), sorted(countries)


def cr_authors(rec: dict[str, Any]) -> list[dict[str, str]]:
    out = []
    for a in rec.get("author") or []:
        name = " ".join(x for x in [clean(a.get("given")), clean(a.get("family"))] if x)
        if name:
            out.append({"name": name, "id": "", "orcid": norm_orcid(a.get("ORCID"))})
    return out


def classify(authors: list[dict[str, str]], source: dict[str, Any], full_list: bool) -> tuple[str, str, str, str]:
    src_authors = source.get("authors") or []
    src_ids = {clean(a.get("id")): clean(a.get("name")) for a in src_authors if clean(a.get("id"))}
    src_orcids = {norm_orcid(a.get("orcid")): clean(a.get("name")) for a in src_authors if norm_orcid(a.get("orcid"))}
    src_names = {norm_name(a.get("name")): clean(a.get("name")) for a in src_authors if norm_name(a.get("name"))}
    verified_overlap, name_overlap = [], []
    strict_self = False
    for a in authors:
        n = norm_name(a.get("name"))
        aid = clean(a.get("id"))
        oid = norm_orcid(a.get("orcid"))
        if n == HAOYU_NAME or oid == HAOYU_ORCID:
            strict_self = True
        if aid and aid in src_ids and norm_name(src_ids[aid]) != HAOYU_NAME:
            verified_overlap.append(src_ids[aid])
        elif oid and oid in src_orcids and oid != HAOYU_ORCID:
            verified_overlap.append(src_orcids[oid])
        elif n and n in src_names and n != HAOYU_NAME:
            name_overlap.append(src_names[n])
    if strict_self:
        return "Strict self-citation", "Verified", "Haoyu Wu is confirmed in the full citing-author list.", "Haoyu Wu"
    if verified_overlap:
        names = sorted(set(verified_overlap))
        return "Coauthor-network citation", "Verified", "Source coauthor overlap confirmed by OpenAlex author ID or ORCID: " + "; ".join(names), "; ".join(names)
    if name_overlap:
        names = sorted(set(name_overlap))
        return "Coauthor-network citation", "Probable", "Exact full-name overlap with source coauthor(s), but no shared persistent identifier: " + "; ".join(names), "; ".join(names)
    if authors and full_list:
        return "Fully independent", "Verified", "No source-author overlap found in a full OpenAlex/Crossref author list.", ""
    return "Author metadata incomplete", "Unresolved", "No sufficiently complete author list was retrieved.", ""


def audit_one(row: dict[str, str], sources: dict[str, Any]) -> dict[str, Any]:
    original_class = clean(row.get("Independence_Class") or row.get("Phase2_Independence_Class"))
    original_conf = clean(row.get("Independence_Confidence") or row.get("Phase2_Independence_Confidence"))
    if original_conf != "Unresolved":
        out = dict(row)
        out.update({
            "Phase3_Fast_Attempted": "No", "Phase3_Fast_Match_Source": "Not needed",
            "Phase3_Fast_Title_Score": "", "Phase3_Fast_OpenAlex_ID": "",
            "Phase3_Fast_Canonical_DOI": norm_doi(row.get("Citing_DOI")),
            "Phase3_Fast_Authors": clean(row.get("Enriched_Authors")) or clean(row.get("Citing_Authors_Raw")),
            "Phase3_Fast_Institutions": clean(row.get("Institutions_Available")),
            "Phase3_Fast_Countries": clean(row.get("Countries_Available")),
            "Final_Independence_Class": original_class,
            "Final_Independence_Confidence": original_conf,
            "Final_Independence_Note": clean(row.get("Independence_Note")),
            "Overlapping_Source_Authors": "", "Phase3_Fast_Changed": "No",
        })
        return out

    oa, score, match_source = fetch_openalex(row)
    authors, insts, countries = [], [], []
    full_list = False
    oa_id = ""
    canonical_doi = norm_doi(row.get("Citing_DOI"))
    if oa:
        authors, insts, countries = oa_authors(oa)
        full_list = bool(authors)
        oa_id = clean(oa.get("id"))
        canonical_doi = norm_doi(oa.get("doi")) or canonical_doi
    if not authors:
        cr = fetch_crossref(row)
        if cr:
            authors = cr_authors(cr)
            full_list = bool(authors)
            match_source = (match_source + "; Crossref DOI fallback").strip("; ")
    final_class, final_conf, note, overlaps = classify(authors, sources[row["Source_ID"]], full_list)
    out = dict(row)
    out.update({
        "Phase3_Fast_Attempted": "Yes", "Phase3_Fast_Match_Source": match_source,
        "Phase3_Fast_Title_Score": round(score, 1) if score else "",
        "Phase3_Fast_OpenAlex_ID": oa_id,
        "Phase3_Fast_Canonical_DOI": canonical_doi,
        "Phase3_Fast_Authors": "; ".join(a["name"] for a in authors),
        "Phase3_Fast_Institutions": "; ".join(insts),
        "Phase3_Fast_Countries": "; ".join(countries),
        "Final_Independence_Class": final_class,
        "Final_Independence_Confidence": final_conf,
        "Final_Independence_Note": note,
        "Overlapping_Source_Authors": overlaps,
        "Phase3_Fast_Changed": "Yes" if (final_class, final_conf) != (original_class, original_conf) else "No",
    })
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields, seen = [], set()
    for row in rows:
        for k in row:
            if k not in seen:
                seen.add(k); fields.append(k)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def main() -> None:
    rows, _ = decode_existing_seed()
    sources = extract_sources(decode_existing_script(), rows)
    print("Recovered source-author counts:", {pid: len(sources[pid]["authors"]) for pid in PIDS}, flush=True)
    unresolved_before = sum(clean(r.get("Independence_Confidence") or r.get("Phase2_Independence_Confidence")) == "Unresolved" for r in rows)
    print(f"Unresolved rows before fast audit: {unresolved_before}", flush=True)
    results: list[dict[str, Any] | None] = [None] * len(rows)
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(audit_one, row, sources): i for i, row in enumerate(rows)}
        for done, fut in enumerate(as_completed(futs), 1):
            i = futs[fut]
            try:
                results[i] = fut.result()
            except Exception as e:
                row = dict(rows[i])
                original_class = clean(row.get("Independence_Class") or row.get("Phase2_Independence_Class"))
                original_conf = clean(row.get("Independence_Confidence") or row.get("Phase2_Independence_Confidence"))
                row.update({
                    "Phase3_Fast_Attempted": "Yes", "Phase3_Fast_Match_Source": "Error",
                    "Phase3_Fast_Title_Score": "", "Phase3_Fast_OpenAlex_ID": "",
                    "Phase3_Fast_Canonical_DOI": norm_doi(row.get("Citing_DOI")),
                    "Phase3_Fast_Authors": "", "Phase3_Fast_Institutions": "", "Phase3_Fast_Countries": "",
                    "Final_Independence_Class": original_class,
                    "Final_Independence_Confidence": original_conf,
                    "Final_Independence_Note": f"Fast metadata audit error: {type(e).__name__}: {e}",
                    "Overlapping_Source_Authors": "", "Phase3_Fast_Changed": "No",
                })
                results[i] = row
            if done % 25 == 0:
                print(f"processed {done}/299", flush=True)

    final = [r for r in results if r is not None]
    unresolved = [r for r in final if r.get("Final_Independence_Confidence") == "Unresolved"]
    changes = [r for r in final if r.get("Phase3_Fast_Changed") == "Yes"]
    write_csv(OUTPUT / "phase3_metadata_fast_master.csv", final)
    write_csv(OUTPUT / "phase3_metadata_fast_changes.csv", changes)
    write_csv(OUTPUT / "phase3_metadata_fast_unresolved.csv", unresolved)
    counts = Counter((r.get("Final_Independence_Class", ""), r.get("Final_Independence_Confidence", "")) for r in final)
    summary = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "exact_google_scholar_edges": len(final),
        "unresolved_before": unresolved_before,
        "unresolved_after": len(unresolved),
        "rows_changed": len(changes),
        "resolved_in_this_run": unresolved_before - len(unresolved),
        "final_independence": {f"{k[0]} | {k[1]}": v for k, v in sorted(counts.items())},
        "methodology": [
            "Only Phase 2 unresolved rows were sent to OpenAlex/Crossref; settled rows were retained.",
            "Verified coauthor-network overlap requires an OpenAlex author ID or ORCID match; full-name-only overlap is probable.",
            "A full API author list with no source-author overlap is classified as verified independent.",
            "This is a research audit, not legal advice.",
        ],
    }
    (OUTPUT / "phase3_metadata_fast_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
