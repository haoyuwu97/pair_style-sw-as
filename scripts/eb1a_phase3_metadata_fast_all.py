#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import runpy
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "scripts" / "eb1a_phase3_metadata_fast.py"
ns = runpy.run_path(str(BASE_SCRIPT), run_name="phase3_fast_base")
OUTPUT = ns["OUTPUT"]


def main() -> None:
    rows, _ = ns["decode_existing_seed"]()
    sources = ns["extract_sources"](ns["decode_existing_script"](), rows)
    print("Recovered source-author counts:", {pid: len(sources[pid]["authors"]) for pid in ns["PIDS"]}, flush=True)

    # The compact bundle omits the Phase 2 independence columns. Force every exact
    # Scholar edge through the API resolver, then merge the result with Phase 2 by
    # Edge_ID when building the final workbook.
    for row in rows:
        row["Independence_Class"] = "Author metadata incomplete"
        row["Independence_Confidence"] = "Unresolved"
        row["Independence_Note"] = "Compact seed: full API author audit required."

    results = [None] * len(rows)
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(ns["audit_one"], row, sources): i for i, row in enumerate(rows)}
        for done, fut in enumerate(as_completed(futs), 1):
            i = futs[fut]
            try:
                results[i] = fut.result()
            except Exception as e:
                row = dict(rows[i])
                row.update({
                    "Phase3_Fast_Attempted": "Yes",
                    "Phase3_Fast_Match_Source": "Error",
                    "Phase3_Fast_Title_Score": "",
                    "Phase3_Fast_OpenAlex_ID": "",
                    "Phase3_Fast_Canonical_DOI": ns["norm_doi"](row.get("Citing_DOI")),
                    "Phase3_Fast_Authors": "",
                    "Phase3_Fast_Institutions": "",
                    "Phase3_Fast_Countries": "",
                    "Final_Independence_Class": "Author metadata incomplete",
                    "Final_Independence_Confidence": "Unresolved",
                    "Final_Independence_Note": f"Fast metadata audit error: {type(e).__name__}: {e}",
                    "Overlapping_Source_Authors": "",
                    "Phase3_Fast_Changed": "No",
                })
                results[i] = row
            if done % 25 == 0:
                print(f"processed {done}/299", flush=True)

    final = [r for r in results if r is not None]
    unresolved = [r for r in final if r.get("Final_Independence_Confidence") == "Unresolved"]
    changes = [r for r in final if r.get("Final_Independence_Confidence") != "Unresolved"]
    ns["write_csv"](OUTPUT / "phase3_metadata_fast_master.csv", final)
    ns["write_csv"](OUTPUT / "phase3_metadata_fast_changes.csv", changes)
    ns["write_csv"](OUTPUT / "phase3_metadata_fast_unresolved.csv", unresolved)

    counts = Counter((r.get("Final_Independence_Class", ""), r.get("Final_Independence_Confidence", "")) for r in final)
    summary = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "exact_google_scholar_edges": len(final),
        "api_audit_targets": len(rows),
        "unresolved_after": len(unresolved),
        "resolved_by_api": len(rows) - len(unresolved),
        "final_independence": {f"{k[0]} | {k[1]}": v for k, v in sorted(counts.items())},
        "methodology": [
            "All 299 exact Google Scholar/PoP edges were resolved independently through OpenAlex with Crossref fallback.",
            "Verified coauthor-network overlap requires an OpenAlex author ID or ORCID match; full-name-only overlap is probable.",
            "A full API author list with no source-author overlap is classified as verified independent.",
            "Results are merged with the Phase 2 master by Edge_ID; unresolved API rows retain the conservative Phase 2 status.",
            "This is a research audit, not legal advice.",
        ],
    }
    (OUTPUT / "phase3_metadata_fast_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
