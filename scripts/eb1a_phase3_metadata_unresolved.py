#!/usr/bin/env python3
from __future__ import annotations

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

UNRESOLVED_EDGE_IDS = {
    "E006","E007","E011","E013","E015","E016","E018","E020","E025","E026",
    "E032","E033","E036","E039","E042","E048","E052","E054","E055","E057",
    "E058","E062","E063","E070","E071","E073","E074","E075","E082","E085",
    "E087","E091","E093","E095","E103","E104","E105","E107","E110","E112",
    "E114","E115","E117","E118","E120","E121","E122","E123","E125","E126",
    "E129","E132","E137","E148","E152","E156","E160","E189","E191","E192",
    "E194","E195","E218","E226","E243","E244","E246","E257","E259","E269",
    "E275","E282","E290",
}


def main() -> None:
    rows, _ = ns["decode_existing_seed"]()
    sources = ns["extract_sources"](ns["decode_existing_script"](), rows)
    targets = [r for r in rows if r.get("Edge_ID") in UNRESOLVED_EDGE_IDS]
    if len(targets) != 73:
        raise RuntimeError(f"Expected 73 unresolved targets, found {len(targets)}")
    print("Recovered source-author counts:", {pid: len(sources[pid]["authors"]) for pid in ns["PIDS"]}, flush=True)

    for row in targets:
        row["Independence_Class"] = "Author metadata incomplete"
        row["Independence_Confidence"] = "Unresolved"
        row["Independence_Note"] = "Phase 2 unresolved record: full API author audit required."

    results = [None] * len(targets)
    with ThreadPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(ns["audit_one"], row, sources): i for i, row in enumerate(targets)}
        for done, fut in enumerate(as_completed(futs), 1):
            i = futs[fut]
            try:
                results[i] = fut.result()
            except Exception as e:
                row = dict(targets[i])
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
            if done % 10 == 0:
                print(f"processed {done}/73", flush=True)

    final = [r for r in results if r is not None]
    unresolved = [r for r in final if r.get("Final_Independence_Confidence") == "Unresolved"]
    resolved = [r for r in final if r.get("Final_Independence_Confidence") != "Unresolved"]
    ns["write_csv"](OUTPUT / "phase3_metadata_fast_unresolved_targets.csv", final)
    ns["write_csv"](OUTPUT / "phase3_metadata_fast_resolved.csv", resolved)
    ns["write_csv"](OUTPUT / "phase3_metadata_fast_still_unresolved.csv", unresolved)

    counts = Counter((r.get("Final_Independence_Class", ""), r.get("Final_Independence_Confidence", "")) for r in final)
    summary = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "phase2_unresolved_targets": len(targets),
        "resolved_by_api": len(resolved),
        "still_unresolved": len(unresolved),
        "target_results": {f"{k[0]} | {k[1]}": v for k, v in sorted(counts.items())},
        "methodology": [
            "The 73 Phase 2 unresolved Edge_IDs were sent to OpenAlex with Crossref DOI fallback.",
            "Verified coauthor-network overlap requires an OpenAlex author ID or ORCID match; full-name-only overlap is probable.",
            "A full API author list with no source-author overlap is classified as verified independent.",
            "Target results are merged into the immutable 299-edge Phase 2 master by Edge_ID.",
            "This is a research audit, not legal advice.",
        ],
    }
    (OUTPUT / "phase3_metadata_fast_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
