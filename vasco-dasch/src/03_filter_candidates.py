"""Stage 3: Classify VASCO sources based on Harvard coverage + modern catalog check.

Classification flags:
  TRULY_ABSENT_WITH_COVERAGE  — no APASS counterpart + Harvard plates in 1949-1957
                                 → PRIMARY CANDIDATE: download FITS, run source extraction
  MODERN_MATCH_WITH_COVERAGE  — APASS source found + Harvard coverage
                                 → likely persistent star (VASCO false positive);
                                   worth checking Harvard plates to confirm persistence
  TRULY_ABSENT_NO_COVERAGE    — absent in modern catalogs but no Harvard coverage
                                 → interesting but can't cross-check with Harvard
  MODERN_MATCH_NO_COVERAGE    — APASS match + no Harvard coverage
                                 → lowest priority; probable false positive, can't verify
  NO_DATA                     — query failed or position not reached yet

Output: data/results/candidates.csv

Usage:
    poetry run python src/03_filter_candidates.py
"""

import sys
import json
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.database import get_conn, init_db

RESULTS_DIR = Path(__file__).parent.parent / "data" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def classify(n_window: int, has_modern_match: bool) -> str:
    if n_window > 0 and not has_modern_match:
        return "TRULY_ABSENT_WITH_COVERAGE"
    elif n_window > 0 and has_modern_match:
        return "MODERN_MATCH_WITH_COVERAGE"
    elif n_window == 0 and not has_modern_match:
        return "TRULY_ABSENT_NO_COVERAGE"
    else:
        return "MODERN_MATCH_NO_COVERAGE"


def main():
    init_db()

    with get_conn() as conn:
        coverage = {
            row["vasco_id"]: dict(row)
            for row in conn.execute(
                "SELECT vasco_id, ra, dec, n_plates, n_window FROM plate_coverage"
            ).fetchall()
        }
        refcat = {
            row["vasco_id"]: row["gsc_bin_index"]  # -1 = no match (TRULY_ABSENT)
            for row in conn.execute(
                "SELECT vasco_id, gsc_bin_index FROM refcat_lookup"
            ).fetchall()
        }

    print(f"Classifying {len(coverage)} sources...")
    print(f"  Refcat lookup available for: {len(refcat)} positions")

    records = []
    counts = {
        "TRULY_ABSENT_WITH_COVERAGE": 0,
        "MODERN_MATCH_WITH_COVERAGE": 0,
        "TRULY_ABSENT_NO_COVERAGE": 0,
        "MODERN_MATCH_NO_COVERAGE": 0,
        "NO_DATA": 0,
    }

    for vid, cov in coverage.items():
        if vid not in refcat:
            flag = "NO_DATA"
        else:
            has_modern = refcat[vid] != -1
            flag = classify(cov["n_window"], has_modern)

        counts[flag] += 1
        records.append({
            "vasco_id": vid,
            "ra": cov["ra"],
            "dec": cov["dec"],
            "flag": flag,
            "n_harvard_plates": cov["n_plates"],
            "n_window_plates": cov["n_window"],
            "has_modern_match": flag.startswith("MODERN_MATCH"),
        })

    df = pd.DataFrame(records)
    out = RESULTS_DIR / "candidates.csv"
    df.to_csv(out, index=False)

    print(f"\n=== Stage 3 Complete ===")
    print(f"Results: {out}")
    print(f"\nClassification summary:")
    total = len(records)
    for flag, count in counts.items():
        pct = 100 * count / max(total, 1)
        marker = " ← PRIMARY TARGETS" if flag == "TRULY_ABSENT_WITH_COVERAGE" else ""
        print(f"  {flag:<40s}: {count:5d} ({pct:.1f}%){marker}")

    n_primary = counts["TRULY_ABSENT_WITH_COVERAGE"]
    print(f"\n{n_primary} primary candidates for FITS inspection (Stage 4).")


if __name__ == "__main__":
    main()
