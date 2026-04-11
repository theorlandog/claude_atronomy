"""Utility for building consecutive plate pairs from the same telescope series.

For transient detection: a source on plate A but absent on plate B
(same series, consecutive dates, overlapping field) is a transient candidate.
"""

import json
from collections import defaultdict
from datetime import datetime

from .database import get_conn

DEEP_SERIES = {"mc", "mf", "rb", "b"}
WINDOW_START = "1949-11-01"
WINDOW_END = "1957-10-04"
DEFAULT_MAX_GAP_DAYS = 30


def _gap_days(plate_a: dict, plate_b: dict) -> int:
    da = datetime.fromisoformat(plate_a["expdate"][:10])
    db = datetime.fromisoformat(plate_b["expdate"][:10])
    return abs((db - da).days)


def build_consecutive_pairs(vasco_id: str, max_gap_days: int = DEFAULT_MAX_GAP_DAYS) -> list[tuple[dict, dict]]:
    """Build consecutive plate pairs for one VASCO position.

    Filters to deep series (mc, mf, rb, b) within the 1949-1957 window,
    groups by series, sorts by expdate, and returns consecutive pairs.
    """
    with get_conn() as conn:
        row = conn.execute(
            "SELECT plates_json FROM plate_coverage WHERE vasco_id = ?",
            (vasco_id,),
        ).fetchone()
    if not row:
        return []

    plates = json.loads(row["plates_json"])

    # Filter to deep series + window dates
    filtered = []
    for p in plates:
        series = p.get("series", "")
        expdate = p.get("expdate", "")[:10]
        if series in DEEP_SERIES and WINDOW_START <= expdate <= WINDOW_END:
            filtered.append(p)

    # Group by series, sort by expdate, form consecutive pairs
    by_series = defaultdict(list)
    for p in filtered:
        by_series[p["series"]].append(p)

    pairs = []
    for series, group in by_series.items():
        group.sort(key=lambda p: p.get("expdate", ""))
        for i in range(len(group) - 1):
            if _gap_days(group[i], group[i + 1]) <= max_gap_days:
                pairs.append((group[i], group[i + 1]))

    return pairs


def get_all_unique_pairs(max_gap_days: int = DEFAULT_MAX_GAP_DAYS) -> list[tuple[str, dict, dict]]:
    """Get all unique consecutive plate pairs across all VASCO positions.

    Deduplicates by (plate_id_A, plate_id_B) since multiple positions
    may share the same plates. Returns (vasco_id, plate_A, plate_B).

    Iterates one row at a time to avoid loading all plates_json into memory.
    """
    conn = get_conn()
    cursor = conn.execute(
        "SELECT vasco_id FROM plate_coverage WHERE n_window > 0 OR n_plates > 0"
    )
    vasco_ids = [row["vasco_id"] for row in cursor]
    conn.close()

    seen = set()
    unique_pairs = []

    for vasco_id in vasco_ids:
        pairs = build_consecutive_pairs(vasco_id, max_gap_days=max_gap_days)
        for a, b in pairs:
            pid_a = a.get("plate_id") or f"{a['series']}{int(a['platenum']):05d}"
            pid_b = b.get("plate_id") or f"{b['series']}{int(b['platenum']):05d}"
            pair_key = (pid_a, pid_b)
            if pair_key not in seen:
                seen.add(pair_key)
                unique_pairs.append((vasco_id, a, b))

    return unique_pairs
