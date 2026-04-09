"""Stage 1: Query DASCH plate coverage for each VASCO source.

For each source in the catalog, queries DASCH for all exposures covering
that sky position, then filters to plates within the 1949-1957 window.

Results stored in SQLite (data/pipeline.db), table: plate_coverage
Checkpoints every N queries; safe to restart — skips already-queried positions.

Usage:
    poetry run python src/01_query_plate_coverage.py [--catalog vetted|full|test]
"""

import sys
import yaml
import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent))
from utils.dasch_api import DASCHClient, parse_csv_response
from utils.database import init_db, coverage_already_queried, save_coverage

CONFIG_PATH = Path(__file__).parent.parent / "config.yaml"
CATALOG_DIR = Path(__file__).parent.parent / "data" / "vasco_catalog"


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def load_catalog(which: str) -> pd.DataFrame:
    paths = {
        "vetted": CATALOG_DIR / "vetted_5399.csv",
        "full": CATALOG_DIR / "full_107k.csv",
        "test": CATALOG_DIR / "test_200.csv",
    }
    if which not in paths:
        raise ValueError(f"Unknown catalog: {which}. Choose vetted, full, or test.")
    p = paths[which]
    if not p.exists():
        raise FileNotFoundError(f"Catalog not found: {p}\nRun src/00_fetch_vasco_catalog.py first.")
    df = pd.read_csv(p)
    print(f"Loaded {len(df)} sources from {p.name}")
    return df


def filter_to_window(plates: list, date_start: str, date_end: str) -> list:
    """Keep only plates whose obs date falls in the survey window."""
    in_window = []
    for p in plates:
        obs = p.get("expdate", "")
        if not obs:
            continue
        # expdate format: "1949-11-19T00:41:44.7Z" or "1950-03-12"
        date_part = obs[:10]
        if date_start <= date_part <= date_end:
            in_window.append(p)
    return in_window


def query_one(client: DASCHClient, vasco_id: str, ra: float, dec: float,
              date_start: str, date_end: str) -> tuple[list, int]:
    """Query coverage for one source. Returns (all_plates, n_in_window)."""
    raw = client.query_exposures(ra_deg=ra, dec_deg=dec)
    if not isinstance(raw, list) or len(raw) < 2:
        return [], 0
    plates = parse_csv_response(raw)
    # Add normalized plate_id in daschlab format: {series}{platenum:05d}
    for p in plates:
        try:
            p["plate_id"] = f"{p['series']}{int(p['platenum']):05d}"
        except (KeyError, ValueError):
            p["plate_id"] = ""
    in_window = filter_to_window(plates, date_start, date_end)
    return plates, len(in_window)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", default=None,
                        choices=["vetted", "full", "test"],
                        help="Which catalog to use (default: auto-detect)")
    args = parser.parse_args()

    cfg = load_config()
    date_start = cfg["pipeline"]["date_start"]
    date_end = cfg["pipeline"]["date_end"]
    checkpoint_n = cfg["pipeline"]["checkpoint_interval"]

    # Auto-detect catalog
    if args.catalog is None:
        if (CATALOG_DIR / "vetted_5399.csv").exists():
            catalog = "vetted"
        elif (CATALOG_DIR / "full_107k.csv").exists():
            catalog = "full"
        else:
            catalog = "test"
        print(f"Auto-selected catalog: {catalog}")
    else:
        catalog = args.catalog

    df = load_catalog(catalog)
    init_db()
    client = DASCHClient()

    # Stats tracking
    n_queried = 0
    n_skipped = 0
    n_with_coverage = 0
    n_errors = 0

    sources = list(df.itertuples(index=False))
    total = len(sources)

    print(f"\nQuerying DASCH plate coverage for {total} sources")
    print(f"Window: {date_start} – {date_end}")
    print(f"Checkpoint every {checkpoint_n} queries\n")

    with tqdm(total=total, unit="src") as pbar:
        for row in sources:
            vasco_id = str(row.source_id)
            ra = float(row.ra)
            dec = float(row.dec)

            # Skip near-polar sources: |Dec| > 88° causes API timeouts
            # because thousands of overlapping plates exceed the Lambda response limit
            if abs(dec) > 88.0:
                tqdm.write(f"  SKIP polar source {vasco_id} (Dec={dec:.1f}°)")
                save_coverage(vasco_id, ra, dec, [], 0)
                n_errors += 1
                pbar.update(1)
                pbar.set_postfix(queried=n_queried, skipped=n_skipped,
                                 coverage=n_with_coverage, errors=n_errors)
                continue

            if coverage_already_queried(vasco_id):
                n_skipped += 1
                pbar.update(1)
                pbar.set_postfix(queried=n_queried, skipped=n_skipped,
                                 coverage=n_with_coverage, errors=n_errors)
                continue

            try:
                plates, n_window = query_one(client, vasco_id, ra, dec, date_start, date_end)
                save_coverage(vasco_id, ra, dec, plates, n_window)
                n_queried += 1
                if n_window > 0:
                    n_with_coverage += 1
            except Exception as e:
                tqdm.write(f"  ERROR {vasco_id} ({ra:.4f},{dec:.4f}): {e}")
                n_errors += 1
                # Save empty coverage so we don't retry on resume
                save_coverage(vasco_id, ra, dec, [], 0)

            pbar.update(1)
            pbar.set_postfix(queried=n_queried, skipped=n_skipped,
                             coverage=n_with_coverage, errors=n_errors)

    print(f"\n=== Stage 1 Complete ===")
    print(f"Newly queried : {n_queried}")
    print(f"Skipped (done): {n_skipped}")
    print(f"With window coverage: {n_with_coverage} / {n_queried} "
          f"({100*n_with_coverage/max(n_queried,1):.1f}%)")
    print(f"Errors: {n_errors}")

    if n_queried > 0:
        pct = 100 * n_with_coverage / n_queried
        print(f"\nGo/No-Go gate: {pct:.1f}% of sources have 1949-1957 coverage")
        if pct < 20:
            print("WARNING: < 20% coverage — review project scope before Stage 2.")
        else:
            print("Coverage is sufficient. Proceed to Stage 2.")


if __name__ == "__main__":
    main()
