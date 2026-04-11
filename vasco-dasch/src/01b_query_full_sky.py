"""Stage 1b: Query DASCH plate coverage across the full sky.

Queries a grid of positions to build a complete inventory of deep-series
plates in the 1949-1957 window, not limited to VASCO positions.

This gives an unbiased sample for rate comparison and shadow deficit tests.

Usage:
    poetry run python src/01b_query_full_sky.py [--step 5] [--workers 2]
"""

import sys
import argparse
import numpy as np
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

sys.path.insert(0, str(Path(__file__).parent))
from utils.dasch_api import DASCHClient, parse_csv_response
from utils.database import init_db, get_conn, save_coverage

DEEP_SERIES = {"mc", "mf", "rb", "b"}
WINDOW_START = "1949-11-01"
WINDOW_END = "1957-10-04"


def build_grid(step_deg: float) -> list[tuple[str, float, float]]:
    """Build a sky grid, spacing positions evenly in RA/Dec."""
    positions = []
    dec_steps = np.arange(-85, 86, step_deg)
    for dec in dec_steps:
        # Adjust RA spacing for cos(dec) to keep roughly uniform sky density
        cos_dec = max(np.cos(np.radians(dec)), 0.1)
        ra_step = step_deg / cos_dec
        ra_steps = np.arange(0, 360, ra_step)
        for ra in ra_steps:
            grid_id = f"grid_{ra:07.3f}_{dec:+07.3f}"
            positions.append((grid_id, float(ra), float(dec)))
    return positions


def query_one(client: DASCHClient, ra: float, dec: float) -> list[dict]:
    """Query exposures and filter to deep series in window."""
    raw = client.query_exposures(ra_deg=ra, dec_deg=dec)
    if not isinstance(raw, list) or len(raw) < 2:
        return []
    plates = parse_csv_response(raw)
    filtered = []
    for p in plates:
        series = p.get("series", "")
        expdate = p.get("expdate", "")[:10]
        if series in DEEP_SERIES and WINDOW_START <= expdate <= WINDOW_END:
            try:
                p["plate_id"] = f"{series}{int(p['platenum']):05d}"
            except (KeyError, ValueError):
                p["plate_id"] = ""
            filtered.append(p)
    return filtered


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", type=float, default=5.0,
                        help="Grid spacing in degrees (default 5)")
    parser.add_argument("--workers", type=int, default=2,
                        help="Parallel query threads (default 2)")
    args = parser.parse_args()

    init_db()
    grid = build_grid(args.step)
    print(f"Stage 1b: Full-sky plate inventory")
    print(f"  Grid step: {args.step}°")
    print(f"  Grid positions: {len(grid)}")
    print(f"  Workers: {args.workers}")

    # Skip already-queried grid positions
    with get_conn() as conn:
        done = {row["vasco_id"] for row in
                conn.execute("SELECT vasco_id FROM plate_coverage WHERE vasco_id LIKE 'grid_%'").fetchall()}

    to_query = [(gid, ra, dec) for gid, ra, dec in grid if gid not in done]
    print(f"  Already queried: {len(done)}")
    print(f"  Remaining: {len(to_query)}\n")

    if not to_query:
        print("All grid positions already queried.")
        return

    _thread_clients = {}
    _lock = threading.Lock()

    def get_client():
        tid = threading.current_thread().ident
        if tid not in _thread_clients:
            with _lock:
                if tid not in _thread_clients:
                    _thread_clients[tid] = DASCHClient()
        return _thread_clients[tid]

    def worker(item):
        gid, ra, dec = item
        c = get_client()
        try:
            plates = query_one(c, ra, dec)
            return (gid, ra, dec, plates, None)
        except Exception as e:
            return (gid, ra, dec, [], e)

    n_queried = 0
    n_errors = 0
    n_with_plates = 0
    batch_size = args.workers * 2

    with tqdm(total=len(to_query), unit="pos") as pbar:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            for batch_start in range(0, len(to_query), batch_size):
                batch = to_query[batch_start:batch_start + batch_size]
                futures = {pool.submit(worker, item): item for item in batch}
                for future in as_completed(futures):
                    gid, ra, dec, plates, err = future.result()
                    if err:
                        tqdm.write(f"  ERROR {gid}: {err}")
                        n_errors += 1
                        save_coverage(gid, ra, dec, [], 0)
                    else:
                        save_coverage(gid, ra, dec, plates, len(plates))
                        n_queried += 1
                        if len(plates) > 0:
                            n_with_plates += 1
                    pbar.update(1)
                    pbar.set_postfix(queried=n_queried, hits=n_with_plates,
                                     errors=n_errors)
                del futures

    print(f"\n=== Stage 1b Complete ===")
    print(f"Grid positions queried: {n_queried}")
    print(f"With deep window plates: {n_with_plates}")
    print(f"Errors: {n_errors}")


if __name__ == "__main__":
    main()
