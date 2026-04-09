"""Stage 4: Download FITS plate mosaics for SINGLE_DETECTION candidates.

Downloads low-res (bin_factor=16) first for triage.
For sources flagged as interesting after visual inspection, re-run with
--full to download bin_factor=01 (full resolution).

Requires authenticated API key (DASCHLAB_API_KEY in .env).

Output: data/fits_cutouts/{plate_id}_bin16.fits

Usage:
    poetry run python src/04_download_fits.py [--full] [--limit N]
"""

import sys
import json
import argparse
import sqlite3
import pandas as pd
from pathlib import Path
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from utils.dasch_api import DASCHClient
from utils.database import get_conn

CONFIG_PATH = Path(__file__).parent.parent / "config.yaml"
FITS_DIR = Path(__file__).parent.parent / "data" / "fits_cutouts"
RESULTS_DIR = Path(__file__).parent.parent / "data" / "results"

DOWNLOAD_LOG = FITS_DIR / "downloaded.json"


def load_downloaded() -> set:
    if DOWNLOAD_LOG.exists():
        return set(json.loads(DOWNLOAD_LOG.read_text()))
    return set()


def save_downloaded(done: set):
    DOWNLOAD_LOG.write_text(json.dumps(sorted(done)))


def get_plate_ids_for_candidates() -> list[tuple[str, str]]:
    """Return (vasco_id, plate_id) pairs for SINGLE_DETECTION candidates."""
    with get_conn() as conn:
        # Get candidates
        cands = conn.execute(
            "SELECT vasco_id FROM candidates WHERE flag = 'SINGLE_DETECTION'"
        ).fetchall()
        if not cands:
            # Fall back to candidates.csv
            candidates_csv = RESULTS_DIR / "candidates.csv"
            if candidates_csv.exists():
                df = pd.read_csv(candidates_csv)
                cand_ids = set(df[df["flag"] == "SINGLE_DETECTION"]["vasco_id"])
            else:
                return []
        else:
            cand_ids = {r["vasco_id"] for r in cands}

        # Get plate IDs from coverage
        pairs = []
        for row in conn.execute(
            "SELECT vasco_id, plates_json FROM plate_coverage"
        ).fetchall():
            if row["vasco_id"] not in cand_ids:
                continue
            plates = json.loads(row["plates_json"])
            for p in plates:
                plate_id = p.get("plate_id") or f"{p.get('series','')}{int(p.get('platenum',0)):05d}"
                if plate_id:
                    pairs.append((row["vasco_id"], plate_id))
    return pairs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true",
                        help="Download full resolution (bin_factor=1) instead of preview")
    parser.add_argument("--limit", type=int, default=None,
                        help="Maximum number of plates to download")
    args = parser.parse_args()

    bin_factor = 1 if args.full else 16
    FITS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Try to populate candidates table from CSV if empty
    candidates_csv = RESULTS_DIR / "candidates.csv"
    if candidates_csv.exists():
        df_cands = pd.read_csv(candidates_csv)
        with get_conn() as conn:
            for _, row in df_cands.iterrows():
                conn.execute(
                    """INSERT OR IGNORE INTO candidates (vasco_id, ra, dec, flag, n_detections, n_window)
                       VALUES (?, ?, ?, ?, ?, ?)""",
                    (row["vasco_id"], row["ra"], row["dec"],
                     row["flag"], row.get("n_detections", 0), row.get("n_window", 0))
                )

    pairs = get_plate_ids_for_candidates()
    if not pairs:
        print("No SINGLE_DETECTION candidates found. Run Stage 3 first.")
        return

    if args.limit:
        pairs = pairs[:args.limit]

    done = load_downloaded()
    to_download = [(v, p) for v, p in pairs if f"{p}_bin{bin_factor:02d}" not in done]

    print(f"FITS download: {len(to_download)} plates to fetch "
          f"(bin_factor={bin_factor}, {len(done)} already done)")

    client = DASCHClient()
    n_ok = 0
    n_err = 0

    with tqdm(total=len(to_download), unit="plate") as pbar:
        for vasco_id, plate_id in to_download:
            dest = FITS_DIR / f"{plate_id}_bin{bin_factor:02d}.fits"
            key = f"{plate_id}_bin{bin_factor:02d}"
            try:
                client.download_mosaic(plate_id, binning=bin_factor, dest_path=dest)
                done.add(key)
                n_ok += 1
            except Exception as e:
                tqdm.write(f"  ERROR {plate_id}: {e}")
                n_err += 1
            pbar.update(1)
            pbar.set_postfix(ok=n_ok, errors=n_err)
            # Save checkpoint every 10 downloads
            if (n_ok + n_err) % 10 == 0:
                save_downloaded(done)

    save_downloaded(done)
    print(f"\n=== Stage 4 Complete ===")
    print(f"Downloaded: {n_ok}, Errors: {n_err}")
    print(f"FITS files in: {FITS_DIR}")
    print("\nNEXT: Visually spot-check a sample, then run Stage 5 (source extraction).")


if __name__ == "__main__":
    main()
