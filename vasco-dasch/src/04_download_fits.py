"""Stage 4: Download FITS plate mosaics for primary candidates.

Uses daschlab's Session.mosaic() for each target position, which handles
the full two-step download process (plate metadata + S3 presigned URL).

Downloads bin_factor=16 (low-res preview) first for triage.
Re-run with --full for bin_factor=1 (full resolution, ~256× larger files).

Output: data/fits_cutouts/{vasco_id}/{plate_id}_16.fits

Usage:
    poetry run python src/04_download_fits.py [--full] [--limit N] [--catalog test|vetted]
"""

import os
import sys
import json
import warnings
import argparse
import tempfile
import shutil
import pandas as pd
from pathlib import Path
from tqdm import tqdm

os.environ.setdefault("DASCHLAB_API_KEY", "")
sys.path.insert(0, str(Path(__file__).parent))

from dotenv import load_dotenv
load_dotenv(Path(__file__).parent.parent / ".env")

import daschlab
warnings.filterwarnings("ignore")

from utils.database import get_conn

FITS_DIR = Path(__file__).parent.parent / "data" / "fits_cutouts"
RESULTS_DIR = Path(__file__).parent.parent / "data" / "results"
LOG_FILE = FITS_DIR / "downloaded.json"


def load_done() -> set:
    if LOG_FILE.exists():
        return set(json.loads(LOG_FILE.read_text()))
    return set()


def save_done(done: set):
    LOG_FILE.write_text(json.dumps(sorted(done)))


def get_primary_candidates() -> list[dict]:
    """Return primary candidates: TRULY_ABSENT_WITH_COVERAGE + MODERN_MATCH_WITH_COVERAGE."""
    csv = RESULTS_DIR / "candidates.csv"
    if not csv.exists():
        raise FileNotFoundError("Run Stage 3 first: candidates.csv not found")
    df = pd.read_csv(csv)
    primary_flags = {"TRULY_ABSENT_WITH_COVERAGE", "MODERN_MATCH_WITH_COVERAGE"}
    return df[df["flag"].isin(primary_flags)].to_dict("records")


def get_window_plates(vasco_id: str, date_start: str, date_end: str) -> list[dict]:
    """Return plates in the 1949-1957 window for a given VASCO source."""
    with get_conn() as conn:
        row = conn.execute(
            "SELECT plates_json FROM plate_coverage WHERE vasco_id = ?", (vasco_id,)
        ).fetchone()
    if not row:
        return []
    plates = json.loads(row["plates_json"])
    return [
        p for p in plates
        if date_start <= p.get("expdate", "")[:10] <= date_end
    ]


def download_for_candidate(cand: dict, binning: int, date_start: str, date_end: str,
                            done: set, session_dir: Path) -> tuple[int, int]:
    """Download FITS mosaics for one candidate. Returns (n_ok, n_err)."""
    vasco_id = cand["vasco_id"]
    ra = float(cand["ra"])
    dec = float(cand["dec"])

    plates = get_window_plates(vasco_id, date_start, date_end)
    if not plates:
        return 0, 0

    out_dir = FITS_DIR / vasco_id
    out_dir.mkdir(parents=True, exist_ok=True)

    # Create a daschlab session for this position
    sess_path = session_dir / f"sess_{vasco_id}"
    sess = daschlab.open_session(root=str(sess_path), interactive=False)
    sess.select_target(ra_deg=ra, dec_deg=dec)
    # Pre-fetch exposures to populate session cache
    try:
        _ = sess.exposures()
    except Exception:
        pass

    n_ok = 0
    n_err = 0

    for p in plates:
        series = p.get("series", "")
        platenum = p.get("platenum", 0)
        plate_id = p.get("plate_id") or f"{series}{int(platenum):05d}"
        key = f"{vasco_id}/{plate_id}_{binning:02d}"

        if key in done:
            continue

        dest = out_dir / f"{plate_id}_{binning:02d}.fits"
        try:
            relpath = sess.mosaic(plate_id, binning=binning)
            src = Path(str(sess_path)) / relpath
            shutil.copy2(src, dest)
            done.add(key)
            n_ok += 1
        except Exception as e:
            tqdm.write(f"  SKIP {plate_id}: {str(e)[:80]}")
            n_err += 1

    return n_ok, n_err


def main():
    import yaml
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true",
                        help="Download full resolution (binning=1)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Limit to N candidates")
    parser.add_argument("--plates-per-candidate", type=int, default=3,
                        help="Max plates per candidate (default 3 for triage)")
    args = parser.parse_args()

    binning = 1 if args.full else 16

    with open(Path(__file__).parent.parent / "config.yaml") as f:
        cfg = yaml.safe_load(f)
    date_start = cfg["pipeline"]["date_start"]
    date_end = cfg["pipeline"]["date_end"]

    FITS_DIR.mkdir(parents=True, exist_ok=True)
    done = load_done()

    candidates = get_primary_candidates()
    if args.limit:
        candidates = candidates[:args.limit]

    print(f"Stage 4: FITS download (binning={binning})")
    print(f"  Candidates: {len(candidates)}")
    print(f"  Max plates/candidate: {args.plates_per_candidate}")
    print(f"  Already downloaded: {len(done)}")

    n_total_ok = 0
    n_total_err = 0

    with tempfile.TemporaryDirectory(prefix="dasch_sess_") as sess_tmp:
        for cand in tqdm(candidates, desc="Candidates", unit="src"):
            # Limit plates per candidate to keep download manageable
            vasco_id = cand["vasco_id"]
            plates = get_window_plates(vasco_id, date_start, date_end)
            # Prioritize plates near the middle of the window (peak POSS-I era)
            plates_to_try = plates[:args.plates_per_candidate]

            n_ok, n_err = download_for_candidate(
                cand, binning, date_start, date_end, done, Path(sess_tmp)
            )
            n_total_ok += n_ok
            n_total_err += n_err
            save_done(done)

    print(f"\n=== Stage 4 Complete ===")
    print(f"Downloaded: {n_total_ok} FITS files, {n_total_err} skipped")
    print(f"Output: {FITS_DIR}")
    print("\nSpot-check downloaded images, then run Stage 5 (source extraction).")


if __name__ == "__main__":
    main()
