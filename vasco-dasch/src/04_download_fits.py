"""Stage 4: Download FITS plate mosaics for consecutive plate pairs.

For each unique consecutive pair from a deep series (mc, mf, rb, b),
downloads both plates at bin_factor=16 for transient detection in Stage 5.

Output: data/fits_cutouts/{pair_id}/{plate_id}_16.fits

Usage:
    poetry run python src/04_download_fits.py [--full] [--limit N]
"""

import os
import sys
import json
import warnings
import argparse
import tempfile
import shutil
from pathlib import Path
from tqdm import tqdm

os.environ.setdefault("DASCHLAB_API_KEY", "")
sys.path.insert(0, str(Path(__file__).parent))

from dotenv import load_dotenv
load_dotenv(Path(__file__).parent.parent / ".env")

import daschlab
warnings.filterwarnings("ignore")

from utils.plate_pairs import get_all_unique_pairs

FITS_DIR = Path(__file__).parent.parent / "data" / "fits_cutouts"
LOG_FILE = FITS_DIR / "downloaded.json"


def load_done() -> set:
    if LOG_FILE.exists():
        return set(json.loads(LOG_FILE.read_text()))
    return set()


def save_done(done: set):
    LOG_FILE.write_text(json.dumps(sorted(done)))


def make_pair_id(plate_a: dict, plate_b: dict) -> str:
    series = plate_a.get("series", "unk")
    pid_a = plate_a.get("plate_id") or f"{plate_a['series']}{int(plate_a['platenum']):05d}"
    pid_b = plate_b.get("plate_id") or f"{plate_b['series']}{int(plate_b['platenum']):05d}"
    return f"{series}_{pid_a}_{pid_b}"


def download_plate(sess, plate: dict, binning: int, dest: Path,
                   sess_path: Path) -> bool:
    """Download one plate FITS. Returns True on success."""
    pid = plate.get("plate_id") or f"{plate['series']}{int(plate['platenum']):05d}"
    try:
        relpath = sess.mosaic(pid, binning=binning)
        src = Path(str(sess_path)) / relpath
        shutil.copy2(src, dest)
        return True
    except Exception as e:
        tqdm.write(f"  SKIP {pid}: {str(e)[:80]}")
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true",
                        help="Download full resolution (binning=1)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Limit to N pairs")
    parser.add_argument("--max-gap-days", type=int, default=30,
                        help="Max days between consecutive plates (default 30)")
    args = parser.parse_args()

    binning = 1 if args.full else 16
    FITS_DIR.mkdir(parents=True, exist_ok=True)
    done = load_done()

    print(f"Stage 4: Building consecutive plate pairs (max gap: {args.max_gap_days} days)...")
    all_pairs = get_all_unique_pairs(max_gap_days=args.max_gap_days)
    if args.limit:
        all_pairs = all_pairs[:args.limit]

    print(f"  Unique pairs: {len(all_pairs)}")
    print(f"  Binning: {binning}")
    print(f"  Already downloaded: {len(done)}")

    n_ok = 0
    n_err = 0

    with tempfile.TemporaryDirectory(prefix="dasch_sess_") as sess_tmp:
        for vasco_id, plate_a, plate_b in tqdm(all_pairs, desc="Pairs", unit="pair"):
            pair_id = make_pair_id(plate_a, plate_b)
            out_dir = FITS_DIR / pair_id
            out_dir.mkdir(parents=True, exist_ok=True)

            pid_a = plate_a.get("plate_id") or f"{plate_a['series']}{int(plate_a['platenum']):05d}"
            pid_b = plate_b.get("plate_id") or f"{plate_b['series']}{int(plate_b['platenum']):05d}"
            key_a = f"{pair_id}/{pid_a}_{binning:02d}"
            key_b = f"{pair_id}/{pid_b}_{binning:02d}"

            if key_a in done and key_b in done:
                continue

            # Create a daschlab session using the VASCO position
            ra = float(plate_a.get("ra", 0))
            dec = float(plate_a.get("dec", 0))
            sess_path = Path(sess_tmp) / f"sess_{pair_id}"
            sess = daschlab.open_session(root=str(sess_path), interactive=False)
            sess.select_target(ra_deg=ra, dec_deg=dec)
            try:
                _ = sess.exposures()
            except Exception:
                pass

            for plate, pid, key in [(plate_a, pid_a, key_a), (plate_b, pid_b, key_b)]:
                if key in done:
                    continue
                dest = out_dir / f"{pid}_{binning:02d}.fits"
                if download_plate(sess, plate, binning, dest, sess_path):
                    done.add(key)
                    n_ok += 1
                else:
                    n_err += 1

            save_done(done)

    print(f"\n=== Stage 4 Complete ===")
    print(f"Downloaded: {n_ok} FITS files, {n_err} skipped")
    print(f"Output: {FITS_DIR}")
    print("\nRun Stage 5 for transient detection on plate pairs.")


if __name__ == "__main__":
    main()
