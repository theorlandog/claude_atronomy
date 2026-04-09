"""Stage 2: Modern catalog vetting for VASCO positions.

VASCO sources are objects that VANISHED — by definition they should NOT be in
modern reference catalogs like APASS. This stage checks each VASCO position
against APASS to categorize it:

  MODERN_MATCH     — APASS source within search radius
                     → likely a persistent star; VASCO detection may be false positive
  TRULY_ABSENT     — no APASS source found
                     → position is empty in modern sky; consistent with genuine transient
                     → HIGH PRIORITY for Harvard plate inspection (Stage 4+)

This information feeds directly into Stage 3 classification:
  TRULY_ABSENT + Harvard coverage → primary candidates for FITS download
  MODERN_MATCH → lower priority (but still worth checking Harvard plates)

Results stored in SQLite: refcat_lookup table.
Safe to resume: skips already-queried positions.

Usage:
    poetry run python src/02_retrieve_lightcurves.py
"""

import sys
import math
import yaml
from pathlib import Path
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from utils.dasch_api import DASCHClient, parse_csv_response
from utils.database import (
    init_db, get_positions_with_window_coverage,
    refcat_already_queried, save_refcat,
)

CONFIG_PATH = Path(__file__).parent.parent / "config.yaml"


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def angular_sep_arcsec(ra1, dec1, ra2, dec2) -> float:
    """Great-circle separation in arcseconds."""
    d2r = math.pi / 180
    dra = (ra2 - ra1) * d2r
    d1 = dec1 * d2r
    d2 = dec2 * d2r
    a = (math.sin(dra / 2) ** 2
         + math.cos(d1) * math.cos(d2) * math.sin((d2 - d1) / 2) ** 2)
    return 2 * math.asin(min(1, math.sqrt(a))) * (180 / math.pi) * 3600


def find_nearest_refcat(rows: list[dict], ra: float, dec: float) -> dict | None:
    best = None
    best_sep = float("inf")
    for r in rows:
        try:
            sep = angular_sep_arcsec(ra, dec, float(r["ra_deg"]), float(r["dec_deg"]))
            if sep < best_sep:
                best_sep = sep
                best = dict(r)
                best["_sep_arcsec"] = sep
        except (KeyError, ValueError):
            continue
    return best


def main():
    cfg = load_config()
    radius = cfg["pipeline"]["search_radius_arcsec"]
    refcat = cfg["pipeline"].get("refcat", "apass")

    init_db()
    client = DASCHClient()

    positions = get_positions_with_window_coverage()
    print(f"Stage 2: Modern catalog check (APASS) for {len(positions)} positions")
    print(f"Search radius: {radius} arcsec | Refcat: {refcat}")
    print(f"Science goal: TRULY_ABSENT = no APASS match = genuine transient candidate\n")

    n_modern = 0
    n_absent = 0
    n_errors = 0

    with tqdm(total=len(positions), desc="APASS lookup", unit="src") as pbar:
        for vasco_id, ra, dec in positions:
            if refcat_already_queried(vasco_id):
                pbar.update(1)
                continue
            try:
                raw = client.query_refcat(ra_deg=ra, dec_deg=dec,
                                          radius_arcsec=radius, refcat=refcat)
                rows = parse_csv_response(raw) if isinstance(raw, list) else []
                nearest = find_nearest_refcat(rows, ra, dec)
                if nearest:
                    save_refcat(
                        vasco_id, ra, dec,
                        gsc_bin_index=int(nearest["gsc_bin_index"]),
                        ref_number=int(nearest["ref_number"]),
                        sep_arcsec=nearest.get("_sep_arcsec", 0.0),
                        refcat=refcat,
                    )
                    n_modern += 1
                else:
                    # No modern match — save sentinel (-1) indicating TRULY_ABSENT
                    save_refcat(vasco_id, ra, dec, -1, -1, -1.0, refcat)
                    n_absent += 1
            except Exception as e:
                tqdm.write(f"  ERROR {vasco_id}: {e}")
                save_refcat(vasco_id, ra, dec, -1, -1, -1.0, refcat)
                n_errors += 1
            pbar.update(1)
            pbar.set_postfix(modern=n_modern, absent=n_absent, errors=n_errors)

    total = n_modern + n_absent
    print(f"\n=== Stage 2 Complete ===")
    print(f"MODERN_MATCH  (APASS found):  {n_modern} ({100*n_modern/max(total,1):.1f}%)")
    print(f"TRULY_ABSENT  (no APASS):     {n_absent} ({100*n_absent/max(total,1):.1f}%)")
    print(f"Errors: {n_errors}")
    print(f"\nTRULY_ABSENT sources with Harvard coverage are primary targets for Stage 4.")


if __name__ == "__main__":
    main()
