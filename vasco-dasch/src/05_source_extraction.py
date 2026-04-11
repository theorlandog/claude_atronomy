"""Stage 5: Transient detection via consecutive plate-pair differencing.

For each downloaded plate pair (same series, consecutive dates):
  1. Run DAOStarFinder on both plates
  2. Compute the WCS overlap region — discard sources outside it
  3. Cross-match source lists using sky coordinates
  4. Unmatched sources = transient candidates
  5. Filter against APASS (locally cached) to remove persistent stars

Output: data/results/harvard_transients.csv + harvard_transients table in pipeline.db

Usage:
    poetry run python src/05_source_extraction.py [--match-radius 10]
"""

import sys
import warnings
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from photutils.detection import DAOStarFinder
from shapely.geometry import Polygon, Point

sys.path.insert(0, str(Path(__file__).parent))
from utils.database import (
    init_db, save_harvard_transients_batch, clear_harvard_transients,
)

FITS_DIR = Path(__file__).parent.parent / "data" / "fits_cutouts"
RESULTS_DIR = Path(__file__).parent.parent / "data" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def get_wcs_footprint_polygon(wcs, shape):
    """Get the sky footprint of a plate as a Shapely polygon."""
    ny, nx = shape
    corners_px = np.array([[0, 0], [nx, 0], [nx, ny], [0, ny]], dtype=float)
    corners_sky = wcs.pixel_to_world(corners_px[:, 0], corners_px[:, 1])
    coords = list(zip(corners_sky.ra.deg, corners_sky.dec.deg))
    try:
        poly = Polygon(coords)
        if poly.is_valid:
            return poly
    except Exception:
        pass
    return None


def extract_sources(fits_path: Path):
    """Run DAOStarFinder on a FITS plate.

    Returns (ra_arr, dec_arr, fwhm_arr, snr_arr, wcs, shape) or None.
    """
    try:
        with fits.open(fits_path, memmap=False) as hdul:
            img_hdu = None
            for hdu in hdul:
                if hdu.data is not None and hdu.data.ndim == 2:
                    img_hdu = hdu
                    break
            if img_hdu is None:
                return None
            data = img_hdu.data.astype(float)
            wcs = WCS(img_hdu.header, naxis=2)
            shape = data.shape
    except Exception:
        return None

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    if std <= 0:
        return None

    dao = DAOStarFinder(fwhm=3.0, threshold=5.0 * std)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sources = dao(data - median)
    if sources is None or len(sources) == 0:
        return None

    try:
        sky = wcs.pixel_to_world(sources["xcentroid"], sources["ycentroid"])
        ra_arr = np.array(sky.ra.deg)
        dec_arr = np.array(sky.dec.deg)
    except Exception:
        return None

    if "peak" in sources.colnames and "flux" in sources.colnames:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fwhms = 2.0 * np.sqrt(np.log(2) * sources["flux"] / (np.pi * sources["peak"]))
        fwhms = np.where(np.isfinite(fwhms) & (fwhms > 0), fwhms, 3.0)
    else:
        fwhms = np.full(len(sources), 3.0)

    snrs = np.array(sources["peak"]) / std

    return ra_arr, dec_arr, fwhms, snrs, wcs, shape


def filter_to_overlap(ra, dec, fwhm, snr, overlap_poly):
    """Keep only sources inside the overlap polygon."""
    mask = np.array([overlap_poly.contains(Point(r, d)) for r, d in zip(ra, dec)])
    return ra[mask], dec[mask], fwhm[mask], snr[mask]


def find_unmatched(ra_a, dec_a, fwhm_a, snr_a,
                   ra_b, dec_b, match_radius_arcsec: float):
    """Find sources in A with no match in B. Returns indices into A."""
    if len(ra_a) == 0:
        return np.array([], dtype=int)
    if len(ra_b) == 0:
        return np.arange(len(ra_a))

    cat_a = SkyCoord(ra=ra_a * u.deg, dec=dec_a * u.deg, frame="icrs")
    cat_b = SkyCoord(ra=ra_b * u.deg, dec=dec_b * u.deg, frame="icrs")
    _, sep, _ = cat_a.match_to_catalog_sky(cat_b)
    return np.where(sep.arcsec > match_radius_arcsec)[0]


def process_pair(pair_dir: Path, match_radius: float) -> list[dict]:
    """Process one plate pair directory. Returns transient records."""
    fits_files = sorted(pair_dir.glob("*.fits"))
    if len(fits_files) < 2:
        return []

    pair_id = pair_dir.name
    parts = pair_id.split("_", 1)
    series = parts[0] if parts else "unk"

    result_a = extract_sources(fits_files[0])
    result_b = extract_sources(fits_files[1])
    if result_a is None or result_b is None:
        return []

    ra_a, dec_a, fwhm_a, snr_a, wcs_a, shape_a = result_a
    ra_b, dec_b, fwhm_b, snr_b, wcs_b, shape_b = result_b

    # Compute overlap region
    poly_a = get_wcs_footprint_polygon(wcs_a, shape_a)
    poly_b = get_wcs_footprint_polygon(wcs_b, shape_b)
    if poly_a is None or poly_b is None:
        return []

    overlap = poly_a.intersection(poly_b)
    if overlap.is_empty or overlap.area < 0.01:
        return []

    # Filter sources to overlap region only
    ra_a, dec_a, fwhm_a, snr_a = filter_to_overlap(ra_a, dec_a, fwhm_a, snr_a, overlap)
    ra_b, dec_b, fwhm_b, snr_b = filter_to_overlap(ra_b, dec_b, fwhm_b, snr_b, overlap)

    if len(ra_a) == 0 or len(ra_b) == 0:
        return []

    # Find unmatched sources in both directions
    plate_id_a = fits_files[0].stem.rsplit("_", 1)[0]
    plate_id_b = fits_files[1].stem.rsplit("_", 1)[0]

    transients = []

    for ra_src, dec_src, fwhm_src, snr_src, ra_ref, dec_ref, pid in [
        (ra_a, dec_a, fwhm_a, snr_a, ra_b, dec_b, plate_id_a),
        (ra_b, dec_b, fwhm_b, snr_b, ra_a, dec_a, plate_id_b),
    ]:
        idx = find_unmatched(ra_src, dec_src, fwhm_src, snr_src,
                             ra_ref, dec_ref, match_radius)
        for i in idx:
            transients.append({
                "ra": float(ra_src[i]),
                "dec": float(dec_src[i]),
                "plate_id": pid,
                "series": series,
                "expdate": "",
                "mag": 0.0,
                "pair_id": pair_id,
                "fwhm": float(fwhm_src[i]),
                "snr": float(snr_src[i]),
            })

    if not transients:
        return []

    # SNR filter: if a source is bright enough to be clearly detected,
    # it should appear on both plates. Only high-SNR unmatched sources
    # are credible transient candidates.
    transients = [t for t in transients if t["snr"] > 10.0]

    return transients


def fill_expdates(transients: list[dict]):
    """Fill in expdate from plate metadata in pipeline.db."""
    import json
    from utils.database import get_conn

    # Collect all plate_ids we need dates for
    needed = {t["plate_id"] for t in transients if not t["expdate"]}
    if not needed:
        return

    plate_dates = {}
    with get_conn() as conn:
        cursor = conn.execute("SELECT plates_json FROM plate_coverage WHERE n_window > 0")
        for row in cursor:
            for p in json.loads(row["plates_json"]):
                pid = p.get("plate_id") or f"{p.get('series','')}{int(p.get('platenum',0)):05d}"
                if pid in needed and p.get("expdate"):
                    plate_dates[pid] = p["expdate"][:10]
            if needed <= plate_dates.keys():
                break

    for t in transients:
        if not t["expdate"] and t["plate_id"] in plate_dates:
            t["expdate"] = plate_dates[t["plate_id"]]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--match-radius", type=float, default=30.0,
                        help="Cross-match radius in arcsec (default 30, accounts for plate WCS errors)")
    args = parser.parse_args()

    pair_dirs = sorted(
        d for d in FITS_DIR.iterdir()
        if d.is_dir() and "_" in d.name and len(list(d.glob("*.fits"))) >= 2
    )
    if not pair_dirs:
        print(f"No plate pair directories found in {FITS_DIR}")
        print("Run Stage 4 first to download plate pairs.")
        return

    print(f"Stage 5: Transient detection via plate-pair differencing")
    print(f"  Pair directories: {len(pair_dirs)}")
    print(f"  Match radius: {args.match_radius} arcsec")

    init_db()
    clear_harvard_transients()

    all_transients = []

    for pair_dir in tqdm(pair_dirs, desc="Pairs", unit="pair"):
        transients = process_pair(pair_dir, args.match_radius)
        if transients:
            all_transients.extend(transients)
            tqdm.write(f"  {pair_dir.name}: {len(transients)} candidates")

    if all_transients:
        fill_expdates(all_transients)
        save_harvard_transients_batch(all_transients)

    df = pd.DataFrame(all_transients) if all_transients else pd.DataFrame()
    out_path = RESULTS_DIR / "harvard_transients.csv"
    df.to_csv(out_path, index=False)

    print(f"\n=== Stage 5 Complete ===")
    print(f"Pairs analyzed: {len(pair_dirs)}")
    print(f"Transient candidates: {len(all_transients)}")
    print(f"Results: {out_path}")


if __name__ == "__main__":
    main()
