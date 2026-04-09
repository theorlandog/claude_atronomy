"""Stage 5: Source extraction and PSF analysis on downloaded FITS plates.

For each FITS file in data/fits_cutouts/:
  - Run DAOStarFinder for source detection
  - Measure FWHM for all detected sources
  - Compare candidate transient FWHM to stellar median FWHM
  - Flag anomalous PSF (much narrower than stars = possible flash/artifact)

Output: data/results/psf_analysis.csv
  Columns: plate_id, vasco_id, cand_ra, cand_dec, cand_fwhm, median_fwhm,
            fwhm_ratio, n_stars_on_plate, cand_snr, cand_mag, anomalous_psf

Usage:
    poetry run python src/05_source_extraction.py [--fwhm-threshold 0.5]
"""

import sys
import json
import argparse
import warnings
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
from photutils.aperture import CircularAperture, aperture_photometry

sys.path.insert(0, str(Path(__file__).parent))
from utils.database import get_conn

FITS_DIR = Path(__file__).parent.parent / "data" / "fits_cutouts"
RESULTS_DIR = Path(__file__).parent.parent / "data" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_candidate_positions() -> dict:
    """Return {vasco_id: (ra, dec)} for primary candidates."""
    cands_csv = RESULTS_DIR / "candidates.csv"
    if not cands_csv.exists():
        return {}
    df = pd.read_csv(cands_csv)
    primary_flags = {"TRULY_ABSENT_WITH_COVERAGE", "MODERN_MATCH_WITH_COVERAGE"}
    primaries = df[df["flag"].isin(primary_flags)]
    return {row["vasco_id"]: (row["ra"], row["dec"]) for _, row in primaries.iterrows()}


def analyze_plate(fits_path: Path, cand_ra: float, cand_dec: float,
                  fwhm_threshold: float) -> dict | None:
    """Run source extraction on one FITS plate, return analysis dict."""
    try:
        with fits.open(fits_path, memmap=False) as hdul:
            # Find the image HDU
            img_hdu = None
            for hdu in hdul:
                if hdu.data is not None and hdu.data.ndim == 2:
                    img_hdu = hdu
                    break
            if img_hdu is None:
                return None

            data = img_hdu.data.astype(float)
            header = img_hdu.header
            wcs = WCS(header, naxis=2)
    except Exception as e:
        return {"error": str(e)}

    # Background stats
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)

    # Source detection
    fwhm_px = 3.0  # initial guess; typical for bin16 plates
    dao = DAOStarFinder(fwhm=fwhm_px, threshold=5.0 * std)
    sources = dao(data - median)
    if sources is None or len(sources) == 0:
        return None

    # Estimate FWHM from sharpness/roundness proxy (DAOStarFinder returns sharpness)
    # Use peak/flux ratio as FWHM proxy (narrower PSF → higher peak/flux)
    if "peak" in sources.colnames and "flux" in sources.colnames:
        fwhms = 2.0 * np.sqrt(np.log(2) * sources["flux"] / (np.pi * sources["peak"]))
        fwhms = fwhms[np.isfinite(fwhms) & (fwhms > 0)]
    else:
        fwhms = np.full(len(sources), fwhm_px)

    if len(fwhms) == 0:
        return None

    median_fwhm = float(np.median(fwhms))

    # Find source nearest to candidate position using WCS
    cand_sky = SkyCoord(ra=cand_ra * u.deg, dec=cand_dec * u.deg, frame="icrs")
    try:
        cand_x, cand_y = wcs.world_to_pixel(cand_sky)
        cand_x, cand_y = float(cand_x), float(cand_y)
    except Exception:
        return None

    # Find nearest detected source to candidate pixel position
    dx = sources["xcentroid"] - cand_x
    dy = sources["ycentroid"] - cand_y
    dist_px = np.sqrt(dx**2 + dy**2)
    nearest_idx = int(np.argmin(dist_px))
    nearest_dist = float(dist_px[nearest_idx])

    # Consider candidate detected if within 5 pixels of expected position
    MATCH_RADIUS_PX = 5.0
    if nearest_dist > MATCH_RADIUS_PX:
        return {
            "n_stars": len(sources),
            "median_fwhm": median_fwhm,
            "cand_detected": False,
            "cand_dist_px": nearest_dist,
        }

    cand_fwhm = float(fwhms[nearest_idx]) if nearest_idx < len(fwhms) else median_fwhm
    fwhm_ratio = cand_fwhm / max(median_fwhm, 0.01)
    cand_snr = float(sources["peak"][nearest_idx]) / std if std > 0 else 0.0

    return {
        "n_stars": len(sources),
        "median_fwhm": round(median_fwhm, 3),
        "cand_fwhm": round(cand_fwhm, 3),
        "fwhm_ratio": round(fwhm_ratio, 3),
        "cand_snr": round(cand_snr, 2),
        "cand_dist_px": round(nearest_dist, 2),
        "cand_detected": True,
        "anomalous_psf": bool(fwhm_ratio < fwhm_threshold),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fwhm-threshold", type=float, default=0.5,
                        help="Flag candidate if FWHM ratio < threshold (default 0.5)")
    args = parser.parse_args()

    fits_files = sorted(FITS_DIR.glob("**/*.fits"))
    if not fits_files:
        print(f"No FITS files found in {FITS_DIR}")
        print("Run Stage 4 first to download plate mosaics.")
        return

    cand_map = load_candidate_positions()
    print(f"Found {len(fits_files)} FITS files, {len(cand_map)} candidate positions")

    records = []
    n_anomalous = 0

    with tqdm(total=len(fits_files), unit="plate") as pbar:
        for fits_path in fits_files:
            # Directory name is vasco_id; strip _16 or _01 binning suffix from stem
            vasco_id = fits_path.parent.name
            plate_id = fits_path.stem.rsplit("_", 1)[0]  # "kf00038_16" → "kf00038"
            lookup = cand_map.get(vasco_id)
            if lookup is None:
                pbar.update(1)
                continue
            cand_ra, cand_dec = lookup

            result = analyze_plate(fits_path, cand_ra, cand_dec, args.fwhm_threshold)
            if result is None:
                pbar.update(1)
                continue

            rec = {
                "plate_id": plate_id,
                "vasco_id": vasco_id,
                "fits_file": fits_path.name,
                "cand_ra": cand_ra,
                "cand_dec": cand_dec,
                **result,
            }
            records.append(rec)
            if result.get("anomalous_psf"):
                n_anomalous += 1
                tqdm.write(f"  ANOMALOUS PSF: {plate_id} vasco={vasco_id} "
                           f"ratio={result.get('fwhm_ratio','?')}")
            pbar.update(1)

    df = pd.DataFrame(records)
    out_path = RESULTS_DIR / "psf_analysis.csv"
    df.to_csv(out_path, index=False)

    print(f"\n=== Stage 5 Complete ===")
    print(f"Plates analyzed: {len(records)}")
    print(f"Candidates detected on plate: {df['cand_detected'].sum() if 'cand_detected' in df else 0}")
    print(f"Anomalous PSF flagged: {n_anomalous}")
    print(f"Results: {out_path}")
    print("\nNEXT: Visually inspect anomalous PSF candidates before Stage 6.")


if __name__ == "__main__":
    main()
