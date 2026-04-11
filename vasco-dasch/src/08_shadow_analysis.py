"""Stage 8: Earth shadow deficit analysis on Harvard transients.

Tests whether Harvard transients (from Stage 5) show the same shadow
deficit as VASCO — i.e., fewer transients in the antisolar direction
(Earth's umbral shadow cone, ~1.4% of sky).

This is an independent test: different telescope, different plates,
different source extraction.

Usage:
    poetry run python src/08_shadow_analysis.py
"""

import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from astropy.time import Time
from astropy.coordinates import get_sun, SkyCoord
import astropy.units as u

sys.path.insert(0, str(Path(__file__).parent))
from utils.database import get_harvard_transients
from utils.statistics import earth_shadow_test

RESULTS_DIR = Path(__file__).parent.parent / "data" / "results"


def main():
    transients = get_harvard_transients()

    print(f"Earth shadow analysis (Harvard transients)")
    print(f"  Total transients: {len(transients)}")

    # Filter to those with valid expdates
    with_dates = [t for t in transients if t.get("expdate")]
    print(f"  With observation dates: {len(with_dates)}")

    if len(with_dates) == 0:
        print("  No transients with dates — cannot run shadow test.")
        return

    ra_list = np.array([t["ra"] for t in with_dates])
    dec_list = np.array([t["dec"] for t in with_dates])
    dates = [t["expdate"] for t in with_dates]

    # Compute solar elongations
    elongations = []
    for ra, dec, date_str in zip(ra_list, dec_list, dates):
        try:
            t = Time(date_str)
            sun = get_sun(t)
            src = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
            elongations.append(float(sun.separation(src).deg))
        except Exception:
            elongations.append(np.nan)

    elongations = np.array(elongations)
    valid = np.isfinite(elongations)
    print(f"  Valid elongations: {valid.sum()}")

    if valid.sum() == 0:
        print("  Insufficient data for shadow test.")
        return

    # Shadow cone threshold
    SHADOW_THRESHOLD = 179.73
    in_shadow = elongations[valid] > SHADOW_THRESHOLD
    n_shadow = int(in_shadow.sum())
    n_total = int(valid.sum())
    frac_obs = n_shadow / n_total
    frac_exp = 0.014

    print(f"\nShadow cone analysis:")
    print(f"  In shadow: {n_shadow} / {n_total} = {frac_obs:.4f}")
    print(f"  Expected (geometric): {frac_exp:.3f}")

    result = earth_shadow_test(
        obs_times_utc=list(np.array(dates)[valid]),
        ra_list=ra_list[valid],
        dec_list=dec_list[valid],
    )
    print(f"  Chi-square: {result.get('chi2', '?')}, p={result.get('p_value', '?')}")
    print(f"  Significant deficit: {result.get('significant', False)}")

    out = RESULTS_DIR / "shadow_analysis.json"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, default=str))

    # Plot
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(elongations[valid], bins=36, range=(0, 180),
            color="steelblue", edgecolor="white")
    ax.axvline(SHADOW_THRESHOLD, color="red", linestyle="--",
               label=f"Shadow threshold ({SHADOW_THRESHOLD}°)")
    ax.set_xlabel("Solar elongation (degrees)")
    ax.set_ylabel("Number of transients")
    ax.set_title("Solar Elongation Distribution — Harvard Transients")
    ax.legend()
    fig.tight_layout()

    fig_dir = RESULTS_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / "fig_shadow_elongation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"\nResults saved to {out}")
    print(f"Figure saved to {fig_dir}/fig_shadow_elongation.png")


if __name__ == "__main__":
    main()
