"""Stage 7: Spatial correlation — Harvard transients vs VASCO positions.

Tests whether independently-detected Harvard transients cluster near
VASCO transient positions more than expected by chance.

Usage:
    poetry run python src/07_spatial_correlation.py [--radius-arcsec 10] [--n-mc 10000]
"""

import sys
import argparse
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.database import get_harvard_transients
from utils.statistics import spatial_correlation_mc

RESULTS_DIR = Path(__file__).parent.parent / "data" / "results"
CATALOG_DIR = Path(__file__).parent.parent / "data" / "vasco_catalog"


def load_vasco_positions() -> pd.DataFrame:
    for name in ["vetted_5399.csv", "full_107k.csv", "test_200.csv"]:
        p = CATALOG_DIR / name
        if p.exists():
            return pd.read_csv(p)
    raise FileNotFoundError("No VASCO catalog found")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius-arcsec", type=float, default=10.0)
    parser.add_argument("--n-mc", type=int, default=5000)
    args = parser.parse_args()

    transients = get_harvard_transients()
    vasco = load_vasco_positions()

    print(f"Spatial correlation test")
    print(f"  Harvard transients: {len(transients)}")
    print(f"  VASCO positions: {len(vasco)}")
    print(f"  Match radius: {args.radius_arcsec} arcsec")
    print(f"  Monte Carlo shuffles: {args.n_mc}")

    if len(transients) == 0:
        print("No Harvard transients. Cannot run test.")
        return

    import numpy as np
    test_ra = np.array([t["ra"] for t in transients])
    test_dec = np.array([t["dec"] for t in transients])

    result = spatial_correlation_mc(
        test_ra=test_ra,
        test_dec=test_dec,
        ref_ra=vasco["ra"].values,
        ref_dec=vasco["dec"].values,
        n_mc=args.n_mc,
        radius_arcsec=args.radius_arcsec,
    )

    print(f"\nResults:")
    print(f"  Observed matches: {result['n_obs_matches']}")
    print(f"  MC expected:      {result['mc_mean']:.1f} ± {result['mc_std']:.1f}")
    print(f"  z-score:          {result['z_score']:.3f}")
    print(f"  p-value:          {result['p_value']:.4f}")
    print(f"  Significant:      {result['significant']}")

    if result["significant"]:
        if result["z_score"] > 0:
            print("\nINTERPRETATION: Harvard transients cluster near VASCO positions")
            print("  → Consistent with a real astrophysical transient population")
        else:
            print("\nINTERPRETATION: Harvard transients ANTI-cluster near VASCO positions")
    else:
        print("\nINTERPRETATION: No significant spatial correlation found")
        print("  → Harvard data neither confirms nor refutes VASCO positions")


if __name__ == "__main__":
    main()
