"""Stage 6: Transient rate comparison — Harvard vs POSS-I.

Compares the transient rate found on Harvard plates (from Stage 5)
against the POSS-I rate reported by VASCO.

POSS-I reference:
  - 936 fields, each ~6.6° × 6.6° = ~43.6 sq deg
  - 107,875 transients (Bruehl & Villarroel 2025)
  - Two exposures per field (red + blue), so 936 plate pairs

Usage:
    poetry run python src/06_rate_comparison.py
"""

import sys
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.database import get_harvard_transients
from utils.statistics import poisson_rate_test

RESULTS_DIR = Path(__file__).parent.parent / "data" / "results"
FITS_DIR = Path(__file__).parent.parent / "data" / "fits_cutouts"

# POSS-I parameters
POSSI_N_FIELDS = 936
POSSI_FIELD_AREA_SQDEG = 6.6 * 6.6  # ~43.6 sq deg
POSSI_N_TRANSIENTS = 107_875


def count_pairs() -> int:
    """Count plate pair directories that have 2 FITS files."""
    n = 0
    if FITS_DIR.exists():
        for d in FITS_DIR.iterdir():
            if d.is_dir() and "_" in d.name:
                fits_files = list(d.glob("*.fits"))
                if len(fits_files) >= 2:
                    n += 1
    return n


def estimate_pair_area_sqdeg(series: str) -> float:
    """Rough field-of-view area for deep DASCH series.

    These are approximate values based on DASCH documentation.
    """
    areas = {
        "mc": 25.0,   # Metcalf triplet, ~5° field
        "mf": 25.0,   # similar to mc
        "rb": 16.0,   # Ross-Baker, ~4° field
        "b":  25.0,   # Bache doublet, ~5° field
    }
    return areas.get(series, 20.0)


def main():
    transients = get_harvard_transients()
    n_pairs = count_pairs()

    print("Stage 6: Transient rate comparison — Harvard vs POSS-I")
    print(f"  Harvard transients detected: {len(transients)}")
    print(f"  Harvard plate pairs analyzed: {n_pairs}")

    if n_pairs == 0:
        print("  No plate pairs found. Run Stages 4–5 first.")
        return

    # Count by series for weighted area estimate
    series_counts = {}
    for t in transients:
        s = t.get("series", "unk")
        series_counts[s] = series_counts.get(s, 0) + 1

    mean_area = 20.0  # default
    if series_counts:
        total_area = sum(estimate_pair_area_sqdeg(s) * c for s, c in series_counts.items())
        mean_area = total_area / sum(series_counts.values())

    # Rates: transients per plate-pair per sq deg
    harvard_rate = len(transients) / n_pairs / mean_area if n_pairs > 0 else 0
    possi_rate = POSSI_N_TRANSIENTS / POSSI_N_FIELDS / POSSI_FIELD_AREA_SQDEG

    print(f"\n  Harvard rate: {harvard_rate:.4f} transients / pair / sq deg")
    print(f"  POSS-I rate:  {possi_rate:.4f} transients / pair / sq deg")

    # Expected Harvard count under POSS-I rate
    expected = possi_rate * n_pairs * mean_area
    result = poisson_rate_test(k_obs=len(transients), k_exp=expected)

    print(f"\n  Expected (if POSS-I rate): {expected:.1f}")
    print(f"  Observed: {len(transients)}")
    print(f"  Rate ratio: {result['rate_ratio']}")
    print(f"  p-value: {result['p_value']}")
    print(f"  Significant: {result['significant']}")

    output = {
        "harvard_transients": len(transients),
        "harvard_pairs": n_pairs,
        "harvard_rate_per_pair_sqdeg": round(harvard_rate, 6),
        "possi_rate_per_pair_sqdeg": round(possi_rate, 6),
        **result,
    }

    out_path = RESULTS_DIR / "rate_comparison.json"
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(output, indent=2, default=str))
    print(f"\n  Results saved to {out_path}")


if __name__ == "__main__":
    main()
