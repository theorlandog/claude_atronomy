"""Stage 6: Statistical analysis of DASCH cross-match results.

Runs four tests:
  1. Transient rate comparison (Poisson): Harvard rate vs POSS-I rate
  2. Spatial correlation (Monte Carlo): Are Harvard transients near VASCO positions?
  3. Earth shadow deficit: fraction of transients in antisolar direction
  4. Nuclear test temporal correlation: chi-square contingency table

Output: data/results/statistical_tests.json

Usage:
    poetry run python src/06_statistical_analysis.py
"""

import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils.statistics import (
    poisson_rate_test,
    spatial_correlation_mc,
    earth_shadow_test,
    contingency_test,
)
from utils.database import get_conn

RESULTS_DIR = Path(__file__).parent.parent / "data" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_candidates() -> pd.DataFrame:
    p = RESULTS_DIR / "candidates.csv"
    if not p.exists():
        raise FileNotFoundError("Run Stage 3 first: candidates.csv not found")
    return pd.read_csv(p)


def load_nuclear_dataset() -> pd.DataFrame:
    """Load the Bruehl nuclear test correlation dataset if available."""
    p = Path(__file__).parent.parent / "data" / "vasco_catalog" / "raw_nuclear_dataset.xlsx"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_excel(p)
    if "Date" in df.columns:
        df["date"] = pd.to_datetime(df["Date"])
    return df


def test1_rate_comparison(df_cands: pd.DataFrame) -> dict:
    """Test 1: Compare Harvard single-detection rate to POSS-I expectation.

    POSS-I had ~170 transients per plate per 36 sq deg (from VASCO paper).
    We compare this to our observed Harvard rate.
    """
    print("Test 1: Transient rate comparison")

    n_single = int((df_cands["flag"] == "SINGLE_DETECTION").sum())
    n_sources_with_coverage = int((df_cands["flag"] != "NO_DETECTION").sum())

    # Get total Harvard plates in window
    with get_conn() as conn:
        rows = conn.execute(
            "SELECT SUM(n_window) as total_plates FROM plate_coverage"
        ).fetchone()
        total_plates = int(rows["total_plates"] or 0)

    if total_plates == 0:
        return {"error": "No plate coverage data"}

    # POSS-I reference: ~170 transients / plate / 36 sq deg
    # Rough estimate for Harvard plates (smaller fields, typically ~1 sq deg)
    poss_rate_per_plate = 170.0 / 36.0  # per sq deg per plate
    # Harvard plate field: assume ~1 sq deg average for A-series
    harvard_field_sqdeg = 1.0
    k_expected = poss_rate_per_plate * harvard_field_sqdeg * total_plates

    result = poisson_rate_test(n_single, k_expected)
    result.update({
        "test": "rate_comparison",
        "n_single_detections": n_single,
        "total_harvard_plates": total_plates,
        "poss_rate_per_plate_per_sqdeg": poss_rate_per_plate,
        "harvard_field_sqdeg": harvard_field_sqdeg,
    })
    print(f"  n_single={n_single}, expected={k_expected:.1f}, "
          f"ratio={result['rate_ratio']:.2f}, p={result['p_value']:.4f}")
    return result


def test2_spatial_correlation(df_cands: pd.DataFrame) -> dict:
    """Test 2: Are Harvard single detections spatially correlated with VASCO positions?"""
    print("Test 2: Spatial correlation")

    singles = df_cands[df_cands["flag"] == "SINGLE_DETECTION"]
    if len(singles) == 0:
        return {"error": "No single-detection candidates"}

    # All VASCO positions are the reference set
    all_ra = df_cands["ra"].values
    all_dec = df_cands["dec"].values

    result = spatial_correlation_mc(
        test_ra=singles["ra"].values,
        test_dec=singles["dec"].values,
        ref_ra=all_ra,
        ref_dec=all_dec,
        n_mc=1000,
        radius_arcsec=10.0,
    )
    result["test"] = "spatial_correlation"
    print(f"  matches={result['n_obs_matches']}, MC_mean={result['mc_mean']:.1f}, "
          f"z={result['z_score']:.2f}, p={result['p_value']:.4f}")
    return result


def test3_earth_shadow(df_cands: pd.DataFrame) -> dict:
    """Test 3: Earth shadow deficit analysis."""
    print("Test 3: Earth shadow deficit")

    singles = df_cands[df_cands["flag"] == "SINGLE_DETECTION"]
    if len(singles) == 0 or "utc_date" not in df_cands.columns:
        return {"error": "No single-detection candidates or no UTC dates"}

    singles_with_date = singles.dropna(subset=["utc_date"]) if "utc_date" in singles.columns else singles.iloc[0:0]
    if len(singles_with_date) == 0:
        return {"error": "No observation dates available for shadow test"}

    result = earth_shadow_test(
        obs_times_utc=singles_with_date["utc_date"].tolist(),
        ra_list=singles_with_date["ra"].values,
        dec_list=singles_with_date["dec"].values,
    )
    result["test"] = "earth_shadow"
    print(f"  n_shadow={result.get('n_shadow','?')}/{result.get('n_total','?')}, "
          f"frac={result.get('frac_shadow','?'):.3f} vs expected {result.get('expected_frac','?'):.3f}, "
          f"p={result.get('p_value','?')}")
    return result


def test4_nuclear_correlation(df_nuclear: pd.DataFrame) -> dict:
    """Test 4: Nuclear test temporal correlation (replicates Bruehl & Villarroel 2025)."""
    print("Test 4: Nuclear test temporal correlation")

    if df_nuclear.empty:
        return {"error": "Nuclear dataset not available — download raw_nuclear_dataset.xlsx"}

    if "Nuclear_Testing_YN" not in df_nuclear.columns:
        return {"error": "Expected column 'Nuclear_Testing_YN' not found"}

    on_test = df_nuclear[df_nuclear["Nuclear_Testing_YN"] == 1]
    off_test = df_nuclear[df_nuclear["Nuclear_Testing_YN"] == 0]

    n_on = int(on_test["Transient_Positive"].sum())
    n_off = int(off_test["Transient_Positive"].sum())

    result = contingency_test(
        n_transients_on_test_days=n_on,
        n_test_days=len(on_test),
        n_transients_off_test_days=n_off,
        n_non_test_days=len(off_test),
    )
    result["test"] = "nuclear_correlation"
    print(f"  rate_on={result['rate_on_test_days']:.4f}, "
          f"rate_off={result['rate_off_test_days']:.4f}, "
          f"ratio={result['rate_ratio']:.2f}, p={result['p_value']:.4f}")
    return result


def main():
    print("=== Stage 6: Statistical Analysis ===\n")

    df_cands = load_candidates()
    df_nuclear = load_nuclear_dataset()

    results = {}
    results["test1_rate"] = test1_rate_comparison(df_cands)
    print()
    results["test2_spatial"] = test2_spatial_correlation(df_cands)
    print()
    results["test3_shadow"] = test3_earth_shadow(df_cands)
    print()
    results["test4_nuclear"] = test4_nuclear_correlation(df_nuclear)

    out_path = RESULTS_DIR / "statistical_tests.json"
    # Convert numpy/python bools to plain Python types for JSON
    def _jsonify(obj):
        if isinstance(obj, dict):
            return {k: _jsonify(v) for k, v in obj.items()}
        if isinstance(obj, (bool,)):
            return bool(obj)
        if hasattr(obj, 'item'):  # numpy scalar
            return obj.item()
        return obj

    out_path.write_text(json.dumps(_jsonify(results), indent=2))

    print(f"\n=== Stage 6 Complete ===")
    print(f"Results written to {out_path}")
    print("\nSummary:")
    for name, res in results.items():
        if "error" in res:
            print(f"  {name}: SKIPPED ({res['error']})")
        elif "p_value" in res:
            sig = "SIGNIFICANT" if res.get("significant") else "not significant"
            print(f"  {name}: p={res['p_value']:.4f} ({sig})")


if __name__ == "__main__":
    main()
