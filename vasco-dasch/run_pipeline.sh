#!/bin/bash
# VASCO × DASCH Pipeline — Master Run Script
# Runs all stages in sequence. Safe to re-run: each stage skips already-done work.
#
# Usage:
#   ./run_pipeline.sh              # auto-detect catalog (vetted > full > test)
#   ./run_pipeline.sh --catalog test
#   ./run_pipeline.sh --catalog vetted
#   ./run_pipeline.sh --stage 3   # start from a specific stage

set -e

CATALOG="auto"
START_STAGE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --catalog) CATALOG="$2"; shift 2 ;;
        --stage)   START_STAGE="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

echo "================================================"
echo " VASCO × DASCH Cross-Match Pipeline"
echo " Catalog: $CATALOG | Starting at stage: $START_STAGE"
echo "================================================"
echo ""

run_stage() {
    local stage=$1
    local desc=$2
    local cmd=$3

    if [ "$stage" -lt "$START_STAGE" ]; then
        echo "[SKIP] Stage $stage: $desc (before start stage)"
        return
    fi

    echo ""
    echo "--- Stage $stage: $desc ---"
    eval "$cmd"
    echo "--- Stage $stage done ---"
}

run_stage 0 "Fetch/validate VASCO catalog" \
    "poetry run python src/00_fetch_vasco_catalog.py"

if [ "$CATALOG" = "auto" ]; then
    CATALOG_FLAG=""
else
    CATALOG_FLAG="--catalog $CATALOG"
fi

run_stage 1 "Query DASCH plate coverage" \
    "poetry run python src/01_query_plate_coverage.py $CATALOG_FLAG --workers 2"

run_stage 2 "Retrieve lightcurves" \
    "poetry run python src/02_retrieve_lightcurves.py"

run_stage 3 "Filter and classify candidates" \
    "poetry run python src/03_filter_candidates.py"

echo ""
echo "=== Go/No-Go Check ==="
poetry run python -c "
import pandas as pd
from pathlib import Path
csv = Path('data/results/candidates.csv')
if csv.exists():
    df = pd.read_csv(csv)
    n_single = (df['flag'] == 'SINGLE_DETECTION').sum()
    n_total = len(df)
    print(f'Single-detection candidates: {n_single} / {n_total}')
    if n_single == 0:
        print('No candidates for FITS followup. Check lightcurve data.')
    else:
        print(f'Proceed to Stage 4 to download {n_single} FITS cutouts.')
else:
    print('candidates.csv not found — check Stage 3')
"

echo ""
echo "NOTE: Stages 4-8 are not run automatically."
echo "Review candidates, then run:"
echo "  poetry run python src/04_download_fits.py [--limit 50]   # download plate pairs"
echo "  poetry run python src/05_source_extraction.py            # detect transients via pair diff"
echo "  poetry run python src/06_rate_comparison.py              # Harvard vs POSS-I rate"
echo "  poetry run python src/07_spatial_correlation.py          # Harvard ↔ VASCO clustering"
echo "  poetry run python src/08_shadow_analysis.py              # Earth shadow deficit"
echo "  poetry run python src/09_generate_figures.py"
echo ""
echo "================================================"
echo " Pipeline stages 0-3 complete."
echo "================================================"
