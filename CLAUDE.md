# VASCO × DASCH Cross-Match Pipeline

## Project Goal

Independently test VASCO transient claims by cross-matching the published catalog of
107,875 pre-Sputnik transients against Harvard's DASCH DR7 plate archive.

**Key question:** Do Harvard plates from 1949–1957 show transient phenomena at the
same sky coordinates as Palomar's POSS-I transients?

**Outcome:** Publishable result regardless of finding (confirmation, constraint, or refutation).

## Data Sources

- **VASCO catalog** — Bruehl & Villarroel (2025), DOI: 10.1038/s41598-025-21620-3
  - 107,875 sources with J2000 RA/Dec and UTC timestamps
  - Vetted 5,399-source subset from Solano et al. (2022) via Spanish VO
- **DASCH DR7** — Harvard plate archive, fully free/open-access
  - API: no auth for metadata; registered API key required for FITS downloads
  - Registration: starglass.cfa.harvard.edu

## Project Structure

```
vasco-dasch/
├── config.yaml                  # API keys, paths, parameters
├── requirements.txt
├── data/
│   ├── vasco_catalog/           # full_107k.csv, vetted_5399.csv
│   ├── dasch_metadata/          # plate coverage results
│   ├── lightcurves/             # retrieved DASCH lightcurves
│   ├── fits_cutouts/            # downloaded plate images
│   └── results/                 # analysis outputs
├── src/
│   ├── 00_fetch_vasco_catalog.py
│   ├── 01_query_plate_coverage.py
│   ├── 02_retrieve_lightcurves.py
│   ├── 03_filter_candidates.py
│   ├── 04_download_fits.py
│   ├── 05_source_extraction.py
│   ├── 06_rate_comparison.py
│   ├── 07_spatial_correlation.py
│   ├── 08_shadow_analysis.py
│   ├── 09_generate_figures.py
│   └── utils/
│       ├── dasch_api.py         # API wrapper with rate limiting
│       ├── coordinates.py       # astropy coordinate helpers
│       ├── database.py          # SQLite storage layer
│       └── statistics.py        # statistical test functions
├── notebooks/
│   ├── exploration.ipynb
│   ├── coverage_map.ipynb
│   └── results_summary.ipynb
└── tests/
    ├── test_api_connection.py
    ├── test_coordinates.py
    └── test_on_known_source.py
```

## Python Dependencies

Managed via poetry (`pyproject.toml`). Run `poetry install` to set up the venv.

```
daschlab        # official DASCH Python toolkit
astropy         # coordinates, FITS, tables
photutils       # source extraction, PSF fitting
requests        # API calls
pandas / numpy / scipy / matplotlib
python-dotenv   # loads DASCHLAB_API_KEY from .env
openpyxl        # reads the nuclear dataset Excel file
```

## DASCH API — Actual Endpoints (discovered during development)

**Base URL:** `https://api.starglass.cfa.harvard.edu/full` (authenticated)
**Method:** POST with JSON body (NOT GET params)
**Auth:** `x-api-key` header, loaded from `.env` via `DASCHLAB_API_KEY`
**Response:** CSV-like list of strings; row 0 = column headers

| Endpoint | Payload keys |
|---|---|
| `/dasch/dr7/queryexps` | `ra_deg`, `dec_deg` |
| `/dasch/dr7/querycat` | `refcat`, `ra_deg`, `dec_deg`, `radius_arcsec` |
| `/dasch/dr7/lightcurve` | `refcat`, `gsc_bin_index`, `ref_number` |

**Lightcurve queries need refcat IDs** — must call `querycat` first to get
`gsc_bin_index` + `ref_number`, then call `lightcurve` with those.

**Real-world performance (measured):**
- `queryexps`: 15–20 sec/position (large buffered Lambda response)
- Near-polar sources (|Dec| > 88°) time out — excluded automatically
- Config: 120s timeout, 3 retries, 0.5 req/sec, exponential backoff
- Observed coverage: ~96% of test positions have 1949–1957 Harvard plates

**Timing estimates:**
- Vetted 5,399 sources: ~24 hours for Stage 1 (run overnight)
- Full 107,875 sources: ~500 hours — requires batching across days

## Pipeline Stages

### Stage 0 — Setup & Validation
- Install dependencies, register StarGlass account, download VASCO catalog
- Test DASCH API against 5–10 known coordinates
- **Claude role:** write setup script and test harness
- **Human role:** register account, verify credentials

### Stage 1 — Plate Coverage Query (automated)
- For each VASCO (RA, Dec): POST `/dasch/dr7/queryexps`, filter to 1949-11-01–1957-10-04
- Store plate_id, obs_date, exposure_time, series, limiting_mag in SQLite
- Rate limiting: 0.5 req/sec with exponential backoff; 3 retries
- Checkpoint every 100 queries; resume skips already-queried positions
- Sources at |Dec| > 88° are skipped (API timeout due to plate density)
- **Actual timing:** vetted set ~24 hours, full set ~500 hours (multi-day batch)

### Stage 2 — Lightcurve Retrieval (automated)
- For each position with coverage: POST `/dasch/dr7/lightcurve` with APASS catalog
- Store raw lightcurves in SQLite; same rate limiting and checkpointing pattern

### Stage 3 — Candidate Filtering (automated)
- Classify each VASCO source: `SINGLE_DETECTION`, `MULTI_DETECTION`, or `NO_DETECTION`
- Apply quality cuts (reject plate-edge detections, bad photometry)
- Output: `candidate_list.csv`

### Stage 4 — FITS Download (automated, authenticated)
- Download `bin_factor=16` first for triage; `bin_factor=01` for interesting candidates
- Requires `x-api-key` header
- Estimated: 100–500 cutouts for single-detection candidates

### Stage 5 — Source Extraction & PSF Analysis
- `photutils` DAOStarFinder for detection; measure FWHM vs stellar median
- Flag anomalous PSF (narrower than stars = flash-like)
- **Human role:** visually inspect flagged candidates (critical step)

### Stage 6 — Statistical Analysis (automated)
1. **Transient rate comparison** — chi-square/Poisson test vs POSS-I rate
2. **Spatial correlation** — nearest-neighbor vs Monte Carlo shuffle
3. **Earth shadow deficit** — classify by umbral cone, compare to ~1.4% geometric expectation
4. **Nuclear test temporal correlation** — cross-ref with DOE/NV-209, contingency table

### Stage 7 — Figures & Paper Draft
- Sky coverage map, plate coverage histogram, lightcurve examples, PSF comparison plots
- **Claude role:** generate figure scripts, assist with methods section
- **Human role:** write interpretation, discussion, conclusions

## Key Engineering Patterns

All API query loops must include:
- Rate limiting at 1–2 req/sec
- Exponential backoff on errors
- Checkpoint/resume logic (skip already-queried records)
- SQLite as the storage layer (not flat files) for query-ability

## Go/No-Go Gate

After Stage 1: if fewer than 20% of VASCO positions have any DASCH coverage in
1949–1957, project scope needs rethinking. This is a 1-evening check.

## Scale & Storage

| Resource           | Vetted (5,399) | Full (107,875) |
|--------------------|----------------|----------------|
| API query time     | ~1.5 hours     | ~30 hours      |
| Metadata storage   | ~160 MB        | ~3 GB          |
| Lightcurve storage | ~500 MB        | ~10 GB         |
| FITS downloads     | ~20–60 GB      | ~200–1200 GB   |
| RAM needed         | 8 GB           | 16 GB ideal    |

## Target Journals

- **RNAAS** — 1–2 page result, fast turnaround, no page charges (good for null result)
- **PASP** — where Villarroel's aligned transients paper was published
- **MNRAS** — where the triple transient paper appeared
- **Scientific Reports** — where the nuclear test correlation paper appeared

## Standing Instructions for Claude

At the start of each conversation, review all files in `/notes/` and summarize any new research notes, TODOs, or decisions that are relevant to the current pipeline stage.

## Catalog Status

The VASCO positional catalog is **not auto-downloadable** — it must be obtained manually:
- Go to https://doi.org/10.1038/s41598-025-21620-3 → Supplementary Information
- Save RA/Dec table as `vasco-dasch/data/vasco_catalog/vetted_5399.csv`
- Required columns: `ra` (deg J2000), `dec` (deg J2000)
- A 200-source synthetic test catalog is at `data/vasco_catalog/test_200.csv`
- The nuclear test correlation dataset is at `data/vasco_catalog/raw_nuclear_dataset.xlsx`

## Running the Pipeline

```bash
cd vasco-dasch/
poetry install                          # first time only
poetry run python src/00_fetch_vasco_catalog.py   # validate catalog
./run_pipeline.sh --catalog vetted      # stages 0-3 overnight
# Then manually review candidates before:
poetry run python src/04_download_fits.py --limit 50
poetry run python src/05_source_extraction.py
poetry run python src/06_statistical_analysis.py
poetry run python src/09_generate_figures.py
```

## First Evening Checklist (✅ = already done)

- [x] Register at starglass.cfa.harvard.edu, get API key → in `.env`
- [x] `poetry install` — all dependencies installed
- [x] `poetry run python tests/test_api_connection.py` — PASSES
- [x] `poetry run python tests/test_coordinates.py` — PASSES
- [x] Stage 1 tested on 200-source synthetic catalog — 96% coverage
- [ ] Download VASCO supplementary data from Scientific Reports DOI above
- [ ] Join DASCH email list (dasch@gaggle.email) — optional but helpful
- [ ] Run `./run_pipeline.sh --catalog vetted` overnight
