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

```
daschlab        # official DASCH Python toolkit
astropy         # coordinates, FITS, tables
photutils       # source extraction, PSF fitting
requests        # API calls
pandas
numpy
scipy
matplotlib
```

## Pipeline Stages

### Stage 0 — Setup & Validation
- Install dependencies, register StarGlass account, download VASCO catalog
- Test DASCH API against 5–10 known coordinates
- **Claude role:** write setup script and test harness
- **Human role:** register account, verify credentials

### Stage 1 — Plate Coverage Query (automated)
- For each VASCO (RA, Dec): POST `/dasch/dr7/queryexps`, filter to 1949-11-01–1957-10-04
- Store plate_id, obs_date, exposure_time, series, limiting_mag in SQLite
- Rate limiting: 1–2 req/sec with exponential backoff
- Checkpoint every 100 queries; resume skips already-queried positions
- Vetted set: ~45 min. Full set: ~15 hours.

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

## First Evening Checklist

- [ ] `pip install daschlab astropy photutils requests pandas`
- [ ] Register at starglass.cfa.harvard.edu, get API key
- [ ] Join DASCH email list (dasch@gaggle.email)
- [ ] Download VASCO supplementary data from Scientific Reports DOI above
- [ ] Run: `python -c "import daschlab; print('ready')"`
- [ ] Test: query DASCH for a known bright star's lightcurve
- [ ] Test: query DASCH for a VASCO transient coordinate
- [ ] Verify: does the position have Harvard plate coverage 1949–1957?
