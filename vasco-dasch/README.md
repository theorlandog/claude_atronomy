# VASCO × DASCH Cross-Match Pipeline

An independent test of the VASCO transient claims, cross-matching the published
catalog of pre-Sputnik flashes against Harvard's DASCH DR7 plate archive.

## Background and motivation

In 2025, Bruehl & Villarroel published a catalog of **107,875 transient sources**
detected on Palomar Sky Survey plates from **1949–1957** — brief flashes that
appear on a single plate and never again, all pre-dating the launch of Sputnik
in October 1957. The VASCO group reported two striking statistical signatures:

- A **22σ deficit** of transients inside Earth's umbral shadow, consistent with
  reflective objects in geosynchronous orbit before any human-made satellite
  could plausibly exist there.
- A **45% increase** in transient rate within 24 hours of atmospheric nuclear tests.

The claims are contested. Watters et al. (2026) argued the signals are explained
by plate defects and observation-schedule artifacts; Villarroel rebutted in
February 2026. Independent replications by Doherty (April 2026) and Busko (March
2026, using Hamburg plates via APPLAUSE) recovered both signatures on different
data, but **no one has yet cross-matched VASCO against Harvard's DASCH archive**
— the largest digitized photographic plate collection in the world (430,000+
plates, 1885–1992), made fully public in December 2024.

That's the gap this project fills. Whatever the result — confirmation,
constraint, or refutation — it is publishable, because it is the first
cross-instrument check of VASCO using a fully independent telescope archive.

We explicitly are **not** trying to catch the same flash on a Harvard plate
(probability ~0% given exposure timing). The science questions are statistical:

1. Do Harvard plates show a comparable transient rate at VASCO positions?
2. Is the Earth-shadow deficit reproduced in DASCH detections?
3. Does the nuclear-test temporal correlation hold up on a second telescope?
4. Do VASCO positions cluster spatially in a way consistent with reflective
   objects, or is the distribution consistent with noise?

## Approach

For each VASCO source, ask DASCH: *what plates from 1949–1957 cover this sky
position, and is anything detected there?* Sources with Harvard coverage but
no modern (APASS) catalog match at that location are the primary candidates —
they represent positions that were genuinely empty in modern data, so any
detection on a 1949–1957 Harvard plate is interesting.

The pipeline is built around the DASCH StarGlass API
(`api.starglass.cfa.harvard.edu/full`), with checkpointed SQLite storage and
resumable stages so multi-day runs can be interrupted safely.

## Data pipeline

```mermaid
flowchart TD
    A[VASCO catalog<br/>vetted_5399.csv] --> S0[Stage 0<br/>Validate catalog]
    S0 --> S1[Stage 1<br/>queryexps per position]
    S1 -->|plate metadata| DB[(pipeline.db<br/>SQLite)]
    S1 --> S2[Stage 2<br/>APASS modern-catalog check]
    S2 -->|match / no match| DB
    DB --> S3[Stage 3<br/>Classify candidates]
    S3 --> CAND[candidates.csv<br/>TRULY_ABSENT_WITH_COVERAGE]
    CAND --> GATE{Go / No-Go<br/>≥20% coverage?}
    GATE -->|yes| S4[Stage 4<br/>FITS cutout download<br/>mc+mf+rb+b series]
    GATE -->|no| STOP[Rethink scope]
    S4 --> FITS[(data/fits_cutouts/)]
    FITS --> S5[Stage 5<br/>photutils source extraction<br/>+ PSF analysis]
    S5 --> HUMAN[Human visual<br/>inspection of flags]
    HUMAN --> S6[Stage 6<br/>Statistical tests]
    S6 --> R1[Rate comparison<br/>vs POSS-I]
    S6 --> R2[Spatial correlation<br/>nearest-neighbor MC]
    S6 --> R3[Earth shadow<br/>deficit test]
    S6 --> R4[Nuclear test<br/>temporal correlation]
    R1 --> S9[Stage 9<br/>Figures + draft]
    R2 --> S9
    R3 --> S9
    R4 --> S9
```

### Stage interactions with the DASCH API

```mermaid
sequenceDiagram
    participant P as Pipeline worker
    participant DB as pipeline.db
    participant API as StarGlass API
    P->>DB: already queried?
    DB-->>P: no
    P->>API: POST /dasch/dr7/queryexps {ra, dec}
    Note right of API: ~15–20 s/position<br/>Lambda buffered
    API-->>P: CSV rows (plate_id, date, exp, lim_mag, ...)
    P->>P: filter to 1949-11-01 .. 1957-10-04
    P->>DB: insert plate rows + checkpoint
    P->>API: POST /dasch/dr7/querycat {refcat=APASS, ra, dec, r}
    API-->>P: APASS matches (or empty)
    P->>DB: insert match flag
```

Four worker threads run `queryexps` in parallel; the API response time (not
our 0.5 req/s rate limit) is the bottleneck, so concurrency gives a clean ~4x
speedup. Checkpoints land every 100 queries so a crash loses at most a minute
of work.

### Candidate classification

```mermaid
flowchart LR
    V[VASCO source] --> Q{Harvard<br/>coverage<br/>1949–1957?}
    Q -->|no| N1[NO_COVERAGE]
    Q -->|yes| M{APASS match<br/>at position?}
    M -->|yes| FP[MODERN_MATCH<br/>likely false positive]
    M -->|no| TGT[TRULY_ABSENT<br/>WITH_COVERAGE<br/><b>primary target</b>]
    TGT --> F4[Stage 4 FITS]
```

## Repository layout

```
vasco-dasch/
├── config.yaml                  # API config, paths, parameters
├── pyproject.toml               # poetry-managed dependencies
├── run_pipeline.sh              # master script, runs stages 0–3
├── data/
│   ├── vasco_catalog/           # input catalogs
│   ├── pipeline.db              # SQLite: plates, matches, candidates
│   ├── fits_cutouts/            # downloaded plate images
│   └── results/                 # candidate lists, figures
├── src/
│   ├── 00_fetch_vasco_catalog.py
│   ├── 01_query_plate_coverage.py     # Stage 1, parallelized
│   ├── 02_retrieve_lightcurves.py     # APASS check
│   ├── 03_filter_candidates.py
│   ├── 04_download_fits.py
│   ├── 05_source_extraction.py        # photutils + PSF
│   ├── 06_rate_comparison.py
│   ├── 07_spatial_correlation.py
│   ├── 08_shadow_analysis.py
│   ├── 09_generate_figures.py
│   └── utils/                          # API client, coords, db, stats
└── tests/
```

## Running the pipeline

```bash
poetry install                                  # first time only
poetry run python tests/test_api_connection.py  # smoke test

./run_pipeline.sh --catalog vetted              # stages 0–3, overnight
# review candidates, then:
poetry run python src/04_download_fits.py --limit 50
poetry run python src/05_source_extraction.py
poetry run python src/06_rate_comparison.py
poetry run python src/09_generate_figures.py
```

All long-running stages are resumable: re-running `run_pipeline.sh` skips
positions already in `pipeline.db`. Sources at `|Dec| > 88°` are auto-excluded
because the API times out at the poles.

## Plate-series depth note

Not all DASCH series are deep enough for VASCO transients (mag 15–16). Patrol
cameras (`ai`, `fa`, `ka`) bottom out at mag ~12 and are useful only for
brighter rate comparisons. The series with limiting magnitudes that actually
reach VASCO depths are `mc`, `mf`, `rb`, and `b`. Stage 4 filters to those.

## Status

Stage 1 is currently running on the vetted 5,399-source catalog with 4 parallel
workers; ETA ~5 hours. The full 107,875-source catalog is not publicly available
and must be requested from the authors.

## References

- Bruehl & Villarroel (2025), *Sci. Rep.*, [10.1038/s41598-025-21620-3](https://doi.org/10.1038/s41598-025-21620-3)
- Solano et al. (2022), *MNRAS* 515, 1380 — vetted 5,399 subset
- Watters et al. (2026) — critique
- Villarroel (Feb 2026), arXiv:2602.15171 — rebuttal
- DASCH DR7: https://dasch.cfa.harvard.edu/dr7/
- StarGlass: https://starglass.cfa.harvard.edu
