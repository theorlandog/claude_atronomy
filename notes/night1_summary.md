# VASCO × DASCH Project — Night 1 Summary
## April 8, 2026

---

## What we learned

### The VASCO debate
- The VASCO project (led by Beatriz Villarroel, Nordita/Stockholm) found 107,875
  transient sources on 1949–1957 Palomar Sky Survey plates — brief flashes that
  appear once and never again, all pre-dating Sputnik
- Key findings: 22σ deficit of transients inside Earth's shadow (consistent with
  sunlight reflections off orbiting objects at ~GEO altitude), 45% increased
  transient rate within 24 hours of nuclear tests
- Watters et al. (2026) published a 30-page critique arguing plate defects and
  observation schedule artifacts
- Villarroel rebutted (Feb 2026, arXiv:2602.15171)
- Independent replications: Doherty (April 2026) confirmed shadow deficit and
  nuclear correlation; Busko (March 2026) found similar transients on Hamburg
  Observatory plates via APPLAUSE — first independent telescope confirmation

### The gap we identified
- Nobody has cross-matched VASCO transients against Harvard's DASCH plate archive
- DASCH DR7 became fully public December 2024 — 430,000+ plates, 1885–1992
- The Hamburg APPLAUSE paper explicitly called for cross-instrument validation
- This is a novel, publishable project regardless of outcome

### What we ruled out
- Solar transit detection: ran the math, objects at GEO altitude produce
  0.09 ppm dimming against the solar disk — 5 orders of magnitude below
  photographic plate sensitivity. Dead end.
- Catching the identical flash on a Harvard plate: probability ~0%.
  Different approach needed (rate comparison, spatial correlation)

---

## What we accomplished tonight

### Live API queries against DASCH/StarGlass
- Confirmed the StarGlass search API works at:
  `POST https://api.starglass.cfa.harvard.edu/public/plates/search`
- Learned that "medium" precision search uses WCS bounds (can miss wide-field
  plates and produce false positives)
- "Low" precision search uses sky bins and returns comprehensive results

### Coverage confirmed
- Queried the nine-transient field (RA=83.5°, Dec=+28.9°) for all scanned
  Harvard plates 1949–1957
- Result: **~1,400 plates** potentially covering the field
- Multiple telescope series: ai, fa, ka, kb, kc, ke, kf, kg, b, bi, rh, and others

### Depth assessment
- Patrol cameras (ai, fa series): 1.5-inch aperture, wide field, 5–51 min
  exposures. Too shallow for VASCO transients (mag 15–16). Useful for rate
  comparison of brighter transients only.
- **8-inch Bache Doublet (b series): 30–45 min exposures. Comparable depth
  to Palomar's 48-inch Schmidt for the relevant magnitude range.**

### The 11 key plates identified
All 8-inch Bache Doublet, Bloemfontein, South Africa:

| Plate   | Date           | RA     | Dec    | Exp  | Offset from VASCO field |
|---------|----------------|--------|--------|------|------------------------|
| b74433  | 1949-02-03     | 78.7°  | +25.5° | 30m  | ~5.5°                  |
| b74449  | 1949-02-17     | 79.0°  | +25.7° | 45m  | ~5.3°                  |
| b74510  | 1949-03-22     | 85.5°  | +26.4° | 45m  | **~3.0°** ⭐           |
| b75044  | 1949-08-30     | 85.8°  | +26.4° | 45m  | **~3.0°** ⭐           |
| b75125  | 1949-10-21     | 79.4°  | +25.7° | 45m  | ~5.0°                  |
| b75567  | 1950-09-18     | 79.6°  | +25.7° | 45m  | ~5.0°                  |
| b75652  | 1950-11-12     | 79.7°  | +25.8° | 45m  | ~4.9°                  |
| b76004  | 1951-11-26     | 79.6°  | +25.6° | 45m  | ~5.1°                  |
| b76115  | 1952-03-18     | 85.5°  | +26.3° | 45m  | **~3.0°** ⭐           |
| b76610  | 1953-02-05     | 79.6°  | +25.6° | 45m  | ~5.0°                  |
| b76664  | 1953-03-08     | 79.4°  | +25.9° | 45m  | ~4.8°                  |

**Three best plates (⭐): b74510, b75044, b76115** — centered ~3° from the
nine-transient position, near plate center where image quality is best.

### Timeline proximity to VASCO events
- Nine-transient event: April 12, 1950
  - b74510: March 22, 1949 (13 months before)
  - b75044: August 30, 1949 (8 months before)
  - b75567: September 18, 1950 (5 months after)
- Triple transient event: July 19, 1952
  - b76115: March 18, 1952 (4 months before)

---

## Later progress (continued)

### Completed
- [x] Sign up at starglass.cfa.harvard.edu/signin
- [x] Get API key from user profile
- [x] Install daschlab via poetry (`pyproject.toml`)
- [x] Full 10-stage pipeline built and tested on 200-source synthetic catalog
- [x] VASCO vetted catalog (5,399 sources) downloaded from SVO ConeSearch
- [x] Stage 1 running on vetted catalog (~24 hour ETA)
- [x] Stage 4 test FITS downloads running (~32 hour ETA)
- [x] Stage 5 source extraction tested and bugs fixed

### VASCO catalog obtained
The vetted 5,399-source catalog was downloaded from the Spanish Virtual Observatory:
- **URL:** `http://svocats.cab.inta-csic.es/vanish-possi/`
- **Method:** ConeSearch VO protocol (`cs.php?RA=180&DEC=0&SR=180&VERB=2`)
- **Saved to:** `data/vasco_catalog/vetted_5399.csv`
- **Source:** Solano et al. (2022), MNRAS 515, 1380

The full 107,875-source catalog (Bruehl & Villarroel 2025) is **not publicly available** —
must be requested from the authors (stephen.bruehl@vumc.org or Beatriz Villarroel).

### DASCH API corrections (from Night 1 plan)
The plan had the API wrong. Actual API:
- **Base URL:** `https://api.starglass.cfa.harvard.edu/full`
- **Method:** POST with JSON body (NOT GET with query params)
- **Auth:** `x-api-key` header from `DASCHLAB_API_KEY` env var (set in `.env`)
- **Response:** CSV-like list of strings; row 0 = column headers
- **Performance:** 15–20 sec/position for `queryexps` (Lambda buffered response)
- Near-polar sources (|Dec| > 88°) time out — excluded automatically

### Science redesign: Stage 2
The original plan said "retrieve DASCH lightcurves for VASCO positions." This was wrong —
VASCO sources vanished and are NOT in modern APASS. Lightcurves require catalog IDs.

Stage 2 redesigned as **APASS modern catalog check**:
- APASS match → persistent star → VASCO false positive
- No APASS match → genuinely absent → genuine transient candidate

Classification (Stage 3):
- `TRULY_ABSENT_WITH_COVERAGE` — no APASS + Harvard plates → **PRIMARY TARGETS**
- `MODERN_MATCH_WITH_COVERAGE` — APASS match + Harvard plates → likely false positive
- `TRULY_ABSENT_NO_COVERAGE` / `MODERN_MATCH_NO_COVERAGE` — lower priority

### Stage 5 bugs fixed
1. `glob("*.fits")` missed files in subdirectories → `glob("**/*.fits")`
2. Filename stripping (`_bin16`) didn't match actual format (`_16`) → `rsplit("_", 1)[0]`
3. Filtered on `SINGLE_DETECTION` (a Stage 5 output) → now uses Stage 3 primary flags

### Pipeline test results (200-source synthetic catalog)
| Metric | Value |
|---|---|
| Harvard coverage (1949-1957) | **88%** of queried positions (59/200 so far) |
| Plates per position (median) | **~1,455** in the 8-year window |
| APASS matches | 0/48 with coverage (0%) |
| Primary candidates | **45/54 = 83%** TRULY_ABSENT_WITH_COVERAGE |
| Nuclear test correlation | rate ratio 1.37, p=0.20 (not significant alone) |

### Plate series depth analysis (critical finding)

**The b-series is NOT the best plate for the 1949–1957 window.**

Measured limiting magnitudes from DASCH `limMagApass` field for window plates:

| Series | Median lim mag | VASCO detectable? | Window plates (52 pos) |
|--------|---------------|-------------------|------------------------|
| mc (Metcalf) | 17.0 | **Yes — best** | 233 |
| mf | 16.4 | **Yes** | 389 |
| rb | 15.8 | **Marginal** | 585 |
| b (Bache) | 15.3 | **Marginal** | 325 |
| rh | 15.0 | **Barely** | 1,772 |
| ac | 14.1 | No | 7,419 |
| ai/fa/ka | 11.9–12.1 | No | ~33,000 each |

The b-series peaked in 1885–1920 and was winding down by 1949:

| Decade | b-series plates |
|--------|----------------|
| 1880s | 1,561 |
| 1890s | 4,116 |
| 1900s | 4,923 |
| 1910s | 2,867 |
| 1920s | 34 |
| 1930s | 1,338 |
| 1940s | 1,366 |
| 1950s | 308 |

**For FITS downloads, prioritize mc + mf + rb + b series.** These are the only
series with limiting magnitudes reaching VASCO transient brightness (mag 15–16).

### FITS download strategy

Disk: 293 GB free, 30 GB reserved = 263 GB usable.

Full vetted catalog at all window plates: **~20 TB** — impossible.
Deep plates only (mc+mf+rb+b): ~9.3 per position × 4,400 primaries × 3 MB = **~120 GB** — fits.
These plate counts are already position-filtered (queryexps returns only plates covering each RA/Dec).

Plan:
1. Let test-catalog Stage 4 finish (~196 GB, synthetic — useful for testing only)
2. Delete test FITS to free space
3. Run vetted Stage 4 filtered to mc+mf+rb+b series only (~120 GB)

### Processes running
- **Stage 1 (test):** PID 84697, ~59/200 sources, ~1 hour remaining
- **Stage 1 (vetted):** PID 108198, 5,399 sources, ~24 hours
- **Stage 4 (test):** PID 98367, downloading FITS, ~32 hours

---

## Next steps

### After Stage 1 (vetted) completes (~April 10)
- [ ] Run Stages 2→3 on vetted catalog to classify all 5,399 sources
- [ ] Delete test FITS data to free ~196 GB
- [ ] Add series filter to Stage 4 (mc+mf+rb+b only)
- [ ] Run Stage 4 on vetted primaries (~120 GB, ~14 hours)
- [ ] Run Stage 5 PSF analysis on deep plates
- [ ] Visual inspection of flagged candidates (critical human step)
- [ ] Re-run Stage 6 with real detections

### Remaining code work
- Stage 4 needs series filter (currently downloads all window plates)
- Statistical tests 2 (spatial) and 3 (shadow) need Stage 5 detections
- Figure 4 (lightcurve examples) not yet implemented

### Longer term
- [ ] Email Bruehl/Villarroel for full 107k catalog
- [ ] Join DASCH mailing list (dasch@gaggle.email)

---

## Key URLs
- StarGlass: https://starglass.cfa.harvard.edu
- StarGlass API docs: https://starglass.cfa.harvard.edu/docs/api/
- DASCH DR7: https://dasch.cfa.harvard.edu/dr7/
- daschlab docs: https://daschlab.readthedocs.io
- VASCO vetted catalog: http://svocats.cab.inta-csic.es/vanish-possi/
- VASCO paper: DOI 10.1038/s41598-025-21620-3
- DASCH mailing list: dasch@gaggle.email
