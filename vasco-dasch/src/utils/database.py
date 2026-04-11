"""SQLite storage layer for the pipeline."""

import sqlite3
import json
from pathlib import Path
import yaml

CONFIG_PATH = Path(__file__).parent.parent.parent / "config.yaml"


def db_path() -> Path:
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)
    p = Path(__file__).parent.parent.parent / cfg["paths"]["db"]
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(db_path())
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    return conn


def init_db():
    with get_conn() as conn:
        conn.executescript("""
        CREATE TABLE IF NOT EXISTS plate_coverage (
            vasco_id    TEXT PRIMARY KEY,
            ra          REAL,
            dec         REAL,
            queried_at  TEXT,
            n_plates    INTEGER,
            n_window    INTEGER,
            plates_json TEXT
        );

        CREATE TABLE IF NOT EXISTS refcat_lookup (
            vasco_id        TEXT PRIMARY KEY,
            ra              REAL,
            dec             REAL,
            queried_at      TEXT,
            gsc_bin_index   INTEGER,
            ref_number      INTEGER,
            sep_arcsec      REAL,
            refcat          TEXT
        );

        CREATE TABLE IF NOT EXISTS lightcurves (
            vasco_id    TEXT PRIMARY KEY,
            ra          REAL,
            dec         REAL,
            queried_at  TEXT,
            lc_json     TEXT
        );

        CREATE TABLE IF NOT EXISTS candidates (
            vasco_id        TEXT PRIMARY KEY,
            ra              REAL,
            dec             REAL,
            flag            TEXT,
            n_detections    INTEGER,
            n_window        INTEGER,
            notes           TEXT
        );

        CREATE TABLE IF NOT EXISTS harvard_transients (
            transient_id    INTEGER PRIMARY KEY AUTOINCREMENT,
            ra              REAL,
            dec             REAL,
            plate_id        TEXT,
            series          TEXT,
            expdate         TEXT,
            mag             REAL,
            pair_id         TEXT,
            fwhm            REAL,
            snr             REAL
        );
        """)


def coverage_already_queried(vasco_id: str) -> bool:
    with get_conn() as conn:
        return conn.execute(
            "SELECT 1 FROM plate_coverage WHERE vasco_id = ?", (vasco_id,)
        ).fetchone() is not None


def save_coverage(vasco_id: str, ra: float, dec: float,
                  plates: list, n_window: int):
    with get_conn() as conn:
        conn.execute(
            """INSERT OR REPLACE INTO plate_coverage
               (vasco_id, ra, dec, queried_at, n_plates, n_window, plates_json)
               VALUES (?, ?, ?, datetime('now'), ?, ?, ?)""",
            (vasco_id, ra, dec, len(plates), n_window, json.dumps(plates)),
        )


def refcat_already_queried(vasco_id: str) -> bool:
    with get_conn() as conn:
        return conn.execute(
            "SELECT 1 FROM refcat_lookup WHERE vasco_id = ?", (vasco_id,)
        ).fetchone() is not None


def save_refcat(vasco_id: str, ra: float, dec: float,
                gsc_bin_index: int, ref_number: int,
                sep_arcsec: float, refcat: str):
    with get_conn() as conn:
        conn.execute(
            """INSERT OR REPLACE INTO refcat_lookup
               (vasco_id, ra, dec, queried_at, gsc_bin_index, ref_number, sep_arcsec, refcat)
               VALUES (?, ?, ?, datetime('now'), ?, ?, ?, ?)""",
            (vasco_id, ra, dec, gsc_bin_index, ref_number, sep_arcsec, refcat),
        )


def lightcurve_already_queried(vasco_id: str) -> bool:
    with get_conn() as conn:
        return conn.execute(
            "SELECT 1 FROM lightcurves WHERE vasco_id = ?", (vasco_id,)
        ).fetchone() is not None


def save_lightcurve(vasco_id: str, ra: float, dec: float, lc: list):
    with get_conn() as conn:
        conn.execute(
            """INSERT OR REPLACE INTO lightcurves
               (vasco_id, ra, dec, queried_at, lc_json)
               VALUES (?, ?, ?, datetime('now'), ?)""",
            (vasco_id, ra, dec, json.dumps(lc)),
        )


def get_positions_with_window_coverage() -> list:
    """Return (vasco_id, ra, dec) for positions with plates in the 1949-1957 window."""
    with get_conn() as conn:
        rows = conn.execute(
            "SELECT vasco_id, ra, dec FROM plate_coverage WHERE n_window > 0"
        ).fetchall()
    return [(r["vasco_id"], r["ra"], r["dec"]) for r in rows]


def save_harvard_transient(ra: float, dec: float, plate_id: str, series: str,
                           expdate: str, mag: float, pair_id: str,
                           fwhm: float, snr: float):
    with get_conn() as conn:
        conn.execute(
            """INSERT INTO harvard_transients
               (ra, dec, plate_id, series, expdate, mag, pair_id, fwhm, snr)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (ra, dec, plate_id, series, expdate, mag, pair_id, fwhm, snr),
        )


def save_harvard_transients_batch(records: list[dict]):
    with get_conn() as conn:
        conn.executemany(
            """INSERT INTO harvard_transients
               (ra, dec, plate_id, series, expdate, mag, pair_id, fwhm, snr)
               VALUES (:ra, :dec, :plate_id, :series, :expdate, :mag, :pair_id, :fwhm, :snr)""",
            records,
        )


def get_harvard_transients() -> list[dict]:
    with get_conn() as conn:
        rows = conn.execute("SELECT * FROM harvard_transients").fetchall()
    return [dict(r) for r in rows]


def clear_harvard_transients():
    with get_conn() as conn:
        conn.execute("DELETE FROM harvard_transients")


def get_refcat_for_lightcurve() -> list:
    """Return rows needed to query lightcurves: (vasco_id, ra, dec, gsc_bin_index, ref_number, refcat)."""
    with get_conn() as conn:
        rows = conn.execute("""
            SELECT r.vasco_id, r.ra, r.dec, r.gsc_bin_index, r.ref_number, r.refcat
            FROM refcat_lookup r
            JOIN plate_coverage p ON r.vasco_id = p.vasco_id
            WHERE p.n_window > 0
        """).fetchall()
    return [dict(r) for r in rows]
