"""DASCH StarGlass API wrapper with rate limiting and retry logic.

All endpoints use POST with JSON payloads.
Base URL: https://api.starglass.cfa.harvard.edu/full  (authenticated)
"""

import os
import time
import yaml
import requests
from pathlib import Path
from dotenv import load_dotenv

CONFIG_PATH = Path(__file__).parent.parent.parent / "config.yaml"
ENV_PATH = Path(__file__).parent.parent.parent / ".env"
_API_ROOT = "https://api.starglass.cfa.harvard.edu"


def load_config():
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


class DASCHClient:
    def __init__(self):
        load_dotenv(ENV_PATH)
        cfg = load_config()
        self.api_key = os.environ.get("DASCHLAB_API_KEY", "")
        if not self.api_key:
            raise RuntimeError("DASCHLAB_API_KEY not set. Add it to vasco-dasch/.env")
        self.base_url = f"{_API_ROOT}/full"  # authenticated endpoint
        self.rate_limit = cfg["dasch"]["rate_limit_rps"]
        self.max_retries = cfg["dasch"]["max_retries"]
        self.backoff_base = cfg["dasch"]["backoff_base"]
        self.timeout = cfg["dasch"].get("request_timeout_sec", 120)
        self._last_request = 0.0

    def _headers(self):
        return {
            "x-api-key": self.api_key,
            "Accept": "application/json",
            "User-Agent": "vasco-dasch-pipeline",
        }

    def _throttle(self):
        elapsed = time.monotonic() - self._last_request
        gap = 1.0 / self.rate_limit
        if elapsed < gap:
            time.sleep(gap - elapsed)
        self._last_request = time.monotonic()

    def post(self, path: str, payload: dict) -> list:
        """POST to a DASCH API endpoint; returns the parsed JSON response."""
        url = f"{self.base_url}{path}"
        for attempt in range(self.max_retries):
            self._throttle()
            try:
                resp = requests.post(
                    url, json=payload, headers=self._headers(),
                    timeout=self.timeout, allow_redirects=False
                )
                if resp.status_code == 429:
                    wait = self.backoff_base ** (attempt + 1)
                    print(f"  Rate limited; sleeping {wait:.1f}s")
                    time.sleep(wait)
                    continue
                resp.raise_for_status()
                return resp.json()
            except requests.RequestException as e:
                wait = self.backoff_base ** (attempt + 1)
                print(f"  Request error ({e}); retry {attempt+1}/{self.max_retries} in {wait:.1f}s")
                time.sleep(wait)
        raise RuntimeError(f"DASCH API failed after {self.max_retries} retries: {url}")

    def query_exposures(self, ra_deg: float, dec_deg: float) -> list:
        """Query all exposures covering a sky position.

        Returns a CSV-like list of strings; first row is column names.
        """
        return self.post("/dasch/dr7/queryexps", {
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
        })

    def query_refcat(self, ra_deg: float, dec_deg: float,
                     radius_arcsec: float = 10.0, refcat: str = "apass") -> list:
        """Query reference catalog sources near a position."""
        return self.post("/dasch/dr7/querycat", {
            "refcat": refcat,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "radius_arcsec": radius_arcsec,
        })

    def get_lightcurve(self, gsc_bin_index: int, ref_number: int,
                       refcat: str = "apass") -> list:
        """Retrieve lightcurve by reference catalog identifiers."""
        return self.post("/dasch/dr7/lightcurve", {
            "refcat": refcat,
            "gsc_bin_index": int(gsc_bin_index),
            "ref_number": int(ref_number),
        })

    def get_mosaic_url(self, plate_id: str, binning: int = 16) -> str:
        """Get presigned S3 URL for a plate mosaic via mosaic_package API."""
        raw = self.post("/dasch/dr7/mosaic_package", {
            "plate_id": plate_id,
            "binning": binning,
        })
        if not isinstance(raw, dict) or "baseFitsUrl" not in raw:
            raise RuntimeError(f"Unexpected mosaic_package response for {plate_id}: {raw!r}")
        return raw["baseFitsUrl"]

    def download_mosaic(self, plate_id: str, binning: int = 16,
                        dest_path: Path = None) -> Path:
        """Download a plate FITS mosaic.

        plate_id format: '{series}{platenum:05d}'  e.g. 'ai11826' or 'a00123'
        binning: 1 (full res) or 16 (preview, ~1/256 data volume)
        """
        s3_url = self.get_mosaic_url(plate_id, binning)
        if dest_path is None:
            dest_path = Path(f"{plate_id}_{binning:02d}.fits")
        self._throttle()
        resp = requests.get(s3_url, stream=True, timeout=self.timeout)
        resp.raise_for_status()
        with open(dest_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1 << 20):
                f.write(chunk)
        return dest_path


def parse_csv_response(rows: list) -> list[dict]:
    """Parse daschlab CSV-style response (first row = headers) into list of dicts."""
    if not rows:
        return []
    headers = rows[0].split(",")
    result = []
    for row in rows[1:]:
        pieces = row.split(",")
        result.append(dict(zip(headers, pieces)))
    return result
