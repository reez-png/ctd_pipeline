"""Cast identity, positions and sun time, shared by the steps so they cannot drift apart.

Nothing here assumes a cruise: cast files are matched by their trailing number (any prefix, any zero padding), and
the daylight test uses local SOLAR time, UTC + longitude / 15 h, because non-photochemical quenching follows the
sun, not the ship's clock.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd

_TRAILING = re.compile(r"(\d+)\s*$")
POSITION_FILES = ("{cruise}_CTD_cast_summary.csv", "stations.csv", "{cruise}_stations.csv")


def cast_number(value) -> int | None:
    """The cast number in a cell: 7, "7", "7.0", "P45_06_CTD_07", "ZZ12-03_CTD_107" -> int; blank or NaN -> None."""
    if value is None:
        return None
    s = str(value).strip()
    if not s or s.lower() in ("nan", "none", "nat"):
        return None
    try:
        f = float(s)
        return int(f) if np.isfinite(f) and f == int(f) else None
    except ValueError:
        m = _TRAILING.search(s)
        return int(m.group(1)) if m else None


def match_cast(value, cast_ids: Iterable[str]) -> str | None:
    """The cast id among `cast_ids` that `value` names: the id itself, or the one id whose trailing number equals the
    cast number in `value`. None if there is no such id or more than one (ambiguous)."""
    ids = list(cast_ids)
    s = "" if value is None else str(value).strip()
    if s in ids:
        return s
    n = cast_number(value)
    if n is None:
        return None
    hits = [c for c in ids if cast_number(c) == n]
    return hits[0] if len(hits) == 1 else None


def cnv_start_time(path) -> pd.Timestamp | None:
    """The `# start_time` of an SBE .cnv header, as UTC (the CTD clock is taken to be UTC, as step02 does)."""
    try:
        with open(path, encoding="latin-1", errors="replace") as fh:
            for ln in fh:
                m = re.match(r"#\s*start_time\s*=\s*([A-Za-z]{3}\s+\d{1,2}\s+\d{4}\s+\d{2}:\d{2}:\d{2})", ln)
                if m:
                    return pd.to_datetime(m.group(1), format="%b %d %Y %H:%M:%S", utc=True)
                if ln.strip().upper().startswith("*END*"):
                    break
    except OSError:
        pass
    return None


def _num(v):
    v = pd.to_numeric(v, errors="coerce")
    return np.nan if pd.isna(v) or v == 0.0 else float(v)   # an exact 0 is a dropped fix; 0.6 N is a real station


def read_cast_positions(meta_dir, cruise_id: str, cast_ids: Iterable[str] | None = None) -> pd.DataFrame:
    """Start-of-cast position and UTC start time per cast, from the cruise's cast summary (or stations.csv).

    Index: cast_id. Columns: lat, lon (decimal degrees; an exact 0 is read as a dropped fix, NaN), start_utc
    (tz-aware UTC or NaT), source (the file). A `cast` number is matched to `cast_ids` by trailing number; with no
    `cast_ids` it becomes <cruise_id>_CTD_<NN>. Empty frame when no position file exists."""
    meta_dir = Path(meta_dir)
    path = next((meta_dir / n.format(cruise=cruise_id) for n in POSITION_FILES
                 if (meta_dir / n.format(cruise=cruise_id)).exists()), None)
    cols = ["lat", "lon", "start_utc", "source"]
    if path is None:
        return pd.DataFrame(columns=cols, index=pd.Index([], name="cast_id"))
    s = pd.read_csv(path)
    s.columns = [c.strip().lower() for c in s.columns]
    ids = list(cast_ids) if cast_ids is not None else None
    if "cast_id" in s.columns:
        raw = s["cast_id"]
    elif "cast" in s.columns:
        raw = s["cast"]
    else:
        return pd.DataFrame(columns=cols, index=pd.Index([], name="cast_id"))
    if ids is not None:
        cid = [match_cast(v, ids) for v in raw]
    else:
        cid = [(str(v).strip() if cast_number(v) is None or not re.fullmatch(r"[\d.]+", str(v).strip())
                else f"{cruise_id}_CTD_{cast_number(v):02d}") for v in raw]
    latc = next((c for c in ("latitude_start", "lat_deg", "lat", "latitude") if c in s.columns), None)
    lonc = next((c for c in ("longitude_start", "lon_deg", "lon", "longitude") if c in s.columns), None)
    tc = next((c for c in ("date_time_utc_start", "datetime_utc_start", "datetime_start", "date_start") if c in s.columns),
              next((c for c in s.columns if "start" in c and ("date" in c or "time" in c)), None))
    out = pd.DataFrame({
        "cast_id": cid,
        "lat": [_num(v) for v in s[latc]] if latc else np.nan,
        "lon": [_num(v) for v in s[lonc]] if lonc else np.nan,
        "start_utc": pd.to_datetime(s[tc], errors="coerce", utc=True) if tc else pd.NaT,
        "source": path.name,
    })
    out = out[out["cast_id"].notna()].drop_duplicates("cast_id").set_index("cast_id")
    return out[cols]


def solar_hour(start_utc, lon) -> float | None:
    """Local mean solar hour (0 to 24) of a UTC time at longitude `lon` (degrees east): UTC + lon / 15 h.
    None if either is missing. (The equation of time, up to about 16 min, is ignored: the daylight window is hours.)"""
    if start_utc is None or pd.isna(start_utc) or lon is None or not np.isfinite(lon):
        return None
    t = pd.Timestamp(start_utc) + pd.Timedelta(hours=float(lon) / 15.0)
    return t.hour + t.minute / 60.0 + t.second / 3600.0


def npq_hour(cast_id: str, positions: pd.DataFrame, ctd_files: Mapping[str, Path] | None = None,
             fallback_time=None) -> tuple[float | None, str]:
    """Hour for the NPQ daylight test, and which rule gave it.

    "solar": solar hour at the cast start (start from the cast summary, else the CTD file header) at the cast's
    longitude. "clock": the hour of `fallback_time` (e.g. a bottle's logged time) when the cast has no longitude
    or start. (None, "none") when there is nothing to go on."""
    lon = positions.at[cast_id, "lon"] if cast_id in positions.index else np.nan
    start = positions.at[cast_id, "start_utc"] if cast_id in positions.index else pd.NaT
    if (start is None or pd.isna(start)) and ctd_files and cast_id in ctd_files:
        start = cnv_start_time(ctd_files[cast_id])
    h = solar_hour(start, lon)
    if h is not None:
        return h, "solar"
    if fallback_time is not None and pd.notna(fallback_time):
        t = pd.Timestamp(fallback_time)
        return t.hour + t.minute / 60.0, "clock"
    return None, "none"
