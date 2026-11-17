"""Per-cruise SUNA clock offset: cruises/<CRUISE>/metadata/suna_time_offsets.csv.

The SUNA logs with its own internal clock, set by hand before each cruise, so its time stamps can be off from the
CTD's by a constant. step04_suna_time_lag measures it from when the two instruments enter and leave the water
(independent of the nitrate profile); tools/diagnose_suna_lag.py reports it with its interval and drift test.

The file states the correction once per cruise, like suna_calibration.json for the calibration:

    cast_id,offset_s,source
    P45_07_CTD_01,-16.0,water entry (step04): +16.1 s, 95% CI +15.5 to +17.9, 27 casts, no drift

offset_s is ADDED to the SUNA time, so it is minus the measured lag (SUNA stamps 16 s late -> offset -16).
step02_merge applies it whenever the file exists (TIME_OFFSET_MODE "auto"); step03 carries it into the frames;
step02b applies it to the depths of the uncorrected co-located product. Cruise files are not published (*.csv is
git-ignored everywhere).
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

FILE_NAME = "suna_time_offsets.csv"


def offsets_path(cruise_dir) -> Path:
    return Path(cruise_dir) / "metadata" / FILE_NAME


def load_offsets(cruise_dir) -> dict:
    """{cast_id: offset_s} from the cruise's file; {} when there is none. Raises ValueError on a malformed file
    (missing columns, a non-numeric or non-finite offset, a cast listed twice)."""
    path = offsets_path(cruise_dir)
    if not path.exists():
        return {}
    t = pd.read_csv(path, encoding="utf-8-sig")
    t.columns = [c.strip() for c in t.columns]
    if not {"cast_id", "offset_s"} <= set(t.columns):
        raise ValueError(f"{path} needs columns cast_id and offset_s; found {list(t.columns)}")
    t = t.dropna(subset=["cast_id"])
    cast = t["cast_id"].astype(str).str.strip()
    if cast.duplicated().any():
        raise ValueError(f"{path}: casts listed twice: {sorted(set(cast[cast.duplicated()]))}")
    off = pd.to_numeric(t["offset_s"], errors="coerce")
    bad = ~off.map(lambda v: isinstance(v, float) and math.isfinite(v))
    if bad.any():
        raise ValueError(f"{path}: offset_s not a finite number for {list(cast[bad])}")
    return dict(zip(cast, off.astype(float)))


def describe(cruise_dir) -> str:
    """One line for the run log: what the file says, or that there is none."""
    path = offsets_path(cruise_dir)
    o = load_offsets(cruise_dir)
    if not o:
        return f"no {FILE_NAME}: SUNA times used as logged"
    vals = sorted(set(o.values()))
    src = ""
    try:
        s = pd.read_csv(path, encoding="utf-8-sig")
        if "source" in s.columns and s["source"].notna().any():
            src = f" ({str(s['source'].dropna().iloc[0])})"
    except Exception:
        pass
    span = f"{vals[0]:+g} s" if len(vals) == 1 else f"{vals[0]:+g} to {vals[-1]:+g} s"
    return f"{FILE_NAME}: {len(o)} casts, offset {span} added to SUNA time{src}"


def write_offsets(cruise_dir, casts, offset_s: float, source: str, overwrite: bool = False) -> Path:
    """Write one offset for every cast. Refuses to replace an existing file unless overwrite; needs a source."""
    if not source or not str(source).strip():
        raise ValueError("a source is required: say where the offset comes from (e.g. step04 water entry, n, CI)")
    if not math.isfinite(float(offset_s)):
        raise ValueError(f"offset must be finite, got {offset_s!r}")
    casts = sorted({str(c).strip() for c in casts if str(c).strip()})
    if not casts:
        raise ValueError("no casts to write")
    path = offsets_path(cruise_dir)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists:\n{path.read_text(encoding='utf-8-sig').strip()[:400]}\nuse --force to replace it")
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cast_id": casts, "offset_s": float(offset_s), "source": str(source).strip()}
                 ).to_csv(path, index=False, encoding="utf-8")
    return path
