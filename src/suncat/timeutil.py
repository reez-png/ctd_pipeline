"""Times as numbers, independent of how pandas stores them.

Why this exists: `pd.to_datetime(x).astype("int64") / 1e9` was used to turn times into seconds. That is only
right when the times are stored in NANOSECONDS. pandas 3 stores times parsed from text, or built from Python
datetimes, in MICROSECONDS (datetime64[us]), so the same expression returns values 1000 times too small: a 1 s
interval becomes 0.001 s. Anything that compares an interval with a threshold in seconds (a matching gap, a
clip, a time constant) is then silently wrong, while anything that only compares two such numbers with each
other still works, which is why the error can go unnoticed. These helpers work in any resolution.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

_EPOCH = pd.Timestamp("1970-01-01", tz="UTC")


def to_utc(values, fmt: str | None = "ISO8601") -> pd.Series:
    """Times as a tz-aware UTC Series. Text is parsed with format="ISO8601" by default, which reads mixed
    forms (with and without fractional seconds, with or without an offset) in one column; without a format,
    pandas guesses one from the first value and turns every value written the other way into NaT. Tz-naive
    times are taken as UTC. Unreadable values become NaT."""
    s = values if isinstance(values, pd.Series) else pd.Series(values)
    if pd.api.types.is_datetime64_any_dtype(s):
        return s.dt.tz_localize("UTC") if s.dt.tz is None else s.dt.tz_convert("UTC")
    return pd.to_datetime(s, utc=True, errors="coerce", format=fmt)


def epoch_seconds(values) -> np.ndarray:
    """Seconds since 1970-01-01 UTC as float64 (NaN for NaT), whatever the stored resolution (s, ms, us, ns)."""
    return (to_utc(values) - _EPOCH).dt.total_seconds().to_numpy(float)


def elapsed_seconds(values) -> np.ndarray:
    """Seconds since the first valid time, as float64 (NaN for NaT): smaller numbers than epoch_seconds, so
    differences keep full precision."""
    t = to_utc(values)
    t0 = t.dropna().min() if t.notna().any() else _EPOCH
    return (t - t0).dt.total_seconds().to_numpy(float)
