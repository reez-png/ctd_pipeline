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


def to_utc(values, fmt: str | None = "ISO8601", what: str | None = None) -> pd.Series:
    """Times as a tz-aware UTC Series. Tz-naive times are taken as UTC.

    Text is parsed with format="ISO8601" first, which reads mixed ISO forms in one column (with and without
    fractional seconds, with or without an offset). Without a format, pandas guesses one from the FIRST value and
    then either raises on, or (errors="coerce") silently turns into NaT, every value written the other way: a file
    whose times are "11:27:56+00:00" on some rows and "11:27:56.500000+00:00" on others loses rows. Values that
    are not ISO at all (e.g. "08/04/2026 11:27") get a second, per-value pass (format="mixed"). What is still
    unreadable becomes NaT; with `what` given, the count is printed, so a loss is never silent."""
    s = values if isinstance(values, pd.Series) else pd.Series(values)
    if pd.api.types.is_datetime64_any_dtype(s):
        return s.dt.tz_localize("UTC") if s.dt.tz is None else s.dt.tz_convert("UTC")
    t = pd.to_datetime(s, utc=True, errors="coerce", format=fmt)
    text = s.notna() & s.astype(str).str.strip().ne("") & ~s.astype(str).str.strip().str.lower().isin(["nan", "nat"])
    miss = t.isna() & text
    if fmt is not None and miss.any():
        t = t.copy()
        t[miss] = pd.to_datetime(s[miss], utc=True, errors="coerce", format="mixed")
        miss = t.isna() & text
    if what and miss.any():
        print(f"  WARNING: {int(miss.sum())} of {len(s)} {what} times unreadable (NaT), e.g. {s[miss].iloc[0]!r}")
    return t


def epoch_seconds(values) -> np.ndarray:
    """Seconds since 1970-01-01 UTC as float64 (NaN for NaT), whatever the stored resolution (s, ms, us, ns)."""
    return (to_utc(values) - _EPOCH).dt.total_seconds().to_numpy(float)


def elapsed_seconds(values) -> np.ndarray:
    """Seconds since the first valid time, as float64 (NaN for NaT): smaller numbers than epoch_seconds, so
    differences keep full precision."""
    t = to_utc(values)
    t0 = t.dropna().min() if t.notna().any() else _EPOCH
    return (t - t0).dt.total_seconds().to_numpy(float)
