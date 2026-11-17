"""Per-cruise SUNA response flags: cruises/<CRUISE>/metadata/suna_response.csv.

step02b fits a first-order response time tau to every cast (both legs, no bottles). Most casts share one tau,
but a few can respond far more slowly (on P45_07 casts 01, 09, 10 and 14 at 20 to 45 s against a cruise median
near 6 s). Their down-cast profiles are smeared across the nitracline more than the rest, which the bottles
cannot detect. The file states which casts those are, once per cruise, like suna_calibration.json and
suna_time_offsets.csv:

    cast_id,tau_s,slow_response,source
    P45_07_CTD_01,44.4,True,"step02b tau per cast; slow = tau >= 3 x cruise median 5.9 s (17.7 s)"

tools/set_suna_slow_casts.py writes it from step02b's table. step08 copies tau_s and slow_response into the L3
nitrate (suna_tau_s, nitrate_slow_response) and step06 marks the slow casts on the nitrate sections. With no
file, the L3 columns are empty: not assessed, which is not the same as not slow. Cruise files are not
published (*.csv is git-ignored everywhere).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

FILE_NAME = "suna_response.csv"
FACTOR = 3.0            # slow = tau at least this many times the cruise median


def response_path(cruise_dir) -> Path:
    return Path(cruise_dir) / "metadata" / FILE_NAME


def _truthy(v) -> bool:
    return str(v).strip().lower() in ("true", "1", "1.0", "yes")


def load_response(cruise_dir) -> dict:
    """{cast_id: (tau_s, slow)} from the cruise's file; {} when there is none. tau_s is NaN where step02b gave
    none (an edge or a skipped cast). Raises ValueError on a malformed file."""
    path = response_path(cruise_dir)
    if not path.exists():
        return {}
    t = pd.read_csv(path, encoding="utf-8-sig")
    t.columns = [c.strip() for c in t.columns]
    if not {"cast_id", "slow_response"} <= set(t.columns):
        raise ValueError(f"{path} needs columns cast_id and slow_response; found {list(t.columns)}")
    t = t.dropna(subset=["cast_id"])
    cast = t["cast_id"].astype(str).str.strip()
    if cast.duplicated().any():
        raise ValueError(f"{path}: casts listed twice: {sorted(set(cast[cast.duplicated()]))}")
    bad = ~t["slow_response"].astype(str).str.strip().str.lower().isin(["true", "false", "1", "0", "1.0", "0.0"])
    if bad.any():
        raise ValueError(f"{path}: slow_response must be True or False; check {list(cast[bad])}")
    tau = pd.to_numeric(t["tau_s"], errors="coerce") if "tau_s" in t else pd.Series(np.nan, index=t.index)
    return {c: (float(x), _truthy(s)) for c, x, s in zip(cast, tau, t["slow_response"])}


def describe(cruise_dir) -> str:
    """One line for the run log."""
    path = response_path(cruise_dir)
    r = load_response(cruise_dir)
    if not r:
        return f"no {FILE_NAME}: response not assessed per cast (L3 suna_tau_s and nitrate_slow_response empty)"
    slow = sorted(c for c, (_, s) in r.items() if s)
    src = ""
    try:
        s = pd.read_csv(path, encoding="utf-8-sig")
        if "source" in s.columns and s["source"].notna().any():
            src = f" ({str(s['source'].dropna().iloc[0])})"
    except Exception:
        pass
    names = ", ".join(c.split("_")[-1] for c in slow) if slow else "none"
    return f"{FILE_NAME}: {len(r)} casts, {len(slow)} slow: {names}{src}"


def flag_slow(table: pd.DataFrame, cruise_median: float, factor: float = FACTOR) -> pd.DataFrame:
    """From step02b's per-cast table (cast, tau, edge): tau_s, and slow = tau >= factor x cruise_median, or tau at
    the grid's upper edge (slower than the grid could measure). A skipped cast gets tau NaN and slow False."""
    if not (math.isfinite(cruise_median) and cruise_median > 0):
        raise ValueError(f"cruise median tau must be a positive number, got {cruise_median!r}")
    tau = pd.to_numeric(table["tau"], errors="coerce")
    edge = table["edge"].map(_truthy) if "edge" in table else pd.Series(False, index=table.index)
    slow = (tau >= factor * cruise_median) | edge
    return pd.DataFrame({"cast_id": table["cast"].astype(str), "tau_s": tau.round(1), "slow_response": slow.astype(bool)})


def write_response(cruise_dir, flags: pd.DataFrame, source: str, overwrite: bool = False) -> Path:
    """Write the flags (cast_id, tau_s, slow_response). Refuses to replace an existing file unless overwrite;
    needs a source."""
    if not source or not str(source).strip():
        raise ValueError("a source is required: say where the taus and the rule come from")
    if flags.empty:
        raise ValueError("no casts to write")
    path = response_path(cruise_dir)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists:\n{path.read_text(encoding='utf-8-sig').strip()[:400]}\n"
                              "use --force to replace it")
    path.parent.mkdir(parents=True, exist_ok=True)
    out = flags[["cast_id", "tau_s", "slow_response"]].copy()
    out["source"] = str(source).strip()
    out.sort_values("cast_id").to_csv(path, index=False, encoding="utf-8")
    return path
