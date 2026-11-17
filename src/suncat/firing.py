"""Bottle depths from the firing record, not the sample sheet.

A sample sheet gives each Niskin a nominal depth (the target: 30 m, 15 m, ...). The bottle closed wherever the
rosette was when it fired, which on a cast fired on the fly can be metres away from the target: on P45_07 cast
02, Niskin 2 (nominal 30 m) closed at 27.5 m. The firing record is metadata/bottle_nitrate/<CRUISE>_bottle_ctd.csv
(build_bottle_nitrate_table.py, from the .btl files): cast_id, bottle_no, depth_m, ...

fired_depths() puts that depth on every sample row that names its cast and Niskin, keeps the sheet's value as
depth_nominal_m, and says per row which one it carries (depth_source "fired" or "nominal"). Rows it cannot match
keep the nominal depth, and the report counts them, so nothing changes silently.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from suncat.casts import match_cast


def firings_path(cruise_dir, cruise_id: str) -> Path:
    return Path(cruise_dir) / "metadata" / "bottle_nitrate" / f"{cruise_id}_bottle_ctd.csv"


def load_firings(cruise_dir, cruise_id: str) -> pd.DataFrame:
    """cast_id, bottle_no, depth_m from the firing record; empty when there is none."""
    f = firings_path(cruise_dir, cruise_id)
    if not f.exists():
        return pd.DataFrame(columns=["cast_id", "bottle_no", "depth_m"])
    b = pd.read_csv(f, encoding="utf-8-sig")
    b.columns = [c.strip().lower() for c in b.columns]
    out = pd.DataFrame({"cast_id": b["cast_id"].astype(str).str.strip(),
                        "bottle_no": pd.to_numeric(b["bottle_no"], errors="coerce"),
                        "depth_m": pd.to_numeric(b["depth_m"], errors="coerce")}).dropna()
    out["bottle_no"] = out["bottle_no"].astype(int)
    return out.drop_duplicates(["cast_id", "bottle_no"])


def fired_depths(samples: pd.DataFrame, firings: pd.DataFrame, cast_col: str = "cast_id",
                 niskin_col: str = "niskin", depth_col: str = "depth_m") -> tuple[pd.DataFrame, dict]:
    """samples with depth_col replaced by the fired depth where (cast, Niskin) is in the firing record.
    Adds depth_nominal_m and depth_source ("fired" or "nominal"). Returns (table, report): n, n_fired,
    n_nominal, and the nominal minus fired differences (median, min, max, in metres)."""
    s = samples.copy()
    # sheets often give whole-metre depths, read as integers: make the column float before fired depths go in
    # (pandas 3 refuses to put 27.48 into an int64 column)
    s[depth_col] = pd.to_numeric(s[depth_col], errors="coerce").astype(float)
    s["depth_nominal_m"] = s[depth_col].copy()
    s["depth_source"] = "nominal"
    rep = dict(n=len(s), n_fired=0, n_nominal=len(s), d_median=np.nan, d_min=np.nan, d_max=np.nan,
               has_record=not firings.empty, has_niskin=niskin_col in s.columns)
    if firings.empty or niskin_col not in s.columns:
        return s, rep
    casts = firings["cast_id"].unique()
    key = {(c, int(n)): z for c, n, z in zip(firings["cast_id"], firings["bottle_no"], firings["depth_m"])}
    nk = pd.to_numeric(s[niskin_col], errors="coerce")
    fired = []
    for c, n in zip(s[cast_col], nk):
        cid = match_cast(c, casts) if pd.notna(c) else None
        fired.append(key.get((cid, int(n))) if cid is not None and pd.notna(n) else None)
    fz = pd.Series([np.nan if v is None else float(v) for v in fired], index=s.index)
    hit = fz.notna()
    s.loc[hit, depth_col] = fz[hit]
    s.loc[hit, "depth_source"] = "fired"
    d = (s.loc[hit, "depth_nominal_m"] - fz[hit]).dropna()
    rep.update(n_fired=int(hit.sum()), n_nominal=int((~hit).sum()))
    if len(d):
        rep.update(d_median=float(d.median()), d_min=float(d.min()), d_max=float(d.max()))
    return s, rep


def describe(rep: dict, label: str = "samples") -> str:
    """One line for the run log."""
    if not rep["has_record"]:
        return (f"{label}: no firing record (<CRUISE>_bottle_ctd.csv): depths are the sample sheet's nominal depths "
                f"({rep['n']} rows); build it with build_bottle_nitrate_table.py")
    if not rep["has_niskin"]:
        return f"{label}: no Niskin column: depths are the sample sheet's nominal depths ({rep['n']} rows)"
    s = f"{label}: {rep['n_fired']} of {rep['n']} rows at their fired depth (firing record)"
    if rep["n_fired"]:
        s += (f"; nominal minus fired: median {rep['d_median']:+.1f} m, range {rep['d_min']:+.1f} to "
              f"{rep['d_max']:+.1f} m")
    if rep["n_nominal"]:
        s += f"; {rep['n_nominal']} left at the nominal depth (cast or Niskin not in the record)"
    return s
