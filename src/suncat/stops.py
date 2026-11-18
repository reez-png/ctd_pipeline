"""Where the rosette stopped: bottle depths from the CTD's own pressure record, for a cruise without a firing record.

On a cruise whose Niskins were fired without an SBE 55 auto-fire log (P45_06: no .afm, so no .btl and no
<CRUISE>_bottle_ctd.csv), the sample sheet's depth is the target, not where the bottle closed. The rosette still had
to stop for each bottle, and the stop is in the CTD record: a stretch of the up-cast where the pressure holds still.

find_stops()    stretches where the pressure, smoothed over smooth_s seconds (so ship heave neither breaks a stop up
                nor fakes one), changes by less than v_max dbar/s for at least min_s seconds, with the package in
                the water; the hold is that stretch widened by half the window at each end, its pressure the
                smoothed pressure over the stretch. Each stop is labelled down, bottom or up by where it falls
                against the deepest point, and an up-cast stop shallower than SURFACE_M as surface (the wait before
                recovery: on P45_06 it sat 0.6 to 2.1 m deep and took every 5 m sheet depth until Nov 2026). On synthetic casts (1 and 4 Hz; swell 0 to 1.5 dbar, 6 to 11 s; winch
                0.3 to 1 m/s; holds of 20 to 60 s) every hold is found, within 0.6 dbar. In a heavy short swell (1 to
                1.5 dbar at 6 s) on a 0.3 m/s up-cast the package itself stops with every swell, and the up-cast
                also shows brief false holds (about 22 s) every few metres: so match_depths takes the longest hold.
match_depths()  each sample depth to the longest up-cast or bottom stop within max(tol_m, tol_frac x depth) of it
                (the surface wait only for a sheet depth shallower than SURFACE_M);
                two sample depths that pick the same stop both keep their sheet depth, and say so.
cast_stops()    both, from one cast's .cnv (timeS, pressure, temperature, conductivity when present).

A stop shows the package was held there; it does not show a bottle was fired there. On P45_06 (Nov 2026) the
up-casts held only at the bottom and at the surface (and once at 12 m): no hold at any mid-depth sheet depth, as
expected of bottles fired while the rosette was moving, so build_bottle_table.py reports the stops there and keeps
the sheet depths (STOP_DEPTHS = "report"). The match is by depth, so a
sample whose sheet depth is far from every stop keeps its sheet depth, and a stop no sample matches is listed.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from suncat.calc.depth import pressure_to_depth
from suncat.io.cnv import read_cnv

SMOOTH_S = 16.0        # running mean over this many seconds before the velocity (two 8 s swell periods)
V_MAX = 0.20           # dbar/s: slower than this is held still (winch speeds are about 0.5 to 1 m/s)
MIN_S = 4.0            # the smoothed pressure holds still at least this long (a hold of about MIN_S + SMOOTH_S)
MERGE_S = 16.0         # still stretches this close in time and within MERGE_DBAR are one stop (heave splits holds)
MERGE_DBAR = 2.0
P_MIN = 0.3            # dbar: shallower than this is the surface, not a stop
TOL_M, TOL_FRAC = 5.0, 0.15
UPCAST_MIN_S = 30.0    # a record with less than this after its deepest point has no up-cast
SURFACE_M = 3.0        # an up-cast stop shallower than this is the wait at the surface before recovery

STOP_COLS = ["t0", "t1", "dur_s", "p_dbar", "p_sd", "depth_m", "temp_c", "leg"]


def find_stops(t, p, temp=None, wet=None, lat=5.0, smooth_s=SMOOTH_S, v_max=V_MAX, min_s=MIN_S, merge_s=MERGE_S,
               p_min=P_MIN) -> pd.DataFrame:
    """Stops in one cast. t (s), p (dbar), optional temp and wet (bool, e.g. conductivity above 1)."""
    t, p = np.asarray(t, float), np.asarray(p, float)
    ok = np.isfinite(t) & np.isfinite(p)
    if wet is not None:
        ok &= np.asarray(wet, bool)
    if ok.sum() < 10:
        return pd.DataFrame(columns=STOP_COLS)
    t, p = t[ok], p[ok]
    tc = None if temp is None else np.asarray(temp, float)[ok]
    o = np.argsort(t); t, p = t[o], p[o]
    if tc is not None:
        tc = tc[o]
    ps = pd.Series(p, index=pd.to_datetime(t, unit="s")).rolling(f"{smooth_s}s", center=True, min_periods=1).mean()
    ps = ps.to_numpy()
    v = np.gradient(ps, t)
    still = (np.abs(v) < v_max) & (p > p_min)
    edges = np.flatnonzero(np.diff(np.r_[0, still.astype(int), 0]))
    runs = [[s, e - 1] for s, e in zip(edges[::2], edges[1::2])]
    merged = []
    for s, e in runs:
        if merged and t[s] - t[merged[-1][1]] <= merge_s and abs(np.median(ps[s:e + 1]) - np.median(
                ps[merged[-1][0]:merged[-1][1] + 1])) < MERGE_DBAR:
            merged[-1][1] = e
        else:
            merged.append([s, e])
    imax = int(np.argmax(ps)); pmax = float(ps[imax])
    half = smooth_s / 2
    rows = []
    for s, e in merged:
        if t[e] - t[s] < min_s:
            continue
        # the smoothing shrinks a hold by about half a window at each end: the hold is the still core widened by
        # that much, its pressure the smoothed pressure over the core (heave averaged out)
        a, b = t[s] - half, t[e] + half
        held = (t >= a) & (t <= b)
        pm = float(np.median(ps[s:e + 1]))
        if s <= imax <= e or (abs(pm - pmax) < 2.0 and abs(t[s] - t[imax]) < 120):
            leg = "bottom"
        else:
            leg = "down" if e < imax else "up"
        if leg == "up" and pressure_to_depth(pm, lat) < SURFACE_M:
            leg = "surface"
        rows.append(dict(t0=a, t1=b, dur_s=b - a, p_dbar=pm, p_sd=float(np.std(p[held])),
                         depth_m=float(pressure_to_depth(pm, lat)),
                         temp_c=float(np.nanmedian(tc[held])) if tc is not None else np.nan, leg=leg))
    return pd.DataFrame(rows, columns=STOP_COLS)


def has_upcast(t, p) -> bool:
    t, p = np.asarray(t, float), np.asarray(p, float)
    ok = np.isfinite(t) & np.isfinite(p)
    if ok.sum() < 10:
        return False
    t, p = t[ok], p[ok]
    return t[-1] - t[int(np.argmax(p))] >= UPCAST_MIN_S


def match_depths(depths, stops: pd.DataFrame, tol_m=TOL_M, tol_frac=TOL_FRAC) -> pd.DataFrame:
    """One row per sample depth: nominal_m, stop (index into stops or -1), stop_depth_m, dz_m (stop minus sheet),
    dur_s, temp_c, note."""
    out = []
    for z in depths:
        r = dict(nominal_m=float(z), stop=-1, stop_depth_m=np.nan, dz_m=np.nan, dur_s=np.nan, temp_c=np.nan, note="")
        # the surface wait counts only for a sheet depth that is itself at the surface
        cand = stops[stops["leg"].isin(["up", "bottom"] + (["surface"] if z < SURFACE_M else []))]
        if cand.empty or not np.isfinite(z):
            r["note"] = "no up-cast stop in the record"
        else:
            d = (cand["depth_m"] - z).abs()
            near = d[d <= max(tol_m, tol_frac * z)]
            # the longest hold within tolerance: a bottle stop outlasts the brief holds heave can make on a slow up-cast
            i = cand.loc[near.index, "dur_s"].idxmax() if len(near) else d.idxmin()
            if len(near):
                r.update(stop=int(i), stop_depth_m=float(cand.at[i, "depth_m"]), dz_m=float(cand.at[i, "depth_m"] - z),
                         dur_s=float(cand.at[i, "dur_s"]), temp_c=float(cand.at[i, "temp_c"]))
            else:
                r["note"] = f"nearest up-cast stop {cand.at[i, 'depth_m']:.1f} m is {d[i]:.1f} m away"
        out.append(r)
    out = pd.DataFrame(out, columns=["nominal_m", "stop", "stop_depth_m", "dz_m", "dur_s", "temp_c", "note"])
    shared = out["stop"].ge(0) & out["stop"].duplicated(keep=False)
    if shared.any():
        for i in out.index[shared]:
            others = out.loc[shared & (out["stop"] == out.at[i, "stop"]) & (out.index != i), "nominal_m"].tolist()
            out.loc[i, "note"] = (f"shares the {out.at[i, 'stop_depth_m']:.1f} m stop with sheet depth(s) "
                                  f"{', '.join(f'{x:g}' for x in others)} m: sheet depth kept")
        out.loc[shared, ["stop", "stop_depth_m", "dz_m", "dur_s", "temp_c"]] = [-1, np.nan, np.nan, np.nan, np.nan]
    return out


def _col(d, *starts):
    return next((c for c in d.columns if c.lower().startswith(starts)), None)


def cast_stops(paths, lat=5.0, **kw) -> tuple[pd.DataFrame, str]:
    """Stops from the first .cnv in paths that has an up-cast. Returns (stops, the file used or why none)."""
    why = []
    for path in paths:
        path = Path(path)
        if not path.exists():
            continue
        d = read_cnv(path, mask_badflag=True)
        if d is None or d.empty:
            why.append(f"{path.name}: no data"); continue
        tcol, pcol = _col(d, "times"), _col(d, "prdm", "prm", "pr")
        if tcol is None or pcol is None:
            why.append(f"{path.name}: no time or pressure column"); continue
        if not has_upcast(d[tcol], d[pcol]):
            why.append(f"{path.name}: no up-cast"); continue
        tmp, cnd = _col(d, "t090", "tv290", "t068"), _col(d, "c0", "cond")
        wet = (d[cnd].to_numpy(float) > 1.0) if cnd else None
        return find_stops(d[tcol], d[pcol], d[tmp] if tmp else None, wet, lat, **kw), path.name
    return pd.DataFrame(columns=STOP_COLS), "; ".join(why) or "no .cnv found"
