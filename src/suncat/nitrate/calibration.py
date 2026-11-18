"""Homogeneous-water-only bottle calibration of the SUNA (step03_suna_bottle_calibration).

A SUNA with a slow response reads the water it came from: on the up-cast through a nitracline it reads high, on the
down-cast low. Where the water is vertically homogeneous that does not matter, so the calibration uses only the
bottles whose depth shows the SUNA's down-cast and up-cast agreeing (|down - up| <= hyst_tol, in the SUNA's own
units), and fits SUNA = offset + gain x bottle through them.

The agreement test needs both legs. Until Nov 2026 it needed min_in_win clean frames in the window in total, so a
bottle could pass on one frame of one leg against several of the other: one frame cannot show the water is
homogeneous, and the window's pooled median is then the other leg's, lag and all. Now each leg needs min_per_leg
frames. On P45_06 that removed 5 of 18 bottles (all at 0 uM, with 1 or 2 down-cast or 1 up-cast frame): gain 0.877
to 0.862, offset -0.75 to -0.52. On P45_07 it changed nothing (6 bottles either way).

A fit through few bottles can rest on one or two of them: leave_one_out() gives the gain and offset without each
bottle. On P45_07 the 6 homogeneous bottles hold only 2 above 2 uM (casts 08 and 15, both near 16 uM); alone they
imply gains of 0.89 and 1.03, so the fit is 0.957 with cast 15 and 0.891 (offset -0.40) without it.

leg_medians()       down, up (median, count) and the pooled median of the clean frames within +-win dbar of a depth
homogeneous_pairs() one row per bottle: the leg medians and counts, |down - up|, the SUNA value, used, and why not
fit_line()          SUNA = offset + gain x bottle (or offset only), with R2, the residual SD and the covariance
leave_one_out()     the fit without each bottle in turn

The SUNA value at a bottle is the pooled median of the window's frames (value="pooled", what P45_07's adopted
calibration used) or the mean of the two leg medians (value="legs"), which cancels a symmetric lag to first order.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

DEPTH_WIN = 2.5        # dbar window around the bottle depth (bottle depth in m taken as dbar: within 0.6 % on the shelf)
HYST_TOL = 2.0         # uM: |down - up| at most this for the water to count as homogeneous
MIN_PER_LEG = 3        # clean frames needed on EACH leg inside the window

COLS = ["cast", "station", "depth_m", "bottle", "suna_down", "n_down", "suna_up", "n_up", "hyst", "suna", "n",
        "used", "note"]


def _clean(frames: pd.DataFrame) -> pd.DataFrame:
    ok = frames["L1_clean"].astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])
    return frames[ok & frames["pressure_db"].notna() & frames["nitrate_uM"].notna()]


def leg_medians(clean: pd.DataFrame, cast: str, depth: float, win: float = DEPTH_WIN) -> dict:
    d = clean[(clean["cast"] == cast) & (clean["pressure_db"] >= depth - win) & (clean["pressure_db"] <= depth + win)]
    dn = d.loc[d["leg"] == "down", "nitrate_uM"]
    up = d.loc[d["leg"] == "up", "nitrate_uM"]
    return dict(down=float(dn.median()) if len(dn) else np.nan, n_down=len(dn),
                up=float(up.median()) if len(up) else np.nan, n_up=len(up),
                pooled=float(d["nitrate_uM"].median()) if len(d) else np.nan, n=len(d))


def homogeneous_pairs(frames: pd.DataFrame, bottles: pd.DataFrame, win: float = DEPTH_WIN, hyst_tol: float = HYST_TOL,
                      min_per_leg: int = MIN_PER_LEG, include_flagged: bool = False, value: str = "pooled",
                      min_total: int | None = None) -> pd.DataFrame:
    """frames: the L2 co-located product (cast, leg, pressure_db, nitrate_uM, L1_clean). bottles: cast_id, depth_m,
    nitrate_uM, optional station and flag. min_total (the rule before Nov 2026: frames in the window in total,
    either leg) replaces the per-leg rule when given; kept to compare the two."""
    clean = _clean(frames)
    flag = bottles["flag"].fillna("").astype(str).str.strip() if "flag" in bottles else pd.Series("", bottles.index)
    rows = []
    for i, r in bottles.iterrows():
        row = dict(cast=r["cast_id"], station=r.get("station", ""), depth_m=r["depth_m"], bottle=r["nitrate_uM"],
                   suna_down=np.nan, n_down=0, suna_up=np.nan, n_up=0, hyst=np.nan, suna=np.nan, n=0, used=False,
                   note="")
        if flag[i] and not include_flagged:
            row["note"] = "flagged bottle excluded"
            rows.append(row); continue
        m = leg_medians(clean, r["cast_id"], float(r["depth_m"]), win)
        hyst = abs(m["down"] - m["up"]) if np.isfinite(m["down"]) and np.isfinite(m["up"]) else np.nan
        val = m["pooled"] if value == "pooled" else (np.nanmean([m["down"], m["up"]]) if m["n"] else np.nan)
        row.update(suna_down=m["down"], n_down=m["n_down"], suna_up=m["up"], n_up=m["n_up"], hyst=hyst, suna=val,
                   n=m["n"])
        if min_total is not None:
            enough, short = m["n"] >= min_total, f"<{min_total} frames"
        else:
            enough = m["n_down"] >= min_per_leg and m["n_up"] >= min_per_leg
            short = f"down {m['n_down']} / up {m['n_up']} frames (<{min_per_leg} on a leg)"
        if m["n"] == 0:
            row["note"] = "no SUNA match"
        elif not np.isfinite(hyst):
            row["note"] = "one leg only"
        elif not enough:
            row["note"] = short
        elif hyst > hyst_tol:
            row["note"] = f"stratified (down/up {hyst:.1f} uM)"
        else:
            row["used"] = True
        rows.append(row)
    return pd.DataFrame(rows, columns=COLS)


def fit_line(bottle, suna, fit_slope: bool = True) -> dict:
    b, s = np.asarray(bottle, float), np.asarray(suna, float)
    out = dict(n=len(b), gain=np.nan, offset=np.nan, R2=np.nan, resid_sd=np.nan, cov=None)
    if len(b) < 2:
        return out
    if fit_slope and np.ptp(b) > 0 and len(b) >= 3:
        (g, o), cov = np.polyfit(b, s, 1, cov=True)
    elif fit_slope and np.ptp(b) > 0:
        (g, o), cov = np.polyfit(b, s, 1), None
    else:
        g, o, cov = 1.0, float(np.mean(s - b)), None
    pred = o + g * b
    sst = float(np.sum((s - s.mean()) ** 2))
    out.update(gain=float(g), offset=float(o), cov=cov, resid_sd=float(np.std(s - pred)),
               R2=1 - float(np.sum((s - pred) ** 2)) / sst if sst > 0 else np.nan)
    return out


def leave_one_out(fit: pd.DataFrame, fit_slope: bool = True) -> pd.DataFrame:
    """fit: the used rows of homogeneous_pairs(). One row per bottle: cast, depth_m, bottle, and the gain and offset
    of the fit without it."""
    rows = []
    for i in fit.index:
        rest = fit.drop(index=i)
        r = fit_line(rest["bottle"], rest["suna"], fit_slope)
        rows.append(dict(cast=fit.at[i, "cast"], depth_m=fit.at[i, "depth_m"], bottle=fit.at[i, "bottle"],
                         gain=r["gain"], offset=r["offset"]))
    return pd.DataFrame(rows, columns=["cast", "depth_m", "bottle", "gain", "offset"])
