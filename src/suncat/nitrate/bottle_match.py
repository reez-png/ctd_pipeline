"""Pairing bottles with SUNA profiles, one rule for every tool that asks "which estimate matches the bottles?".

step02b_suna_response_correction and tools/diagnose_suna_lag.py (through estimate_response_time --bottles) both
compare the down-cast with lag-corrected profiles against the bottles. Until Nov 2026 each had its own bottle
list, down-cast rule and inclusion rule, so on P45_07 they compared 59 and 46 bottles and disagreed on R2 (0.768
against 0.820) for the same estimator. Here is the one rule:

  bottles   metadata/bottle_nitrate/<CRUISE>_bottle_nitrate.csv: rows with no flag and a finite nitrate and
            depth, after step03's replicate QC (a replicate group spanning more than REPLICATE_MAX_DIFF umol/L
            is dropped). Each remaining row is one bottle, as in step03's match-ups; a cast given as a number
            names the cast id with that number (suncat.casts.match_cast).
  down      step07's rule on the down-cast: 1 m depth cells centred on the bottle depth, within +-WIN_M; the
            mean of the frames in each cell, then the median of those cell means. Only down-cast frames with a
            value whose package was MOVING (|w| >= MOVING_W m/s), so a soak or a heave stop on the way down
            does not count.
  profiles  each fitted profile (N@0, N@tau, ...) read at the bottle depth by linear interpolation, only where
            the bottle lies within that cast's profile grid (no extrapolation).
  included  a bottle enters the comparison when EVERY required estimator has a value for it; each bottle that
            does not gets a reason (no frames on this cast, no profile for this cast, outside the profile's
            depth range, no moving down-cast frame within +-WIN_M), so a count can be traced.
  scores    sensor = offset + gain * bottle by OLS (gain with a t interval), R2, LOOCV RMSE in bottle units,
            and paired LOOCV squared errors against the reference (exact two-sided sign test).

Which frames are paired stays each tool's own choice and is printed beside the scores: step02b fits the L2
co-located product (SUNA frame times, step02's clock applied to the depths), the diagnosis fits step03's
frames (merged onto the CTD 1 s record, down-cast limited to scans Loop Edit kept).
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from suncat.casts import match_cast

WIN_M = 2.0                 # down-cast cells within +-WIN_M of the bottle (step07)
MOVING_W = 0.2              # m/s: slower frames are a soak or a stop, not the down-cast profile
REPLICATE_MAX_DIFF = 0.5    # umol/L: step03's replicate QC (Zheng)

REASONS = ("no frames on this cast", "no profile for this cast", "outside the profile's depth range",
           f"no moving down-cast frame within +-{WIN_M:g} m")


# ------------------------------------------------------------------ bottles
def load_bottles(path, replicate_max_diff: float = REPLICATE_MAX_DIFF) -> pd.DataFrame:
    """cast, depth_m, bottle from a bottle nitrate table (cast_id, depth_m, nitrate_uM; optional flag and
    replicate_id). Empty frame when the file is missing."""
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=["cast", "depth_m", "bottle"])
    b = pd.read_csv(path, encoding="utf-8-sig")
    b.columns = [c.strip() for c in b.columns]
    b["nitrate_uM"] = pd.to_numeric(b["nitrate_uM"], errors="coerce")
    b["depth_m"] = pd.to_numeric(b["depth_m"], errors="coerce")
    if "flag" in b.columns:
        b = b[b["flag"].fillna("").astype(str).str.strip().eq("")]
    b = b.dropna(subset=["cast_id", "nitrate_uM", "depth_m"])
    if "replicate_id" in b.columns and b["replicate_id"].notna().any():
        spread = b.groupby(["cast_id", "depth_m", "replicate_id"], dropna=False)["nitrate_uM"].transform(
            lambda v: v.max() - v.min() if v.notna().sum() >= 2 else 0.0)
        b = b[~(spread > replicate_max_diff)]                     # file order kept
    return pd.DataFrame({"cast": b["cast_id"].astype(str).str.strip(), "depth_m": b["depth_m"].astype(float),
                         "bottle": b["nitrate_uM"].astype(float)}).reset_index(drop=True)


# ------------------------------------------------------------------ the down-cast estimator
def moving(t_s, z, w_min: float = MOVING_W) -> np.ndarray:
    """|dz/dt| >= w_min, from a centred difference over the neighbouring frames (one-sided at the ends); t in
    seconds, sorted. False where the speed cannot be computed."""
    t, z = np.asarray(t_s, float), np.asarray(z, float)
    if len(t) < 2:
        return np.zeros(len(t), bool)
    i0 = np.clip(np.arange(len(t)) - 1, 0, len(t) - 1)
    i1 = np.clip(np.arange(len(t)) + 1, 0, len(t) - 1)
    dt = t[i1] - t[i0]
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.where(dt > 0, (z[i1] - z[i0]) / dt, np.nan)
    return np.isfinite(w) & (np.abs(w) >= w_min)


def down_value(z, y, use, z0: float, win: float = WIN_M) -> tuple[float, int]:
    """Median of 1 m cell means within +-win of z0, cells centred on z0, over the frames in `use` with a value.
    Returns (value, cells); (NaN, 0) when no frame is in range."""
    z, y, use = np.asarray(z, float), np.asarray(y, float), np.asarray(use, bool)
    ok = use & np.isfinite(z) & np.isfinite(y)
    c = np.floor(z[ok] - z0 + 0.5)
    m = np.abs(c) <= win
    if not m.any():
        return np.nan, 0
    means = pd.Series(y[ok][m]).groupby(c[m]).mean()
    return float(np.median(means)), int(means.size)


# ------------------------------------------------------------------ pairing
def pair(bottles: pd.DataFrame, frames: dict, profiles: dict, names, win: float = WIN_M) -> pd.DataFrame:
    """One row per bottle: cast, depth_m, bottle, down, one column per profile name, and `reason` (empty when
    every estimator has a value, else why the bottle is left out).
    frames:   {cast: dict(z, y, down)} with down already limited to the moving down-cast frames
    profiles: {cast: dict(grid, <name>: values on grid, ...)}; a name whose values are None is missing"""
    rows = []
    for _, b in bottles.iterrows():
        cast = match_cast(b["cast"], frames.keys()) or str(b["cast"])     # a cast number names its cast id
        zb = float(b["depth_m"])
        rec = dict(cast=cast, depth_m=zb, bottle=float(b["bottle"]), down=np.nan, **{n: np.nan for n in names},
                   reason="")
        F, P = frames.get(cast), profiles.get(cast)
        if F is None:
            rec["reason"] = REASONS[0]
            rows.append(rec)
            continue
        rec["down"] = down_value(F["z"], F["y"], F["down"], zb, win)[0]
        if P is None:
            rec["reason"] = REASONS[1]
        elif not (P["grid"][0] <= zb <= P["grid"][-1]):
            rec["reason"] = REASONS[2]
        else:
            for n in names:
                if P.get(n) is not None:
                    rec[n] = float(np.interp(zb, P["grid"], P[n]))
            if any(not np.isfinite(rec[n]) for n in names):
                rec["reason"] = REASONS[1]
        if not rec["reason"] and not np.isfinite(rec["down"]):
            rec["reason"] = f"no moving down-cast frame within +-{win:g} m"
        rows.append(rec)
    return pd.DataFrame(rows, columns=["cast", "depth_m", "bottle", "down", *names, "reason"])


# ------------------------------------------------------------------ scores
def _t975(df: int) -> float:
    table = [12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145,
             2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048,
             2.045, 2.042]
    return table[df - 1] if 1 <= df <= len(table) else 1.96


def ols(b, s) -> dict:
    """sensor = offset + gain * bottle: gain with its 95% t interval, offset, R2."""
    b, s = np.asarray(b, float), np.asarray(s, float)
    (g, o), cov = np.polyfit(b, s, 1, cov=True)
    se, tc = float(np.sqrt(cov[0, 0])), _t975(len(b) - 2)
    return dict(gain=float(g), lo=float(g - tc * se), hi=float(g + tc * se), offset=float(o),
                R2=float(1 - np.sum((s - o - g * b) ** 2) / np.sum((s - s.mean()) ** 2)))


def loo_errors(b, s) -> np.ndarray:
    """Leave-one-out prediction errors in bottle units: (s_i - offset_-i) / gain_-i - b_i."""
    b, s = np.asarray(b, float), np.asarray(s, float)
    e = np.full(len(b), np.nan)
    for i in range(len(b)):
        k = np.arange(len(b)) != i
        g, o = np.polyfit(b[k], s[k], 1)
        if g != 0:
            e[i] = (s[i] - o) / g - b[i]
    return e


def sign_test(d) -> float:
    """Exact two-sided sign test on paired differences (zeros dropped); NaN if none differ."""
    d = np.asarray(d, float)
    d = d[np.isfinite(d) & (d != 0)]
    n, k = len(d), int((d > 0).sum())
    if n == 0:
        return float("nan")
    return min(1.0, 2 * sum(math.comb(n, j) for j in range(min(k, n - k) + 1)) / 2 ** n)


def compare(table: pd.DataFrame, names, reference: str = "down", min_n: int = 5) -> dict:
    """Scores on the bottles every estimator has a value for (reason empty). Returns n, the counts left out by
    reason, and per estimator: gain, lo, hi, offset, R2, loocv; and for each non-reference estimator its paired
    LOOCV squared errors against the reference: median_d, better, n, p_sign."""
    est = [reference, *[n for n in names if n != reference]]
    both = table[table["reason"].eq("")].dropna(subset=est)
    out = dict(n=len(both), n_bottles=len(table), left_out=table.loc[table["reason"].ne(""), "reason"].value_counts().to_dict(),
               fits={}, both=both)
    if len(both) < max(min_n, 3):                                  # a gain interval needs three bottles
        return out
    bb = both["bottle"].to_numpy(float)
    loo = {}
    for e in est:
        s = both[e].to_numpy(float)
        f = ols(bb, s)
        loo[e] = loo_errors(bb, s)
        f["loocv"] = float(np.sqrt(np.nanmean(loo[e] ** 2)))
        out["fits"][e] = f
    for e in est[1:]:
        d = loo[e] ** 2 - loo[reference] ** 2
        out["fits"][e].update(median_d=float(np.nanmedian(d)), better=int((d < 0).sum()), n=int(len(d)),
                              p_sign=sign_test(d))
    return out


def report(cmp: dict, frames_label: str, unit: str = "umol/L", say=print) -> None:
    """Print the comparison: frames used, bottles kept and left out (by reason), the fits, the paired tests."""
    left = "; ".join(f"{v} {k}" for k, v in cmp["left_out"].items())
    say(f"  bottles: {cmp['n_bottles']} after flag and replicate QC; {cmp['n']} compared (every estimator has a value)"
        + (f"; left out: {left}" if left else ""))
    say(f"  frames : {frames_label}")
    if not cmp["fits"]:
        say("  too few bottles for a comparison")
        return
    names = list(cmp["fits"])
    say(f"  {'estimator':13s} {'gain (95% CI)':>27s} {'offset':>8s} {'R2':>6s}  LOOCV RMSE ({unit})")
    for e in names:
        f = cmp["fits"][e]
        say(f"  {e:13s} {f['gain']:7.4f} ({f['lo']:.4f} to {f['hi']:.4f}) {f['offset']:+8.3f} {f['R2']:6.3f}  {f['loocv']:.3f}")
    say(f"  paired LOOCV squared errors vs {names[0]} (negative = predicts the bottles better):")
    for e in names[1:]:
        f = cmp["fits"][e]
        say(f"    {e:13s} median {f['median_d']:+.3f}, better on {f['better']} of {f['n']}, sign test p = {f['p_sign']:.3f}")
