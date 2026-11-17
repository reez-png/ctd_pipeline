"""Sensor response time (tau) from each cast's own profiles: a first-order lag fitted to the down-cast and the up-cast
together, without bottles.

    python tools/estimate_response_time.py P45_07                      SUNA nitrate (default)
    python tools/estimate_response_time.py P45_07 --var chl            chlorophyll fluorescence (4 Hz CTD scans)
    python tools/estimate_response_time.py P45_07 --col no3_stage1     another nitrate column
    python tools/estimate_response_time.py P45_07 --tau-max 120        a longer grid

Why not the offset sweep (compare_match_estimators.py --sweep)? A lag both delays and smooths the profile. A
depth shift can undo the delay but not the smoothing, so the sweep's tau saturates: on synthetic first-order lags
of 12 to 25 s it gave 6.7 to 9.5 s. This tool models the smoothing too.

Model. A first-order sensor obeys dy/dt = (x(t) - y(t)) / tau: it relaxes towards the water it is in, x(t), with
time constant tau (the time to cover 63% of a step). The water is a profile N(z), the same on both legs, so
x(t) = N(z(t)) along the package's path z(t), stops and heave included. For a given tau the reading y is LINEAR in
N: integrating the equation exactly between frames, with x varying linearly over each interval,
    y_i = a y_(i-1) + x_i [1 - r] + x_(i-1) [r - a],   a = exp(-dt/tau),   r = (tau/dt)(1 - a),
and x_i = N(z_i) interpolated on a depth grid. So for each tau, N is found by least squares (with a small penalty
on its curvature, which keeps cells with few frames from fitting noise), and the residual sum of squares RSS(tau)
says how well one profile seen through that lag explains both legs. The least-squares tau is the minimum of RSS
over a grid, refined by a parabola through the three lowest points. With tau = 0 the model is "both legs see the
same profile", so the drop from RSS(0) to RSS(tau) measures the hysteresis the lag explains. The frames run in
time order through the whole cast; only clean frames with a value on a leg (down-cast scans Loop Edit kept, up-cast)
enter the RSS, but every frame with a time and a depth carries the sensor's state forward (soak, stops).

Uncertainty. Frames are serially correlated, so RSS-based (likelihood) intervals for one cast would be too
narrow. Each cast is an independent estimate instead: the cruise value is the median over casts, with a bootstrap
over casts (90% interval) and the interquartile range. A pooled estimate (the minimum of RSS summed over casts)
is given too; it weights long casts more.

Tests (synthetic cruise, 1 Hz frames, heave on the down-cast, stops on the up-cast, lag integrated on a fine time
step): true tau 0, 2, 5, 16, 20, 30, 60 s gave 0.0, 2.0, 5.1, 16.1, 20.3, 30.1, 60.1 s (median over 4 casts), with
noise of 0.3 to 1.0 umol/L. Unlike the sweep, it does not saturate.

Assumptions, and what breaks them: the water column does not change between the down-cast and the up-cast. If the
up-cast sees the profile displaced by d metres (internal waves, drift), tau is biased by about d / (v_down + v_up)
(synthetic: +1 m gave +0.5 s at 1.8 m/s combined, +3 m gave +1.6 to +2.8 s); random displacements average out
over casts, a systematic one does not; the lag is first order (a sensor
with a dead time plus a lag, or a flow cell that mixes, shows a curved RSS minimum and a tau between the two); the
time stamps are right. A cast needs both legs (a single leg cannot separate a lag from the profile).

Output: a per-cast table on screen and in <cruise>/_audit/<cruise>__response_time_<var>.csv, with the RSS curve
in <cruise>/_audit/<cruise>__response_time_<var>_curve.csv. Reads only; changes no product.
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("compare_match_estimators", _HERE / "compare_match_estimators.py")
CME = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CME)
pp = CME.pp

DEFAULTS = {"nitrate": dict(tau_max=90.0, step=1.0, dz=1.0), "chl": dict(tau_max=12.0, step=0.25, dz=0.5)}
SMOOTH = 0.1          # curvature penalty weight, in (sensor units per m^2)^-2 relative to one frame's squared error
MIN_FRAMES = 30       # per leg, with a value


# ------------------------------------------------------------------ the model
def _interp(grid: np.ndarray, z: np.ndarray):
    k = np.clip(np.searchsorted(grid, z, side="right") - 1, 0, len(grid) - 2)
    f = np.clip((z - grid[k]) / (grid[k + 1] - grid[k]), 0.0, 1.0)
    return k, f


def lag_design(t: np.ndarray, z: np.ndarray, grid: np.ndarray, tau: float) -> np.ndarray:
    """Matrix A with y = A @ N: the first-order sensor's reading at each frame for the profile N on grid."""
    n, G = len(t), len(grid)
    k, f = _interp(grid, z)
    A = np.zeros((n, G))
    rows = np.arange(n)
    A[rows, k] += 1 - f
    A[rows, k + 1] += f
    if tau <= 0:
        return A
    M = A.copy()                                           # x_i = M[i] @ N
    dt = np.diff(t)
    a = np.exp(-dt / tau)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(dt > 0, tau / dt * (1 - a), 1.0)
    w1, w0 = 1 - r, r - a                                  # weights on x_i and x_(i-1)
    for i in range(1, n):
        A[i] = a[i - 1] * A[i - 1] + w1[i - 1] * M[i] + w0[i - 1] * M[i - 1]
    return A


def fit_profile(t, z, y, grid, tau, smooth=SMOOTH):
    """Least-squares profile N for one tau; returns (rss, N, residuals at the frames with a value)."""
    A = lag_design(t, z, grid, tau)
    ok = np.isfinite(y)
    G = len(grid)
    D = np.zeros((G - 2, G))
    i = np.arange(G - 2)
    D[i, i], D[i, i + 1], D[i, i + 2] = 1.0, -2.0, 1.0
    dz = float(np.median(np.diff(grid)))
    lhs = np.vstack([A[ok], np.sqrt(smooth) * D / dz ** 2])
    rhs = np.concatenate([y[ok], np.zeros(G - 2)])
    N, *_ = np.linalg.lstsq(lhs, rhs, rcond=None)
    res = y[ok] - A[ok] @ N
    return float(np.sum(res ** 2)), N, res


def refine_min(taus: np.ndarray, rss: np.ndarray):
    """Grid minimum refined by a parabola through it and its neighbours. edge=True when the minimum is the LAST
    grid value (no minimum within the grid). A minimum at tau = 0 is an estimate (no lag resolved), not an edge."""
    j = int(np.nanargmin(rss))
    if j == len(taus) - 1:
        return float(taus[j]), True
    if j == 0:
        return float(taus[0]), False
    x, yv = taus[j - 1:j + 2], rss[j - 1:j + 2]
    den = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
    A_ = (x[2] * (yv[1] - yv[0]) + x[1] * (yv[0] - yv[2]) + x[0] * (yv[2] - yv[1])) / den
    B_ = (x[2] ** 2 * (yv[0] - yv[1]) + x[1] ** 2 * (yv[2] - yv[0]) + x[0] ** 2 * (yv[1] - yv[2])) / den
    if A_ <= 0:
        return float(taus[j]), False
    return float(np.clip(-B_ / (2 * A_), x[0], x[2])), False


def cast_tau(s: pd.DataFrame, taus: np.ndarray, dz: float, smooth: float = SMOOTH) -> dict:
    """tau for one cast's samples (t, z, v, down, up as compare_match_estimators reads them)."""
    s = s[s["t"].notna() & s["z"].notna()].sort_values("t", kind="stable").copy()
    # residuals only on the legs (down: Loop Edit kept; up): surface soak and flagged scans carry the sensor's
    # state forward in time but are not fitted
    s["v"] = s["v"].where(s["down"].astype(bool) | s["up"].astype(bool))
    out = dict(n=int(s["v"].notna().sum()), n_down=int((s["v"].notna() & s["down"]).sum()),
               n_up=int((s["v"].notna() & s["up"]).sum()), tau=np.nan, edge=False, rss0=np.nan, rss_min=np.nan,
               rms0=np.nan, rms=np.nan, hyst0=np.nan, hyst=np.nan, curve=None)
    if min(out["n_down"], out["n_up"]) < MIN_FRAMES:
        out["skip"] = f"too few frames with a value on a leg (down {out['n_down']}, up {out['n_up']})"
        return out
    t = (s["t"] - s["t"].iloc[0]).dt.total_seconds().to_numpy(float)
    z = s["z"].to_numpy(float)
    y = s["v"].to_numpy(float)
    grid = np.arange(np.floor(z.min() / dz) * dz, np.ceil(z.max() / dz) * dz + dz / 2, dz)
    if len(grid) < 4:
        out["skip"] = "profile shorter than four grid cells"
        return out
    rss = np.array([fit_profile(t, z, y, grid, tau, smooth)[0] for tau in taus])
    tau_hat, edge = refine_min(taus, rss)
    ok = np.isfinite(y)
    is_up, is_dn = s["up"].to_numpy(bool)[ok], s["down"].to_numpy(bool)[ok]

    def hyst(res):                                # mean misfit up minus down: the hysteresis left unexplained
        return float(np.mean(res[is_up]) - np.mean(res[is_dn])) if is_up.any() and is_dn.any() else np.nan

    r0 = fit_profile(t, z, y, grid, 0.0, smooth)
    r1 = fit_profile(t, z, y, grid, tau_hat, smooth)
    n = int(ok.sum())
    out.update(tau=np.nan if edge else tau_hat, tau_grid=tau_hat, edge=edge, rss0=r0[0], rss_min=r1[0],
               rms0=float(np.sqrt(r0[0] / n)), rms=float(np.sqrt(r1[0] / n)), hyst0=hyst(r0[2]), hyst=hyst(r1[2]),
               curve=rss)
    return out


# ------------------------------------------------------------------ a cruise
def read_cast(cruise: str, var: str, cast_file: Path, col: str | None):
    if var == "chl":
        return CME.read_ctd_scans(cast_file)
    return CME.read_suna_frames_df(pd.read_csv(cast_file), vcol=col)


def estimate(cruise: str, var: str = "nitrate", col: str | None = None, tau_max: float | None = None,
             step: float | None = None, dz: float | None = None, smooth: float = SMOOTH, n_boot: int = 2000,
             verbose: bool = True, write: bool = True) -> dict:
    say = print if verbose else (lambda *a, **k: None)
    d = DEFAULTS[var]
    tau_max, step, dz = tau_max or d["tau_max"], step or d["step"], dz or d["dz"]
    taus = np.arange(0.0, tau_max + step / 2, step)
    folder, pattern = (pp.product_dir(cruise, "ctd_derived"), "*_der.cnv") if var == "chl" else \
                      (pp.product_dir(cruise, "nitrate_qc"), "*_nitrate_QC.csv")
    files = sorted(folder.glob(pattern))
    if not files:
        raise FileNotFoundError(f"no {pattern} in {folder}")
    unit = "mg/m3" if var == "chl" else "umol/L"
    say(f"{cruise} {var}: first-order lag fitted to each cast's down and up profiles together (no bottles)")
    say(f"  tau grid 0 to {tau_max:g} s, step {step:g} s; profile grid {dz:g} m; curvature penalty {smooth:g}"
        + (f"; value column {col}" if col else ""))
    rows, curves = [], []
    for f in files:
        cast = f.name.replace("_der.cnv", "").replace("_nitrate_QC.csv", "")
        r = cast_tau(read_cast(cruise, var, f, col), taus, dz, smooth)
        rows.append({k: v for k, v in dict(cast=cast, **r).items() if k != "curve"})
        if r["curve"] is not None:
            curves.append(pd.DataFrame({"cast": cast, "tau_s": taus, "rss": r["curve"]}))
    t = pd.DataFrame(rows)
    good = t[t["tau"].notna()]
    say(f"\n  {'cast':16s} {'frames':>6s} {'down/up':>9s} {'tau (s)':>8s}   RMS misfit {unit}: tau=0 -> tau   "
        f"up-minus-down misfit: tau=0 -> tau")
    for _, r in t.iterrows():
        if "skip" in r and isinstance(r.get("skip"), str):
            say(f"  {r['cast']:16s} skipped: {r['skip']}")
            continue
        tau_s = "  EDGE" if r["edge"] else f"{r['tau']:8.1f}"
        say(f"  {r['cast']:16s} {r['n']:6d} {r['n_down']:4d}/{r['n_up']:<4d} {tau_s:>8s}   "
            f"{r['rms0']:8.3f} -> {r['rms']:6.3f}            {r['hyst0']:+8.3f} -> {r['hyst']:+7.3f}")
    res = dict(table=t, taus=taus)
    n_edge = int(t["edge"].sum())
    if n_edge:
        say(f"\n  {n_edge} cast(s) at the end of the grid (EDGE): no minimum below {tau_max:g} s; left out. "
            "If many, widen --tau-max.")
    if len(good) < 3:
        say("  fewer than 3 casts with an estimate: no cruise value")
        return res
    v = good["tau"].to_numpy(float)
    rng = np.random.default_rng(0)
    boot = np.median(rng.choice(v, (n_boot, len(v)), replace=True), axis=1)
    res.update(tau_median=float(np.median(v)), tau_iqr=(float(np.percentile(v, 25)), float(np.percentile(v, 75))),
               tau_ci90=(float(np.percentile(boot, 5)), float(np.percentile(boot, 95))), n_casts=len(v))
    if curves:
        cv = pd.concat(curves)
        cv = cv[cv["cast"].isin(good["cast"])]
        tot = cv.groupby("tau_s")["rss"].sum()
        res["tau_pooled"], _ = refine_min(tot.index.to_numpy(float), tot.to_numpy(float))
    share = 1 - good["rss_min"].sum() / good["rss0"].sum()
    res["rss_share_explained"] = float(share)
    say(f"\n  tau over {len(v)} casts: median {res['tau_median']:.1f} s (IQR {res['tau_iqr'][0]:.1f} to "
        f"{res['tau_iqr'][1]:.1f} s); bootstrap over casts, 90% interval {res['tau_ci90'][0]:.1f} to "
        f"{res['tau_ci90'][1]:.1f} s")
    if "tau_pooled" in res:
        say(f"  pooled (RSS summed over casts): {res['tau_pooled']:.1f} s")
    say(f"  the lag explains {share:.0%} of the misfit of a lag-free model (RSS, summed over these casts); "
        f"median up-minus-down misfit {good['hyst0'].median():+.3f} -> {good['hyst'].median():+.3f} {unit}")
    say("  tau is the e-folding time: the sensor covers 63% of a step change in tau seconds, 95% in 3 tau.")
    if write:
        out = pp.cruise_dir(cruise) / "_audit" / f"{cruise}__response_time_{var}.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        t.to_csv(out, index=False, encoding="utf-8-sig")
        if curves:
            pd.concat(curves).to_csv(out.with_name(out.stem + "_curve.csv"), index=False, encoding="utf-8-sig")
        say(f"\n  per-cast table -> {out}")
        res["out"] = out
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cruise")
    ap.add_argument("--var", choices=["nitrate", "chl"], default="nitrate")
    ap.add_argument("--col", default=None, help="nitrate column (default the latest stage: no3_stage2, no3_stage1, no3_qc)")
    ap.add_argument("--tau-max", type=float, default=None, help="longest tau tried, s (default 90 nitrate, 12 chl)")
    ap.add_argument("--step", type=float, default=None, help="tau grid step, s (default 1 nitrate, 0.25 chl)")
    ap.add_argument("--dz", type=float, default=None, help="profile grid, m (default 1 nitrate, 0.5 chl)")
    ap.add_argument("--smooth", type=float, default=SMOOTH, help=f"curvature penalty (default {SMOOTH:g})")
    a = ap.parse_args()
    estimate(a.cruise, a.var, a.col, a.tau_max, a.step, a.dz, a.smooth)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
