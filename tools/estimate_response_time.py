"""Sensor response time (tau) from each cast's own profiles: a first-order lag fitted to the down-cast and the up-cast
together, without bottles.

    python tools/estimate_response_time.py P45_07                      SUNA nitrate (default)
    python tools/estimate_response_time.py P45_07 --var chl            chlorophyll fluorescence (4 Hz CTD scans)
    python tools/estimate_response_time.py P45_07 --col no3_stage1     another nitrate column
    python tools/estimate_response_time.py P45_07 --tau-max 120        a longer grid
    python tools/estimate_response_time.py P45_07 --bottles            also: do lag-corrected profiles match the
                                                                       bottles better than the down-cast?

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
import os
import sys
from pathlib import Path

# these matrices are small (hundreds of columns): a multi-threaded BLAS spends more time coordinating threads than
# computing, badly so on many-core Windows machines. One thread, unless the user has set otherwise. (Takes effect
# only when this tool is the program; under pytest numpy is already loaded.)
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("compare_match_estimators", _HERE / "compare_match_estimators.py")
CME = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CME)
pp = CME.pp
sys.path.insert(0, str(CME.REPO / "src"))                 # suncat, also when it is not installed
import suncat.nitrate.hysteresis as H  # noqa: E402
from suncat.nitrate.hysteresis import (SMOOTH, MIN_FRAMES, lag_design, fit_profile, refine_min,  # noqa: E402,F401
                                       search_tau, summarise)
from suncat.timeutil import elapsed_seconds  # noqa: E402

DEFAULTS = {"nitrate": dict(tau_max=90.0, step=1.0, dz=1.0), "chl": dict(tau_max=12.0, step=0.25, dz=0.5)}


# ------------------------------------------------------------------ the model (suncat.nitrate.hysteresis)
def cast_tau(s: pd.DataFrame, taus: np.ndarray, dz: float, smooth: float = SMOOTH) -> dict:
    """tau for one cast's samples (t, z, v, down, up as compare_match_estimators reads them). Residuals only on
    the legs (down: Loop Edit kept; up): surface soak and flagged scans carry the sensor's state forward in time
    but are not fitted."""
    s = s[s["t"].notna() & s["z"].notna()]
    return H.cast_tau(elapsed_seconds(s["t"]), s["z"].to_numpy(float), s["v"].to_numpy(float),
                      s["down"].to_numpy(bool), s["up"].to_numpy(bool), taus, dz, smooth)


# ------------------------------------------------------------------ against the bottles
def _adopted(cruise):
    """The cruise's stated SUNA calibration (suna_calibration.json), or None."""
    try:
        from suncat.suna_cal import load_suna_calibration
        return load_suna_calibration(pp.cruise_dir(cruise))
    except Exception:
        return None


def bottle_check(cruise, var, res, samples, data, smooth=SMOOTH, say=print):
    """Does correcting the lag make the sensor match the bottles better? Per bottle, on the same bottles:
      down          the production estimator (step07's down-cast cell median, compare_match_estimators 'down')
      N@0           the profile fitted with no lag (both legs averaged, in effect)
      N@tau_cruise  the lag-corrected profile with the cruise's median tau (one constant for the sensor)
      N@tau_cast    the lag-corrected profile with each cast's own tau
    each read at the bottle depth, then sensor = offset + gain * bottle by OLS, LOOCV RMSE in bottle units, and
    paired LOOCV squared errors against 'down' (exact sign test; Wilcoxon when scipy is present).
    Then the depth offset (offset_sweep): the lag-corrected profile carries no lag, so reading it at the bottle depth
    + D and scoring each D measures any offset between where the bottle's water came from and the depth recorded
    for it (rosette geometry, firing timing, bottle closure), free of the lag that confounds the same sweep on raw
    legs."""
    unit = "mg/m3" if var == "chl" else "umol/L"
    bottles = CME.load_bottles(cruise, var)
    tau_c = res["tau_median"]
    per_cast = res["table"].set_index("cast")["tau"]
    prof = {}
    for cast, (t, z, y, grid) in data.items():
        tc = per_cast.get(cast, np.nan)
        prof[cast] = dict(grid=grid, n0=fit_profile(t, z, y, grid, 0.0, smooth)[1],
                          nc=fit_profile(t, z, y, grid, tau_c, smooth)[1],
                          nk=fit_profile(t, z, y, grid, tc, smooth)[1] if np.isfinite(tc) else None)
    rows = []
    for _, b in bottles.iterrows():
        cast, zb = b["cast"], float(b["depth_m"])
        rec = dict(cast=cast, depth_m=zb, bottle=float(b["bottle"]))
        if cast in prof:
            P = prof[cast]
            inside = P["grid"][0] <= zb <= P["grid"][-1]
            s = samples[cast]
            rec["down"] = CME.cell_median(s[s["down"].astype(bool) & s["moving"].astype(bool)], zb)[0]
            rec["N@0"] = float(np.interp(zb, P["grid"], P["n0"])) if inside else np.nan
            rec["N@tau_cruise"] = float(np.interp(zb, P["grid"], P["nc"])) if inside else np.nan
            rec["N@tau_cast"] = float(np.interp(zb, P["grid"], P["nk"])) if inside and P["nk"] is not None else np.nan
        rows.append(rec)
    T = pd.DataFrame(rows)
    est = [e for e in ("down", "N@0", "N@tau_cruise", "N@tau_cast") if e in T and T[e].notna().sum() >= 5]
    both = T.dropna(subset=est)
    out = dict(table=T, fits={}, loocv_rmse={}, paired={})
    say(f"\n  against the bottles: {len(T)} bottles, {len(both)} with every estimator; tau_cruise = {tau_c:.1f} s")
    if len(both) < 5 or "down" not in est:
        say("  too few bottles for a comparison")
        return out
    bb = both["bottle"].to_numpy(float)
    say(f"  {'estimator':13s} {'gain (95% CI)':>27s} {'offset':>8s} {'R2':>6s}  LOOCV RMSE ({unit})")
    loo = {}
    for e in est:
        f = CME.ols(bb, both[e].to_numpy(float))
        loo[e] = CME.loo_errors(bb, both[e].to_numpy(float))
        out["fits"][e], out["loocv_rmse"][e] = f, float(np.sqrt(np.nanmean(loo[e] ** 2)))
        say(f"  {e:13s} {f['gain']:7.4f} ({f['lo']:.4f} to {f['hi']:.4f}) {f['offset']:+8.3f} {f['R2']:6.3f}  "
            f"{out['loocv_rmse'][e]:.3f}")
    try:
        from scipy.stats import wilcoxon
    except ImportError:
        wilcoxon = None
    say("  paired LOOCV squared errors vs down (negative = predicts the bottles better):")
    for e in est:
        if e == "down":
            continue
        d = loo[e] ** 2 - loo["down"] ** 2
        p_s = CME.sign_test(d)
        p_w = float(wilcoxon(d).pvalue) if (wilcoxon is not None and np.any(d != 0)) else float("nan")
        out["paired"][e] = dict(median=float(np.nanmedian(d)), better=int((d < 0).sum()), n=len(d), p_sign=p_s,
                                p_wilcoxon=p_w)
        say(f"    {e:13s} median {np.nanmedian(d):+.3f}, better on {int((d < 0).sum())} of {len(d)}, sign test "
            f"p = {p_s:.3f}" + (f", Wilcoxon p = {p_w:.3f}" if np.isfinite(p_w) else ""))
    if "N@tau_cruise" in est:
        out["offset"] = offset_sweep(both, prof, unit, say=say)
    cal = _adopted(cruise) if var == "nitrate" else None
    if cal and cal.get("mode") == "adopted":
        f = out["fits"].get("N@tau_cruise")
        say(f"  the cruise's adopted calibration: gain {cal['a1']:g}, offset {cal['a0']:+g} ({cal['source']})")
        if f:
            inside = f["lo"] <= cal["a1"] <= f["hi"]
            say(f"    the lag-corrected fit's gain {f['gain']:.3f} (95% CI {f['lo']:.3f} to {f['hi']:.3f}) "
                + ("includes" if inside else "excludes") + " it")
    return out


OFFSETS = np.arange(-6.0, 6.001, 0.5)


def offset_sweep(both: pd.DataFrame, prof: dict, unit: str, deltas=OFFSETS, n_boot: int = 500,
                 say=print) -> dict:
    """Read each cast's lag-corrected profile (cruise tau) at bottle depth + D (D > 0 = deeper; beyond the profile's
    ends the end value is held, the mixed layer or the bottom water) and fit sensor = offset + gain * bottle at each
    D. Score: cross-validated R2, Q2 = 1 - PRESS / SS (scale-free, as in compare_match_estimators --sweep), with the
    LOOCV RMSE in bottle units alongside. The best D comes with the D values whose 1 - Q2 is within 10% of its own and a
    bootstrap over bottles (best D by in-sample R2 in each resample; 90% interval). The best D is chosen on the
    same bottles it is scored on, so its LOOCV error is optimistic; the interval is the honest part."""
    bb = both["bottle"].to_numpy(float)
    V = np.full((len(both), len(deltas)), np.nan)
    for i, (cast, zb) in enumerate(zip(both["cast"], both["depth_m"])):
        P = prof[cast]
        V[i] = np.interp(float(zb) + deltas, P["grid"], P["nc"])
    rows = []
    for j, D in enumerate(deltas):
        v = V[:, j]
        e = CME.loo_forward_errors(bb, v)
        q2 = 1.0 - float(np.sum(e ** 2)) / float(np.sum((v - v.mean()) ** 2))
        rows.append(dict(D=float(D), q2=q2, R2=CME.ols(bb, v)["R2"],
                         loocv=float(np.sqrt(np.nanmean(CME.loo_errors(bb, v) ** 2)))))
    T = pd.DataFrame(rows)
    i = int(T["q2"].idxmax())
    near = T.loc[(1 - T["q2"]) <= 1.1 * (1 - T.at[i, "q2"]), "D"]     # unexplained share within 10% of the best
    rng = np.random.default_rng(0)
    best = []
    for _ in range(n_boot):
        k = rng.integers(0, len(bb), len(bb))
        b = bb[k]
        if np.ptp(b) == 0:
            continue
        r2 = []
        for j in range(len(deltas)):
            v = V[k, j]
            g, o = np.polyfit(b, v, 1)
            r2.append(1 - np.sum((v - (o + g * b)) ** 2) / max(np.sum((v - v.mean()) ** 2), 1e-12))
        best.append(deltas[int(np.argmax(r2))])
    ci = (float(np.percentile(best, 5)), float(np.percentile(best, 95))) if best else (np.nan, np.nan)
    j0 = int(np.argmin(np.abs(deltas)))
    say(f"\n  depth offset, lag-corrected profile (cruise tau) read at bottle depth + D (D > 0 = deeper), "
        f"{len(bb)} bottles:")
    say(f"     D      Q2      R2   LOOCV RMSE ({unit})")
    for _, r in T.iterrows():
        say(f"  {r['D']:+5.1f}  {r['q2']:+.3f}  {r['R2']:.3f}  {r['loocv']:.3f}" + ("   <- best" if r.name == i else ""))
    say(f"  best D {T.at[i, 'D']:+.1f} m (Q2 {T.at[i, 'q2']:+.3f}; at D = 0: {T.at[j0, 'q2']:+.3f}); unexplained share "
        f"1 - Q2 within 10% of the best: {near.min():+.1f} to {near.max():+.1f} m; bootstrap over bottles, 90% interval {ci[0]:+.1f} to "
        f"{ci[1]:+.1f} m")
    say("  (D < 0: the bottles match water shallower than their recorded depth. With no lag left in the profile this")
    say("   is the bottle-depth offset itself: rosette geometry, firing timing or bottle closure, not the sensor.)")
    if T.at[i, "D"] in (deltas.min(), deltas.max()):
        say("  NOTE: the best D sits at the end of the range: a limit, not an estimate")
    return dict(table=T, best_D=float(T.at[i, "D"]), q2_best=float(T.at[i, "q2"]), q2_at_0=float(T.at[j0, "q2"]),
                near=(float(near.min()), float(near.max())), ci90=ci)


# ------------------------------------------------------------------ a cruise
def read_cast(cruise: str, var: str, cast_file: Path, col: str | None):
    if var == "chl":
        return CME.read_ctd_scans(cast_file)
    return CME.read_suna_frames_df(pd.read_csv(cast_file), vcol=col)


def estimate(cruise: str, var: str = "nitrate", col: str | None = None, tau_max: float | None = None,
             step: float | None = None, dz: float | None = None, smooth: float = SMOOTH, n_boot: int = 2000,
             verbose: bool = True, write: bool = True, bottles: bool = False) -> dict:
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
    rows, curves, data, samples, full = [], [], {}, {}, {}
    say(f"\n  {'cast':16s} {'frames':>6s} {'down/up':>9s} {'tau (s)':>8s}   RMS misfit {unit}: tau=0 -> tau   "
        f"up-minus-down misfit: tau=0 -> tau")
    for f in files:                                       # each row printed as it is done
        cast = f.name.replace("_der.cnv", "").replace("_nitrate_QC.csv", "")
        samples[cast] = read_cast(cruise, var, f, col)
        r = cast_tau(samples[cast], taus, dz, smooth)
        rows.append({k: v for k, v in dict(cast=cast, **r).items() if k not in ("curve", "data")})
        if r["curve"] is not None:
            curves.append(pd.DataFrame({"cast": cast, "tau_s": taus, "rss": r["curve"]}))
            data[cast] = r["data"]
            full[cast] = r
        if isinstance(r.get("skip"), str):
            say(f"  {cast:16s} skipped: {r['skip']}", flush=True)
            continue
        tau_s = "  EDGE" if r["edge"] else f"{r['tau']:8.1f}"
        say(f"  {cast:16s} {r['n']:6d} {r['n_down']:4d}/{r['n_up']:<4d} {tau_s:>8s}   "
            f"{r['rms0']:8.3f} -> {r['rms']:6.3f}            {r['hyst0']:+8.3f} -> {r['hyst']:+7.3f}", flush=True)
    t = pd.DataFrame(rows)
    good = t[t["tau"].notna()]
    res = dict(table=t, taus=taus)
    n_edge = int(t["edge"].sum())
    if n_edge:
        say(f"\n  {n_edge} cast(s) at the end of the grid (EDGE): no minimum below {tau_max:g} s; left out. "
            "If many, widen --tau-max.")
    if len(good) < 3:
        say("  fewer than 3 casts with an estimate: no cruise value")
        return res
    res.update(summarise(full, taus, smooth, n_boot))
    share = res["rss_share_explained"]
    say(f"\n  tau over {res['n_casts']} casts: median {res['tau_median']:.1f} s (IQR {res['tau_iqr'][0]:.1f} to "
        f"{res['tau_iqr'][1]:.1f} s); bootstrap over casts, 90% interval {res['tau_ci90'][0]:.1f} to "
        f"{res['tau_ci90'][1]:.1f} s")
    if "tau_pooled" in res:
        say(f"  pooled (RSS summed over casts): {res['tau_pooled']:.1f} s")
    say(f"  the lag explains {share:.0%} of the misfit of a lag-free model (RSS, summed over these casts); "
        f"median up-minus-down misfit {good['hyst0'].median():+.3f} -> {good['hyst'].median():+.3f} {unit}")
    say("  tau is the e-folding time: the sensor covers 63% of a step change in tau seconds, 95% in 3 tau.")
    if bottles:
        res["bottles"] = bottle_check(cruise, var, res, samples, data, smooth, say)
    if write:
        out = pp.cruise_dir(cruise) / "_audit" / f"{cruise}__response_time_{var}.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        t.to_csv(out, index=False, encoding="utf-8-sig")
        if curves:
            pd.concat(curves).to_csv(out.with_name(out.stem + "_curve.csv"), index=False, encoding="utf-8-sig")
        if bottles and "bottles" in res:
            res["bottles"]["table"].to_csv(out.with_name(out.stem + "_bottles.csv"), index=False,
                                           encoding="utf-8-sig")
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
    ap.add_argument("--bottles", action="store_true",
                    help="also test the lag-corrected profiles against the bottles (vs the down-cast estimator)")
    a = ap.parse_args()
    estimate(a.cruise, a.var, a.col, a.tau_max, a.step, a.dz, a.smooth, bottles=a.bottles)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
