"""Sensor response time (tau) from each cast's own profiles: a first-order lag fitted to the down-cast and the up-cast
together, without bottles.

    python tools/estimate_response_time.py P45_07                      SUNA nitrate (default)
    python tools/estimate_response_time.py P45_07 --var chl            chlorophyll fluorescence (4 Hz CTD scans)
    python tools/estimate_response_time.py P45_07 --col no3_stage1     another nitrate column
    python tools/estimate_response_time.py P45_07 --tau-max 120        a longer grid
    python tools/estimate_response_time.py P45_07 --bottles            also: do lag-corrected profiles match the
                                                                       bottles better than the down-cast?
    python tools/estimate_response_time.py P45_07 --delay 16           a known pure delay first (e.g. step04's
                                                                       water-entry clock), then tau
    python tools/estimate_response_time.py P45_07 --fit-delay          delay and tau fitted together

Delay and response. A pure delay d (a clock offset, or water travelling through tubing) shifts the record; a
first-order response also smooths it. The first-order model alone absorbs a delay into tau (synthetic: 16 s delay
+ 8 s response fitted as tau of 25 to 49 s). --delay applies a known d first (synthetic: d = 16 given, tau 8.1 s
recovered; a d wrong by 2 s moves tau by about 2 s the other way, so delay + tau is conserved). --fit-delay
searches both: on the synthetic cruise (noise 0.5 umol/L, heave and stops) 16 s + 8 s came back as 16.3 s and
7.6 s; on plain casts with no stops, a response under about half the delay can be fitted as tau = 0 and a longer
delay, with delay + tau still right. The split then needs an independent delay (step04's water-entry clock).

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
import suncat.nitrate.bottle_match as BM  # noqa: E402
from suncat.nitrate.hysteresis import (SMOOTH, MIN_FRAMES, lag_design, fit_profile, refine_min,  # noqa: E402,F401
                                       search_tau, summarise)
from suncat.timeutil import elapsed_seconds  # noqa: E402

DEFAULTS = {"nitrate": dict(tau_max=90.0, step=1.0, dz=1.0), "chl": dict(tau_max=12.0, step=0.25, dz=0.5)}


# ------------------------------------------------------------------ the model (suncat.nitrate.hysteresis)
def cast_tau(s: pd.DataFrame, taus: np.ndarray, dz: float, smooth: float = SMOOTH, delay: float = 0.0,
             delays=None, mode: str = "time") -> dict:
    """tau for one cast's samples (t, z, v, down, up as compare_match_estimators reads them), with a fixed pure
    delay (default none) or, given a grid of `delays`, the delay and tau fitted together. Residuals only on the
    legs (down: Loop Edit kept; up): surface soak and flagged scans carry the sensor's state forward in time but
    are not fitted."""
    s = s[s["t"].notna() & s["z"].notna()]
    a = (elapsed_seconds(s["t"]), s["z"].to_numpy(float), s["v"].to_numpy(float),
         s["down"].to_numpy(bool), s["up"].to_numpy(bool))
    if delays is not None:
        return H.cast_delay_tau(*a, delays=delays, taus=taus, dz=dz, smooth=smooth, mode=mode)
    return H.cast_tau(*a, taus, dz, smooth, delay=delay, mode=mode)


# ------------------------------------------------------------------ against the bottles
def _adopted(cruise):
    """The cruise's stated SUNA calibration (suna_calibration.json), or None."""
    try:
        from suncat.suna_cal import load_suna_calibration
        return load_suna_calibration(pp.cruise_dir(cruise))
    except Exception:
        return None


def bottle_check(cruise, var, res, samples, data, smooth=SMOOTH, say=print):
    """Does correcting the lag make the sensor match the bottles better? On the same bottles, paired by the one
    rule in suncat.nitrate.bottle_match (shared with step02b):
      down          step07's down-cast cell median over the moving down-cast frames
      N@0           the profile fitted with no lag (both legs averaged, in effect)
      N@tau_cruise  the lag-corrected profile with the cruise's median tau (one constant for the sensor)
    each read at the bottle depth, then sensor = offset + gain * bottle by OLS, LOOCV RMSE in bottle units, and
    paired LOOCV squared errors against 'down' (exact sign test). N@tau_cast, each cast's own tau, is scored on
    its own line over the bottles it has a value for.
    Then the depth offset (offset_sweep): the lag-corrected profile carries no lag, so reading it at the bottle depth
    + D and scoring each D measures any offset between where the bottle's water came from and the depth recorded
    for it (rosette geometry, firing timing, bottle closure), free of the lag that confounds the same sweep on raw
    legs."""
    unit = "mg/m3" if var == "chl" else "umol/L"
    if var == "chl":
        bottles = CME.load_bottles(cruise, var)
        src = "step07's chlorophyll calibration table"
    else:
        bottles = BM.load_bottles(pp.cruise_dir(cruise) / "metadata" / "bottle_nitrate" / f"{cruise}_bottle_nitrate.csv")
        src = f"{cruise}_bottle_nitrate.csv"
    tau_c = res["tau_median"]
    d_c = res.get("delay_median", res.get("delay_fixed", 0.0))
    per_cast = res["table"].set_index("cast")["tau"]
    md = res.get("mode", "time")
    prof, frames = {}, {}
    for cast, (t, zd, y, grid, z_raw) in data.items():
        tc = per_cast.get(cast, np.nan)
        prof[cast] = {"grid": grid, "N@0": fit_profile(t, z_raw, y, grid, 0.0, smooth)[1],
                      "N@tau_cruise": fit_profile(t, z_raw, y, grid, tau_c, smooth, d_c, md)[1],
                      "N@tau_cast": fit_profile(t, zd, y, grid, tc, smooth, mode=md)[1] if np.isfinite(tc) else None}
        prof[cast]["nc"] = prof[cast]["N@tau_cruise"]                       # offset_sweep reads it
    for cast, s in samples.items():
        s = s[s["t"].notna() & s["z"].notna()].sort_values("t")
        t_s = elapsed_seconds(s["t"])
        frames[cast] = dict(z=s["z"].to_numpy(float), y=s["v"].to_numpy(float),
                            down=s["down"].to_numpy(bool) & BM.moving(t_s, s["z"].to_numpy(float)))
    names = ["N@0", "N@tau_cruise"]
    T = BM.pair(bottles, frames, prof, names)
    say(f"\n  against the bottles ({src}); tau_cruise = {tau_c:.1f} "
        f"{'m (flushing length)' if md == 'distance' else 's'}" + (f", delay {d_c:+.1f} s" if d_c else "")
        + " (N@0: no delay, no lag)")
    cmp = BM.compare(T, names)
    BM.report(cmp, "step03's nitrate frames (merged on the CTD 1 s record; down-cast: scans Loop Edit kept, moving)",
              unit, say)
    out = dict(table=T, fits=cmp["fits"], loocv_rmse={e: f["loocv"] for e, f in cmp["fits"].items()},
               paired={e: dict(median=f["median_d"], better=f["better"], n=f["n"], p_sign=f["p_sign"])
                       for e, f in cmp["fits"].items() if "p_sign" in f}, n=cmp["n"], left_out=cmp["left_out"])
    if not cmp["fits"]:
        return out
    Tk = BM.pair(bottles, frames, prof, ["N@tau_cast"])
    ck = BM.compare(Tk, ["N@tau_cast"])
    if "N@tau_cast" in ck["fits"]:
        f = ck["fits"]["N@tau_cast"]
        out["fits"]["N@tau_cast"], out["loocv_rmse"]["N@tau_cast"] = f, f["loocv"]
        say(f"  {'N@tau_cast':13s} {f['gain']:7.4f} ({f['lo']:.4f} to {f['hi']:.4f}) {f['offset']:+8.3f} {f['R2']:6.3f}  "
            f"{f['loocv']:.3f}   (each cast's own tau; {ck['n']} bottles, against down on the same: better on "
            f"{f['better']} of {f['n']}, sign test p = {f['p_sign']:.3f})")
    out["offset"] = offset_sweep(cmp["both"], prof, unit, say=say)
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
    if np.isfinite(ci[0]) and (ci[1] - ci[0]) >= 0.5 * (deltas.max() - deltas.min()):
        say("  NOT DETERMINED: the bootstrap interval covers half the range or more, so these bottles do not fix an")
        say("  offset (a curve with two optima, or one carried by a few bottles). Read the best D as no estimate.")
    return dict(table=T, determined=not (np.isfinite(ci[0]) and (ci[1] - ci[0]) >= 0.5 * (deltas.max() - deltas.min())),
                best_D=float(T.at[i, "D"]), q2_best=float(T.at[i, "q2"]), q2_at_0=float(T.at[j0, "q2"]),
                near=(float(near.min()), float(near.max())), ci90=ci)


# ------------------------------------------------------------------ a cruise
def applied_clock(cruise: str) -> dict:
    """{cast: offset_s} that step02/step03 already applied to the nitrate frames' times (the provenance column
    suna_time_offset_s, written only when an offset was applied). A delay given here on top of it would count the
    clock twice."""
    out = {}
    for f in sorted(pp.product_dir(cruise, "nitrate_qc").glob("*_nitrate_QC.csv")):
        try:
            head = pd.read_csv(f, nrows=1, encoding="utf-8-sig")
        except Exception:
            continue
        if "suna_time_offset_s" in head.columns and pd.notna(head["suna_time_offset_s"].iloc[0]):
            out[f.name.replace("_nitrate_QC.csv", "")] = float(head["suna_time_offset_s"].iloc[0])
    return out


def read_cast(cruise: str, var: str, cast_file: Path, col: str | None):
    if var == "chl":
        return CME.read_ctd_scans(cast_file)
    return CME.read_suna_frames_df(pd.read_csv(cast_file), vcol=col)


def estimate(cruise: str, var: str = "nitrate", col: str | None = None, tau_max: float | None = None,
             step: float | None = None, dz: float | None = None, smooth: float = SMOOTH, n_boot: int = 2000,
             verbose: bool = True, write: bool = True, bottles: bool = False, delay: float = 0.0,
             fit_delay: bool = False, delay_min: float = -10.0, delay_max: float = 40.0, mode: str = "time") -> dict:
    say = print if verbose else (lambda *a, **k: None)
    d = dict(DEFAULTS[var])
    if mode == "distance":                                 # tau is a flushing length in metres
        d.update(tau_max=40.0, step=0.5)
    tu = "m" if mode == "distance" else "s"
    tau_max, step, dz = tau_max or d["tau_max"], step or d["step"], dz or d["dz"]
    taus = np.arange(0.0, tau_max + step / 2, step)
    folder, pattern = (pp.product_dir(cruise, "ctd_derived"), "*_der.cnv") if var == "chl" else \
                      (pp.product_dir(cruise, "nitrate_qc"), "*_nitrate_QC.csv")
    files = sorted(folder.glob(pattern))
    if not files:
        raise FileNotFoundError(f"no {pattern} in {folder}")
    unit = "mg/m3" if var == "chl" else "umol/L"
    say(f"{cruise} {var}: first-order lag fitted to each cast's down and up profiles together (no bottles)")
    say(f"  tau grid 0 to {tau_max:g} {tu}, step {step:g} {tu}; profile grid {dz:g} m; curvature penalty {smooth:g}"
        + (f"; value column {col}" if col else ""))
    if mode == "distance":
        say("  RESPONSE PER METRE TRAVELLED: tau is a flushing length L (m); the sensor relaxes as the package moves "
            "and is frozen while it stops (effective time constant L / |w|)")
    delays = np.arange(delay_min, delay_max + 1.0 / 2, 1.0) if fit_delay else None
    applied = applied_clock(cruise) if var == "nitrate" else {}
    if applied:
        vals = sorted(set(applied.values()))
        say(f"  NOTE the frames already carry a clock correction from step02 ({vals[0]:+g}"
            + (f" to {vals[-1]:+g}" if len(vals) > 1 else "") + f" s, {len(applied)} casts): any delay given or fitted "
            "here is what remains on top of it")
        if delay:
            say(f"  WARNING --delay {delay:+g} s on corrected frames counts the clock twice unless it is a residual")
    if fit_delay:
        say(f"  DELAY AND RESPONSE fitted together: delay grid {delay_min:+g} to {delay_max:+g} s (a pure delay shifts "
            "the record; a response also smooths it)")
    elif delay:
        say(f"  FIXED DELAY {delay:+g} s applied first (frames carry the water from {delay:g} s earlier), then tau fitted")
    rows, curves, data, samples, full = [], [], {}, {}, {}
    dcol = f" {'delay':>6s}" if fit_delay else ""
    say(f"\n  {'cast':16s} {'frames':>6s} {'down/up':>9s}{dcol} {'tau (' + tu + ')':>8s}   RMS misfit {unit}: none -> fitted   "
        f"up-minus-down misfit: none -> fitted" + ("   [tau with no delay]" if fit_delay else ""))
    for f in files:                                       # each row printed as it is done
        cast = f.name.replace("_der.cnv", "").replace("_nitrate_QC.csv", "")
        samples[cast] = read_cast(cruise, var, f, col)
        r = cast_tau(samples[cast], taus, dz, smooth, delay=delay, delays=delays, mode=mode)
        rows.append({k: v for k, v in dict(cast=cast, **r).items()
                     if k not in ("curve", "data", "z_raw", "curve_delay")})
        if r["curve"] is not None:
            curves.append(pd.DataFrame({"cast": cast, "tau_s": taus, "rss": r["curve"]}))
            data[cast] = (*r["data"], r["z_raw"])
            full[cast] = r
        if isinstance(r.get("skip"), str):
            say(f"  {cast:16s} skipped: {r['skip']}", flush=True)
            continue
        tau_s = "  EDGE" if r["edge"] else f"{r['tau']:8.1f}"
        d_s = ("" if not fit_delay else (f" {'EDGE':>6s}" if r.get("delay_edge") else f" {r['delay']:+6.1f}"))
        say(f"  {cast:16s} {r['n']:6d} {r['n_down']:4d}/{r['n_up']:<4d}{d_s} {tau_s:>8s}   "
            f"{r['rms0']:8.3f} -> {r['rms']:6.3f}            {r['hyst0']:+8.3f} -> {r['hyst']:+7.3f}"
            + (f"   [{r['tau_only']:.1f}]" if fit_delay else ""), flush=True)
    t = pd.DataFrame(rows)
    good = t[t["tau"].notna()]
    res = dict(table=t, taus=taus, mode=mode, unit_tau=tu, samples=samples, data=data, full=full,
               applied_clock=applied)
    n_edge = int(t["edge"].sum())
    if n_edge:
        say(f"\n  {n_edge} cast(s) at the end of the grid (EDGE): no minimum below {tau_max:g} {tu}; left out. "
            "If many, widen --tau-max.")
    if len(good) < 3:
        say("  fewer than 3 casts with an estimate: no cruise value")
        return res
    if fit_delay:                                      # a delay at the edge of its grid is no estimate either
        full = {c: r for c, r in full.items() if not r.get("delay_edge")}
    res.update(summarise(full, taus, smooth, n_boot))
    if delay and not fit_delay:
        res["delay_fixed"] = float(delay)
    share = res["rss_share_explained"]
    say(f"\n  tau over {res['n_casts']} casts: median {res['tau_median']:.1f} {tu} (IQR {res['tau_iqr'][0]:.1f} to "
        f"{res['tau_iqr'][1]:.1f} {tu}); bootstrap over casts, 90% interval {res['tau_ci90'][0]:.1f} to "
        f"{res['tau_ci90'][1]:.1f} {tu}")
    if "tau_pooled" in res:
        say(f"  pooled (RSS summed over casts): {res['tau_pooled']:.1f} {tu}")
    if "delay_fixed" in res:
        say(f"  delay fixed at {res['delay_fixed']:+.1f} s (given, not estimated)"
            + (f"; delay + tau: median {res['delay_fixed'] + res['tau_median']:.1f} s" if mode == "time" else ""))
    elif "delay_median" in res:
        say(f"  delay over these casts: median {res['delay_median']:+.1f} s (IQR {res['delay_iqr'][0]:+.1f} to "
            f"{res['delay_iqr'][1]:+.1f} s); bootstrap over casts, 90% interval {res['delay_ci90'][0]:+.1f} to "
            f"{res['delay_ci90'][1]:+.1f} s;  delay + tau: median {res['sum_median']:.1f} s")
        if res.get("n_tau_zero"):
            say(f"  NOTE {res['n_tau_zero']} of {res['n_casts']} casts put tau at 0: there the fit could not tell a short "
                "response from a slightly longer delay (on synthetic casts this happens when tau is under about half")
            say("  the delay at noise ~0.5 umol/L; the sum delay + tau is then still about right, the split is not).")
    say(f"  the {'delay and lag explain' if (fit_delay or delay) else 'lag explains'} {share:.0%} of the misfit of a "
        f"model with neither (RSS, summed over these casts); "
        f"median up-minus-down misfit {good['hyst0'].median():+.3f} -> {good['hyst'].median():+.3f} {unit}")
    say("  tau is the e-folding " + ("distance: the sensor covers 63% of a step change in tau metres of travel, 95% in 3 tau."
        if mode == "distance" else "time: the sensor covers 63% of a step change in tau seconds, 95% in 3 tau."))
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
    ap.add_argument("--delay", type=float, default=0.0,
                    help="fixed pure delay, s, applied before tau is fitted (e.g. the water-entry clock from step04)")
    ap.add_argument("--fit-delay", action="store_true", help="fit a pure delay and tau together")
    ap.add_argument("--response", choices=["time", "distance"], default="time",
                    help="relax per second (time constant, s) or per metre travelled (flushing length, m)")
    ap.add_argument("--delay-range", type=float, nargs=2, default=(-10.0, 40.0), metavar=("MIN", "MAX"),
                    help="delay grid for --fit-delay, s (default -10 40)")
    ap.add_argument("--bottles", action="store_true",
                    help="also test the lag-corrected profiles against the bottles (vs the down-cast estimator)")
    a = ap.parse_args()
    estimate(a.cruise, a.var, a.col, a.tau_max, a.step, a.dz, a.smooth, bottles=a.bottles, delay=a.delay,
             fit_delay=a.fit_delay, delay_min=a.delay_range[0], delay_max=a.delay_range[1], mode=a.response)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
