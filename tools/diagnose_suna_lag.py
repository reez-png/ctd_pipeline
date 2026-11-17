"""Diagnose a cruise's SUNA timing and response in one run: the clock, the kind of response, and the bottles.

    python tools/diagnose_suna_lag.py P45_07              needs step04 run first (its water entry/exit events)
    python tools/diagnose_suna_lag.py P45_07 --delay 16   a delay of your own instead of step04's
    python tools/diagnose_suna_lag.py P45_07 --no-bottles

Reads only; writes <cruise>/_audit/suna/<cruise>__suna_lag_diagnosis.csv.

1. CLOCK. step04 times when the SUNA and the CTD enter and leave the water. ENTRY is used for the delay: when the
   path fills on immersion the SUNA cannot be late for any reason but its clock. EXIT is reported beside it,
   because an optical path that keeps water after leaving the sea reads "wet" for longer and biases exit. The
   delay is the median entry offset (bootstrap 95% CI of the median, cast-to-cast drift tested).

2. RESPONSE. With that delay applied first, two first-order responses are fitted to every cast
   (suncat.nitrate.hysteresis):
     time      the sensor relaxes per SECOND: a time constant tau (s), whatever the package does;
     distance  it relaxes per METRE travelled: a flushing length L (m). Water in the optical path is replaced
               only as the package moves (an open path turned horizontal, with no pump, is one way to get this),
               so the response freezes when the package stops and its effective time constant is L / |w|.
   On one cast at one steady speed the two are the same model (tau = L / |w|); within a cast they differ only where
   the speed changes (heave, stops, different down and up speeds), so the misfit is the weaker test. What
   separates them cleanly is how they vary BETWEEN casts profiled at different speeds: a time
   response keeps tau constant and L grows with speed; a flushing response keeps L constant and tau shrinks as
   1 / speed. The test is the slope of log(tau) on log(speed) across casts (OLS, 95% t interval): 0 for a time
   response, -1 for a flushing response. The slope of log(L) is the same number plus 1, so one regression decides.
   Synthetic check (tests/test_diagnose_suna_lag.py, 4 casts at 0.33 to 1.25 m/s): a true 8 s time response gave
   tau 8.0 to 8.1 s and slope -0.00 (95% CI -0.04 to +0.04); a true 4 m flushing length gave L 4.0 m on every
   cast, tau 3.4 to 12.6 s and slope -1.00 (-1.04 to -0.95).

3. BOTTLES. The better-supported response (by the slope; by the summed misfit if the slope cannot decide) is
   tested against the bottles as in estimate_response_time --bottles.

Limits: the speeds of P45_07's casts span about a factor 2, so the slope's interval is wide; casts with tau at 0
(no response resolved) or at the grid's edge are left out of the regression; the long, deep casts are listed
separately because they have behaved differently in every fit so far.
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("estimate_response_time", _HERE / "estimate_response_time.py")
ERT = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ERT)
CME, pp = ERT.CME, ERT.pp


def _boot_median_ci(v, n_boot=5000, seed=0):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) < 3:
        return (np.nan, np.nan)
    b = np.median(np.random.default_rng(seed).choice(v, (n_boot, len(v))), axis=1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def _t_p(t, df):
    """Two-sided p of a t statistic (scipy when present, else the normal approximation, noted)."""
    try:
        from scipy import stats
        return float(2 * stats.t.sf(abs(t), df)), ""
    except ImportError:
        return float(math.erfc(abs(t) / math.sqrt(2))), " (normal approx.; scipy not installed)"


def _slope(x, y):
    """OLS y = a + b x: b, its 95% t interval, se, n."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    k = np.isfinite(x) & np.isfinite(y)
    x, y = x[k], y[k]
    if len(x) < 4 or np.ptp(x) == 0:
        return dict(b=np.nan, lo=np.nan, hi=np.nan, se=np.nan, n=len(x))
    b, a = np.polyfit(x, y, 1)
    res = y - (a + b * x)
    se = float(np.sqrt(np.sum(res ** 2) / (len(x) - 2) / np.sum((x - x.mean()) ** 2)))
    tc = CME.t975(len(x) - 2)
    return dict(b=float(b), lo=float(b - tc * se), hi=float(b + tc * se), se=se, n=len(x))


# ------------------------------------------------------------------ 1. clock
def clock(cruise: str, say=print) -> dict:
    tl = pp.cruise_dir(cruise) / "_audit" / "suna" / "time_lag"
    ev_f = tl / f"{cruise}__suna_clock_events.csv"
    if not ev_f.exists():
        say(f"  no {ev_f.name}: run step04_suna_time_lag for {cruise} first (or give --delay)")
        return {}
    ev = pd.read_csv(ev_f)
    ent, ext = ev["entry_offset_s"].astype(float), ev["exit_offset_s"].astype(float)
    out = dict(entry_median=float(ent.median()), entry_ci=_boot_median_ci(ent), n_entry=int(ent.notna().sum()),
               exit_median=float(ext.median()), exit_ci=_boot_median_ci(ext), n_exit=int(ext.notna().sum()))
    both = ev[ent.notna() & ext.notna()]
    d = (both["exit_offset_s"] - both["entry_offset_s"]).to_numpy(float)
    out.update(n_both=len(d), exit_minus_entry=float(np.median(d)) if len(d) else np.nan,
               p_sign=CME.sign_test(d) if len(d) else np.nan)
    say(f"  water entry: SUNA minus CTD median {out['entry_median']:+.1f} s (95% CI {out['entry_ci'][0]:+.1f} to "
        f"{out['entry_ci'][1]:+.1f} s, {out['n_entry']} casts)")
    say(f"  water exit : median {out['exit_median']:+.1f} s (95% CI {out['exit_ci'][0]:+.1f} to {out['exit_ci'][1]:+.1f} s, "
        f"{out['n_exit']} casts)")
    if len(d):
        say(f"  exit minus entry, same cast: median {out['exit_minus_entry']:+.1f} s over {len(d)} casts, sign test "
            f"p = {out['p_sign']:.3f}  (a path that holds water after leaving the sea makes exit later than entry)")
    per = tl / f"{cruise}__suna_time_lag_per_cast.csv"
    if per.exists():
        st = pd.read_csv(per)[["cast", "start_utc"]].merge(ev[["cast", "entry_offset_s"]], on="cast")
        h = (pd.to_datetime(st["start_utc"], format="ISO8601", utc=True)
             - pd.to_datetime(st["start_utc"], format="ISO8601", utc=True).min()).dt.total_seconds() / 3600
        r = _slope(h, st["entry_offset_s"])
        if np.isfinite(r["b"]):
            p, note = _t_p(r["b"] / r["se"], r["n"] - 2)
            out.update(drift_s_per_h=r["b"], drift_p=p)
            say(f"  drift of the entry offset through the cruise: {r['b']:+.4f} s/h (95% CI {r['lo']:+.4f} to "
                f"{r['hi']:+.4f}), p = {p:.3f}{note}, over {h.max():.0f} h")
    return out


# ------------------------------------------------------------------ 2. response
def cast_speeds(samples: dict) -> pd.DataFrame:
    rows = []
    for cast, s in samples.items():
        m = s["moving"].astype(bool) & (s["down"].astype(bool) | s["up"].astype(bool)) & s["w"].notna()
        rows.append(dict(cast=cast, speed=float(s.loc[m, "w"].abs().median()) if m.any() else np.nan,
                         max_depth=float(s["z"].max()) if s["z"].notna().any() else np.nan))
    return pd.DataFrame(rows)


def response(cruise: str, delay: float, say=print) -> dict:
    say(f"  fitting both responses with the delay {delay:+.1f} s applied first ...", flush=True)
    rt = ERT.estimate(cruise, "nitrate", delay=delay, mode="time", verbose=False, write=False)
    rd = ERT.estimate(cruise, "nitrate", delay=delay, mode="distance", verbose=False, write=False)
    T = (rt["table"][["cast", "n", "tau", "edge", "rss_min", "rss0"]]
         .rename(columns={"tau": "tau_s", "edge": "tau_edge", "rss_min": "rss_time"})
         .merge(rd["table"][["cast", "tau", "edge", "rss_min"]]
                .rename(columns={"tau": "L_m", "edge": "L_edge", "rss_min": "rss_distance"}), on="cast")
         .merge(cast_speeds(rt["samples"]), on="cast", how="left"))
    say(f"\n  {'cast':16s} {'frames':>6s} {'speed':>6s} {'depth':>6s} {'tau (s)':>8s} {'L (m)':>7s} {'L/speed':>8s}"
        f"   misfit time / distance      [speed: median |w| on the legs, m/s; L/speed: L's time constant, s]")
    for _, r in T.iterrows():
        if not np.isfinite(r["rss_time"]):
            continue
        tau_s = "EDGE" if r["tau_edge"] else f"{r['tau_s']:.1f}"
        L_s = "EDGE" if r["L_edge"] else f"{r['L_m']:.1f}"
        say(f"  {r['cast']:16s} {int(r['n']):6d} {r['speed']:6.2f} {r['max_depth']:6.0f} {tau_s:>8s} {L_s:>7s} "
            f"{r['L_m'] / r['speed'] if r['speed'] > 0 else np.nan:8.1f}   {r['rss_time'] / r['rss_distance']:.3f}")
    ok = T["tau_s"].gt(0.5) & T["L_m"].gt(0.1) & T["speed"].gt(0) & ~T["tau_edge"] & ~T["L_edge"]
    sl = _slope(np.log(T.loc[ok, "speed"]), np.log(T.loc[ok, "tau_s"]))
    both = T["rss_time"].notna() & T["rss_distance"].notna()
    rss_ratio = float(T.loc[both, "rss_time"].sum() / T.loc[both, "rss_distance"].sum())
    iqr = lambda v: (np.percentile(v, 75) - np.percentile(v, 25)) / np.median(v)
    spread_t, spread_L = float(iqr(T.loc[ok, "tau_s"])), float(iqr(T.loc[ok, "L_m"]))
    say(f"\n  slope of log(tau) on log(speed) across {sl['n']} casts: {sl['b']:+.2f} (95% CI {sl['lo']:+.2f} to "
        f"{sl['hi']:+.2f}); a time response gives 0, a flushing (per metre) response gives -1")
    say(f"  spread across casts (IQR / median): tau {spread_t:.2f}, L {spread_L:.2f}  (the right model has the "
        f"smaller spread); casts left out (tau or L at 0 or at the grid's edge): {int((~ok & both).sum())}")
    n_dist = int((T.loc[both, "rss_time"] > T.loc[both, "rss_distance"]).sum())
    say(f"  summed misfit, time / distance: {rss_ratio:.3f} (above 1 favours distance); distance fits better on "
        f"{n_dist} of {int(both.sum())} casts. Within one cast the two differ only where the speed changes (heave,")
    say("  stops, different down and up speeds), so this is the weaker test; the slope across casts is the main one.")
    has0, has1 = sl["lo"] <= 0 <= sl["hi"], sl["lo"] <= -1 <= sl["hi"]
    if has1 and not has0:
        verdict, mode = "FLUSHING: the response scales with 1 / speed (per metre travelled), not with time", "distance"
    elif has0 and not has1:
        verdict, mode = "TIME: the response does not scale with speed; a time constant describes it", "time"
    else:
        mode = "distance" if rss_ratio > 1 else "time"
        verdict = (f"UNDETERMINED by the speeds: the slope's interval holds {'both 0 and -1' if (has0 and has1) else 'neither 0 nor -1'};"
                   f" the summed misfit leans to {mode}")
    say(f"  -> {verdict}")
    med_t = float(np.median(T.loc[ok, "tau_s"])) if ok.any() else np.nan
    long_ = T[ok & (T["tau_s"] > 2 * med_t)]
    if len(long_):
        say(f"  casts with tau above twice the median ({med_t:.1f} s): "
            + ", ".join(f"{c.split('_')[-1]} ({t:.0f} s, {d:.0f} m, {v:.2f} m/s)"
                        for c, t, d, v in zip(long_["cast"], long_["tau_s"], long_["max_depth"], long_["speed"])))
        say("  (their depth, not their speed, is what they share if the speeds above are ordinary)")
    return dict(table=T, slope=sl, rss_ratio=rss_ratio, spread_tau=spread_t, spread_L=spread_L, mode=mode,
                verdict=verdict, res={"time": rt, "distance": rd})


def diagnose(cruise: str, delay: float | None = None, bottles: bool = True, say=print) -> dict:
    say(f"{cruise}: SUNA timing and response diagnosis\n\n1. CLOCK (step04's water entry and exit)")
    ck = clock(cruise, say)
    if delay is None:
        if not ck:
            raise SystemExit("no delay: run step04 for this cruise, or pass --delay")
        delay = round(ck["entry_median"] * 2) / 2
        say(f"  -> delay used: {delay:+.1f} s (median water-entry offset, to 0.5 s)")
    else:
        say(f"  -> delay used: {delay:+.1f} s (given)")
    say("\n2. RESPONSE (time constant or flushing length?)")
    rp = response(cruise, delay, say)
    out = dict(clock=ck, delay=delay, response=rp)
    if bottles:
        say(f"\n3. BOTTLES, with the {rp['mode']} response and the {delay:+.1f} s delay")
        r = rp["res"][rp["mode"]]
        r["delay_fixed"] = delay
        out["bottles"] = ERT.bottle_check(cruise, "nitrate", r, r["samples"], r["data"], ERT.SMOOTH, say)
    f = pp.cruise_dir(cruise) / "_audit" / "suna" / f"{cruise}__suna_lag_diagnosis.csv"
    f.parent.mkdir(parents=True, exist_ok=True)
    rp["table"].assign(delay_s=delay).to_csv(f, index=False, encoding="utf-8-sig")
    say(f"\nper-cast table -> {f}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cruise")
    ap.add_argument("--delay", type=float, default=None, help="delay, s (default: step04's median water-entry offset)")
    ap.add_argument("--no-bottles", action="store_true")
    a = ap.parse_args()
    diagnose(a.cruise, a.delay, not a.no_bottles)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
