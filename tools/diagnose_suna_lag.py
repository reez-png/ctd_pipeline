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

3. START OF CAST. Is the long tau of a few casts the sensor being slow throughout, or a transient at
   immersion? Each cast's time out of the water before it (from the end of the previous cast's frames) and
   whether it was the day's first cast are set against tau. Then each cast is refitted with a transient
   a * exp(-t'/T0) beside the lag (t' from the first frame in the water; T0 of 30, 60, 120 or 240 s, the best
   kept). A transient (fresh water or air left in the path, warm-up after power-on) decays from immersion and
   falls mostly on the early down-cast; a lag shows on both legs all through the cast, so the two separate.
   Leaving the first minutes out instead does not work: it removes the down-cast's crossing of a shallow
   nitracline, the part of the cast that carries the lag. Synthetic check (tests/test_diagnose_suna_lag.py): a
   6 s sensor with a -6 umol/L transient over 60 s fits as a 10.4 s lag alone and as 6.0 s and -5.9 with the
   transient; a true 25 s sensor keeps 25 s with the transient allowed.

4. BOTTLES. The better-supported response (by the slope; by the summed misfit if the slope cannot decide) is
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


def cast_conditions(cruise: str) -> pd.DataFrame:
    """Per cast, from the nitrate frame files: temperature range on the legs (5th to 95th percentile, degC, so one
    bad frame does not set it), coldest temperature, and how long the legs took (s)."""
    rows = []
    for f in sorted(pp.product_dir(cruise, "nitrate_qc").glob("*_nitrate_QC.csv")):
        d = pd.read_csv(f, encoding="utf-8-sig")
        leg = d["leg"].astype(str).str.lower().isin(["down", "up"]) if "leg" in d else pd.Series(True, index=d.index)
        tc = pd.to_numeric(d.get("temp_c"), errors="coerce")[leg].dropna()
        tt = CME.parse_utc(d["utc_time"], "frame utc_time")[leg].dropna()
        rows.append(dict(cast=f.name.replace("_nitrate_QC.csv", ""),
                         t_range=float(np.percentile(tc, 95) - np.percentile(tc, 5)) if len(tc) > 5 else np.nan,
                         t_min=float(np.percentile(tc, 5)) if len(tc) > 5 else np.nan,
                         duration=float((tt.max() - tt.min()).total_seconds()) if len(tt) > 1 else np.nan))
    return pd.DataFrame(rows)


def spearman(x, y, n_perm=10000, seed=0):
    """Spearman rank correlation and a two-sided permutation p-value (no scipy needed)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    k = np.isfinite(x) & np.isfinite(y)
    x, y = x[k], y[k]
    if len(x) < 5:
        return np.nan, np.nan, len(x)
    rx, ry = pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy()
    rho = float(np.corrcoef(rx, ry)[0, 1])
    rng = np.random.default_rng(seed)
    perm = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(n_perm)])
    return rho, float((np.sum(np.abs(perm) >= abs(rho)) + 1) / (n_perm + 1)), len(x)


def conditions_test(T: pd.DataFrame, say=print) -> dict:
    """Which cast property goes with a long tau? Spearman correlations of tau (time model) with speed, depth,
    temperature range, coldest temperature and duration across casts, and the correlations among those
    properties, because they travel together (deep casts are long, cold and, on P45_07, slow). Observational: a
    correlation here says which candidate the data favour, not what causes the long tau."""
    ok = T["tau_s"].gt(0.5) & ~T["tau_edge"]
    if ok.sum() < 5:
        say(f"\n  cast properties vs tau: too few casts with a tau ({int(ok.sum())}; 5 needed)")
        return {}
    props = [("speed", "speed (m/s)"), ("max_depth", "max depth (m)"), ("t_range", "temperature range (degC)"),
             ("t_min", "coldest temperature (degC)"), ("duration", "duration of the legs (s)"),
             ("idle_h", "time out of the water (h)")]
    out = {}
    say(f"\n  which cast property goes with a long tau? Spearman rank correlation with tau across {int(ok.sum())} casts"
        " (permutation p):")
    for c, lab in props:
        if c in T:
            rho, p, n = spearman(T.loc[ok, c], T.loc[ok, "tau_s"])
            out[c] = dict(rho=rho, p=p, n=n)
            say(f"    {lab:28s} rho {rho:+.2f}  p = {p:.3f}")
    pairs = [("max_depth", "t_range"), ("max_depth", "duration"), ("max_depth", "speed"), ("t_range", "duration"),
             ("max_depth", "idle_h")]
    say("  how those properties go together (rho): " + ", ".join(
        f"{a} / {b} {spearman(T.loc[ok, a], T.loc[ok, b], n_perm=0)[0]:+.2f}" for a, b in pairs if a in T and b in T))
    say("  (observational: these move together across casts, so the strongest correlate is a lead, not a cause)")
    return out


def idle_times(samples: dict) -> pd.DataFrame:
    """Per cast: hours since the previous cast's frames ended (NaN for the first cast) and whether it is the
    first cast of its UTC day."""
    rows = []
    for cast, s in samples.items():
        t = pd.to_datetime(s["t"], utc=True).dropna()
        if len(t):
            rows.append(dict(cast=cast, start=t.min(), end=t.max()))
    if not rows:
        return pd.DataFrame(columns=["cast", "idle_h", "first_of_day"])
    I = pd.DataFrame(rows).sort_values("start").reset_index(drop=True)
    I["idle_h"] = (I["start"] - I["end"].shift()).dt.total_seconds() / 3600
    I["first_of_day"] = I["start"].dt.date.ne(I["start"].dt.date.shift())
    return I[["cast", "idle_h", "first_of_day"]]


TRANSIENT_T0 = (30.0, 60.0, 120.0, 240.0)   # s: e-folding times tried for a start-of-cast transient


def _fit_transient(t, zd, y, grid, tau, e, smooth):
    """The lag model plus one extra column e (a transient, amplitude a): least squares with the curvature penalty
    on the profile only. Returns (rss, N, a). e = None fits the lag model alone."""
    A = ERT.H.lag_design(t, zd, grid, tau)
    ok = np.isfinite(y)
    G = len(grid)
    X = A[ok] if e is None else np.hstack([A[ok], e[ok, None]])
    P = np.zeros((X.shape[1], X.shape[1]))
    D = np.zeros((G - 2, G))
    i = np.arange(G - 2)
    D[i, i], D[i, i + 1], D[i, i + 2] = 1.0, -2.0, 1.0
    dz = float(np.median(np.diff(grid)))
    P[:G, :G] = (smooth / dz ** 4) * (D.T @ D)
    M = X.T @ X + P
    M[np.diag_indices(len(M))] += 1e-9 * max(float(np.trace(M)) / len(M), 1e-12)
    b = np.linalg.solve(M, X.T @ y[ok])
    res = y[ok] - X @ b
    return float(np.sum(res ** 2)), b[:G], (float(b[G]) if e is not None else 0.0)


def start_of_cast(s: pd.DataFrame, taus, dz: float = 1.0, delay: float = 0.0, t0s=TRANSIENT_T0,
                  smooth: float = ERT.SMOOTH) -> dict:
    """tau fitted with a start-of-cast transient beside the lag: the reading carries an extra a * exp(-t'/T0),
    t' the time since the first frame with a value on a leg. A transient (water or air left in the path,
    warm-up after power-on) decays with time from immersion and touches mainly the early down-cast; a lag shows
    on both legs all through the cast, so the two separate. Returns tau_tr (s), the transient's amplitude a at
    immersion (value units) and e-folding T0 (s), and the RMS misfit without and with it."""
    s = s[s["t"].notna() & s["z"].notna()].sort_values("t")
    out = dict(tau_lag=np.nan, tau_tr=np.nan, tr_amp=np.nan, tr_t0=np.nan, rms_lag=np.nan, rms_tr=np.nan)
    t = ERT.elapsed_seconds(s["t"])
    z, y = s["z"].to_numpy(float), s["v"].to_numpy(float)
    dn, up = s["down"].to_numpy(bool), s["up"].to_numpy(bool)
    base = ERT.H.cast_tau(t, z, y, dn, up, taus, dz, smooth, delay=delay)
    if "skip" in base or not np.isfinite(base.get("tau_grid", np.nan)):
        return out
    tt, zd, yk, grid = base["data"]
    on = np.isfinite(yk)
    since = np.clip(tt - tt[on][0], 0.0, None)
    coarse = taus[::3] if len(taus) > 30 else taus
    best = (np.inf, None, None)
    for T0 in t0s:
        e = np.exp(-since / T0)
        for tau in coarse:
            r = _fit_transient(tt, zd, yk, grid, tau, e, smooth)[0]
            if r < best[0]:
                best = (r, tau, T0)
    _, tau_c, T0 = best
    e = np.exp(-since / T0)
    step = float(np.median(np.diff(taus)))
    fine = np.arange(max(0.0, tau_c - 3 * step), min(taus[-1], tau_c + 3 * step) + step / 4, step / 2)
    rss = np.array([_fit_transient(tt, zd, yk, grid, tau, e, smooth)[0] for tau in fine])
    tau_tr, _ = ERT.H.refine_min(fine, rss)
    r_tr, _, a = _fit_transient(tt, zd, yk, grid, tau_tr, e, smooth)
    n = int(on.sum())
    out.update(tau_lag=float(base["tau"]), tau_tr=float(tau_tr), tr_amp=a, tr_t0=float(T0), rms_lag=float(base["rms"]),
               rms_tr=float(np.sqrt(r_tr / n)))
    return out


def start_test(T: pd.DataFrame, samples: dict, delay: float, say=print) -> dict:
    """Section 3: time out of the water, first cast of the day, and tau refitted with a start-of-cast transient."""
    taus = np.arange(0.0, 90.5, 1.0)
    say("  refitting each cast with a transient at immersion beside the lag ...", flush=True)
    S = pd.DataFrame([dict(cast=c, **start_of_cast(samples[c], taus, delay=delay)) for c in T["cast"] if c in samples])
    T = T.drop(columns=[c for c in S.columns if c != "cast" and c in T], errors="ignore").merge(S, on="cast", how="left")
    say(f"\n  {'cast':16s} {'idle h':>7s} {'1st':>4s} {'tau (s)':>8s} {'tau_tr':>7s} {'amp':>7s} {'T0 s':>5s}"
        f"   RMS lag -> +transient   [tau_tr: tau fitted with a transient amp * exp(-t / T0) from immersion; amp in "
        f"the value's units]")
    for _, r in T.iterrows():
        if not np.isfinite(r.get("rss_time", np.nan)) or not np.isfinite(r.get("tau_tr", np.nan)):
            continue
        say(f"  {r['cast']:16s} {r['idle_h']:7.1f} {'yes' if r['first_of_day'] else '':>4s} {r['tau_s']:8.1f} "
            f"{r['tau_tr']:7.1f} {r['tr_amp']:+7.2f} {r['tr_t0']:5.0f}   {r['rms_lag']:.3f} -> {r['rms_tr']:.3f}")
    ok = T["tau_s"].gt(0.5) & ~T["tau_edge"]
    out = dict(table=T)
    rho, p, n = spearman(T.loc[ok, "idle_h"], T.loc[ok, "tau_s"])
    out["idle_tau"] = dict(rho=rho, p=p, n=n)
    say(f"\n  time out of the water vs tau: " + (f"Spearman rho {rho:+.2f}, p = {p:.3f}, {n} casts" if np.isfinite(rho)
                                                else f"too few casts with an idle time ({n}; 5 needed)"))
    # depth travels with tau and, through the day's station order, with idle time: hold it fixed
    k = ok & T["max_depth"].gt(0) & T["idle_h"].notna()
    if k.sum() >= 5:
        lx, ly = np.log(T.loc[k, "max_depth"]), np.log(T.loc[k, "tau_s"])
        b, a = np.polyfit(lx, ly, 1)
        rho_p, p_p, n_p = spearman(T.loc[k, "idle_h"], ly - (a + b * lx))
        rho_d = spearman(T.loc[k, "idle_h"], T.loc[k, "max_depth"], n_perm=0)[0]
        out["idle_tau_given_depth"] = dict(rho=rho_p, p=p_p, n=n_p)
        say(f"  ... with depth held fixed (log tau on log depth, slope {b:+.2f}; rho with the remainder): rho "
            f"{rho_p:+.2f}, p = {p_p:.3f}  (idle time vs depth: rho {rho_d:+.2f})")
    fd = T["first_of_day"].fillna(False).astype(bool)
    for lab, col in (("tau", "tau_s"), ("tau with a transient", "tau_tr"), ("transient amplitude", "tr_amp")):
        a_, b_ = T.loc[ok & fd, col].dropna(), T.loc[ok & ~fd, col].dropna()
        if len(a_) and len(b_):
            say(f"  {lab:22s} median: first cast of the day {a_.median():+6.2f} ({len(a_)} casts), others "
                f"{b_.median():+6.2f} ({len(b_)} casts)")
    ea = T.loc[ok, ["idle_h", "tr_amp"]].dropna()
    if len(ea) >= 5:
        rho_e, p_e, n_e = spearman(ea["idle_h"], ea["tr_amp"])
        out["idle_amp"] = dict(rho=rho_e, p=p_e, n=n_e)
        say(f"  transient amplitude vs time out of the water: rho {rho_e:+.2f}, p = {p_e:.3f}, {n_e} casts")
    med = float(T.loc[ok, "tau_s"].median()) if ok.any() else np.nan
    long_ = T[ok & (T["tau_s"] > 2 * med) & T["tau_tr"].notna()]
    if len(long_):
        say(f"  the long casts (tau above {2 * med:.1f} s), with a transient allowed:")
        for _, r in long_.iterrows():
            say(f"    {r['cast'].split('_')[-1]}: tau {r['tau_s']:.1f} -> {r['tau_tr']:.1f} s, transient {r['tr_amp']:+.2f} "
                f"(T0 {r['tr_t0']:.0f} s), RMS {r['rms_lag']:.3f} -> {r['rms_tr']:.3f}  (idle {r['idle_h']:.1f} h"
                f"{', first of the day' if r['first_of_day'] else ''})")
        ratio = float((long_["tau_tr"] / long_["tau_s"]).median())
        out["long_ratio"] = ratio
        if ratio <= 0.5:
            v = (f"TRANSIENT AT IMMERSION: on the long casts a start transient takes over most of the lag (tau falls to "
                 f"{ratio:.0%} of itself, median). A one-off transient fits them better than a slow sensor")
        elif ratio >= 0.75:
            v = (f"NOT A START TRANSIENT: on the long casts tau keeps {ratio:.0%} of itself (median) with a transient "
                 "allowed; the slow response lasts through the cast")
        else:
            v = f"PARTLY: on the long casts tau keeps {ratio:.0%} of itself (median) with a transient allowed"
        out["verdict"] = v
        say(f"  -> {v}")
    say("  (observational: idle time is a proxy for powering off, rinsing and drying between casts, and on this "
        "cruise the deep stations were often the day's first)")
    return out


def response(cruise: str, delay: float, say=print) -> dict:
    say(f"  fitting both responses with the delay {delay:+.1f} s applied first ...", flush=True)
    rt = ERT.estimate(cruise, "nitrate", delay=delay, mode="time", verbose=False, write=False)
    rd = ERT.estimate(cruise, "nitrate", delay=delay, mode="distance", verbose=False, write=False)
    T = (rt["table"][["cast", "n", "tau", "edge", "rss_min", "rss0"]]
         .rename(columns={"tau": "tau_s", "edge": "tau_edge", "rss_min": "rss_time"})
         .merge(rd["table"][["cast", "tau", "edge", "rss_min"]]
                .rename(columns={"tau": "L_m", "edge": "L_edge", "rss_min": "rss_distance"}), on="cast")
         .merge(cast_speeds(rt["samples"]), on="cast", how="left")
         .merge(cast_conditions(cruise), on="cast", how="left")
         .merge(idle_times(rt["samples"]), on="cast", how="left"))
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
    cond = conditions_test(T, say)
    return dict(table=T, slope=sl, rss_ratio=rss_ratio, spread_tau=spread_t, spread_L=spread_L, mode=mode,
                verdict=verdict, res={"time": rt, "distance": rd}, conditions=cond)


def diagnose(cruise: str, delay: float | None = None, bottles: bool = True, say=print) -> dict:
    say(f"{cruise}: SUNA timing and response diagnosis\n\n1. CLOCK (step04's water entry and exit)")
    ck = clock(cruise, say)
    applied = ERT.applied_clock(cruise)
    a_off = float(np.median(list(applied.values()))) if applied else 0.0
    if delay is None:
        if not ck:
            raise SystemExit("no delay: run step04 for this cruise, or pass --delay")
        # step04's events are measured on the raw (uncorrected) co-located product, so they are the RAW clock;
        # if step02/step03 already shifted the frames, only the remainder is applied here
        delay = round((ck["entry_median"] + a_off) * 2) / 2
        if applied:
            say(f"  the nitrate frames already carry step02's clock correction ({a_off:+.1f} s, {len(applied)} casts);"
                f" remaining delay {ck['entry_median']:+.1f} {a_off:+.1f} = {ck['entry_median'] + a_off:+.1f} s")
        say(f"  -> delay used: {delay:+.1f} s ({'remainder after the applied correction' if applied else 'median water-entry offset'}, to 0.5 s)")
    else:
        say(f"  -> delay used: {delay:+.1f} s (given)"
            + (f"; NOTE the frames already carry {a_off:+.1f} s from step02, so this is on top of it" if applied else ""))
    say("\n2. RESPONSE (time constant or flushing length?)")
    rp = response(cruise, delay, say)
    say("\n3. START OF CAST (a transient at immersion, or a slow sensor?)")
    st = start_test(rp["table"], rp["res"]["time"]["samples"], delay, say)
    rp["table"] = st["table"]
    out = dict(clock=ck, delay=delay, response=rp, start=st)
    if bottles:
        say(f"\n4. BOTTLES, with the {rp['mode']} response and the {delay:+.1f} s delay")
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
