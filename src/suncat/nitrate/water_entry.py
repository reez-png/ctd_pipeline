"""When did the SUNA enter the water, by its own optics, compared with when the CTD did? (step04's clock check)

The CTD's entry is its first wet scan (conductivity above 1), from suncat.nitrate.spectral_qc.ctd_wet_intervals. The
SUNA's entry is read from its peak counts (and abs254 where present), not from step01's f_outofwater flag: on deck
the SUNA sits in one of two states, a dry optical path (peak about 26,000 counts on P45_07) or a wet one (about
38,000 to 42,000), and seawater reads about 45,000. Until Nov 2026 step04 used the flag, which a file's own window
could not set right when deck time dominated the file (P45_06, P45_07 file A0000027), so the check failed there.

Per cast, with times on each instrument's own clock (SUNA minus CTD, so a positive offset means the SUNA stamps late):
  water level   median peak (and abs254) of the SUNA frames well inside the cast: the CTD's wet interval minus
                inner_margin_s at each end, and away from any surface break (spectral_qc.ctd_wet_intervals)
  water-like    a frame within the step01 tolerance of that level: peak at least (1 - frac) of it, abs254 at
                most tol above it (suncat.nitrate.spectral_qc uses the same numbers)
  stable water  the first water-like frame, within search_s of the CTD's entry, after which at least min_good of
                the next run frames (itself included) are water-like; bracketed by the frame before it. THIS IS THE
                CLOCK ESTIMATE: it needs no deck reference, and on P45_07 it reproduces the flag-based water entry
                behind the applied -16 s (stable water +15.8 s, 95% CI +15.0 to +18.2, 30 casts, Nov 2026).
  first change  walking back from stable water, the earliest frame of the unbroken run that has departed from
                the deck level just before it (peak more than depart_frac from the median of the non-water frames
                20 to 60 s before stable water), at most contact_max_s back and never across a gap in the record
                longer than max_gap_s; it is stable water itself when the frame just before is still at the deck
                level. Bracketed by the frame before it. Reported, not used for the clock: splash and a part-filled
                path make it a few seconds early (on P45_07 a median 1.4 s before stable water). A first version
                (Nov 2026) took the deck level 20 to 120 s back with a 3 % threshold: the lamp warming after the
                SUNA was switched on moved the deck level more than that, and first change ran up to 34 s early
                (P45_07 casts 02, 05, 24).
  exit          the last water-like frame, near the CTD's exit, before which at least min_good of the previous run
                frames are water-like; bracketed by the frame after it. Exit is biased late by water held in the
                optical path, so it is reported, not used for the clock.
Each event is the midpoint of its bracket, with the half-width as its timing resolution (about 0.5 s at 1 Hz,
about 6 s for P45_06's bursts 12 s apart). A cast is left out, with a reason, when the CTD record starts or ends
wet (less than min_dry_s of dry record), when the SUNA has no frames on both sides of the entry, or when no deck
or water level can be formed.

summarize() gives the cruise estimate for one definition: the median of the per-cast midpoints with a bootstrap
95% interval over casts, the IQR, and a verdict, MEASURED when at least min_casts casts give it and the interval is
at most max_ci_s wide, else UNDETERMINED with the reason.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

FRAC, TOL = 0.10, 0.05          # water-like: as step01's CTD-derived gate
DEPART_FRAC = 0.05              # first change: peak more than 5 % away from the deck level just before
CONTACT_MAX_S = 20.0            # ... at most 20 s before stable water
MAX_GAP_S = 5.0                 # ... and never across a gap in the SUNA record longer than this
RUN, MIN_GOOD = 6, 5            # stable water: 5 of the next 6 frames water-like
SEARCH_S = 90.0                 # look for the SUNA's entry within +-90 s of the CTD's
INNER_MARGIN_S = 60.0           # water level from frames this far inside the cast
MIN_DRY_S = 10.0                # the CTD record must hold this much dry time before entry / after exit
MIN_CASTS, MAX_CI_S = 5, 4.0    # verdict


def _water_like(peak, a254, wpk, wab, frac=FRAC, tol=TOL):
    ok = np.asarray(peak, float) >= (1 - frac) * wpk
    if a254 is not None and np.isfinite(wab):
        ok &= np.asarray(a254, float) <= wab + tol
    return ok


def cast_events(t, peak, a254, t_in, t_out, t_first, t_last, search_s=SEARCH_S, inner_margin_s=INNER_MARGIN_S,
                frac=FRAC, tol=TOL, depart_frac=DEPART_FRAC, run=RUN, min_good=MIN_GOOD, min_dry_s=MIN_DRY_S,
                contact_max_s=CONTACT_MAX_S, max_gap_s=MAX_GAP_S, deck_window=(60.0, 20.0), breaks=()) -> dict:
    """Entry (stable water, the clock estimate; first change, contact_*) and exit of one cast. t, peak, a254: the SUNA frames (epoch s, sorted);
    t_in, t_out: the CTD's first and last wet scan; t_first, t_last: the CTD record's ends. Offsets in s, SUNA minus
    CTD; lo/hi bracket each event. note says why an event is missing."""
    out = dict(contact_s=np.nan, contact_lo=np.nan, contact_hi=np.nan, stable_s=np.nan, stable_lo=np.nan,
               stable_hi=np.nan, exit_s=np.nan, exit_lo=np.nan, exit_hi=np.nan, water_peak=np.nan,
               deck_peak=np.nan, note="")
    t, peak = np.asarray(t, float), np.asarray(peak, float)
    a = None if a254 is None else np.asarray(a254, float)
    inner = (t >= t_in + inner_margin_s) & (t <= t_out - inner_margin_s) & np.isfinite(peak)
    for b0, b1 in breaks or ():                     # the package back at the surface inside the cast (spectral_qc)
        inner &= ~((t >= b0 - inner_margin_s) & (t <= b1 + inner_margin_s))
    if inner.sum() < 10:
        out["note"] = f"fewer than 10 SUNA frames more than {inner_margin_s:g} s inside the cast: no water level"
        return out
    wpk = float(np.median(peak[inner]))
    wab = float(np.nanmedian(a[inner])) if a is not None and np.isfinite(a[inner]).any() else np.nan
    out["water_peak"] = wpk
    wl = _water_like(peak, a, wpk, wab, frac, tol)
    notes = []
    # ---- entry
    if t_in - t_first < min_dry_s:
        notes.append("CTD record starts wet")
    else:
        idx = np.flatnonzero((t >= t_in - search_s) & (t <= t_in + search_s))
        if len(idx) < run + 1:
            notes.append("too few SUNA frames around the entry")
        else:
            st = next((i for i in idx if i >= 1 and wl[i] and not wl[i - 1] and len(wl[i:i + run]) == run
                       and wl[i:i + run].sum() >= min_good), None)
            if st is None:
                notes.append("no change from deck to stable water near the CTD's entry")
            else:
                out.update(stable_lo=t[st - 1] - t_in, stable_hi=t[st] - t_in)
                out["stable_s"] = (out["stable_lo"] + out["stable_hi"]) / 2
                deck = (t >= t[st] - deck_window[0]) & (t <= t[st] - deck_window[1]) & ~wl & np.isfinite(peak)
                if deck.sum() < 3:
                    notes.append("no deck level before the entry: first change taken as stable water")
                    c = st
                else:
                    dpk = float(np.median(peak[deck]))
                    out["deck_peak"] = dpk
                    c = st
                    while (c - 1 >= 1 and t[st] - t[c - 1] <= contact_max_s and t[c] - t[c - 1] <= max_gap_s
                           and abs(peak[c - 1] - dpk) > depart_frac * dpk):
                        c -= 1
                out.update(contact_lo=t[c - 1] - t_in, contact_hi=t[c] - t_in)
                out["contact_s"] = (out["contact_lo"] + out["contact_hi"]) / 2
    # ---- exit
    if t_last - t_out < min_dry_s:
        notes.append("CTD record ends wet")
    else:
        idx = np.flatnonzero((t >= t_out - search_s) & (t <= t_out + search_s))
        if len(idx) >= run + 1:
            ex = next((i for i in idx[::-1] if wl[i] and i + 1 < len(t) and not wl[i + 1] and i - run + 1 >= 0
                       and wl[i - run + 1:i + 1].sum() >= min_good), None)
            if ex is not None:
                out.update(exit_lo=t[ex] - t_out, exit_hi=t[ex + 1] - t_out)
                out["exit_s"] = (out["exit_lo"] + out["exit_hi"]) / 2
    out["note"] = "; ".join(notes)
    return out


def entry_table(frames: pd.DataFrame, intervals: pd.DataFrame, **kw) -> pd.DataFrame:
    """Per cast: frames (t, peak_counts, optional abs254; epoch s), intervals (cast, t_in, t_out, t_first, t_last,
    from spectral_qc.ctd_wet_intervals). One row per cast with cast_events' fields."""
    f = frames.dropna(subset=["t", "peak_counts"]).sort_values("t").reset_index(drop=True)
    a = f["abs254"].to_numpy(float) if "abs254" in f else None
    rows = []
    for _, r in intervals.iterrows():
        br = r.get("breaks", ())
        e = cast_events(f["t"].to_numpy(float), f["peak_counts"].to_numpy(float), a, r["t_in"], r["t_out"],
                        r.get("t_first", -np.inf), r.get("t_last", np.inf),
                        breaks=br if isinstance(br, (list, tuple)) else (), **kw)
        rows.append(dict(cast=r["cast"], **e))
    return pd.DataFrame(rows)


def _boot_ci(v, n_boot=5000, seed=0):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 3:
        return (np.nan, np.nan)
    b = np.median(np.random.default_rng(seed).choice(v, (n_boot, len(v))), axis=1)
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


def summarize(table: pd.DataFrame, col: str, min_casts: int = MIN_CASTS, max_ci_s: float = MAX_CI_S) -> dict:
    """Cruise estimate from the per-cast midpoints in col (contact_s, stable_s or exit_s)."""
    v = table[col].to_numpy(float); v = v[np.isfinite(v)]
    half = ((table[col.replace("_s", "_hi")] - table[col.replace("_s", "_lo")]) / 2).to_numpy(float)
    half = half[np.isfinite(half)]
    out = dict(n=len(v), median=float(np.median(v)) if len(v) else np.nan, ci=_boot_ci(v),
               iqr=(float(np.percentile(v, 25)), float(np.percentile(v, 75))) if len(v) else (np.nan, np.nan),
               half_width=float(np.median(half)) if len(half) else np.nan)
    reasons = []
    if len(v) < min_casts:
        reasons.append(f"only {len(v)} cast(s) (at least {min_casts} needed)")
    elif not (out["ci"][1] - out["ci"][0] <= max_ci_s):
        reasons.append(f"the 95% interval is {out['ci'][1] - out['ci'][0]:.1f} s wide (at most {max_ci_s:g} s)")
    out["verdict"] = "MEASURED" if not reasons else "UNDETERMINED: " + "; ".join(reasons)
    return out


def describe(s: dict, label: str) -> str:
    if not s["n"]:
        return f"  {label}: no cast gives one"
    return (f"  {label}: SUNA minus CTD median {s['median']:+.1f} s (95% CI {s['ci'][0]:+.1f} to {s['ci'][1]:+.1f} s, "
            f"IQR {s['iqr'][0]:+.1f} to {s['iqr'][1]:+.1f} s, {s['n']} casts, typical timing resolution "
            f"+-{s['half_width']:.1f} s) -> {s['verdict']}")
