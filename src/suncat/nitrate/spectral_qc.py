"""SUNA per-frame spectral QC (suna/step01_suna_spectral_qc), with the in-water reference taken from the CTD.

Two quantities in step01 rest on knowing what "in water" looks like for a file:
  reference   the median spectrum every frame's absorbance is computed against (long-wave baseline, shape)
  window      the optical out-of-water gate: peak counts, abs254 and spec_avg inside median +- k robust sigma

Until Nov 2026 both were taken from the file's own frames, which is right only when most of a file is in water.
On P45_06 the SUNA logged for hours on deck between casts, with a wet optical path, so in several files the deck
frames were the majority: the window and the reference described the deck, and real seawater frames were flagged
out of water or as a baseline anomaly (A0000017: in-water peak about 44,000 counts against a window of 28,065 to
34,532; casts 09 and 10 kept 68 and 5 usable frames). P45_07's two-day file A0000027 failed the other way: its
window (peak 10,382 to 60,000) admitted everything.

Here both come from frames taken while the CTD was in the water: inside a CTD cast's wet interval (conductivity
above 1, the first to the last wet scan) shrunk by margin_s at each end. The margin is wider than any plausible
SUNA clock error (16 s on P45_07), so no clock correction is needed and the step stays ahead of co-location. A file
with too few such frames, or with no CTD record, keeps the file's own population (the earlier rule); with too few
of those, the fixed fallback window. Which source each file used is reported.

With the window taken from seawater, two more rules follow from what out-of-water frames look like on P45_06 and
P45_07 (inspect_gate_losses, Nov 2026):
  one-sided   air, a wet optical path on deck and bubbles all give LESS light and MORE apparent absorbance than
              seawater: on both cruises the deck state read 15 % lower peak counts and abs254 about 0.08 higher.
              Seawater frames that read brighter (the first frame of each P45_06 burst, the lamp still warming at
              a soak) are not out of water, so a CTD-derived window bounds only the dark side.
  minimum     a window from stable seawater is narrow (peak +- 1 to 3 %), and real seawater varies more than
  width       that: near-bottom turbidity (P45_06 cast 04 at 90 to 93 dbar, peak 4 to 9 % low) and the settling
              of the optical path during the soak. The dark-side bound is therefore never closer to the seawater
              median than oow_dark_frac (10 %) for peak counts and spec_avg, or oow_abs254_tol (0.05) for abs254:
              about half to two-thirds of the deck-to-seawater separation.
A window from the file's own frames or the fixed fallback stays two-sided, as before.

Two refinements (Nov 2026):
  surface breaks  a CTD's conductivity cell stays full out of the water, so its wet interval runs from the first
                  to the last wet scan even when the package comes back to the surface between the soak and the
                  profile. P45_07 cast 20 sat at the surface for 194 s after its soak, and its deck-state frames
                  were among those A0000030's window and reference were built from: leaving the break out took
                  that set from 663 to 531 frames (the window's peak bound moved 0.3 %). Breaks found (Nov 2026):
                  P45_07 casts 14 and 20, P45_06 cast 07.
                  ctd_wet_intervals now lists each stretch of at least min_break_s with pressure below p_surface
                  inside the wet interval (breaks), and in_cast leaves each out, with the same margin.
  per-cast        a file can span many casts and days (P45_07's A0000027: two days), and the seawater level the
  windows         window is built around can move between casts (water type, the lamp, the optical path). With
                  oow_per_cast, the frames of a cast with at least oow_min_n_ctd in-cast frames are judged
                  against a window built from that cast's own in-cast frames, by the same one-sided rule and
                  minimum width; the other frames keep the file's window. The reference spectrum stays per file.
                  How far each cast's seawater level sits from the file's is reported, so a cast whose own window
                  rests on an unusual level is seen. In A0000027 casts 11 and 12 read 1 to 2 % above the file's
                  seawater level and casts 13 and 14 read 9 and 4 % below it, and the darker, more variable
                  frames of cast 13 had widened the file's window to 25 % below its median; each cast's own
                  window restores the 10 % minimum width there (33 more frames flagged in casts 11, 12 and 14).

Timing gives a second, separate flag, f_offcast: a frame more than margin_s away from every CTD wet interval was
not taken during a cast (deck time, a deck test), whatever its optics. It is fatal (the frame has no CTD record
to place it) but kept apart from f_outofwater, which stays purely optical so the water-entry check in step04
still reads the instrument, not the CTD.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from suncat.io.cnv import read_cnv


@dataclass
class QCConfig:
    sat_counts: float = 60000
    low_counts: float = 2000
    rmse_hi: float = 1.0e-3
    hum_hi: float = 10.0
    warmup_s: float = 60.0
    fit_lo: float = 217.0
    fit_hi: float = 240.0
    longwave: tuple = (310.0, 360.0)
    longwave_a_hi: float = 0.10
    shape: tuple = (250.0, 360.0)
    shape_dev_hi: float = 0.01
    oow_adaptive: bool = True
    oow_k: float = 6.0
    oow_min_n: int = 50            # frames needed to derive from the file's own population
    oow_min_n_ctd: int = 30        # frames needed to derive from the in-cast population
    ref_min_n: int = 5             # frames needed for a reference spectrum
    fallback_peak: tuple = (20000, 58000)       # P45_06 values, used only when nothing can be derived
    fallback_abs254: tuple = (-0.60, -0.15)
    fallback_specavg: tuple = (8000, 30000)
    oow_dark_frac: float = 0.10    # CTD-derived window: dark-side bound at least this far below the median (peak, spec_avg)
    oow_abs254_tol: float = 0.05   # CTD-derived window: abs254 bound at least this far above the median
    oow_per_cast: bool = True      # a cast with enough in-cast frames gets its own window (needs cast_of)
    peak_hard: tuple = (2000, 60000)
    abs254_hard: tuple = (-5.0, 5.0)
    specavg_hard: tuple = (0.0, 60000)


L0_FLAGS = ["f_warmup", "f_saturated", "f_lowsignal", "f_highrmse", "f_humid"]
OTHER_FATAL = ["f_warmup", "f_saturated", "f_lowsignal", "f_highrmse", "f_satfit", "f_baseline", "f_shape"]
FATAL = ["f_warmup", "f_saturated", "f_lowsignal", "f_highrmse", "f_outofwater", "f_satfit", "f_baseline", "f_shape",
         "f_offcast"]


# ------------------------------------------------------------------ the CTD's time in the water
def _start_time(path) -> datetime | None:
    for ln in Path(path).read_text(encoding="latin-1", errors="replace").splitlines():
        if ln.startswith("# start_time"):
            s = " ".join(ln.split("=", 1)[1].split("[")[0].split())
            try:
                return datetime.strptime(s, "%b %d %Y %H:%M:%S")
            except ValueError:
                return None
        if ln.strip().upper().startswith("*END*"):
            break
    return None


def surface_breaks(ts, p, i0: int, i1: int, p_surface: float = 0.5, min_break_s: float = 30.0) -> list:
    """(start, end) elapsed times of each stretch of at least min_break_s, between scan indices i0 and i1, with
    pressure below p_surface dbar: the package at or above the surface inside a wet interval."""
    ts, p = np.asarray(ts, float)[i0:i1 + 1], np.asarray(p, float)[i0:i1 + 1]
    low = np.isfinite(p) & (p < p_surface)
    edges = np.flatnonzero(np.diff(np.r_[0, low.astype(int), 0]))
    return [(float(ts[a]), float(ts[b - 1])) for a, b in zip(edges[::2], edges[1::2])
            if ts[b - 1] - ts[a] >= min_break_s]


def ctd_wet_intervals(paths, cond_in: float = 1.0, p_in: float = 1.0, p_surface: float = 0.5,
                      min_break_s: float = 30.0) -> tuple[pd.DataFrame, list]:
    """(cast, t_in, t_out, by, t_first, t_last, breaks) per CTD .cnv, epoch seconds UTC: the first and last scan
    with conductivity above cond_in (1 separates air from seawater in S/m or mS/cm), or, in a file with no
    conductivity column, with pressure above p_in dbar (by says which), and the record's first and last scan.
    breaks: [(start, end)] of each surface break inside the wet interval (surface_breaks; empty without a pressure
    column). Also returns [(file name, reason)] for files left out."""
    rows, skipped = [], []
    for p in sorted(Path(x) for x in paths):
        start = _start_time(p)
        d = read_cnv(p, mask_badflag=True)
        if start is None or d is None or d.empty:
            skipped.append((p.name, "no start_time or no data")); continue
        tcol = next((c for c in d.columns if c.lower().startswith("times")), None)
        ccol = next((c for c in d.columns if c.lower().startswith(("c0", "cond"))), None)
        pcol = next((c for c in d.columns if c.lower().startswith(("prdm", "prm", "pr"))), None)
        if tcol is None or (ccol is None and pcol is None):
            skipped.append((p.name, f"no elapsed time, conductivity or pressure column ({list(d.columns)[:8]})"))
            continue
        col, thr, by = (ccol, cond_in, "conductivity") if ccol else (pcol, p_in, "pressure")
        c = d[col].to_numpy(float)
        wet = np.flatnonzero(np.isfinite(c) & (c > thr))
        if not len(wet):
            skipped.append((p.name, "never wet")); continue
        t0 = (pd.Timestamp(start) - pd.Timestamp("1970-01-01")).total_seconds()
        ts = d[tcol].to_numpy(float)
        br = ([(t0 + a, t0 + b) for a, b in surface_breaks(ts, d[pcol], wet[0], wet[-1], p_surface, min_break_s)]
              if pcol else [])
        rows.append(dict(cast=p.stem, t_in=t0 + ts[wet[0]], t_out=t0 + ts[wet[-1]], by=by,
                         t_first=t0 + float(np.nanmin(ts)), t_last=t0 + float(np.nanmax(ts)), breaks=br))
    return pd.DataFrame(rows, columns=["cast", "t_in", "t_out", "by", "t_first", "t_last", "breaks"]), skipped


def in_cast(t_epoch, intervals: pd.DataFrame, margin_s: float = 60.0) -> np.ndarray:
    """True where t lies inside a cast's wet interval shrunk by margin_s at each end, and more than margin_s from
    each of its surface breaks (the breaks column, when present)."""
    t = np.asarray(t_epoch, float)
    m = np.zeros(t.shape, bool)
    for a, b in zip(intervals["t_in"], intervals["t_out"]):
        if b - a > 2 * margin_s:
            m |= (t >= a + margin_s) & (t <= b - margin_s)
    if "breaks" in intervals:
        for br in intervals["breaks"]:
            for a, b in (br if isinstance(br, (list, tuple)) else []):
                m &= ~((t >= a - margin_s) & (t <= b + margin_s))
    return m


def off_cast(t_epoch, intervals: pd.DataFrame, margin_s: float = 60.0) -> np.ndarray:
    """True where t lies more than margin_s outside every cast's wet interval (False everywhere with no casts)."""
    t = np.asarray(t_epoch, float)
    if intervals is None or not len(intervals):
        return np.zeros(t.shape, bool)
    near = np.zeros(t.shape, bool)
    for a, b in zip(intervals["t_in"], intervals["t_out"]):
        near |= (t >= a - margin_s) & (t <= b + margin_s)
    return ~near


def which_cast(t_epoch, intervals: pd.DataFrame) -> np.ndarray:
    """The cast whose wet interval (no margin) holds each time; '' between casts."""
    t = np.asarray(t_epoch, float)
    out = np.full(t.shape, "", dtype=object)
    for c, a, b in zip(intervals["cast"], intervals["t_in"], intervals["t_out"]):
        out[(t >= a) & (t <= b)] = c
    return out


# ------------------------------------------------------------------ the gate
def robust_window(x, k: float, hard, min_n: int):
    """median +- k * MAD-sigma, clipped to hard physical bounds; None if fewer than min_n finite values."""
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if x.size < min_n:
        return None
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    if not np.isfinite(mad) or mad <= 0:
        s = x.std(); mad = s if s > 0 else max(abs(med) * 0.05, 1.0)
    return (float(max(med - k * mad, hard[0])), float(min(med + k * mad, hard[1])))


def _dark_side(x, k: float, hard, min_n: int, frac: float | None = None, tol: float | None = None, high_is_dark=False):
    """One-sided window from seawater frames: (lower bound, hard upper) for peak and spec_avg, (hard lower, upper
    bound) for abs254 (high_is_dark). The bound is median -+ k * MAD-sigma, but never closer to the median than frac
    (a fraction of it) or tol (an absolute amount). None if fewer than min_n finite values."""
    w = robust_window(x, k, hard, min_n)
    if w is None:
        return None
    xx = np.asarray(x, float); med = float(np.median(xx[np.isfinite(xx)]))
    if high_is_dark:
        return (float(hard[0]), float(min(max(w[1], med + tol), hard[1])))
    return (float(max(min(w[0], med - abs(med) * frac), hard[0])), float(hard[1]))


def derive_window(df: pd.DataFrame, good, cfg: QCConfig, incast=None) -> tuple[dict, str, int]:
    """The out-of-water window, where it came from ("ctd", "file" or "fixed") and how many frames it rests on.
    "ctd": good frames inside a CTD cast (at least oow_min_n_ctd), bounded on the dark side only and never closer
    to the median than oow_dark_frac / oow_abs254_tol; "file": all good frames (at least oow_min_n), two-sided;
    "fixed": the fallback values."""
    fallback = {"peak": cfg.fallback_peak, "abs254": cfg.fallback_abs254, "spec_avg": cfg.fallback_specavg}
    good = np.asarray(good, bool)
    if not cfg.oow_adaptive:
        return fallback, "fixed", 0
    if incast is not None and (good & incast).sum() >= cfg.oow_min_n_ctd:
        use, src, mn = good & incast, "ctd", cfg.oow_min_n_ctd
    elif good.sum() >= cfg.oow_min_n:
        use, src, mn = good, "file", cfg.oow_min_n
    else:
        return fallback, "fixed", 0
    g = df.loc[use]
    if src == "ctd":
        win = {"peak": _dark_side(g["peak_counts"], cfg.oow_k, cfg.peak_hard, mn, frac=cfg.oow_dark_frac)
                       or fallback["peak"],
               "abs254": _dark_side(g["abs254"], cfg.oow_k, cfg.abs254_hard, mn, tol=cfg.oow_abs254_tol,
                                    high_is_dark=True) or fallback["abs254"],
               "spec_avg": _dark_side(g["spec_avg"], cfg.oow_k, cfg.specavg_hard, mn, frac=cfg.oow_dark_frac)
                           or fallback["spec_avg"]}
        return win, src, int(use.sum())
    win = {"peak": robust_window(g["peak_counts"], cfg.oow_k, cfg.peak_hard, mn) or fallback["peak"],
           "abs254": robust_window(g["abs254"], cfg.oow_k, cfg.abs254_hard, mn) or fallback["abs254"],
           "spec_avg": robust_window(g["spec_avg"], cfg.oow_k, cfg.specavg_hard, mn) or fallback["spec_avg"]}
    return win, src, int(use.sum())


def qc_frames(df: pd.DataFrame, spec: np.ndarray, wl: np.ndarray, cfg: QCConfig | None = None, incast=None,
              offcast=None, cast_of=None):
    """L0 integrity + L1 spectral QC for one file. incast (bool per frame, from in_cast) makes the reference
    spectrum and the out-of-water window come from frames taken while the CTD was in the water; None keeps the
    file's own population (the rule before Nov 2026). offcast (bool per frame, from off_cast) sets f_offcast;
    None leaves it False. cast_of (the cast each frame belongs to, '' between casts, from which_cast) with
    cfg.oow_per_cast gives each cast with at least oow_min_n_ctd in-cast frames its own window. Returns (df with
    flags, the file's window, info) where info holds window_source, window_n, ref_source, ref_n, ref (the reference
    spectrum) and cast_windows (one dict per cast: cast, source "cast" or "file", n, the window, and the cast's
    median peak, abs254 and spec_avg against the file's in-cast medians)."""
    cfg = cfg or QCConfig()
    if df.empty:
        return df, None, dict(window_source="none", window_n=0, ref_source="none", ref_n=0)
    df = df.copy()
    t0 = df["time_utc"].iloc[0]
    secs = (df["time_utc"] - t0).dt.total_seconds().to_numpy()
    peak = spec.max(axis=1)
    df["peak_counts"] = peak
    df["f_warmup"] = secs < cfg.warmup_s
    df["f_saturated"] = peak > cfg.sat_counts
    df["f_lowsignal"] = peak < cfg.low_counts
    df["f_highrmse"] = df["fit_rmse"] > cfg.rmse_hi
    df["f_humid"] = df["humidity"] > cfg.hum_hi
    l0ok = ~df[L0_FLAGS].any(axis=1).to_numpy()
    inc = None if incast is None else np.asarray(incast, bool)
    # reference spectrum: in-cast frames when there are enough, else the file's L0-good frames, else all frames
    if inc is not None and (l0ok & inc).sum() >= max(cfg.ref_min_n, cfg.oow_min_n_ctd):
        ref, ref_src, ref_n = np.median(spec[l0ok & inc], axis=0), "ctd", int((l0ok & inc).sum())
    elif l0ok.sum() >= cfg.ref_min_n:
        ref, ref_src, ref_n = np.median(spec[l0ok], axis=0), "file", int(l0ok.sum())
    else:
        ref, ref_src, ref_n = np.median(spec, axis=0), "file (all frames)", len(spec)
    dk = float(np.median(df["dark"]))
    A = -np.log10(np.clip(spec - dk, 1, None) / np.clip(ref - dk, 1, None))
    lw = (wl >= cfg.longwave[0]) & (wl <= cfg.longwave[1])
    sb = (wl >= cfg.shape[0]) & (wl <= cfg.shape[1])
    fitb = (wl >= cfg.fit_lo) & (wl <= cfg.fit_hi)
    df["longwave_A"] = np.mean(np.abs(A[:, lw]), axis=1)
    rb = ref[sb] - ref[sb].mean()
    sbn = spec[:, sb] - spec[:, sb].mean(axis=1, keepdims=True)
    corr = (sbn @ rb) / (np.sqrt((sbn ** 2).sum(1)) * np.sqrt((rb ** 2).sum()) + 1e-9)
    df["shape_dev"] = 1.0 - corr
    df["f_satfit"] = spec[:, fitb].max(axis=1) >= cfg.sat_counts
    df["f_baseline"] = df["longwave_A"] > cfg.longwave_a_hi
    df["f_shape"] = df["shape_dev"] > cfg.shape_dev_hi
    good0 = ~df[OTHER_FATAL].any(axis=1).to_numpy()
    win, src, n = derive_window(df, good0, cfg, inc)

    def _inside(w, m):
        return (df.loc[m, "peak_counts"].between(*w["peak"]) & df.loc[m, "abs254"].between(*w["abs254"])
                & df.loc[m, "spec_avg"].between(*w["spec_avg"]))

    allm = np.ones(len(df), bool)
    ok = _inside(win, allm).to_numpy().copy()
    cast_windows = []
    if cast_of is not None and inc is not None:
        co = np.asarray(cast_of, dtype=object)
        base = good0 & inc
        fmed = {q: float(np.median(df.loc[base, q])) if base.any() else np.nan
                for q in ("peak_counts", "abs254", "spec_avg")}
        for c in sorted(x for x in set(co) if x):
            m = co == c
            use = base & m
            row = dict(cast=c, source="file", n=int(use.sum()), window=win)
            for q in ("peak_counts", "abs254", "spec_avg"):
                row[f"{q}_median"] = float(np.median(df.loc[use, q])) if use.any() else np.nan
                row[f"{q}_file_median"] = fmed[q]
            if cfg.oow_per_cast and src == "ctd" and use.sum() >= cfg.oow_min_n_ctd:
                wc, _, _ = derive_window(df, use, cfg, use)
                ok[m] = _inside(wc, m).to_numpy()
                row.update(source="cast", window=wc)
            cast_windows.append(row)
    df["f_outofwater"] = ~ok
    df["f_offcast"] = np.zeros(len(df), bool) if offcast is None else np.asarray(offcast, bool)
    df["L1_clean"] = ~df[FATAL].any(axis=1)
    return df, win, dict(window_source=src, window_n=n, ref_source=ref_src, ref_n=ref_n, ref=ref,
                         cast_windows=cast_windows)
