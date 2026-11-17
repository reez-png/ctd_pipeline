"""Which sensor value best represents the water a bottle closed in? Compares bottle-matching estimators.

    python tools/compare_match_estimators.py P45_07                 chlorophyll fluorescence (default)
    python tools/compare_match_estimators.py P45_07 --var nitrate   SUNA nitrate
    python tools/compare_match_estimators.py P45_07 --var nitrate --flush-s 15
    python tools/compare_match_estimators.py P45_07 --var nitrate --sweep   also the offset sweep (lag vs geometry)
    python tools/compare_match_estimators.py P45_07 --var nitrate --sweep --range 30   a wider sweep

Run from the repo root, after step07 (chl) or step03 (nitrate). A diagnostic: it reads the products and writes a
per-bottle table to the cruise's _audit; it changes nothing in the pipeline.

Bottles. chl: step07's calibration table (_audit/productivity/<CRUISE>__chl_calibration.csv), without NPQ-flagged
bottles and dead-fluorometer casts. nitrate: step03's match-up table (_audit/step03_bottle_matchups.csv), else
metadata/bottle_nitrate/<CRUISE>_bottle_nitrate.csv.

Samples. chl: the cast's native-rate (4 Hz) Derived scans (L1/CTD/Derived/<cast>_der.cnv), sensor = the first of
flECO-AFL, flECO-AFL1, wetCDOM, flC; blanks (<= -0.1) dropped; down-cast = Loop Edit-kept scans up to the maximum
kept pressure. nitrate: step03's frames (L2/nitrate/<cast>_nitrate_QC.csv) that passed SUNA spectral QC
(L1_clean), value = no3_stage2 (nitrate before the bottle calibration; else no3_stage1, else no3_qc), leg as step03
labelled it, down-cast frames limited to those whose concurrent CTD scan Loop Edit kept when step03 recorded it.

Estimators, per bottle:
  (down and up use samples moving faster than MOVING_W = 0.2 m/s; the stationary ones are the stop's)
  down       median of 1 m depth-cell means within +-WIN of the bottle depth, down-cast samples (step07's rule; the
             cells are laid out centred on the depth matched, which for a whole-metre depth is step07's grid)
  up         the same on the up-cast samples
  down@fire  as down, centred on the bottle's FIRING depth (bottle tables carry nominal depths), so down -> down@fire
             isolates the depth effect and down@fire -> stop the leg-and-timing effect
  up@fire    as up, centred on the firing depth. When bottles are fired on the fly (the rosette still hauling, no
             stop), this is the sensor's view of the bottle's water: the same water, passing the sensor at about
             the firing time
  stop       median of the STATIONARY samples (|w| < MOVING_W) while the rosette sits at the bottle: from
             STOP_BEFORE_S before to STOP_AFTER_S after the firing time, within STOP_DP dbar of the firing pressure.
             Samples passing that pressure on the move are not a stop and do not count: a cast fired on the fly has
             no stop estimator, and the speed at each firing (fire_w) is reported. For nitrate the first --flush-s
             seconds of that stationary period are discarded (the flow cell is still flushing with the water the
             package arrived from); the early-minus-late difference is reported as the flushing signature
Firing time, pressure and depth come from metadata/bottle_nitrate/<CRUISE>_bottle_ctd.csv (build_bottle_nitrate_
table.py, from the Sea-Bird .btl files; CTD clock, UTC); a bottle is paired with the firing on its cast nearest its
depth, within PAIR_DZ m. Without that file only down and up are compared.

Statistics, on the SAME bottles (those with every estimator):
  * OLS sensor = offset + gain * bottle (step07's fit), gain with 95% t interval, R2;
  * leave-one-out cross-validation (LOOCV): each bottle predicted, (sensor - offset) / gain, from the fit made
    without it; RMSE in the bottle's units. It scores how well an estimator calibrates water it has not seen,
    so the estimator that fits its own bottles best is not automatically favoured (Hastie et al., The Elements
    of Statistical Learning, ch. 7);
  * paired LOOCV squared errors against `down`, bottle by bottle: median difference, exact sign test and, when
    scipy is installed, Wilcoxon signed-rank. With 15 to 20 bottles read p-values as weak evidence; the same
    direction on two cruises carries more weight than either p-value.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import project_paths as pp  # noqa: E402

WIN, BLANK = 2.0, -0.1
STOP_BEFORE_S, STOP_AFTER_S, STOP_DP = 30.0, 3.0, 1.0
PAIR_DZ = 3.0
SENSOR_COLS = ["flECO-AFL", "flECO-AFL1", "wetCDOM", "flC"]
ESTIMATORS = ("down", "up", "down@fire", "up@fire", "stop")
FIRE_W_S = 3.0       # s either side of a firing over which the rosette's speed at the firing (fire_w) is taken


def t975(df: int) -> float:
    table = [12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145,
             2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048,
             2.045, 2.042]
    return table[df - 1] if 1 <= df <= len(table) else 1.96


def _truthy(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])


def _first_col(names, *patterns):
    for pat in patterns:
        for c in names:
            if re.fullmatch(pat, str(c).strip(), re.I):
                return c
    return None


def parse_utc(col: pd.Series, what: str = "time") -> pd.Series:
    """Timestamps read as UTC with format="ISO8601", which accepts whole and fractional seconds in one column.
    Without it pandas guesses the format from the FIRST value and, with errors="coerce", silently turns every value
    written the other way into NaT (step03 writes frame times both ways). Warns if any value still fails."""
    t = pd.to_datetime(col, utc=True, errors="coerce", format="ISO8601")
    bad = int(t.isna().sum() - col.isna().sum())
    if bad > 0:
        print(f"  WARNING: {bad} of {len(col)} {what} values could not be read as times")
    return t


# ------------------------------------------------------------------ samples: t, p, z, v, down (bool), up (bool)
def read_ctd_scans(path: Path) -> pd.DataFrame:
    """Derived 4 Hz scans of one cast as samples (chlorophyll fluorescence)."""
    lines = path.read_text(encoding="latin-1", errors="replace").splitlines()
    names, start, i0 = [], None, None
    for i, ln in enumerate(lines):
        if ln.startswith("# name"):
            names.append(ln.split("=", 1)[1].split(":")[0].strip())
        m = re.match(r"#\s*start_time\s*=\s*([A-Za-z]{3}\s+\d{1,2}\s+\d{4}\s+\d{2}:\d{2}:\d{2})", ln)
        if m:
            start = pd.to_datetime(m.group(1), format="%b %d %Y %H:%M:%S", utc=True)
        if ln.strip().upper().startswith("*END*"):
            i0 = i + 1
            break
    rows = []
    for ln in lines[i0:]:
        parts = ln.split()
        if len(parts) >= len(names):
            try:
                rows.append([float(x) for x in parts[:len(names)]])
            except ValueError:
                pass
    d = np.array(rows, float)
    d[(d <= -8e-29) & (d >= -1.2e-28)] = np.nan                       # Loop Edit's bad flag -> NaN
    df = pd.DataFrame(d, columns=names)
    pcol = _first_col(names, r"pr[dDsS]?M", r"pressure.*")
    zcol = _first_col(names, r"depSM", r"depFM", r"depth.*") or pcol
    scol = next((c for c in SENSOR_COLS if c in names), None)
    kept = df["flag"].eq(0) if "flag" in df else pd.Series(True, index=df.index)
    imax = int(np.nanargmax(df[pcol].where(kept).to_numpy()))
    v = pd.to_numeric(df[scol], errors="coerce") if scol else pd.Series(np.nan, index=df.index)
    t = start + pd.to_timedelta(df["timeS"], unit="s") if (start is not None and "timeS" in df) else pd.NaT
    idx = pd.Series(np.arange(len(df)), index=df.index)
    return add_motion(pd.DataFrame({"t": t, "p": df[pcol], "z": df[zcol], "v": v.where(~(v <= BLANK)),
                                    "down": kept & (idx <= imax), "up": idx > imax}))


def read_suna_frames(path: Path) -> pd.DataFrame:
    """step03 frames of one cast as samples (nitrate before the bottle calibration), spectrally clean only."""
    return read_suna_frames_df(pd.read_csv(path))


def read_suna_frames_df(f: pd.DataFrame, vcol: str | None = None) -> pd.DataFrame:
    """Frames as samples. vcol: the nitrate column (default the latest stage present: no3_stage2, no3_stage1,
    no3_qc)."""
    if vcol is not None and vcol not in f.columns:
        raise KeyError(f"no column {vcol!r} in the frame file (columns: {', '.join(f.columns)})")
    vcol = vcol or next((c for c in ("no3_stage2", "no3_stage1", "no3_qc") if c in f.columns), None)
    clean = _truthy(f["L1_clean"]) if "L1_clean" in f else pd.Series(True, index=f.index)
    leg = f["leg"].astype(str).str.lower() if "leg" in f else pd.Series("down", index=f.index)
    kept = _truthy(f["ctd_scan_kept"]) if "ctd_scan_kept" in f else pd.Series(True, index=f.index)
    # pressure for the stop test: the frames' own CTD pressure (1 s record). NOT ctd_scan_pressure_db: Loop Edit
    # rejects scans slower than its minimum velocity, so at a bottle stop every scan is flagged and that column is
    # empty there
    pcol = "pressure_db" if "pressure_db" in f else "ctd_scan_pressure_db"
    v = pd.to_numeric(f[vcol], errors="coerce") if vcol else pd.Series(np.nan, index=f.index)
    return pd.DataFrame({"t": parse_utc(f["utc_time"], "frame utc_time"),
                         "p": pd.to_numeric(f[pcol], errors="coerce"),
                         "z": pd.to_numeric(f["depth_m"], errors="coerce"),
                         "v": v.where(clean), "down": leg.eq("down") & kept, "up": leg.eq("up")}).pipe(add_motion)


MOVING_W = 0.2       # m/s: slower samples are "stationary" (a bottle stop), the 'stop' estimator's, not a leg's


def add_motion(s: pd.DataFrame) -> pd.DataFrame:
    """Vertical speed w (m/s, positive down) from depth over a centred difference (+-1 s, at least one sample each
    side, at most 6 s), and `moving` = |w| >= MOVING_W.
    down and up use moving samples only; the stationary ones (bottle stops) belong to the 'stop' estimator."""
    s = s.copy()
    if s["t"].isna().all():
        s["w"], s["moving"] = np.nan, True
        return s
    # time order, rows with a time and a depth only: a frame file need not be sorted and may lack some times
    g = s[s["t"].notna() & s["z"].notna()].sort_values("t", kind="stable")
    tt = (g["t"] - g["t"].iloc[0]).dt.total_seconds().to_numpy(float)
    z = g["z"].to_numpy(float)
    n = len(tt)
    w = np.full(n, np.nan)
    if n >= 2:
        # +-1 s, but always reaching at least the previous and the next sample: SUNA frames are often more than
        # 1 s apart, and a frame with no neighbour inside +-1 s must still get a speed. Spans over 6 s (a gap, or
        # a stop's edge) give no speed.
        i = np.arange(n)
        j0 = np.clip(np.minimum(np.searchsorted(tt, tt - 1.0), i - 1), 0, n - 1)
        j1 = np.clip(np.maximum(np.searchsorted(tt, tt + 1.0, side="right") - 1, i + 1), 0, n - 1)
        span = tt[j1] - tt[j0]
        ok = (span > 0) & (span <= 6.0)
        w[ok] = (z[j1] - z[j0])[ok] / span[ok]
    s["w"] = np.nan
    s.loc[g.index, "w"] = w
    s["moving"] = s["w"].abs() >= MOVING_W                         # no speed (no time) -> not moving -> no leg
    return s


def net_speed(samples, leg: str) -> float:
    """Net profiling speed of a leg (m/s): depth travelled over time while moving, summed over casts. Ship heave
    averages out of this; it is the speed that converts a response time into a lag distance."""
    dz, dt = 0.0, 0.0
    for s in samples:
        if s is None or s["t"].isna().all():
            continue
        g = s[s[leg] & s["moving"]].dropna(subset=["t", "z"]).sort_values("t")
        if len(g) < 2:
            continue
        step_t = np.diff((g["t"] - g["t"].iloc[0]).dt.total_seconds().to_numpy())
        step_z = np.diff(g["z"].to_numpy(float))
        keep = (step_t > 0) & (step_t <= 2.5)                    # consecutive samples only, not across a stop
        dz += float(np.sum(step_z[keep])) * (1.0 if leg == "down" else -1.0)
        dt += float(np.sum(step_t[keep]))
    return dz / dt if dt > 0 else float("nan")


# ------------------------------------------------------------------ estimators
def cell_median(s: pd.DataFrame, z0: float) -> tuple[float, int]:
    """Median of 1 m depth-cell means within +-WIN of z0, the cells laid out CENTRED ON z0 (z0 + k m, k integer).
    For a whole-metre z0 this is step07's rule exactly (cells on whole metres). Centring on z0 matters for a
    non-integer depth such as a firing depth: on a fixed whole-metre grid a shift under 0.5 m selects the same
    cells and cannot change the value."""
    ok = s["z"].notna() & s["v"].notna()
    c = np.floor(s.loc[ok, "z"].to_numpy(float) - z0 + 0.5)                 # cell index relative to z0
    m = np.abs(c) <= WIN
    if not m.any():
        return np.nan, 0
    means = pd.Series(s.loc[ok, "v"].to_numpy(float)[m]).groupby(c[m]).mean()
    return float(np.median(means)), int(means.size)


def stop_value(s: pd.DataFrame, t_fire, p_fire: float, flush_s: float) -> dict:
    """Median of the stationary samples at the bottle stop, after discarding the first flush_s seconds of it; also
    fire_w, the rosette's median |speed| within FIRE_W_S s of the firing (m/s): above MOVING_W, the bottle was fired
    on the fly and there is no stop to sample."""
    dt = (s["t"] - t_fire).dt.total_seconds()
    still = ~s["moving"].astype(bool) if "moving" in s else pd.Series(True, index=s.index)
    w = s[(dt >= -STOP_BEFORE_S) & (dt <= STOP_AFTER_S) & ((s["p"] - p_fire).abs() <= STOP_DP) & s["v"].notna()
          & still]
    out = dict(stop=np.nan, stop_n=0, stop_early_minus_late=np.nan, stop_dt_at_p=np.nan, stop_dp_at_t=np.nan,
               fire_w=np.nan)
    if "w" in s:
        near = s.loc[dt.abs() <= FIRE_W_S, "w"].abs().dropna()
        if len(near):
            out["fire_w"] = float(near.median())
    # diagnostics, kept for every bottle: how far in time from the firing the stationary samples AT the firing
    # pressure are (all samples there if none is stationary),
    # and how far in pressure from the firing pressure the samples nearest the firing TIME are
    sv = s[s["v"].notna() & s["t"].notna() & s["p"].notna()]
    if len(sv):
        at_p = sv[(sv["p"] - p_fire).abs() <= STOP_DP]
        if "moving" in sv and (~at_p["moving"].astype(bool)).any():
            at_p = at_p[~at_p["moving"].astype(bool)]          # the stop itself, not a frame passing that pressure
        if len(at_p):
            d = (at_p["t"] - t_fire).dt.total_seconds()
            out["stop_dt_at_p"] = float(d.iloc[int(np.argmin(np.abs(d.to_numpy())))])
        k = int(np.argmin(np.abs((sv["t"] - t_fire).dt.total_seconds().to_numpy())))
        out["stop_dp_at_t"] = float(sv["p"].iloc[k] - p_fire)
    if w.empty:
        return out
    t0 = w["t"].min()
    late = w[(w["t"] - t0).dt.total_seconds() >= flush_s]
    early = w[(w["t"] - t0).dt.total_seconds() < flush_s]
    if len(late):
        out.update(stop=float(late["v"].median()), stop_n=int(len(late)))
    if len(late) and len(early):
        out["stop_early_minus_late"] = float(early["v"].median() - late["v"].median())
    return out


# ------------------------------------------------------------------ statistics
def ols(b, s) -> dict:
    b, s = np.asarray(b, float), np.asarray(s, float)
    (g, o), cov = np.polyfit(b, s, 1, cov=True)
    r = s - (o + g * b)
    se = float(np.sqrt(cov[0, 0]))
    tc = t975(len(b) - 2)
    return dict(gain=g, lo=g - tc * se, hi=g + tc * se, offset=o, R2=1 - np.sum(r ** 2) / np.sum((s - s.mean()) ** 2))


def loo_errors(b, s) -> np.ndarray:
    """Leave-one-out prediction errors, in bottle units: (s_i - offset_-i) / gain_-i - b_i."""
    b, s = np.asarray(b, float), np.asarray(s, float)
    err = np.full(len(b), np.nan)
    for i in range(len(b)):
        k = np.arange(len(b)) != i
        g, o = np.polyfit(b[k], s[k], 1)
        if g != 0:
            err[i] = (s[i] - o) / g - b[i]
    return err


def loo_forward_errors(b, s) -> np.ndarray:
    """Leave-one-out errors of the forward fit, in SENSOR units: s_i - (offset_-i + gain_-i * b_i). Used by the
    sweep: it stays finite when a shifted window drives the gain towards zero, where the inverse form explodes."""
    b, s = np.asarray(b, float), np.asarray(s, float)
    err = np.full(len(b), np.nan)
    for i in range(len(b)):
        k = np.arange(len(b)) != i
        g, o = np.polyfit(b[k], s[k], 1)
        err[i] = s[i] - (o + g * b[i])
    return err


def sign_test(d) -> float:
    """Exact two-sided sign test on paired differences (zeros dropped); NaN if none differ."""
    d = np.asarray(d, float)
    d = d[np.isfinite(d) & (d != 0)]
    n, k = len(d), int((d > 0).sum())
    if n == 0:
        return float("nan")
    return min(1.0, 2 * sum(math.comb(n, j) for j in range(0, min(k, n - k) + 1)) / 2 ** n)


# ------------------------------------------------------------------ bottles and firings
def load_bottles(cruise: str, var: str) -> pd.DataFrame:
    cdir = pp.cruise_dir(cruise)
    if var == "chl":
        cal = pd.read_csv(cdir / "_audit" / "productivity" / f"{cruise}__chl_calibration.csv", encoding="utf-8-sig")
        note = cal["note"].astype(str) if "note" in cal else pd.Series("", index=cal.index)
        npq = _truthy(cal["npq"]) if "npq" in cal else pd.Series(False, index=cal.index)
        cal = cal[cal["bottle"].notna() & ~npq & ~note.str.contains("dead", case=False)]
        return pd.DataFrame({"cast": cal["cast"].astype(str), "depth_m": cal["depth_m"].astype(float),
                             "bottle": cal["bottle"].astype(float)})
    mu = cdir / "_audit" / "step03_bottle_matchups.csv"
    if mu.exists():
        m = pd.read_csv(mu, encoding="utf-8-sig")
        return pd.DataFrame({"cast": m["cast_id"].astype(str), "depth_m": m["depth_m"].astype(float),
                             "bottle": m["bottle_no3"].astype(float)}).dropna()
    b = pd.read_csv(cdir / "metadata" / "bottle_nitrate" / f"{cruise}_bottle_nitrate.csv", encoding="utf-8-sig")
    return pd.DataFrame({"cast": b["cast_id"].astype(str), "depth_m": b["depth_m"].astype(float),
                         "bottle": b["nitrate_uM"].astype(float)}).dropna()


def load_firings(cruise: str) -> pd.DataFrame | None:
    f = pp.cruise_dir(cruise) / "metadata" / "bottle_nitrate" / f"{cruise}_bottle_ctd.csv"
    if not f.exists():
        return None
    b = pd.read_csv(f, encoding="utf-8-sig")
    b.columns = [c.strip().lower() for c in b.columns]
    pc = next((c for c in ("pressure_dbar", "pressure_db", "prdm") if c in b.columns), None)
    if pc is None or "time_utc" not in b or "depth_m" not in b:
        return None
    return pd.DataFrame({"cast": b["cast_id"].astype(str), "depth_m": pd.to_numeric(b["depth_m"], errors="coerce"),
                         "p": pd.to_numeric(b[pc], errors="coerce"),
                         "t": parse_utc(b["time_utc"], "bottle time_utc")})


# ------------------------------------------------------------------ the comparison
def compare(cruise: str, var: str = "chl", flush_s: float | None = None, verbose: bool = True) -> dict:
    flush = (10.0 if var == "nitrate" else 0.0) if flush_s is None else float(flush_s)
    bottles = load_bottles(cruise, var)
    firings = load_firings(cruise)
    say = print if verbose else (lambda *a, **k: None)
    if firings is None:
        say(f"no {cruise}_bottle_ctd.csv with time_utc and pressure: only 'down' and 'up' "
            "(make it with build_bottle_nitrate_table.py from the .btl files)")
    src = (pp.product_dir(cruise, "ctd_derived"), "{}_der.cnv", read_ctd_scans) if var == "chl" else \
          (pp.product_dir(cruise, "nitrate_qc"), "{}_nitrate_QC.csv", read_suna_frames)
    cache, rows = {}, []
    for _, r in bottles.iterrows():
        cast, z0 = r["cast"], float(r["depth_m"])
        if cast not in cache:
            f = src[0] / src[1].format(cast)
            cache[cast] = src[2](f) if f.exists() else None
        s = cache[cast]
        rec = dict(cast=cast, depth_m=z0, bottle=float(r["bottle"]))
        if s is None:
            rows.append(rec)
            continue
        dn, up = s[s["down"] & s["moving"]], s[s["up"] & s["moving"]]
        rec["down"], rec["down_cells"] = cell_median(dn, z0)
        rec["up"], rec["up_cells"] = cell_median(up, z0)
        rec.update({"down@fire": np.nan, "up@fire": np.nan, "fire_dz": np.nan, "stop": np.nan, "stop_n": 0,
                    "stop_early_minus_late": np.nan, "stop_dt_at_p": np.nan, "stop_dp_at_t": np.nan,
                    "fire_w": np.nan})
        if firings is not None:
            fb = firings[firings["cast"] == cast]
            if len(fb):
                j = (fb["depth_m"] - z0).abs().idxmin()
                dz = float(fb.at[j, "depth_m"] - z0)
                rec["fire_dz"] = dz
                if abs(dz) <= PAIR_DZ:
                    rec["down@fire"], _ = cell_median(dn, z0 + dz)
                    rec["up@fire"], _ = cell_median(up, z0 + dz)
                    if pd.notna(fb.at[j, "t"]) and s["t"].notna().any():
                        rec.update(stop_value(s, fb.at[j, "t"], float(fb.at[j, "p"]), flush))
        rows.append(rec)
    t = pd.DataFrame(rows)
    out_dir = pp.cruise_dir(cruise) / "_audit" / ("productivity" if var == "chl" else "")
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{cruise}__match_estimators_{var}.csv"
    t.to_csv(out, index=False, encoding="utf-8-sig")
    est = [e for e in ESTIMATORS if e in t and t[e].notna().sum() >= 4]
    both = t.dropna(subset=est)
    unit = "mg/m3" if var == "chl" else "umol/L"
    say(f"{cruise} {var}: {len(t)} bottles; {len(both)} have all of {est}"
        + (f"; stop discards the first {flush:g} s (flushing)" if var == "nitrate" else ""))
    if "stop" in est:
        say(f"  firings paired within +-{PAIR_DZ:g} m: {int((t['fire_dz'].abs() <= PAIR_DZ).sum())}; samples per stop "
            f"(after flushing) median {t.loc[t['stop_n'] > 0, 'stop_n'].median():g}")
        e_l = t["stop_early_minus_late"].dropna()
        if var == "nitrate" and len(e_l):
            say(f"  flushing signature, early minus late part of the stop: median {e_l.median():+.3f} {unit} "
                f"(n = {len(e_l)}, sign test p = {sign_test(e_l):.3f})")
    paired = t[(t["fire_dz"].abs() <= PAIR_DZ)] if "fire_dz" in t else t.iloc[0:0]
    fw = paired["fire_w"].dropna() if len(paired) else pd.Series(dtype=float)
    if len(fw):
        fly = int((fw >= MOVING_W).sum())
        say(f"  rosette speed at the firings (+-{FIRE_W_S:g} s): median {fw.median():.2f} m/s; {fly} of {len(fw)} "
            f"fired while moving (|w| >= {MOVING_W:g} m/s)")
        if fly > len(fw) / 2:
            say("    -> bottles fired ON THE FLY: there is no stop to sample, so no stop estimator. The bottle's water")
            say("       passed the sensors at about the firing time, on the up-cast: up@fire is the matching estimator.")
    empty = paired[(paired["stop_n"] == 0) & ~(paired["fire_w"] >= MOVING_W)] if len(paired) else paired
    if len(empty):
        dtp, dpt = empty["stop_dt_at_p"].dropna(), empty["stop_dp_at_t"].dropna()
        say(f"  {len(empty)} of {len(paired)} paired bottles have no samples at their stop. Why, from the samples:")
        if len(dtp):
            say(f"    samples at the firing pressure (+-{STOP_DP:g} dbar) are nearest the firing at a median "
                f"{dtp.median():+.0f} s (range {dtp.min():+.0f} to {dtp.max():+.0f} s); the stop window is "
                f"-{STOP_BEFORE_S:g} to +{STOP_AFTER_S:g} s")
        if len(dpt):
            say(f"    the samples nearest the firing time are a median {dpt.median():+.1f} dbar from the firing pressure "
                f"(range {dpt.min():+.1f} to {dpt.max():+.1f})")
        if len(dtp) >= 3 and dtp.abs().median() > 2 * STOP_BEFORE_S and dtp.std() < 0.25 * dtp.abs().median():
            say(f"    -> a consistent time offset: these samples' clock looks {dtp.median():+.0f} s off the firing times")
    for a, b in (("up", "down"), ("down@fire", "down"), ("up@fire", "up"), ("up@fire", "down@fire"),
                 ("stop", "down@fire"), ("stop", "up@fire")):
        if a in est and b in est:
            if var == "chl":
                d = (both[a] - both[b]) / both[b].abs().clip(lower=0.05)
                say(f"  {a} vs {b}: median relative difference {d.median():+.1%}, median |.| {d.abs().median():.1%}")
            else:                                         # nitrate reaches zero: relative differences mean nothing
                d = both[a] - both[b]
                say(f"  {a} vs {b}: median difference {d.median():+.3f} {unit}, median |.| {d.abs().median():.3f} {unit}")
    result = dict(table=t, estimators=est, fits={}, loocv_rmse={}, paired={}, out=out)
    if len(both) < 5:
        say("  fewer than 5 bottles with every estimator: no calibration comparison")
        return result
    bb = both["bottle"].to_numpy(float)
    lev = 1 / len(bb) + (bb - bb.mean()) ** 2 / np.sum((bb - bb.mean()) ** 2)
    if lev.max() > 0.5:
        i = int(np.argmax(lev))
        say(f"  WARNING: one bottle carries the fit (leverage {lev[i]:.2f} > 0.5: {both['cast'].iloc[i]} "
            f"{both['depth_m'].iloc[i]:g} m, {bb[i]:g} {unit}); leaving it out changes the fit drastically, so the "
            "LOOCV numbers below are not a fair comparison")
    result["max_leverage"] = float(lev.max())
    say(f"\ncalibration on the same {len(both)} bottles (sensor = offset + gain * bottle)")
    say(f"estimator   gain   (95% CI)             offset      R2    LOOCV RMSE ({unit})")
    loo = {}
    for e in est:
        f = ols(both["bottle"], both[e])
        loo[e] = loo_errors(both["bottle"], both[e])
        result["fits"][e] = f
        result["loocv_rmse"][e] = float(np.sqrt(np.nanmean(loo[e] ** 2)))
        say(f"{e:10s}{f['gain']:8.4f} ({f['lo']:.4f} to {f['hi']:.4f})  {f['offset']:+8.4f}  {f['R2']:.3f}   "
            f"{result['loocv_rmse'][e]:.4f}")
    try:
        from scipy.stats import wilcoxon
    except ImportError:
        wilcoxon = None
    say("\npaired LOOCV squared errors vs down (negative = this estimator predicts the bottles better)")
    for e in est:
        if e == "down":
            continue
        d = loo[e] ** 2 - loo["down"] ** 2
        p_sign = sign_test(d)
        p_w = float(wilcoxon(d).pvalue) if (wilcoxon is not None and np.any(d != 0)) else float("nan")
        result["paired"][e] = dict(median=float(np.nanmedian(d)), better=int((d < 0).sum()), n=int(len(d)),
                                   p_sign=p_sign, p_wilcoxon=p_w)
        say(f"  {e:9s}: median {np.nanmedian(d):+.4f}, better on {int((d < 0).sum())} of {len(d)} bottles, "
            f"sign test p = {p_sign:.3f}" + (f", Wilcoxon p = {p_w:.3f}" if np.isfinite(p_w) else ""))
    ref = both["stop"] if "stop" in est else both["up"]
    big = both.reindex((ref - both["down"]).abs().sort_values(ascending=False).index).head(8)
    say("\nlargest differences from down:")
    say(big[[c for c in ("cast", "depth_m", "bottle", "down", "up", "down@fire", "up@fire", "stop", "down_cells",
                         "up_cells", "stop_n", "fire_dz", "fire_w") if c in big]].round(4).to_string(index=False))
    say(f"\nper-bottle table -> {out}")
    return result


# ------------------------------------------------------------------ offset sweep: lag versus rosette geometry
def sweep_deltas(var: str, half_range: float | None = None) -> np.ndarray:
    """Offsets to try: +-3 m in 0.25 m steps for chl, +-10 m in 0.5 m for nitrate; a wider --range in 1 m steps."""
    r = half_range if half_range is not None else (3.0 if var == "chl" else 10.0)
    step = 0.25 if r <= 3 else (0.5 if r <= 10 else 1.0)
    return np.arange(-r, r + step / 1000, step)


def sweep(cruise: str, var: str = "chl", deltas=None, verbose: bool = True, n_boot: int = 500) -> dict:
    """Slide each leg's window by an offset D (positive = deeper than the bottle) and score the match by LOOCV.

    The window is centred on the bottle's firing depth when a firing is paired (else its nominal depth) plus D,
    for D over a range, on the down-cast and the up-cast separately; at each D the down/up cell medians are
    fitted against the bottles as in compare() (forward fit, sensor = offset + gain * bottle) and scored by the
    cross-validated R2, Q2 = 1 - PRESS / SS, where PRESS is the sum of squared leave-one-out errors and SS the sum
    of squares of the sensor values about their mean. Q2 is scale-free on purpose. An error in sensor units (LOOCV
    RMSE) shrinks wherever the shifted window sits in a flat layer, because the sensor barely varies there, and so
    rewards windows that carry no match at all; the inverse form (bottle units) divides by the gain and explodes in
    the same place. Q2 asks what fraction of the sensor's variation across bottles the bottles predict out of
    sample: about 0 or below in a flat layer, whatever its noise. Moving samples only: the
    stationary ones at the stop are flushed and carry no lag. Two physical effects move the best D:
      * a sensor response lag (the flow cell still holding water from earlier in the profile) shows, on the way
        down, water from ABOVE the sensor and, on the way up, water from BELOW it: it pushes the best D down-cast
        deeper (D > 0) and up-cast shallower (D < 0) by about the same distance, speed x response time;
      * the rosette geometry (bottles mounted above the sensors; the firing depth is the sensor's) puts the bottle's
        water shallower than the firing depth on BOTH legs: it shifts both best D the same way.
    With h the geometric offset (m) and tau the response time (s): D_down = h + v_down tau, D_up = h - v_up tau, where
    v_down and v_up are the net profiling speeds, measured from the samples (depth travelled over time while
    moving). Hence tau = (D_down - D_up) / (v_down + v_up) and h = D_down - v_down tau. (Averaging the two D
    would leak lag into h whenever the legs run at different speeds.) A first-order lag also smooths the profile,
    which no shift can undo, so tau comes out LOW: on synthetic first-order lags of 3, 5 and 8 s the sweep gave
    2.2, 3.9 and 5.9 s (about three quarters), with |h| under 1 m. Read tau as a lower bound. Each best D is
    reported with the range of D whose Q2 is within 0.05 of its maximum, and the bottles are resampled
    (bootstrap, 500 draws) to give 90% intervals for both best D, tau and h: with 15 to 50 bottles and a jagged
    curve those intervals, not the single best values, are the honest answer. The bootstrap covers sampling noise,
    not the method's own bias: because tau comes out low, h = D_down - v_down tau comes out high when the lag is
    several seconds (synthetic, true h = 0: +0.21, -0.13, +0.81 m for tau = 3, 5, 8 s; with noise +0.6 m, 90%
    interval +0.15 to +1.49 m). For a lagged sensor, read h to about +-0.8 m at best. Longer lags are worse: a
    shift cannot undo the smoothing, which grows with tau, so tau saturates (true 12, 16, 20, 25 s gave 6.7, 7.3,
    9.5, 9.5 s on the synthetic cruise, +-30 m sweep) and h's bias reaches +2.2 m (true 0). Beyond about 8 s the
    sweep gives a floor for tau and no usable h."""
    if deltas is None:
        deltas = sweep_deltas(var)
    say = print if verbose else (lambda *a, **k: None)
    bottles, firings = load_bottles(cruise, var), load_firings(cruise)
    src = (pp.product_dir(cruise, "ctd_derived"), "{}_der.cnv", read_ctd_scans) if var == "chl" else \
          (pp.product_dir(cruise, "nitrate_qc"), "{}_nitrate_QC.csv", read_suna_frames)
    cache, centres = {}, []
    for _, r in bottles.iterrows():
        cast, z0 = r["cast"], float(r["depth_m"])
        if cast not in cache:
            f = src[0] / src[1].format(cast)
            cache[cast] = src[2](f) if f.exists() else None
        zc = z0
        if firings is not None:
            fb = firings[firings["cast"] == cast]
            if len(fb):
                j = (fb["depth_m"] - z0).abs().idxmin()
                if abs(float(fb.at[j, "depth_m"]) - z0) <= PAIR_DZ:
                    zc = float(fb.at[j, "depth_m"])
        centres.append((cast, zc, float(r["bottle"])))
    speed = {leg: net_speed(cache.values(), leg) for leg in ("down", "up")}
    deltas = np.asarray(deltas, float)
    B = np.array([c[2] for c in centres], float)                   # bottle values
    V = {}                                                          # leg -> (bottle, D) matched sensor values
    for leg in ("down", "up"):
        V[leg] = np.full((len(centres), len(deltas)), np.nan)
        for i, (cast, zc, _) in enumerate(centres):
            s = cache.get(cast)
            if s is None:
                continue
            m = s[s[leg] & s["moving"]]
            for j, D in enumerate(deltas):
                V[leg][i, j] = cell_median(m, zc + D)[0]
    rows = []
    for leg in ("down", "up"):
        for j, D in enumerate(deltas):
            k = np.isfinite(V[leg][:, j])
            b, v = B[k], V[leg][k, j]
            if len(b) >= 5 and np.ptp(b) > 0 and np.ptp(v) > 0:
                e = loo_forward_errors(b, v)
                press, ss = float(np.nansum(e ** 2)), float(np.sum((v - v.mean()) ** 2))
                rows.append(dict(leg=leg, D=float(D), n=len(b), q2=1.0 - press / ss,
                                 loocv_rmse=float(np.sqrt(np.nanmean(e ** 2))), R2=ols(b, v)["R2"]))
    t = pd.DataFrame(rows)
    res = dict(table=t, best={}, speed=speed)
    if t.empty:
        say("  sweep: too few bottles")
        return res
    unit = "mg/m3" if var == "chl" else "umol/L"
    say(f"\noffset sweep, {cruise} {var}: window centred at the firing depth (else nominal) + D; D > 0 = deeper")
    wide = t.pivot(index="D", columns="leg", values="q2")
    nn = t.pivot(index="D", columns="leg", values="n")
    say("      D     down Q2 (n)     up Q2 (n)    [cross-validated R2 of the forward fit; 1 = perfect, <= 0 = no skill]")
    for D in wide.index:
        cells = []
        for leg in ("down", "up"):
            if leg in wide and np.isfinite(wide.at[D, leg]):
                cells.append(f"{wide.at[D, leg]:+8.3f} ({int(nn.at[D, leg]):2d})")
            else:
                cells.append(" " * 13)
        say(f"  {D:+6.2f}  " + "   ".join(cells))
    for leg in ("down", "up"):
        tl = t[t["leg"] == leg]
        if tl.empty:
            continue
        tl = tl[tl["n"] >= 0.8 * tl["n"].max()]                  # compare D values on nearly the same bottles
        i = tl["q2"].idxmax()
        near = tl[tl["q2"] >= tl.at[i, "q2"] - 0.05]["D"]
        res["best"][leg] = dict(D=float(tl.at[i, "D"]), q2=float(tl.at[i, "q2"]), rmse=float(tl.at[i, "loocv_rmse"]),
                                lo=float(near.min()), hi=float(near.max()))
        say(f"  best D {leg:4s}: {tl.at[i, 'D']:+.2f} m (Q2 {tl.at[i, 'q2']:+.3f}, n {int(tl.at[i, 'n'])}); "
            f"Q2 within 0.05 of that: {near.min():+.2f} to {near.max():+.2f} m")
    edge = [leg for leg, b in res["best"].items() if b["D"] in (float(deltas.min()), float(deltas.max()))]
    res["edge"] = edge
    if edge:
        say(f"  NO INTERIOR OPTIMUM for {', '.join(edge)}: the best D sits at the end of the range "
            f"({deltas.min():+g} to {deltas.max():+g} m), so it is a limit, not an optimum. tau and h are not computed;")
        say("  widen the range (--range) to see the whole curve. If there is still no optimum, the simple model (a lag")
        say("  of a few seconds plus a fixed offset) does not describe this sensor.")
    elif {"down", "up"} <= set(res["best"]) and np.isfinite(speed["down"]) and np.isfinite(speed["up"]):
        dd, du = res["best"]["down"]["D"], res["best"]["up"]["D"]
        vd, vu = speed["down"], speed["up"]
        tau = (dd - du) / (vd + vu)
        h = dd - vd * tau
        res["tau_s"], res["geometry_m"] = tau, h
        say(f"  profiling speed: down {vd:.2f} m/s, up {vu:.2f} m/s (net, while moving)")
        say(f"  -> response time tau = (D_down - D_up)/(v_down + v_up) = {tau:+.1f} s "
            f"(positive: the sensor shows water it has already passed; a lower bound, about 3/4 of a first-order "
            f"lag's true time constant in synthetic tests)")
        say(f"  -> geometric offset h = D_down - v_down tau = {h:+.2f} m (negative: the matching water is shallower "
            f"than the firing depth, as for bottles mounted above the sensors)")
        if tau > 8:
            say(f"  NOTE tau > 8 s: past what the shift can measure. On synthetic first-order lags of 12, 16, 20 and 25 s")
            say("  the sweep gave 6.7, 7.3, 9.5 and 9.5 s, and h = +0.4 to +2.2 m with no geometric offset at all. Read")
            say("  tau as a floor (the true time constant is longer, by an unknown amount) and h as unresolved.")
        if tau < 0:
            say("  NOTE tau < 0: the up-cast's best window is deeper than the down-cast's. A response lag cannot do that")
            say("  (it delays, it does not anticipate), so unless the interval below reaches well above 0 these are not")
            say("  physical estimates: the lag-plus-offset model does not describe this sensor.")
        ci = _bootstrap_sweep(B, V, deltas, vd, vu, n_boot=n_boot)
        res["bootstrap"] = ci
        if ci:
            say(f"  bootstrap over bottles ({ci['n_ok']} of {n_boot} resamples usable), 90% intervals:")
            say(f"     best D down {ci['down'][0]:+.2f} to {ci['down'][1]:+.2f} m, up {ci['up'][0]:+.2f} to {ci['up'][1]:+.2f} m;"
                f" tau {ci['tau'][0]:+.1f} to {ci['tau'][1]:+.1f} s; h {ci['h'][0]:+.2f} to {ci['h'][1]:+.2f} m")
        say("     These intervals, not the single best values, are what the data support. They cover sampling noise only:")
        say("     with a lag of several seconds h is also biased high (up to about +0.8 m in synthetic tests).")
    return res


def _bootstrap_sweep(B, V, deltas, vd, vu, n_boot=500, seed=0) -> dict:
    """Resample the bottles with replacement, find each leg's best D again (largest in-sample R2 of the forward fit,
    scale-free like the sweep's Q2, on D values with at least 80% of the bottles), and turn the pair into tau and h;
    returns 5th and 95th percentiles. The resampling asks how much the best D would move with a different draw of bottles from the same
    water: with 15 to 50 bottles and a jagged curve, that is the main uncertainty."""
    rng = np.random.default_rng(seed)
    nb = len(B)
    out = {"down": [], "up": [], "tau": [], "h": []}
    for _ in range(n_boot):
        idx = rng.integers(0, nb, nb)
        best = {}
        for leg in ("down", "up"):
            Vi, Bi = V[leg][idx], B[idx]
            n = np.isfinite(Vi).sum(axis=0)
            if n.max() < 5:
                break
            score = np.full(len(deltas), -np.inf)
            for j in np.where(n >= max(5, 0.8 * n.max()))[0]:
                k = np.isfinite(Vi[:, j])
                b, v = Bi[k], Vi[k, j]
                if np.ptp(b) == 0 or np.ptp(v) == 0:
                    continue
                g, o = np.polyfit(b, v, 1)
                score[j] = 1.0 - np.sum((v - (o + g * b)) ** 2) / np.sum((v - v.mean()) ** 2)
            if not np.isfinite(score).any():
                break
            best[leg] = deltas[int(np.argmax(score))]
        if len(best) < 2:
            continue
        tau = (best["down"] - best["up"]) / (vd + vu)
        out["down"].append(best["down"]); out["up"].append(best["up"])
        out["tau"].append(tau); out["h"].append(best["down"] - vd * tau)
    if len(out["tau"]) < 0.5 * n_boot:
        return {}
    ci = {k: (float(np.percentile(v, 5)), float(np.percentile(v, 95))) for k, v in out.items()}
    ci["n_ok"] = len(out["tau"])
    return ci


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cruise")
    ap.add_argument("--var", choices=["chl", "nitrate"], default="chl")
    ap.add_argument("--flush-s", type=float, default=None,
                    help="nitrate: seconds discarded at the start of each stop (default 10; chl default 0)")
    ap.add_argument("--sweep", action="store_true",
                    help="also slide each leg's window by an offset and score it (separates lag from rosette geometry)")
    ap.add_argument("--range", type=float, default=None,
                    help="sweep: half-width in metres (default 3 for chl, 10 for nitrate; wider ranges step 1 m)")
    a = ap.parse_args()
    compare(a.cruise, a.var, a.flush_s)
    if a.sweep:
        sweep(a.cruise, a.var, deltas=sweep_deltas(a.var, a.range))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
