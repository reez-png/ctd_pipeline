"""suncat.nitrate.water_entry: the SUNA's water entry from its own peak counts, against the CTD's first wet scan.
Synthetic frames shaped like P45_06 and P45_07 (Nov 2026)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from suncat.nitrate import water_entry as W

T_IN, T_OUT = 1000.0, 1600.0


def _cast(deck=38000, water=45000, clock=16.0, dt=1.0, gaps=None, extra=(), seed=0, abs_deck=0.08, abs_water=0.0):
    """SUNA frames from 800 to 1800 s (SUNA clock = CTD clock + clock): deck level until the SUNA enters at
    T_IN + clock, seawater until it leaves at T_OUT + clock + 8 (water held in the path), deck after."""
    rng = np.random.default_rng(seed)
    t = np.arange(800.0, 1800.0, dt)
    if gaps:                                         # bursts of `gaps[0]` frames every `gaps[1]` s
        n, every = gaps
        t = np.concatenate([b + np.arange(n) * dt for b in np.arange(800.0, 1800.0, every)])
    wet = (t >= T_IN + clock) & (t <= T_OUT + clock + 8)
    peak = np.where(wet, water, deck) * (1 + rng.normal(0, 0.003, t.size))
    a254 = np.where(wet, abs_water, abs_deck) + rng.normal(0, 0.002, t.size)
    for at, val in extra:                           # splash frames
        i = int(np.argmin(np.abs(t - at)))
        peak[i] = val
    return pd.DataFrame(dict(t=t, peak_counts=peak, abs254=a254))


IV = pd.DataFrame(dict(cast=["X"], t_in=[T_IN], t_out=[T_OUT], t_first=[T_IN - 40], t_last=[T_OUT + 40]))


def test_dry_path_entry_at_one_hertz():
    e = W.entry_table(_cast(deck=25800), IV).iloc[0]
    assert e["stable_lo"] < 16 <= e["stable_hi"] and e["stable_hi"] - e["stable_lo"] == pytest.approx(1.0)
    assert e["contact_s"] == pytest.approx(e["stable_s"])            # no splash: contact is stable water
    assert 20 < e["exit_s"] < 26                                      # held water: exit late by about 8 s


def test_wet_path_with_a_splash_before_stable_water():
    e = W.entry_table(_cast(extra=[(T_IN + 14.0, 2568), (T_IN + 15.0, 65535)]), IV).iloc[0]
    assert e["contact_lo"] < 14 <= e["contact_hi"]                   # first contact: the splash
    assert e["stable_lo"] < 16 <= e["stable_hi"]                     # stable water once in
    assert e["deck_peak"] == pytest.approx(38000, rel=0.01)


def test_a_warming_lamp_does_not_drag_first_change_back():
    """The deck level rises about 8 % in the 2 minutes before entry (lamp warming after switch-on, P45_07 casts 02
    and 05): first change must stay at the entry, not run back through the drifting deck frames."""
    d = _cast(deck=25800)
    pre = d["t"] < T_IN + 16
    d.loc[pre, "peak_counts"] *= 1 + 0.08 * np.clip((d.loc[pre, "t"] - (T_IN - 104)) / 120, 0, 1)
    e = W.entry_table(d, IV).iloc[0]
    assert e["stable_lo"] < 16 <= e["stable_hi"] and e["stable_s"] - e["contact_s"] <= 2
    old = W.entry_table(d, IV, depart_frac=0.03, contact_max_s=60, max_gap_s=1e9, deck_window=(120.0, 20.0)).iloc[0]
    assert old["stable_s"] - old["contact_s"] > 15                 # the first version ran back through the drift


def test_bursts_twelve_seconds_apart_give_a_wide_bracket():
    e = W.entry_table(_cast(clock=4.0, gaps=(1, 12.0)), IV).iloc[0]
    assert e["stable_lo"] < 4 <= e["stable_hi"] and 10 <= e["stable_hi"] - e["stable_lo"] <= 13


def test_a_ctd_record_that_starts_wet_or_too_few_frames_is_left_out():
    iv = IV.assign(t_first=T_IN - 3)
    e = W.entry_table(_cast(), iv).iloc[0]
    assert np.isnan(e["contact_s"]) and "starts wet" in e["note"]
    e = W.entry_table(_cast().query("t > 1100"), IV).iloc[0]         # SUNA switched on after the entry
    assert np.isnan(e["stable_s"]) and e["note"]


def test_summary_measured_and_undetermined():
    rows = []
    for i, clk in enumerate([15.2, 16.1, 15.8, 16.4, 15.9, 16.0, 16.7]):
        e = W.entry_table(_cast(clock=clk, seed=i), IV).iloc[0]
        rows.append(e)
    s = W.summarize(pd.DataFrame(rows), "stable_s")
    assert s["verdict"] == "MEASURED" and abs(s["median"] - 16.0) < 1.0
    rows = []
    for i, clk in enumerate([-8.0, 0.5, 10.0, 6.0, 21.0, -3.0]):
        rows.append(W.entry_table(_cast(clock=clk, seed=i, gaps=(1, 12.0)), IV).iloc[0])
    s = W.summarize(pd.DataFrame(rows), "stable_s")
    assert s["verdict"].startswith("UNDETERMINED") and "wide" in s["verdict"]
    s = W.summarize(pd.DataFrame(rows[:3]), "stable_s")
    assert "only 3 cast" in s["verdict"]


def test_step04_runs_on_the_example_cruise(tmp_path):
    """step01 -> step02_suna_colocation -> step04 on the shipped example (its CTD records start wet, so the entry
    check must say so and propose no clock)."""
    import os, shutil, subprocess, sys
    from pathlib import Path
    pytest.importorskip("nbclient"); pytest.importorskip("ipykernel")
    repo = Path(__file__).resolve().parents[1]
    shutil.copytree(repo / "cruises" / "_EXAMPLE", tmp_path / "cruises" / "_EXAMPLE")
    (tmp_path / "metadata" / "calibration").mkdir(parents=True)
    shutil.copy(repo / "metadata" / "calibration" / "SNA2503C.cal", tmp_path / "metadata" / "calibration")
    env = dict(os.environ, CTD_PIPELINE_ROOT=str(tmp_path), MPLBACKEND="Agg")
    for nb in ("step01_suna_spectral_qc", "step02_suna_colocation", "step04_suna_time_lag"):
        r = subprocess.run([sys.executable, str(repo / "tools" / "run_notebook.py"),
                            str(repo / "processing_scripts" / "notebook" / "suna" / f"{nb}.ipynb"),
                            "--set", "CRUISE_ID=_EXAMPLE", "--quiet", "--timeout", "300"],
                           env=env, capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, f"{nb}:\n" + r.stdout[-3000:] + r.stderr[-3000:]
    ev = pd.read_csv(tmp_path / "cruises" / "_EXAMPLE" / "_audit" / "suna" / "time_lag" / "_EXAMPLE__suna_clock_events.csv")
    assert len(ev) == 3 and ev["entry_offset_s"].isna().all() and ev["note"].str.contains("starts wet").all()


def test_a_surface_break_is_kept_out_of_the_water_level():
    """P45_07 cast 20: back at the surface for 2.5 min inside the CTD's wet interval. Without the break the deck
    frames (38,000) would be over half the frames the water level is taken from."""
    d = _cast(deck=38000)
    brk = (d["t"] > T_IN + 100) & (d["t"] < T_IN + 420)
    d.loc[brk, "peak_counts"] = 38000.0
    br = [(T_IN + 84, T_IN + 404)]                                # CTD times (the SUNA stamps 16 s late here)
    iv = IV.copy(); iv["breaks"] = pd.Series([br], dtype=object)
    e = W.entry_table(d, iv).iloc[0]
    assert e["water_peak"] == pytest.approx(45000, rel=0.01)
    assert e["stable_lo"] < 16 <= e["stable_hi"]
    bad = W.entry_table(d, IV).iloc[0]
    assert bad["water_peak"] < 40000                              # the deck state took over the median
