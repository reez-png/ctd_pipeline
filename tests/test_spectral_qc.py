"""suncat.nitrate.spectral_qc: the in-water population for step01's out-of-water gate and reference spectrum comes
from frames taken while the CTD was in the water. Synthetic data only."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from suncat.nitrate import spectral_qc as Q

REPO = Path(__file__).resolve().parents[1]
WL = np.linspace(190.0, 395.0, 256)


def _cnv(path, start="Apr 16 2026 11:00:00", cond=True, n=2400, wet=(240, 2160)):
    """0.25 s scans; wet (conductivity 5 S/m, pressure > 1) between scan indices wet[0] and wet[1]."""
    t = np.arange(n) * 0.25
    w = (np.arange(n) >= wet[0]) & (np.arange(n) < wet[1])
    cols = [("timeS", t), ("prdM", np.where(w, 20.0, 0.1))] + ([("c0S/m", np.where(w, 5.0, 0.002))] if cond else [])
    head = [f"# name {i} = {c}: {c}" for i, (c, _) in enumerate(cols)]
    body = [" ".join(f"{v:.4f}" for v in row) for row in zip(*[v for _, v in cols])]
    Path(path).write_text("\n".join([f"# start_time = {start} [Instrument's time stamp, header]", *head, "*END*", *body]))


def test_wet_intervals_from_conductivity_or_pressure(tmp_path):
    _cnv(tmp_path / "ZZ_CTD_01.cnv")
    _cnv(tmp_path / "ZZ_CTD_02.cnv", start="Apr 16 2026 12:00:00", cond=False)
    (tmp_path / "ZZ_CTD_03.cnv").write_text("# name 0 = timeS: t\n*END*\n1.0\n")
    iv, skipped = Q.ctd_wet_intervals(sorted(tmp_path.glob("*.cnv")))
    t0 = (pd.Timestamp("2026-04-16 11:00:00") - pd.Timestamp("1970-01-01")).total_seconds()
    r = iv.set_index("cast")
    assert r.loc["ZZ_CTD_01", "t_in"] == pytest.approx(t0 + 60.0) and r.loc["ZZ_CTD_01", "t_out"] == pytest.approx(t0 + 539.75)
    assert r.loc["ZZ_CTD_01", "by"] == "conductivity" and r.loc["ZZ_CTD_02", "by"] == "pressure"
    assert r.loc["ZZ_CTD_02", "t_in"] == pytest.approx(t0 + 3600 + 60.0)
    assert [s[0] for s in skipped] == ["ZZ_CTD_03.cnv"]


def test_in_cast_keeps_a_margin_and_which_cast_does_not():
    iv = pd.DataFrame(dict(cast=["A", "B"], t_in=[100.0, 1000.0], t_out=[400.0, 1090.0]))
    t = np.array([30, 90, 120, 170, 330, 390, 470, 1010, 1050])
    assert Q.in_cast(t, iv, 60.0).tolist() == [False, False, False, True, True, False, False, False, False]  # B too short
    assert Q.which_cast(t, iv).tolist() == ["", "", "A", "A", "A", "A", "", "B", "B"]
    assert Q.off_cast(t, iv, 60.0).tolist() == [True, False, False, False, False, False, True, False, False]
    assert not Q.off_cast(t, iv.iloc[:0], 60.0).any()


def _file(n_deck=420, n_water=180, seed=3):
    """A file that is mostly deck time with a wet optical path (peak ~38,000, tilted spectrum), then one cast in
    seawater (peak ~45,000). Frame times 1 s apart from 2026-04-16 11:00; the water frames are the last n_water."""
    rng = np.random.default_rng(seed)
    n = n_deck + n_water
    water = np.arange(n) >= n_deck
    lamp = np.exp(-((WL - 260.0) / 70.0) ** 2)
    tilt = 1.0 + 0.35 * (WL - 300.0) / 100.0
    spec = np.where(water[:, None], 45000 * lamp, 38000 * lamp * tilt) * (1 + rng.normal(0, 0.002, (n, WL.size)))
    df = pd.DataFrame(dict(time_utc=pd.Timestamp("2026-04-16 11:00:00") + pd.to_timedelta(np.arange(n), unit="s"),
                           nitrate_uM=1.0, abs254=np.where(water, -0.43, -0.35) + rng.normal(0, 0.003, n),
                           abs350=0.0, bromide=0.0, spec_avg=np.where(water, 18500.0, 17000.0) + rng.normal(0, 40, n),
                           dark=600.0, fit_rmse=2e-4, humidity=2.0))
    t = (df["time_utc"] - pd.Timestamp("1970-01-01")).dt.total_seconds().to_numpy()
    iv = pd.DataFrame(dict(cast=["C1"], t_in=[t[n_deck]], t_out=[t[-1]]))
    return df, spec, water, Q.in_cast(t, iv, 30.0)


def test_deck_majority_inverts_the_old_rule_and_not_the_new_one():
    df, spec, water, inc = _file()
    cfg = Q.QCConfig()
    old, _, info_old = Q.qc_frames(df, spec, WL, cfg)
    new, win, info = Q.qc_frames(df, spec, WL, cfg, inc)
    assert not new["f_offcast"].any() and not old["f_offcast"].any()
    assert info_old["window_source"] == "file" and info["window_source"] == "ctd" and info["ref_source"] == "ctd"
    assert old.loc[water, "L1_clean"].mean() < 0.1                     # the old rule rejects the seawater frames
    assert new.loc[water, "L1_clean"].mean() > 0.95                    # the new one keeps them
    deck = ~water & ~new["f_warmup"].to_numpy()
    assert new.loc[deck, "f_outofwater"].mean() > 0.95                 # and flags the deck frames
    med = float(np.median(new.loc[water, "peak_counts"]))
    assert 38000 < win["peak"][0] <= 0.9 * med + 1 and win["peak"][1] == cfg.peak_hard[1]   # dark side only


def test_too_few_in_cast_frames_keep_the_file_rule():
    df, spec, water, inc = _file(n_deck=100, n_water=400)
    few = inc & (np.cumsum(inc) <= 10)                                  # only 10 frames inside a cast
    cfg = Q.QCConfig()
    a, wa, ia = Q.qc_frames(df, spec, WL, cfg)
    b, wb, ib = Q.qc_frames(df, spec, WL, cfg, few)
    assert ib["window_source"] == "file" and ib["ref_source"] == "file" and wa == wb
    assert a["L1_clean"].equals(b["L1_clean"])
    _, _, ic = Q.qc_frames(df, spec, WL, Q.QCConfig(oow_adaptive=False), inc)
    assert ic["window_source"] == "fixed"


def test_step01_runs_on_the_example_cruise(tmp_path):
    pytest.importorskip("nbclient"); pytest.importorskip("ipykernel")
    shutil.copytree(REPO / "cruises" / "_EXAMPLE", tmp_path / "cruises" / "_EXAMPLE")
    (tmp_path / "metadata" / "calibration").mkdir(parents=True)
    shutil.copy(REPO / "metadata" / "calibration" / "SNA2503C.cal", tmp_path / "metadata" / "calibration")
    env = dict(os.environ, CTD_PIPELINE_ROOT=str(tmp_path), MPLBACKEND="Agg")
    nb = REPO / "processing_scripts" / "notebook" / "suna" / "step01_suna_spectral_qc.ipynb"
    r = subprocess.run([sys.executable, str(REPO / "tools" / "run_notebook.py"), str(nb), "--set", "CRUISE_ID=_EXAMPLE",
                        "--quiet", "--timeout", "300"], env=env, capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    audit = tmp_path / "cruises" / "_EXAMPLE" / "_audit" / "suna"
    qc = pd.read_csv(audit / "_EXAMPLE__suna_L1_qc_summary.csv", encoding="utf-8-sig")
    assert int(qc["clean"].sum()) == 543 and int(qc["n"].sum()) == 633           # README_EXAMPLE: unchanged
    cmp = pd.read_csv(audit / "_EXAMPLE__suna_oow_gate_comparison.csv", encoding="utf-8-sig")
    assert set(cmp["source"]) == {"ctd"} and (cmp["clean_old"] - cmp["clean_new"]).abs().max() <= 10
    win = pd.read_csv(audit / "_EXAMPLE__suna_oow_windows.csv", encoding="utf-8-sig")
    assert {"source", "n_frames", "reference"} <= set(win.columns)


def test_ctd_window_is_one_sided_with_a_minimum_width():
    """Seawater frames that read brighter, or a few per cent darker (turbidity, a soak), stay in water; the deck
    state (15 % darker, abs254 + 0.08) does not."""
    df, spec, water, inc = _file(n_deck=0, n_water=600)
    rows = {300: (1.06, 0.0), 310: (0.94, 0.02), 320: (0.85, 0.0), 330: (1.0, 0.08), 340: (0.95, 0.04)}
    for i, (f, da) in rows.items():
        df.loc[i, "spec_avg"] *= f; df.loc[i, "abs254"] += da; spec[i] *= f
    out, win, info = Q.qc_frames(df, spec, WL, Q.QCConfig(), inc | True)
    oow = out["f_outofwater"]
    assert info["window_source"] == "ctd"
    assert not oow[300] and not oow[310] and not oow[340]          # brighter; 6 % darker; 5 % darker and +0.04
    assert oow[320] and oow[330]                                   # 15 % darker; abs254 + 0.08
    med = float(np.median(out["peak_counts"]))
    assert win["peak"][0] <= 0.9 * med + 1 and win["abs254"][1] >= float(np.median(out["abs254"])) + 0.05 - 1e-9


def test_offcast_frames_are_not_clean_and_not_out_of_water():
    df, spec, water, inc = _file(n_deck=0, n_water=300)
    off = np.zeros(len(df), bool); off[100:120] = True
    out, _, _ = Q.qc_frames(df, spec, WL, Q.QCConfig(), inc | True, off)
    assert out.loc[100:119, "f_offcast"].all() and not out.loc[100:119, "L1_clean"].any()
    assert not out.loc[100:119, "f_outofwater"].any()              # the optical flag is untouched


def _cnv_break(path, start="Apr 16 2026 11:00:00", n=2400):
    """0.25 s scans; conductivity wet from scan 240 to 2160 throughout, but the package back at the surface
    (pressure 0.0) from scan 600 to 1200 (150 s): P45_07 cast 20's shape."""
    i = np.arange(n); t = i * 0.25
    wet = (i >= 240) & (i < 2160)
    p = np.where(wet, 20.0, 0.1); p[(i >= 600) & (i < 1200)] = 0.0
    cols = [("timeS", t), ("prdM", p), ("c0S/m", np.where(wet, 5.0, 0.002))]
    head = [f"# name {k} = {c}: {c}" for k, (c, _) in enumerate(cols)]
    body = [" ".join(f"{v:.4f}" for v in row) for row in zip(*[v for _, v in cols])]
    Path(path).write_text("\n".join([f"# start_time = {start} [Instrument's time stamp, header]", *head, "*END*", *body]))


def test_a_surface_break_inside_the_wet_interval_is_listed_and_left_out(tmp_path):
    _cnv_break(tmp_path / "ZZ_CTD_20.cnv")
    _cnv(tmp_path / "ZZ_CTD_19.cnv", start="Apr 16 2026 12:00:00")
    iv, _ = Q.ctd_wet_intervals(sorted(tmp_path.glob("*.cnv")))
    r = iv.set_index("cast")
    t0 = (pd.Timestamp("2026-04-16 11:00:00") - pd.Timestamp("1970-01-01")).total_seconds()
    assert r.loc["ZZ_CTD_20", "t_in"] == pytest.approx(t0 + 60) and r.loc["ZZ_CTD_20", "t_out"] == pytest.approx(t0 + 539.75)
    (a, b), = r.loc["ZZ_CTD_20", "breaks"]
    assert a == pytest.approx(t0 + 150) and b == pytest.approx(t0 + 299.75)
    assert r.loc["ZZ_CTD_19", "breaks"] == []
    t = t0 + np.array([130, 200, 320, 400, 470])
    assert Q.in_cast(t, iv, 30.0).tolist() == [False, False, False, True, True]  # 130 and 320: within 30 s of the break
    assert Q.in_cast(t, iv.drop(columns="breaks"), 30.0).tolist() == [True, True, True, True, True]
    assert Q.surface_breaks(np.arange(100.0), np.r_[np.full(50, 5.0), np.full(20, 0.1), np.full(30, 5.0)], 0, 99) == []


def _two_casts(level2=0.85, n1=300, n2=200, seed=5):
    """One file, two casts: cast A in seawater at peak ~45,000; cast B's seawater reads level2 x that (another
    water mass, the lamp, the optical path), with a few deck-state frames (15 % darker than B's own level)."""
    rng = np.random.default_rng(seed)
    n = n1 + n2
    lamp = np.exp(-((WL - 260.0) / 70.0) ** 2)
    lvl = np.r_[np.ones(n1), np.full(n2, level2)]
    deck = np.zeros(n, bool); deck[n1 + 50:n1 + 60] = True
    lvl[deck] *= 0.85
    spec = (45000 * lvl)[:, None] * lamp * (1 + rng.normal(0, 0.002, (n, WL.size)))
    df = pd.DataFrame(dict(time_utc=pd.Timestamp("2026-08-03 08:00:00") + pd.to_timedelta(np.arange(n), unit="s"),
                           nitrate_uM=1.0, abs254=np.where(deck, -0.35, -0.43) + rng.normal(0, 0.003, n),
                           abs350=0.0, bromide=0.0, spec_avg=18500.0 * lvl + rng.normal(0, 40, n),
                           dark=600.0, fit_rmse=2e-4, humidity=2.0))
    cast_of = np.r_[np.full(n1, "A"), np.full(n2, "B")].astype(object)
    return df, spec, deck, cast_of


def test_a_cast_reading_darker_seawater_keeps_its_frames_with_its_own_window():
    df, spec, deck, co = _two_casts()
    inc = np.ones(len(df), bool)
    file_only, _, _ = Q.qc_frames(df, spec, WL, Q.QCConfig(oow_per_cast=False), inc, None, co)
    per_cast, win, info = Q.qc_frames(df, spec, WL, Q.QCConfig(), inc, None, co)
    b = (co == "B") & ~deck
    assert file_only.loc[b, "f_outofwater"].mean() > 0.9          # 15 % below the file's level: out of water
    assert not per_cast.loc[b, "f_outofwater"].any()                # against its own level: in water
    assert per_cast.loc[deck, "f_outofwater"].all()                 # B's deck frames are still caught
    cw = {w["cast"]: w for w in info["cast_windows"]}
    assert cw["A"]["source"] == cw["B"]["source"] == "cast"
    assert cw["B"]["peak_counts_median"] / cw["B"]["peak_counts_file_median"] == pytest.approx(0.85, abs=0.03)
    assert per_cast.loc[co == "A", "f_outofwater"].sum() == file_only.loc[co == "A", "f_outofwater"].sum()


def test_a_cast_with_too_few_frames_keeps_the_file_window():
    df, spec, deck, co = _two_casts(level2=1.0, n2=20)
    df2, _, info = Q.qc_frames(df, spec, WL, Q.QCConfig(), np.ones(len(df), bool), None, co)
    cw = {w["cast"]: w for w in info["cast_windows"]}
    assert cw["A"]["source"] == "cast" and cw["B"]["source"] == "file" and cw["B"]["n"] == 20
