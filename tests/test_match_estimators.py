"""tools/compare_match_estimators.py on a synthetic cruise whose water is the same on both legs and at the stops:
every estimator must recover the known calibration, and the statistics must behave on known inputs."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


GEN = _load("make_synthetic_cruise", "tools/make_synthetic_cruise.py")
CME = _load("compare_match_estimators", "tools/compare_match_estimators.py")
CRUISE = "ZZ12_03"


@pytest.fixture(scope="module")
def cruise(tmp_path_factory):
    root = tmp_path_factory.mktemp("match_root")
    C = GEN.make_cruise(root, CRUISE)
    fires = pd.read_csv(C / "metadata/bottle_nitrate" / f"{CRUISE}_bottle_ctd.csv")
    fires = fires[fires["depth_m"] > 3]                                   # the 2 m stop is above the profiles
    rng = np.random.default_rng(11)
    nominal = fires["depth_m"].round()                                    # bottle tables carry nominal depths
    (C / "_audit/productivity").mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cast": fires["cast_id"], "station": "T", "depth_m": nominal,
                  "bottle": GEN.chl_true(fires["depth_m"]) * (1 + rng.normal(0, 0.03, len(fires))),
                  "npq": False, "used": True, "note": ""}
                 ).to_csv(C / "_audit/productivity" / f"{CRUISE}__chl_calibration.csv", index=False)
    pd.DataFrame({"cast_id": fires["cast_id"], "depth_m": nominal,
                  "bottle_no3": GEN.no3_true(fires["depth_m"]) + rng.normal(0, 0.2, len(fires))}
                 ).to_csv(C / "_audit" / "step03_bottle_matchups.csv", index=False)
    return root


@pytest.mark.parametrize("var, truth_gain, tol, max_rmse",
                         [("chl", GEN.GAIN, 0.05, 0.05), ("nitrate", 1.0, 0.05, 0.8)])
def test_every_estimator_recovers_the_known_calibration(cruise, monkeypatch, var, truth_gain, tol, max_rmse):
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(cruise))
    r = CME.compare(CRUISE, var, verbose=False)
    assert r["estimators"] == ["down", "up", "down@fire", "up@fire", "stop"]
    t = r["table"]
    assert (t["fire_dz"].abs() <= CME.PAIR_DZ).all() and (t["stop_n"] > 0).all()
    for e in r["estimators"]:
        assert abs(r["fits"][e]["gain"] - truth_gain) < tol * truth_gain, (e, r["fits"][e])
    # LOOCV error bounded in bottle units (mg/m3, umol/L). Not equal across estimators: 1 Hz SUNA frames leave 1 m
    # cells sparse, so in the nitracline the down-cast cell median is less precise than the stop at the closing depth
    assert all(v < max_rmse for v in r["loocv_rmse"].values()), r["loocv_rmse"]
    assert r["out"].exists()


def test_nitrate_stop_discards_the_flushing_seconds(cruise, monkeypatch):
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(cruise))
    n0 = CME.compare(CRUISE, "nitrate", flush_s=0, verbose=False)["table"]["stop_n"]
    n10 = CME.compare(CRUISE, "nitrate", flush_s=10, verbose=False)["table"]["stop_n"]
    assert (n10 < n0).all() and (n10 > 0).all()            # stationary samples of 20 to 30 s stops at 1 Hz: ~10 fewer


def test_without_firing_times_only_down_and_up(cruise, monkeypatch, tmp_path):
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(cruise))
    f = cruise / "cruises" / CRUISE / "metadata/bottle_nitrate" / f"{CRUISE}_bottle_ctd.csv"
    f.rename(tmp_path / f.name)
    try:
        assert CME.compare(CRUISE, "chl", verbose=False)["estimators"] == ["down", "up"]
    finally:
        (tmp_path / f.name).rename(f)


def test_stop_survives_loop_edit_rejecting_the_stationary_scans(cruise, tmp_path):
    # Loop Edit flags scans slower than its minimum velocity, so at a stop step03's ctd_scan_pressure_db is empty;
    # the stop must use the frames' own CTD pressure
    C = cruise / "cruises" / CRUISE
    f = pd.read_csv(C / "L2/nitrate" / f"{CRUISE}_CTD_001_nitrate_QC.csv")
    up = f["leg"].eq("up")
    f.loc[up, "ctd_scan_pressure_db"] = np.nan
    f.loc[up, "ctd_scan_kept"] = False
    g = tmp_path / "frames.csv"
    f.to_csv(g, index=False)
    s = CME.read_suna_frames(g)
    fire = pd.read_csv(C / "metadata/bottle_nitrate" / f"{CRUISE}_bottle_ctd.csv").iloc[1]
    r = CME.stop_value(s, pd.Timestamp(fire["time_utc"], tz="UTC"), float(fire["pressure_dbar"]), 10.0)
    assert r["stop_n"] > 0 and np.isfinite(r["stop"])


@pytest.mark.parametrize("var", ["chl", "nitrate"])
def test_sweep_finds_no_offset_where_there_is_none(cruise, monkeypatch, var):
    # synthetic water has no response lag and the sensors sit at the firing depth: best D near 0 on both legs
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(cruise))
    r = CME.sweep(CRUISE, var, deltas=np.arange(-3.0, 3.001, 0.25), verbose=False)
    for leg in ("down", "up"):
        b = r["best"][leg]
        assert b["lo"] <= 0.0 <= b["hi"] or abs(b["D"]) <= 0.5, (leg, b)
    assert abs(r["tau_s"]) <= 1.0 and abs(r["geometry_m"]) <= 0.75
    lo, hi = r["bootstrap"]["h"]
    assert lo <= 0.0 <= hi or max(abs(lo), abs(hi)) <= 0.75                # no offset: 0 inside (or next to) the CI


def test_mixed_timestamp_formats_are_all_read():
    # step03 writes frame times both with and without fractional seconds; pandas' format guess would drop half
    col = pd.Series(["2026-07-29 12:26:01+00:00", "2026-07-29 12:26:01.500000+00:00",
                     "2026-07-29 12:26:02+00:00", "2026-07-29 12:26:02.250000+00:00"])
    t = CME.parse_utc(col)
    assert t.notna().all() and (t.diff().dropna().dt.total_seconds() > 0).all()


def test_sweep_recovers_a_known_response_time(tmp_path_factory, monkeypatch):
    # a first-order sensor lag, tau = 5 s, applied to the synthetic nitrate frames in time order (both legs); bottles
    # spread through the nitracline (20 to 60 m, where the profile changes): the sweep must find tau and no geometry
    root = tmp_path_factory.mktemp("lag_root")
    C = GEN.make_cruise(root, CRUISE)
    tau = 5.0
    for f in sorted((C / "L2/nitrate").glob("*_nitrate_QC.csv")):
        d = pd.read_csv(f)
        tt = pd.to_datetime(d["utc_time"])
        t = (tt - tt.iloc[0]).dt.total_seconds().to_numpy()        # not astype(int64): its unit differs by pandas version
        x = GEN.no3_true(d["depth_m"].to_numpy(float))
        y = np.empty_like(x); y[0] = x[0]
        for i in range(1, len(x)):
            y[i] = y[i - 1] + (x[i] - y[i - 1]) * (1 - np.exp(-(t[i] - t[i - 1]) / tau))
        d["no3_qc"] = y
        d["L1_clean"] = True
        d.to_csv(f, index=False)
    (C / "metadata/bottle_nitrate" / f"{CRUISE}_bottle_ctd.csv").unlink()          # centre on the bottle depths
    rows = [(f"{CRUISE}_CTD_{n:03d}", z) for n in (2, 3, 4) for z in range(20, 61, 4)]
    pd.DataFrame({"cast_id": [r[0] for r in rows], "depth_m": [r[1] for r in rows],
                  "bottle_no3": GEN.no3_true(np.array([r[1] for r in rows], float))}
                 ).to_csv(C / "_audit" / "step03_bottle_matchups.csv", index=False)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(root))
    r = CME.sweep(CRUISE, "nitrate", deltas=np.arange(-10.0, 10.001, 0.5), verbose=False)
    assert r["best"]["down"]["D"] > 0 > r["best"]["up"]["D"]
    assert abs(r["tau_s"] - tau) <= 1.5, (r["tau_s"], r["best"], r["speed"])
    assert abs(r["geometry_m"]) <= 1.0, (r["geometry_m"], r["best"])


def test_sweep_is_not_drawn_to_a_flat_layer(tmp_path_factory, monkeypatch):
    # no lag, no geometry, but a nitracline between two flat layers: a window slid into a flat layer sees almost no
    # sensor variation, so an error in sensor units shrinks there although the match is worthless. The score must
    # be scale-free (cross-validated R2) and keep the best D near 0 on both legs
    root = tmp_path_factory.mktemp("flat_root")
    C = GEN.make_cruise(root, CRUISE)
    rng = np.random.default_rng(5)
    for f in sorted((C / "L2/nitrate").glob("*_nitrate_QC.csv")):
        d = pd.read_csv(f)
        d["no3_qc"] = GEN.no3_true(d["depth_m"].to_numpy(float)) + rng.normal(0, 0.3, len(d))
        d["L1_clean"] = True
        d.to_csv(f, index=False)
    (C / "metadata/bottle_nitrate" / f"{CRUISE}_bottle_ctd.csv").unlink()
    rows = [(f"{CRUISE}_CTD_{n:03d}", z) for n in (2, 3, 4) for z in range(28, 53, 3)]
    zb = np.array([r[1] for r in rows], float)
    pd.DataFrame({"cast_id": [r[0] for r in rows], "depth_m": zb,                  # imperfect bottles, as in real data
                  "bottle_no3": GEN.no3_true(zb) + rng.normal(0, 2.0, len(zb))}
                 ).to_csv(C / "_audit" / "step03_bottle_matchups.csv", index=False)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(root))
    r = CME.sweep(CRUISE, "nitrate", deltas=CME.sweep_deltas("nitrate", 20.0), verbose=False, n_boot=50)
    for leg in ("down", "up"):
        assert abs(r["best"][leg]["D"]) <= 1.0, (leg, r["best"][leg])


def test_motion_copes_with_unsorted_frames_and_missing_times(cruise):
    # real frame files need not be in time order and can lack a time: the speed must still come out per frame
    f = pd.read_csv(cruise / "cruises" / CRUISE / "L2/nitrate" / f"{CRUISE}_CTD_002_nitrate_QC.csv")
    ref = CME.read_suna_frames_df(f)
    g = f.sample(frac=1.0, random_state=1).reset_index(drop=True)
    g.loc[[3, 40], "utc_time"] = None
    s = CME.read_suna_frames_df(g)
    assert s["w"].isna().sum() >= 2 and s.loc[[3, 40], "moving"].eq(False).all()
    m = pd.merge(ref[["t", "w"]].dropna(), s[["t", "w"]].dropna(), on="t", suffixes=("_ref", "_shuffled"))
    assert len(m) > 0.9 * len(ref) and np.allclose(m["w_ref"], m["w_shuffled"], equal_nan=True, atol=0.3)


def test_motion_with_frames_more_than_a_second_apart(cruise):
    # real SUNA frames are often 1 to 2 s apart: every moving frame must still get a speed (a +-1 s window alone
    # left most of them without one, and both legs empty)
    f = pd.read_csv(cruise / "cruises" / CRUISE / "L2/nitrate" / f"{CRUISE}_CTD_003_nitrate_QC.csv")
    rng = np.random.default_rng(4)
    g = f.iloc[np.sort(rng.choice(len(f), size=int(len(f) * 0.6), replace=False))].copy()   # gaps of 1 to 4 s
    jitter = pd.to_timedelta(rng.uniform(0, 0.4, len(g)), unit="s")
    g["utc_time"] = (pd.to_datetime(g["utc_time"], format="ISO8601", utc=True) + jitter).astype(str)
    full, thin = CME.read_suna_frames_df(f), CME.read_suna_frames_df(g)
    for leg in ("down", "up"):
        a = (full[leg] & full["moving"]).sum() / full[leg].sum()
        b = (thin[leg] & thin["moving"]).sum() / thin[leg].sum()
        assert b > 0.8 * a, (leg, a, b)


def test_empty_stop_reports_a_clock_offset(cruise, tmp_path):
    # frames whose clock is 120 s late against the firing times: the stop is empty and the diagnostic says why
    C = cruise / "cruises" / CRUISE
    f = pd.read_csv(C / "L2/nitrate" / f"{CRUISE}_CTD_002_nitrate_QC.csv")
    f["utc_time"] = (pd.to_datetime(f["utc_time"], format="ISO8601", utc=True) + pd.Timedelta(seconds=120)).astype(str)
    s = CME.read_suna_frames_df(f)
    fire = pd.read_csv(C / "metadata/bottle_nitrate" / f"{CRUISE}_bottle_ctd.csv")
    fire = fire[fire["cast_id"] == f"{CRUISE}_CTD_002"].iloc[1]
    r = CME.stop_value(s, pd.Timestamp(fire["time_utc"], tz="UTC"), float(fire["pressure_dbar"]), 10.0)
    assert r["stop_n"] == 0
    assert 85 <= r["stop_dt_at_p"] <= 125                                 # stationary samples begin ~90 s after


def test_bottles_fired_on_the_fly_have_no_stop(cruise, monkeypatch, tmp_path):
    # move every firing to a moment when the rosette is hauling (the up-cast, between stops): samples pass the
    # firing pressure but none is stationary, so there is no stop; fire_w says why and up@fire takes over
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(cruise))
    C = cruise / "cruises" / CRUISE
    f = C / "metadata/bottle_nitrate" / f"{CRUISE}_bottle_ctd.csv"
    orig = pd.read_csv(f)
    rows = []
    try:
        for cast, g in orig.groupby("cast_id"):
            s = CME.read_suna_frames_df(pd.read_csv(C / "L2/nitrate" / f"{cast}_nitrate_QC.csv"))
            m = s[s["up"] & s["moving"]].dropna(subset=["t", "p"]).sort_values("t")
            for z in g["depth_m"]:
                k = (m["z"] - (z - 2.0)).abs().idxmin()     # 2 m above the old stop: hauling, still pairs (PAIR_DZ 3)
                rows.append(dict(cast_id=cast, bottle_no=len(rows) + 1, pressure_dbar=float(m.at[k, "p"]),
                                 depth_m=float(m.at[k, "z"]), time_utc=str(m.at[k, "t"])))
        pd.DataFrame(rows).to_csv(f, index=False)
        r = CME.compare(CRUISE, "nitrate", verbose=False)
    finally:
        orig.to_csv(f, index=False)
    t = r["table"]
    assert "stop" not in r["estimators"] and "up@fire" in r["estimators"]
    assert (t["stop_n"] == 0).all()
    assert (t["fire_w"].dropna() >= CME.MOVING_W).all()


def test_sweep_refuses_tau_and_h_at_the_edge_of_the_range(tmp_path_factory, monkeypatch):
    # a 20 s lag needs offsets of 15 to 20 m; a +-2 m sweep can only hit its edge: no tau, no h
    root = tmp_path_factory.mktemp("edge_root")
    C = GEN.make_cruise(root, CRUISE)
    for f in sorted((C / "L2/nitrate").glob("*_nitrate_QC.csv")):
        d = pd.read_csv(f)
        tt = pd.to_datetime(d["utc_time"], format="ISO8601")
        t = (tt - tt.iloc[0]).dt.total_seconds().to_numpy()
        x = GEN.no3_true(d["depth_m"].to_numpy(float))
        y = np.empty_like(x); y[0] = x[0]
        for i in range(1, len(x)):
            y[i] = y[i - 1] + (x[i] - y[i - 1]) * (1 - np.exp(-(t[i] - t[i - 1]) / 20.0))
        d["no3_qc"] = y
        d["L1_clean"] = True
        d.to_csv(f, index=False)
    (C / "metadata/bottle_nitrate" / f"{CRUISE}_bottle_ctd.csv").unlink()
    rows = [(f"{CRUISE}_CTD_{n:03d}", z) for n in (2, 3, 4) for z in range(20, 61, 4)]
    pd.DataFrame({"cast_id": [r[0] for r in rows], "depth_m": [r[1] for r in rows],
                  "bottle_no3": GEN.no3_true(np.array([r[1] for r in rows], float))}
                 ).to_csv(C / "_audit" / "step03_bottle_matchups.csv", index=False)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(root))
    r = CME.sweep(CRUISE, "nitrate", deltas=CME.sweep_deltas("nitrate", 2.0), verbose=False, n_boot=20)
    assert r["edge"] and "tau_s" not in r and "geometry_m" not in r


def test_cell_window_follows_a_non_integer_depth():
    z = np.arange(20.0, 40.0, 0.05)
    s = pd.DataFrame({"z": z, "v": z})                                    # value = depth: the window's centre shows
    assert CME.cell_median(s, 30.0)[0] == pytest.approx(30.0, abs=0.03)
    assert CME.cell_median(s, 30.3)[0] == pytest.approx(30.3, abs=0.03)   # a fixed whole-metre grid would give 30.0
    assert CME.cell_median(s, 30.3)[1] == 5


def test_statistics_on_known_inputs():
    assert CME.sign_test([1, 1, 1, 1, 1, 1]) == pytest.approx(2 / 64)
    assert CME.sign_test([1, -1]) == 1.0 and np.isnan(CME.sign_test([0, 0]))
    b = np.array([0.1, 0.4, 0.8, 1.2, 2.0, 3.0])
    assert np.allclose(CME.loo_errors(b, 0.04 + 1.6 * b), 0, atol=1e-9)  # a perfect line predicts every bottle
    f = CME.ols(b, 0.04 + 1.6 * b)
    assert f["gain"] == pytest.approx(1.6) and f["offset"] == pytest.approx(0.04) and f["R2"] == pytest.approx(1.0)
