"""tools/estimate_response_time.py: a first-order lag fitted to both legs must recover a known response time from
synthetic casts (heave, stops, noise), including lags far longer than the offset sweep can measure."""
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
ERT = _load("estimate_response_time", "tools/estimate_response_time.py")
CRUISE = "ZZ12_03"


def lagged(t, z, tau, f, sub=20):
    """The true continuous first-order response to f(z(t)), integrated on a step 20 times finer than the frames
    (not the estimator's own discrete formula, so the test does not grade the model against itself)."""
    tf = np.concatenate([np.linspace(t[i], t[i + 1], sub, endpoint=False) for i in range(len(t) - 1)] + [t[-1:]])
    zf = np.interp(tf, t, z)
    x = f(zf)
    if tau <= 0:
        return np.interp(t, tf, x)
    y = np.empty_like(x)
    y[0] = x[0]
    for i in range(1, len(x)):
        y[i] = y[i - 1] + (x[i] - y[i - 1]) * (1 - np.exp(-(tf[i] - tf[i - 1]) / tau))
    return np.interp(t, tf, y)


def _lagged_cruise(root, tau, noise=0.5, shift_up=0.0, seed=1):
    C = GEN.make_cruise(root, CRUISE)
    rng = np.random.default_rng(seed)
    for f in sorted((C / "L2/nitrate").glob("*_nitrate_QC.csv")):
        d = pd.read_csv(f)
        tt = pd.to_datetime(d["utc_time"], format="ISO8601")
        t = (tt - tt.iloc[0]).dt.total_seconds().to_numpy(float)
        z = d["depth_m"].to_numpy(float) + np.where(d["leg"].astype(str).eq("up"), shift_up, 0.0)
        d["no3_qc"] = lagged(t, z, tau, GEN.no3_true) + rng.normal(0, noise, len(d))
        d["L1_clean"] = True
        d.to_csv(f, index=False)
    return C


@pytest.mark.parametrize("tau", [0.0, 5.0, 20.0, 45.0])
def test_recovers_a_known_response_time(tmp_path, monkeypatch, tau):
    _lagged_cruise(tmp_path, tau)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    r = ERT.estimate(CRUISE, "nitrate", verbose=False, write=False)
    assert abs(r["tau_median"] - tau) <= max(0.5, 0.05 * tau), (tau, r["table"]["tau"].tolist())
    lo, hi = r["tau_ci90"]
    assert lo - 0.5 <= tau <= hi + 0.5 or abs(r["tau_median"] - tau) <= 0.3
    assert abs(r["tau_pooled"] - tau) <= max(0.5, 0.05 * tau)


def test_lag_explains_the_hysteresis(tmp_path, monkeypatch):
    # with a 20 s lag the up-cast reads systematically differently from the down-cast; the fitted lag removes it
    _lagged_cruise(tmp_path, 20.0)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    r = ERT.estimate(CRUISE, "nitrate", verbose=False, write=False)
    t = r["table"]
    assert (t["rms"] < 0.6 * t["rms0"]).all()
    assert r["rss_share_explained"] > 0.5
    assert (t["hyst"].abs() < 0.25 * t["hyst0"].abs()).all()


def test_a_displaced_up_cast_biases_tau_by_displacement_over_speed(tmp_path, monkeypatch):
    # the documented failure: water displaced by d between the legs looks like a lag of about d / (v_down + v_up)
    _lagged_cruise(tmp_path, 16.0, shift_up=2.0)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    r = ERT.estimate(CRUISE, "nitrate", verbose=False, write=False)
    bias = r["tau_median"] - 16.0
    assert 0.3 <= bias <= 3.0, bias                     # 2 m / (0.8 + 1.0 m/s) = 1.1 s


def test_too_long_a_lag_is_an_edge_not_an_estimate(tmp_path, monkeypatch):
    _lagged_cruise(tmp_path, 40.0)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    r = ERT.estimate(CRUISE, "nitrate", tau_max=20.0, verbose=False, write=False)
    assert r["table"]["edge"].all() and "tau_median" not in r


def test_writes_the_table_and_the_curve(tmp_path, monkeypatch):
    _lagged_cruise(tmp_path, 5.0)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    r = ERT.estimate(CRUISE, "nitrate", verbose=False)
    assert r["out"].exists() and r["out"].with_name(r["out"].stem + "_curve.csv").exists()
    assert len(pd.read_csv(r["out"])) == 4


def test_chlorophyll_on_ctd_scans_without_a_lag(tmp_path, monkeypatch):
    GEN.make_cruise(tmp_path, CRUISE)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    r = ERT.estimate(CRUISE, "chl", tau_max=4.0, step=0.5, verbose=False, write=False)
    assert r["tau_median"] <= 0.5


def test_design_matrix_matches_the_continuous_lag():
    # the exact interval formula against a fine-step integration, on a path with a stop and a reversal
    t = np.arange(0, 200, 1.0)
    z = np.concatenate([np.linspace(2, 80, 100), np.full(20, 80.0), np.linspace(80, 2, 80)])
    grid = np.arange(0, 82, 1.0)
    N = GEN.no3_true(grid)
    for tau in (0.0, 3.0, 25.0):
        y = ERT.lag_design(t, z, grid, tau) @ N
        ref = lagged(t, z, tau, lambda zz: np.interp(zz, grid, N), sub=50)
        assert np.max(np.abs(y - ref)) < 0.02, tau


def test_lag_correction_matches_the_bottles_better_than_the_down_cast(tmp_path, monkeypatch):
    # a 20 s lag: the down-cast reads water from above the bottle; the lag-corrected profile reads the bottle's own
    # water, so it calibrates against the bottles with gain 1 and a smaller LOOCV error than down or N@0
    C = _lagged_cruise(tmp_path, 20.0, noise=0.3)
    rows = [(f"{CRUISE}_CTD_{n:03d}", z) for n in (2, 3, 4) for z in range(20, 61, 4)]
    zb = np.array([r[1] for r in rows], float)
    pd.DataFrame({"cast_id": [r[0] for r in rows], "depth_m": zb, "bottle_no3": GEN.no3_true(zb)}
                 ).to_csv(C / "_audit" / "step03_bottle_matchups.csv", index=False)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    r = ERT.estimate(CRUISE, "nitrate", verbose=False, write=False, bottles=True)["bottles"]
    L = r["loocv_rmse"]
    assert L["N@tau_cruise"] < 0.5 * L["down"] and L["N@tau_cruise"] < 0.5 * L["N@0"], L
    f = r["fits"]["N@tau_cruise"]
    assert f["lo"] <= 1.0 <= f["hi"] or abs(f["gain"] - 1) < 0.03, f
    assert r["paired"]["N@tau_cruise"]["better"] >= 0.75 * r["paired"]["N@tau_cruise"]["n"]


def test_offset_sweep_finds_a_bottle_depth_offset(tmp_path, monkeypatch):
    # bottles whose water came from 2 m above their recorded depth (a firing or rosette offset), on a lagged sensor:
    # the lag-corrected profile has no lag left, so the offset sweep must find D = -2 m
    C = _lagged_cruise(tmp_path, 20.0, noise=0.3)
    rows = [(f"{CRUISE}_CTD_{n:03d}", z) for n in (2, 3, 4) for z in range(24, 61, 3)]
    zb = np.array([r[1] for r in rows], float)
    pd.DataFrame({"cast_id": [r[0] for r in rows], "depth_m": zb, "bottle_no3": GEN.no3_true(zb - 2.0)}
                 ).to_csv(C / "_audit" / "step03_bottle_matchups.csv", index=False)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    o = ERT.estimate(CRUISE, "nitrate", verbose=False, write=False, bottles=True)["bottles"]["offset"]
    assert abs(o["best_D"] + 2.0) <= 0.5, o["best_D"]
    assert o["ci90"][0] - 0.5 <= -2.0 <= o["ci90"][1] + 0.5
    assert o["q2_best"] > o["q2_at_0"] and o["determined"]
