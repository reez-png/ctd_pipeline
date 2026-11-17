"""tools/diagnose_suna_lag.py: tell a time-constant response from a flushing (per metre travelled) response by
how each scales across casts profiled at different speeds, and read the clock from step04's water entry."""
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
DSL = _load("diagnose_suna_lag", "tools/diagnose_suna_lag.py")
CRUISE = "ZZ12_03"
STRETCH = (0.8, 1.25, 2.0, 3.0)        # cast k's clock runs this much slower: speeds from about 1.1 to 0.3 m/s


def _response(t, z, kind, par, delay, f, sub=20):
    """The sensor's reading: water from z(t - delay), seen through a first-order response per second (kind
    "time", par = tau s) or per metre travelled (kind "distance", par = L m), integrated 20 times finer."""
    tf = np.concatenate([np.linspace(t[i], t[i + 1], sub, endpoint=False) for i in range(len(t) - 1)] + [t[-1:]])
    zf = np.interp(tf - delay, t, z)
    x = f(zf)
    sf = np.concatenate([[0.0], np.cumsum(np.abs(np.diff(zf)))])
    y = np.empty_like(x)
    y[0] = x[0]
    for i in range(1, len(x)):
        step = tf[i] - tf[i - 1] if kind == "time" else sf[i] - sf[i - 1]
        y[i] = y[i - 1] + (x[i] - y[i - 1]) * (1 - np.exp(-step / par))
    return np.interp(t, tf, y)


def _cruise(root, kind, par, delay=16.0, noise=0.2):
    C = GEN.make_cruise(root, CRUISE)
    rng = np.random.default_rng(2)
    for k, f in enumerate(sorted((C / "L2/nitrate").glob("*_nitrate_QC.csv"))):
        d = pd.read_csv(f)
        tt = pd.to_datetime(d["utc_time"], format="ISO8601")
        t = (tt - tt.iloc[0]).dt.total_seconds().to_numpy(float) * STRETCH[k]   # same path, slower or faster
        d["utc_time"] = (tt.iloc[0] + pd.to_timedelta(t, unit="s")).astype(str)
        d["no3_qc"] = _response(t, d["depth_m"].to_numpy(float), kind, par, delay, GEN.no3_true) \
            + rng.normal(0, noise, len(d))
        d["L1_clean"] = True
        d.to_csv(f, index=False)
    return C


@pytest.mark.parametrize("kind, par", [("time", 8.0), ("distance", 4.0)])
def test_the_speed_dependence_tells_time_from_flushing(tmp_path, monkeypatch, kind, par):
    _cruise(tmp_path, kind, par)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    r = DSL.response(CRUISE, 16.0, say=lambda *a, **k: None)
    T = r["table"]
    assert r["mode"] == kind, (r["verdict"], r["slope"])
    if kind == "time":
        assert np.allclose(T["tau_s"], par, atol=1.5), T["tau_s"].tolist()
        assert abs(r["slope"]["b"]) < 0.3
    else:
        assert np.allclose(T["L_m"], par, atol=0.75), T["L_m"].tolist()
        assert abs(r["slope"]["b"] + 1) < 0.3


def test_clock_uses_water_entry_and_reports_exit(tmp_path, monkeypatch):
    C = GEN.make_cruise(tmp_path, CRUISE)
    tl = C / "_audit" / "suna" / "time_lag"
    tl.mkdir(parents=True)
    n = 12
    pd.DataFrame({"cast": [f"{CRUISE}_CTD_{i:03d}" for i in range(n)],
                  "entry_offset_s": 16 + np.random.default_rng(0).normal(0, 1, n),
                  "exit_offset_s": 22 + np.random.default_rng(1).normal(0, 4, n)}       # wet path lingers at exit
                 ).to_csv(tl / f"{CRUISE}__suna_clock_events.csv", index=False)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    ck = DSL.clock(CRUISE, say=lambda *a, **k: None)
    assert abs(ck["entry_median"] - 16) < 1 and ck["entry_ci"][0] < 16 < ck["entry_ci"][1]
    assert ck["exit_minus_entry"] > 3 and ck["p_sign"] < 0.05


def test_conditions_test_names_the_property_that_goes_with_tau():
    rng = np.random.default_rng(1)          # a seed whose random speeds are not correlated with tau by chance
    n = 20
    t_range = rng.uniform(1, 15, n)
    T = pd.DataFrame({"tau_s": 3 + 2.5 * t_range + rng.normal(0, 1, n), "tau_edge": False,
                      "speed": rng.uniform(0.4, 0.7, n), "max_depth": rng.uniform(30, 300, n), "t_range": t_range,
                      "t_min": 25 - t_range, "duration": rng.uniform(100, 1200, n)})
    c = DSL.conditions_test(T, say=lambda *a, **k: None)
    assert c["t_range"]["rho"] > 0.8 and c["t_range"]["p"] < 0.01
    assert c["t_min"]["rho"] < -0.8
    assert abs(c["speed"]["rho"]) < 0.5


def test_a_clock_already_applied_is_not_counted_twice(tmp_path, monkeypatch):
    # frames written by step03 after a -16 s clock correction carry suna_time_offset_s; step04's events still
    # report the raw +16 s, so the delay the diagnosis applies must be the remainder, 0
    C = _cruise(tmp_path, "time", 8.0, delay=0.0)
    for f in sorted((C / "L2/nitrate").glob("*_nitrate_QC.csv")):
        d = pd.read_csv(f)
        d["suna_time_offset_s"] = -16.0
        d.to_csv(f, index=False)
    tl = C / "_audit" / "suna" / "time_lag"
    tl.mkdir(parents=True)
    pd.DataFrame({"cast": [f"{CRUISE}_CTD_{i:03d}" for i in range(1, 5)], "entry_offset_s": [15.9, 16.1, 16.2, 16.0],
                  "exit_offset_s": [16.5, 15.5, 16.3, 16.1]}).to_csv(tl / f"{CRUISE}__suna_clock_events.csv", index=False)
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    out = DSL.diagnose(CRUISE, bottles=False, say=lambda *a, **k: None)
    assert out["delay"] == 0.0
    assert np.allclose(out["response"]["table"]["tau_s"], 8.0, atol=1.5)
