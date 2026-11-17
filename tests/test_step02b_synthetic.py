"""step02b (SUNA response time) on synthetic cruises with a known lag: the real notebook, run with
tools/run_notebook.py against a temporary data root. It must recover tau with no bottles, say the lag explains the
hysteresis, and decide from the bottle statistics (usable when the bottles carry no other error, not usable when
they do). Skipped where the notebook runner (nbclient, ipykernel) is not installed."""
from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("nbclient")
pytest.importorskip("ipykernel")

REPO = Path(__file__).resolve().parents[1]
NB = REPO / "processing_scripts" / "notebook" / "suna" / "step02b_suna_response_correction.ipynb"
CRUISE = "ZZ12_03"
TAU = 25.0


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


GEN = _load("make_synthetic_cruise", "tools/make_synthetic_cruise.py")
TRT = _load("test_response_time", "tests/test_response_time.py")       # its continuous-lag integrator


def _cruise(root: Path, bottle_noise: float) -> Path:
    """Synthetic cruise with step02's L2 co-located frames (raw SUNA = truth seen through a 25 s lag, + noise),
    bottles on the cruise's stated scale (adopted a1 = 0.89, a0 = -0.39), and that calibration file."""
    C = GEN.make_cruise(root, CRUISE)
    rng = np.random.default_rng(3)
    a1, a0 = 0.89, -0.39
    rows, bots = [], []
    for f in sorted((C / "L2/nitrate").glob("*_nitrate_QC.csv")):
        d = pd.read_csv(f)
        cast = f.name.replace("_nitrate_QC.csv", "")
        tt = pd.to_datetime(d["utc_time"], format="ISO8601")
        t = (tt - tt.iloc[0]).dt.total_seconds().to_numpy(float)
        z = d["depth_m"].to_numpy(float)
        truth = lambda zz: a0 + a1 * GEN.no3_true(zz)                     # the sensor's raw scale
        raw = TRT.lagged(t, z, TAU, truth) + rng.normal(0, 0.3, len(d))
        rows.append(pd.DataFrame({"time_utc": d["utc_time"], "nitrate_uM": raw, "L1_clean": True,
                                  "pressure_db": d["pressure_db"], "cast": cast, "in_cast": True, "leg": d["leg"]}))
        for zb in range(20, 61, 5):
            if zb < z.max() - 2:
                bots.append(dict(cast_id=cast, depth_m=float(zb), flag="",
                                 nitrate_uM=float(GEN.no3_true(zb)) + rng.normal(0, bottle_noise)))
    l2 = C / "L2" / "SUNA"
    l2.mkdir(parents=True, exist_ok=True)
    pd.concat(rows).to_csv(l2 / f"{CRUISE}__suna_L2_colocated.csv", index=False)
    (C / "metadata/bottle_nitrate").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(bots).to_csv(C / "metadata/bottle_nitrate" / f"{CRUISE}_bottle_nitrate.csv", index=False)
    (C / "metadata" / "suna_calibration.json").write_text(json.dumps(
        {"mode": "adopted", "a1": a1, "a0": a0, "source": "synthetic"}))
    return C


def _run(root: Path):
    env = dict(os.environ, CTD_PIPELINE_ROOT=str(root), MPLBACKEND="Agg")
    cmd = [sys.executable, str(REPO / "tools" / "run_notebook.py"), str(NB), "--set", f"CRUISE_ID={CRUISE}",
           "--set", "LAT_DEG=0.0", "--quiet", "--timeout", "600"]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, f"step02b failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}"
    return r


def _summary(C: Path) -> dict:
    return json.loads((C / "_audit" / "suna" / f"{CRUISE}__suna_response_summary.json").read_text())


def test_step02b_recovers_tau_and_finds_usable_pairs_when_bottles_are_clean(tmp_path):
    C = _cruise(tmp_path, bottle_noise=0.05)
    _run(tmp_path)
    s = _summary(C)
    assert abs(s["tau_median"] - TAU) <= 1.5, s["tau_median"]
    assert s["lag_explained"] and s["rss_share_explained"] > 0.8
    assert s["pairs_usable"]
    f = s["fits"]["N@tau"]
    assert f["lo"] <= 0.89 <= f["hi"] or abs(f["gain"] - 0.89) < 0.03         # the scale the bottles were put on
    assert s["fits"]["down"]["gain"] < f["gain"]                                # the lag compresses the down-cast
    assert (C / "_audit" / "suna" / f"{CRUISE}__suna_response_time.csv").exists()


def test_step02b_keeps_the_calibration_when_the_bottles_scatter(tmp_path):
    # the lag is still found and explained, but bottles with 3 uM of their own error cannot carry a depth-matched
    # calibration: the verdict must say so, from R2, not from the lag
    C = _cruise(tmp_path, bottle_noise=3.0)
    _run(tmp_path)
    s = _summary(C)
    assert abs(s["tau_median"] - TAU) <= 1.5 and s["lag_explained"]
    assert not s["pairs_usable"] and s["fits"]["N@tau"]["R2"] < 0.90
