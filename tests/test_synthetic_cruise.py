"""The phase-4 path (step07 -> step08 -> step06) on synthetic cruises that differ from P45 wherever the code could
assume it (tools/make_synthetic_cruise.py: other cruise id, three-digit cast files, prDM, a cast without Loop Edit
flags, mean-pressure bins, southern hemisphere east of Greenwich, an equatorial line). Runs the real notebooks with
tools/run_notebook.py against a temporary data root (CTD_PIPELINE_ROOT); nothing is written into the repo.
Skipped where the notebook runner (nbclient, ipykernel) is not installed.
"""
from __future__ import annotations

import importlib.util
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
NB = REPO / "processing_scripts" / "notebook"
SOUTH, EQ = "ZZ12_03", "EQ01_01"


def _gen():
    spec = importlib.util.spec_from_file_location("make_synthetic_cruise", REPO / "tools" / "make_synthetic_cruise.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


GEN = _gen()


@pytest.fixture(scope="module")
def root(tmp_path_factory):
    r = tmp_path_factory.mktemp("synthetic_root")
    GEN.make_cruise(r, SOUTH, lat0=-12.4, lon0=40.6, utc_offset=3.0)
    GEN.make_cruise(r, EQ, lat0=0.45, lon0=-3.0, utc_offset=0.0)
    return r


def run(root, notebook, cruise, *sets):
    env = dict(os.environ, CTD_PIPELINE_ROOT=str(root), MPLBACKEND="Agg")
    cmd = [sys.executable, str(REPO / "tools" / "run_notebook.py"), str(notebook), "--set", f"CRUISE_ID={cruise}",
           "--quiet", "--timeout", "300"]
    for s in sets:
        cmd += ["--set", s]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, f"{notebook.name} failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}"


@pytest.fixture(scope="module")
def south(root):
    run(root, NB / "productivity" / "step07_productivity_calibration.ipynb", SOUTH, "INTERACTIVE_ZOOM=False")
    run(root, NB / "step08_grid_l3.ipynb", SOUTH)
    return root / "cruises" / SOUTH


def test_step07_resolves_every_cast_and_recovers_the_calibration(south):
    cal = pd.read_csv(south / "_audit" / "productivity" / f"{SOUTH}__chl_calibration.csv", encoding="utf-8-sig")
    assert set(cal["cast"]) == {f"{SOUTH}_CTD_{n:03d}" for n in (1, 2, 3, 4)}   # numbers AND full hyphenated ids
    u = cal[cal["used"].astype(str).str.lower() == "true"]
    gain, offset = np.polyfit(u["bottle"].astype(float), u["sensor"].astype(float), 1)
    assert abs(gain - GEN.GAIN) < 0.06 and abs(offset - GEN.OFFSET) < 0.03


def test_step07_npq_follows_the_sun(south):
    cal = pd.read_csv(south / "_audit" / "productivity" / f"{SOUTH}__chl_calibration.csv", encoding="utf-8-sig")
    shallow = cal[cal["depth_m"] <= 10].set_index("cast")["npq"].astype(str).str.lower().eq("true")
    # solar start hours at ~41 E: 08.7, 11.7, 14.7, 20.7 -> the first three casts are daylit
    assert shallow.to_dict() == {f"{SOUTH}_CTD_001": True, f"{SOUTH}_CTD_002": True,
                                 f"{SOUTH}_CTD_003": True, f"{SOUTH}_CTD_004": False}


def test_step08_bins_carry_the_ctd_scan_counts_and_the_truth(south):
    for n, daylit in ((1, True), (2, True), (3, True), (4, False)):
        d = pd.read_csv(south / "L3" / "gridded" / f"{SOUTH}_CTD_{n:03d}_L3_1m_down.csv", encoding="utf-8-sig")
        live = d["chl_n"] > 0
        assert live.all() and (d.loc[live, "chl_n"] == d.loc[live, "ctd_n_scans"]).all()
        assert abs(np.median(d["chl_mg_m3"] - GEN.chl_true(d["depSM"]))) < 0.02
        assert abs(np.nanmedian(d["nitrate_uM"] - GEN.no3_true(d["depSM"]))) < 0.1
        risk = d["chl_npq_risk"].astype(str).str.lower().eq("true")
        assert risk.any() == daylit and not risk[d["pressure_db"] > 10].any()


def test_step06_keeps_equatorial_stations(root):
    run(root, NB / "step06_hydrographic_sections.ipynb", EQ, "GROUP_BY=all")
    m = pd.read_csv(root / "cruises" / EQ / "_audit" / "sections" / f"{EQ}__sampling_map.csv", encoding="utf-8-sig")
    assert len(m) == 4 and m["lat"].between(0.0, 1.0).all()                  # 0.10 to 0.45 N, none dropped
    assert (root / "cruises" / EQ / "_audit" / "sections" / f"{EQ}__section__temperature.csv").exists()
