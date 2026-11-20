"""The README's quick start: the production path runs on cruises/_EXAMPLE and recovers what was planted in it.

Planted (cruises/_EXAMPLE/README_EXAMPLE.md): the SUNA reads SUNA = -1.10 + 0.879 x true nitrate; the fluorometer
reads 1.14 x true chlorophyll. Until Nov 2026 only the SUNA QA path ran on the example: step02 found no SUNA files on
a case-sensitive file system (A9000001.CSV against "*.csv"), step03 had no suna_calibration.json, step07 no Derived
files, and step08 no bin-type header. Runs the real notebooks against a copy (CTD_PIPELINE_ROOT); nothing is written
into the repo. Skipped where the notebook runner (nbclient, ipykernel) is not installed.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("nbclient")
pytest.importorskip("ipykernel")

REPO = Path(__file__).resolve().parents[1]
NB = REPO / "processing_scripts" / "notebook"
STEPS = ["suna/step01_suna_spectral_qc", "step02_merge_ctd_suna", "step03_nitrate_qc",
         "productivity/step07_productivity_calibration", "step08_grid_l3", "step04_plot_profiles",
         "step06_hydrographic_sections"]


@pytest.fixture(scope="module")
def example(tmp_path_factory):
    root = tmp_path_factory.mktemp("example_root")
    shutil.copytree(REPO / "cruises" / "_EXAMPLE", root / "cruises" / "_EXAMPLE")
    shutil.copytree(REPO / "metadata" / "calibration", root / "metadata" / "calibration")
    env = dict(os.environ, CTD_PIPELINE_ROOT=str(root), MPLBACKEND="Agg", SUNCAT_HEADLESS="1")
    out = {}
    for step in STEPS:
        r = subprocess.run([sys.executable, str(REPO / "tools" / "run_notebook.py"), str(NB / f"{step}.ipynb"),
                            "--set", "CRUISE_ID=_EXAMPLE", "--timeout", "300"],
                           env=env, capture_output=True, text=True, timeout=900)
        assert r.returncode == 0, f"{step} failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}"
        out[step] = r.stdout
    return root / "cruises" / "_EXAMPLE", out


def test_step02_pairs_every_cast(example):
    _, out = example
    assert "Paired casts: 3   Unpaired CTD 1 s files: 0" in out["step02_merge_ctd_suna"]


def test_step03_recovers_the_planted_suna_calibration(example):
    cruise, _ = example
    fit = pd.read_csv(cruise / "_audit" / "step03_fit_coefficients.csv").set_index("stage")
    a1, a0 = float(fit.loc["stage3_applied_calibration", "a1"]), float(fit.loc["stage3_applied_calibration", "a0"])
    assert a1 == pytest.approx(0.879, abs=0.01) and a0 == pytest.approx(-1.10, abs=0.1)
    assert len(list((cruise / "L2" / "nitrate").glob("*_nitrate_QC.csv"))) == 3


def test_step07_recovers_the_planted_chlorophyll_gain(example):
    _, out = example
    line = next(ln for ln in out["productivity/step07_productivity_calibration"].splitlines()
                if ln.strip().startswith("gain "))
    assert float(line.split()[1]) == pytest.approx(1.14, abs=0.03)


def test_step08_grids_nitrate_and_chlorophyll_and_step04_plots_the_grid(example):
    cruise, out = example
    tables = sorted((cruise / "L3" / "gridded").glob("*_L3_1m_down.csv"))
    assert len(tables) == 3
    for t in tables:
        d = pd.read_csv(t)
        assert d["nitrate_uM"].notna().any() and d["chl_mg_m3"].notna().any()
        assert d["chl_npq_risk"].fillna(False).astype(bool).any()        # daylight casts: the top bins at risk
    assert out["step04_plot_profiles"].count("nitrate: gridded (step08)") == 3
    readable = pd.read_csv(cruise / "L3" / "CTD" / "profiles_readable" / "_EXAMPLE_CTD_03_readable.csv")
    assert {"Nitrate_corrected_uM", "Chl_calibrated_mg_m3", "Nitrate_downup_over_2uM"} <= set(readable.columns)
