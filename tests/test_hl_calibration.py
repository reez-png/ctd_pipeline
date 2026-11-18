"""suncat.nitrate.calibration: the homogeneous-layer bottle calibration needs clean SUNA frames on BOTH legs at a
bottle depth (P45_06: 5 bottles passed the old 3-in-total rule on 1 or 2 frames of one leg)."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from suncat.nitrate import calibration as HC

REPO = Path(__file__).resolve().parents[1]


def _frames(cast, depth, down, up, n_down, n_up):
    p = np.r_[np.linspace(depth - 2, depth + 2, n_down), np.linspace(depth - 2, depth + 2, n_up)]
    return pd.DataFrame(dict(cast=cast, pressure_db=p, leg=["down"] * n_down + ["up"] * n_up,
                             nitrate_uM=np.r_[np.full(n_down, down), np.full(n_up, up)], L1_clean=True))


FR = pd.concat([_frames("A", 55, 14.0, 15.2, 1, 12),       # one down-cast frame: cannot show homogeneity
                _frames("B", 80, 15.0, 15.4, 8, 9),        # homogeneous, both legs well sampled
                _frames("C", 30, 6.0, 11.0, 6, 6),         # stratified
                _frames("D", 5, 0.1, 0.0, 0, 9),           # up-cast only
                _frames("E", 85, 13.0, 13.5, 5, 5)], ignore_index=True)
BOT = pd.DataFrame(dict(cast_id=list("ABCDE"), depth_m=[55, 80, 30, 5, 85], nitrate_uM=[16.0, 17.0, 9.0, 0.0, 15.0],
                        flag=["", "", "", "", "suspect"]))


def test_each_leg_needs_its_own_frames():
    c = HC.homogeneous_pairs(FR, BOT).set_index("cast")
    assert not c.loc["A", "used"] and "down 1 / up 12" in c.loc["A", "note"]
    assert c.loc["B", "used"] and c.loc["B", "n_down"] == 8 and c.loc["B", "n_up"] == 9
    assert not c.loc["C", "used"] and c.loc["C", "note"].startswith("stratified")
    assert c.loc["D", "note"] == "one leg only"
    assert c.loc["E", "note"] == "flagged bottle excluded"


def test_the_old_rule_is_reproduced_when_asked():
    c = HC.homogeneous_pairs(FR, BOT, min_total=3).set_index("cast")
    assert c.loc["A", "used"] and c.loc["B", "used"] and not c.loc["C", "used"]
    assert c.loc["A", "suna"] == pytest.approx(15.2)                 # the pooled median is the up-cast's


def test_leg_mean_value_and_the_fit():
    c = HC.homogeneous_pairs(FR, BOT, value="legs").set_index("cast")
    assert c.loc["B", "suna"] == pytest.approx(15.2)                 # (15.0 + 15.4) / 2
    r = HC.fit_line([0, 5, 10, 17], [-0.4, 4.05, 8.5, 14.73])
    assert r["gain"] == pytest.approx(0.89, abs=1e-3) and r["offset"] == pytest.approx(-0.4, abs=1e-2)
    assert HC.fit_line([1.0], [1.0])["n"] == 1 and np.isnan(HC.fit_line([1.0], [1.0])["gain"])


def test_leave_one_out_shows_a_fit_resting_on_one_bottle():
    """P45_07's shape: four surface bottles and two near 16 uM that disagree; dropping either moves the gain."""
    fit = pd.DataFrame(dict(cast=list("abcdef"), depth_m=[5, 5, 5, 6, 55, 55],
                            bottle=[6.3, 3.0, 1.0, 0.1, 16.27, 15.89], suna=[5.6, 2.4, 0.6, -0.7, 13.9, 15.75]))
    loo = HC.leave_one_out(fit).set_index("cast")
    assert loo.loc["f", "gain"] < 0.92 < HC.fit_line(fit["bottle"], fit["suna"])["gain"]
    assert loo.loc["e", "gain"] > 1.0


def test_the_notebook_runs_on_the_example_cruise(tmp_path):
    pytest.importorskip("nbclient"); pytest.importorskip("ipykernel")
    shutil.copytree(REPO / "cruises" / "_EXAMPLE", tmp_path / "cruises" / "_EXAMPLE")
    (tmp_path / "metadata" / "calibration").mkdir(parents=True)
    shutil.copy(REPO / "metadata" / "calibration" / "SNA2503C.cal", tmp_path / "metadata" / "calibration")
    env = dict(os.environ, CTD_PIPELINE_ROOT=str(tmp_path), MPLBACKEND="Agg", SUNCAT_HEADLESS="1")
    for nb in ("step01_suna_spectral_qc", "step02_suna_colocation", "step03_suna_bottle_calibration"):
        r = subprocess.run([sys.executable, str(REPO / "tools" / "run_notebook.py"),
                            str(REPO / "processing_scripts" / "notebook" / "suna" / f"{nb}.ipynb"),
                            "--set", "CRUISE_ID=_EXAMPLE", "--timeout", "300"],
                           env=env, capture_output=True, text=True, timeout=600)
        assert r.returncode == 0, f"{nb}:\n" + r.stdout[-3000:] + r.stderr[-3000:]
    assert "check: the fit under each rule" in r.stdout
    cal = pd.read_csv(tmp_path / "cruises" / "_EXAMPLE" / "_audit" / "suna" / "_EXAMPLE__suna_bottle_calibration.csv")
    assert {"n_down", "n_up", "used", "note"} <= set(cal.columns)
    assert (cal.loc[cal["used"], ["n_down", "n_up"]] >= HC.MIN_PER_LEG).all().all()
