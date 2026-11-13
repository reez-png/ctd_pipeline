"""suncat.suna_cal: a cruise's SUNA calibration is stated in a file, validated, and never defaulted."""
from __future__ import annotations

import json

import pytest

from suncat import suna_cal as sc


def test_missing_file_is_an_error_that_says_how_to_fix_it(tmp_path):
    with pytest.raises(FileNotFoundError, match="set_suna_calibration.py"):
        sc.load_suna_calibration(tmp_path / "XX_01")


def test_fit_and_adopted_round_trip(tmp_path):
    c1, c2 = tmp_path / "XX_01", tmp_path / "XX_02"
    sc.write_suna_calibration(c1, "fit")
    r = sc.load_suna_calibration(c1)
    assert r["mode"] == "fit" and r["a1"] is None and r["a0"] is None
    sc.write_suna_calibration(c2, "adopted", a1=0.89, a0=-0.39, source="homogeneous-layer fit")
    r = sc.load_suna_calibration(c2)
    assert (r["mode"], r["a1"], r["a0"]) == ("adopted", 0.89, -0.39) and "homogeneous" in r["describe"]


@pytest.mark.parametrize("raw, msg", [
    ({"mode": "guess"}, "mode"),
    ({}, "mode"),
    ({"mode": "fit", "a1": 0.9}, "remove"),                                   # ambiguous: fit with numbers
    ({"mode": "adopted", "a1": 0.89, "source": "x"}, "a0"),
    ({"mode": "adopted", "a1": "abc", "a0": 0, "source": "x"}, "a1"),
    ({"mode": "adopted", "a1": 0.0, "a0": 0, "source": "x"}, "positive"),
    ({"mode": "adopted", "a1": float("nan"), "a0": 0, "source": "x"}, "finite"),
    ({"mode": "adopted", "a1": 0.89, "a0": -0.39}, "source"),                 # provenance is required
])
def test_ambiguous_or_incomplete_records_are_refused(tmp_path, raw, msg):
    p = sc.calibration_path(tmp_path / "XX_01"); p.parent.mkdir(parents=True)
    p.write_text(json.dumps(raw))
    with pytest.raises(sc.SunaCalibrationError, match=msg):
        sc.load_suna_calibration(tmp_path / "XX_01")


def test_existing_file_is_not_replaced_without_overwrite(tmp_path):
    c = tmp_path / "XX_01"
    sc.write_suna_calibration(c, "fit")
    with pytest.raises(FileExistsError):
        sc.write_suna_calibration(c, "adopted", a1=1.0, a0=0.0, source="x")
    sc.write_suna_calibration(c, "adopted", a1=1.0, a0=0.0, source="x", overwrite=True)
    assert sc.load_suna_calibration(c)["mode"] == "adopted"


def test_invalid_json_is_reported(tmp_path):
    p = sc.calibration_path(tmp_path / "XX_01"); p.parent.mkdir(parents=True)
    p.write_text("{mode: fit")
    with pytest.raises(sc.SunaCalibrationError, match="JSON"):
        sc.load_suna_calibration(tmp_path / "XX_01")
