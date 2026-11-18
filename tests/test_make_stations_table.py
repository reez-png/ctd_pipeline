"""tools/make_stations_table.py: step06's position table from the .cnv headers and the nutrient file. Synthetic only."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("make_stations_table", REPO / "tools" / "make_stations_table.py")
M = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(M)


def _cnv(path, start, nmea=None):
    head = [f"# start_time = {start} [NMEA time, header]"]
    if nmea:
        head = [f"* NMEA Latitude = {nmea[0]}", f"* NMEA Longitude = {nmea[1]}"] + head
    Path(path).write_text("\n".join(["* Sea-Bird SBE 19plus V2", *head, "# name 0 = prdM: p", "*END*", "1.0"]))


@pytest.mark.parametrize("text, want", [("05 29.37 N", 5 + 29.37 / 60), ("000 00.99 W", -0.99 / 60),
                                        ("5.48945", 5.48945), ("-0.016483", -0.016483), ("12 30 00 S", -12.5)])
def test_parse_coord(text, want):
    assert M.parse_coord(text) == pytest.approx(want)


def test_header_fix_first_then_nutrient_file_then_blank(tmp_path):
    _cnv(tmp_path / "ZZ_CTD_01.cnv", "Apr 16 2026 09:00:00", ("05 29.37 N", "000 01.00 W"))
    _cnv(tmp_path / "ZZ_CTD_02.cnv", "Apr 16 2026 11:00:00")
    _cnv(tmp_path / "ZZ_CTD_03.cnv", "Apr 17 2026 08:00:00")
    _cnv(tmp_path / "ZZ_CTD_04.cnv", "Apr 17 2026 10:00:00", ("00 00.00 N", "000 00.00 E"))   # GPS with no fix
    nut = pd.DataFrame(dict(cruise=["ZZ"] * 5 + ["YY"], station=["A", "A", "B", "B", "B", "B"],
                            latitude=[5.4900, 5.4900, 5.2000, 5.2001, 5.1999, 9.0],
                            longitude=[-0.0170, -0.0170, 0.1000, 0.1000, 0.1000, 9.0]))
    audit = pd.DataFrame(dict(station=["A", "B"], matched_cast=["ZZ_CTD_01", "ZZ_CTD_02"],
                              confidence=["high", "medium"], time_gap_min=[3.0, 59.2]))
    t = M.build(M.cast_headers(sorted(tmp_path.glob("*.cnv"))), M.nutrient_positions(nut, "ZZ", audit)).set_index("cast_id")
    assert t.loc["ZZ_CTD_01", "position_source"] == "cnv header NMEA fix"
    assert t.loc["ZZ_CTD_01", "lat_deg"] == pytest.approx(5 + 29.37 / 60)
    assert 0 < t.loc["ZZ_CTD_01", "header_vs_nutrient_km"] < 0.2                   # header and samples agree
    assert t.loc["ZZ_CTD_02", "lat_deg"] == pytest.approx(5.2) and t.loc["ZZ_CTD_02", "lon_deg"] == pytest.approx(0.1)
    assert "station B" in t.loc["ZZ_CTD_02", "position_source"] and "medium" in t.loc["ZZ_CTD_02", "position_source"]
    assert np.isnan(t.loc["ZZ_CTD_03", "lat_deg"]) and np.isnan(t.loc["ZZ_CTD_04", "lat_deg"])
    assert t.loc["ZZ_CTD_03", "date_time_utc_start"] == "2026-04-17 08:00:00"


def test_two_stations_on_one_cast_are_left_blank(tmp_path):
    _cnv(tmp_path / "ZZ_CTD_01.cnv", "Apr 16 2026 09:00:00")
    nut = pd.DataFrame(dict(cruise="ZZ", station=["A", "B"], latitude=[5.0, 5.3], longitude=[0.0, 0.1]))
    audit = pd.DataFrame(dict(station=["A", "B"], matched_cast=["ZZ_CTD_01"] * 2, confidence=["high"] * 2,
                              time_gap_min=[1.0, 2.0]))
    t = M.build(M.cast_headers(sorted(tmp_path.glob("*.cnv"))), M.nutrient_positions(nut, "ZZ", audit))
    assert len(t) == 1 and np.isnan(t.loc[0, "lat_deg"]) and t.loc[0, "station"] == "A+B"
    assert "by hand" in t.loc[0, "position_source"]
