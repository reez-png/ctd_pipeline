"""suncat.casts and the gridder's daylight test: no assumption about the cruise name, cast padding or longitude."""
from __future__ import annotations

import numpy as np
import pandas as pd

from suncat import casts as c
from suncat import gridding as g

IDS = ["ZZ12_03_CTD_001", "ZZ12_03_CTD_002", "ZZ12_03_CTD_010"]


def test_cast_number_reads_numbers_and_ids():
    assert c.cast_number(7) == 7 and c.cast_number("7") == 7 and c.cast_number("7.0") == 7
    assert c.cast_number("P45_06_CTD_07") == 7 and c.cast_number("ZZ12-03_CTD_107") == 107
    assert c.cast_number("") is None and c.cast_number(np.nan) is None and c.cast_number("nan") is None
    assert c.cast_number("7.5") is None


def test_match_cast_any_prefix_any_padding():
    assert c.match_cast("ZZ12_03_CTD_002", IDS) == "ZZ12_03_CTD_002"     # the id itself
    assert c.match_cast(2, IDS) == "ZZ12_03_CTD_002"                      # a number, three-digit files
    assert c.match_cast("2.0", IDS) == "ZZ12_03_CTD_002"
    assert c.match_cast("ZZ12-03_CTD_002", IDS) == "ZZ12_03_CTD_002"      # hyphenated cruise in a bottle log
    assert c.match_cast(10, IDS) == "ZZ12_03_CTD_010"
    assert c.match_cast(3, IDS) is None                                   # no such cast: not invented
    assert c.match_cast(1, IDS + ["OTHER_CTD_1"]) is None                 # ambiguous: refused


def test_positions_keep_the_equator_and_drop_only_zeroed_fixes(tmp_path):
    pd.DataFrame({"cast": [1, 2, 10], "latitude_start": [0.6, 0.0, -0.3], "longitude_start": [-3.0, -3.1, -3.2],
                  "date_time_utc_start": ["2026-03-10 06:00:00", "2026-03-11 09:00:00", ""]}
                 ).to_csv(tmp_path / "ZZ12_03_CTD_cast_summary.csv", index=False)
    pos = c.read_cast_positions(tmp_path, "ZZ12_03", IDS)
    assert list(pos.index) == IDS
    assert pos.at["ZZ12_03_CTD_001", "lat"] == 0.6                      # 0.6 N is a real station
    assert np.isnan(pos.at["ZZ12_03_CTD_002", "lat"])                    # an exact 0 is a dropped fix
    assert pos.at["ZZ12_03_CTD_010", "lat"] == -0.3
    assert pd.isna(pos.at["ZZ12_03_CTD_010", "start_utc"])
    assert c.read_cast_positions(tmp_path / "none", "ZZ12_03", IDS).empty


def test_solar_hour_follows_longitude():
    t = pd.Timestamp("2026-03-10 06:00", tz="UTC")
    assert c.solar_hour(t, 0.0) == 6.0
    assert abs(c.solar_hour(t, 45.0) - 9.0) < 1e-9                       # 45 E: 3 h ahead of UTC
    assert abs(c.solar_hour(t, -90.0) - 0.0) < 1e-9
    assert c.solar_hour(t, np.nan) is None and c.solar_hour(pd.NaT, 10.0) is None


def test_npq_hour_prefers_the_sun_and_says_when_it_falls_back(tmp_path):
    pos = pd.DataFrame({"lat": [-12.4], "lon": [40.6], "start_utc": [pd.Timestamp("2026-03-10 06:00", tz="UTC")],
                        "source": "x"}, index=pd.Index(["A_CTD_001"], name="cast_id"))
    h, rule = c.npq_hour("A_CTD_001", pos, fallback_time=pd.Timestamp("2026-03-10 23:00"))
    assert rule == "solar" and abs(h - (6 + 40.6 / 15)) < 1e-6
    h, rule = c.npq_hour("A_CTD_002", pos, fallback_time=pd.Timestamp("2026-03-10 13:30"))
    assert rule == "clock" and h == 13.5
    assert c.npq_hour("A_CTD_002", pos) == (None, "none")


def test_gridder_daylight_uses_solar_time_when_the_longitude_is_known():
    cfg = g.GridConfig()
    t = pd.Timestamp("2026-03-10 06:00", tz="UTC")                       # 09:00 local at 45 E, before 07:00 at 0 E
    p = np.array([5.0, 20.0])
    assert list(g.npq_risk(p, t, cfg, lon=45.0)) == [True, False]
    assert list(g.npq_risk(p, t, cfg, lon=None)) == [False, False]       # no longitude: clock time, offset 0
    assert list(g.npq_risk(p, t, g.GridConfig(utc_offset_hours=3.0), lon=None)) == [True, False]
    assert list(g.npq_risk(p, t, g.GridConfig(solar_time=False, utc_offset_hours=0.0), lon=45.0)) == [False, False]
