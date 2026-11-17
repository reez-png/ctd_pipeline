"""suncat.firing: discrete samples at the depth their Niskin was fired, not the sheet's nominal depth."""
from __future__ import annotations

import numpy as np
import pandas as pd

from suncat import firing as F


def _record(tmp_path):
    f = F.firings_path(tmp_path, "ZZ")
    f.parent.mkdir(parents=True)
    pd.DataFrame({"cast_id": ["ZZ_CTD_02"] * 3, "bottle_no": [1, 2, 4], "pressure_dbar": [28.9, 27.6, 14.7],
                  "depth_m": [28.74, 27.48, 14.6], "time_utc": ["2026-07-29T11:48:21"] * 3}).to_csv(f, index=False)


def test_fired_depth_replaces_the_nominal_one(tmp_path):
    _record(tmp_path)
    s = pd.DataFrame({"cast_id": ["ZZ_CTD_02", "ZZ_CTD_02", "ZZ_CTD_02", "ZZ_CTD_02"],
                      "niskin": [2, 2, 4, 6], "depth_m": [30.0, 30.0, 15.0, 5.0], "value": [0.55, 0.55, 1.12, 0.2]})
    out, rep = F.fired_depths(s, F.load_firings(tmp_path, "ZZ"))
    assert out["depth_m"].tolist() == [27.48, 27.48, 14.6, 5.0]                 # replicates share their Niskin
    assert out["depth_source"].tolist() == ["fired", "fired", "fired", "nominal"]
    assert out["depth_nominal_m"].tolist() == [30.0, 30.0, 15.0, 5.0]
    assert rep["n_fired"] == 3 and rep["n_nominal"] == 1 and abs(rep["d_max"] - 2.52) < 1e-9
    line = F.describe(rep, "chl.csv")
    assert "3 of 4 rows at their fired depth" in line and "1 left at the nominal depth" in line


def test_no_record_or_no_niskin_keeps_the_sheet_and_says_so(tmp_path):
    s = pd.DataFrame({"cast_id": ["ZZ_CTD_02"], "niskin": [2], "depth_m": [30.0]})
    out, rep = F.fired_depths(s, F.load_firings(tmp_path, "ZZ"))
    assert out["depth_m"].tolist() == [30.0] and "no firing record" in F.describe(rep)
    _record(tmp_path)
    out, rep = F.fired_depths(s.drop(columns="niskin"), F.load_firings(tmp_path, "ZZ"))
    assert out["depth_m"].tolist() == [30.0] and "no Niskin column" in F.describe(rep)


def test_a_cast_number_finds_its_cast(tmp_path):
    _record(tmp_path)
    out, rep = F.fired_depths(pd.DataFrame({"cast_id": ["2"], "niskin": [1], "depth_m": [31.0]}),
                              F.load_firings(tmp_path, "ZZ"))
    assert out["depth_m"].tolist() == [28.74] and np.isclose(rep["d_median"], 2.26)
