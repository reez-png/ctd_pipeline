"""suncat.nitrate.bottle_match: the one bottle-pairing rule step02b and the lag diagnosis share."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from suncat.nitrate import bottle_match as BM


def test_bottles_flag_replicate_qc_and_cast_numbers(tmp_path):
    f = tmp_path / "b.csv"
    pd.DataFrame({"cast_id": ["X_CTD_01", "X_CTD_01", "X_CTD_01", "X_CTD_01", "2", "X_CTD_03"],
                  "depth_m": [10, 30, 30, 50, 20, 5],
                  "nitrate_uM": [1.0, 8.0, 8.9, 12.0, 3.0, np.nan],
                  "replicate_id": [np.nan, "a", "a", np.nan, np.nan, np.nan],
                  "flag": ["", "", "", "BDL", "", ""]}).to_csv(f, index=False)
    b = BM.load_bottles(f)
    # the 30 m replicates differ by 0.9 > 0.5 umol/L: both dropped; the BDL row and the empty value too
    assert list(zip(b["cast"], b["depth_m"])) == [("X_CTD_01", 10.0), ("2", 20.0)]
    assert BM.load_bottles(tmp_path / "missing.csv").empty


def test_moving_and_the_down_cast_cell_median():
    t = np.arange(10.0)
    z = np.array([0, 0, 0, 1, 2, 3, 4, 4, 4, 5], float)        # a soak, 1 m/s, a stop, on again
    m = BM.moving(t, z)
    # centred difference: the soak and the middle of the stop are still, their edges count as moving
    assert m.tolist() == [False, False, True, True, True, True, True, False, True, True]
    y = z.copy()
    y[7] = 99.0                                                  # a value read while stopped
    v_all, n_all = BM.down_value(z, y, np.ones(10, bool), 2.0)
    v, n = BM.down_value(z, y, m, 2.0)
    # cells centred on 2 m within +-2 m hold 0..4 m; the stopped frame lifts the 4 m cell only when it counts
    assert n == 5 and v == pytest.approx(2.0) and n_all == 5 and v_all == pytest.approx(2.0)
    v, n = BM.down_value(z, y, m, 4.0)
    # cells 2, 3, 4, 5 m: moving frames give means 2, 3, 4, 5 (median 3.5); with the stopped frame the 4 m cell is 35.7
    assert v == pytest.approx(3.5) and BM.down_value(z, y, np.ones(10, bool), 4.0)[0] == pytest.approx(4.0)
    assert np.isnan(BM.down_value(z, y, m, 20.0)[0])


def test_pair_gives_every_bottle_a_value_or_a_reason():
    grid = np.arange(0.0, 51.0)
    frames = {"C_01": dict(z=grid, y=grid / 10, down=np.ones_like(grid, bool)),
              "C_02": dict(z=grid, y=grid / 10, down=np.zeros_like(grid, bool))}
    prof = {"C_01": {"grid": grid, "A": grid / 10}, "C_02": {"grid": grid, "A": grid / 10}}
    b = pd.DataFrame({"cast": ["1", "C_01", "C_02", "C_09", "C_01"], "depth_m": [20.0, 70.0, 20.0, 5.0, 40.0],
                      "bottle": [2.0, 7.0, 2.0, 0.5, 4.0]})
    T = BM.pair(b, frames, prof, ["A"])
    assert T.loc[0, "cast"] == "C_01" and T.loc[0, "reason"] == "" and T.loc[0, "A"] == pytest.approx(2.0)
    assert T.loc[1, "reason"] == "outside the profile's depth range"
    assert T.loc[2, "reason"].startswith("no moving down-cast frame")
    assert T.loc[3, "reason"] == "no frames on this cast"
    c = BM.compare(T, ["A"])
    assert not c["fits"] and c["n"] == 2 and c["n_bottles"] == 5 and sum(c["left_out"].values()) == 3


def test_compare_scores_and_paired_test():
    rng = np.random.default_rng(0)
    b = np.linspace(0, 20, 40)
    T = pd.DataFrame({"cast": "C", "depth_m": b, "bottle": b, "down": 0.9 * b + rng.normal(0, 2, 40),
                      "good": 0.9 * b + rng.normal(0, 0.2, 40), "reason": ""})
    c = BM.compare(T, ["good"])
    f = c["fits"]["good"]
    assert abs(f["gain"] - 0.9) < 0.02 and f["R2"] > 0.99 and f["loocv"] < c["fits"]["down"]["loocv"]
    assert f["better"] > 30 and f["p_sign"] < 0.001
