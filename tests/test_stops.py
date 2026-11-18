"""suncat.stops: bottle depths from where the rosette was held on the up-cast. Synthetic 4 Hz casts with ship heave."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from suncat import stops as S


def trace(plan, hz=4.0, heave=0.4, seed=0):
    """plan: [(target dbar, hold s)] visited in order at 0.8 dbar/s; heave: amplitude of an 8 s swell (dbar)."""
    rng = np.random.default_rng(seed)
    p, cur = [], 0.0
    for target, hold in plan:
        n = int(abs(target - cur) / 0.8 * hz)
        p += list(np.linspace(cur, target, max(n, 1), endpoint=False))
        p += [target] * int(hold * hz)
        cur = target
    p = np.asarray(p, float)
    t = np.arange(p.size) / hz
    p = np.clip(p + heave * np.sin(2 * np.pi * t / 8.0) * (p > 1) + rng.normal(0, 0.02, p.size), 0, None)
    temp = 29.0 - 0.13 * p
    return t, p, temp


PLAN = [(0.0, 30), (10.0, 60), (0.5, 10), (85.0, 40), (60.0, 30), (45.0, 40), (25.0, 40), (5.0, 40), (0.0, 20)]


def test_stops_are_found_and_labelled_through_heave():
    t, p, temp = trace(PLAN)
    st = S.find_stops(t, p, temp, lat=5.4)
    up = st[st["leg"] != "down"]
    assert list(up["leg"]) == ["bottom", "up", "up", "up", "up"]          # 5 m is a hold, not the surface wait
    assert np.allclose(up["p_dbar"], [85, 60, 45, 25, 5], atol=0.3)
    assert (st.loc[st["leg"] == "down", "p_dbar"].round() == 10).any()          # the soak, on the way down
    assert np.allclose(up["temp_c"], 29.0 - 0.13 * up["p_dbar"], atol=0.1)
    assert up["depth_m"].iloc[0] == pytest.approx(84.5, abs=0.3)                  # UNESCO: 85 dbar is about 84.5 m


def test_sheet_depths_take_the_nearest_upcast_stop_or_keep_the_sheet_depth():
    t, p, temp = trace(PLAN)
    st = S.find_stops(t, p, temp, lat=5.4)
    m = S.match_depths([85, 47, 25, 5, 15, 10], st).set_index("nominal_m")
    assert m.loc[47, "stop_depth_m"] == pytest.approx(44.7, abs=0.4)             # sheet 47 m, held at 45 dbar
    assert m.loc[5, "stop"] >= 0 and m.loc[85, "stop"] >= 0
    assert m.loc[15, "stop"] == -1 and "away" in m.loc[15, "note"]                # 10 m from the nearest stop
    assert m.loc[10, "stop"] == -1                                                # the soak is a down-cast stop


def test_two_sheet_depths_on_one_stop_both_keep_the_sheet_depth():
    t, p, temp = trace(PLAN)
    m = S.match_depths([44, 47], S.find_stops(t, p, temp))
    assert (m["stop"] == -1).all() and m["note"].str.contains("shares").all()


def test_a_record_without_an_upcast_gives_no_stops(tmp_path):
    t, p, temp = trace([(0.0, 30), (10.0, 60), (85.0, 5)])
    f = tmp_path / "ZZ_CTD_01.cnv"
    _write_cnv(f, t, p, temp)
    st, why = S.cast_stops([f, tmp_path / "missing.cnv"])
    assert st.empty and "no up-cast" in why


def _write_cnv(path, t, p, temp, start="Apr 16 2026 14:33:22"):
    cols = [("timeS", t), ("prdM", p), ("t090C", temp), ("c0S/m", np.where(p > 0.2, 5.5, 0.001))]
    head = [f"# name {i} = {c}: {c}" for i, (c, _) in enumerate(cols)]
    body = [" ".join(f"{v:.4f}" for v in row) for row in zip(*[v for _, v in cols])]
    Path(path).write_text("\n".join([f"# start_time = {start} [NMEA time, header]", *head, "*END*", *body]))


def test_cast_stops_reads_a_cnv(tmp_path):
    t, p, temp = trace(PLAN)
    f = tmp_path / "ZZ_CTD_06.cnv"
    _write_cnv(f, t, p, temp)
    st, used = S.cast_stops([f])
    assert used == "ZZ_CTD_06.cnv" and (st["leg"] == "up").sum() == 4


def test_the_surface_wait_takes_only_a_surface_sheet_depth():
    """P45_06: the up-casts held at the bottom and at about 1 m before recovery, never at 5 m; the 5 m sheet depth must
    not take the 1 m wait (it did in the first version), but a 1 m sheet depth may."""
    t, p, temp = trace([(0.0, 30), (10.0, 60), (0.5, 10), (60.0, 30), (1.2, 40), (0.0, 20)])
    st = S.find_stops(t, p, temp)
    assert list(st.loc[st["leg"] != "down", "leg"]) == ["bottom", "surface"]
    m = S.match_depths([5, 1, 25], st).set_index("nominal_m")
    assert m.loc[5, "stop"] == -1 and m.loc[25, "stop"] == -1
    assert m.loc[1, "stop_depth_m"] == pytest.approx(1.2, abs=0.3)
