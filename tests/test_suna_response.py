"""suncat.suna_response and tools/set_suna_slow_casts.py: per-cast SUNA response flags stated once per cruise."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from suncat import suna_response as R

REPO = Path(__file__).resolve().parents[1]


def _tool():
    spec = importlib.util.spec_from_file_location("set_suna_slow_casts", REPO / "tools" / "set_suna_slow_casts.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


TAB = pd.DataFrame({"cast": ["C_01", "C_02", "C_03", "C_04", "C_05"],
                    "tau": [44.4, 6.2, 12.0, np.nan, np.nan],
                    "edge": [False, False, False, True, False]})


def test_rule_flags_slow_and_edge_casts_only():
    f = R.flag_slow(TAB, 5.9).set_index("cast_id")
    # 3 x 5.9 = 17.7 s: 44.4 is slow, 12.0 is not; a tau at the grid's edge is slower than measurable
    assert f.loc["C_01", "slow_response"] and not f.loc["C_03", "slow_response"]
    assert f.loc["C_04", "slow_response"] and not f.loc["C_05", "slow_response"]
    assert R.flag_slow(TAB, 5.9, factor=2.0).set_index("cast_id").loc["C_03", "slow_response"]
    with pytest.raises(ValueError):
        R.flag_slow(TAB, float("nan"))


def test_write_load_describe_and_no_silent_overwrite(tmp_path):
    assert R.load_response(tmp_path) == {} and "not assessed" in R.describe(tmp_path)
    with pytest.raises(ValueError):
        R.write_response(tmp_path, R.flag_slow(TAB, 5.9), "")
    R.write_response(tmp_path, R.flag_slow(TAB, 5.9), "synthetic rule")
    r = R.load_response(tmp_path)
    assert r["C_01"] == (44.4, True) and r["C_02"] == (6.2, False) and np.isnan(r["C_05"][0])
    assert "2 slow: 01, 04" in R.describe(tmp_path) and "synthetic rule" in R.describe(tmp_path)
    with pytest.raises(FileExistsError):
        R.write_response(tmp_path, R.flag_slow(TAB, 5.9), "again")


def test_malformed_file_is_refused(tmp_path):
    (tmp_path / "metadata").mkdir()
    R.response_path(tmp_path).write_text("cast_id,tau_s,slow_response\nC_01,44,maybe\n")
    with pytest.raises(ValueError, match="True or False"):
        R.load_response(tmp_path)


def test_tool_writes_from_step02b_and_hand_lists_need_a_source(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    a = tmp_path / "cruises" / "ZZ" / "_audit" / "suna"
    a.mkdir(parents=True)
    TAB.to_csv(a / "ZZ__suna_response_time.csv", index=False)
    (a / "ZZ__suna_response_summary.json").write_text(json.dumps({"tau_median": 5.9}))
    T = _tool()
    assert T.main(["ZZ"]) == 0
    out = capsys.readouterr().out
    assert "threshold 17.7 s" in out and "slowest unflagged 12.0 s, fastest flagged 44.4 s" in out
    C = tmp_path / "cruises" / "ZZ"
    assert {c for c, (_, s) in R.load_response(C).items() if s} == {"C_01", "C_04"}
    assert T.main(["ZZ", "--casts", "C_02"]) == 1                       # no source
    assert T.main(["ZZ", "--casts", "C_02", "--source", "by eye", "--force"]) == 0
    assert {c for c, (_, s) in R.load_response(C).items() if s} == {"C_02"}
