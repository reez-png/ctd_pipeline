"""build_bottle_table.py (P45_06 style: nutrient samples paired to casts by time). A medium-confidence pairing is
kept with a pairing_note; flag stays empty, because step03 and step02b leave out any row with a flag."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import build_bottle_table as BT  # noqa: E402


def test_medium_pairing_is_kept_with_a_note_not_a_flag(tmp_path, monkeypatch):
    ctd, l2, bd = tmp_path / "ctd", tmp_path / "l2", tmp_path / "bottles"
    ctd.mkdir(); l2.mkdir()
    for n, start in (("03", "Apr 16 2026 11:16:12"), ("06", "Apr 16 2026 14:33:22"), ("07", "Apr 17 2026 12:04:27")):
        (ctd / f"ZZ_CTD_{n}.cnv").write_text(f"# start_time = {start} [Instrument's time stamp, header]\n*END*\n")
    for n in ("03", "06"):
        (l2 / f"ZZ_CTD_{n}_SUNA_1s.csv").write_text("x\n")
    nut = tmp_path / "nut.csv"
    pd.DataFrame({"cruise": ["Z-1"] * 4, "station": ["P1", "P1", "P3", "J9"],
                  "date_time": ["2026-04-16 12:15", "2026-04-16 12:15", "2026-04-16 14:40", "2026-04-17 12:10"],
                  "depth": [5, 30, 5, 5], "nitrate_nitrite": [0.1, 0.2, 0.3, 0.4], "replicate": [1, 1, 1, 1]}
                 ).to_csv(nut, index=False)
    for k, v in dict(CRUISE_ID="ZZ", NUTRIENT_FILE=nut, NUTRIENT_CRUISE_LABEL="Z-1", L1_CTD_ROOT=ctd,
                     L2_SUNA_ROOT=l2, BOTTLE_DIR=bd).items():
        monkeypatch.setattr(BT, k, v)
    BT.main()
    b = pd.read_csv(bd / "ZZ_bottle_nitrate.csv", keep_default_na=False)
    assert set(b["station"]) == {"P1", "P3"}                      # J9's cast has no SUNA: excluded
    p1 = b[b["station"] == "P1"]
    assert (p1["cast_id"] == "ZZ_CTD_03").all() and (p1["flag"] == "").all()
    assert p1["pairing_note"].str.startswith("REVIEW_PAIRING_59min").all()
    assert (b.loc[b["station"] == "P3", "pairing_note"] == "").all()   # 7 min: high confidence, no note


def test_a_hand_flag_survives_a_rebuild(tmp_path, monkeypatch, capsys):
    test_medium_pairing_is_kept_with_a_note_not_a_flag(tmp_path, monkeypatch)
    path = tmp_path / "bottles" / "ZZ_bottle_nitrate.csv"
    b = pd.read_csv(path, keep_default_na=False)
    b.loc[(b["station"] == "P1") & (b["depth_m"] == 30), "flag"] = "suspect: duplicates another depth"
    b.to_csv(path, index=False)
    BT.main()
    b = pd.read_csv(path, keep_default_na=False)
    assert b.loc[(b["station"] == "P1") & (b["depth_m"] == 30), "flag"].tolist() == ["suspect: duplicates another depth"]
    assert (b.loc[b["depth_m"] != 30, "flag"] == "").all()
    assert "hand flag kept: P1 30" in capsys.readouterr().out
