"""suncat.gridding: L2 products onto the SBE 1 m down-cast pressure bins (synthetic data only)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from suncat import gridding as g

BIN_CNV = (
    "* synthetic Bin Average output\n"
    "# name 0 = prdM: Pressure, Strain Gauge [db]\n"
    "# name 1 = tv290C: Temperature [ITS-90, deg C]\n"
    "# name 2 = depSM: Depth [salt water, m], lat = 5\n"
    "# name 3 = nbin: number of scans per bin\n"
    "# name 4 = flag: flag\n"
    "# name 5 = sigma-\u00e900: Density [sigma-theta, kg/m^3]\n"
    "# start_time = Apr 16 2026 11:16:12 [Instrument's time stamp, header]\n"
    "# binavg_bintype = decibars\n"
    "# binavg_binsize = 1\n"
    "*END*\n"
    "      8.000    28.90      7.956          5 0.0000e+00   22.30\n"
    "      9.000    28.80      8.950          6 0.0000e+00   22.32\n"
    "     10.000    28.70      9.945          6 0.0000e+00   22.34\n"
    "     11.000    28.60     10.939          5 0.0000e+00   22.36\n"
)


@pytest.fixture
def ctd(tmp_path):
    f = tmp_path / "C1_1m_down.cnv"
    f.write_text(BIN_CNV, encoding="latin-1")
    return g.read_sbe_bins(f)


def frames(rows, gap=0.3):
    f = pd.DataFrame(rows, columns=["pressure_db", "leg", "L1_clean", "no3_qc", "lag_limited"])
    f.insert(2, "ctd_gap_s", gap)
    return f


def test_reader_gets_sbe_bin_metadata(ctd):
    df, meta = ctd
    assert meta["bintype"] == "decibars" and meta["binsize"] == 1.0
    assert meta["start_time"] == pd.Timestamp("2026-04-16 11:16:12", tz="UTC")
    assert list(df["prdM"]) == [8.0, 9.0, 10.0, 11.0]


def test_rows_are_the_sbe_bins_and_ctd_columns_are_untouched(ctd):
    df, meta = ctd
    out, st = g.grid_cast("C1", df, meta)
    assert list(out["pressure_db"]) == [8.0, 9.0, 10.0, 11.0]
    assert list(out["tv290C"]) == [28.9, 28.8, 28.7, 28.6]
    assert list(out["ctd_n_scans"]) == [5, 6, 6, 5]
    assert "flag" not in out and "prdM" not in out and list(out["ctd_pressure_db"]) == [8.0, 9.0, 10.0, 11.0]
    assert list(out["sigma-theta00"]) == [22.30, 22.32, 22.34, 22.36]
    assert st.ctd_bins == 4 and (out["nitrate_n"] == 0).all()


def test_nitrate_mean_n_sd_and_every_filter(ctd):
    df, meta = ctd
    nf = frames([
        (8.40, "down", True, 1.0, True),     # bin 8
        (8.49, "down", True, 3.0, True),     # bin 8 -> mean 2, n 2, sd sqrt(2)
        (8.60, "down", True, 5.0, True),     # bin 9 (8.6 rounds to 9)
        (9.20, "down", False, 99.0, True),   # failed spectral QC
        (9.30, "up", True, 77.0, True),      # upcast
        (10.1, "down", True, np.nan, True),  # no calibrated value
        (30.0, "down", True, 9.0, True),     # below the CTD's deepest bin
    ])
    out, st = g.grid_cast("C1", df, meta, nf)
    r = out.set_index("pressure_db")
    assert r.loc[8.0, "nitrate_uM"] == 2.0 and r.loc[8.0, "nitrate_n"] == 2
    assert r.loc[8.0, "nitrate_sd"] == pytest.approx(np.sqrt(2.0))
    assert r.loc[9.0, "nitrate_uM"] == 5.0 and r.loc[9.0, "nitrate_n"] == 1 and np.isnan(r.loc[9.0, "nitrate_sd"])
    assert np.isnan(r.loc[10.0, "nitrate_uM"]) and r.loc[10.0, "nitrate_n"] == 0
    assert st.nitrate_used == 3
    assert st.nitrate_dropped == {"upcast": 1, "ctd_not_kept": 0, "not_clean": 1, "no_value": 1, "outside_ctd_bins": 1}
    assert st.nitrate_n_max == 2
    assert r.loc[8.0, "nitrate_lag_limited"] == True and pd.isna(r.loc[10.0, "nitrate_lag_limited"])  # noqa: E712


def test_bin_edges(ctd):
    df, meta = ctd
    # bin b spans [b - 0.5, b + 0.5): an exact edge goes up (7.5 -> bin 8, 8.5 -> bin 9, 10.4999 -> bin 10)
    nf = frames([(7.5, "down", True, 1.0, False), (8.4999, "down", True, 3.0, False),
                 (8.5, "down", True, 5.0, False), (10.4999, "down", True, 7.0, False)])
    out, _ = g.grid_cast("C1", df, meta, nf)
    r = out.set_index("pressure_db")
    assert r.loc[8.0, "nitrate_n"] == 2 and r.loc[9.0, "nitrate_n"] == 1 and r.loc[10.0, "nitrate_n"] == 1


def test_all_frames_unclean_gives_no_nitrate_and_a_note(ctd):
    df, meta = ctd
    out, st = g.grid_cast("C1", df, meta, frames([(9.0, "down", False, 4.0, True)] * 5))
    assert out["nitrate_uM"].isna().all() and st.nitrate_used == 0
    assert any("no nitrate frame passed" in n for n in st.notes)


def test_frames_where_the_ctd_kept_no_data_are_left_out(ctd):
    # a 2-minute hold at 10 dbar: Loop Edit flagged the CTD scans, so the 1 s record has no seconds there
    # and the frames borrowed values from up to 60 s away; only the frame on the CTD record is used
    df, meta = ctd
    hold = frames([(10.0, "down", True, 3.0 + 0.01 * i, True) for i in range(120)], gap=np.r_[0.4, np.linspace(1.2, 60, 119)])
    out, st = g.grid_cast("C1", df, meta, hold)
    r = out.set_index("pressure_db")
    assert r.loc[10.0, "nitrate_n"] == 1 and r.loc[10.0, "nitrate_uM"] == 3.0
    assert st.nitrate_dropped["ctd_not_kept"] == 119


def test_scan_level_rule_follows_loop_edit_and_bins_by_the_scan_pressure(ctd):
    # P45_06 CTD_04's bottom stop: the 1 s record has every second (gap small), but Loop Edit kept no scan
    # during the heave; only frames whose concurrent 4 Hz scan was kept enter the bin, at that scan's pressure
    df, meta = ctd
    nf = frames([(10.4, "down", True, 1.0, True)] * 6, gap=0.2)
    nf["ctd_scan_pressure_db"] = [10.40, 10.45, 10.47, 10.46, 10.52, 9.40]
    nf["ctd_scan_kept"] = [True, False, False, False, True, True]
    nf["no3_qc"] = [1.0, 50.0, 50.0, 50.0, 3.0, 7.0]
    out, st = g.grid_cast("C1", df, meta, nf)
    r = out.set_index("pressure_db")
    assert r.loc[10.0, "nitrate_n"] == 1 and r.loc[10.0, "nitrate_uM"] == 1.0     # 10.40 kept
    assert r.loc[11.0, "nitrate_n"] == 1 and r.loc[11.0, "nitrate_uM"] == 3.0     # binned by scan p 10.52
    assert r.loc[9.0, "nitrate_uM"] == 7.0 and st.nitrate_dropped["ctd_not_kept"] == 3
    assert not any("gap rule" in n for n in st.notes)


def test_without_scan_status_the_1s_gap_rule_is_used_and_noted(ctd):
    df, meta = ctd
    nf = frames([(9.0, "down", True, 2.0, False)], gap=0.3)
    nf["ctd_scan_pressure_db"], nf["ctd_scan_kept"] = np.nan, pd.NA
    out, st = g.grid_cast("C1", df, meta, nf)
    assert out.set_index("pressure_db").loc[9.0, "nitrate_n"] == 1
    assert any("gap rule" in n for n in st.notes)


def test_mean_pressure_bins_are_keyed_by_centre_and_ambiguous_files_refused(ctd):
    df, meta = ctd
    shifted = df.assign(prdM=[8.247, 9.023, 10.002, 10.957])           # P45_07: mean pressure per bin
    out, st = g.grid_cast("C1", shifted, meta, frames([(9.1, "down", True, 2.0, False)]))
    assert list(out["pressure_db"]) == [8.0, 9.0, 10.0, 11.0]
    assert list(out["ctd_pressure_db"]) == [8.247, 9.023, 10.002, 10.957]
    assert st.ctd_bin_offset_max == pytest.approx(0.247)
    with pytest.raises(ValueError, match="one-to-one"):
        g.grid_cast("C1", df.assign(prdM=[8.0, 8.3, 10.0, 11.0]), meta)   # two rows in bin 8


def test_old_nitrate_product_is_refused(ctd):
    df, meta = ctd
    old = pd.DataFrame({"pressure_db": [9.0], "leg": ["down"], "L1_clean": [True], "no3_qc": [1.0]})
    with pytest.raises(ValueError, match="re-run step03"):
        g.grid_cast("C1", df, meta, old)


def test_chlorophyll_binned_from_its_own_rows_skipping_blanks_and_other_casts(ctd):
    df, meta = ctd
    chl = pd.DataFrame({"pressure_db": [8.0, 9.0, 9.0, 10.0, 8.0],
                        "sensor_blank": [False, False, True, False, False],
                        "cast": ["C1", "C1", "C1", "C1", "C2"],
                        "chl_corrected": [0.5, 0.7, 9.9, np.nan, 5.0]})
    out, st = g.grid_cast("C1", df, meta, chl=chl)
    r = out.set_index("pressure_db")
    assert r.loc[8.0, "chl_mg_m3"] == 0.5 and r.loc[9.0, "chl_mg_m3"] == 0.7
    assert r.loc[9.0, "chl_n"] == 1 and r.loc[10.0, "chl_n"] == 0 and st.chl_used == 2


def test_npq_risk_follows_step07_rule(ctd):
    df, meta = ctd
    out, _ = g.grid_cast("C1", df, meta)                         # 11:16 UTC = local daylight
    assert list(out["chl_npq_risk"]) == [True, True, True, False]  # <= 10 dbar
    night = dict(meta, start_time=pd.Timestamp("2026-04-16 22:00", tz="UTC"))
    out2, _ = g.grid_cast("C1", df, night)
    assert not out2["chl_npq_risk"].any()


def test_depth_binned_file_is_refused(ctd):
    df, meta = ctd
    with pytest.raises(ValueError, match="not decibars"):
        g.grid_cast("C1", df, dict(meta, bintype="meters"))


def test_grid_cruise_writes_one_table_per_cast(tmp_path, ctd):
    ctd_dir = tmp_path / "ctd"
    ctd_dir.mkdir()
    (ctd_dir / "C1_1m_down.cnv").write_text(BIN_CNV, encoding="latin-1")
    nit = tmp_path / "nit"
    nit.mkdir()
    frames([(9.0, "down", True, 4.0, True)]).to_csv(nit / "C1_nitrate_QC.csv", index=False)
    tables, stats = g.grid_cruise("X", ctd_dir, nit, None, None, tmp_path / "out")
    t = pd.read_csv(tmp_path / "out" / "C1_L3_1m_down.csv")
    assert list(t.columns[:3]) == ["cruise_id", "cast_id", "pressure_db"]
    assert t.set_index("pressure_db").loc[9.0, "nitrate_uM"] == 4.0
    assert stats.loc[0, "cast_id"] == "C1" and stats.loc[0, "nitrate_used"] == 1
