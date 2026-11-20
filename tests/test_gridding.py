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


def test_slow_response_travels_with_the_nitrate(ctd):
    df, meta = ctd
    f = frames([(8.0, "down", True, 1.0, True), (9.0, "down", True, 2.0, True)])
    out, st = g.grid_cast("C1", df, meta, f, response=(44.4, True))
    r = out.set_index("pressure_db")
    assert r.loc[8.0, "nitrate_slow_response"] == True and pd.isna(r.loc[10.0, "nitrate_slow_response"])  # noqa: E712
    assert (out["suna_tau_s"] == 44.4).all() and any("SLOW" in n for n in st.notes)
    out, _ = g.grid_cast("C1", df, meta, f, response=(6.2, False))
    assert out.set_index("pressure_db").loc[8.0, "nitrate_slow_response"] == False  # noqa: E712
    # no response file for the cruise: not assessed, so empty, not False
    out, _ = g.grid_cast("C1", df, meta, f)
    assert out["nitrate_slow_response"].isna().all() and out["suna_tau_s"].isna().all()


def test_updown_difference_needs_both_legs():
    d = g.updown_difference([10, 20], p_down=[9.5, 10, 10.5, 20], v_down=[1, 1, 1, 5], p_up=[9.6, 10.2, 10.4, 20, 20.1],
                            v_up=[4, 4, 4, 5, 5], win=1.0, min_per_leg=3)
    assert d[0] == pytest.approx(3.0) and np.isnan(d[1])          # at 20 dbar a leg has fewer than 3 values


def test_bins_where_down_and_up_disagree_are_flagged(ctd):
    """A lag, or water held in the optical path: the up-cast reads high below 9.5 dbar. The down-cast value stays
    the grid's value; the bin is flagged where |up - down| > 2 uM, left empty where a leg is too thin."""
    df, meta = ctd
    rows = []
    for z in (8.0, 9.0, 10.0, 11.0):
        rows += [(z - 0.2, "down", True, 1.0, True), (z, "down", True, 1.0, True), (z + 0.2, "down", True, 1.0, True)]
    for z, v, n in ((8.0, 1.3, 3), (9.0, 1.3, 3), (10.0, 4.5, 3), (11.0, 4.5, 2)):
        rows += [(z + 0.1 * k, "up", True, v, True) for k in range(n)]
    out, st = g.grid_cast("C1", df, meta, frames(rows), cfg=g.GridConfig(hyst_win_dbar=0.45))
    r = out.set_index("pressure_db")
    assert r.loc[8.0, "nitrate_updown_uM"] == pytest.approx(0.3) and r.loc[8.0, "nitrate_hyst_flag"] == False  # noqa: E712
    assert r.loc[10.0, "nitrate_updown_uM"] == pytest.approx(3.5) and r.loc[10.0, "nitrate_hyst_flag"] == True  # noqa: E712
    assert pd.isna(r.loc[11.0, "nitrate_hyst_flag"]) and np.isnan(r.loc[11.0, "nitrate_updown_uM"])
    assert r.loc[10.0, "nitrate_uM"] == 1.0                                     # the value is still the down-cast's
    assert st.nitrate_hyst_assessed == 3 and st.nitrate_hyst_flagged == 1
    assert any("differ by more than 2 uM in 1 of 3 assessed bins (10 to 10 dbar" in n for n in st.notes)
    row = st.row()
    assert row["nitrate_hyst_flagged"] == 1 and row["nitrate_hyst_assessed"] == 3


def test_no_upcast_frames_leave_the_flag_empty(ctd):
    df, meta = ctd
    out, st = g.grid_cast("C1", df, meta, frames([(9.0, "down", True, 2.0, False)] * 4))
    assert out["nitrate_hyst_flag"].isna().all() and st.nitrate_hyst_assessed == 0


def test_a_gradient_with_frames_at_different_depths_is_not_a_disagreement():
    """No lag, nitrate rising 1 uM/dbar; the up-cast's frames crowd at a bottle stop 2 dbar below the centre.
    Window medians would read 2 uM of disagreement; the legs read at the centre agree."""
    pd_ = np.linspace(47.6, 52.4, 12); pu = np.r_[np.full(10, 52.0), 49.0, 51.0]
    d = g.updown_difference([50.0], pd_, pd_ - 40.0, pu, pu - 40.0, win=2.5, min_per_leg=3)
    assert abs(d[0]) < 0.05
    pu = np.full(8, 52.0)                                         # every up frame at the stop: too far to read 50
    assert np.isnan(g.updown_difference([50.0], pd_, pd_ - 40.0, pu, pu - 40.0, 2.5, 3)[0])


def test_chl_above_the_calibrations_range_is_flagged(ctd):
    df, meta = ctd
    chl = pd.DataFrame({"pressure_db": [8.0, 9.0, 10.0], "cast": "C1", "chl_corrected": [0.4, 2.6, 1.9],
                        "chl_cal_max": 2.0})
    out, st = g.grid_cast("C1", df, meta, chl=chl)
    r = out.set_index("pressure_db")
    assert r.loc[9.0, "chl_above_cal"] == True and r.loc[8.0, "chl_above_cal"] == False  # noqa: E712
    assert pd.isna(r.loc[11.0, "chl_above_cal"]) and st.chl_above_cal_bins == 1
    assert any("above the calibration's range (2 mg/m3" in n for n in st.notes)
    out2, _ = g.grid_cast("C1", df, meta, chl=chl.drop(columns="chl_cal_max"))     # an older product: not assessed
    assert out2["chl_above_cal"].isna().all()


def test_npq_risk_by_par_beside_the_daylight_rule(ctd):
    """Daylight cast (11:16 UTC): the clock rule marks 8 to 10 dbar. A cloudy day's PAR reaches the threshold only at
    8 dbar, so the two rules disagree at 9 and 10 dbar, and that is reported, not resolved."""
    df, meta = ctd
    df = df.assign(**{"par/sat/log": [250.0, 80.0, 40.0, 20.0]})
    out, st = g.grid_cast("C1", df, meta)
    assert list(out["chl_npq_risk_par"]) == [True, False, False, False]
    assert list(out["chl_npq_risk"]) == [True, True, True, False]               # unchanged
    assert (st.npq_bins_clock, st.npq_bins_par, st.npq_bins_agree, st.npq_bins_both_known) == (3, 1, 2, 4)
    assert any("they disagree on 2 bins" in n for n in st.notes)
    out2, st2 = g.grid_cast("C1", ctd[0], meta)                                  # no PAR sensor: not assessed
    assert out2["chl_npq_risk_par"].isna().all()
    assert st2.npq_bins_clock == 3 and st2.npq_bins_both_known == 0             # the daylight rule still counted
