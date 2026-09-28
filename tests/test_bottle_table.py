"""
Tests for build_bottle_nitrate_table.parse_btl against the real SBE Bottle Summary
.btl layout (no '# name N =' lines; a 'Bottle ... Date ...' title row; a derived
Salinity column whose sdev line is blank). Built from an actual P45_07 .btl so a
future change to the parser can't silently drop bottles again.

Run from the project root:  pytest -q
"""

from pathlib import Path

import build_bottle_nitrate_table as B


# A faithful slice of a real SBE .btl: raw sensor-XML header (no '# name' lines),
# the two-row title header, then avg/sdev pairs. Sal00 is derived -> no sdev value.
_REAL_BTL = (
    "* Sea-Bird SBE19plus  Data File:\r\n"
    "* FileName = C:\\CTD_SUNA_Operations\\03_Cruises\\P45-07\\test1.hex\r\n"
    "# interval = seconds: 0.25\r\n"
    "# datcnv_ox_tau_correction = no\r\n"
    "# datcnv_bottle_scan_range_source = AFM file\r\n"
    "# bottlesum_date = Sep 28 2026 18:10:02, 7.26.7.129\r\n"
    "    Bottle        Date      Sal00       PrdM     Tv290C      C0S/m    Sbeox0V\r\n"
    "  Position        Time                                                       \r\n"
    "      1    Aug 07 2026    35.7661     36.448    19.4747   4.832433     1.5703 (avg)\r\n"
    "              15:23:16                 0.426     0.0040   0.000323     0.0071 (sdev)\r\n"
    "      2    Aug 07 2026    35.8043     30.116    19.6104   4.850774     1.6351 (avg)\r\n"
    "              15:23:26                 0.303     0.0470   0.007811     0.0060 (sdev)\r\n"
)


def _write(tmp_path, text):
    p = tmp_path / "P45_07_CTD_31.btl"
    p.write_bytes(text.encode("latin-1"))
    return p


def test_column_names_from_title_row(tmp_path):
    names, _ = B.parse_btl(_write(tmp_path, _REAL_BTL))
    assert names == ["Sal00", "PrdM", "Tv290C", "C0S/m", "Sbeox0V"]


def test_parses_all_bottles(tmp_path):
    _, recs = B.parse_btl(_write(tmp_path, _REAL_BTL))
    assert [r["bottle_no"] for r in recs] == [1, 2]


def test_values_and_depth(tmp_path):
    _, recs = B.parse_btl(_write(tmp_path, _REAL_BTL))
    r = recs[0]
    assert r["pressure_dbar"] == 36.448
    assert r["temp_c"] == 19.4747
    assert r["salinity"] == 35.7661
    assert r["salinity_source"] == "btl"     # taken from Sal00, not recomputed
    # depth from pressure at the tool's LATITUDE, close to the pressure on the shelf
    assert abs(r["depth_m"] - 36.2) < 0.5


def test_time_from_sdev_line(tmp_path):
    _, recs = B.parse_btl(_write(tmp_path, _REAL_BTL))
    assert recs[0]["time_utc"] == "2026-08-07T15:23:16"
    assert recs[1]["time_utc"] == "2026-08-07T15:23:26"


def _btl_file(dirpath, cast, bottles):
    """Write a minimal real-format .btl. bottles: {no: (pressure, temp, salinity)}."""
    hdr = ("* Sea-Bird SBE19plus  Data File:\r\n# interval = seconds: 0.25\r\n"
           "    Bottle        Date      Sal00       PrdM     Tv290C      C0S/m    Sbeox0V\r\n"
           "  Position        Time\r\n")
    body = ""
    for n, (p, t, s) in bottles.items():
        body += f"      {n}    Aug 07 2026    {s:.4f}     {p:.3f}    {t:.4f}   4.800000     1.5000 (avg)\r\n"
        body += f"              15:2{n}:16                 0.400     0.0040   0.000300     0.0070 (sdev)\r\n"
    (dirpath / f"{cast}.btl").write_bytes((hdr + body).encode("latin-1"))


def test_lab_join_matches_and_reports_unmatched(tmp_path, monkeypatch):
    """The lab-join branch: matched bottles carry nitrate; unmatched (both directions)
    are blanked in the product and listed in the UNMATCHED report."""
    import pandas as pd

    bdir = tmp_path / "Bottle"; bdir.mkdir()
    _btl_file(bdir, "P45_07_CTD_02", {1: (36.4, 19.5, 35.77), 2: (30.1, 19.6, 35.80)})
    _btl_file(bdir, "P45_07_CTD_03", {1: (15.7, 20.9, 35.79), 2: (14.6, 21.0, 35.82)})

    lab = pd.DataFrame({
        "cast_id":    ["P45_07_CTD_02", "P45_07_CTD_02", "P45_07_CTD_03", "P45_07_CTD_99"],
        "bottle_no":  [1, 2, 1, 5],          # last row: lab sample with no matching bottle
        "nitrate_uM": [0.42, 1.10, 3.55, 9.9],
        "station":    ["S1", "S1", "S2", "SX"],
        "flag":       [2, 2, 2, 2],
    })
    labf = tmp_path / "lab.csv"; lab.to_csv(labf, index=False)
    out = tmp_path / "out"; out.mkdir()

    monkeypatch.setattr(B, "BOTTLE_DIR", bdir)
    monkeypatch.setattr(B, "OUTPUT_DIR", out)
    monkeypatch.setattr(B, "LAB_FILE", labf)
    monkeypatch.setattr(B.sys, "argv", ["build_bottle_nitrate_table.py"])  # no CLI arg

    assert B.main() == 0

    nit = pd.read_csv(out / "P45_07_bottle_nitrate.csv")
    assert len(nit) == 4
    assert int(nit["nitrate_uM"].notna().sum()) == 3          # 3 of 4 matched
    # the unmatched bottle (CTD_03, bottle 2) has no nitrate
    unmatched_row = nit[nit["nitrate_uM"].isna()]
    assert len(unmatched_row) == 1

    report = (out / "P45_07_bottle_nitrate_UNMATCHED.csv").read_text()
    assert "P45_07_CTD_03,2" in report      # bottle with no lab
    assert "P45_07_CTD_99,5" in report      # lab row with no bottle


def test_lab_template_emitted_without_lab_file(tmp_path, monkeypatch):
    """With no LAB_FILE, a pre-keyed lab template is written for every bottle."""
    import pandas as pd
    bdir = tmp_path / "Bottle"; bdir.mkdir()
    _btl_file(bdir, "P45_07_CTD_02", {1: (36.4, 19.5, 35.77), 2: (30.1, 19.6, 35.80)})
    out = tmp_path / "out"; out.mkdir()
    monkeypatch.setattr(B, "BOTTLE_DIR", bdir)
    monkeypatch.setattr(B, "OUTPUT_DIR", out)
    monkeypatch.setattr(B, "LAB_FILE", None)
    monkeypatch.setattr(B.sys, "argv", ["build_bottle_nitrate_table.py"])

    assert B.main() == 0
    tmpl = pd.read_csv(out / "P45_07_lab_nitrate_TEMPLATE.csv")
    assert list(tmpl.columns) == ["cast_id", "bottle_no", "depth_m", "time_utc",
                                  "nitrate_uM", "replicate_id", "station", "flag"]
    assert len(tmpl) == 2
    assert list(tmpl["bottle_no"]) == [1, 2]


def test_salinity_computed_when_absent(tmp_path):
    # A .btl WITHOUT a salinity column but with C, T, P -> salinity is computed (PSS-78)
    no_sal = _REAL_BTL.replace(
        "    Bottle        Date      Sal00       PrdM     Tv290C      C0S/m    Sbeox0V\r\n",
        "    Bottle        Date       PrdM     Tv290C      C0S/m    Sbeox0V\r\n",
    )
    # drop the Sal00 value from each avg line so column counts line up
    no_sal = (no_sal
              .replace("35.7661     36.448", "36.448")
              .replace("35.8043     30.116", "30.116"))
    p = tmp_path / "nosal.btl"
    p.write_bytes(no_sal.encode("latin-1"))
    _, recs = B.parse_btl(p)
    assert recs[0]["salinity_source"] == "computed"
    assert 30 < recs[0]["salinity"] < 40      # plausible open-ocean value
