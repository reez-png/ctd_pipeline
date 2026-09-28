"""
Tests for ctd_lib.py — the pipeline's core parsers and seawater functions.

Run from the project root:
    pip install pytest
    pytest -q

These guard the functions most likely to break silently and corrupt data:
the .cnv parser, the soak+dip Section step (including CRLF preservation), and the
PSS-78 salinity / potential-temperature / depth conversions (checked against
published reference values).
"""

import numpy as np
import pytest

import ctd_lib as L


# --------------------------------------------------------------------------- #
# helpers to build synthetic SBE files (CRLF, like real SBE output)
# --------------------------------------------------------------------------- #
def _write_cnv(path, pressures, temps=None, conds=None, interval=1 / 24.0, nl="\r\n"):
    n = len(pressures)
    temps = [15.0] * n if temps is None else temps
    conds = [4.0] * n if conds is None else conds
    hdr = [
        "* Sea-Bird SBE 19plus Data File:",
        "# nquan = 4",
        f"# nvalues = {n}",
        f"# interval = seconds: {interval:.6f}",
        "# name 0 = scan: Scan Count",
        "# name 1 = prDM: Pressure, Strain Gauge [db]",
        "# name 2 = t090C: Temperature [ITS-90, deg C]",
        "# name 3 = c0S/m: Conductivity [S/m]",
        "*END*",
    ]
    rows = [f"{i+1:8d}{pressures[i]:11.3f}{temps[i]:11.4f}{conds[i]:11.5f}" for i in range(n)]
    path.write_bytes((nl.join(hdr) + nl + nl.join(rows) + nl).encode("latin-1"))


# --------------------------------------------------------------------------- #
# .cnv parsing
# --------------------------------------------------------------------------- #
def test_read_cnv_names_and_shape(tmp_path):
    f = tmp_path / "c.cnv"
    _write_cnv(f, [1.0, 2.0, 3.0])
    names, data = L.read_cnv(f)
    assert names == ["scan", "prDM", "t090C", "c0S/m"]
    assert data.shape == (3, 4)
    assert L.cnv_pressure_column(names) == 1
    np.testing.assert_allclose(data[:, 1], [1.0, 2.0, 3.0])


def test_cnv_sample_rate():
    assert L.cnv_sample_rate(["# interval = seconds: 0.041667"]) == pytest.approx(24.0, abs=0.01)
    assert L.cnv_sample_rate(["# nothing here"]) == 24.0  # fallback


def test_read_cnv_no_end_raises(tmp_path):
    f = tmp_path / "bad.cnv"
    f.write_bytes(b"* header only\r\n# name 0 = prDM: Pressure\r\n")
    with pytest.raises(ValueError):
        L.read_cnv(f)


# --------------------------------------------------------------------------- #
# descent onset + sectioning
# --------------------------------------------------------------------------- #
def test_descent_onset_removes_soak_and_dip():
    fs = 24.0
    p = np.concatenate([
        np.full(200, 1.0),                       # surface soak
        np.linspace(1, 11, 150),                 # equilibration dip down
        np.linspace(11, 2, 150),                 # dip back up
        np.linspace(2, 60, 1500),                # real descent
        np.linspace(60, 0, 2000),                # upcast
    ])
    start = L.descent_onset(p, fs)
    # onset should land at the start of the real descent (~2 db), after the 500 soak+dip scans
    assert 480 <= start <= 520
    assert p[start] < 5.0


def test_section_cnv_preserves_crlf_and_trims(tmp_path):
    n_soak = 300
    pressures = list(np.concatenate([
        np.full(n_soak, 1.0),                    # soak
        np.linspace(1, 40, 1200),                # descent
        np.linspace(40, 0, 1200),                # upcast
    ]))
    src = tmp_path / "in.cnv"
    out = tmp_path / "in_sec.cnv"
    _write_cnv(src, pressures)
    info = L.section_cnv(src, out)

    raw = out.read_bytes()
    assert b"\r\n" in raw                                  # CRLF preserved
    assert raw.replace(b"\r\n", b"").count(b"\n") == 0     # no bare LF introduced
    # nvalues header updated to the kept count
    assert f"# nvalues = {info['rows_kept']}".encode() in raw
    # kept rows are a verbatim tail of the original data rows
    src_rows = src.read_bytes().split(b"*END*\r\n", 1)[1].splitlines(keepends=True)
    src_rows = [r for r in src_rows if r.strip()]
    out_rows = raw.split(b"*END*\r\n", 1)[1].splitlines(keepends=True)
    out_rows = [r for r in out_rows if r.strip()]
    assert out_rows == src_rows[info["onset_scan"]:]
    assert info["rows_removed"] >= n_soak - 24            # roughly the soak was dropped


# --------------------------------------------------------------------------- #
# seawater reference values
# --------------------------------------------------------------------------- #
def test_sal78_standard_seawater():
    # By definition C(35,15,0) = C_STD, so R = 1 at T68=15, P=0 -> S = 35.000
    s = L.sal78(1.0, 15.0, 0.0)
    assert s == pytest.approx(35.0, abs=1e-3)


def test_practical_salinity_from_conductivity():
    # C = C_STD at T90~15, P=0 should give ~35 (T90->T68 tiny correction)
    s = L.practical_salinity(L.C_STD, 15.0, 0.0)
    assert s == pytest.approx(35.0, abs=0.05)


def test_potential_temperature_reference():
    # Fofonoff & Millard (1983) check value is defined in IPTS-68:
    #   S=40, T68=40, P=10000 -> theta68(ref 0) = 36.89073
    # ctd_lib works in ITS-90 (input * 1.00024 -> T68 internally, result / 1.00024 -> T90).
    # So feed T90 = 40/1.00024 (internal T68 = 40) and expect the T90-scaled result.
    t90_in = 40.0 / 1.00024
    theta90 = L.potential_temperature(40.0, t90_in, 10000.0)
    assert theta90 == pytest.approx(36.89073 / 1.00024, abs=2e-3)


def test_pressure_to_depth_reference():
    # UNESCO 1983 check: p=10000 dbar at lat 30 deg -> depth ~ 9712.653 m
    d = L.pressure_to_depth(10000.0, 30.0)
    assert d == pytest.approx(9712.653, abs=0.5)
    # shallow shelf: pressure ~ depth within ~0.5 %
    d50 = L.pressure_to_depth(50.0, 5.0)
    assert abs(d50 - 50.0) < 0.5


def test_salinity_monotonic_in_conductivity():
    lo = L.practical_salinity(3.8, 15.0, 0.0)
    hi = L.practical_salinity(4.6, 15.0, 0.0)
    assert hi > lo
