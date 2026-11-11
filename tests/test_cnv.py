import numpy as np
from suncat.io import cnv

# A tiny synthetic Sea-Bird .cnv: a header (with a sampling interval), three scans, and on the last
# scan a bad-flag sentinel in a DATA column (tv290C) and in the flag column. No real cruise data.
SYNTH = (
    "* Sea-Bird SBE 19plus V2 synthetic test file\n"
    "# interval = seconds: 0.25\n"
    "# name 0 = prdM: Pressure, Strain Gauge [db]\n"
    "# name 1 = tv290C: Temperature [ITS-90, deg C]\n"
    "# name 2 = flag: flag\n"
    "*END*\n"
    "   1.0  28.00  0.0000e+00\n"
    "   2.0  27.50  0.0000e+00\n"
    "   3.0  -9.990e-29  -9.990e-29\n"
)


def _write(tmp_path):
    p = tmp_path / "synthetic.cnv"
    p.write_text(SYNTH, encoding="latin-1")
    return p


def test_columns_and_scans(tmp_path):
    df = cnv.read_cnv(_write(tmp_path))
    assert list(df.columns) == ["prdM", "tv290C", "flag"]
    assert len(df) == 3
    assert df.loc[0, "prdM"] == 1.0
    assert df.loc[1, "tv290C"] == 27.5


def test_badflag_masking(tmp_path):
    p = _write(tmp_path)
    raw = cnv.read_cnv(p, mask_badflag=False)
    masked = cnv.read_cnv(p, mask_badflag=True)
    assert raw.loc[2, "tv290C"] < -8e-29              # unmasked keeps the sentinel
    assert np.isnan(masked.loc[2, "tv290C"])          # masked -> NaN in a data column
    assert masked.loc[2, "flag"] < -8e-29             # but the flag column keeps its sentinel


def test_read_cnv_full_unpacks_like_the_old_tuple(tmp_path):
    header, cols, data, rows_raw = cnv.read_cnv_full(_write(tmp_path))
    assert any("*END*" in h.upper() for h in header)
    assert cols == ["prdM", "tv290C", "flag"]
    assert len(data) == 3
    assert len(rows_raw) == 3


def test_sampling_metadata_in_attrs(tmp_path):
    # step00c's resolve_sample_rate reads these from df.attrs
    df = cnv.read_cnv(_write(tmp_path))
    assert df.attrs["interval_kind"] == "seconds"
    assert df.attrs["interval_value"] == 0.25
    assert df.attrs["sample_rate_hz"] == 4.0          # 1 / 0.25
    assert df.attrs["source_path"].endswith("synthetic.cnv")
