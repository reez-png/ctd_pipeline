"""suncat.timeutil: seconds from times must not depend on how pandas stores them (the astype("int64") bug)."""
from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import pytest

from suncat.timeutil import elapsed_seconds, epoch_seconds, to_utc

TEXT = ["2026-08-04 11:27:56", "2026-08-04 11:27:57.100000", "2026-08-04 11:27:59+00:00"]
EXPECT = np.array([0.0, 1.1, 3.0])


@pytest.mark.parametrize("unit", ["s", "ms", "us", "ns"])
def test_epoch_seconds_is_the_same_in_every_resolution(unit):
    t = pd.Series(pd.to_datetime(["2026-08-04 11:27:56", "2026-08-04 11:27:57.5"], format="ISO8601")).astype(f"datetime64[{unit}]")
    e = epoch_seconds(t)
    assert e[0] == pd.Timestamp("2026-08-04 11:27:56", tz="UTC").timestamp()
    assert abs((e[1] - e[0]) - (1.0 if unit == "s" else 1.5)) < 1e-6


def test_the_old_expression_is_wrong_in_microseconds():
    # the bug this module replaces: astype("int64") / 1e9 is seconds only for nanosecond storage
    t = pd.Series(pd.to_datetime(["2026-08-04 11:27:56", "2026-08-04 11:27:57"])).astype("datetime64[us]")
    old = t.astype("int64").to_numpy() / 1e9
    assert np.diff(old)[0] == pytest.approx(0.001)
    assert np.diff(epoch_seconds(t))[0] == pytest.approx(1.0)


def test_mixed_text_formats_are_all_read():
    assert np.allclose(elapsed_seconds(TEXT), EXPECT)


def test_python_datetimes_and_tz():
    base = datetime(2026, 8, 4, 11, 27, 56)
    py = [base + timedelta(seconds=s) for s in (0, 1.1, 3)]
    assert np.allclose(elapsed_seconds(py), EXPECT)
    aware = pd.Series(pd.to_datetime(TEXT, format="ISO8601", utc=True)).dt.tz_convert("Africa/Accra")
    assert np.allclose(epoch_seconds(aware), epoch_seconds(TEXT))


def test_unreadable_values_become_nan():
    e = epoch_seconds(["2026-08-04 11:27:56", "not a time", None])
    assert np.isfinite(e[0]) and np.isnan(e[1:]).all()
    assert to_utc(["x"]).isna().all()


def test_non_iso_text_is_read_by_the_second_pass():
    t = to_utc(["2026-08-04 11:27:56", "08/04/2026 11:27:57", "Aug 04 2026 11:27:59"])
    assert t.notna().all()
    assert np.allclose(elapsed_seconds(t), [0.0, 1.0, 3.0])


def test_unreadable_times_are_counted_not_dropped_silently(capsys):
    to_utc(["2026-08-04 11:27:56", "garbage", "", None], what="merged 1 s")
    out = capsys.readouterr().out
    assert "1 of 4 merged 1 s times unreadable" in out and "'garbage'" in out
    to_utc(["2026-08-04 11:27:56"], what="clean")
    assert capsys.readouterr().out == ""


def test_the_old_guess_loses_rows_that_iso8601_keeps():
    # pandas guesses the format from the first value; coerce then turns the other form into NaT
    col = pd.Series(["2026-08-04 11:27:56+00:00", "2026-08-04 11:27:56.500000+00:00", "2026-08-04 11:27:57+00:00"])
    old = pd.to_datetime(col, utc=True, errors="coerce")
    assert old.isna().sum() == 1
    assert to_utc(col).notna().all()
