"""suncat.suna_clock: the per-cruise SUNA clock offset file."""
from __future__ import annotations

import pytest

from suncat.suna_clock import describe, load_offsets, offsets_path, write_offsets


def test_absent_file_means_no_offset(tmp_path):
    assert load_offsets(tmp_path) == {}
    assert "no suna_time_offsets.csv" in describe(tmp_path)


def test_write_then_load(tmp_path):
    p = write_offsets(tmp_path, ["C_02", "C_01"], -16.0, "water entry (step04): +16.1 s")
    assert p == offsets_path(tmp_path)
    assert load_offsets(tmp_path) == {"C_01": -16.0, "C_02": -16.0}
    d = describe(tmp_path)
    assert "2 casts" in d and "-16 s" in d and "water entry" in d


def test_refuses_to_overwrite_and_needs_a_source(tmp_path):
    write_offsets(tmp_path, ["C_01"], -16.0, "x")
    with pytest.raises(FileExistsError):
        write_offsets(tmp_path, ["C_01"], -10.0, "y")
    write_offsets(tmp_path, ["C_01"], -10.0, "y", overwrite=True)
    assert load_offsets(tmp_path) == {"C_01": -10.0}
    with pytest.raises(ValueError):
        write_offsets(tmp_path / "other", ["C_01"], -10.0, "  ")


def test_malformed_files_are_errors(tmp_path):
    p = offsets_path(tmp_path)
    p.parent.mkdir(parents=True)
    p.write_text("cast_id,offset_s\nC_01,abc\n")
    with pytest.raises(ValueError, match="finite"):
        load_offsets(tmp_path)
    p.write_text("cast_id,offset_s\nC_01,-16\nC_01,-15\n")
    with pytest.raises(ValueError, match="twice"):
        load_offsets(tmp_path)
    p.write_text("cast,offset\nC_01,-16\n")
    with pytest.raises(ValueError, match="columns"):
        load_offsets(tmp_path)
