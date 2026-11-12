"""Tests for suncat.sbe.resolve: executable + PSA resolution, and the pure output-filename builder.

Filesystem cases use pytest's tmp_path, so they are isolated and run anywhere (including CI).
"""
from __future__ import annotations
import pytest

from suncat.sbe.resolve import find_executable, resolve_psa, module_output_file


# --- module_output_file (pure) ---------------------------------------------

def test_output_file_explicit_template():
    m = {"name": "01_datcnv", "output_file": "{cast_id}.cnv"}
    assert module_output_file(m, "P45_06_CTD_03") == "P45_06_CTD_03.cnv"


def test_output_file_from_suffix():
    m = {"name": "08b_binavg_1m_down", "output_suffix": "_1m_down.cnv"}
    assert module_output_file(m, "P45_06_CTD_03") == "P45_06_CTD_03_1m_down.cnv"


def test_output_file_requires_one_or_the_other():
    with pytest.raises(ValueError):
        module_output_file({"name": "x"}, "P45_06_CTD_03")


# --- find_executable -------------------------------------------------------

def test_find_executable_direct(tmp_path):
    (tmp_path / "DatCnvW.exe").write_text("")
    assert find_executable("DatCnvW.exe", tmp_path) == tmp_path / "DatCnvW.exe"


def test_find_executable_recursive(tmp_path):
    sub = tmp_path / "bin"
    sub.mkdir()
    (sub / "FilterW.exe").write_text("")
    assert find_executable("FilterW.exe", tmp_path) == sub / "FilterW.exe"


def test_find_executable_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        find_executable("Nope.exe", tmp_path)


def test_find_executable_bad_root_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        find_executable("x.exe", tmp_path / "does_not_exist")


# --- resolve_psa -----------------------------------------------------------

def test_resolve_psa_direct_candidate(tmp_path):
    (tmp_path / "01_datcnv.psa").write_text("")
    m = {"name": "01_datcnv", "psa_candidates": ["01_datcnv.psa"]}
    assert resolve_psa(m, tmp_path) == tmp_path / "01_datcnv.psa"


def test_resolve_psa_canonical_is_tried(tmp_path):
    (tmp_path / "05_filter.psa").write_text("")
    m = {"name": "05_filter", "psa_canonical": "05_filter.psa"}
    assert resolve_psa(m, tmp_path) == tmp_path / "05_filter.psa"


def test_resolve_psa_recursive_match(tmp_path):
    sub = tmp_path / "nested"
    sub.mkdir()
    (sub / "07_derive.psa").write_text("")
    m = {"name": "07_derive", "psa_candidates": ["07_derive.psa"]}
    assert resolve_psa(m, tmp_path) == sub / "07_derive.psa"


def test_resolve_psa_token_fallback(tmp_path):
    # No candidate matches, but a .psa filename contains the module token ("loopedit").
    (tmp_path / "my_loopedit_custom.psa").write_text("")
    m = {"name": "06_loopedit", "psa_candidates": []}
    assert resolve_psa(m, tmp_path) == tmp_path / "my_loopedit_custom.psa"


def test_resolve_psa_not_found_raises(tmp_path):
    m = {"name": "99_missing", "psa_candidates": ["none.psa"]}
    with pytest.raises(FileNotFoundError):
        resolve_psa(m, tmp_path)
