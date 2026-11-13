"""Pins the data-level layout registry in project_paths (see docs/DATA_LEVELS.md).

v1 must equal the folders the steps hardcode today, so moving a step onto product_dir() under v1 cannot
change where it reads or writes. v2 must follow the decisions agreed with Amala and Edem. Every
location in either layout must be git-ignored, because these folders hold cruise data.
"""
from __future__ import annotations

import pytest

import project_paths as pp

# The folders the notebooks and scripts build by hand today (inventory, Nov 2026).
V1_TODAY = {
    "ctd_raw": "L0/CTD", "suna_raw": "L0/SUNA", "ctd_datcnv": "L1/CTD",
    "ctd_filter": "L2/CTD/Filter", "ctd_align": "L2/CTD/Align_CTD", "ctd_celltm": "L2/CTD/Cell_tm",
    "ctd_section": "L2/CTD/Section", "ctd_loopedit": "L2/CTD/Loop_Edit", "ctd_derived": "L2/CTD/Derived",
    "ctd_1s": "L2/CTD/CTD_1s", "suna_l1": "L1/SUNA", "suna_l2": "L2/SUNA", "ctd_bottle": "L2/CTD/Bottle",
    "nitrate_qc": "L3", "suna_calibrated": "L3/SUNA", "chl_calibrated": "L3/productivity",
    "chl_calibrated_noNPQ": "L3/productivity_noNPQ", "ctd_grid_1m_down": "L2/CTD/CTD_1m_down",
    "grid_1m_down": "L3/gridded", "qc_plots": "L2/CTD/qc_plots", "profile_plots": "L2/CTD/profile_plots",
    "profiles_table": "L2/CTD",
}
SBE_STAGES = ("ctd_filter", "ctd_align", "ctd_celltm", "ctd_section", "ctd_loopedit", "ctd_derived")


def test_layouts_cover_the_same_products():
    assert set(pp.LAYOUTS["v1"]) == set(pp.LAYOUTS["v2"]) == set(pp.MIGRATE_GLOB)


def test_v1_is_exactly_todays_folders():
    assert pp.LAYOUTS["v1"] == V1_TODAY


def test_v2_follows_the_agreed_decisions():
    v2 = pp.LAYOUTS["v2"]
    for k in SBE_STAGES + ("ctd_1s",):                     # decisions 1 and 2: SBE stages and 1 s base in L1
        assert v2[k].startswith("L1/CTD/"), k
    for k in ("suna_l2", "nitrate_qc", "suna_calibrated", "chl_calibrated", "chl_calibrated_noNPQ"):
        assert v2[k].startswith("L2/"), k                  # merged, calibrated and flagged products in L2
    gridded = {k for k, v in v2.items() if v.startswith("L3")}
    assert gridded == {"ctd_grid_1m_down", "grid_1m_down", "profiles_table"}  # decisions 3 and 4: only down-cast grids in L3
    assert not any(v.startswith("L2/CTD/") and k in SBE_STAGES for k, v in v2.items())


@pytest.mark.parametrize("lay", ["v1", "v2"])
def test_no_two_products_share_a_folder(lay):
    dirs = list(pp.LAYOUTS[lay].values())
    assert len(dirs) == len(set(dirs)), f"{lay}: a folder is used by two products"


def test_product_dir_resolves_both_layouts(tmp_path, monkeypatch):
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    base = tmp_path / "cruises" / "P45_06"
    assert pp.product_dir("P45_06", "nitrate_qc", layout="v1") == base / "L3"
    assert pp.product_dir("P45_06", "nitrate_qc", layout="v2") == base / "L2" / "nitrate"
    made = pp.product_dir("P45_06", "ctd_1s", layout="v2", create=True)
    assert made == base / "L1" / "CTD" / "CTD_1s" and made.is_dir()
    monkeypatch.setenv("CTD_PIPELINE_LAYOUT", "v2")
    assert pp.product_dir("P45_06", "ctd_grid_1m_down") == base / "L3" / "CTD" / "CTD_1m_down"


def test_default_layout_is_v1_until_the_flip(monkeypatch):
    monkeypatch.delenv("CTD_PIPELINE_LAYOUT", raising=False)
    assert pp.active_layout() == "v1"


def test_unknown_names_fail_loudly(monkeypatch):
    monkeypatch.delenv("CTD_PIPELINE_LAYOUT", raising=False)
    with pytest.raises(KeyError, match="known products"):
        pp.product_dir("P45_06", "nitrate")
    with pytest.raises(ValueError):
        pp.product_dir("P45_06", "nitrate_qc", layout="v3")
    monkeypatch.setenv("CTD_PIPELINE_LAYOUT", "v9")
    with pytest.raises(ValueError):
        pp.active_layout()


def test_every_location_is_git_ignored(monkeypatch):
    monkeypatch.delenv("CTD_PIPELINE_ROOT", raising=False)
    rules = (pp.project_root() / ".gitignore").read_text(encoding="utf-8").splitlines()
    rules = {r.strip() for r in rules}
    for lay, table in pp.LAYOUTS.items():
        for key, rel in table.items():
            top = rel.split("/")[0]
            assert f"cruises/*/{top}/" in rules, f"{lay}:{key} lives in cruises/*/{top}/, which .gitignore does not cover"
