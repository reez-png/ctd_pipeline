"""tools/migrate_layout.py and the cruise-layout guard in project_paths (phase 3, docs/DATA_LEVELS.md).

Builds a small fake cruise in a temporary folder (no real data), migrates it v1 -> v2 and back, and checks
that every byte arrives where the registry says, that conflicts are refused before anything moves, that a
repeated or interrupted run is safe, and that product_dir refuses to mix layouts.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

import project_paths as pp

_spec = importlib.util.spec_from_file_location(
    "migrate_layout", Path(__file__).resolve().parents[1] / "tools" / "migrate_layout.py")
ml = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ml)

CRUISE = "P99_01"
# cruise-relative v1 paths of a small but representative cruise
V1_FILES = [
    "L0/CTD/P99_01_CTD_01.hex",                        # raw: never moves
    "L0/SUNA/A0000001.CSV",
    "L1/CTD/P99_01_CTD_01.cnv",                        # DatCnv: stays in L1/CTD
    "L1/SUNA/A0000001__L1.csv",
    "L2/CTD/Filter/P99_01_CTD_01_filt.cnv",            # SBE stages: L2/CTD/<stage> -> L1/CTD/<stage>
    "L2/CTD/Derived/P99_01_CTD_01_der.cnv",
    "L2/CTD/CTD_1s/P99_01_CTD_01_1s.cnv",
    "L2/CTD/CTD_1m_down/P99_01_CTD_01_1m_down.cnv",    # 1 m grid -> L3/CTD/CTD_1m_down, sub-folder too
    "L2/CTD/CTD_1m_down/asc_semicolon/P99_01_CTD_01_1m_down_ascii.asc",
    "L2/CTD/Bottle/P99_01_CTD_01.btl",                 # stays
    "L2/CTD/combined_all_casts_final_1m_down_profiles.csv",
    "L2/CTD/qc_plots/P99_01_CTD_01_filter.png",
    "L2/SUNA/P99_01_CTD_01_SUNA_1s.csv",               # stays
    "L3/P99_01_CTD_01_nitrate_QC.csv",                 # -> L2/nitrate
    "L3/SUNA/P99_01__suna_L3_corrected.csv",           # -> L2/SUNA_calibrated
    "L3/productivity/P99_01__chl_L3_corrected.csv",    # -> L2/productivity
]
UNCLAIMED = ["L3/notes_by_hand.csv", "L2/CTD/stray.txt"]


def run(repo: Path, *args: str) -> int:
    """migrate_layout.main, always pointed at a throwaway baseline folder (never the real one)."""
    if "--baseline-dir" not in args:
        args += ("--baseline-dir", str(repo / "no_baseline"))
    return ml.main(list(args))


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    monkeypatch.delenv("CTD_PIPELINE_LAYOUT", raising=False)
    cdir = tmp_path / "cruises" / CRUISE
    for i, rel in enumerate(V1_FILES + UNCLAIMED):
        p = cdir / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(f"{rel}\n{i}\n".encode() * (i + 1))
    return tmp_path


def _digest_tree(cdir: Path) -> dict:
    return {p.relative_to(cdir).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in cdir.rglob("*") if p.is_file() and "_audit" not in p.parts}


def _expected_v2(tree_v1: dict) -> dict:
    out = {}
    for rel, h in tree_v1.items():
        _key, new = ml.relocate(rel, "v1", "v2")
        out[new or rel] = h
    return out


# ------------------------------------------------------------------ registry properties the migration relies on
def test_every_moving_folder_is_distinct_across_layouts():
    v1_only = {pp.LAYOUTS["v1"][k] for k in pp.LAYOUTS["v1"] if pp.LAYOUTS["v1"][k] != pp.LAYOUTS["v2"][k]}
    v2_only = {pp.LAYOUTS["v2"][k] for k in pp.LAYOUTS["v2"] if pp.LAYOUTS["v1"][k] != pp.LAYOUTS["v2"][k]}
    assert not v1_only & set(pp.LAYOUTS["v2"].values())
    assert not v2_only & set(pp.LAYOUTS["v1"].values())


@pytest.mark.parametrize("lay", ["v1", "v2"])
def test_each_product_folder_classifies_to_itself(lay):
    for key, folder in pp.LAYOUTS[lay].items():
        glob = pp.MIGRATE_GLOB[key]
        name = "x.csv" if glob == "*" else glob.replace("*", "P99_01_CTD_01")
        assert ml.classify(f"{folder}/{name}", lay) == key


# ------------------------------------------------------------------ the migration itself
def test_dry_run_moves_nothing(repo):
    before = _digest_tree(repo / "cruises" / CRUISE)
    assert run(repo) == 0
    assert _digest_tree(repo / "cruises" / CRUISE) == before
    assert pp.cruise_layout(CRUISE) == "v1"


def test_v1_to_v2_moves_every_byte_where_the_registry_says(repo, capsys):
    cdir = repo / "cruises" / CRUISE
    before = _digest_tree(cdir)
    assert run(repo, "--apply") == 0
    after = _digest_tree(cdir)
    assert after == _expected_v2(before)
    assert "L1/CTD/Filter/P99_01_CTD_01_filt.cnv" in after
    assert "L3/CTD/CTD_1m_down/asc_semicolon/P99_01_CTD_01_1m_down_ascii.asc" in after
    assert "L2/nitrate/P99_01_CTD_01_nitrate_QC.csv" in after
    assert "_figures/qc/P99_01_CTD_01_filter.png" in after
    out = capsys.readouterr().out
    for rel in UNCLAIMED:                                       # never touched, but reported
        assert rel in after and rel in out
    assert not (cdir / "L2" / "CTD" / "Filter").exists()        # emptied source folders are removed
    assert (cdir / "L2" / "CTD" / "Bottle").is_dir()            # unchanged products stay
    assert pp.cruise_layout(CRUISE) == "v2"
    logs = list((cdir / "_audit").glob("layout_migration_*.csv"))
    assert len(logs) == 1
    rows = logs[0].read_text(encoding="utf-8").strip().splitlines()[1:]
    assert len(rows) == sum(1 for r in before if ml.relocate(r, "v1", "v2")[1]) and all(",moved" in r for r in rows)


def test_second_run_is_a_no_op(repo):
    assert run(repo, "--apply") == 0
    tree = _digest_tree(repo / "cruises" / CRUISE)
    assert run(repo, "--apply") == 0
    assert _digest_tree(repo / "cruises" / CRUISE) == tree


def test_round_trip_restores_the_original_tree(repo):
    cdir = repo / "cruises" / CRUISE
    before = _digest_tree(cdir)
    assert run(repo, "--apply") == 0
    assert run(repo, "--to", "v1", "--apply") == 0
    assert _digest_tree(cdir) == before
    assert pp.cruise_layout(CRUISE) == "v1"


def test_conflict_refuses_the_whole_cruise(repo):
    cdir = repo / "cruises" / CRUISE
    clash = cdir / "L1" / "CTD" / "Filter" / "P99_01_CTD_01_filt.cnv"
    clash.parent.mkdir(parents=True)
    clash.write_text("different content")
    before = _digest_tree(cdir)
    assert run(repo, "--apply") == 1
    assert _digest_tree(cdir) == before                         # nothing moved, not even unrelated files
    assert not (cdir / pp.LAYOUT_MARKER).exists()


def test_interrupted_run_resumes(repo):
    cdir = repo / "cruises" / CRUISE
    src = cdir / "L3" / "P99_01_CTD_01_nitrate_QC.csv"
    dst = cdir / "L2" / "nitrate" / src.name
    dst.parent.mkdir(parents=True)
    dst.write_bytes(src.read_bytes())                           # copied but source not yet removed
    before = _digest_tree(cdir)
    before.pop("L2/nitrate/" + src.name)
    assert run(repo, "--apply") == 0
    assert _digest_tree(cdir) == _expected_v2(before)


# ------------------------------------------------------------------ the guard in product_dir
def test_product_dir_refuses_to_mix_layouts(repo, monkeypatch):
    assert pp.product_dir(CRUISE, "nitrate_qc") == repo / "cruises" / CRUISE / "L3"   # v1 data, v1 code
    monkeypatch.setenv("CTD_PIPELINE_LAYOUT", "v2")
    with pytest.raises(RuntimeError, match="stored in layout v1"):
        pp.product_dir(CRUISE, "nitrate_qc")                    # v2 code must not write beside v1 data
    assert run(repo, "--apply") == 0
    assert pp.product_dir(CRUISE, "nitrate_qc") == repo / "cruises" / CRUISE / "L2" / "nitrate"
    monkeypatch.setenv("CTD_PIPELINE_LAYOUT", "v1")
    with pytest.raises(RuntimeError, match="stored in layout v2"):
        pp.product_dir(CRUISE, "nitrate_qc")
    assert pp.product_dir(CRUISE, "nitrate_qc", layout="v1").name == "L3"   # explicit layout: no guard


def test_fresh_cruise_has_no_layout(repo):
    (repo / "cruises" / "P99_02" / "L0" / "CTD").mkdir(parents=True)
    assert pp.cruise_layout("P99_02") is None
    assert pp.product_dir("P99_02", "ctd_raw").name == "CTD"


# ------------------------------------------------------------------ frozen baseline
def test_baseline_is_remapped_with_its_cruise(repo):
    cdir = repo / "cruises" / CRUISE
    bdir = repo / "baseline"
    files = []
    for rel in ("L3/P99_01_CTD_01_nitrate_QC.csv", "L2/SUNA/P99_01_CTD_01_SUNA_1s.csv"):
        (bdir / rel).parent.mkdir(parents=True, exist_ok=True)
        (bdir / rel).write_bytes((cdir / rel).read_bytes())
        files.append({"path": rel, "sha256": hashlib.sha256((bdir / rel).read_bytes()).hexdigest()})
    (bdir / "manifest.json").write_text(json.dumps({"cruise": CRUISE, "files": files, "products": []}))
    assert run(repo, "--apply", "--baseline-dir", str(bdir)) == 0
    man = json.loads((bdir / "manifest.json").read_text())
    assert man["layout"] == "v2"
    paths = [f["path"] for f in man["files"]]
    assert paths == ["L2/nitrate/P99_01_CTD_01_nitrate_QC.csv", "L2/SUNA/P99_01_CTD_01_SUNA_1s.csv"]
    for f in man["files"]:                                      # frozen copy and live product line up again
        assert (bdir / f["path"]).read_bytes() == (cdir / f["path"]).read_bytes()
        assert hashlib.sha256((bdir / f["path"]).read_bytes()).hexdigest() == f["sha256"]
    assert (bdir / "manifest.v1.json").is_file()
