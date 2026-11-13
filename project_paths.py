"""project_paths.py

One place to resolve shared project inputs, so no script hardcodes a path.

Any step, notebook, or helper can do:

    from project_paths import find_calibration, find_xmlcon, find_suna_cal, cruise_dir

and stop writing literal paths. Calibration files resolve to the shared
``metadata/calibration/`` folder first, then fall back to a cruise's own
``Instrument_Calibration_Files/`` if a shared copy is not there. That fallback
makes the transition safe: you can delete the per-cruise duplicates and every
caller keeps working as long as the shared copy exists, and a cruise that truly
needs its own calibration can still keep a local copy.

The project root is found by walking up from this file until a marker is seen,
so it does not matter whether a notebook runs from the repo root or from
``processing_scripts/notebook/``. Set the environment variable
``CTD_PIPELINE_ROOT`` to override.

Pure standard library, no third-party imports.
"""
from __future__ import annotations

import json
import os
import warnings
from pathlib import Path

# Files/dirs that only exist at the repo root, used to locate it.
_ROOT_MARKERS = ("ctd_lib.py", ".git", "requirements.txt", "metadata")

# Discrete-sample subfolders under a cruise's metadata/ (see README).
_DISCRETE_KINDS = ("bottle_nitrate", "chla", "winkler_oxygen")


def project_root(start: str | os.PathLike | None = None) -> Path:
    """Return the repository root.

    Order: the ``CTD_PIPELINE_ROOT`` environment variable if set, then a walk up
    from ``start`` (or this file) until a directory containing a root marker is
    found, then the current working directory as a last resort.
    """
    env = os.environ.get("CTD_PIPELINE_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    here = Path(start).resolve() if start is not None else Path(__file__).resolve()
    for cand in (here, *here.parents):
        if cand.is_dir() and any((cand / m).exists() for m in _ROOT_MARKERS):
            return cand
    return Path.cwd().resolve()


def calibration_dir() -> Path:
    """The shared instrument-reference folder: ``metadata/calibration/``."""
    return project_root() / "metadata" / "calibration"


def cruise_dir(cruise_id: str) -> Path:
    """``cruises/<cruise_id>/``."""
    return project_root() / "cruises" / cruise_id


def cruise_metadata_dir(cruise_id: str) -> Path:
    """``cruises/<cruise_id>/metadata/``."""
    return cruise_dir(cruise_id) / "metadata"


def find_calibration(
    filename: str | None = None,
    *,
    cruise_id: str | None = None,
    pattern: str | None = None,
) -> Path:
    """Resolve an instrument calibration or config file.

    Pass either an exact ``filename`` (for example ``"SNA2503C.cal"``) or a glob
    ``pattern`` (for example ``"*.xmlcon"``). Search order:

    1. shared ``metadata/calibration/``
    2. the cruise's ``Instrument_Calibration_Files/`` (only if ``cruise_id`` given)

    Returns the first match. Raises ``FileNotFoundError`` listing where it looked.
    """
    if not filename and not pattern:
        raise ValueError("find_calibration needs filename= or pattern=")

    search_dirs = [calibration_dir()]
    if cruise_id:
        search_dirs.append(cruise_dir(cruise_id) / "Instrument_Calibration_Files")

    for d in search_dirs:
        if not d.is_dir():
            continue
        if filename and (d / filename).exists():
            return d / filename
        if pattern:
            hits = sorted(d.glob(pattern))
            if hits:
                return hits[0]

    looked = "; ".join(str(d) for d in search_dirs)
    target = filename if filename else pattern
    raise FileNotFoundError(f"calibration not found ({target}). Looked in: {looked}")


def find_xmlcon(cruise_id: str | None = None, filename: str | None = None) -> Path:
    """The CTD config (``*.xmlcon``). Exact ``filename`` wins, else first ``*.xmlcon``."""
    if filename:
        return find_calibration(filename, cruise_id=cruise_id)
    return find_calibration(pattern="*.xmlcon", cruise_id=cruise_id)


def find_suna_cal(cruise_id: str | None = None, filename: str | None = None) -> Path:
    """The SUNA calibration (``*.cal``). Exact ``filename`` wins, else first ``*.cal``."""
    if filename:
        return find_calibration(filename, cruise_id=cruise_id)
    return find_calibration(pattern="*.cal", cruise_id=cruise_id)


def find_discrete(cruise_id: str, kind: str, pattern: str | None = None) -> Path:
    """Resolve a discrete-sample table under ``cruises/<cruise_id>/metadata/<kind>/``.

    ``kind`` is one of ``bottle_nitrate``, ``chla``, ``winkler_oxygen``. Template
    files (``*_TEMPLATE.csv``) are skipped. Give ``pattern`` to narrow the match;
    otherwise the first non-template ``*.csv`` is returned.
    """
    if kind not in _DISCRETE_KINDS:
        raise ValueError(f"kind must be one of {_DISCRETE_KINDS}, got {kind!r}")
    d = cruise_metadata_dir(cruise_id) / kind
    if not d.is_dir():
        raise FileNotFoundError(f"no {kind} folder for {cruise_id}: {d}")
    hits = [
        p
        for p in sorted(d.glob(pattern or "*.csv"))
        if not p.name.upper().endswith("_TEMPLATE.CSV")
    ]
    if not hits:
        raise FileNotFoundError(f"no {kind} table (non-template .csv) in {d}")
    return hits[0]


def find_station_coords(cruise_id: str) -> Path:
    """The step06 station coordinate table ``cruises/<cruise_id>/metadata/station_coords.csv``."""
    p = cruise_metadata_dir(cruise_id) / "station_coords.csv"
    if not p.exists():
        raise FileNotFoundError(f"station_coords.csv not found for {cruise_id}: {p}")
    return p


def suna_cast_map(cruise_id: str) -> Path:
    """``cruises/<cruise_id>/metadata/suna_cast_map.csv``."""
    return cruise_metadata_dir(cruise_id) / "suna_cast_map.csv"


# ---------------------------------------------------------------------------
# Data-level layout (see docs/DATA_LEVELS.md)
# ---------------------------------------------------------------------------
# Every step asks for a product folder by NAME, product_dir(cruise_id, "nitrate_qc"), instead of
# building "L2"/"L3" paths itself, so the level scheme lives in exactly one place.
#
#   v1  the layout the pipeline was built on (binned at L2, calibrated at L3)
#   v2  the layout agreed with Amala and Edem, Nov 2026: calibrate and flag BEFORE binning
#       L1 = cleaned per cast (native resolution), L2 = merged, calibrated and flagged (native
#       resolution), L3 = gridded (1 m down-cast bins of calibrated, flagged data)
#
# ACTIVE_LAYOUT selects one; the environment variable CTD_PIPELINE_LAYOUT overrides it. The migration
# glob says which files belong to a product when its folder is shared or moved.

#   key                     v1 (built on)              v2 (agreed)                migrate glob
PRODUCTS = (
    ("ctd_raw",              "L0/CTD",                 "L0/CTD",                  "*"),
    ("suna_raw",             "L0/SUNA",                "L0/SUNA",                 "*"),
    ("ctd_datcnv",           "L1/CTD",                 "L1/CTD",                  "*"),
    ("ctd_filter",           "L2/CTD/Filter",          "L1/CTD/Filter",           "*"),
    ("ctd_align",            "L2/CTD/Align_CTD",       "L1/CTD/Align_CTD",        "*"),
    ("ctd_celltm",           "L2/CTD/Cell_tm",         "L1/CTD/Cell_tm",          "*"),
    ("ctd_section",          "L2/CTD/Section",         "L1/CTD/Section",          "*"),
    ("ctd_loopedit",         "L2/CTD/Loop_Edit",       "L1/CTD/Loop_Edit",        "*"),
    ("ctd_derived",          "L2/CTD/Derived",         "L1/CTD/Derived",          "*"),
    ("ctd_1s",               "L2/CTD/CTD_1s",          "L1/CTD/CTD_1s",           "*"),
    ("suna_l1",              "L1/SUNA",                "L1/SUNA",                 "*"),
    ("suna_l2",              "L2/SUNA",                "L2/SUNA",                 "*"),
    ("ctd_bottle",           "L2/CTD/Bottle",          "L2/CTD/Bottle",           "*"),
    ("nitrate_qc",           "L3",                     "L2/nitrate",              "*_nitrate_QC.csv"),
    ("suna_calibrated",      "L3/SUNA",                "L2/SUNA_calibrated",      "*"),
    ("chl_calibrated",       "L3/productivity",        "L2/productivity",         "*"),
    ("chl_calibrated_noNPQ", "L3/productivity_noNPQ",  "L2/productivity_noNPQ",   "*"),
    ("ctd_grid_1m_down",     "L2/CTD/CTD_1m_down",     "L3/CTD/CTD_1m_down",      "*"),
    ("grid_1m_down",         "L3/gridded",             "L3/gridded",              "*"),
    ("qc_plots",             "L2/CTD/qc_plots",        "_figures/qc",             "*"),
    ("profile_plots",        "L2/CTD/profile_plots",   "_figures/profiles",       "*"),
    ("profiles_table",       "L2/CTD",                 "L3/CTD",                  "combined_all_casts_*.csv"),
)
LAYOUTS = {
    "v1": {key: v1 for key, v1, _v2, _g in PRODUCTS},
    "v2": {key: v2 for key, _v1, v2, _g in PRODUCTS},
}
MIGRATE_GLOB = {key: g for key, _v1, _v2, g in PRODUCTS}

ACTIVE_LAYOUT = "v1"   # flips to "v2" only after every step uses product_dir (docs/DATA_LEVELS.md)


def active_layout() -> str:
    """The layout in force: ``CTD_PIPELINE_LAYOUT`` if set, else ``ACTIVE_LAYOUT``."""
    lay = os.environ.get("CTD_PIPELINE_LAYOUT", ACTIVE_LAYOUT)
    if lay not in LAYOUTS:
        raise ValueError(f"unknown layout {lay!r}; choose one of {sorted(LAYOUTS)}")
    return lay


def product_dir(cruise_id: str, key: str, *, layout: str | None = None, create: bool = False) -> Path:
    """Folder that holds one named product for a cruise, under the active (or given) layout.

    >>> product_dir("P45_06", "nitrate_qc")            # doctest: +SKIP
    ...cruises/P45_06/L3                (v1)   or   ...cruises/P45_06/L2/nitrate   (v2)
    """
    lay = layout or active_layout()
    if lay not in LAYOUTS:
        raise ValueError(f"unknown layout {lay!r}; choose one of {sorted(LAYOUTS)}")
    try:
        rel = LAYOUTS[lay][key]
    except KeyError:
        raise KeyError(f"unknown product {key!r}; known products: {sorted(LAYOUTS[lay])}") from None
    if layout is None:
        _check_cruise_layout(cruise_id, lay)
    path = cruise_dir(cruise_id) / rel
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


# ---------------------------------------------------------------------------
# Which layout a cruise's files are actually in
# ---------------------------------------------------------------------------
# tools/migrate_layout.py writes this marker after moving a cruise. Without it, the layout is inferred
# from which layout's own folders hold files. product_dir() refuses to hand out a folder in one layout
# for a cruise stored in the other, so a step can never scatter new outputs across both.
LAYOUT_MARKER = Path("_audit") / "layout.json"


def _has_files(folder: Path, glob: str) -> bool:
    if not folder.is_dir():
        return False
    found = folder.rglob("*") if glob == "*" else folder.glob(glob)
    return any(p.is_file() for p in found)


def cruise_layout(cruise_id: str) -> str | None:
    """Layout a cruise's files are stored in: "v1", "v2", "mixed", or None (nothing layout-specific yet).

    The migration marker wins; otherwise each product that moves between layouts is checked for files in
    its v1 folder and in its v2 folder."""
    cdir = cruise_dir(cruise_id)
    marker = cdir / LAYOUT_MARKER
    if marker.is_file():
        try:
            return json.loads(marker.read_text(encoding="utf-8"))["layout"]
        except (ValueError, KeyError, OSError):
            pass                                    # unreadable marker: fall back to inference
    found = set()
    for _key, v1, v2, glob in PRODUCTS:
        if v1 == v2:
            continue
        if _has_files(cdir / v1, glob):
            found.add("v1")
        if _has_files(cdir / v2, glob):
            found.add("v2")
    if len(found) > 1:
        return "mixed"
    return found.pop() if found else None


def _check_cruise_layout(cruise_id: str, lay: str) -> None:
    stored = cruise_layout(cruise_id)
    if stored is None or stored == lay:
        return
    fix = f"run  python tools/migrate_layout.py --cruise {cruise_id}  (dry run first; docs/DATA_LEVELS.md)"
    if stored == "mixed":
        msg = f"cruise {cruise_id} has files in both the v1 and v2 folders; {fix}"
        if lay == "v1":                             # the status quo before the flip: warn, keep working
            warnings.warn(msg, stacklevel=3)
            return
        raise RuntimeError(msg)
    raise RuntimeError(
        f"cruise {cruise_id} is stored in layout {stored} but the code is set to {lay}. Either migrate it "
        f"({fix.replace('--cruise', f'--to {lay} --cruise')}) or read it as it is by setting the environment "
        f"variable CTD_PIPELINE_LAYOUT={stored}.")


if __name__ == "__main__":
    # Quick self-report, useful to confirm the root is found correctly.
    print("project_root     :", project_root())
    print("calibration_dir  :", calibration_dir(), "(exists:", calibration_dir().is_dir(), ")")
    try:
        print("SUNA cal (shared):", find_suna_cal())
    except FileNotFoundError as e:
        print("SUNA cal         :", e)
    try:
        print("xmlcon (shared)  :", find_xmlcon())
    except FileNotFoundError as e:
        print("xmlcon           :", e)
    print("data-level layout:", active_layout())
    croot = project_root() / "cruises"
    for c in sorted(p.name for p in croot.iterdir() if p.is_dir()) if croot.is_dir() else []:
        print(f"  {c:12s} stored as {cruise_layout(c)}")
