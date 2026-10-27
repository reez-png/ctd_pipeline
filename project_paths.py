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

import os
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
