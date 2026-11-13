"""Per-cruise SUNA calibration: which Stage 3 calibration (Zheng et al. 2024, N_QC = (N - a0) / a1) a cruise uses.

The choice lives in one small file per cruise, cruises/<CRUISE>/metadata/suna_calibration.json, so it travels with
the cruise's data and cannot be changed by a run command or a notebook default:

    {"mode": "fit"}
        a1, a0 are fitted in step03 from this cruise's own bottle match-ups (P45_06).
    {"mode": "adopted", "a1": 0.89, "a0": -0.39, "source": "..."}
        a1, a0 are taken as given; `source` says where they come from (P45_07: the homogeneous-layer fit, because
        the response lag corrupts depth-matched pairs).

A missing or malformed file is an error, not a default: a cruise must state its calibration before nitrate is
produced. Write the file with tools/set_suna_calibration.py.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

FILE_NAME = "suna_calibration.json"
MODES = ("fit", "adopted")


class SunaCalibrationError(ValueError):
    """The calibration file is missing or does not say unambiguously what to apply."""


def calibration_path(cruise_dir) -> Path:
    return Path(cruise_dir) / "metadata" / FILE_NAME


def validate(raw: dict, where: str = "suna_calibration") -> dict:
    """Check a calibration record and return it normalised: mode, a1, a0 (None for "fit"), source, describe."""
    if not isinstance(raw, dict):
        raise SunaCalibrationError(f"{where}: expected a JSON object, got {type(raw).__name__}")
    mode = str(raw.get("mode", "")).strip().lower()
    if mode not in MODES:
        raise SunaCalibrationError(f"{where}: 'mode' must be one of {MODES}, got {raw.get('mode')!r}")
    source = str(raw.get("source") or "").strip()
    if mode == "fit":
        extra = [k for k in ("a1", "a0") if raw.get(k) is not None]
        if extra:
            raise SunaCalibrationError(f"{where}: mode 'fit' fits a1 and a0 from the bottles; remove {extra} "
                                       "(or use mode 'adopted')")
        return dict(mode="fit", a1=None, a0=None, source=source or "fitted in step03 from this cruise's bottles",
                    describe="FIT from this cruise's bottle match-ups")
    vals = {}
    for k in ("a1", "a0"):
        try:
            vals[k] = float(raw.get(k))
        except (TypeError, ValueError):
            raise SunaCalibrationError(f"{where}: mode 'adopted' needs a numeric {k}, got {raw.get(k)!r}") from None
        if not math.isfinite(vals[k]):
            raise SunaCalibrationError(f"{where}: {k} must be finite, got {raw.get(k)!r}")
    if vals["a1"] <= 0:
        raise SunaCalibrationError(f"{where}: gain a1 must be positive (N_QC = (N - a0) / a1), got {vals['a1']}")
    if not source:
        raise SunaCalibrationError(f"{where}: mode 'adopted' needs a 'source' saying where a1 and a0 come from")
    return dict(mode="adopted", a1=vals["a1"], a0=vals["a0"], source=source,
                describe=f"ADOPTED a1={vals['a1']:g}, a0={vals['a0']:g} ({source})")


def load_suna_calibration(cruise_dir) -> dict:
    """The cruise's calibration record (validate's keys plus `path`). Raises FileNotFoundError if the file is
    missing and SunaCalibrationError if it is unreadable or ambiguous."""
    path = calibration_path(cruise_dir)
    cruise = Path(cruise_dir).name
    if not path.exists():
        raise FileNotFoundError(
            f"No SUNA calibration for {cruise}: {path} is missing. State it before producing nitrate, e.g.\n"
            f"  python tools/set_suna_calibration.py {cruise} fit\n"
            f"  python tools/set_suna_calibration.py {cruise} adopted --a1 0.89 --a0 -0.39 --source \"...\"")
    try:
        raw = json.loads(path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as e:
        raise SunaCalibrationError(f"{path}: not valid JSON ({e})") from None
    rec = validate(raw, where=str(path))
    rec["path"] = path
    return rec


def write_suna_calibration(cruise_dir, mode: str, a1: float | None = None, a0: float | None = None,
                           source: str | None = None, overwrite: bool = False) -> Path:
    """Validate and write the cruise's calibration file. Refuses to replace an existing file unless overwrite."""
    raw = {"mode": mode}
    if a1 is not None:
        raw["a1"] = a1
    if a0 is not None:
        raw["a0"] = a0
    if source:
        raw["source"] = source
    validate(raw)
    path = calibration_path(cruise_dir)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists:\n  {path.read_text(encoding='utf-8-sig').strip()}\n"
                              "use --force to replace it")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(raw, indent=2) + "\n", encoding="utf-8")
    return path
