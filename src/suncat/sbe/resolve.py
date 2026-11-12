"""Locate the Sea-Bird executable and PSA setup file for a processing module (optional backend).

Extracted verbatim from step01. The two directory roots (the SBE install folder and the PSA folder)
were module globals there; here they are explicit parameters, so the functions are self-contained and
testable with temporary folders. step01 keeps thin wrappers that pass SBE_BIN_DIR and PSA_ROOT.

module_output_file is pure: it only formats a filename from the module record.
"""
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict


def _require_dir(path: Path, label: str) -> None:
    """Raise a clear error if a required folder is missing (step01's require_folder)."""
    if not path.exists():
        raise FileNotFoundError(f"{label} folder does not exist:\n{path}")
    if not path.is_dir():
        raise NotADirectoryError(f"{label} path exists but is not a folder:\n{path}")


def find_executable(exe_name: str, sbe_bin: Path) -> Path:
    """Find a Sea-Bird executable inside the SBE Data Processing folder."""
    _require_dir(sbe_bin, "SBE executable")
    direct = sbe_bin / exe_name
    if direct.exists():
        return direct
    matches = list(sbe_bin.rglob(exe_name))
    if matches:
        return matches[0]
    raise FileNotFoundError(
        f"Could not find SeaBird executable:\n{exe_name}\n\n"
        f"Expected under:\n{sbe_bin}\n\n"
        "Check your SBE Data Processing installation path."
    )


def resolve_psa(module: Dict[str, Any], psa_root: Path) -> Path:
    """Resolve the best PSA file for a module using its candidate names."""
    _require_dir(psa_root, "PSA setup")

    candidates = list(module.get("psa_candidates", []))
    canonical = module.get("psa_canonical") or module.get("psa") or module.get("psa_file") or ""
    canonical = Path(canonical).name if canonical else ""
    if canonical and canonical not in candidates:
        candidates.insert(0, canonical)

    candidates = [Path(name).name for name in candidates if str(name).strip()]

    # First try exact names in psa_root.
    for name in candidates:
        direct = psa_root / name
        if direct.exists():
            return direct

    # Then try recursive exact name matches.
    for name in candidates:
        matches = sorted(psa_root.rglob(name))
        if matches:
            return matches[0]

    # Final fallback: look for the logical module token anywhere in a .psa filename.
    token = module["name"].split("_", 1)[-1].lower()
    all_psa = sorted(psa_root.rglob("*.psa"))
    loose = [p for p in all_psa if token in p.name.lower()]
    if loose:
        return loose[0]

    raise FileNotFoundError(
        f"Could not find a PSA setup file for module {module['name']}.\n\n"
        f"Looked under:\n{psa_root}\n\n"
        "Accepted candidate names were:\n"
        + "\n".join(f"  - {name}" for name in candidates)
    )


def module_output_file(module: Dict[str, Any], cast_id: str) -> str:
    """Build the expected output filename for one module and cast. Pure."""
    if "output_file" in module and module["output_file"]:
        return str(module["output_file"]).format(cast_id=cast_id)
    output_suffix = str(module.get("output_suffix", ""))
    if not output_suffix:
        raise ValueError(f"Module {module.get('name', '')} has no output_file or output_suffix.")
    return f"{cast_id}{output_suffix}"
