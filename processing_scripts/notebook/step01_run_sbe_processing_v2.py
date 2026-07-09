"""
step01_run_sbe_processing_v2.py

Script version of the step01 notebook, for PyCharm Community Edition.
Runs the official Sea-Bird SBE Data Processing modules to produce L1 + L2 products.

HOW TO RUN
  1. Edit the USER SETTINGS section below (CRUISE_ID, paths, RUN_SBE_COMMANDS).
  2. Keep RUN_SBE_COMMANDS = False for a dry run first (writes the command plan only,
     launches no Sea-Bird binary).
  3. Right-click in the editor -> Run, or use the green play button.
  4. Inspect printed tables and the files in <cruise>/_audit, then set RUN_SBE_COMMANDS = True
     and run again to actually process.

Flat top-to-bottom script: settings and helpers are at module level, identical scoping to the notebook.
"""


from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd




# ==========================================================================
# Step 1. Run official Sea-Bird SBE Data Processing from Python (v2)
# ==========================================================================


# ==========================================================================
# Shared XMLCON mode for the April 2026 casts
# ==========================================================================


# ==========================================================================
# 1. User settings
# ==========================================================================

from pathlib import Path

# ===========================================================================
# 1. Project + cruise paths  (v2 = cruise-centric L0/L1/L2/L3)
# ===========================================================================
CTD_ROOT     = Path(r"C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd")
PSA_ROOT     = CTD_ROOT / "psa"
CRUISES_ROOT = CTD_ROOT / "cruises"

# The cruise being processed. 05 = December 2025, 06 = April 2026.
CRUISE_ID  = "P45_06"
CRUISE_DIR = CRUISES_ROOT / CRUISE_ID

# --- Level roots for this cruise -------------------------------------------
L0_CTD_ROOT  = CRUISE_DIR / "L0" / "CTD"     # renamed *.hex live directly here
L0_SUNA_ROOT = CRUISE_DIR / "L0" / "SUNA"    # *.bin
L1_CTD_ROOT  = CRUISE_DIR / "L1" / "CTD"     # readable full-cast *.cnv, NO binning
L1_SUNA_ROOT = CRUISE_DIR / "L1" / "SUNA"    # *.csv from SeaBird UCI
L2_CTD_ROOT  = CRUISE_DIR / "L2" / "CTD"
L2_SUNA_ROOT = CRUISE_DIR / "L2" / "SUNA"    # *_SUNA_1s.csv (merge, step02)
L3_ROOT      = CRUISE_DIR / "L3"             # SUNA nitrate QC product (reserved)

# --- L2 process subfolders (names match supervisor slides 6 and 7) ---------
L2_ALIGN   = L2_CTD_ROOT / "Align_CTD"
L2_LOOP    = L2_CTD_ROOT / "Loop_Edit"
L2_DERIVED = L2_CTD_ROOT / "Derived_Parameter"   # unbinned hub
L2_1M_DOWN = L2_CTD_ROOT / "CTD_1m_down"         # science profiles
L2_1S      = L2_CTD_ROOT / "CTD_1s"              # for SUNA merge

# Raw input + audit output for this run.
RAW_INPUT_ROOT = L0_CTD_ROOT
OUTPUT_ROOT    = CRUISE_DIR / "_audit"           # per-cruise audit (was global outputs/)

# Calibration: per-cruise calibration folder (spec slide 3).
CALIBRATION_ROOT     = CRUISE_DIR / "Instrument_Calibration_Files"
USER_SUPPLIED_XMLCON = CALIBRATION_ROOT / "19-8460.xmlcon"

# ===========================================================================
# 2. Sea-Bird installation + no-space working folder
# ===========================================================================
SBE_BIN_DIR = Path(r"C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32")
# 64-bit install alternative:
# SBE_BIN_DIR = Path(r"C:\Program Files\Sea-Bird\SBEDataProcessing-Win32")

# Short, space-free staging folder. Kept inside the project in v2 so the whole
# pipeline is self-contained; still has no spaces, which is what SBE needs.
SBE_WORK_ROOT        = CTD_ROOT / "_sbe_work"
SBE_WORK_PSA_ROOT    = SBE_WORK_ROOT / "sbe_psa"
SBE_WORK_RAW_ROOT    = SBE_WORK_ROOT / "raw"
SBE_WORK_OUTPUT_ROOT = SBE_WORK_ROOT / "sbe_outputs"

# ===========================================================================
# 3. File patterns
# ===========================================================================
RAW_FILE_PATTERNS = ["*.hex", "*.HEX"]
CONFIG_FILE_PATTERNS = ["*.xmlcon", "*.XMLCON", "*.con", "*.CON"]

# ===========================================================================
# 4. Config-matching behaviour (one shared xmlcon for the batch)
# ===========================================================================
USE_ORIGINAL_USER_XMLCON_FOR_COMMAND = True
SEARCH_WHOLE_RAW_TREE_FOR_CONFIG = False
BORROW_MISSING_XMLCON = False
CREATE_CAST_CONFIG_ALIAS = False
XMLCON_VALIDATION_REQUIRED = True

# ===========================================================================
# 5. Gated optional modules
#    Default OFF so the canonical run matches the 5 PSAs shipped in psa/.
#    See decisions memo before turning these on.
# ===========================================================================
APPLY_WILDEDIT = False   # Despike (Sea-Bird manual). Off pending despiking decision.
APPLY_CELLTM   = False   # Cell Thermal Mass. Off pending "Question for Drew" (slide 6).
APPLY_FILTER   = False   # Low-pass pressure filter. Off; not in v2 PSA set.

# ===========================================================================
# 6. Run-mode safety switches
# ===========================================================================
TEST_SINGLE_CAST_ONLY = False
TEST_CAST_ID = "P45_06_CTD_01"

RUN_SBE_COMMANDS = True          # dry run first; flip to True after log looks right
STOP_ON_ERROR = True
CLEAN_WORK_FOLDER_FIRST = True

# ASCII export is optional in v2. The SUNA merge (step02) reads the 1 s CNV
# directly, so ASCII Out is only needed if you also want semicolon .asc copies.
RUN_ASCII_OUT = False

# Run modules from the input file's own folder using local (basename) args.
RUN_FROM_INPUT_FILE_FOLDER = True
USE_LOCAL_INPUT_AND_CONFIG_NAMES = True
COPY_INPUT_FILES_TO_CAST_FOLDER = True

print("CRUISE_ID           :", CRUISE_ID)
print("RAW_INPUT_ROOT (L0) :", RAW_INPUT_ROOT)
print("L1 readable CNV     :", L1_CTD_ROOT)
print("L2 derived (hub)    :", L2_DERIVED)
print("L2 1 m down         :", L2_1M_DOWN)
print("L2 1 s full         :", L2_1S)
print("Audit               :", OUTPUT_ROOT)
print("APPLY_WILDEDIT / CELLTM / FILTER:", APPLY_WILDEDIT, APPLY_CELLTM, APPLY_FILTER)


# ==========================================================================
# 2. Define the Sea-Bird module sequence and branch products
# ==========================================================================

# ===========================================================================
# LINEAR SEQUENCE
#   DatCnv writes the readable full-cast CNV to L1/CTD (no binning).
#   Every later module writes into an L2/CTD/<process>/ folder.
#   WildEdit, CellTM and Filter are GATED (default off) so the canonical run
#   matches the five PSAs in psa/. Their output_folder is L2_ALIGN so that,
#   when enabled, they stay in the alignment stage and feed LoopEdit cleanly.
# ===========================================================================
LINEAR_MODULE_SEQUENCE = [
    {
        "name": "01_datcnv",
        "description": "Data Conversion from raw Sea-Bird hex to engineering/oceanographic units (L1 readable full cast)",
        "exe": "DatCnvW.exe",
        "psa": "01_datcnv.psa",
        "enabled": True,
        "input_kind": "raw",
        "output_folder": L1_CTD_ROOT,                 # <-- L1, not OUTPUT_ROOT
        "output_suffix": ".cnv",                      # readable full cast, plain stem
        "needs_config": True,
        "extra_args": [],
    },
    {
        "name": "02_alignctd",
        "description": "Align CTD sensor timing (reduces salinity/oxygen spiking)",
        "exe": "AlignCTDW.exe",
        "psa": "02_alignctd.psa",
        "enabled": True,
        "input_kind": "previous",
        "output_folder": L2_ALIGN,                    # <-- L2/CTD/Align_CTD
        "output_suffix": "_al.cnv",
        "needs_config": False,
        "extra_args": [],
    },
    {
        "name": "03_wildedit",
        "description": "Despike obvious outliers (GATED: Sea-Bird manual recommends before LoopEdit)",
        "exe": "WildEditW.exe",
        "psa": "03_wildedit.psa",
        "enabled": APPLY_WILDEDIT,                    # default False
        "input_kind": "previous",
        "output_folder": L2_ALIGN,
        "output_suffix": "_al_we.cnv",
        "needs_config": False,
        "extra_args": [],
    },
    {
        "name": "04_celltm",
        "description": "Correct conductivity cell thermal mass (GATED: 'Question for Drew', slide 6)",
        "exe": "CellTMW.exe",
        "psa": "04_celltm.psa",
        "enabled": APPLY_CELLTM,                      # default False
        "input_kind": "previous",
        "output_folder": L2_ALIGN,
        "output_suffix": "_al_ctm.cnv",
        "needs_config": False,
        "extra_args": [],
    },
    {
        "name": "05_filter",
        "description": "Low-pass filter pressure before LoopEdit (GATED: not in v2 PSA set)",
        "exe": "FilterW.exe",
        "psa": "05_filter.psa",
        "enabled": APPLY_FILTER,                      # default False
        "input_kind": "previous",
        "output_folder": L2_ALIGN,
        "output_suffix": "_al_filt.cnv",
        "needs_config": False,
        "extra_args": [],
    },
    {
        "name": "06_loopedit",
        "description": "Mark pressure reversals and slow movement loops (ship heave)",
        "exe": "LoopEditW.exe",
        "psa": "06_loopedit.psa",
        "enabled": True,
        "input_kind": "previous",
        "output_folder": L2_LOOP,                     # <-- L2/CTD/Loop_Edit
        "output_suffix": "_loop.cnv",
        "needs_config": False,
        "extra_args": [],
    },
    {
        "name": "07_derive",
        "description": "Derive oceanographic variables (salinity, density, etc.) - unbinned hub",
        "exe": "DeriveW.exe",
        "psa": "07_derive.psa",
        "enabled": True,
        "input_kind": "previous",
        "output_folder": L2_DERIVED,                  # <-- L2/CTD/Derived_Parameter
        "output_suffix": "_der.cnv",
        "needs_config": False,
        "extra_args": [],
    },
]

# ===========================================================================
# BIN AVERAGE BRANCHES
#   Both read the SAME derived unbinned CNV. Canonical v2 set = 2 products.
# ===========================================================================
BINAVG_BRANCHES = [
    {
        "name": "08b_binavg_1m_down",
        "description": "1 m depth-bin average, downcast only (science profiles, slide 6)",
        "exe": "BinAvgW.exe",
        "psa": "08b_binavg_1m_down.psa",
        "enabled": True,
        "input_kind": "derived_unbinned",
        "output_folder": L2_1M_DOWN,                  # <-- L2/CTD/CTD_1m_down
        "output_suffix": "_1m_down.cnv",
        "needs_config": False,
        "extra_args": [],
    },
    {
        "name": "08d_binavg_1s_full",
        "description": "1 s time-bin average, full cast for SUNA matching (slide 7)",
        "exe": "BinAvgW.exe",
        "psa": "08d_binavg_1s_full.psa",
        "enabled": True,
        "input_kind": "derived_unbinned",
        "output_folder": L2_1S,                       # <-- L2/CTD/CTD_1s
        "output_suffix": "_1s.cnv",
        "needs_config": False,
        "extra_args": [],
    },
]

# ===========================================================================
# ASCII BRANCHES (optional; only if RUN_ASCII_OUT = True)
#   Semicolon-separated, header labels on, time conversion on.
# ===========================================================================
ASCII_BRANCHES = [
    {
        "name": "09b_asciiout_1m_down",
        "description": "ASCII export of 1 m downcast CNV",
        "exe": "ASCII_OutW.exe",
        "psa": "09_asciiout_semicolon.psa",
        "psa_candidates": ["09b_asciiout_1m_down.psa", "09_asciiout_semicolon.psa", "09_asciiout.psa"],
        "enabled": RUN_ASCII_OUT,
        "input_folder": L2_1M_DOWN,
        "output_folder": L2_1M_DOWN / "asc_semicolon",
        "input_suffix": "_1m_down.cnv",
        "output_suffix": "_1m_down_ascii.asc",
        "needs_config": False,
        "extra_args": [],
    },
    {
        "name": "09d_asciiout_1s_full",
        "description": "ASCII export of 1 s full cast CNV for SUNA matching",
        "exe": "ASCII_OutW.exe",
        "psa": "09_asciiout_semicolon.psa",
        "psa_candidates": ["09d_asciiout_1s_full.psa", "09_asciiout_semicolon.psa", "09_asciiout.psa"],
        "enabled": RUN_ASCII_OUT,
        "input_folder": L2_1S,
        "output_folder": L2_1S / "asc_semicolon",
        "input_suffix": "_1s.cnv",
        "output_suffix": "_1s_ascii.asc",
        "needs_config": False,
        "extra_args": [],
    },
]

# MODULE_SEQUENCE is kept for validation, PSA copying, summaries and reports.
MODULE_SEQUENCE = LINEAR_MODULE_SEQUENCE + BINAVG_BRANCHES + ASCII_BRANCHES

FINAL_LINEAR_MODULE_NAME = "07_derive"

# ## 2B. Normalize Sea-Bird module records
#
# This cell makes the module dictionaries compatible with the no-space working
# folder functions. It adds `psa_canonical` automatically from the existing
# `psa` value and keeps branch modules in the same validation table.

from pathlib import Path

def normalize_module_record(module: dict) -> dict:
    """
    Return a safe copy of one Sea-Bird module dictionary.
    """
    m = dict(module)

    if "psa_canonical" not in m:
        if "psa" in m:
            m["psa_canonical"] = Path(m["psa"]).name
        elif "psa_file" in m:
            m["psa_canonical"] = Path(m["psa_file"]).name
        else:
            raise KeyError(
                "Each module must have either 'psa', 'psa_file', or 'psa_canonical'. "
                f"Problem module: {m}"
            )

    if "psa_candidates" not in m:
        m["psa_candidates"] = [m["psa_canonical"]]
    elif m["psa_canonical"] not in m["psa_candidates"]:
        m["psa_candidates"] = [m["psa_canonical"], *list(m["psa_candidates"])]

    if "exe_canonical" not in m:
        if "exe" in m:
            m["exe_canonical"] = Path(m["exe"]).name
        elif "exe_file" in m:
            m["exe_canonical"] = Path(m["exe_file"]).name
        else:
            m["exe_canonical"] = ""

    m.setdefault("enabled", True)
    m.setdefault("needs_config", False)
    m.setdefault("extra_args", [])
    return m


def expected_enabled_modules() -> list:
    """
    Return enabled Sea-Bird modules with normalized key names.
    """
    return [
        normalize_module_record(module)
        for module in MODULE_SEQUENCE
        if module.get("enabled", True)
    ]


def enabled_linear_modules() -> list:
    return [
        normalize_module_record(module)
        for module in LINEAR_MODULE_SEQUENCE
        if module.get("enabled", True)
    ]


def enabled_binavg_branches() -> list:
    return [
        normalize_module_record(module)
        for module in BINAVG_BRANCHES
        if module.get("enabled", True)
    ]


def enabled_ascii_branches() -> list:
    return [
        normalize_module_record(module)
        for module in ASCII_BRANCHES
        if module.get("enabled", True)
    ]


module_check = pd.DataFrame(
    [
        {
            "name": m.get("name", ""),
            "exe": m.get("exe", ""),
            "psa_canonical": m.get("psa_canonical", ""),
            "psa_candidates": "; ".join(m.get("psa_candidates", [])),
            "output_folder": str(m.get("output_folder", "")),
            "output_suffix": m.get("output_suffix", ""),
            "input_suffix": m.get("input_suffix", ""),
            "enabled": m.get("enabled", True),
        }
        for m in expected_enabled_modules()
    ]
)

print(module_check)

assert "psa_canonical" in expected_enabled_modules()[0], "psa_canonical was not created."
print("Module records normalized successfully.")


# ==========================================================================
# 3. Helper functions
# ==========================================================================


def safe_name(value: str) -> str:
    """Make a safe file or folder name."""
    value = str(value)
    for bad in ['<', '>', ':', '"', '/', '\\', '|', '?', '*']:
        value = value.replace(bad, "_")
    value = "_".join(value.split()).strip("_")
    return value or "unnamed"


def require_folder(path: Path, label: str) -> None:
    """Raise a clear error if a required folder is missing."""
    if not path.exists():
        raise FileNotFoundError(f"{label} folder does not exist:\n{path}")
    if not path.is_dir():
        raise NotADirectoryError(f"{label} path exists but is not a folder:\n{path}")


def ensure_folder(path: Path) -> Path:
    """Create a folder if needed and return it."""
    path.mkdir(parents=True, exist_ok=True)
    return path


def ensure_no_spaces(path: Path, label: str) -> None:
    """SeaBird command-line path arguments should not contain spaces."""
    if " " in str(path):
        raise ValueError(
            f"{label} contains a space:\n  {path}\n\n"
            "SeaBird command-line modules often fail with quoted paths. "
            "Use the short working folder under C:\\sbe_work for command arguments."
        )


def find_executable(exe_name: str) -> Path:
    """Find a SeaBird executable inside the SBE Data Processing folder."""
    require_folder(SBE_BIN_DIR, "SBE executable")
    direct = SBE_BIN_DIR / exe_name
    if direct.exists():
        return direct
    matches = list(SBE_BIN_DIR.rglob(exe_name))
    if matches:
        return matches[0]
    raise FileNotFoundError(
        f"Could not find SeaBird executable:\n{exe_name}\n\n"
        f"Expected under:\n{SBE_BIN_DIR}\n\n"
        "Check your SBE Data Processing installation path."
    )


def resolve_psa_for_module(module: Dict[str, Any]) -> Path:
    """Resolve the best PSA file for a module using the candidate names."""
    require_folder(PSA_ROOT, "PSA setup")

    candidates = list(module.get("psa_candidates", []))
    canonical = module.get("psa_canonical") or module.get("psa") or module.get("psa_file") or ""
    canonical = Path(canonical).name if canonical else ""
    if canonical and canonical not in candidates:
        candidates.insert(0, canonical)

    candidates = [Path(name).name for name in candidates if str(name).strip()]

    # First try exact names in PSA_ROOT.
    for name in candidates:
        direct = PSA_ROOT / name
        if direct.exists():
            return direct

    # Then try recursive exact name matches.
    for name in candidates:
        matches = sorted(PSA_ROOT.rglob(name))
        if matches:
            return matches[0]

    # Final fallback: look for the logical module token anywhere in a .psa filename.
    token = module["name"].split("_", 1)[-1].lower()
    all_psa = sorted(PSA_ROOT.rglob("*.psa"))
    loose = [p for p in all_psa if token in p.name.lower()]
    if loose:
        return loose[0]

    raise FileNotFoundError(
        f"Could not find a PSA setup file for module {module['name']}.\n\n"
        f"Looked under:\n{PSA_ROOT}\n\n"
        "Accepted candidate names were:\n"
        + "\n".join(f"  - {name}" for name in candidates)
    )


# expected_enabled_modules() is defined in Section 2B. Do not redefine it here.

def module_output_file_for_cast(module: Dict[str, Any], cast_id: str) -> str:
    """Build the expected output filename for one module and cast."""
    if "output_file" in module and module["output_file"]:
        return str(module["output_file"]).format(cast_id=cast_id)
    output_suffix = str(module.get("output_suffix", ""))
    if not output_suffix:
        raise ValueError(f"Module {module.get('name', '')} has no output_file or output_suffix.")
    return f"{cast_id}{output_suffix}"


def validate_module_sequence(module_sequence: List[Dict[str, Any]]) -> pd.DataFrame:
    """Validate executables, PSA files and output names before processing."""
    records = []

    for step_number, module in enumerate(module_sequence, start=1):
        module = normalize_module_record(module)
        name = module.get("name", "")
        enabled = bool(module.get("enabled", True))
        exe_name = module.get("exe", "")
        canonical_psa = module.get("psa_canonical", "")
        output_folder = module.get("output_folder", "")
        output_suffix = module.get("output_suffix", "")
        input_suffix = module.get("input_suffix", "")

        if not enabled:
            records.append({
                "step": step_number,
                "name": name,
                "enabled": False,
                "exe": exe_name,
                "exe_found": None,
                "exe_path": "",
                "canonical_psa": canonical_psa,
                "selected_psa_found": None,
                "selected_psa_path": "",
                "output_folder": str(output_folder),
                "output_suffix": output_suffix,
                "input_suffix": input_suffix,
                "status": "SKIPPED_DISABLED",
            })
            continue

        status_parts = []
        exe_path: Optional[Path] = None
        selected_psa: Optional[Path] = None

        try:
            exe_path = find_executable(exe_name)
        except Exception:
            status_parts.append("MISSING_EXECUTABLE")

        try:
            selected_psa = resolve_psa_for_module(module)
        except Exception:
            status_parts.append("MISSING_PSA")

        if not output_folder:
            status_parts.append("MISSING_OUTPUT_FOLDER")
        if not output_suffix:
            status_parts.append("MISSING_OUTPUT_SUFFIX")

        status = "OK" if not status_parts else " | ".join(status_parts)

        records.append({
            "step": step_number,
            "name": name,
            "enabled": True,
            "exe": exe_name,
            "exe_found": exe_path is not None,
            "exe_path": str(exe_path) if exe_path else "",
            "canonical_psa": canonical_psa,
            "selected_psa_found": selected_psa is not None,
            "selected_psa_path": str(selected_psa) if selected_psa else "",
            "output_folder": str(output_folder),
            "output_suffix": output_suffix,
            "input_suffix": input_suffix,
            "status": status,
        })

    validation_df = pd.DataFrame(records)
    problems = validation_df[(validation_df["enabled"] == True) & (validation_df["status"] != "OK")]

    if not problems.empty:
        print(validation_df)
        raise FileNotFoundError(
            "One or more Sea-Bird processing requirements are missing.\n\n"
            "Check the validation table above before running the batch processor."
        )

    return validation_df



def validate_processing_setup() -> pd.DataFrame:
    """Run all setup validation checks before batch processing."""
    require_folder(RAW_INPUT_ROOT, "Raw input")
    require_folder(PSA_ROOT, "PSA setup")
    require_folder(SBE_BIN_DIR, "SBE executable")
    ensure_folder(OUTPUT_ROOT)
    ensure_folder(L1_CTD_ROOT)
    ensure_folder(L2_CTD_ROOT)
    ensure_folder(L2_DERIVED)

    for module in expected_enabled_modules():
        output_folder = module.get("output_folder")
        if output_folder:
            ensure_folder(Path(output_folder))

    return validate_module_sequence(MODULE_SEQUENCE)


def find_files_by_patterns(root: Path, patterns: List[str]) -> List[Path]:
    """Find files matching any pattern inside a folder."""
    require_folder(root, "Search root")
    files: List[Path] = []
    for pattern in patterns:
        files.extend(root.rglob(pattern))
    return sorted(set(files))




def is_probably_valid_seabird_config(config_path: Path) -> Tuple[bool, str]:
    """
    Lightweight validation for SeaBird .xmlcon or .con files.

    This does not replace DatCnv's own parser, but it catches common problems before
    DatCnvW opens a blocking popup, such as:
      - wrong file copied with .xmlcon extension;
      - OneDrive placeholder or failed download text;
      - empty or tiny file;
      - plain XML that is not an instrument configuration.
    """
    try:
        if not config_path.exists():
            return False, "file does not exist"
        if not config_path.is_file():
            return False, "path is not a file"

        size = config_path.stat().st_size
        if size < 200:
            return False, f"file is too small to be a normal instrument configuration ({size} bytes)"

        raw = config_path.read_bytes()[:20000]
        text = raw.decode("utf-8", errors="ignore")
        lower = text.lower()

        if "<html" in lower or "<!doctype html" in lower:
            return False, "file looks like HTML, not a SeaBird instrument configuration"
        if "onedrive" in lower and "instrument" not in lower:
            return False, "file may be a cloud placeholder or sync message, not an instrument configuration"
        if "not found" in lower and "instrument" not in lower:
            return False, "file content looks like an error message"

        # SeaBird XMLCON files usually contain instrument/configuration/sensor/calibration tags.
        xmlcon_tokens = [
            "sbe_instrumentconfiguration",
            "instrumentconfiguration",
            "instrument configuration",
            "sensorarray",
            "calibrationcoefficients",
            "seabird",
            "sea-bird",
            "sbe 19",
            "sbe19",
            "sbe 19plus",
            "sbe 25",
            "sbe25",
            "sbe 9",
            "sbe9",
        ]

        # Legacy .con files are text-like and often include instrument and sensor/calibration information.
        con_tokens = [
            "configuration report",
            "instrument type",
            "serial number",
            "temperature sensor",
            "conductivity sensor",
            "pressure sensor",
            "oxygen sensor",
            "calibration",
            "sbe",
            "sea-bird",
        ]

        tokens = xmlcon_tokens + con_tokens
        matched = [tok for tok in tokens if tok in lower]
        if not matched:
            preview = " ".join(text.split()[:40])
            return False, f"no SeaBird configuration markers found. Preview: {preview[:250]}"

        return True, f"passed lightweight validation; matched markers: {', '.join(matched[:5])}"

    except Exception as exc:
        return False, f"could not read config file: {exc!r}"


def validate_or_raise_config(config_path: Path, context: str = "") -> Path:
    """Return config_path if valid enough, otherwise raise a useful error."""
    if not XMLCON_VALIDATION_REQUIRED:
        return config_path

    ok, reason = is_probably_valid_seabird_config(config_path)
    if ok:
        return config_path

    try:
        preview = config_path.read_text(encoding="utf-8", errors="replace")[:1200]
    except Exception:
        preview = "(could not read text preview)"

    message = f"""The selected SeaBird configuration file does not look valid.

Context: {context}
Path: {config_path}
Reason: {reason}

This matches the DatCnv popup pattern: the file exists, but DatCnv does not accept it as a real instrument configuration file.

Fix:
  1. Use the exact .xmlcon or .con exported from SeaBird for this CTD package.
  2. Make sure the file is fully downloaded from OneDrive.
  3. Set USER_SUPPLIED_XMLCON to that valid file if automatic matching picks the wrong one.

Beginning of selected file:
{preview}
"""
    raise ValueError(message)


def auto_find_xmlcon(user_path: Optional[str], hex_path: Path) -> Optional[Path]:
    """
    Find the best SeaBird configuration file for a raw .hex file.

    Safer v10 behavior:
    1. If USER_SUPPLIED_XMLCON or user_path is provided, use only that file.
    2. Prefer same-stem config beside the raw file.
    3. Prefer the only valid config in the raw file folder.
    4. Search the whole raw tree only when SEARCH_WHOLE_RAW_TREE_FOR_CONFIG is True.
    5. Reject files that do not look like real SeaBird .xmlcon or .con files.
    """
    manual = user_path
    if manual is None and "USER_SUPPLIED_XMLCON" in globals() and USER_SUPPLIED_XMLCON is not None:
        manual = str(USER_SUPPLIED_XMLCON)

    if manual:
        p = Path(manual).expanduser()
        if p.exists() and p.suffix.lower() in {".xmlcon", ".con"}:
            return validate_or_raise_config(p.resolve(), context=f"manual config for {hex_path.name}")
        return None

    candidates: List[Path] = []

    # Strongest match: exact same path stem beside the raw file.
    candidates.extend([
        hex_path.with_suffix(".xmlcon"),
        hex_path.with_suffix(".XMLCON"),
        hex_path.with_suffix(".con"),
        hex_path.with_suffix(".CON"),
    ])

    # Same folder configs only. This avoids accidentally taking a bad config from
    # a different cast folder.
    for pattern in ("*.xmlcon", "*.XMLCON", "*.con", "*.CON"):
        candidates.extend(sorted(hex_path.parent.glob(pattern)))

    # Whole raw tree fallback, disabled by default for safety.
    if globals().get("SEARCH_WHOLE_RAW_TREE_FOR_CONFIG", False):
        if "RAW_INPUT_ROOT" in globals() and RAW_INPUT_ROOT.exists():
            for pattern in ("*.xmlcon", "*.XMLCON", "*.con", "*.CON"):
                candidates.extend(sorted(RAW_INPUT_ROOT.rglob(pattern)))

    # Deduplicate while preserving order.
    seen = set()
    unique_candidates = []
    for p in candidates:
        key = str(p).lower()
        if key not in seen:
            seen.add(key)
            unique_candidates.append(p)

    existing = [p for p in unique_candidates if p.exists() and p.suffix.lower() in {".xmlcon", ".con"}]
    if not existing:
        return None

    raw_stem = hex_path.stem.lower()

    valid_candidates: List[Tuple[Path, str]] = []
    invalid_records: List[Dict[str, str]] = []

    for p in existing:
        ok, reason = is_probably_valid_seabird_config(p)
        if ok or not XMLCON_VALIDATION_REQUIRED:
            valid_candidates.append((p.resolve(), reason))
        else:
            invalid_records.append({
                "raw_file": str(hex_path),
                "candidate_config": str(p),
                "reason": reason,
            })

    if invalid_records:
        try:
            pd.DataFrame(invalid_records).to_csv(
                OUTPUT_ROOT / "01_rejected_invalid_config_candidates.csv",
                index=False,
                encoding="utf-8-sig",
            )
        except Exception:
            pass

    if not valid_candidates:
        # Raise on exact same-stem candidate because this is usually what caused
        # the DatCnv popup and should be corrected explicitly.
        exact = [p for p in existing if p.stem.lower() == raw_stem]
        if exact:
            validate_or_raise_config(exact[0].resolve(), context=f"automatic exact-stem match for {hex_path.name}")
        return None

    def score(item: Tuple[Path, str]) -> tuple:
        p, _reason = item
        suffix_score = 0 if p.suffix.lower() == ".xmlcon" else 1
        same_folder_score = 0 if p.parent == hex_path.parent else 1
        same_stem_score = 0 if p.stem.lower() == raw_stem else 1
        loose_stem_score = 0 if (raw_stem in p.stem.lower() or p.stem.lower() in raw_stem) else 1
        return (same_folder_score, same_stem_score, loose_stem_score, suffix_score, len(str(p)))

    return sorted(valid_candidates, key=score)[0][0]

def create_cast_config_alias(source_config: Path, target_raw: Path) -> Path:
    """
    Copy the selected config beside the working raw file using several names.

    This is deliberately redundant because different SeaBird PSA files may remember
    either the original config basename or the cast-specific basename.

    Example:
        C:/sbe_work/raw/J1_15_04_25/J1_15_04_25.hex
        C:/sbe_work/raw/J1_15_04_25/J1_15_04_25.xmlcon
        C:/sbe_work/raw/J1_15_04_25/J1_15_04_25.XMLCON
        C:/sbe_work/raw/J1_15_04_25/J1_15_04_25.con
        C:/sbe_work/raw/J1_15_04_25/J1_15_04_25.CON
        plus the original source filename.
    """
    source_config = validate_or_raise_config(source_config, context=f"before copying config for {target_raw.name}")

    suffix = source_config.suffix.lower()
    if suffix not in {".xmlcon", ".con"}:
        raise ValueError(f"Not a valid SeaBird configuration file suffix: {source_config}")

    target_raw.parent.mkdir(parents=True, exist_ok=True)

    # Preferred config passed to DatCnv. Keep the original suffix.
    alias_config = target_raw.with_suffix(suffix)
    shutil.copy2(source_config, alias_config)
    validate_or_raise_config(alias_config, context=f"after copying alias {alias_config.name}")

    # Also create common extension variants. The content is identical; this handles
    # PSA files that internally ask for .xmlcon, .XMLCON, .con or .CON.
    alias_names = [
        target_raw.with_suffix(".xmlcon"),
        target_raw.with_suffix(".XMLCON"),
        target_raw.with_suffix(".con"),
        target_raw.with_suffix(".CON"),
        target_raw.parent / source_config.name,
    ]

    for alias in alias_names:
        if alias.name.lower() == alias_config.name.lower():
            continue
        shutil.copy2(source_config, alias)

    return alias_config

def choose_config_for_raw(raw_file: Path, config_files: List[Path]) -> Optional[Path]:
    """
    Choose the configuration file for one raw file.

    In the fixed XMLCON workflow, USER_SUPPLIED_XMLCON takes priority for every cast.
    That prevents accidental use of S1_18_04_26.XML, S2_18_04_26.XML or other sidecar XML files.
    """
    found = auto_find_xmlcon(user_path=None, hex_path=raw_file)
    if found is not None:
        return found

    if not config_files:
        return None

    raw_stem = raw_file.stem.lower()

    def priority(path: Path) -> tuple:
        suffix_score = 0 if path.suffix.lower() == ".xmlcon" else 1
        same_folder_score = 0 if path.parent == raw_file.parent else 1
        same_stem_score = 0 if path.stem.lower() == raw_stem else 1
        loose_stem_score = 0 if (raw_stem in path.stem.lower() or path.stem.lower() in raw_stem) else 1
        return (same_folder_score, same_stem_score, loose_stem_score, suffix_score, len(str(path)))

    valid = [c for c in config_files if c.exists() and c.suffix.lower() in {".xmlcon", ".con"}]
    if not valid:
        return None

    return validate_or_raise_config(sorted(valid, key=priority)[0], context=f"fallback match for {raw_file.name}")


def borrow_missing_xmlcon_files(raw_root: Path) -> pd.DataFrame:
    """
    Copy a sister cast's .xmlcon next to each .hex that has no neighbouring .xmlcon or .con.

    This refuses to act if existing .xmlcon files have different SHA-256 hashes.
    That prevents silently using the wrong calibration if the sensor package changed mid-cruise.
    """
    xmlcon_files = sorted(set(raw_root.rglob("*.xmlcon")) | set(raw_root.rglob("*.XMLCON")))

    if not xmlcon_files:
        raise FileNotFoundError(
            f"No .xmlcon files exist anywhere under:\n  {raw_root}\n\n"
            "Cannot borrow without at least one source config. Recreate an .xmlcon first."
        )

    hashes: Dict[str, List[Path]] = {}
    for p in xmlcon_files:
        h = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
        hashes.setdefault(h, []).append(p)

    if len(hashes) != 1:
        msg_lines = [
            f"Found {len(hashes)} distinct .xmlcon contents under {raw_root}.",
            "Refusing to auto-borrow because the sensor package may have changed.",
            "Inspect manually:",
        ]
        for h, paths in hashes.items():
            msg_lines.append(f"  hash {h}: {len(paths)} file(s)")
            for p in paths[:5]:
                msg_lines.append(f"    {p}")
        raise ValueError("\n".join(msg_lines))

    source = xmlcon_files[0]
    hex_files = sorted(set(raw_root.rglob("*.hex")) | set(raw_root.rglob("*.HEX")))

    records = []
    for hex_path in hex_files:
        sibling_candidates = [
            hex_path.with_suffix(".xmlcon"),
            hex_path.with_suffix(".XMLCON"),
            hex_path.with_suffix(".con"),
            hex_path.with_suffix(".CON"),
        ]
        existing = next((s for s in sibling_candidates if s.exists()), None)

        if existing is not None:
            records.append({
                "hex_file": str(hex_path),
                "action": "skipped_already_has_config",
                "target": str(existing),
                "source": "",
            })
            continue

        target = hex_path.with_suffix(".xmlcon")
        shutil.copy2(source, target)
        records.append({
            "hex_file": str(hex_path),
            "action": "borrowed_xmlcon",
            "target": str(target),
            "source": str(source),
        })

    log = pd.DataFrame(records)
    log.to_csv(OUTPUT_ROOT / "01_xmlcon_borrow_log.csv", index=False, encoding="utf-8-sig")
    return log


def normalise_config_name_for_sbe(config_file: Path) -> str:
    """Return a SeaBird-friendly configuration filename."""
    suffix = config_file.suffix.lower()
    if suffix == ".xml":
        raise ValueError(
            f"Refusing to use a .xml file as a configuration:\n  {config_file}\n\n"
            "Use a real .xmlcon or .con file instead."
        )
    return config_file.name


def prepare_no_space_work_folder(inventory: pd.DataFrame, clean_existing_work_folder: bool) -> pd.DataFrame:
    """Copy PSA and raw files to C:\\sbe_work, while preserving the selected command XMLCON."""
    ensure_no_spaces(SBE_WORK_ROOT, "SBE_WORK_ROOT")

    if clean_existing_work_folder and SBE_WORK_ROOT.exists():
        shutil.rmtree(SBE_WORK_ROOT)

    ensure_folder(SBE_WORK_ROOT)
    ensure_folder(SBE_WORK_PSA_ROOT)
    ensure_folder(SBE_WORK_RAW_ROOT)
    ensure_folder(SBE_WORK_OUTPUT_ROOT)

    # Copy selected PSA files into the work folder using canonical names.
    psa_copy_records = []
    for module in expected_enabled_modules():
        source_psa = resolve_psa_for_module(module)

        canonical_name = (
            module.get("psa_canonical")
            or module.get("psa")
            or module.get("psa_file")
        )

        if canonical_name is None:
            raise KeyError(
                "Could not determine the PSA filename for this module. "
                f"Module record was: {module}"
            )

        canonical_name = Path(canonical_name).name

        target_psa = SBE_WORK_PSA_ROOT / canonical_name
        if source_psa.stat().st_size == 0:
            raise IOError(
                f"Source PSA is zero bytes (likely a OneDrive Files-On-Demand "
                f"placeholder):\n{source_psa}\n\n"
                "Right-click the parent folder in Explorer and choose "
                "'Always keep on this device', then rerun."
            )
        shutil.copy2(source_psa, target_psa)
        psa_copy_records.append({
            "module": module["name"],
            "source_psa": str(source_psa),
            "work_psa": str(target_psa),
        })

    pd.DataFrame(psa_copy_records).to_csv(
        OUTPUT_ROOT / "01_psa_copy_log.csv",
        index=False,
        encoding="utf-8-sig",
    )

    records = []
    for _, row in inventory.iterrows():
        cast_id = safe_name(row["cast_id"])
        source_raw = Path(row["raw_file"])
        source_config = Path(row["config_file"])

        if not source_raw.exists():
            raise FileNotFoundError(f"Raw file does not exist:\n{source_raw}")
        if not source_config.exists():
            raise FileNotFoundError(f"Configuration file does not exist:\n{source_config}")

        source_config = validate_or_raise_config(source_config.resolve(), context=f"selected config for {cast_id}")

        cast_work_raw_dir = ensure_folder(SBE_WORK_RAW_ROOT / cast_id)
        target_raw = cast_work_raw_dir / source_raw.name
        if source_raw.stat().st_size == 0:
            raise IOError(
                f"Source raw file is zero bytes (likely a OneDrive "
                f"Files-On-Demand placeholder):\n{source_raw}\n\n"
                "Right-click the parent folder in Explorer and choose "
                "'Always keep on this device', then rerun."
            )
        shutil.copy2(source_raw, target_raw)

        # We still keep one archive copy of the selected config in C:\\sbe_work, but the command
        # can use the original Calibration folder path when USE_ORIGINAL_USER_XMLCON_FOR_COMMAND is True.
        if CREATE_CAST_CONFIG_ALIAS:
            target_config = create_cast_config_alias(source_config=source_config, target_raw=target_raw)
        else:
            target_config = cast_work_raw_dir / normalise_config_name_for_sbe(source_config)
            shutil.copy2(source_config, target_config)
            validate_or_raise_config(target_config, context=f"working copy for {cast_id}")

        if (
            globals().get("USE_ORIGINAL_USER_XMLCON_FOR_COMMAND", False)
            and globals().get("USER_SUPPLIED_XMLCON", None) is not None
        ):
            command_config = validate_or_raise_config(
                Path(globals()["USER_SUPPLIED_XMLCON"]).resolve(),
                context=f"command config for {cast_id}",
            )
        else:
            command_config = target_config

        records.append({
            "cast_id": cast_id,
            "source_raw_file": str(source_raw),
            "source_config_file": str(source_config),
            "work_raw_file": str(target_raw),
            "work_config_file": str(target_config),
            "command_config_file": str(command_config),
            "config_alias_created": CREATE_CAST_CONFIG_ALIAS,
            "uses_original_user_xmlcon_for_command": bool(
                globals().get("USE_ORIGINAL_USER_XMLCON_FOR_COMMAND", False)
                and globals().get("USER_SUPPLIED_XMLCON", None) is not None
            ),
        })

    work_inventory = pd.DataFrame(records)
    work_inventory.to_csv(
        OUTPUT_ROOT / "01_working_copy_inventory.csv",
        index=False,
        encoding="utf-8-sig",
    )
    return work_inventory


def find_work_psa(module: Dict[str, Any]) -> Path:
    """Find the canonical copied PSA file in the no-space working folder.

    Tolerant of un-normalised module dicts: tries 'psa_canonical' first,
    then falls back to 'psa' or 'psa_file' (using only the basename), so
    we don't crash with KeyError if Section 2B was skipped.
    """
    canonical = (
        module.get("psa_canonical")
        or (Path(module["psa"]).name if module.get("psa") else None)
        or (Path(module["psa_file"]).name if module.get("psa_file") else None)
    )
    if canonical is None:
        raise KeyError(
            "Module record has no PSA filename. Expected one of "
            "'psa_canonical', 'psa', or 'psa_file'. "
            f"Got keys: {sorted(module.keys())}"
        )
    psa_path = SBE_WORK_PSA_ROOT / canonical
    if not psa_path.exists():
        raise FileNotFoundError(
            f"Copied PSA file does not exist:\n{psa_path}\n\n"
            "Run Section 5D again."
        )
    return psa_path


def build_sbe_command(
    exe_path: Path,
    input_file: Path,
    output_dir: Path,
    output_file: str,
    psa_file: Path,
    config_file: Optional[Path] = None,
    extra_args: Optional[List[str]] = None,
) -> List[str]:
    """Build a Sea-Bird command line call using no-space arguments."""
    extra_args = extra_args or []
    ensure_folder(output_dir)

    for label, p in [
        ("input_file", input_file),
        ("output_dir", output_dir),
        ("psa_file", psa_file),
        ("config_file", config_file),
    ]:
        if p is None:
            continue
        ensure_no_spaces(p, label)

    use_local = bool(globals().get("RUN_FROM_INPUT_FILE_FOLDER", False)) and bool(
        globals().get("USE_LOCAL_INPUT_AND_CONFIG_NAMES", True)
    )

    input_arg = input_file.name if use_local else str(input_file)
    config_arg = None
    if config_file is not None:
        # Only use local config name when it is beside the input file.
        if use_local and config_file.parent.resolve() == input_file.parent.resolve():
            config_arg = config_file.name
        else:
            config_arg = str(config_file)

    # SBE command line modules expect /f to be the output basename.
    # Give /f without .cnv or .asc to avoid doubled extensions.
    output_stem = Path(output_file).stem
    command = [
        str(exe_path),
        f"/i{input_arg}",
        f"/p{psa_file}",
        f"/f{output_stem}",
        f"/o{output_dir}",
    ]

    if config_arg is not None:
        command.append(f"/c{config_arg}")

    command.extend(extra_args)
    command.append("/s")
    return command

def command_to_text(command: List[str]) -> str:
    """Convert command list to readable Windows command text."""
    return subprocess.list2cmdline(command)


def run_command(
    command: List[str],
    log_path: Path,
    working_dir: Optional[Path] = None,
    timeout_seconds: Optional[int] = None,
) -> Tuple[int, str, str]:
    """Run one SeaBird command and save stdout and stderr."""
    ensure_folder(log_path.parent)
    started = datetime.now().isoformat(timespec="seconds")

    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            shell=False,
            cwd=str(working_dir) if working_dir else None,
            timeout=timeout_seconds,
        )
        return_code = result.returncode
        stdout_text = result.stdout or ""
        stderr_text = result.stderr or ""
    except subprocess.TimeoutExpired as exc:
        return_code = 124
        stdout_text = exc.stdout or ""
        stderr_text = (exc.stderr or "") + "\nCommand timed out."
    except Exception as exc:
        return_code = 1
        stdout_text = ""
        stderr_text = repr(exc)

    finished = datetime.now().isoformat(timespec="seconds")
    log_text = [
        f"Started: {started}",
        f"Finished: {finished}",
        f"Return code: {return_code}",
        f"Working directory: {working_dir if working_dir else ''}",
        "",
        "Command:",
        command_to_text(command),
        "",
        "STDOUT:",
        stdout_text or "(empty)",
        "",
        "STDERR:",
        stderr_text or "(empty)",
        "",
    ]
    log_path.write_text("\n".join(log_text), encoding="utf-8")
    return return_code, stdout_text, stderr_text


def find_expected_output(
    output_dir: Path,
    output_file: str,
    started_at: Optional[datetime] = None,
) -> Optional[Path]:
    """Locate a module output, tolerant of PSA-driven renaming.

    First tries the exact expected filename. If that is missing but the
    module exited cleanly, scan the output directory for files with the same
    extension written since `started_at`. This catches common SBE cases where
    a PSA append flag or output template changes the final filename.
    """
    expected = output_dir / output_file
    if expected.exists():
        return expected
    if not output_dir.exists():
        return None

    suffix = Path(output_file).suffix.lower()
    if suffix:
        patterns = [f"*{suffix}"]
        if suffix == ".asc":
            patterns.extend(["*.txt", "*.csv"])
    else:
        patterns = ["*"]

    candidates: List[Path] = []
    for pattern in patterns:
        candidates.extend(output_dir.glob(pattern))

    # Remove folders and duplicate paths.
    unique_candidates = []
    seen = set()
    for candidate in candidates:
        key = str(candidate).lower()
        if candidate.is_file() and key not in seen:
            unique_candidates.append(candidate)
            seen.add(key)
    candidates = unique_candidates

    if started_at is not None:
        cutoff = started_at.timestamp() - 2  # 2 second grace for clock skew
        candidates = [p for p in candidates if p.stat().st_mtime >= cutoff]
    if not candidates:
        return None
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]



def copy_with_log(source: Path, destination: Path, records: List[Dict[str, Any]], cast_id: str, product_type: str) -> None:
    """Copy one product and append an audit record."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    status = "planned_not_run"
    message = ""

    if RUN_SBE_COMMANDS:
        if source.exists():
            shutil.copy2(source, destination)
            status = "copied"
        else:
            status = "source_missing"
            message = f"Source file does not exist: {source}"

    records.append({
        "cast_id": cast_id,
        "product_type": product_type,
        "source_path": str(source),
        "destination_path": str(destination),
        "run_requested": RUN_SBE_COMMANDS,
        "status": status,
        "message": message,
    })


def df_to_markdown_safe(df: pd.DataFrame) -> str:
    """Create a markdown table without requiring tabulate."""
    try:
        return df.to_markdown(index=False)
    except Exception:
        return "```text\n" + df.to_csv(index=False) + "\n```"


def write_validation_report(output_root: Path, module_df: pd.DataFrame) -> None:
    """Write validation results to CSV files for audit."""
    ensure_folder(output_root)
    module_df.to_csv(output_root / "00_preflight_module_check.csv", index=False, encoding="utf-8-sig")


def write_markdown_readme(output_root: Path, module_validation_df: Optional[pd.DataFrame] = None) -> None:
    """Write a README for the processing run."""
    ensure_folder(output_root)
    module_table = df_to_markdown_safe(pd.DataFrame(MODULE_SEQUENCE))
    validation_table = df_to_markdown_safe(module_validation_df) if module_validation_df is not None else "Validation table was not supplied."

    text = f"""# Step 1 SBE terminal processing

Generated on: {datetime.now().isoformat(timespec="seconds")}

## Purpose

This folder contains the internal audit trail for official Sea-Bird SBE Data Processing runs controlled from Python.

Python is used to copy files, build commands, run official Sea-Bird modules, record logs and save summaries.

## Processing design

The notebook runs the linear CTD sequence first:

```text
DatCnv > AlignCTD > WildEdit > CellTM > Filter > LoopEdit > Derive
```

After Derive, the workflow branches into four Bin Average products from the same unbinned derived CNV:

```text
1 m full cast
1 m downcast only
1 m upcast only
1 s full cast for SUNA matching
```

## Why a no-space working folder is used

The Sea-Bird command line modules can fail when input, output, configuration or PSA paths contain spaces.

This version copies working inputs to:

```text
{SBE_WORK_ROOT}
```

## Active raw input folder

```text
{RAW_INPUT_ROOT}
```

## PSA setup folder

```text
{PSA_ROOT}
```

## SBE executable folder

```text
{SBE_BIN_DIR}
```

## Internal output audit folder

```text
{OUTPUT_ROOT}
```

## L1 CTD deliverable folder

```text
{L1_CTD_ROOT}
```

## Module sequence and branches

{module_table}

## Module validation

{validation_table}
"""
    (output_root / "README_sbe_terminal_processing.md").write_text(text, encoding="utf-8")


# ==========================================================================
# 4. Preflight checks
# ==========================================================================


OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

try:
    module_validation_df = validate_processing_setup()
    preflight = module_validation_df.copy()

    write_validation_report(
        output_root=OUTPUT_ROOT,
        module_df=module_validation_df,
    )

    write_markdown_readme(
        output_root=OUTPUT_ROOT,
        module_validation_df=module_validation_df,
    )

    print(preflight)

    print("Preflight check passed.")
    print("Original source folders, SeaBird executables and usable PSA files were found.")
    print(f"Preflight files written to: {OUTPUT_ROOT}")

except Exception as error:
    print("Preflight check failed.")
    print("\nProblem:")
    print(error)
    print("\nFix this before setting RUN_SBE_COMMANDS = True.")
    raise


# ==========================================================================
# 5. Find raw files and selected configuration file
# ==========================================================================

raw_files = find_files_by_patterns(
    root=RAW_INPUT_ROOT,
    patterns=RAW_FILE_PATTERNS,
)

if TEST_SINGLE_CAST_ONLY:
    raw_files = [p for p in raw_files if p.stem == TEST_CAST_ID]
    if not raw_files:
        raise FileNotFoundError(
            f"TEST_SINGLE_CAST_ONLY is True, but no raw file named {TEST_CAST_ID}.hex was found under:\n"
            f"{RAW_INPUT_ROOT}\n\n"
            "Either check TEST_CAST_ID or set TEST_SINGLE_CAST_ONLY = False."
        )

if not raw_files:
    raise FileNotFoundError(
        f"No raw files found under:\n{RAW_INPUT_ROOT}\n\n"
        f"Patterns used:\n{RAW_FILE_PATTERNS}"
    )

sidecar_xml_files = sorted(set(RAW_INPUT_ROOT.rglob("*.XML")) | set(RAW_INPUT_ROOT.rglob("*.xml")))

if USER_SUPPLIED_XMLCON is not None:
    selected_xmlcon = validate_or_raise_config(
        Path(USER_SUPPLIED_XMLCON).expanduser().resolve(),
        context="USER_SUPPLIED_XMLCON from settings cell",
    )
    config_files = [selected_xmlcon]
else:
    config_files = find_files_by_patterns(
        root=RAW_INPUT_ROOT,
        patterns=CONFIG_FILE_PATTERNS,
    )
    if not config_files and not BORROW_MISSING_XMLCON:
        raise FileNotFoundError(
            "No valid .xmlcon or .con configuration files were found and BORROW_MISSING_XMLCON is False."
        )

raw_inventory_preview = pd.DataFrame(
    {
        "raw_file": [str(p) for p in raw_files],
        "cast_id": [p.stem for p in raw_files],
        "size_bytes": [p.stat().st_size for p in raw_files],
    }
)

print(raw_inventory_preview)
print(f"Raw files selected for this run: {len(raw_files)}")
print(f"TEST_SINGLE_CAST_ONLY: {TEST_SINGLE_CAST_ONLY}")
print(f"TEST_CAST_ID: {TEST_CAST_ID}")


# ==========================================================================
# 5A. Configuration inventory
# ==========================================================================


if USER_SUPPLIED_XMLCON is not None:
    print("USER_SUPPLIED_XMLCON is active. Skipping XMLCON borrowing and raw-tree config search.")
    config_files = [
        validate_or_raise_config(
            Path(USER_SUPPLIED_XMLCON).expanduser().resolve(),
            context="USER_SUPPLIED_XMLCON inventory refresh",
        )
    ]
else:
    if BORROW_MISSING_XMLCON:
        print("BORROW_MISSING_XMLCON = True, checking whether safe borrowing is possible...")
        xmlcon_borrow_log = borrow_missing_xmlcon_files(RAW_INPUT_ROOT)
        print(xmlcon_borrow_log)
    else:
        print("BORROW_MISSING_XMLCON = False, no configuration files were copied.")

    config_files = find_files_by_patterns(
        root=RAW_INPUT_ROOT,
        patterns=CONFIG_FILE_PATTERNS,
    )

if not config_files:
    raise FileNotFoundError(
        "No .con or .xmlcon files are available. Set USER_SUPPLIED_XMLCON to the correct SeaBird configuration file."
    )

config_inventory_records = []
for config_file in config_files:
    config_file = validate_or_raise_config(Path(config_file).resolve(), context="config inventory")
    stat = config_file.stat()
    config_inventory_records.append({
        "config_stem": config_file.stem,
        "config_suffix": config_file.suffix.lower(),
        "config_file": str(config_file),
        "config_size_kb": round(stat.st_size / 1024, 2),
        "config_date_modified": datetime.fromtimestamp(stat.st_mtime).isoformat(timespec="seconds"),
        "config_source": "USER_SUPPLIED_XMLCON" if USER_SUPPLIED_XMLCON is not None else "auto_search",
    })

config_inventory = pd.DataFrame(config_inventory_records)
config_inventory.to_csv(OUTPUT_ROOT / "01_config_file_inventory.csv", index=False, encoding="utf-8-sig")

print(config_inventory)
print(f"Configuration files available after inventory step: {len(config_inventory)}")


# ==========================================================================
# 5B. Match raw files to the selected configuration file
# ==========================================================================


inventory_records = []

for raw_file in raw_files:
    config_file = choose_config_for_raw(raw_file=raw_file, config_files=config_files)
    raw_stat = raw_file.stat()

    if config_file is not None:
        config_match_status = "matched"
        config_suffix = config_file.suffix.lower()
        config_stem = config_file.stem
    else:
        config_match_status = "not_matched_manual_review_needed"
        config_suffix = ""
        config_stem = ""

    inventory_records.append({
        "cast_id": raw_file.stem,
        "raw_suffix": raw_file.suffix.lower(),
        "raw_file": str(raw_file),
        "config_match_status": config_match_status,
        "config_stem": config_stem,
        "config_suffix": config_suffix,
        "config_file": str(config_file) if config_file else "",
        "raw_size_kb": round(raw_stat.st_size / 1024, 2),
        "raw_date_modified": datetime.fromtimestamp(raw_stat.st_mtime).isoformat(timespec="seconds"),
    })

inventory = pd.DataFrame(inventory_records)
inventory.to_csv(OUTPUT_ROOT / "01_raw_file_inventory.csv", index=False, encoding="utf-8-sig")

print(inventory)
print(inventory["config_match_status"].value_counts(dropna=False))


# ==========================================================================
# 5C. Stop if raw/config matching is unsafe
# ==========================================================================


unmatched = inventory[inventory["config_match_status"] != "matched"].copy()
duplicate_cast_ids = inventory[inventory.duplicated(subset=["cast_id"], keep=False)].copy()

if not unmatched.empty:
    unmatched.to_csv(OUTPUT_ROOT / "01_unmatched_raw_files.csv", index=False, encoding="utf-8-sig")
    print(unmatched)
    raise FileNotFoundError(
        "Some raw files could not be matched to a configuration file.\n\n"
        "Do not run official SeaBird processing yet.\n"
        "Fix the raw to configuration file matching first.\n\n"
        f"Review this file:\n{OUTPUT_ROOT / '01_unmatched_raw_files.csv'}"
    )

if not duplicate_cast_ids.empty:
    duplicate_cast_ids.to_csv(OUTPUT_ROOT / "01_duplicate_cast_ids.csv", index=False, encoding="utf-8-sig")
    print(duplicate_cast_ids)
    raise ValueError(
        "Duplicate cast IDs were found.\n\n"
        "Rename files or place them in separate processing batches so outputs do not overwrite each other.\n\n"
        f"Review this file:\n{OUTPUT_ROOT / '01_duplicate_cast_ids.csv'}"
    )

print("Raw/config matching passed.")
print(f"Matched raw files: {len(inventory)}")


# ==========================================================================
# 5D. Prepare the no-space SeaBird working folder
# ==========================================================================


work_inventory = prepare_no_space_work_folder(
    inventory=inventory,
    clean_existing_work_folder=CLEAN_WORK_FOLDER_FIRST,
)

print(work_inventory)

print("No-space SeaBird working folder is ready.")
print(f"Working root: {SBE_WORK_ROOT}")
print(f"Working PSA folder: {SBE_WORK_PSA_ROOT}")
print(f"Working raw folder: {SBE_WORK_RAW_ROOT}")
print(f"Working output folder: {SBE_WORK_OUTPUT_ROOT}")
print(f"CLEAN_WORK_FOLDER_FIRST: {CLEAN_WORK_FOLDER_FIRST}")


# ==========================================================================
# 5E. Validate selected and working XMLCON files before running SeaBird
# ==========================================================================


xmlcon_check_records = []

for _, row in work_inventory.iterrows():
    cast_id = row["cast_id"]
    paths_to_check = [
        ("source_config_file", Path(row["source_config_file"])),
        ("work_config_file", Path(row["work_config_file"])),
        ("command_config_file", Path(row.get("command_config_file", row["work_config_file"]))),
    ]

    seen = set()
    for role, candidate in paths_to_check:
        key = (role, str(candidate).lower())
        if key in seen:
            continue
        seen.add(key)

        exists = candidate.exists()
        size_bytes = candidate.stat().st_size if exists else None
        ok, reason = is_probably_valid_seabird_config(candidate) if exists else (False, "missing")

        try:
            preview = " ".join(candidate.read_text(encoding="utf-8", errors="replace").split()[:30]) if exists else ""
        except Exception:
            preview = ""

        xmlcon_check_records.append({
            "cast_id": cast_id,
            "role": role,
            "candidate": str(candidate),
            "exists": exists,
            "size_bytes": size_bytes,
            "looks_like_valid_seabird_config": ok,
            "reason": reason,
            "preview_first_words": preview[:250],
        })

xmlcon_alias_check = pd.DataFrame(xmlcon_check_records)
xmlcon_alias_check.to_csv(OUTPUT_ROOT / "01_xmlcon_selected_config_check.csv", index=False, encoding="utf-8-sig")

print(xmlcon_alias_check)

bad_xmlcon_files = xmlcon_alias_check[
    (xmlcon_alias_check["exists"] == False)
    | (xmlcon_alias_check["looks_like_valid_seabird_config"] == False)
].copy()

if not bad_xmlcon_files.empty:
    print("A selected configuration file is missing or does not look like a valid SeaBird instrument configuration file.")
    print(f"Review: {OUTPUT_ROOT / '01_xmlcon_selected_config_check.csv'}")
    print(bad_xmlcon_files)
    raise ValueError(
        "Invalid or missing XMLCON/CON content detected before running DatCnv.\n\n"
        "Use the exact .xmlcon or .con file exported from SeaBird for this CTD package."
    )

print("Selected XMLCON/CON check completed.")
print(f"Review file: {OUTPUT_ROOT / '01_xmlcon_selected_config_check.csv'}")


# ==========================================================================
# 5F. Diagnose PSA contents for hardcoded paths and append flags
# ==========================================================================

import xml.etree.ElementTree as ET

def inspect_psa(psa_path: Path) -> Dict[str, Any]:
    """Pull the fields most likely to make a PSA hostile to CLI overrides."""
    info: Dict[str, Any] = {
        "psa": psa_path.name,
        "OutputDir": "",
        "OutputFile": "",
        "AppendOutputFile": "",
        "InputDir": "",
        "warnings": "",
    }
    try:
        tree = ET.parse(psa_path)
        root = tree.getroot()
        for elem in root.iter():
            tag = elem.tag
            val = elem.attrib.get("value", elem.text or "")
            if tag in info and info[tag] == "":
                info[tag] = str(val).strip()
    except Exception as exc:
        info["warnings"] = f"parse_error: {exc!r}"
        return info

    warnings = []
    if info["OutputDir"] and info["OutputDir"] not in {".", ""}:
        warnings.append("OutputDir is hardcoded - /o may be ignored")
    if info["OutputFile"]:
        warnings.append("OutputFile is set - /f may be ignored or appended")
    if info["AppendOutputFile"] not in {"", "0"}:
        warnings.append("AppendOutputFile=1 - SBE will append a suffix to the INPUT stem")
    info["warnings"] = "; ".join(warnings)
    return info

psa_diag_records = []
for psa_path in sorted(SBE_WORK_PSA_ROOT.glob("*.psa")):
    psa_diag_records.append(inspect_psa(psa_path))

psa_diag = pd.DataFrame(psa_diag_records)
psa_diag.to_csv(OUTPUT_ROOT / "01_psa_diagnostic.csv", index=False, encoding="utf-8-sig")
print(psa_diag)

problem_psas = psa_diag[psa_diag["warnings"] != ""]
if not problem_psas.empty:
    print()
    print("One or more PSA files have hardcoded settings that may override the CLI:")
    print("Open each one in SBE Data Processing GUI, set 'Output file' to use the input")
    print("directory and base name (no append, no hardcoded path), and save.")
    print()
    print("This is the most common cause of return_code=0 + missing-output failures")
    print("(see SBE Data Processing User's Manual, 'Running modules from the command line').")
else:
    print("No problematic hardcoded settings detected in PSAs.")


# ==========================================================================
# 6. Run SBE processing sequence
# ==========================================================================

processing_records: List[Dict[str, Any]] = []
deliverable_records: List[Dict[str, Any]] = []

OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
SBE_WORK_OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
L1_CTD_ROOT.mkdir(parents=True, exist_ok=True)
L2_CTD_ROOT.mkdir(parents=True, exist_ok=True)
L2_DERIVED.mkdir(parents=True, exist_ok=True)

# Guard FIRST. If Section 5D was not run, the no-space PSA folder will not
# contain the copied PSA files and the command builder cannot be trusted.
if "work_inventory" not in globals() or work_inventory.empty:
    raise RuntimeError(
        "work_inventory does not exist or is empty. "
        "Run Section 5D before running Section 6."
    )

linear_modules = enabled_linear_modules()
binavg_modules = enabled_binavg_branches()
ascii_modules = enabled_ascii_branches()

all_modules_for_lookup = linear_modules + binavg_modules + ascii_modules

exe_lookup: Dict[str, Path] = {}
psa_lookup: Dict[str, Path] = {}

for module in all_modules_for_lookup:
    exe_lookup[module["name"]] = find_executable(module["exe"])
    psa_lookup[module["name"]] = find_work_psa(module)

unmatched_inventory = inventory[inventory["config_match_status"] != "matched"].copy()
if not unmatched_inventory.empty:
    print(unmatched_inventory)
    raise FileNotFoundError("Section 6 stopped because some raw files do not have matched configuration files.")

stop_entire_batch = False

for _, file_row in inventory.iterrows():
    if stop_entire_batch:
        break

    cast_id = safe_name(file_row["cast_id"])
    match = work_inventory.loc[work_inventory["cast_id"] == cast_id]

    if match.empty:
        processing_records.append({
            "cast_id": cast_id,
            "module": "all",
            "description": "All modules skipped because no working copy was available",
            "phase": "setup",
            "input_file": "",
            "output_dir": "",
            "expected_output_file": "",
            "actual_output_file": "",
            "psa_file": "",
            "config_file": "",
            "command": "",
            "run_requested": RUN_SBE_COMMANDS,
            "status": "skipped_no_working_copy",
            "return_code": "",
            "log_path": "",
        })
        if STOP_ON_ERROR:
            stop_entire_batch = True
        continue

    work_row = match.iloc[0]
    raw_file = Path(work_row["work_raw_file"])
    work_config_file = Path(work_row["work_config_file"])
    config_file = Path(work_row.get("command_config_file", work_row["work_config_file"]))

    if not raw_file.exists() or not config_file.exists():
        processing_records.append({
            "cast_id": cast_id,
            "module": "all",
            "description": "All modules skipped because working raw or configuration file was missing",
            "phase": "setup",
            "input_file": str(raw_file),
            "output_dir": "",
            "expected_output_file": "",
            "actual_output_file": "",
            "psa_file": "",
            "config_file": str(config_file),
            "command": "",
            "run_requested": RUN_SBE_COMMANDS,
            "status": "skipped_missing_working_input",
            "return_code": "",
            "log_path": "",
        })
        if STOP_ON_ERROR:
            stop_entire_batch = True
        continue

    cast_output_dir = SBE_WORK_OUTPUT_ROOT / "per_cast_sbe_outputs" / cast_id
    stage_dirs = {
        "raw_copy": cast_output_dir / "00_raw_copy",
        "logs": cast_output_dir / "logs",
        "commands": cast_output_dir / "commands",
    }
    for folder in stage_dirs.values():
        folder.mkdir(parents=True, exist_ok=True)

    if COPY_INPUT_FILES_TO_CAST_FOLDER:
        shutil.copy2(raw_file, stage_dirs["raw_copy"] / raw_file.name)
        shutil.copy2(config_file, stage_dirs["raw_copy"] / config_file.name)
        if work_config_file.exists() and work_config_file.resolve() != config_file.resolve():
            shutil.copy2(work_config_file, stage_dirs["raw_copy"] / work_config_file.name)

    cast_commands: List[str] = []
    cast_failed = False

    # ------------------------------------------------------------------
    # A. Linear sequence: DatCnv > AlignCTD > [WildEdit] > [CellTM] > [Filter] > LoopEdit > Derive
    #    Bracketed modules are gated by APPLY_* flags and skipped when disabled.
    # ------------------------------------------------------------------
    current_input = raw_file
    derived_unbinned_output: Optional[Path] = None

    for module in linear_modules:
        module_name = module["name"]
        exe_path = exe_lookup[module_name]
        psa_file = psa_lookup[module_name]
        module_output_dir = Path(module["output_folder"])
        module_output_dir.mkdir(parents=True, exist_ok=True)

        output_file = module_output_file_for_cast(module, cast_id)
        expected_output_path = module_output_dir / output_file

        if RUN_SBE_COMMANDS and expected_output_path.exists():
            expected_output_path.unlink()

        command = build_sbe_command(
            exe_path=exe_path,
            input_file=current_input,
            output_dir=module_output_dir,
            output_file=output_file,
            psa_file=psa_file,
            config_file=config_file if module.get("needs_config", False) else None,
            extra_args=module.get("extra_args", []),
        )

        command_text = command_to_text(command)
        cast_commands.append(command_text)
        log_path = stage_dirs["logs"] / f"{module_name}.log"

        record = {
            "cast_id": cast_id,
            "module": module_name,
            "description": module["description"],
            "phase": "linear",
            "input_file": str(current_input),
            "output_dir": str(module_output_dir),
            "expected_output_file": str(expected_output_path),
            "actual_output_file": "",
            "psa_file": str(psa_file),
            "config_file": str(config_file) if module.get("needs_config", False) else "",
            "command": command_text,
            "command_working_dir": str(current_input.parent if RUN_FROM_INPUT_FILE_FOLDER else exe_path.parent),
            "run_requested": RUN_SBE_COMMANDS,
            "status": "planned_not_run",
            "return_code": "",
            "log_path": str(log_path),
        }

        if RUN_SBE_COMMANDS:
            print(f"Running {module_name} for {cast_id}")
            command_working_dir = current_input.parent if RUN_FROM_INPUT_FILE_FOLDER else exe_path.parent
            module_started_at = datetime.now()

            return_code, stdout, stderr = run_command(
                command=command,
                log_path=log_path,
                working_dir=command_working_dir,
            )
            record["return_code"] = return_code

            expected_output = find_expected_output(
                output_dir=module_output_dir,
                output_file=output_file,
                started_at=module_started_at,
            )

            if (
                return_code == 0
                and expected_output is not None
                and expected_output.name != output_file
            ):
                canonical_path = module_output_dir / output_file
                try:
                    expected_output.rename(canonical_path)
                    record["output_renamed_from"] = expected_output.name
                    expected_output = canonical_path
                except OSError as rename_err:
                    record["output_renamed_from"] = (
                        f"FAILED_TO_RENAME:{expected_output.name}:{rename_err}"
                    )

            if return_code == 0 and expected_output is not None:
                record["status"] = "success"
                record["actual_output_file"] = str(expected_output)
                current_input = expected_output
                if module_name == FINAL_LINEAR_MODULE_NAME:
                    derived_unbinned_output = expected_output
            else:
                record["status"] = "failed_or_output_missing"
                cast_failed = True
                processing_records.append(record)

                print(f"Problem while running {module_name} for {cast_id}.")
                print(f"  Return code: {return_code}")
                print(f"  Expected output: {expected_output_path}")
                print(f"  Log file: {log_path}")
                if stderr.strip():
                    print("  STDERR excerpt:")
                    for line in stderr.strip().split("\n")[:8]:
                        print(f"    {line[:160]}")
                if stdout.strip() and not stderr.strip():
                    print("  STDOUT excerpt:")
                    for line in stdout.strip().split("\n")[:8]:
                        print(f"    {line[:160]}")

                if STOP_ON_ERROR:
                    stop_entire_batch = True
                    print("Stopping because STOP_ON_ERROR is True.")
                break
        else:
            current_input = expected_output_path
            record["actual_output_file"] = str(expected_output_path)
            if module_name == FINAL_LINEAR_MODULE_NAME:
                derived_unbinned_output = expected_output_path

        processing_records.append(record)

    if cast_failed:
        command_file = stage_dirs["commands"] / f"{cast_id}_sbe_commands.txt"
        command_file.write_text("\n\n".join(cast_commands), encoding="utf-8")
        continue

    if derived_unbinned_output is None:
        raise RuntimeError(
            f"Could not determine the derived unbinned CNV for {cast_id}. "
            f"Expected final linear module: {FINAL_LINEAR_MODULE_NAME}"
        )

    # In v2 the derived unbinned CNV already lives in its canonical L2 home
    # (L2/CTD/Derived_Parameter). The L1 deliverable is DatCnv's readable
    # full-cast CNV, which 01_datcnv wrote directly to L1/CTD. So there is
    # nothing to mirror here. derived_unbinned_output feeds the bin branches.

    # ------------------------------------------------------------------
    # B. Bin Average branches. Each branch starts from the same derived unbinned CNV.
    # ------------------------------------------------------------------
    binavg_outputs_by_suffix: Dict[str, Path] = {}

    for module in binavg_modules:
        module_name = module["name"]
        exe_path = exe_lookup[module_name]
        psa_file = psa_lookup[module_name]
        module_output_dir = Path(module["output_folder"])
        module_output_dir.mkdir(parents=True, exist_ok=True)

        output_file = module_output_file_for_cast(module, cast_id)
        expected_output_path = module_output_dir / output_file

        if RUN_SBE_COMMANDS and expected_output_path.exists():
            expected_output_path.unlink()

        command = build_sbe_command(
            exe_path=exe_path,
            input_file=derived_unbinned_output,
            output_dir=module_output_dir,
            output_file=output_file,
            psa_file=psa_file,
            config_file=None,
            extra_args=module.get("extra_args", []),
        )

        command_text = command_to_text(command)
        cast_commands.append(command_text)
        log_path = stage_dirs["logs"] / f"{module_name}.log"

        record = {
            "cast_id": cast_id,
            "module": module_name,
            "description": module["description"],
            "phase": "binavg_branch",
            "input_file": str(derived_unbinned_output),
            "output_dir": str(module_output_dir),
            "expected_output_file": str(expected_output_path),
            "actual_output_file": "",
            "psa_file": str(psa_file),
            "config_file": "",
            "command": command_text,
            "command_working_dir": str(derived_unbinned_output.parent if RUN_FROM_INPUT_FILE_FOLDER else exe_path.parent),
            "run_requested": RUN_SBE_COMMANDS,
            "status": "planned_not_run",
            "return_code": "",
            "log_path": str(log_path),
        }

        if RUN_SBE_COMMANDS:
            print(f"Running {module_name} for {cast_id}")
            command_working_dir = derived_unbinned_output.parent if RUN_FROM_INPUT_FILE_FOLDER else exe_path.parent
            module_started_at = datetime.now()

            return_code, stdout, stderr = run_command(
                command=command,
                log_path=log_path,
                working_dir=command_working_dir,
            )
            record["return_code"] = return_code

            expected_output = find_expected_output(
                output_dir=module_output_dir,
                output_file=output_file,
                started_at=module_started_at,
            )

            if (
                return_code == 0
                and expected_output is not None
                and expected_output.name != output_file
            ):
                canonical_path = module_output_dir / output_file
                try:
                    expected_output.rename(canonical_path)
                    record["output_renamed_from"] = expected_output.name
                    expected_output = canonical_path
                except OSError as rename_err:
                    record["output_renamed_from"] = (
                        f"FAILED_TO_RENAME:{expected_output.name}:{rename_err}"
                    )

            if return_code == 0 and expected_output is not None:
                record["status"] = "success"
                record["actual_output_file"] = str(expected_output)
                binavg_outputs_by_suffix[module["output_suffix"]] = expected_output

                audit_binavg_path = OUTPUT_ROOT / "08_binavg_cnv" / module_name / output_file
                copy_with_log(
                    source=expected_output,
                    destination=audit_binavg_path,
                    records=deliverable_records,
                    cast_id=cast_id,
                    product_type=f"audit_copy_{module_name}",
                )
            else:
                record["status"] = "failed_or_output_missing"
                cast_failed = True
                processing_records.append(record)

                print(f"Problem while running {module_name} for {cast_id}.")
                print(f"  Return code: {return_code}")
                print(f"  Expected output: {expected_output_path}")
                print(f"  Log file: {log_path}")
                if stderr.strip():
                    print("  STDERR excerpt:")
                    for line in stderr.strip().split("\n")[:8]:
                        print(f"    {line[:160]}")
                if stdout.strip() and not stderr.strip():
                    print("  STDOUT excerpt:")
                    for line in stdout.strip().split("\n")[:8]:
                        print(f"    {line[:160]}")

                if STOP_ON_ERROR:
                    stop_entire_batch = True
                    print("Stopping because STOP_ON_ERROR is True.")
                break
        else:
            record["actual_output_file"] = str(expected_output_path)
            binavg_outputs_by_suffix[module["output_suffix"]] = expected_output_path

        processing_records.append(record)

    if cast_failed:
        command_file = stage_dirs["commands"] / f"{cast_id}_sbe_commands.txt"
        command_file.write_text("\n\n".join(cast_commands), encoding="utf-8")
        continue

    # ------------------------------------------------------------------
    # C. ASCII Out branches. Each ASCII branch reads the matching binned CNV.
    # ------------------------------------------------------------------
    for module in ascii_modules:
        module_name = module["name"]
        exe_path = exe_lookup[module_name]
        psa_file = psa_lookup[module_name]
        module_output_dir = Path(module["output_folder"])
        module_output_dir.mkdir(parents=True, exist_ok=True)

        input_file = Path(module["input_folder"]) / f"{cast_id}{module['input_suffix']}"
        output_file = module_output_file_for_cast(module, cast_id)
        expected_output_path = module_output_dir / output_file

        if RUN_SBE_COMMANDS and expected_output_path.exists():
            expected_output_path.unlink()

        command = build_sbe_command(
            exe_path=exe_path,
            input_file=input_file,
            output_dir=module_output_dir,
            output_file=output_file,
            psa_file=psa_file,
            config_file=None,
            extra_args=module.get("extra_args", []),
        )

        command_text = command_to_text(command)
        cast_commands.append(command_text)
        log_path = stage_dirs["logs"] / f"{module_name}.log"

        record = {
            "cast_id": cast_id,
            "module": module_name,
            "description": module["description"],
            "phase": "ascii_branch",
            "input_file": str(input_file),
            "output_dir": str(module_output_dir),
            "expected_output_file": str(expected_output_path),
            "actual_output_file": "",
            "psa_file": str(psa_file),
            "config_file": "",
            "command": command_text,
            "command_working_dir": str(input_file.parent if RUN_FROM_INPUT_FILE_FOLDER else exe_path.parent),
            "run_requested": RUN_SBE_COMMANDS,
            "status": "planned_not_run",
            "return_code": "",
            "log_path": str(log_path),
        }

        if RUN_SBE_COMMANDS:
            if not input_file.exists():
                record["status"] = "skipped_missing_ascii_input"
                record["return_code"] = ""
                cast_failed = True
                processing_records.append(record)
                print(f"Skipping {module_name} because input does not exist: {input_file}")
                if STOP_ON_ERROR:
                    stop_entire_batch = True
                break

            print(f"Running {module_name} for {cast_id}")
            command_working_dir = input_file.parent if RUN_FROM_INPUT_FILE_FOLDER else exe_path.parent
            module_started_at = datetime.now()

            return_code, stdout, stderr = run_command(
                command=command,
                log_path=log_path,
                working_dir=command_working_dir,
            )
            record["return_code"] = return_code

            expected_output = find_expected_output(
                output_dir=module_output_dir,
                output_file=output_file,
                started_at=module_started_at,
            )

            if (
                return_code == 0
                and expected_output is not None
                and expected_output.name != output_file
            ):
                canonical_path = module_output_dir / output_file
                try:
                    expected_output.rename(canonical_path)
                    record["output_renamed_from"] = expected_output.name
                    expected_output = canonical_path
                except OSError as rename_err:
                    record["output_renamed_from"] = (
                        f"FAILED_TO_RENAME:{expected_output.name}:{rename_err}"
                    )

            if return_code == 0 and expected_output is not None:
                record["status"] = "success"
                record["actual_output_file"] = str(expected_output)

                audit_ascii_path = OUTPUT_ROOT / "09_ascii_csv" / module_name / output_file
                copy_with_log(
                    source=expected_output,
                    destination=audit_ascii_path,
                    records=deliverable_records,
                    cast_id=cast_id,
                    product_type=f"audit_copy_{module_name}",
                )
            else:
                record["status"] = "failed_or_output_missing"
                cast_failed = True
                processing_records.append(record)

                print(f"Problem while running {module_name} for {cast_id}.")
                print(f"  Return code: {return_code}")
                print(f"  Expected output: {expected_output_path}")
                print(f"  Log file: {log_path}")
                if stderr.strip():
                    print("  STDERR excerpt:")
                    for line in stderr.strip().split("\n")[:8]:
                        print(f"    {line[:160]}")
                if stdout.strip() and not stderr.strip():
                    print("  STDOUT excerpt:")
                    for line in stdout.strip().split("\n")[:8]:
                        print(f"    {line[:160]}")

                if STOP_ON_ERROR:
                    stop_entire_batch = True
                    print("Stopping because STOP_ON_ERROR is True.")
                break
        else:
            record["actual_output_file"] = str(expected_output_path)

        processing_records.append(record)

    command_file = stage_dirs["commands"] / f"{cast_id}_sbe_commands.txt"
    command_file.write_text("\n\n".join(cast_commands), encoding="utf-8")

    cast_summary = {
        "cast_id": cast_id,
        "source_raw_file": str(file_row["raw_file"]),
        "source_config_file": str(file_row["config_file"]),
        "work_raw_file": str(raw_file),
        "work_config_file": str(work_config_file),
        "command_config_file": str(config_file),
        "cast_output_dir": str(cast_output_dir),
        "commands_file": str(command_file),
        "run_requested": RUN_SBE_COMMANDS,
        "cast_status": "failed" if cast_failed else "planned" if not RUN_SBE_COMMANDS else "completed",
    }
    pd.DataFrame([cast_summary]).to_csv(
        cast_output_dir / f"{cast_id}_cast_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

processing_log = pd.DataFrame(processing_records)
processing_log.to_csv(OUTPUT_ROOT / "02_sbe_processing_command_log.csv", index=False, encoding="utf-8-sig")

deliverables_log = pd.DataFrame(deliverable_records)
deliverables_log.to_csv(OUTPUT_ROOT / "02_l1_deliverables_copy_log.csv", index=False, encoding="utf-8-sig")

print(processing_log)

if not deliverables_log.empty:
    print(deliverables_log)

if RUN_SBE_COMMANDS:
    failed_steps = processing_log[
        processing_log["status"].isin(["failed_or_output_missing", "skipped_missing_ascii_input"])
    ].copy()
    if failed_steps.empty:
        print("SBE processing completed with no failed module records.")
    else:
        failed_steps.to_csv(OUTPUT_ROOT / "02_failed_sbe_processing_steps.csv", index=False, encoding="utf-8-sig")
        print(failed_steps)
        print("Some SBE processing steps failed or did not produce the expected output.")
        print(f"Review: {OUTPUT_ROOT / '02_failed_sbe_processing_steps.csv'}")
else:
    print("Dry run complete. No Sea-Bird commands were executed.")
    print("Review the command log before setting RUN_SBE_COMMANDS = True.")

print(f"Processing log written to: {OUTPUT_ROOT / '02_sbe_processing_command_log.csv'}")
print(f"Deliverables copy log written to: {OUTPUT_ROOT / '02_l1_deliverables_copy_log.csv'}")


# ==========================================================================
# 7. Review L1 deliverables
# ==========================================================================

deliverable_index_records: List[Dict[str, Any]] = []

# v2 product map: L1 readable CNV + the two L2 binned products (+ optional ASCII).
expected_products = [
    ("L1_cnv_full",    L1_CTD_ROOT,            ".cnv"),
    ("L2_cnv_1m_down", L2_1M_DOWN,             "_1m_down.cnv"),
    ("L2_cnv_1s_full", L2_1S,                  "_1s.cnv"),
]
if RUN_ASCII_OUT:
    expected_products += [
        ("L2_ascii_1m_down", L2_1M_DOWN / "asc_semicolon", "_1m_down_ascii.asc"),
        ("L2_ascii_1s_full", L2_1S / "asc_semicolon",      "_1s_ascii.asc"),
    ]

for cast_id in inventory["cast_id"].map(safe_name).tolist():
    for product_name, product_folder, product_suffix in expected_products:
        expected_path = product_folder / f"{cast_id}{product_suffix}"
        deliverable_index_records.append({
            "cast_id": cast_id,
            "product_name": product_name,
            "expected_path": str(expected_path),
            "exists": expected_path.exists(),
            "size_bytes": expected_path.stat().st_size if expected_path.exists() else None,
            "run_requested": RUN_SBE_COMMANDS,
        })

deliverable_index = pd.DataFrame(deliverable_index_records)
deliverable_index.to_csv(
    OUTPUT_ROOT / "03_l1_l2_deliverable_index.csv",
    index=False,
    encoding="utf-8-sig",
)

print(deliverable_index)

if RUN_SBE_COMMANDS:
    missing = deliverable_index[deliverable_index["exists"] == False]
    if missing.empty:
        print("All expected L1/L2 deliverables are present.")
    else:
        missing.to_csv(OUTPUT_ROOT / "03_missing_deliverables.csv", index=False, encoding="utf-8-sig")
        print(missing)
        print(f"Some deliverables are missing. Review: {OUTPUT_ROOT / '03_missing_deliverables.csv'}")


# ==========================================================================
# 8. Save summary workbook and README
# ==========================================================================

summary_xlsx = OUTPUT_ROOT / "sbe_terminal_processing_summary.xlsx"
readme_path = OUTPUT_ROOT / "README_sbe_terminal_processing.md"
OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)


def variable_exists(name: str) -> bool:
    return name in globals() and globals()[name] is not None


def get_dataframe_if_exists(name: str) -> Optional[pd.DataFrame]:
    if not variable_exists(name):
        return None
    value = globals()[name]
    return value if isinstance(value, pd.DataFrame) else None


def add_sheet_if_available(writer: pd.ExcelWriter, df_name: str, sheet_name: str) -> bool:
    df = get_dataframe_if_exists(df_name)
    if df is None:
        return False
    df.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    return True


write_markdown_readme(
    output_root=OUTPUT_ROOT,
    module_validation_df=module_validation_df if variable_exists("module_validation_df") else None,
)

written_sheets: List[str] = []

try:
    with pd.ExcelWriter(summary_xlsx, engine="openpyxl") as writer:
        sheet_plan = [
            ("module_validation_df", "module_preflight"),
            ("preflight", "preflight"),
            ("config_inventory", "config_inventory"),
            ("inventory", "raw_inventory"),
            ("work_inventory", "working_copy_inventory"),
            ("processing_log", "processing_log"),
            ("deliverables_log", "deliverables_copy_log"),
            ("deliverable_index", "l1_deliverable_index"),
            ("failed_steps", "failed_steps"),
        ]

        for df_name, sheet_name in sheet_plan:
            was_written = add_sheet_if_available(writer=writer, df_name=df_name, sheet_name=sheet_name)
            if was_written:
                written_sheets.append(sheet_name[:31])

        run_summary_records = [
            {"item": "generated_on", "value": datetime.now().isoformat(timespec="seconds")},
            {"item": "ctd_root", "value": str(CTD_ROOT)},
            {"item": "raw_input_root", "value": str(RAW_INPUT_ROOT)},
            {"item": "cruise_id", "value": CRUISE_ID},
            {"item": "l1_ctd_root", "value": str(L1_CTD_ROOT)},
            {"item": "l2_derived", "value": str(L2_DERIVED)},
            {"item": "l2_1m_down", "value": str(L2_1M_DOWN)},
            {"item": "l2_1s", "value": str(L2_1S)},
            {"item": "psa_root", "value": str(PSA_ROOT)},
            {"item": "sbe_bin_dir", "value": str(SBE_BIN_DIR)},
            {"item": "sbe_work_root", "value": str(SBE_WORK_ROOT)},
            {"item": "sbe_work_output_root", "value": str(SBE_WORK_OUTPUT_ROOT)},
            {"item": "output_root", "value": str(OUTPUT_ROOT)},
            {"item": "test_single_cast_only", "value": str(TEST_SINGLE_CAST_ONLY)},
            {"item": "test_cast_id", "value": str(TEST_CAST_ID)},
            {"item": "run_ascii_out", "value": str(RUN_ASCII_OUT)},
            {"item": "run_sbe_commands", "value": str(RUN_SBE_COMMANDS)},
            {"item": "user_supplied_xmlcon", "value": str(USER_SUPPLIED_XMLCON)},
            {"item": "use_original_user_xmlcon_for_command", "value": str(USE_ORIGINAL_USER_XMLCON_FOR_COMMAND)},
            {"item": "create_cast_config_alias", "value": str(CREATE_CAST_CONFIG_ALIAS)},
            {"item": "use_local_input_and_config_names", "value": str(USE_LOCAL_INPUT_AND_CONFIG_NAMES)},
            {"item": "copy_input_files_to_cast_folder", "value": str(COPY_INPUT_FILES_TO_CAST_FOLDER)},
            {"item": "clean_work_folder_first", "value": str(CLEAN_WORK_FOLDER_FIRST)},
            {"item": "stop_on_error", "value": str(STOP_ON_ERROR)},
            {"item": "final_linear_module_name", "value": str(FINAL_LINEAR_MODULE_NAME)},
            {"item": "summary_workbook", "value": str(summary_xlsx)},
            {"item": "readme", "value": str(readme_path)},
            {"item": "sheets_written", "value": ", ".join(written_sheets)},
        ]
        pd.DataFrame(run_summary_records).to_excel(writer, sheet_name="run_summary", index=False)

    print("Summary workbook saved:")
    print(summary_xlsx)

    print("\nSheets written:")
    for sheet in written_sheets:
        print(f"  {sheet}")
    print("  run_summary")

except Exception as error:
    print("Could not write the Excel summary workbook.")
    print(error)
    print("CSV logs and README were still written to OUTPUT_ROOT.")

print("\nREADME saved:")
print(readme_path)
print("\nDone.")

print("\n[step01] Done. Review the tables above and the _audit folder.")
