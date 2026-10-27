"""
combine_readable_profiles.py

Stack the per-cast lab-friendly profile CSVs that step04 writes
(`<cast>_readable.csv`) into ONE combined table:

    combined_all_casts_final_1m_down_profiles.csv

Each row is one depth bin of one cast; Cruise_ID and Cast_ID columns identify
where each row came from. This is the cruise-wide table for a merge tool,
section plots, etc.

It uses the UNION of columns across casts, so casts that lack a variable
(e.g. no SUNA -> no Nitrate_uM) simply have blanks in that column instead of
breaking the concatenation.

USAGE (PyCharm Community: edit CRUISE_ID, then Run; or from the terminal):
    python combine_readable_profiles.py
Optional override of the cruise from the command line:
    python combine_readable_profiles.py P45_05
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
CTD_ROOT   = Path(r"C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd")
CRUISE_ID  = "P45_06"                       # overridden by a command-line arg if given
if len(sys.argv) > 1:
    CRUISE_ID = sys.argv[1]

CRUISE_DIR = CTD_ROOT / "cruises" / CRUISE_ID
PLOT_ROOT  = CRUISE_DIR / "L2" / "CTD" / "profile_plots"   # where step04 wrote _readable.csv
OUT_DIR    = CRUISE_DIR / "L2" / "CTD"                     # combined file lands in L2/CTD
OUT_NAME   = "combined_all_casts_final_1m_down_profiles.csv"

# Preferred left-to-right column order (any extras are appended after these).
PREFERRED_ORDER = [
    "Cruise_ID", "Cast_ID",
    "Depth_m", "Pressure_db",
    "Temperature_C", "PotentialTemp_C",
    "Salinity_PSU", "SigmaT_kg_m3",
    "Oxygen_umol_kg", "Oxygen_umol_L",
    "Fluorescence_mg_m3", "BeamAttenuation_1_m", "PAR",
    "UTC_time", "Nitrate_uM",
]


def cast_id_from_readable(path: Path) -> str:
    name = path.name
    return name[: -len("_readable.csv")] if name.endswith("_readable.csv") else path.stem


def order_columns(cols) -> list:
    front = [c for c in PREFERRED_ORDER if c in cols]
    rest = [c for c in cols if c not in front]
    return front + rest


def main() -> None:
    if not PLOT_ROOT.exists():
        raise SystemExit(f"profile_plots folder not found: {PLOT_ROOT}\n"
                         f"Run step04 first so the *_readable.csv files exist.")

    files = sorted(PLOT_ROOT.glob("*_readable.csv"))
    if not files:
        raise SystemExit(f"No *_readable.csv files in {PLOT_ROOT}. Run step04 first.")

    frames = []
    for f in files:
        cast_id = cast_id_from_readable(f)
        df = pd.read_csv(f)
        df.insert(0, "Cast_ID", cast_id)
        df.insert(0, "Cruise_ID", CRUISE_ID)
        frames.append(df)
        print(f"  {cast_id}: {len(df)} rows, {df.shape[1]} cols")

    combined = pd.concat(frames, ignore_index=True, sort=False)
    combined = combined[order_columns(list(combined.columns))]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / OUT_NAME
    combined.to_csv(out_path, index=False, encoding="utf-8-sig")

    print(f"\nCombined {len(files)} casts -> {len(combined)} rows, {combined.shape[1]} columns")
    print(f"Columns: {list(combined.columns)}")
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()
