"""
build_bottle_nitrate_table.py  -  Turn SBE Bottle Summary (.btl) files (+ optional lab
nitrate) into the discrete bottle products.

This closes the last hand-step in the bottle path: instead of reading firing depth/time
out of each .btl by eye, it parses them all into one tidy table, and (if you supply the
lab nitrate keyed by cast + bottle number) joins nitrate on exactly, producing the
step03-ready table.

Outputs (into OUTPUT_DIR):
  <CRUISE>_bottle_ctd.csv       always — the CTD side of every bottle firing:
                                cast_id, bottle_no, pressure_dbar, depth_m, time_utc,
                                temp_c, salinity
  <CRUISE>_bottle_nitrate.csv   only if a lab file is supplied and joined — the step03
                                schema (cruise_id, cast_id, station, depth_m, time_utc,
                                nitrate_uM, replicate_id, flag), plus an *_UNMATCHED.csv
                                listing bottles/lab rows that did not pair.

HOW TO RUN
  Edit the SETTINGS block, then:  python build_bottle_nitrate_table.py
  (BOTTLE_DIR can also be passed as the first command-line argument.)

Needs ctd_lib.py beside this file. pandas is required; openpyxl only for .xlsx lab files.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

import ctd_lib as L

# =============================== SETTINGS ==================================
CRUISE_ID  = "P45_07"
BOTTLE_DIR = Path(r"C:\Projects\ctd_pipeline\cruises\P45_07\L2\CTD\Bottle")
OUTPUT_DIR = Path(r"C:\Projects\ctd_pipeline\cruises\P45_07\metadata\bottle_nitrate")
LATITUDE   = 5.0        # for pressure -> depth (Gulf of Guinea ~5 N); refine if you have per-cast lat

# --- Optional lab nitrate join (leave LAB_FILE = None to only build the CTD table) ---
LAB_FILE   = None       # e.g. Path(r"C:\...\P45_07_lab_nitrate.csv")  (.csv or .xlsx)
LAB_SHEET  = 0          # sheet name/index if .xlsx
# Column names IN THE LAB FILE that identify each sample. The join is on (cast_id, bottle_no).
LAB_CAST_COL   = "cast_id"     # must hold values like "P45_07_CTD_02"
LAB_BOTTLE_COL = "bottle_no"   # integer bottle/Niskin position matching the .btl
LAB_NITRATE_COL = "nitrate_uM"
LAB_STATION_COL = "station"    # optional; copied through if present
LAB_FLAG_COL    = "flag"       # optional
# ==========================================================================

_MONTHS = {m: i for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"], start=1)}
_TIME_RE = re.compile(r"\b(\d{1,2}):(\d{2}):(\d{2})\b")


def _find_col(names: List[str], *cands: str) -> Optional[int]:
    low = [n.strip().lower() for n in names]
    for c in cands:
        if c in low:
            return low.index(c)
    for c in cands:
        for i, n in enumerate(low):
            if n.startswith(c):
                return i
    return None


def _btl_column_names(raw_lines: List[str]) -> List[str]:
    """Variable names from the .btl two-row title header.

    SBE Bottle Summary writes columns as a title row, not '# name N =' lines, e.g.:
        '    Bottle        Date      Sal00       PrDM     Tv290C      C0S/m    Sbeox0V'
        '  Position        Time'
    The variable names are the tokens after 'Bottle' and 'Date'."""
    for ln in raw_lines:
        toks = ln.split()
        if len(toks) >= 3 and toks[0].lower() == "bottle" and toks[1].lower() == "date":
            return toks[2:]
    return []


def parse_btl(path: Path) -> Tuple[List[str], List[dict]]:
    """Parse one SBE Bottle Summary .btl into (column_names, list-of-bottle-dicts).

    Handles the standard 2-line-per-bottle layout: the mean line begins with the bottle
    number followed by 'Mon DD YYYY' and the per-variable means; the next line carries
    'HH:MM:SS' and the standard deviations. Column identities come from the 'Bottle ...
    Date ...' title row (SBE Bottle Summary format), falling back to any '# name N ='
    header lines a different SBE build might inherit from the source .cnv."""
    raw = path.read_text(encoding="latin-1", errors="replace").splitlines()
    names = _btl_column_names(raw)
    if not names:
        header = [ln for ln in raw if ln.strip().startswith("#") or ln.strip().startswith("*")]
        names = L.cnv_column_names(header)
    if not names:
        raise ValueError(f"{path.name}: could not identify columns "
                         "(no 'Bottle ... Date' title row and no '# name N =' lines).")
    nq = len(names)
    pcol = L.cnv_pressure_column(names)
    tcol = _find_col(names, "t090c", "tv290c", "t090", "tv290", "t68", "t090cm")
    scol = _find_col(names, "sal00", "sal11", "sal78", "sal")
    ccol = _find_col(names, "c0s/m", "c1s/m", "c0ms/cm", "c1ms/cm", "cond0s/m", "cond", "c0")
    # conductivity unit: c0S/m is S/m (what practical_salinity wants); c0mS/cm is mS/cm (/10 -> S/m)
    c_to_Sm = 0.1 if (ccol is not None and "ms/cm" in names[ccol].strip().lower()) else 1.0

    records: List[dict] = []
    for idx, ln in enumerate(raw):
        s = ln.strip()
        if not s or s.startswith("#") or s.startswith("*"):
            continue
        toks = s.split()
        if not toks[0].lstrip("-").isdigit():
            continue
        bottle_no = int(toks[0])
        # detect an optional 'Mon DD YYYY' date immediately after the bottle number
        ndate = 0
        date_iso = ""
        if len(toks) >= 4 and toks[1][:3].lower() in _MONTHS and toks[2].isdigit() and toks[3].isdigit():
            ndate = 3
            mon = _MONTHS[toks[1][:3].lower()]
            day, year = int(toks[2]), int(toks[3])
            date_iso = f"{year:04d}-{mon:02d}-{day:02d}"
        means = toks[1 + ndate: 1 + ndate + nq]
        if len(means) < nq:
            continue
        try:
            vals = [float(v) for v in means]
        except ValueError:
            continue
        # time from the following (std-dev) line, if present
        time_str = ""
        if idx + 1 < len(raw):
            m = _TIME_RE.search(raw[idx + 1])
            if m:
                time_str = m.group(0)
        time_utc = f"{date_iso}T{time_str}" if (date_iso and time_str) else ""

        press = vals[pcol] if pcol is not None else float("nan")
        temp = vals[tcol] if tcol is not None else float("nan")
        # salinity: use the .btl column if present; otherwise compute from averaged C, T, P
        # (PSS-78 via ctd_lib) so a Bottle Summary without a derived-salinity column still works
        if scol is not None:
            sal = vals[scol]
            sal_src = "btl"
        elif ccol is not None and tcol is not None and pcol is not None and press == press and temp == temp:
            sal = float(L.practical_salinity(vals[ccol] * c_to_Sm, temp, press))
            sal_src = "computed"
        else:
            sal, sal_src = float("nan"), ""
        rec = {
            "bottle_no": bottle_no,
            "pressure_dbar": press,
            "depth_m": round(float(L.pressure_to_depth(press, LATITUDE)), 2) if press == press else "",
            "time_utc": time_utc,
            "temp_c": round(temp, 4) if temp == temp else "",
            "salinity": round(sal, 4) if sal == sal else "",
            "salinity_source": sal_src,
        }
        records.append(rec)
    return names, records


def build_ctd_bottle_table(bottle_dir: Path) -> pd.DataFrame:
    btls = sorted(p for p in bottle_dir.glob("*.btl") if p.is_file())
    if not btls:
        btls = sorted(p for p in bottle_dir.glob("*") if p.is_file() and p.suffix.lower() == ".btl")
    if not btls:
        raise FileNotFoundError(f"No .btl files in {bottle_dir}")
    rows = []
    for b in btls:
        cast_id = b.stem
        try:
            _, recs = parse_btl(b)
        except Exception as exc:
            print(f"  {b.name}: SKIPPED ({exc})")
            continue
        for r in recs:
            r = {"cast_id": cast_id, **r}
            rows.append(r)
        print(f"  {b.name}: {len(recs)} bottle(s)")
    df = pd.DataFrame(rows, columns=["cast_id", "bottle_no", "pressure_dbar", "depth_m",
                                     "time_utc", "temp_c", "salinity", "salinity_source"])
    return df


def read_lab(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        return pd.read_excel(path, sheet_name=LAB_SHEET)
    return pd.read_csv(path)


def main() -> int:
    global BOTTLE_DIR
    if len(sys.argv) > 1:
        BOTTLE_DIR = Path(sys.argv[1])
    if not BOTTLE_DIR.exists():
        print(f"BOTTLE_DIR does not exist:\n  {BOTTLE_DIR}")
        return 1
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Parsing .btl files in: {BOTTLE_DIR}")
    ctd = build_ctd_bottle_table(BOTTLE_DIR)
    ctd_path = OUTPUT_DIR / f"{CRUISE_ID}_bottle_ctd.csv"
    ctd.to_csv(ctd_path, index=False, encoding="utf-8-sig")
    print(f"\nCTD bottle table: {len(ctd)} rows -> {ctd_path}")

    if LAB_FILE is None:
        # Emit a lab-nitrate template pre-keyed with every bottle, so whoever runs the
        # nitrate fills values into a sheet whose (cast_id, bottle_no) keys already match
        # the .btl exactly -> the join later matches 100%. depth_m/time_utc are carried
        # through as read-only context to help identify each sample; only nitrate_uM (and
        # optionally replicate_id/station/flag) need filling.
        template = ctd[["cast_id", "bottle_no", "depth_m", "time_utc"]].copy()
        template["nitrate_uM"] = ""
        template["replicate_id"] = ""
        template["station"] = ""
        template["flag"] = ""
        tmpl_path = OUTPUT_DIR / f"{CRUISE_ID}_lab_nitrate_TEMPLATE.csv"
        template.to_csv(tmpl_path, index=False, encoding="utf-8-sig")
        print(f"\nLab-nitrate template ({len(template)} bottles) -> {tmpl_path}")
        print("Fill the nitrate_uM column (and replicate_id/station/flag if used), save it, then set "
              "LAB_FILE to that file and rerun to build the step03 bottle_nitrate table.")
        print("\nNo LAB_FILE set — CTD-side table only for now.")
        return 0

    lab_path = Path(LAB_FILE)
    if not lab_path.exists():
        print(f"LAB_FILE does not exist:\n  {lab_path}")
        return 1
    lab = read_lab(lab_path)
    missing = [c for c in (LAB_CAST_COL, LAB_BOTTLE_COL, LAB_NITRATE_COL) if c not in lab.columns]
    if missing:
        print(f"Lab file is missing required column(s): {missing}")
        print(f"Columns present: {list(lab.columns)}")
        print("Set LAB_CAST_COL / LAB_BOTTLE_COL / LAB_NITRATE_COL to match, or add a cast_id "
              "+ bottle_no key to the lab sheet (station/layer-keyed sheets need your crosswalk first).")
        return 1

    lab = lab.copy()
    lab["_cast"] = lab[LAB_CAST_COL].astype(str).str.strip()
    lab["_bot"] = pd.to_numeric(lab[LAB_BOTTLE_COL], errors="coerce").astype("Int64")
    ctd["_cast"] = ctd["cast_id"].astype(str).str.strip()
    ctd["_bot"] = pd.to_numeric(ctd["bottle_no"], errors="coerce").astype("Int64")

    merged = ctd.merge(lab, on=["_cast", "_bot"], how="left", suffixes=("", "_lab"))
    out = pd.DataFrame({
        "cruise_id": CRUISE_ID,
        "cast_id": merged["cast_id"],
        "station": merged[LAB_STATION_COL] if LAB_STATION_COL in merged.columns else "",
        "depth_m": merged["depth_m"],
        "time_utc": merged["time_utc"],
        "nitrate_uM": merged[LAB_NITRATE_COL],
        "replicate_id": "",
        "flag": merged[LAB_FLAG_COL] if LAB_FLAG_COL in merged.columns else "",
    })
    matched = out["nitrate_uM"].notna()
    nit_path = OUTPUT_DIR / f"{CRUISE_ID}_bottle_nitrate.csv"
    out.to_csv(nit_path, index=False, encoding="utf-8-sig")
    print(f"Bottle nitrate table (step03): {int(matched.sum())} of {len(out)} bottles matched to lab "
          f"-> {nit_path}")

    # unmatched report (both directions)
    lab_keys = set(zip(lab["_cast"], lab["_bot"].astype("object")))
    ctd_keys = set(zip(ctd["_cast"], ctd["_bot"].astype("object")))
    unmatched_bottles = ctd[~ctd.apply(lambda r: (r["_cast"], r["_bot"]) in lab_keys, axis=1)]
    unmatched_lab = lab[~lab.apply(lambda r: (r["_cast"], r["_bot"]) in ctd_keys, axis=1)]
    if len(unmatched_bottles) or len(unmatched_lab):
        rep = OUTPUT_DIR / f"{CRUISE_ID}_bottle_nitrate_UNMATCHED.csv"
        with open(rep, "w", encoding="utf-8-sig") as fh:
            fh.write("# Bottles with no lab nitrate (cast_id,bottle_no):\n")
            for _, r in unmatched_bottles.iterrows():
                fh.write(f"{r['cast_id']},{r['bottle_no']}\n")
            fh.write("# Lab rows with no matching .btl bottle (cast,bottle):\n")
            for _, r in unmatched_lab.iterrows():
                fh.write(f"{r['_cast']},{r['_bot']}\n")
        print(f"Unmatched report -> {rep}  "
              f"({len(unmatched_bottles)} bottles without lab, {len(unmatched_lab)} lab rows without bottle)")
    print("\nReview the nitrate table, add replicate_id where you have replicates, then run step03.")
    return 0


if __name__ == "__main__":
    code = main()
    try:
        input("\nPress Enter to close...")
    except (EOFError, KeyboardInterrupt):
        pass
    sys.exit(code)
