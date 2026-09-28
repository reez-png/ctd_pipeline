"""
verify_reproducibility.py

Check that re-running the step01 pipeline on the same .hex files reproduces the
same numbers. Sea-Bird .cnv headers carry run timestamps that ALWAYS differ, so a
byte-for-byte compare is meaningless. This compares the DATA COLUMNS instead:
same variables, same row count, and per-variable maximum absolute difference.

HOW TO USE
  1. Run step01 once. Copy its output tree somewhere as a baseline, e.g.
       C:\\Projects\\ctd_pipeline\\cruises\\P45_06\\L2   ->   ..._baseline
     (or just copy the specific product folders you care about).
  2. Run step01 again (same .hex, same PSAs).
  3. Point BASELINE_DIR and RERUN_DIR below at the two trees and run this file.
     Files are matched by identical filename found under each tree.

  A deterministic pipeline should give max abs diff = 0.0 for every variable.
  Non-zero differences mean something changed (a PSA edit, a setting, a code
  change, or non-determinism) and are listed in the report.

Exit code is 0 when every matched file reproduces within TOLERANCE, else 1.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ===========================================================================
# USER SETTINGS  (edit these two, or pass them on the command line)
# ===========================================================================
BASELINE_DIR = Path(r"C:\Projects\ctd_pipeline\cruises\P45_06\L2_baseline")
RERUN_DIR    = Path(r"C:\Projects\ctd_pipeline\cruises\P45_06\L2")
GLOB_PATTERN = "*.cnv"          # which products to compare (e.g. "*_1m_down.cnv")
TOLERANCE    = 1e-9             # max abs diff allowed before a variable is "changed"
REPORT_CSV   = Path(r"C:\Projects\ctd_pipeline\cruises\P45_06\_audit\reproducibility_report.csv")

# Sea-Bird bad-flag sentinel written into the flag column for marked scans.
SBE_BADFLAG = -9.990e-29


def read_cnv_data(path: Path) -> Tuple[List[str], np.ndarray]:
    """Return (column_names, data_array) from an SBE .cnv, ignoring the header."""
    lines = Path(path).read_text(encoding="latin-1", errors="replace").splitlines()
    name_re = re.compile(r"#\s*name\s+(\d+)\s*=\s*([^:]+?)\s*[:=]", re.I)
    names: Dict[int, str] = {}
    end = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("#") or s.startswith("*"):
            m = name_re.search(ln)
            if m:
                names[int(m.group(1))] = m.group(2).strip()
            if s.upper().startswith("*END*"):
                end = i + 1
                break
    if end is None:
        raise ValueError(f"No *END* marker in {path.name}")
    cols = [names[k] for k in sorted(names)]
    rows = []
    for ln in lines[end:]:
        if not ln.strip():
            continue
        parts = ln.split()
        if len(parts) >= len(cols):
            rows.append(parts[: len(cols)])
    if not rows:
        return cols, np.empty((0, len(cols)))
    data = pd.DataFrame(rows, columns=cols).apply(pd.to_numeric, errors="coerce").to_numpy(float)
    return cols, data


def count_flags(cols: List[str], data: np.ndarray) -> Optional[int]:
    """Number of Sea-Bird bad-flagged scans. Good scans write flag 0.0; a marked
    scan writes the bad-flag sentinel (default -9.99e-29), so any non-zero flag
    counts as bad."""
    flag_idx = next((i for i, c in enumerate(cols) if c.lower() == "flag"), None)
    if flag_idx is None or data.size == 0:
        return None
    col = data[:, flag_idx]
    return int(np.sum(np.nan_to_num(col) != 0.0))


def compare_pair(base: Path, rerun: Path, tol: float = TOLERANCE) -> Dict[str, object]:
    """Compare one matched file pair by data columns."""
    cb, db = read_cnv_data(base)
    cr, dr = read_cnv_data(rerun)

    rec: Dict[str, object] = {
        "file": base.name,
        "rows_baseline": db.shape[0],
        "rows_rerun": dr.shape[0],
        "cols_match": cb == cr,
        "flags_baseline": count_flags(cb, db),
        "flags_rerun": count_flags(cr, dr),
        "max_abs_diff": "",
        "worst_variable": "",
        "identical": False,
        "note": "",
    }

    if cb != cr:
        rec["note"] = "column names differ"
        return rec
    if db.shape != dr.shape:
        rec["note"] = "row count differs"
        return rec
    if db.size == 0:
        rec["identical"] = bool(dr.size == 0)
        rec["max_abs_diff"] = 0.0
        return rec

    diff = np.abs(db - dr)
    # NaNs in the same place are fine; NaN vs number is a real difference.
    both_nan = np.isnan(db) & np.isnan(dr)
    one_nan = np.isnan(db) ^ np.isnan(dr)
    diff[both_nan] = 0.0
    diff[one_nan] = np.inf
    per_col_max = np.nanmax(diff, axis=0) if diff.size else np.array([0.0])
    worst_i = int(np.argmax(per_col_max))
    rec["max_abs_diff"] = float(per_col_max[worst_i])
    rec["worst_variable"] = cb[worst_i]
    rec["identical"] = bool(rec["max_abs_diff"] <= tol)
    if not rec["identical"]:
        rec["note"] = "data differs beyond tolerance"
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description="Compare two step01 output trees by .cnv data columns.")
    ap.add_argument("--baseline", type=Path, default=BASELINE_DIR)
    ap.add_argument("--rerun", type=Path, default=RERUN_DIR)
    ap.add_argument("--glob", default=GLOB_PATTERN)
    ap.add_argument("--tol", type=float, default=TOLERANCE)
    ap.add_argument("--report", type=Path, default=REPORT_CSV)
    args = ap.parse_args()

    for label, d in (("baseline", args.baseline), ("rerun", args.rerun)):
        if not d.exists():
            print(f"ERROR: {label} folder does not exist:\n  {d}")
            return 2

    base_files = {p.name: p for p in sorted(args.baseline.rglob(args.glob))}
    rerun_files = {p.name: p for p in sorted(args.rerun.rglob(args.glob))}
    common = sorted(set(base_files) & set(rerun_files))
    only_base = sorted(set(base_files) - set(rerun_files))
    only_rerun = sorted(set(rerun_files) - set(base_files))

    if not common:
        print("No files with matching names were found under both folders.")
        print(f"  baseline: {len(base_files)} files matching {args.glob}")
        print(f"  rerun:    {len(rerun_files)} files matching {args.glob}")
        return 2

    records = [compare_pair(base_files[name], rerun_files[name], tol=args.tol) for name in common]
    report = pd.DataFrame(records)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.report, index=False, encoding="utf-8-sig")

    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(report.to_string(index=False))

    n_ok = int(report["identical"].sum())
    n_bad = len(report) - n_ok
    print()
    print(f"Matched files compared : {len(common)}")
    print(f"Reproduced (diff <= {args.tol:g}) : {n_ok}")
    print(f"Changed                : {n_bad}")
    if only_base:
        print(f"Only in baseline ({len(only_base)}): {', '.join(only_base[:8])}{' ...' if len(only_base) > 8 else ''}")
    if only_rerun:
        print(f"Only in rerun ({len(only_rerun)}): {', '.join(only_rerun[:8])}{' ...' if len(only_rerun) > 8 else ''}")
    print(f"\nReport written to: {args.report}")

    if n_bad == 0 and not only_base and not only_rerun:
        print("\nREPRODUCIBLE: every matched product is numerically identical.")
        return 0
    print("\nNOT fully reproducible: see the rows above where identical = False.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
