"""Canonical Sea-Bird .cnv reader for suncat.

Replaces the three near-identical parsers copied across the QC notebooks, preserving their exact
behaviour so the frozen regression baseline still holds:

  - step00a  read_cnv        -> read_cnv(path)                     (plain parse, all scans)
  - step00b  parse_cnv_raw   -> read_cnv(path, mask_badflag=True)  (NaN the Sea-Bird bad-flag sentinel)
  - step00b  read_cnv_full   -> read_cnv_full(path)                (header + raw rows, for re-writing)
"""
from __future__ import annotations
import re
from pathlib import Path
from typing import NamedTuple
import pandas as pd

# Column names come from header lines like '# name 0 = prdM: Pressure, Strain Gauge [db]'.
_NAME_RE = re.compile(r"#\s*name\s+(\d+)\s*=\s*([^:]+?)\s*[:=]", re.I)

# Sea-Bird bad-flag sentinel is -9.990e-29, a TINY negative number written into record ends by shift
# operations (Align CTD, Cell TM). Bracket it tightly; real T/C/P/S/O2 values never fall here.
_BADFLAG_LO, _BADFLAG_HI = -1.2e-28, -8e-29
# The loop-edit 'flag' column legitimately uses that value to mark a flagged scan, so never mask it.
_FLAG_COLS = ("flag",)


class CnvFile(NamedTuple):
    """Full parse of a .cnv. Unpacks exactly like the old (header, cols, data, rows_raw) tuple,
    and also exposes the fields by name."""
    header: list
    columns: list
    data: pd.DataFrame
    rows_raw: list


def _read_lines(path):
    return Path(path).read_text(encoding="latin-1", errors="replace").splitlines()


def _scan_header(lines):
    """Return (ordered column names, index of the first data line). Shared by every reader:
    names come from the '# name N = short: long' lines, data starts after '*END*'."""
    names, data_start = {}, None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("#") or s.startswith("*"):
            m = _NAME_RE.search(ln)
            if m:
                names[int(m.group(1))] = m.group(2).strip()
            if s.upper().startswith("*END*"):
                data_start = i + 1
                break
    cols = [names[k] for k in sorted(names)]
    return cols, data_start


def read_cnv(path, mask_badflag: bool = False):
    """Parse a .cnv into a DataFrame, keeping every scan (including loop-edit flagged ones).

    mask_badflag=False reproduces the plain step00a reader.
    mask_badflag=True reproduces step00b parse_cnv_raw: the bad-flag sentinel becomes NaN so it
    never plots or enters a statistic, while the 'flag' column keeps its sentinel.
    """
    if path is None:
        return None
    lines = _read_lines(path)
    cols, ds = _scan_header(lines)
    if ds is None:
        return None
    rows = [l.split()[:len(cols)] for l in lines[ds:] if l.strip() and len(l.split()) >= len(cols)]
    df = pd.DataFrame(rows, columns=cols).apply(pd.to_numeric, errors="coerce").reset_index(drop=True)
    if mask_badflag:
        bad = (df.le(_BADFLAG_HI) & df.ge(_BADFLAG_LO)) | (df.abs() > 1e30)
        for fc in _FLAG_COLS:
            if fc in bad.columns:
                bad[fc] = False
        df = df.mask(bad)
    return df


def read_cnv_full(path):
    """Parse a .cnv keeping the original header lines and raw data lines, so a sectioned copy can be
    written back out with its header intact. Reproduces step00b read_cnv_full. Returns a CnvFile."""
    if path is None:
        return None
    lines = _read_lines(path)
    cols, ds = _scan_header(lines)
    if ds is None:
        return None
    header = lines[:ds]
    rows_raw = [l for l in lines[ds:] if l.strip()]
    data = pd.DataFrame([r.split()[:len(cols)] for r in rows_raw if len(r.split()) >= len(cols)],
                        columns=cols).apply(pd.to_numeric, errors="coerce")
    return CnvFile(header=header, columns=cols, data=data, rows_raw=rows_raw)
