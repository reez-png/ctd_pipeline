"""
ctd_lib.py  -  Tested core functions for the CTD + SUNA pipeline.

Pure, importable, side-effect-free helpers shared by the standalone tools and covered by
tests/test_ctd_lib.py. The seawater routines (sal78, practical_salinity, potential_temperature,
pressure_to_depth) now live in suncat.calc and are RE-EXPORTED here, so every existing
`from ctd_lib import ...` keeps working unchanged. The cnv parsing, descent onset and sectioning
helpers below remain here for now (they are reconciled with suncat.io.cnv / suncat.section in a
later step).

Functions
  read_cnv(path)                 -> (column_names, data ndarray)  parse an SBE .cnv
  cnv_pressure_column(names)     -> int | None                    index of the pressure column
  cnv_sample_rate(header_lines)  -> float                         Hz from '# interval = seconds:'
  descent_onset(pressure, fs)    -> int                           first scan of the real descent
  section_cnv(in_path, out_path) -> dict                          write soak+dip-removed .cnv (CRLF-safe)
  sal78, practical_salinity, potential_temperature, pressure_to_depth   (re-exported from suncat.calc)

Standard library + numpy/pandas only.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Seawater science now lives in suncat.calc and is re-exported so existing imports keep working.
from suncat.calc.salinity import C_STD, sal78, practical_salinity
from suncat.calc.thermo import potential_temperature
from suncat.calc.depth import pressure_to_depth
from suncat.section import descent_start

_FS_DEFAULT = 24.0
_SOAK_MIN_DROP_DB = 0.10
_SURF_SKIP_DB = 8.0

_NAME_RE = re.compile(r"#\s*name\s+(\d+)\s*=\s*([^:]+?)\s*[:=]", re.I)


# --------------------------------------------------------------------------- #
# .cnv parsing
# --------------------------------------------------------------------------- #
def cnv_column_names(header_lines: List[str]) -> List[str]:
    """Ordered short column names from '# name N = short: description' lines."""
    names: Dict[int, str] = {}
    for ln in header_lines:
        m = _NAME_RE.search(ln)
        if m:
            names[int(m.group(1))] = m.group(2).strip()
    return [names[k] for k in sorted(names)]


def cnv_pressure_column(names: List[str]) -> Optional[int]:
    """0-based index of the pressure column (prDM / prM / prSM / pressure)."""
    for i, c in enumerate(names):
        cl = c.strip().lower()
        if cl in ("prdm", "prm", "prsm") or cl.startswith("prdm"):
            return i
    for i, c in enumerate(names):
        cl = c.strip().lower()
        if "pr" in cl and "spar" not in cl and "par" != cl:
            return i
    return None


def cnv_sample_rate(header_lines: List[str]) -> float:
    """Sample rate (Hz) from '# interval = seconds: X'; fallback to 24 Hz."""
    for ln in header_lines:
        low = ln.lower()
        if "interval" in low and "seconds" in low:
            m = re.search(r"seconds\s*:\s*([0-9.eE+-]+)", ln)
            if m:
                try:
                    dt = float(m.group(1))
                    if dt > 0:
                        return 1.0 / dt
                except ValueError:
                    pass
    return _FS_DEFAULT


def read_cnv(path) -> Tuple[List[str], np.ndarray]:
    """Return (column_names, data_array) from an SBE .cnv. Empty array if no data."""
    lines = Path(path).read_text(encoding="latin-1", errors="replace").splitlines()
    end = None
    for i, ln in enumerate(lines):
        if ln.strip().upper().startswith("*END*"):
            end = i + 1
            break
    if end is None:
        raise ValueError(f"No *END* marker in {Path(path).name}; not a readable .cnv")
    names = cnv_column_names(lines[:end])
    rows = []
    for ln in lines[end:]:
        if not ln.strip():
            continue
        parts = ln.split()
        if names and len(parts) >= len(names):
            rows.append(parts[:len(names)])
    if not rows:
        return names, np.empty((0, len(names)))
    data = pd.DataFrame(rows, columns=names).apply(pd.to_numeric, errors="coerce").to_numpy(float)
    return names, data


# --------------------------------------------------------------------------- #
# Descent onset + sectioning (soak + equilibration dip removal)
# --------------------------------------------------------------------------- #
def descent_onset(pressure, fs: float = _FS_DEFAULT) -> int:
    """Index of the last near-surface scan before the single continuous descent.

    The core walk-back is suncat.section.descent_start; this adds the near-bottom degenerate
    fallback (fall back to the surface crossing when the detected descent is too short or shallow).
    Verified byte-identical to the previous inline implementation across 600+ profiles."""
    p = np.asarray(pressure, float)
    if p.size == 0 or np.isnan(p).all():
        return 0
    imax = int(np.nanargmax(p))
    win = max(int(round(fs)), 1)
    if imax < win + 1:
        return 0
    start, _ = descent_start(p, win, _SOAK_MIN_DROP_DB)
    span_db = (p[imax] - p[start]) if imax > start else 0.0
    if (imax - start) < int(5 * fs) or span_db < 10.0:
        below = np.where(p[:imax + 1] < _SURF_SKIP_DB)[0]
        start = int(below[-1]) if below.size else 0
    return start
def section_cnv(input_cnv, output_cnv) -> Dict[str, object]:
    """Write a sectioned copy of an SBE .cnv: keep [descent onset : end], dropping the
    leading soak + equilibration dip. Byte-faithful (CRLF preserved); only '# nvalues'
    is updated. Returns a summary dict."""
    input_cnv, output_cnv = Path(input_cnv), Path(output_cnv)
    output_cnv.parent.mkdir(parents=True, exist_ok=True)
    raw = input_cnv.read_bytes()
    lines = raw.splitlines(keepends=True)   # bytes, original CRLF preserved

    end_idx = None
    for i, ln in enumerate(lines):
        if ln.strip().upper().startswith(b"*END*"):
            end_idx = i
            break
    if end_idx is None:
        raise ValueError(f"No *END* marker in {input_cnv.name}; not a readable .cnv")

    header = lines[:end_idx + 1]
    data_rows = [ln for ln in lines[end_idx + 1:] if ln.strip() != b""]
    n_in = len(data_rows)
    if n_in == 0:
        raise ValueError(f"{input_cnv.name} has no data rows after *END*")

    header_text = [ln.decode("latin-1", "replace").rstrip("\r\n") for ln in header]
    pcol = cnv_pressure_column(cnv_column_names(header_text))
    if pcol is None:
        raise ValueError(f"No pressure column found in {input_cnv.name}")
    fs = cnv_sample_rate(header_text)

    pressure = np.full(n_in, np.nan)
    for k, ln in enumerate(data_rows):
        parts = ln.decode("latin-1", "replace").split()
        if pcol < len(parts):
            try:
                pressure[k] = float(parts[pcol])
            except ValueError:
                pressure[k] = np.nan

    start = descent_onset(pressure, fs)
    kept = data_rows[start:]
    n_out = len(kept)

    nval_re = re.compile(rb"(#\s*nvalues\s*=\s*)(\d+)", re.I)
    new_header = [nval_re.sub(lambda m: m.group(1) + str(n_out).encode("ascii"), ln)
                  if nval_re.search(ln) else ln for ln in header]

    with open(output_cnv, "wb") as fh:
        fh.writelines(new_header)
        fh.writelines(kept)

    return {"pressure_col": pcol, "sample_rate_hz": round(fs, 4),
            "rows_in": n_in, "rows_kept": n_out, "rows_removed": n_in - n_out,
            "onset_scan": start}
