"""
ctd_lib.py  -  Tested core functions for the CTD + SUNA pipeline.

Pure, importable, side-effect-free helpers shared by the standalone tools and
covered by tests/test_ctd_lib.py. Keeping the reference implementations here (rather
than only inside the flat step01 script) is what lets a regression be caught before
it reaches data.

Functions
  read_cnv(path)                 -> (column_names, data ndarray)  parse an SBE .cnv
  cnv_pressure_column(names)     -> int | None                    index of the pressure column
  cnv_sample_rate(header_lines)  -> float                         Hz from '# interval = seconds:'
  descent_onset(pressure, fs)    -> int                           first scan of the real descent
  section_cnv(in_path, out_path) -> dict                          write soak+dip-removed .cnv (CRLF-safe)
  sal78(R, T68, P)               -> practical salinity (PSS-78)
  practical_salinity(C, T90, P)  -> practical salinity from conductivity [S/m]
  potential_temperature(S,T90,P) -> potential temperature (UNESCO/Fofonoff, ref 0 dbar)
  pressure_to_depth(p, lat)      -> depth [m] from pressure [dbar] (UNESCO 1983)

Standard library + numpy/pandas only.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

C_STD = 4.2914   # conductivity of standard seawater C(35,15,0) in S/m
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

    Walk back from the deepest scan while pressure keeps rising over a ~1 s window;
    fall back to the surface crossing if that is degenerate (a near-bottom stall)."""
    p = np.asarray(pressure, float)
    if p.size == 0 or np.isnan(p).all():
        return 0
    imax = int(np.nanargmax(p))
    win = max(int(round(fs)), 1)
    if imax < win + 1:
        return 0
    seg = pd.Series(p[:imax + 1]).rolling(win, center=True, min_periods=1).median().to_numpy()
    fwd = np.full(len(seg), np.nan)
    fwd[:len(seg) - win] = seg[win:] - seg[:len(seg) - win]
    descending = np.where(np.isnan(fwd), True, fwd > _SOAK_MIN_DROP_DB)
    i = imax
    while i > 0 and descending[i - 1]:
        i -= 1
    start = i
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


# --------------------------------------------------------------------------- #
# Seawater: PSS-78 salinity, potential temperature, depth
# --------------------------------------------------------------------------- #
def sal78(R, T68, P):
    """Practical salinity (PSS-78) from conductivity ratio R = C/C(35,15,0)."""
    a = [0.0080, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081]
    b = [0.0005, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144]
    c = [0.6766097, 2.00564e-2, 1.104259e-4, -6.9698e-7, 1.0031e-9]
    d = [3.426e-2, 4.464e-4, 4.215e-1, -3.107e-3]
    e = [2.070e-5, -6.370e-10, 3.989e-15]
    R = np.asarray(R, float)
    Rp = 1 + (P * (e[0] + e[1] * P + e[2] * P * P)) / (1 + d[0] * T68 + d[1] * T68 * T68 + (d[2] + d[3] * T68) * R)
    rt = c[0] + c[1] * T68 + c[2] * T68**2 + c[3] * T68**3 + c[4] * T68**4
    Rt = R / (Rp * rt)
    Rt = np.maximum(Rt, 1e-12)
    s = np.sqrt(Rt)
    dS = (T68 - 15) / (1 + 0.0162 * (T68 - 15)) * (b[0] + b[1] * s + b[2] * Rt + b[3] * Rt * s + b[4] * Rt**2 + b[5] * Rt**2 * s)
    return a[0] + a[1] * s + a[2] * Rt + a[3] * Rt * s + a[4] * Rt**2 + a[5] * Rt**2 * s + dS


def practical_salinity(C_Sm, T90, P):
    """Practical salinity from conductivity [S/m], temperature [ITS-90 degC], pressure [dbar]."""
    return sal78(np.asarray(C_Sm, float) / C_STD, np.asarray(T90, float) * 1.00024, P)


def _adtg(S, T, P):
    a0, a1, a2, a3 = 3.5803e-5, 8.5258e-6, -6.836e-8, 6.6228e-10
    b0, b1 = 1.8932e-6, -4.2393e-8
    c0, c1, c2, c3 = 1.8741e-8, -6.7795e-10, 8.733e-12, -5.4481e-14
    d0, d1 = -1.1351e-10, 2.7759e-12
    e0, e1, e2 = -4.6206e-13, 1.8676e-14, -2.1687e-16
    return (a0 + (a1 + (a2 + a3 * T) * T) * T + (b0 + b1 * T) * (S - 35)
            + ((c0 + (c1 + (c2 + c3 * T) * T) * T) + (d0 + d1 * T) * (S - 35)) * P + (e0 + (e1 + e2 * T) * T) * P * P)


def potential_temperature(S, T90, P, PR=0.0):
    """Potential temperature (UNESCO/Fofonoff RK4) referenced to PR dbar (default 0)."""
    S = np.asarray(S, float)
    T = np.asarray(T90, float) * 1.00024
    P = np.asarray(P, float)
    H = PR - P
    XK = H * _adtg(S, T, P); T = T + 0.5 * XK; Q = XK; P2 = P + 0.5 * H
    XK = H * _adtg(S, T, P2); T = T + 0.29289322 * (XK - Q); Q = 0.58578644 * XK + 0.121320344 * Q
    XK = H * _adtg(S, T, P2); T = T + 1.707106781 * (XK - Q); Q = 3.414213562 * XK - 4.121320344 * Q
    P3 = P2 + 0.5 * H
    XK = H * _adtg(S, T, P3); T = T + (XK - 2.0 * Q) / 6.0
    return T / 1.00024


def pressure_to_depth(p_dbar, lat_deg=0.0):
    """Depth [m] from pressure [dbar] at latitude (UNESCO 1983, Fofonoff & Millard)."""
    p = np.asarray(p_dbar, float)
    x = np.sin(np.deg2rad(lat_deg)) ** 2
    g = 9.780318 * (1.0 + 5.2788e-3 * x + 2.36e-5 * x * x) + 1.092e-6 * p
    num = (((-1.82e-15 * p + 2.279e-10) * p - 2.2512e-5) * p + 9.72659) * p
    return num / g
