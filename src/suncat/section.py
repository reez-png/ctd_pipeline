"""Motion-based cast sectioning for suncat.

Finds the descent onset by motion, not a fixed depth line, so the surface soak, the equilibration
dip and any flat plateau are excluded while the full descent (including its shallow top) is kept.
Consolidates the byte-identical _descent_start / descent_start used in step00b and step00c; each
notebook keeps a thin wrapper that sets the smoothing window its own way (scans directly, or
seconds x sample-rate).
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def descent_start(pv, win, min_drop: float = 0.10):
    """Positional (start, imax) of the true descent for a pressure array `pv`.

    Walking back from the deepest scan, keep every scan still genuinely descending (smoothed pressure
    rising by more than `min_drop` dbar over the next `win` scans), so the soak, dip and any flat
    plateau are dropped while the full descent is kept. `win` is the smoothing window in SCANS.
    start = descent onset; imax = deepest scan.
    """
    win = int(win)
    pv = np.asarray(pv, float)
    if pv.size == 0 or np.isnan(pv).all():
        return 0, 0
    imax = int(np.nanargmax(pv))
    if imax < win + 1:
        return 0, imax
    seg = pv[:imax + 1]                                   # downcast side only (never peek upcast)
    s = pd.Series(seg).rolling(win, center=True, min_periods=1).median().to_numpy()
    fwd = np.full(len(s), np.nan)
    fwd[:len(s) - win] = s[win:] - s[:len(s) - win]       # pressure gain over the next `win` scans
    descending = np.where(np.isnan(fwd), True, fwd > min_drop)
    i = imax
    while i > 0 and descending[i - 1]:
        i -= 1
    return i, imax
