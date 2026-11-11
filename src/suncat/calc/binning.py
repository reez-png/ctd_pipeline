"""Pressure / x binning: mean of y within equal-width bins (bin-center convention)."""
from __future__ import annotations
import numpy as np
import pandas as pd


def bin_mean(x, y, step):
    """Mean of y within x-bins of width `step`, keyed by bin centre (np.rint, round-half-to-even).
    Returns (bin_centres, means) as float arrays, empty arrays when there is no finite data.

    Consolidates step00b's bin_by_pressure and step00c's bin_profile; step00c's binned_xy is the same
    binning plus the integer bin ids.
    """
    d = pd.DataFrame({"x": x, "y": y}).apply(pd.to_numeric, errors="coerce").dropna()
    if d.empty:
        return np.array([]), np.array([])
    ids = np.rint(d["x"].to_numpy(float) / step).astype(int)
    g = d.assign(_b=ids).groupby("_b")["y"].mean()
    return g.index.to_numpy(float) * step, g.to_numpy(float)
