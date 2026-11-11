"""Depth from pressure (UNESCO 1983, Fofonoff & Millard). Moved verbatim from ctd_lib."""
from __future__ import annotations
import numpy as np


def pressure_to_depth(p_dbar, lat_deg=0.0):
    """Depth [m] from pressure [dbar] at latitude (UNESCO 1983, Fofonoff & Millard)."""
    p = np.asarray(p_dbar, float)
    x = np.sin(np.deg2rad(lat_deg)) ** 2
    g = 9.780318 * (1.0 + 5.2788e-3 * x + 2.36e-5 * x * x) + 1.092e-6 * p
    num = (((-1.82e-15 * p + 2.279e-10) * p - 2.2512e-5) * p + 9.72659) * p
    return num / g
