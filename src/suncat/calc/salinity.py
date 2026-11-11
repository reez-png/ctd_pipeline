"""Practical salinity (PSS-78). Moved verbatim from ctd_lib; ctd_lib now re-exports these."""
from __future__ import annotations
import numpy as np

C_STD = 4.2914   # conductivity of standard seawater C(35,15,0) in S/m


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
