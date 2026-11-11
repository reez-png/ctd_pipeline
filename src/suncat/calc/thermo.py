"""Potential temperature (UNESCO/Fofonoff RK4). Moved verbatim from ctd_lib."""
from __future__ import annotations
import numpy as np


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
