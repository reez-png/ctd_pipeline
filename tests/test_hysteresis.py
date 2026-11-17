"""suncat.nitrate.hysteresis on plain arrays (the form step02b uses): a known lag and the cruise summary."""
from __future__ import annotations

import numpy as np

from suncat.nitrate import hysteresis as H


def _cast(tau, seed, dt=1.1, zmax=70.0, v=0.57, noise=0.5):
    """A cast at v m/s down and up with a 20 s bottom stop; a tanh nitracline seen through a first-order lag
    integrated 20 times finer than the frames."""
    nd = int(zmax / v / dt)
    z = np.concatenate([np.linspace(1, zmax, nd), np.full(int(20 / dt), zmax), np.linspace(zmax, 1, nd)])
    t = np.arange(len(z)) * dt
    tf = np.linspace(0, t[-1], 20 * len(t))
    x = 5 + 8 * (1 + np.tanh((np.interp(tf, t, z) - 30) / 6))
    y = x.copy()
    if tau > 0:
        for i in range(1, len(tf)):
            y[i] = y[i - 1] + (x[i] - y[i - 1]) * (1 - np.exp(-(tf[i] - tf[i - 1]) / tau))
    y = np.interp(t, tf, y) + np.random.default_rng(seed).normal(0, noise, len(t))
    down = np.arange(len(t)) < nd
    return t, z, y, down, ~down


def test_cast_tau_and_summary_recover_a_known_lag():
    taus = np.arange(0, 61, 1.0)
    per = {f"c{k}": H.cast_tau(*_cast(29.0, k), taus=taus, dz=1.0) for k in range(5)}
    assert all(abs(r["tau"] - 29.0) <= 2.0 for r in per.values()), [r["tau"] for r in per.values()]
    s = H.summarise(per, taus)
    assert abs(s["tau_median"] - 29.0) <= 1.5 and s["tau_ci90"][0] <= 29.5 and s["tau_ci90"][1] >= 28.5
    assert abs(s["tau_pooled"] - 29.0) <= 1.5 and s["rss_share_explained"] > 0.8


def test_no_lag_gives_zero_and_too_few_frames_is_skipped():
    taus = np.arange(0, 31, 1.0)
    assert H.cast_tau(*_cast(0.0, 1), taus=taus, dz=1.0)["tau"] <= 0.5
    t, z, y, dn, up = _cast(5.0, 2)
    y[up] = np.nan
    assert "skip" in H.cast_tau(t, z, y, dn, up, taus=taus, dz=1.0)
    assert H.summarise({}, taus) == {}
