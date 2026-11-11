import numpy as np
from suncat.section import descent_start


def _profile():
    # 60-scan surface soak (longer than the 24-scan window), clean descent, then upcast
    soak = np.full(60, 0.5)
    down = np.linspace(1.0, 60.0, 120)
    up = np.linspace(60.0, 2.0, 80)
    return np.concatenate([soak, down, up])


def test_descent_start_on_known_profile():
    start, imax = descent_start(_profile(), win=24, min_drop=0.10)
    assert imax == 179            # deepest scan (end of the descent), deterministic
    assert start == 36            # most of the 60-scan soak excluded, deterministic
    assert 0 < start < 60         # onset lands inside/after the soak, never at 0


def test_guards():
    assert descent_start(np.array([]), 24) == (0, 0)
    assert descent_start(np.full(5, np.nan), 24) == (0, 0)


def test_window_parameterisations_agree():
    # step00b passes the window in scans; step00c passes seconds * sample-rate. Equal windows -> equal result.
    pv = _profile()
    by_scans = descent_start(pv, win=24, min_drop=0.10)
    by_seconds = descent_start(pv, win=max(1, int(round(6.0 * 4.0))), min_drop=0.10)
    assert by_scans == by_seconds
