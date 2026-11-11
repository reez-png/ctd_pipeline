import numpy as np
from suncat.calc.binning import bin_mean


def test_bin_mean_centres_and_means():
    x = np.array([0.4, 0.6, 2.4, 2.6])     # step 2 -> bins centred at 0 and 2
    y = np.array([10.0, 20.0, 30.0, 40.0])
    cen, m = bin_mean(x, y, 2.0)
    assert list(cen) == [0.0, 2.0]
    assert list(m) == [15.0, 35.0]


def test_bin_mean_drops_nan_and_empty():
    cen, m = bin_mean(np.array([np.nan, 1.0]), np.array([5.0, np.nan]), 2.0)
    assert cen.size == 0 and m.size == 0
    cen2, m2 = bin_mean(np.array([]), np.array([]), 2.0)
    assert cen2.size == 0 and m2.size == 0
