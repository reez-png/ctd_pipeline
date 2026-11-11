import numpy as np
from suncat.calc.salinity import sal78, practical_salinity, C_STD
from suncat.calc.thermo import potential_temperature
from suncat.calc.depth import pressure_to_depth


def test_pss78_standard_seawater():
    # by definition of PSS-78, S = 35 at R = 1, T68 = 15, P = 0
    assert abs(sal78(1.0, 15.0, 0.0) - 35.0) < 1e-3


def test_practical_salinity_standard_conductivity():
    s = practical_salinity(C_STD, 15.0 / 1.00024, 0.0)   # standard seawater at ITS-68 15 C
    assert abs(float(s) - 35.0) < 1e-2


def test_potential_temperature_surface_identity():
    T = 12.3
    assert abs(float(potential_temperature(35.0, T, 0.0)) - T) < 1e-9


def test_pressure_to_depth_zero_and_monotonic():
    assert pressure_to_depth(0.0, 0.0) == 0.0
    assert pressure_to_depth(1000.0, 30.0) > pressure_to_depth(500.0, 30.0) > 0
