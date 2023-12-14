import numpy as np
from reactive_transport import Analytical, Crank_Nicolson, Forward_Difference


def test_boundary_values():
    # Test boundary values for Crank_Nicolson function
    cn_result = Crank_Nicolson(
        L=100.0, DX=1.0, Tmax=10.0, por=0.4, DT=0.1, Cf=1.0, Cb=0.5, q=0.8, DH=0.2
    )
    assert cn_result[0][-1] == 1.0
    assert cn_result[-1][-1] == 0.5

    # Test boundary values for Analytical function
    ana_result = Analytical(
        L=100.0, DX=1.0, Tmax=10.0, por=0.4, DT=0.1, Cf=1.0, Cb=0.5, q=0.8, DH=0.2
    )
    assert ana_result[0] == 1.0
    assert ana_result[-1] == 0.5

    # Test boundary values for Forward_Difference function
    fd_result = Forward_Difference(
        L=100.0, DX=1.0, Tmax=10.0, por=0.4, DT=0.1, Cf=1.0, Cb=0.5, q=0.8, DH=0.2
    )
    assert fd_result[0][-1] == 1.0
    assert fd_result[-1][-1] == 0.5


def test_values_within_range():
    # Test if calculated values are within the range of Cf and Cb for Crank_Nicolson
    cn_result = Crank_Nicolson(
        L=100.0, DX=1.0, Tmax=10.0, por=0.4, DT=0.1, Cf=1.0, Cb=0.5, q=0.8, DH=0.2
    )[:, -1]
    assert np.all(np.logical_and(cn_result >= 0.5, cn_result <= 1.0))

    # Test if calculated values are within the range of Cf and Cb for Analytical
    ana_result = Analytical(
        L=100.0, DX=1.0, Tmax=10.0, por=0.4, DT=0.1, Cf=1.0, Cb=0.5, q=0.8, DH=0.2
    )
    assert np.all(np.logical_and(ana_result >= 0.5, ana_result <= 1.0))

    # Test if calculated values are within the range of Cf and Cb for Forward_Difference
    fd_result = Forward_Difference(
        L=100.0, DX=1.0, Tmax=10.0, por=0.4, DT=0.1, Cf=1.0, Cb=0.5, q=0.8, DH=0.2
    )[:, -1]
    assert np.all(np.logical_and(fd_result >= 0.5, fd_result <= 1.0))


def test_extreme_cases():
    # Test extreme values for Crank_Nicolson function when Cf = Cb = c
    constant_c = 0.7  # can change this value to any non-negative constant

    cn_result = Crank_Nicolson(
        L=100.0,
        DX=1.0,
        Tmax=10.0,
        por=0.4,
        DT=0.1,
        Cf=constant_c,
        Cb=constant_c,
        q=0.8,
        DH=0.2,
    )[:, -1]
    assert np.allclose(cn_result, constant_c)

    # Test extreme values for Analytical function when Cf = Cb = c
    ana_result = Analytical(
        L=100.0,
        DX=1.0,
        Tmax=10.0,
        por=0.4,
        DT=0.1,
        Cf=constant_c,
        Cb=constant_c,
        q=0.8,
        DH=0.2,
    )
    assert np.allclose(ana_result, constant_c)

    # Test extreme values for Forward_Difference function when Cf = Cb = c
    fd_result = Forward_Difference(
        L=100.0,
        DX=1.0,
        Tmax=10.0,
        por=0.4,
        DT=0.1,
        Cf=constant_c,
        Cb=constant_c,
        q=0.8,
        DH=0.2,
    )[:, -1]
    assert np.allclose(fd_result, constant_c)


def test_close_results():
    # Test that all three functions return same size and approx. value for same inputs
    cn_result = Crank_Nicolson(
        L=100.0, DX=1.0, Tmax=10.0, por=0.4, DT=0.1, Cf=1.0, Cb=0.5, q=0.8, DH=0.2
    )[:, -1]
    ana_result = Analytical(
        L=100.0, DX=1.0, Tmax=10.0, por=0.4, DT=0.1, Cf=1.0, Cb=0.5, q=0.8, DH=0.2
    )
    fd_result = Forward_Difference(
        L=100.0, DX=1.0, Tmax=10.0, por=0.4, DT=0.1, Cf=1.0, Cb=0.5, q=0.8, DH=0.2
    )[:, -1]

    assert len(cn_result) == len(ana_result) == len(fd_result)
    assert np.allclose(cn_result, ana_result)
    assert np.allclose(ana_result, fd_result)
    assert np.allclose(cn_result, fd_result)
