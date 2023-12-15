import numpy as np
from reactive_transport import Analytical, Crank_Nicolson, Forward_Difference


def test_boundary_values():
    # Test boundary values for the Crank_Nicolson function
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
        L=100.0, DX=1.0, Tmax=5.0, por=0.3, DT=0.05, Cf=100.0, Cb=0.0, q=1.5, DH=2.5
    )[:, -1]
    assert np.all(np.logical_and(cn_result >= 0.0, cn_result <= 100.0))

    # Test if calculated values are within the range of Cf and Cb for Analytical
    ana_result = Analytical(
        L=100.0, DX=1.0, Tmax=5.0, por=0.3, DT=0.05, Cf=100.0, Cb=0.0, q=1.5, DH=2.5
    )
    assert np.all(np.logical_and(ana_result >= 0.0, ana_result <= 100.0))

    # Test if calculated values are within the range of Cf and Cb for Forward_Difference
    fd_result = Forward_Difference(
        L=100.0, DX=1.0, Tmax=5.0, por=0.3, DT=0.05, Cf=100.0, Cb=0.0, q=1.5, DH=2.5
    )[:, -1]
    assert np.all(np.logical_and(fd_result >= 0.0, fd_result <= 100.0))


def test_extreme_cases():
    # Test extreme values for Crank_Nicolson function when Cf = Cb = 0
    cn_result = Crank_Nicolson(
        L=100.0,
        DX=1.0,
        Tmax=10.0,
        por=0.3,
        DT=0.05,
        Cf=0.0,
        Cb=0.0,
        q=1.5,
        DH=2.5,
    )[:, -1]
    assert np.allclose(cn_result, 0.0)

    # Test extreme values for Analytical function when Cf = Cb = 0
    ana_result = Analytical(
        L=100.0,
        DX=1.0,
        Tmax=10.0,
        por=0.3,
        DT=0.05,
        Cf=0.0,
        Cb=0.0,
        q=1.5,
        DH=2.5
    )
    assert np.allclose(ana_result, 0.0)

    # Test extreme values for Forward_Difference function when Cf = Cb = 0
    fd_result = Forward_Difference(
        L=100.0,
        DX=1.0,
        Tmax=10.0,
        por=0.3,
        DT=0.05,
        Cf=0.0,
        Cb=0.0,
        q=1.5,
        DH=2.5
    )[:, -1]
    assert np.allclose(fd_result, 0.0)


def test_close_results():
    # Test that all three functions return the same size and approx. value for same inputs
    cn_result = Crank_Nicolson(
        L=100.0, DX=1.0, Tmax=5.0, por=0.3, DT=0.05, Cf=100.0, Cb=0.0, q=1.5, DH=2.5
    )[:, -1]
    ana_result = Analytical(
        L=100.0, DX=1.0, Tmax=5.0, por=0.3, DT=0.05, Cf=100.0, Cb=0.0, q=1.5, DH=2.5
    )
    fd_result = Forward_Difference(
        L=100.0, DX=1.0, Tmax=5.0, por=0.3, DT=0.05, Cf=100.0, Cb=0.0, q=1.5, DH=2.5
    )[:, -1]

    assert len(cn_result) == len(ana_result) == len(fd_result)
    assert np.allclose(cn_result, ana_result, atol = 0.01)
    assert np.allclose(ana_result, fd_result, atol = 0.01)
    assert np.allclose(cn_result, fd_result, atol = 0.01)
