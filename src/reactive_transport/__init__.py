import math

import numpy as np
from numpy.typing import NDArray


def simple_transport(grad, D, v, grad_c):
    """
    Simple 1D change in concentration over time for one time step at one location.
     Inputs:
         grad (float): Slope of confined groundwater table
         D (float): Diffusion coefficient
         v (float): Velocity of groundwater
         grad_c (float): Slope of concentration gradient
    """
    dc_dt = D * grad**2 - v * grad_c
    return dc_dt


def Crank_Nicolson(
    L: float,
    DX: float,
    Tmax: float,
    por: float,
    DT: float,
    Cf: float,
    Cb: float,
    q: float,
    DH: float,
) -> NDArray[np.float64]:
    """
    Computes the numerical solution using the Crank-Nicolson method

    Inputs:
    - L (float): Length of the domain
    - DX (float): Spatial step size
    - Tmax (float): Maximum simulation time
    - por (float): Porosity
    - DT (float): Time step size
    - Cf (float): Concentration at the front boundary
    - Cb (float): Concentration at the back boundary
    - q (float): Flow rate
    - DH (float): Dispersion coefficient

    Returns:
    - NDArray[np.float64]: Concentration profile across the domain and across time
    """
    N = int(L / DX) + 1  # needs to be changed to match Forward Difference formatting

    number_timesteps = int(Tmax / DT)
    Cnew = np.zeros([N, number_timesteps], dtype=np.float64)

    # set coefficients in in FD equations
    p = (DH * DT) / (DX**2)
    r = (q * DT) / (2 * DX)

    a = p - r
    b = 2 * por + 2 * p
    c = p + r

    for i in range(number_timesteps):
        Cnew[0, i] = Cf  # Boundary condition on the left
        Cnew[50, i] = Cb  # Boundary condition on the right

        D = np.zeros(N + 1, dtype=np.float64)
        E = np.zeros(N + 1, dtype=np.float64)
        F = np.zeros(N + 1, dtype=np.float64)

        for k in range(1, N - 1):
            if k == 1:
                D[k] = (
                    a * Cnew[k + 1, i - 1]
                    + (4 * por - b) * Cnew[k, i - 1]
                    + c * Cnew[k - 1, i - 1]
                )
                E[k] = a / b
                F[k] = (D[k] + c * Cnew[k - 1, i - 1]) / b
            else:
                D[k] = (
                    a * Cnew[k + 1, i - 1]
                    + (4 * por - b) * Cnew[k, i - 1]
                    + c * Cnew[k - 1, i - 1]
                )
                E[k] = a / (b - c * E[k - 1])
                F[k] = (D[k] + c * F[k - 1]) / (b - c * E[k - 1])

        for k in range(N - 2, 0, -1):
            Cnew[k, i] = F[k] + E[k] * Cnew[k + 1, i]
    return Cnew


def Analytical(
    L: float,
    DX: float,
    Tmax: float,
    por: float,
    DT: float,
    Cf: float,
    Cb: float,
    q: float,
    DH: float,
) -> NDArray[np.float64]:
    """
    Computes the analytical solution for reactive transport.

    Inputs:
    - L (float): Length of the domain
    - DX (float): Spatial step size
    - Tmax (float): Maximum simulation time
    - por (float): Porosity
    - DT (float): Time step size
    - Cf (float): Concentration at the front boundary
    - Cb (float): Concentration at the back boundary
    - q (float): Flow rate
    - DH (float): Dispersion coefficient

    Returns:
    - NDArray[np.float64]: Concentration profile across the domain at one time
    """
    N = int(L / DX) + 1

    v = q / por
    D = DH / por

    Cnew_analytic = np.zeros(N, dtype=np.float64)

    Cnew_analytic[0] = Cf
    Cnew_analytic[50] = Cb

    for k in np.arange(1, N - 2, DX):
        x = k * DX
        Cnew_analytic[k] = (Cf / 2) * (
            (math.erfc((x - v * DT) / (2 * math.sqrt(D * DT))))
            + (
                math.exp((v * x) / D)
                * math.erfc((x + v * DT) / (2 * math.sqrt(D * DT)))
            )
        )
    return Cnew_analytic


def Forward_Difference(
    L: float,
    DX: float,
    Tmax: float,
    por: float,
    DT: float,
    Cf: float,
    Cb: float,
    q: float,
    DH: float,
) -> NDArray[np.float64]:
    """
    Computes the numerical solution using the Forward Difference method

    Inputs:
    - L (float): Length of the domain
    - DX (float): Spatial step size
    - Tmax (float): Maximum simulation time
    - por (float): Porosity
    - DT (float): Time step size
    - Cf (float): Concentration at the front boundary
    - Cb (float): Concentration at the back boundary
    - q (float): Flow rate
    - DH (float): Dispersion coefficient

    Returns:
    - NDArray[np.float64]: Concentration profile across the domain and across time
    """
    number_nodes = int(L / DX)
    number_timesteps = int(Tmax / DT)
    Cnew = np.zeros([number_nodes, number_timesteps], dtype=np.float64)

    # set coefficients in in FD equations
    D = DH / por
    p = (D * DT) / (DX**2)
    r = (q * DT) / (2 * por * DX)

    a = p - r
    b = 1.0 - 2 * p
    c = p + r

    for i in range(number_timesteps):
        Cnew[0, i] = Cf
        Cnew[number_nodes - 1, i] = Cb

        for k in np.arange(1, number_nodes - 2, DX):
            Cnew[k, i] = (
                a * Cnew[k + 1, i - 1] + b * Cnew[k, i - 1] + c * Cnew[k - 1, i - 1]
            )
    return Cnew
