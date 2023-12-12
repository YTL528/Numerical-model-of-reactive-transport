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


def Crank_Nicolson(L, DX, Tmax, por, DT, Cf, Cb, q, DH):
    """
    Documentation needed
     Inputs:

    """
    N = int(L / DX) + 1  # needs to be changed to match Forward Difference formatting
    # node_locations = np.arange(0, L + DX, DX)
    number_timesteps = int(Tmax / DT)
    # time_increments = np.arange(DT, Tmax + DT, DT)
    Cnew = np.zeros([N, number_timesteps])

    # set coefficients in in FD equations
    D = DH / por
    p = (D * DT) / (DX**2)
    r = (q * DT) / (2 * DX)  # maybe * por?

    a = p - r
    b = 2 * por + 2 * p
    c = p + r

    for i in range(number_timesteps):
        Cnew[0, i] = Cf  # Boundary condition on the left
        Cnew[50, i] = Cb  # Boundary condition on the right

        D = np.zeros(N + 1)
        E = np.zeros(N + 1)
        F = np.zeros(N + 1)

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


def Analytical(L, DX, Tmax, por, DT, Cf, Cb, q, DH):
    N = int(L / DX) + 1

    v = q / por
    D = DH / por

    # node_locations = np.arange(0, L + DX, DX)

    Cnew_analytic = np.zeros(N)

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


def Forward_Difference(L, DX, Tmax, por, DT, Cf, Cb, q, DH):
    number_nodes = int(L / DX)
    # node_locations = np.arange(0, L + DX, DX)
    number_timesteps = int(Tmax / DT)
    # time_increments = np.arange(DT, Tmax + DT, DT)
    Cnew = np.zeros([number_nodes, number_timesteps])

    # set coefficients in in FD equations
    D = DH / por
    p = (D * DT) / (DX**2)
    r = (q * DT) / (2 * por * DX)  # not sure if por should be here

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
