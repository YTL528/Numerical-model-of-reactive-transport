def simple_transport(grad, D, v, grad_c):
    """
    Simple 1D change in concentration over time for one time step at one location.
     Inputs:
         grad (float): Slope of confined groundwater table
         D (float): Diffusion coefficient
         v (float): Velocity of groundwater
         grad_c (float): Slope of concentration gradient
    """
    dc_dt = D * grad ** 2 - v * grad_c
    return dc_dt
