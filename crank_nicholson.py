import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson, odeint

# currently this assumes time independent boundary conditions

# create a function that transforms initial conditions using Cole-Hopf
# transformation

# look into making this function more robust to transform between phi and 
# u, not just for the initial conditions

# FOR NOW THE TRANSFORMED INITIAL CONDITIONS NEED TO BE COMPUTED ANALYTICALLY

# def transform_init(u_init, x_arr, nu):
#     '''
#     Transforms the initial conditions from u(x,0) to phi(x,0)
#     using the Cole-Hopf transformation.
# 
#     Args:
#         nu (float): kinematic viscosity -  must be nonnegative
#             for it to have physical meaning.
#         u_init (array): 
# 
#     Returns:
#         phi_init (array): 1-d array containing the the transformed
#             u(x,0) which is given by phi(x,0)
#     '''
# 
#     def ode(phi, x):
#         i = (x_arr == x)
#         u = u_init[i][0]
#         return -2 * nu * phi * u
# 
#     # come up with better way to set IC
#     phi_init = odeint(ode, y0=1, x)
# #     U = simpson(u_init, dx = .1)
# #     phi_init = np.exp(-simpson( / (2 * nu)) 
# 
#     return phi_init

 
