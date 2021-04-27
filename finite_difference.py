import numpy as np
import matplotlib.pyplot as plt

# this is an explicit finite difference method
# add in a maximum number of iterations
def finite(U0, x, t, nu):

    dx = x[1] - x[0]
    dt = t[1] - t[0]

    u = np.zeros((len(x), len(t)))
    u[:,0] = U0

    for j in range(0,len(t)-1):
        for i in range(1,len(x)-1):

            # calculate the diffusion and advection terms
            diff = u[i+1,j] - 2 * u[i,j] + u[i-1,j]
            adv = -(u[i+1,j]**2 - u[i-1,j]**2) / 4

            u[i,j+1] = u[i,j] + dt * (nu * diff/(dx**2) + adv/dx)

    return u
