# This file contains the finite elements solvers for the Burger's equation and the functions to plot the solutions

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def FiniteElements_Burgers(nu,grid,time,initial_conditions):
    """
    Returns the solution to the Burger's equation given specific initial conditions
    
    Arguments:
    - nu: as a decimal
    - grid: the grid points as an array (must be evenly spaced)
    - time: an array of the time points desired (must also be evenly spaced)
    - initial_conditions: array, must be the same length as grid
    
    Returns:
    - u: as u[t][x]
    
    """
    
    delta_x = grid[1] - grid[0]
    delta_t = time[1] - time[0]
    
    u = [initial_conditions]
    
    dudt = np.empty(len(grid))
    for t in time[1:]:
        dudt[0] = 0
        dudt[-1] = 0
        for n in range(1,len(grid)-1):
            dudt[n] = nu*(u[-1][n+1] - 2*u[-1][n] + u[-1][n-1])/delta_x**2 - u[-1][n]*(u[-1][n+1] - u[-1][n-1])/(2*delta_x)
            
        u.append(u[-1] + dudt*delta_t)
        
    return np.array(u)


def FiniteElements_Burgers_Higher_Order(nu,grid,time,initial_conditions):
    """
    Returns the solution to the Burger's equation given specific initial conditions
    
    Arguments:
    - nu: as a decimal
    - grid: the grid points as an array (must be evenly spaced)
    - time: an array of the time points desired (must also be evenly spaced)
    - initial_conditions: array, must be the same length as grid
    
    Returns:
    - u: as u[t][x]
    
    """
    
    delta_x = grid[1] - grid[0]
    delta_t = time[1] - time[0]
    
    u = [initial_conditions]
    
    dudt = np.empty(len(grid))
    for t in time[1:]:
        dudt[0] = 0
        dudt[-1] = 0
        
        dudt[1] = nu*(15/4*u[-1][1] - 77/6*u[-1][2] + 107/6*u[-1][3] - 13*u[-1][4] + 61/12*u[-1][5] - 5/6*u[-1][6])/delta_x**2 - u[-1][1]*(-25/12*u[-1][1] + 4*u[-1][2] - 3*u[-1][3] + 4/3*u[-1][4] - 1/4*u[-1][5])/delta_x
        
        dudt[-2] = nu*(15/4*u[-1][-2] - 77/6*u[-1][-3] + 107/6*u[-1][-4] - 13*u[-1][-5] + 61/12*u[-1][-6] - 5/6*u[-1][-7])/delta_x**2 - u[-1][-2]*(25/12*u[-1][-2] - 4*u[-1][-3] + 3*u[-1][-4] - 4/3*u[-1][-5] + 1/4*u[-1][-6])/delta_x
        
        for n in range(2,len(grid)-2):
            dudt[n] = nu*(-1/12*u[-1][n+2] + 4/3*u[-1][n+1] - 5/2*u[-1][n] + 4/3*u[-1][n-1] - 1/12*u[-1][n-2])/delta_x**2 - u[-1][n]*(-1/12*n[-1][n+2] + 2/3*u[-1][n+1] - 2/3*u[-1][n-1] + 1/12*u[-1][n-2])/(2*delta_x)
            
        u.append(u[-1] + dudt*delta_t)
        
    return np.array(u)


def Plot_func(grid,time,u):
    """
    Plots u(x,t) as a 3D plot.
    
    Arguments:
    - grid: a 1D array of the grid points
    - time: a 1D array of the time points
    - u: a 2D array as u[t][x]
    
    Returns:
    - a 3D plot
    
    """
    Grid,Time = np.meshgrid(grid,time)
    
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(Grid,Time,u)
    