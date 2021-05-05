# PHYS 5070 Final Project

The project contains three numerical methods to solve Burger's equation, a simple nonlinear partial differential equation. Two are explicit methods, each using different orders of approximation. The third is an implicit method that uses Crank Nicholson time stepping method. The goal is to understand different ways to solve this PDE, and then use the numerical methods to better understand the PDE itself. 

Tests of these methods are included in the [project write-up](./Project_Write-Up.ipynb) file. There the user can find examples of how to obtain a solution for a given problem. 

The Crank-Nicholson solver has the limitation that the boundary conditions are assumed to be 0. The other two solvers only require time independent boundary conditions. 

The file [analytic_sol.py](./analytic_sol.py) contains two analytic solutions as described by Inan and Bahadir. The solutions are both truncated infinite series that require numerical integration to compute the coefficients in the sum. These analytic solutions are then later added to test the accuracy of the numerical schemes. 

In [finite_difference.py](./finite_difference.py) there are two finite difference solvers. The first solver FiniteElements_Burgers takes uses the most basic finite difference scheme to numerically approximate the solution. The second solver FiniteElements_Burgers_Higher_Order uses a higher order finite difference method. Finally, this file also contains the function Plot_func which plots the solutions in three dimensions.

In [crank_nicholson.py](./crank_nicholson.py), there is a wrapper function that is called and then the solution is retured. Within this function, the transformation between the solution space and the space that the problem is solved is performed, given an initial condition in the solution space. The PDE solver solves the diffusion equation since using the Cole-Hopf transformation the Burger's equation is reduced to a diffusion equation. Within this solver, a system of difference equations is solved. This system of equations can be represented using two tridiagonal matrices, so a solver for that is included as well. All of these functions are tested in the main write-up document. 
