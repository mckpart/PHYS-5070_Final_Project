# PHYS 5070 Final Project

The project contains three numerical methods to solve Burger's equation, a simple nonlinear partial differential equation. Two are explicit methods, each using different orders of approximation. The third is an implicit method that uses Crank Nicholson time stepping method. The goal is to understand different ways to solve this PDE, and then use the numerical methods to better understand the PDE itself. 

Tests of these methods are included in the Project_Write-Up.ipynb file. There the user can find examples of how to obtain a solution for a given problem. 

The Crank-Nicholson solver has the limitation that the boundary conditions are assumed to be 0. In other words, that $u(0,t) = u(L,t) = 0$. 