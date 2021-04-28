# This file contains analytic solutions to compare numerical
# solutions against.

import numpy as np
from scipy.integrate import quad

# since this is an infinite sum, use n as a parameter
def IB_analytic_sum(nu, x, t, n_max, problem = 1):
    '''
    Computes the analytic solution for the either of the first 
    two problems outlined in Inan And Bahadir.

    Args:
        nu (float): kinematic viscosity -  must be nonnegative
            for it to have physical meaning.
        x (array): x-values to solve u(x,t) over
        t (array): t-values to solve u(x,t) over
        n_max (int): each point is a Fourier series, and n_max
            is the truncation value for the series.
        problem (int): determines which problem is to be solved
            - 1: IC: u(x,0) = sin(pi*x), BC: u(0,t) = u(1,t) = 0
            - 2: IC: u(x,0 = 4x(1-x),    BC: u(0,t) = u(1,t) = 0

    Returns:
        u (array): 2-d array containing the solution u(x,t). 
            The element u[i,j] can be interpretted as the analytic
            function evaluated as u(x[i],t[j]). 
    
    '''

    # integrate over the interval (a,b)=(0,1) 
    a = 0; b = 1
    
    def a_n_func(x,n):
        
        if problem == 1:
            f = np.exp(-(1 - np.cos(np.pi * x)) / (2 * np.pi * nu)) 
        
        elif problem == 2:
            f = np.exp(-x**2 * (3 - 2 * x) / (3 * nu))

        return f * np.cos(n * np.pi * x)
   
    def a_n(n):
        # quad returns the integrated function and the error,
        # so the first element of the returned tuple is taken
        a_n = quad(a_n_func, a, b, args=(n))[0]
        
        if n != 0:
            a_n = 2 * a_n
        
        return a_n
     
    def sum_over_terms(x_c, t_c):
        
        num_sum = 0
        den_sum = 0

        for n in range(1, n_max):
            # the numerator and denominator are each truncated
            # infinite series
            num_sum += (a_n(n) * np.exp(-(n * np.pi)**2 * nu * t_c)
                    * n * np.sin(n * np.pi * x_c))
            
            den_sum += (a_n(n) * np.exp(-(n * np.pi)**2 * nu * t_c)
                    * np.cos(n * np.pi * x_c))

        den_sum += a_n(0)
        num_sum *= 2 * np.pi * nu

        return num_sum, den_sum

    u = np.zeros((len(x), len(t)))
    for i in range(len(x)):
        for j in range(len(t)):
            num, den = sum_over_terms(x[i], t[j])
            u[i,j] = num / den

    return u

