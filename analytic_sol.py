# This file contains analytic solutions to compare numerical
# solutions against.

import numpy as np
from scipy.integrate import quad

# add to ref list
# this is using an analytic solution from Inan and Bahadir
# since this is an infinite sum, use n as a parameter
def IB_problem1(nu, x, t, n_max):
    
    # integrate over the interval (a,b)=(0,1)
    # not (x[0],x[-1]) since u(x,t) can be solved at a single 
    # point rather than over an interval
    
    a = 0; b = 1
    
    def a_n_func(x,n):
        
        f = np.exp(-(1 - np.cos(np.pi * x)) / (2 * np.pi * nu)) 
        return f * np.cos(n * np.pi * x)
   
    def a_n(n):
        # maybe come back and create finer grid to integrate over
        # rather than using the grid input to the outer function
        a_n = quad(a_n_func, a, b, args=(n))
        
        if n != 0:
            a_n = 2 * a_n
        
        return a_n[0]
     
    def sum_over_terms(x_c, t_c):
        
        num_sum = 0
        den_sum = 0

        for n in range(1, n_max):
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

