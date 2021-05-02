import numpy as np
import matplotlib.pyplot as plt
# from scipy.integrate import simpson, odeint

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


# initialize tridiagonal matrices that form the 
# system of equations
# this function requires that M > 2
def init_matrices(alpha, M):
    
    A = np.zeros((M-1,M-1))
    B = np.zeros((M-1,M-1))
    
    alpha_2 = alpha * .5
    alpha_p1 = alpha + 1
    alpha_m1 = 1 - alpha
    
    # set first and last rows
    A[0][:2] = np.asarray([alpha_p1, -alpha_2])
    A[M-2][-2:] = np.asarray([-alpha_2, alpha_p1])
    
    B[0][:2] = np.asarray([alpha_m1, alpha_2])
    B[M-2][-2:] = np.asarray([alpha_2, alpha_m1])

    for i in range(1,M-2):
        
        A[i][i-1:i+2] = np.asarray([-alpha_2, alpha_p1, -alpha_2])
        B[i][i-1:i+2] = np.asarray([alpha_2, alpha_m1, alpha_2])
                
    return A, B

# might not need to actually compute A due to the high degree 
# of symmetry
# this is following the tridiagonalization procedure for Landau
def solve_tri_diag(alpha, b):
    
    d = alpha + 1
    a = -alpha/2
    c = -alpha/2 

    h = np.zeros_like(b)
    p = np.zeros_like(b)

    len_b = len(b)

    # set first element of p and h
    p[0] = b[0] / d
    h[0] = c / d

    for i in range(1,len_b):
        h[i] = c / (d - a * h[i-1])
        p[i] = (b[i] - a * p[i-1]) / (d - a * h[i-1])
    
    # initialize solution vector
    x = np.zeros_like(b)
    x[-1] = p[-1]
   
    # backwards solve from x[-1] -> x[0]
    for i in range(0,len_b - 1):
        x[-(i+2)] = p[-(i+2)] - h[-(i+2)] * x[-(i+1)]
    
    return x

# solves the heat equation for phi. This will still need to be transformed 
# after solved to get final solution u(x,t)
def solve_PDE(phi_0, x, t, nu):
    
    delta_t = t[1] - t[0]
    delta_x = x[1] - x[0]

    alpha = delta_t / (delta_x)**2 * nu 
    
    # initialize phi(x,t) and set phi(x,0)
    phi_sol = np.zeros((len(x), len(t)))
    

    return 
