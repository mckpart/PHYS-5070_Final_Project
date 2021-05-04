import numpy as np
import matplotlib.pyplot as plt

# currently this assumes time independent boundary conditions

# initialize tridiagonal matrices that form the 
# system of equations
# this function requires that M > 2
def init_matrices(alpha, M):
    '''
    Initializes both tridiagonal matrices that represent the 
    system of difference equations.

    Args:
        alpha (float): alpha = nu * delta t / (delta x)^2
        M (int): each matrix will have dimension MxM

    Returns:
        A (array): 2-d array representing the matrix that acts
            on the phi_(j+1) array
        B (array): 2-d array representing the matrix that acts
            on the phi_(j) array
    '''  

    A = np.zeros((M, M))
    B = np.zeros((M, M))
    
    alpha_2 = alpha * .5
    alpha_p1 = alpha + 1
    alpha_m1 = 1 - alpha
    
    # set first and last rows
    A[0][:2] = np.asarray([alpha_p1, -alpha_2*2])
    A[M-1][-2:] = np.asarray([-alpha_2*2, alpha_p1])
    
    B[0][:2] = np.asarray([alpha_m1, alpha_2*2])
    B[M-1][-2:] = np.asarray([alpha_2*2, alpha_m1])

    for i in range(1,M-1):
        
        A[i][i-1:i+2] = np.asarray([-alpha_2, alpha_p1, -alpha_2])
        B[i][i-1:i+2] = np.asarray([alpha_2, alpha_m1, alpha_2])
                
    return A, B

# might not need to actually compute A due to the high degree 
# of symmetry
# this is following the tridiagonalization procedure for Landau
def solve_tri_diag(A, b):
   
    d = np.diagonal(A)
    a = np.diagonal(A,offset = -1)
    c = np.diagonal(A, offset = 1) 

    len_b = len(b)
    
    h = np.zeros(len_b-1)
    p = np.zeros_like(b)

    # set first element of p and h
    p[0] = b[0] / d[0]
    h[0] = c[0] / d[0]
    
    # compute h and p. The algorithm in the text says a[i] rather
    # than a[i-1], but this is because they begin with a[2]
    for i in range(1,len_b):
        if i < len_b - 1:
            h[i] = c[i] / (d[i] - a[i-1] * h[i-1])
        p[i] = (b[i] - a[i-1] * p[i-1]) / (d[i] - a[i-1] * h[i-1])
    
    # initialize solution vector
    x = np.zeros_like(b)
    x[-1] = p[-1]
   
    # backwards solve from x[-1] -> x[0]
    for i in range(0,len_b - 1):
        x[-(i+2)] = p[-(i+2)] - h[-(i+1)] * x[-(i+1)]
    
    return x


# solves the heat equation for phi. This will still need to be transformed
# after solved to get final solution u(x,t)
def solve_PDE(phi_0, x, t, nu):

    delta_t = t[1] - t[0]
    delta_x = x[1] - x[0]

    alpha = delta_t / (delta_x)**2 * nu
    
    # the matrices A,B will have dim MxM
    M = len(x)

    # initialize phi(x,t) and set phi(x,0)
    phi_sol = np.zeros((len(x), len(t)))

    # set initial conditions. the boundary conditions are assumed
    # to be phi_x(0,t)=phi_x(L,t) = 0
    phi_sol[:,0] = phi_0

    # compute the A and B matrices
    A, B = init_matrices(alpha, M)

    # initial b = B * phi_0. Then the phi at later times
    # is computed through A*w = b
    b = np.matmul(B, phi_0)
    for i in range(len(t)):
        w = solve_tri_diag(A, b)
        phi_sol[:,i] = w
        b = np.matmul(B,phi_sol[:,i])

    return phi_sol 

# transforms phi(x,t) to the final solution u(x,t)
def transform_phi(phi, del_x, nu):
    u = np.zeros_like(phi)
    
    # set the boundary conditions
    u[0,:] = 0
    u[-1,:] = 0
    
    lx, lr = phi.shape
    
    # do not update boundaries
    for i in range(1,lx-1):
        # each step in time must be updated, including t = 0
        for j in range(lr):
            u[i,j]= -nu * (phi[i+1,j]-phi[i-1,j])/(del_x * phi[i,j])
            
    return u
