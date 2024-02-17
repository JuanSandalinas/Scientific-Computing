
import numpy as np
from scipy.sparse import diags
from scipy.sparse import csr_array
from scipy import special as sp
import math



def jacobi(C, stop):
    """
    Jacobi method:
    Inputs:
        - C: Matrix with initial conditions
        - Stop: Stopping criteria
    """
    w = 1/4
    n_rows = A.shape[0]
    n_cols = A.shape[1]
    diagonals = [np.ones(ncol-1),np.ones(ncol), np.ones(ncols-1)]
    M1 = diags(diagonals , [-1, 0, 1])
    M2 = diags([diagonals[0], diagonals[2]], [-1,1]).toarray()
    M2[-2,0] = 1
    M2[1,-1] = 1
    while True:
        C_d = C
        C_ = 1/4 * M1@C@M2
        C[1:-1] = C_[1:-1]
        if np.abs(C-C_b).all() <= stop:
            break
    return C


def gauss_seidel(C, stop):
    """
    Performs Gauss Seidel iterations
    Inputs:
        - A: Matrix with initial conditions/positions
        - Stop: Stop criteria
    """
    w = 1/4
    C_b = np.ones((C.shape[0], C.shape[1]))
    while True:
        for i in range(0,A.shape[0]):
            C[i,-1] = c*(C[i+1,j] + C[i-1,j] + C[i,j+1] + C[i,-1])
            for j in range(1,A.shape[1]-1):
                C[i,j] = w*(C[i+1,j] + C[i-1,j] + C[i,j+1] + C[i,j-1])

            C[i,-1] = C[i,0]
        if np.abs(C-C_b).all() <= stop:
            break
        C_b = C

    return C

def sor(C,w,stop):
    """
    Performs Successive over relaxation.
    Inputs:
        - A: Matrix A with all values
        - w: weight
        - stop: simulation stopper
    """
    
    C_b = np.ones((C.shape[0], C.shape[1]))
    while True:
        for i in range(1,C.shape[0]-1):
            for j in C.size[1]:
                C_k[i,j] = w*(C_k[i+1,j] + C_k[i-1,j] + C_k[i,j+1] + C_k[i,j-1]) + (1-w)*C_k[i,j]
            C[i,-1] = C[i,0]

        if np.abs(C-C_b).all() <= stop:
            break
        A_b = A
    return A


def c_analytical(y, t, D=1, imax=100):
    """
    Returns analytical solution of diffusion equation
    Inputs:
        - y: y axis solution
        - t: time steps
    """

    lst = []
     
    result = 0
    for i in range(100):
        arg1 = math.erfc((1 - y_i + (2*i)) /(2 * np.sqrt(D)* t))
        arg2 = math.erfc((1 + y_i + (2*i)) /(2 * np.sqrt(D)* t))
        result += (arg1-arg2)
        lst.append(result)
    return np.array(lst)



            

