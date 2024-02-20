
import numpy as np
from scipy.sparse import diags
from scipy.sparse import csr_array
from scipy import special as sp
import math


def jacobi(C, stop,store_step = 10, object = None):
    """
    Jacobi method.
    Inputs:
        - C: Matrix with initial conditions
        - Stop: Stopping criteria
        - Boundary: Boundary conditions that wont change, default top and last row
        - Object: If there is any object, default to 0 and -1 rows
    """
    w = 1/4
    n_cols = C.shape[1]

    diagonals = [np.ones(n_cols-1), np.ones(n_cols-1)]
    
    M1 = diags(diagonals , [-1, 1])
    M2 = diags(diagonals, [-1,1]).toarray()

    ## Cyclic boundary conditions
    M2[-2,0] = 1
    M2[1,-1] = 1

    n_count = 0

    while True:
        n_count += 1
        C_b = np.copy(C)
        c1=  M1@C
        c2 = C@M2
        C[1:-1] = (1/4 *(c1+c2))[1:-1]

        if np.allclose(C, C_b, atol=stop):
            yield C
            break
        if n_count%store_step == 0:
            yield C
    


def gauss_seidel(C, stop,store_step = 10, object = None):
    """
    Performs Gauss Seidel iterations
    Inputs:
        - A: Matrix with initial conditions/positions
        - Stop: Stop criteria
    """
    w = 1/4
    n_count = 0
    while True:
        n_count += 1 
        C_b = np.copy(C)
        for i in range(1,C.shape[0]-1):
            C[i,0] = w*(C[i+1,0] + C[i-1,0] + C[i,1] + C[i,-2])
            for j in range(1,C.shape[1]-1):
                C[i,j] = w*(C[i+1,j] + C[i-1,j] + C[i,j+1] + C[i,j-1])
            C[i,-1] = w*(C[i+1,-1] + C[i-1,-1] + C[i,1] + C[i,-2])

        if np.allclose(C, C_b, atol=stop):
            yield C
            break
        if n_count%store_step == 0:
            yield C


def sor(C,w,stop = 0.001,store_step = 10, object = None):
    """
    Performs Successive over relaxation.
    Inputs:
        - A: Matrix A with all values
        - w: weight
        - stop: simulation stopper
    """
    n_count = 0
    while True:
        n_count += 1
        C_b = np.copy(C)
        for i in range(1,C.shape[0]-1):
            C[i,0] = w*(C[i+1,0] + C[i-1,0] + C[i,1] + C[i,-2]) + (1-w)*C[i,0]
            for j in range(1,C.shape[1]-1):
                C[i,j] = (w/4)*(C[i+1,j] + C[i-1,j] + C[i,j+1] + C[i,j-1]) + (1-w)*C[i,j]
            C[i,-1] = w*(C[i+1,-1] + C[i-1,-1] + C[i,1] + C[i,-2]) + (1-w)*C[i,-1]
        
        if np.allclose(C, C_b, atol=stop):
            yield C
            break
        if n_count%store_step == 0:
            yield C
    



            

