
import numpy as np
from scipy.sparse import diags
from scipy.sparse import csr_array
from scipy import special as sp
import math


def jacobi(C,object_, stop,store_step):
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
    non_cte = np.where(object_ == 0)
    while True:
        n_count += 1
        C_b = np.copy(C)
        c1=  M1@C
        c2 = C@M2
        C[non_cte[0],non_cte[1]] = (1/4 *(c1+c2))[non_cte[0],non_cte[1]]

        if np.allclose(C, C_b, atol=stop):
            yield C
            break
        if n_count%store_step == 0:
            yield C
    


def gauss_seidel(C, object_,stop,store_step):
    """
    Performs Gauss Seidel iterations
    Inputs:
        - A: Matrix with initial conditions/positions
        - Stop: Stop criteria
        - object: matrix that specifies which values need to be updated
    """
    w = 1/4
    n_count = 0
    
    non_cte = np.where(object_ == 0)
   
    while True:
        n_count += 1 
        C_b = np.copy(C)
        for i in np.unique(non_cte[0]):
            
            for j in non_cte[1][non_cte[0] == i]:
                if j == 0:
                    C[i,0] = w*(C[i+1,0] + C[i-1,0] + C[i,1] + C[i,-2])
                elif j == (C.shape[0]-1):
                    C[i,-1] = w*(C[i+1,-1] + C[i-1,-1] + C[i,1] + C[i,-2])
                else:
                    C[i,j] = w*(C[i+1,j] + C[i-1,j] + C[i,j+1] + C[i,j-1])
            
            

        if np.allclose(C, C_b, atol=stop):
            yield C
            break
        if n_count%store_step == 0:
            yield C


def sor(C,object_,w,stop,store_step):
    """
    Performs Successive over relaxation.
    Inputs:
        - A: Matrix A with all values
        - w: weight
        - stop: simulation stopper
    """
    n_count = 0
    non_cte = np.where(object_ == 0)
    while True:
        n_count += 1 
        C_b = np.copy(C)
        for i in np.unique(non_cte[0]):
            
            for j in non_cte[1][non_cte[0] == i]:
                if j == 0:
                    C[i,0] = (w/4)*(C[i+1,0] + C[i-1,0] + C[i,1] + C[i,-2]) + (1-w)*C[i,0]
                elif j == (C.shape[0]-1):
                    C[i,-1] = (w/4)*(C[i+1,-1] + C[i-1,-1] + C[i,1] + C[i,-2]) + (1-w)*C[i,-1]
                else:
                    C[i,j] = (w/4)*(C[i+1,j] + C[i-1,j] + C[i,j+1] + C[i,j-1])+ (1-w)*C[i,j]
            

        if np.allclose(C, C_b, atol=stop):
            yield C
            break
        if n_count%store_step == 0:
            yield C
    



            

