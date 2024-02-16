"""DEFINE JACOBI,etc"""


import numpy as np
from scipy.sparse import diags
from scipy.sparse import csr_array
from scipy.special import erfc


def gauss_seidel(A, stop):
    """
    Performs Gauss Seidel iterations
    Inputs:
        - A: Matrix with initial conditions/positions
        - Stop: Stop criteria
    """
    c = 1/4
    A_b = np.ones((A.shape[0], A.shape[1]))
    while True:
        for i in range(1,A.shape[0]-1):
            for j in range(A.shape[1]-1):
                A[i,j] = c*(A[i+1,j] + A[i-1,j] + A[i,j+1] + A[i,j-1])
            A[i,-1] = A[i,0]
        if np.abs(A-A_b).all() <= stop:
            break
        A_b = A

    return A

def sor(A,w,stop):
    c = w
    A_b = np.ones((A.shape[0], A.shape[1]))
    while True:
        for i in range(1,A.shape[0]-1):
            for j in A.size[1]:
                A_k[i,j] = c*(A_k[i+1,j] + A_k[i-1,j] + A_k[i,j+1] + A_k[i,j-1]) + (1-w)*A_k[i,j]
            A[i,-1] = A[i,0]

        if np.abs(A-A_b).all() <= stop:
            break
        A_b = A
    return A


def analytical(time,tao):
    result = 0.0
    for i in range(time):
        arg1 = 1 - y + 2*i
        arg2 = 1 + y + 2*i
        result += erfc(arg1 / (2 * np.sqrt(D * t))) - erfc(arg2 / (2 * np.sqrt(D * t)))
    return result



            

