"""DEFINE JACOBI,etc"""


import numpy as np

def gauss_seidel(A, stop):
    """
    Performs Gauss Seidel iterations
    Inputs:
        - A: Matrix with initial conditions/positions
        - Stop: Stop criteria
    """
    c = 1/4
    while True:
        A_k = A
        A_kb = np.copy(A_k)
        for i in A.size[0]:
            for j in A.size[1]:
                A_k[i,j] = c*(A_k[i+1,j] + A_k[i-1,j] + A_k[i,j+1], A_k[i,j-1])

        if np.abs(A_k-A_kb).all() <= stop:
            break
        A_kb = np.copy(A_k)

    return A_k

def sor(A,w,stop):
    c = w
    while True:
        A_k = A
        A_kb = np.copy(A_k)
        for i in A.size[0]:
            for j in A.size[1]:
                A_k[i,j] = c*(A_k[i+1,j] + A_k[i-1,j] + A_k[i,j+1], A_k[i,j-1]) + (1-w)*A_k[i,j]

        if np.abs(A_k-A_kb).all() <= stop:
            break
        A_kb = np.copy(A_k)



            

