import numpy as np
from random import random
<<<<<<< HEAD


def growth_model(N, position, w, eta, stop):
    # N is the grid size
    # position is the initial place with a sink

    # C is the 3-dimensional matrix for iteration
    C = np.zeros((1000, N + 3, N + 3))

    # Initial condition
    for i in range(len(C)):
        C[i][0][:] = 1

    # cluster is the 3-dimensional matrix for updating the sinks
    cluster = np.zeros((len(C), N + 3, N + 3))

    # initial condition for the sink
    cluster[0][position[0] + 1][position[1] + 1] = 1

    n_count = 0

    # SOR
    while True:

        n_count += 1
        cluster[n_count] = cluster[n_count - 1]
        non_cte = np.where(cluster[n_count] == 0)

        for i in range(len(non_cte[0])):

            if non_cte[0][i] == 0 or non_cte[0][i] == N + 2:
                pass
            elif non_cte[1][i] == N + 2:
                C[n_count][non_cte[0][i]][N + 2] = w * 0.25 * (
                                C[n_count][non_cte[0][i] - 1][N + 2] + C[n_count - 1][non_cte[0][i] + 1][N + 2] +
                                C[n_count - 1][non_cte[0][i]][1] + C[n_count][non_cte[0][i]][N + 1]) + (1 - w) * \
                                                       C[n_count - 1][non_cte[0][i]][N + 2]
            else:
                C[n_count][non_cte[0][i]][non_cte[1][i]] = w * 0.25 * (
                                C[n_count][non_cte[0][i] - 1][non_cte[1][i]] + C[n_count - 1][non_cte[0][i] + 1][
                            non_cte[1][i]] + C[n_count - 1][non_cte[0][i]][non_cte[1][i] + 1] +
                                C[n_count][non_cte[0][i]][non_cte[1][i] - 1]) + (1 - w) * C[n_count - 1][non_cte[0][i]][
                                                                   non_cte[1][i]]

        if np.max(np.abs(C[n_count]-C[n_count-1])) < stop: 
            break
            
        # Finding all the candidates
        sink = np.where(cluster[n_count] == 1)
        candidate_0 = []
        candidate_1 = []
        c_candidate = 0
        tag = np.zeros((N + 3, N + 3))
=======
import matplotlib.pyplot as plt
from numba import jit


@jit(nopython=True, parallel = False)
def sor(C, cluster,w,stop = 0.0001):
    """
    Performs Successive over relaxation.
    Inputs:
        - C: Matrix A with all values
        - cluster: Matrix with all the sink and cluster positions
        - w: weight
        - stop: simulation stopper
    """

    
    A = C.copy()
    n_count = 0
    non_cte = np.where(cluster == 0)

    while True:
        n_count += 1 
        A_b = A.copy()
        for i in np.unique(non_cte[0]):
            for j in non_cte[1][non_cte[0] == i]:
                if j == 1:
                    A[i,1] = (w/4)*(A[i+1,1] + A[i-1,1] + A[i,2] + A[i,-3]) + (1-w)*A[i,1]

                elif j == (A.shape[0]-2):
                    A[i,-2] = (w/4)*(A[i+1,-2] + A[i-1,-2] + A[i,2] + A[i,-3]) + (1-w)*A[i,-2]
                else:
                    A[i,j] = (w/4)*(A[i+1,j] + A[i-1,j] + A[i,j+1] + A[i,j-1])+ (1-w)*A[i,j]
        
        max_diff = np.max(np.abs(A - A_b))
        
        if  max_diff <= stop:
            return A.copy()
            break

@jit(nopython=True, parallel = False)  
def growth_model(N, position,w,eta, grow_steps = 10, D = 1):
    """

    Executes a time_step ing method given a function

    Inputs:

        - N: Number of divisions
        - w: The SOR parameter
        - eta: eta I Imagine
        - grow_steps: Number of total steps it will do, default + 1000
    """
    C = np.zeros((N+3, N+3)) ## Matrix with concentration
    number = 20
    C[1,:] = 1
    
    C[0,:] = number
    C[-1,:] = number
    C[:,0] = number
    C[:,-1] = number

    data = [C.copy()]
    
    
    cluster = np.zeros((N+3, N+3)) ## Matrix with points. 1 is limit, 2 is sink point, 0 is nothing

    cluster[1,:] = 1
    cluster[-2,:] = 1

    cluster[0,:] = number
    cluster[-1,:] = number
    cluster[:,0] = number
    cluster[:,-1] = number

    cluster[-2,position] = 2

    n_count = 0
    for i in range(1,grow_steps):
        C = sor(C,cluster,w)
        
        # Finding all the candidates
        sink = np.where(cluster == 2)
        candidate_0 = [] # ROW
        candidate_1 = [] # Column
        c_candidate = 0
        tag = np.zeros((N + 3, N + 3)) 
>>>>>>> 16952a5514cfb7a02e6bdeef105963577a71f82e

        for k in range(len(sink[0])):
            a = sink[0][k]
            b = sink[1][k]

            if a == 0 or a == N + 2 or b == 0 or b == N + 2:
<<<<<<< HEAD
                pass
            else:
                for i in [a - 1, a + 1]:
                    if cluster[n_count][i][b] == 0 and tag[i][b] == 0:
                        candidate_0.append(i)
                        candidate_1.append(b)
                        tag[i][b] = 1
                        c_candidate += C[n_count][candidate_0[-1]][candidate_1[-1]]
                for j in [b - 1, b + 1]:
                    if cluster[n_count][a][j] == 0 and tag[a][j] == 0:
                        candidate_0.append(a)
                        candidate_1.append(j)
                        tag[a][j] = 1
                        c_candidate += C[n_count][candidate_0[-1]][candidate_1[-1]]

        # calculating the probability of becoming a sink for candidates
        for k in range(len(candidate_0)):

            if (C[n_count][candidate_0[k]][candidate_1[k]] ** eta) / c_candidate > random():
                C[n_count][candidate_0[k]][candidate_1[k]] = 0
                cluster[n_count][candidate_0[k]][candidate_1[k]] = 1

    return [C, cluster, n_count]

=======
                continue
            else:
                for i in [a - 1, a + 1]:
                    if cluster[i][b] == 0 and tag[i][b] == 0:
                        candidate_0.append(i)
                        candidate_1.append(b)
                        tag[i][b] = 1
                        c_candidate += C[candidate_0[-1]][candidate_1[-1]]
                for j in [b - 1, b + 1]:
                    if cluster[a][j] == 0 and tag[a][j] == 0:
                        candidate_0.append(a)
                        candidate_1.append(j)
                        tag[a][j] = 1
                        c_candidate += C[candidate_0[-1]][candidate_1[-1]]

        # calculating the probability of becoming a sink for candidates
        for k in range(len(candidate_0)):
            if (C[candidate_0[k]][candidate_1[k]] ** eta) / c_candidate > np.random.uniform(0,1):
                C[candidate_0[k]][candidate_1[k]] = 0
                cluster[candidate_0[k]][candidate_1[k]] = 2
        data.append(C.copy())
        

    return C, cluster, data

>>>>>>> 16952a5514cfb7a02e6bdeef105963577a71f82e

def merge(C, cluster, n_count):
    # C and cluster are 3-dimensional matrices

    sink = np.where(cluster[n_count] == 1)

    for i in range(len(sink[0])):
        C[n_count][sink[0][i]][sink[1][i]] = 1

    return C[n_count]
<<<<<<< HEAD
=======

if __name__ == "__main__":
    N = 10
    w = 1.9
    eta = 0.8
    result = growth_model(N, 3,w,eta)
    data = result[2]
    
    fig, ax = plt.subplots()
    ax.imshow(result[0][1:-1, 1:-1], cmap='hot', interpolation='nearest', extent=[0, 1, 0, 1])
    plt.show()  


    
>>>>>>> 16952a5514cfb7a02e6bdeef105963577a71f82e
