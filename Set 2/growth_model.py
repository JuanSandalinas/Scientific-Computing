import numpy as np
from random import random


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
                continue
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
        candidate_0 = [] # ROW
        candidate_1 = [] # Coumn
        c_candidate = 0
        tag = np.zeros((N + 3, N + 3)) ##

        for k in range(len(sink[0])):
            a = sink[0][k]
            b = sink[1][k]

            if a == 0 or a == N + 2 or b == 0 or b == N + 2:
                continue
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


def merge(C, cluster, n_count):
    # C and cluster are 3-dimensional matrices

    sink = np.where(cluster[n_count] == 1)

    for i in range(len(sink[0])):
        C[n_count][sink[0][i]][sink[1][i]] = 1

    return C[n_count]
