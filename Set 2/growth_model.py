import numpy as np
from random import random

def growth_model(N, position, w, eta, stop):
    # N is the grid size
    # position is the initial place with a sink


    # C is the 3-dimensional matrix for iteration
    C = np.zeros((N + 3, N + 3, N + 3))

    # Initial condition
    for i in range(N + 3):
        C[i][0][:] = 1

    # cluster is the 3-dimensional matrix for updating the sinks
    cluster = np.zeros((N + 3, N + 3, N + 3))

    # initial condition for the sink
    cluster[0][position[0]+1][position[1]+1] = 1

    n_count = 0
    measure = 1

    # SOR
    while measure > stop::

        n_count += 1
        cluster[n_count] = cluster[n_count - 1]
        non_cte = np.where(cluster[n_count] == 0)

        for i in range(len(non_cte[0])):

            if non_cte[0][i] == 0 or non_cte[0][i] == N+2:
                measure_n = measure
            else:
                if non_cte[1][i] == 0:
                    C[n_count][non_cte[0][i]][0] = w * 0.25 * (C[n_count][non_cte[0][i] - 1][0] + C[n_count - 1][non_cte[0][i] + 1][0] + C[n_count - 1][non_cte[0][i]][1] + C[n_count][non_cte[0][i]][-1]) + (1 - w) * C[n_count - 1][non_cte[0][i]][0]
                    measure_n = C[n_count][non_cte[0][i]][0] - C[n_count - 1][non_cte[0][i]][0]
                elif non_cte[1][i] == N+2:
                    C[n_count][non_cte[0][i]][N+2] = w * 0.25 * (C[n_count][non_cte[0][i] - 1][N+2] + C[n_count - 1][non_cte[0][i] + 1][N+2] + C[n_count - 1][non_cte[0][i]][1] + C[n_count][non_cte[0][i]][N+1]) + (1 - w) * C[n_count - 1][non_cte[0][i]][N+2]
                    measure_n = C[n_count][non_cte[0][i]][N+2] - C[n_count - 1][non_cte[0][i]][N+2]
                else:
                    C[n_count][non_cte[0][i]][non_cte[1][i]] = w * 0.25 * (C[n_count][non_cte[0][i] - 1][non_cte[1][i]] + C[n_count - 1][non_cte[0][i] + 1][non_cte[1][i]] + C[n_count - 1][non_cte[0][i]][non_cte[1][i] + 1] + C[n_count][non_cte[0][i]][non_cte[1][i] - 1]) + (1 - w) * C[n_count - 1][non_cte[0][i]][non_cte[1][i]]
                    measure_n = C[n_count][non_cte[0][i]][non_cte[1][i]] - C[n_count - 1][non_cte[0][i]][non_cte[1][i]]

            if measure_n > measure:
                measure = measure_n


        # Finding all the candidates
        sink = np.where(cluster[n_count] == 1)
        candidate_0 = []
        candidate_1 = []
        c_candidate = 0
        tag = np.zeros((N + 3, N + 3))

        for k in range(len(sink[0])):
            a = sink[0][k]
            b = sink[1][k]

            if a == 0 or a == N+2 or b == 0 or b == N+2:
                pass
            else:
                for i in [a-1, a+1]:
                    if cluster[n_count][i][b] == 0 and tag[i][b] == 0:
                        candidate_0.append(i)
                        candidate_1.append(b)
                        tag[i][b] = 1
                        c_candidate += C[n_count][candidate_0[-1]][candidate_1[-1]]
                for j in [b-1, b+1]:
                    if cluster[n_count][a][j] == 0 and tag[a][j] == 0:
                        candidate_0.append(a)
                        candidate_1.append(j)
                        tag[a][j] = 1
                        c_candidate += C[n_count][candidate_0[-1]][candidate_1[-1]]

        # calculating the probability of becoming a sink for candidates
        for k in range(len(candidate_0)):

            if (C[n_count][candidate_0[k]][candidate_1[k]] ** eta)/ c_candidate > random():
                C[n_count][candidate_0[k]][candidate_1[k]] = 0
                cluster[n_count][candidate_0[k]][candidate_1[k]] = 1

    return [C, cluster, n_count]