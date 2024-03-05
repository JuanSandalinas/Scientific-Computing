import numpy as np
from random import random

def growth_model(N, position, w, eta, stop):
    # N is the grid size
    # position is the initial place with a sink


    # C is the 3-dimensional matrix for iteration
    C = np.zeros((N + 1, N + 1, N + 1))

    # Initial condition
    for i in range(N + 1):
        C[i][0][:] = 1

    # cluster is the 3-dimensional matrix for updating the sinks
    cluster = np.zeros((N + 1, N + 1, N + 1))

    # initial condition for the sink
    for i in range(len(position[0])):
        cluster[0][position[0][i]][position[1][i]] = 1

    n_count = 0
    measure = 1

    # SOR
    while measure > stop:

        n_count += 1
        cluster[n_count] = cluster[n_count - 1]
        non_cte = np.where(cluster[n_count] == 0)

        for i in range(len(non_cte[0])):

            if non_cte[0][i] == 0 or non_cte[0][i] == N:
                measure_n = measure
            else:
                if non_cte[1][i] == 0:
                    C[n_count][non_cte[0][i]][0] = w * 0.25 * (
                                C[n_count][non_cte[0][i] - 1][0] + C[n_count - 1][non_cte[0][i] + 1][0] +
                                C[n_count - 1][non_cte[0][i]][1] + C[n_count][non_cte[0][i]][-1]) + (1 - w) * \
                                                   C[n_count - 1][non_cte[0][i]][0]
                    measure_n = C[n_count][non_cte[0][i]][0] - C[n_count - 1][non_cte[0][i]][0]
                elif non_cte[1][i] == N:
                    C[n_count][non_cte[0][i]][N] = w * 0.25 * (
                                C[n_count][non_cte[0][i] - 1][N] + C[n_count - 1][non_cte[0][i] + 1][N] +
                                C[n_count - 1][non_cte[0][i]][1] + C[n_count][non_cte[0][i]][N - 1]) + (1 - w) * \
                                                   C[n_count - 1][non_cte[0][i]][N]
                    measure_n = C[n_count][non_cte[0][i]][N] - C[n_count - 1][non_cte[0][i]][N]
                else:
                    C[n_count][non_cte[0][i]][non_cte[1][i]] = w * 0.25 * (
                                C[n_count][non_cte[0][i] - 1][non_cte[1][i]] + C[n_count - 1][non_cte[0][i] + 1][
                            non_cte[1][i]] + C[n_count - 1][non_cte[0][i]][non_cte[1][i] + 1] +
                                C[n_count][non_cte[0][i]][non_cte[1][i] - 1]) + (1 - w) * C[n_count - 1][non_cte[0][i]][
                                                                   non_cte[1][i]]
                    measure_n = C[n_count][non_cte[0][i]][non_cte[1][i]] - C[n_count - 1][non_cte[0][i]][non_cte[1][i]]

            if measure_n > measure:
                measure = measure_n

        # changing the concentration of sinks into 0
        for i in range(N):
            for j in range(N):
                if cluster[n_count][i][j] == 1:
                    C[n_count][i][j] = 0

        # Finding all the candidates
        sink = np.where(cluster[n_count] == 1)
        candidate = np.zeros((N + 1, N + 1))
        c_candidate = 0
        tag = np.zeros((N + 1, N + 1))

        for k in range(len(sink[0])):
            a = sink[0][k]
            b = sink[1][k]

            if a == 0 and b == 0:
                if sink[0][1] != 1 and tag[0][1] == 0:
                    candidate[0].append(0)
                    candidate[1].append(b+1)
                    tag[0][1] = 1
                if sink[1][0] != 1 and tag[1][0] == 0:
                    candidate[0].append(a+1)
                    candidate[1].append(b)
                    tag[1][0] = 1
            elif a == N and b == 0:
                if sink[N-1][0] != 1 and tag[N-1][0] == 0:
                    candidate[0].append(N-1)
                    candidate[1].append(0)
                    tag[N-1][0] = 1
                if sink[N][1] != 1 and tag[N][1] == 0:
                    candidate[0].append(N)
                    candidate[1].append(1)
                    tag[N][1] = 1
            elif a == 0 and b == N:
                if sink[0][N-1] != 1 and tag[0][N-1] == 0:
                    candidate[0].append(0)
                    candidate[1].append(N-1)
                    tag[0][N-1] = 1
                if sink[1][N] != 1 and tag[1][N] == 0:
                    candidate[0].append(1)
                    candidate[1].append(N)
                    tag[1][N] = 1
            elif a == N and b == N:
                if sink[N][N-1] != 1 and tag[N][N-1] == 0:
                    candidate[0].append(N)
                    candidate[1].append(N-1)
                    tag[N][N-1] = 1
                if sink[N-1][N] != 1 and tag[N-1][N] == 0:
                    candidate[0].append(N-1)
                    candidate[1].append(N)
                    tag[N-1][N] = 1
            elif a == 0:
                if sink[0][b-1] != 1 and tag[0][b-1] == 0:
                    candidate[0].append(0)
                    candidate[1].append(b-1)
                    tag[1][b-1] = 1
                if sink[0][b+1] != 1 and tag[0][b+1] == 0:
                    candidate[0].append(0)
                    candidate[1].append(b+1)
                    tag[0][b+1] = 1
                if sink[1][b] != 1 and tag[1][b] == 0:
                    candidate[0].append(1)
                    candidate[1].append(b)
                    tag[1][b] = 1
            elif a == N:
                if sink[N][b-1] != 1 and tag[0][b-1] == 0:
                    candidate[0].append(0)
                    candidate[1].append(b-1)
                    tag[0][b-1] = 1
                if sink[N][b+1] != 1 and tag[0][b+1] == 0:
                    candidate[0].append(0)
                    candidate[1].append(b+1)
                    tag[0][b+1] = 1
                if sink[N-1][b] != 1 and tag[N-1][b] == 0:
                    candidate[0].append(N-1)
                    candidate[1].append(b)
                    tag[N-1][b] = 1
            elif b == 0:
                if sink[a-1][0] != 1 and tag[a-1][0] == 0:
                    candidate[0].append(a-1)
                    candidate[1].append(0)
                    tag[a-1][0] = 1
                if sink[a+1][0] != 1 and tag[a+1][0] == 0:
                    candidate[0].append(a+1)
                    candidate[1].append(0)
                    tag[a+1][0] = 1
                if sink[a][1] != 1 and tag[a][1] == 0:
                    candidate[0].append(a)
                    candidate[1].append(1)
                    tag[a][1] = 1
            elif b == N:
                if sink[a][N-1] != 1 and tag[a][N-1] == 0:
                    candidate[0].append(a)
                    candidate[1].append(N-1)
                    tag[a][N-1] = 1
                if sink[a-1][N] != 1 and tag[a-1][N] == 0:
                    candidate[0].append(a-1)
                    candidate[1].append(N)
                    tag[a-1][b-1] = 1
                if sink[a+1][N] != 1 and tag[a+1][N] == 0:
                    candidate[0].append(a+1)
                    candidate[1].append(N)
                    tag[a+1][N] = 1
            else:
                if sink[a-1][b] != 1 and tag[a-1][b] == 0:
                    candidate[0].append(a-1)
                    candidate[1].append(b)
                    tag[a-1][b] = 1
                if sink[a+1][b] != 1 and tag[a+1][b] == 0:
                    candidate[0].append(a+1)
                    candidate[1].append(b)
                    tag[a+1][b] = 1
                if sink[a][b-1] != 1 and tag[a][b-1] == 0:
                    candidate[0].append(a)
                    candidate[1].append(b-1)
                    tag[a][b-1] = 1
                if sink[a][b+1] != 1 and tag[a][b+1] == 0:
                    candidate[0].append(a)
                    candidate[1].append(b+1)
                    tag[a][b+1] = 1


            c_candidate += C[n_count][candidate[0][-1]][candidate[1][-1]]

        # calculating the probability of becoming a sink for candidates
        for k in range(len(candidate[0])):

            if (C[n_count][candidate[0]][candidate[1]] ** eta)/ c_candidate > random():
                C[n_count][candidate[0]][candidate[1]] = 0
                cluster[n_count][candidate[0]][candidate[1]] = 1

        return [C, cluster, n_count]