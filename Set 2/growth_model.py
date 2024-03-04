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

        sink = np.where(cluster[n_count] == 1)

        for k in range(len(sink[0])):
            a = sink[0][k]
            b = sink[1][k]

            # recording the concentration of four candidates
            if a == 0 and b == 0:
                u = 0
                d = C[n_count][a + 1][b] ** eta
                l = 0
                r = C[n_count][a][b + 1] ** eta
            elif a == N and b == 0:
                u = C[n_count][a - 1][b] ** eta
                d = 0
                l = 0
                r = C[n_count][a][b + 1] ** eta
            elif a == 0 and b == N:
                u = 0
                d = C[n_count][a + 1][b] ** eta
                l = C[n_count][a][b - 1] ** eta
                r = 0
            elif a == N and b == N:
                u = C[n_count][a - 1][b] ** eta
                d = 0
                l = C[n_count][a][b - 1] ** eta
                r = 0
            elif a == 0:
                u = 0
                d = C[n_count][a + 1][b] ** eta
                l = C[n_count][a][b - 1] ** eta
                r = C[n_count][a][b + 1] ** eta
            elif a == N:
                u = C[n_count][a - 1][b] ** eta
                d = 0
                l = C[n_count][a][b - 1] ** eta
                r = C[n_count][a][b + 1] ** eta
            elif b == 0:
                u = C[n_count][a - 1][b] ** eta
                d = C[n_count][a + 1][b] ** eta
                l = 0
                r = C[n_count][a][b + 1] ** eta
            elif b == N:
                u = C[n_count][a - 1][b] ** eta
                d = C[n_count][a + 1][b] ** eta
                l = C[n_count][a][b - 1] ** eta
                r = 0
            else:
                u = C[n_count][a - 1][b] ** eta
                d = C[n_count][a + 1][b] ** eta
                l = C[n_count][a][b - 1] ** eta
                r = C[n_count][a][b + 1] ** eta

            # calculating the probability of becoming a sink for candidates
            if (d + u + l + r) != 0:

                for i in [a - 1, a + 1]:
                    if (C[n_count][i][b] ** eta)/ (d + u + l + r) > random():
                        C[n_count][i][b] = 0
                        cluster[n_count][i][b] = 1

                for j in [b - 1, b + 1]:
                    if (C[n_count][a][j] ** eta) / (d + u + l + r) > random():
                        C[n_count][a][j] = 0
                        cluster[n_count][a][j] = 1

        return [C, cluster, n_count]