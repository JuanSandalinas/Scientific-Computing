import numpy as np
import matplotlib.pyplot as plt
import random

def growth_model(N, position, w, eta, stop):
    if (position[0]) >= (N - 1) or (position[1]) >= (N - 1):
        raise Exception("Outside bounds")

    C = np.zeros((N + 1, N + 1, N + 1))

    C[0][0][:] = 1

    n_count = 0

    cluster = np.zeros((N + 1, N + 1, N + 1))

    cluster[0][position[0]][position[1]] = 1
    C[0][position[0]][position[1]] = 1

    measure = 1

    while measure > stop:

        n_count += 1
        cluster[n_count] = cluster[n_count - 1]
        non_cte = np.where(cluster[n_count] == 0)

        for i in range(len(non_cte)):

            if non_cte[0][i] == 0 or non_cte[0][i] == N:
                pass
                measure_n = measure
            elif non_cte[1][i] == 0:
                C[n_count][non_cte[0][i]][non_cte[1][i]] = w * 0.25 * (
                            C[n_count][non_cte[0][i] - 1][0] + C[n_count - 1][non_cte[0][i] + 1][0] +
                            C[n_count - 1][non_cte[0][i]][2] + C[n_count][non_cte[0][i]][non_cte[1][-1]]) + (1 - w) * \
                                                           C[n_count - 1][non_cte[0][i]][0]
                measure_n = C[n_count][non_cte[0][i]][0] - C[n_count - 1][non_cte[0][i]][non_cte[1][i]]
            elif non_cte[1][i] == N:
                C[n_count][non_cte[0][i]][N] = w * 0.25 * (
                            C[n_count][non_cte[0][i] - 1][N] + C[n_count - 1][non_cte[0][i] + 1][N] +
                            C[n_count][non_cte[0][i]][N + 1] + C[n_count - 1][non_cte[0][i]][N - 1]) + (1 - w) * \
                                               C[n_count - 1][non_cte[0][i]][N]
                measure_n = C[n_count][non_cte[0][i]][N] - C[n_count - 1][non_cte[0][i]][N]
            else:
                C[n_count][non_cte[0][i]][non_cte[1][i]] = w * 0.25 * (
                            C[n_count][non_cte[0][i] - 1][non_cte[1][i]] + C[n_count - 1][non_cte[0][i] + 1][
                        non_cte[1][i]] + C[n_count - 1][non_cte[0][i]][non_cte[1][i] + 1] + C[n_count][non_cte[0][i]][
                                non_cte[1][i] - 1]) + (1 - w) * C[n_count - 1][non_cte[0][i]][non_cte[1][i]]
                measure_n = C[n_count][non_cte[0][i]][non_cte[1][i]] - C[n_count - 1][non_cte[0][i]][non_cte[1][i]]

            if measure_n > measure:
                measure = measure_n

        for i in range(N):
            for j in range(N):
                if cluster[n_count][i][j] == 1:
                    C[n_count][i][j] = 0

        sink = np.where(cluster[n_count] == 1)

        for k in range(len(sink[0])):
            a = sink[0][k]
            b = sink[1][k]

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

            if a != 0:
                if cluster[n_count][a - 1][b] == 1:
                    u = 0
            if a != N:
                if cluster[n_count][a + 1][b] == 1:
                    d = 0
            if b != 0:
                if cluster[n_count][a][b - 1] == 1:
                    l = 0
            if b != N:
                if cluster[n_count][a][b + 1] == 1:
                    r = 0

            if (d + u + l + r) != 0:
                for i in [a - 1, a + 1]:
                    for j in [b - 1, b + 1]:
                        if C[n_count][i][j] / (d + u + l + r) > random():
                            C[n_count][i][j] = 0
                            cluster[n_count][i][j] = 1

        return C