# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 13:35:05 2026

@author: wiewe372
"""

from line_profiler import profile

@profile
def matmult_original(X, Y):
    result = [[0] * len(Y[0]) for _ in range(len(X))]

    for i in range(len(X)):
        for j in range(len(Y[0])):
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]

    return result


def main():
    import random
    N = 250

    X = [[random.randint(0, 100) for _ in range(N)] for _ in range(N)]
    Y = [[random.randint(0, 100) for _ in range(N + 1)] for _ in range(N)]

    matmult_original(X, Y)


if __name__ == "__main__":
    main()