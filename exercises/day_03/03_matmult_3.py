# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 14:59:21 2026

@author: wiewe372
"""

# Program to multiply two matrices using nested loops
import random
import numpy as np

N = 3





# NxN matrix
# --------------------------------------------
X = []
for i in range(N):
    X.append([random.randint(0,100) for r in range(N)])
    
print(X)

# --------------------
# with numpy
X1 = np.random.randint(0, 100, (N, N))
print(X1)





# Nx(N+1) matrix
# --------------------------------------------
Y = []
for i in range(N):
    Y.append([random.randint(0,100) for r in range(N+1)])

print(Y)

# result is Nx(N+1)
result = []
for i in range(N):
    result.append([0] * (N+1))

# iterate through rows of X
for i in range(len(X)):
    # iterate through columns of Y
    for j in range(len(Y[0])):
        # iterate through rows of Y
        for k in range(len(Y)):
            print(X[i][k]),
            print(Y[k][j]),
            print(X[i][k] * Y[k][j])
            result[i][j] += X[i][k] * Y[k][j]
            print(result[i][j])

for r in result:
    print(r)
    
    
    
    
    
# --------------------
# with numpy
Y1 = np.random.randint(0, 100, (N, N+1))
print(Y1)

Y2 = np.zeros((N, N+1))
print(Y2)

Z = np.matmul(X, Y, out=Y2)
print(Z)

