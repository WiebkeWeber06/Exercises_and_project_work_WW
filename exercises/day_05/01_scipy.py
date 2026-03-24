# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:03:19 2026

@author: wiewe372
"""

import numpy as np
from scipy import linalg

# a) Define matrix A
A = np.array([
    [1, -2, 3],
    [4,  5, 6],
    [7,  1, 9]
], dtype=float)

# b) Define vector b
b = np.array([1, 2, 3], dtype=float)

print("Matrix A:")
print(A)
print("\nVector b:")
print(b)

# c) Solve A x = b
x = linalg.solve(A, b)

print("\nSolution x of A x = b:")
print(x)

# d) Check solution
check = A @ x
print("\nCheck A @ x:")
print(check)
print("Matches b:", np.allclose(check, b))

# e) Repeat using a random 3x3 matrix B instead of vector b
B = np.random.rand(3, 3)

print("\nRandom 3x3 matrix B:")
print(B)

# Solve A X = B
X = linalg.solve(A, B)

print("\nSolution X of A X = B:")
print(X)

# Check solution
check_B = A @ X
print("\nCheck A @ X:")
print(check_B)
print("Matches B:", np.allclose(check_B, B))

# f) Eigenvalue problem for A
eigenvalues, eigenvectors = linalg.eig(A)

print("\nEigenvalues of A:")
print(eigenvalues)

print("\nEigenvectors of A (columns):")
print(eigenvectors)

# g) Inverse and determinant of A
A_inv = linalg.inv(A)
A_det = linalg.det(A)

print("\nInverse of A:")
print(A_inv)

print("\nDeterminant of A:")
print(A_det)

# h) Norms of A with different orders
print("\nNorms of A:")
print("Frobenius norm:", linalg.norm(A))
print("1-norm:", linalg.norm(A, 1))
print("Infinity norm:", linalg.norm(A, np.inf))
print("2-norm:", linalg.norm(A, 2))