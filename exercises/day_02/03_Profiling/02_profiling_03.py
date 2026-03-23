# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 10:09:06 2026

@author: wiewe372
"""

# a) matmult.py
# I would start optimizing the line
#     result[i][j] += X[i][k] * Y[k][j]
# because it is the innermost operation inside a triple nested loop.
# Profiling would show that most of the runtime is spent there.

# b) euler72.py
# I would start optimizing the repeated call
#     fraq += fast_phi(i, primes)
# because fast_phi is called for every number from 2 to 10000.
# Inside fast_phi, the main cost comes from
#     factors = factorize(n, primes)
# and inside factorize from:
#     for p in primes:
#         while(n % p == 0):
#             n = n / p
#             factors.append(p)
#         if(p > sqrt(n)):
#             break
# In gen_primes, the expensive part is also the repeated loop
# over primes together with repeated sqrt computation:
#     for d in primes:
#         if(d > sqrt(l[j])):
#

# Profiling tools
# ---------------------------------------------

# I profiled the scripts using cProfile, line_profiler, and Scalene.

# cProfile:
# Used to identify which functions consume the most total runtime.
# Command:
#   python -m cProfile -s cumulative matmult.py
# Result:
# The matrix multiplication function dominated runtime.

# line_profiler:
# Used to analyze performance at the line level.
# I added @profile decorators to functions and ran:
#   kernprof -l -v matmult.py
# Result:
# The innermost line
#     result[i][j] += X[i][k] * Y[k][j]
# was identified as the main bottleneck.

# scalene:
# Used to distinguish between Python time and native code time.
# Command:
#   scalene matmult.py
# Result:
# Most time was spent in Python loops in the original version,
# while the NumPy version shifted computation to native code.

# In Jupyter:
# I used:
#   %prun
#   %lprun
# to interactively profile functions.

# -----------------------------------------------

# c: Improve the performance of matmult.py

import random
import time
import numpy as np

N = 250

# -----------------------------
# Generate matrices ONCE
# -----------------------------
X = [[random.randint(0, 100) for _ in range(N)] for _ in range(N)]
Y = [[random.randint(0, 100) for _ in range(N + 1)] for _ in range(N)]

# Convert to numpy once (same data!)
X_np = np.array(X)
Y_np = np.array(Y)


# -----------------------------
# 1. ORIGINAL VERSION
# -----------------------------
def matmult_original(X, Y):
    N = len(X)
    result = [[0] * len(Y[0]) for _ in range(N)]

    for i in range(len(X)):
        for j in range(len(Y[0])):
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]

    return result


# -----------------------------
# 2. OPTIMIZED PURE PYTHON
# -----------------------------
def matmult_optimized(X, Y):
    N = len(X)
    result = [[0] * len(Y[0]) for _ in range(N)]

    YT = list(zip(*Y))  # transpose once

    for i, Xi in enumerate(X):
        ri = result[i]
        for j, Yj in enumerate(YT):
            s = 0
            for a, b in zip(Xi, Yj):
                s += a * b
            ri[j] = s

    return result


# -----------------------------
# 3. NUMPY VERSION
# -----------------------------
def matmult_numpy(X, Y):
    return X @ Y


# -----------------------------
# TIMING FUNCTION
# -----------------------------
def measure(func, *args):
    start = time.perf_counter()
    result = func(*args)
    end = time.perf_counter()
    return result, end - start


# -----------------------------
# RUN BENCHMARKS
# -----------------------------
print("Running benchmarks...\n")

res1, t1 = measure(matmult_original, X, Y)
print(f"Original version:   {t1:.4f} seconds")

res2, t2 = measure(matmult_optimized, X, Y)
print(f"Optimized version:  {t2:.4f} seconds")

res3, t3 = measure(matmult_numpy, X_np, Y_np)
print(f"NumPy version:      {t3:.4f} seconds")


# -----------------------------
# PRINT RESULTS
# -----------------------------
print("\n--- Sample of result (first 3 rows) ---\n")

print("Original:")
for row in res1[:3]:
    print(row[:10])  # only first 10 elements

print("\nOptimized:")
for row in res2[:3]:
    print(row[:10])

print("\nNumPy:")
for row in res3[:3]:
    print(row[:10])


# -----------------------------
# VERIFY CORRECTNESS
# -----------------------------
print("\nChecking correctness...")
print("Original vs optimized:", res1 == res2)
print("Original vs numpy:   ", res1 == res3.tolist())


# Original version:   1.0089 seconds
# Optimized version:  0.5618 seconds
# Best performance achieved for N = 250:
# NumPy version: 0.0063 seconds










