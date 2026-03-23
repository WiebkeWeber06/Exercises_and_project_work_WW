# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 13:32:39 2026

@author: wiewe372
"""

import numpy as np

# a
# --------------------------------------------
a = np.zeros(10)
a[4] = 1
print(a)

# b
# --------------------------------------------
b = np.arange(10, 50, 1)
print(b)

# c
# --------------------------------------------
c = b[::-1]
print(c)

# d
# --------------------------------------------
d = np.arange(0, 9, 1).reshape(3, 3)
print(d)

# e
# --------------------------------------------
e = np.array([1, 2, 0, 0, 4, 0])
e2 = np.nonzero(e)
print(e2)

# f
# --------------------------------------------
# vector of size 30 with random elements
f = np.random.random((30))
print(f)

# check if size is 30
print(len(f))

# get mean value of vector f
f_mean = f.mean()
print(f_mean)

# g
# --------------------------------------------
g = np.zeros((4, 4))
g_padded = np.pad(g, pad_width=1, mode='constant', constant_values = 1)
print(g_padded)

# h & i
# --------------------------------------------
h = np.array((1, 0))
h2 = np.tile(h, (8, 4))
print(h2)

# h different
# --------------------------------------------
i = np.zeros((8, 8))

i[::2, ::2] = 1
i[1::2, 1::2] = 1
    
print(i)    

# j
# --------------------------------------------
j = np.arange(11)
j[4:8:] = -1 * j[4:8:]

print(j)

# k
# --------------------------------------------
k = np.random.random(10)
k.sort()

print(k)

# l
# --------------------------------------------
A = np.random.randint(0,2,5)
B = np.random.randint(0,2,5)

equal = A == B

print(equal)

# m
# --------------------------------------------
m =np.arange(10, dtype=np.int32)
print(m.dtype)
print(m)

np.square(m, m, dtype=np.int32)
print(m.dtype)
print(m)

# n
# --------------------------------------------
n = np.arange(9).reshape(3,3)
print(n)
o = n + 1
print(o)

p = np.dot(n,o)
print(p)

q = p.diagonal()
print(q)





















