# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:05:21 2026

@author: wiewe372
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# -----------------------------
# a) Discrete random variable: Poisson distribution
# -----------------------------
mu = 4  # mean of Poisson distribution
poisson_rv = stats.poisson(mu=mu)

# x values for plotting
x_pois = np.arange(0, 15)

# PMF and CDF
pmf_pois = poisson_rv.pmf(x_pois)
cdf_pois = poisson_rv.cdf(x_pois)

# 1000 random realizations
samples_pois = poisson_rv.rvs(size=1000)

# Plot PMF
plt.figure(figsize=(6, 4))
plt.stem(x_pois, pmf_pois)
plt.title("Poisson Distribution PMF")
plt.xlabel("x")
plt.ylabel("P(X = x)")
plt.grid(True)
plt.show()

# Plot CDF
plt.figure(figsize=(6, 4))
plt.step(x_pois, cdf_pois, where="post")
plt.title("Poisson Distribution CDF")
plt.xlabel("x")
plt.ylabel("P(X ≤ x)")
plt.grid(True)
plt.show()

# Histogram of samples
plt.figure(figsize=(6, 4))
plt.hist(samples_pois, bins=np.arange(-0.5, max(samples_pois)+1.5, 1), density=True, alpha=0.7, edgecolor="black")
plt.title("Histogram of 1000 Poisson Random Samples")
plt.xlabel("x")
plt.ylabel("Relative Frequency")
plt.grid(True)
plt.show()


# -----------------------------
# b) Continuous random variable: Normal distribution
# -----------------------------
mu = 0       # mean
sigma = 1    # standard deviation
normal_rv = stats.norm(loc=mu, scale=sigma)

# x values for plotting
x_norm = np.linspace(-4, 4, 500)

# PDF and CDF
pdf_norm = normal_rv.pdf(x_norm)
cdf_norm = normal_rv.cdf(x_norm)

# 1000 random realizations
samples_norm = normal_rv.rvs(size=1000)

# Plot PDF
plt.figure(figsize=(6, 4))
plt.plot(x_norm, pdf_norm)
plt.title("Normal Distribution PDF")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.show()

# Plot CDF
plt.figure(figsize=(6, 4))
plt.plot(x_norm, cdf_norm)
plt.title("Normal Distribution CDF")
plt.xlabel("x")
plt.ylabel("P(X ≤ x)")
plt.grid(True)
plt.show()

# Histogram of samples
plt.figure(figsize=(6, 4))
plt.hist(samples_norm, bins=30, density=True, alpha=0.7, edgecolor="black")
plt.title("Histogram of 1000 Normal Random Samples")
plt.xlabel("x")
plt.ylabel("Density")
plt.grid(True)
plt.show()


# -----------------------------
# c) Test whether two independent datasets
#    come from the same distribution
# -----------------------------
# Example datasets
data1 = np.random.normal(loc=0, scale=1, size=100)
data2 = np.random.normal(loc=0, scale=1, size=100)

# Independent t-test
t_stat, p_value = stats.ttest_ind(data1, data2)

print("t-statistic:", t_stat)
print("p-value:", p_value)

if p_value < 0.05:
    print("The two datasets are significantly different.")
else:
    print("No significant difference found between the two datasets.")