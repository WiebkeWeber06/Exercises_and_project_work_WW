# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 12:07:37 2026

@author: wiewe372
"""

import numpy as np
import matplotlib.pyplot as plt

# Create grid data
x = np.linspace(-3, 3, 200)
y = np.linspace(-3, 3, 200)
X, Y = np.meshgrid(x, y)

# Create function (2D Gaussian + waves)
Z = np.exp(-(X**2 + Y**2)) * np.sin(3*X) * np.cos(3*Y)

# Create figure
plt.figure(figsize=(8, 6))

# Heatmap
heatmap = plt.imshow(
    Z,
    extent=[x.min(), x.max(), y.min(), y.max()],
    origin='lower',
    cmap='viridis'
)

# Contour lines on top
contours = plt.contour(X, Y, Z, colors='white', linewidths=0.5)

# Colorbar
plt.colorbar(heatmap, label='Function value')

# Labels & title
plt.title("2D Function Visualization (Heatmap + Contours)")
plt.xlabel("X")
plt.ylabel("Y")

# Layout
plt.tight_layout()

# Save (good for submission)
plt.savefig("heatmap_contour.png", dpi=300)

# Show
plt.show()