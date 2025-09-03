# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 18:00:13 2025

@author: zacha
"""

'''FIGURE OUT - how to determine dna concentrations from input proteins to make this 
dependent on input proteins instead of rna conc'''

import numpy as np
import matplotlib.pyplot as plt

# Parameters
k_tl = 0.05
k_deg = 0.0075

# Protein equation
def f(rna1, rna2, p):
    
    return k_tl*(rna1 + rna2) - k_deg*p

# Create a grid in the (rna1, rna2) plane
rna1_values = np.linspace(0, 6, 5)
rna2_values = np.linspace(0, 6, 5)
p_values = np.linspace(0, 80, 4)
R1, R2, P = np.meshgrid(rna1_values, rna2_values, p_values)

# Compute the derivative on the grid
dp = f(R1, R2, P)

# For quiver, we need vectors. We'll show change along R1 and R2 as zeros (since only p changes)
U = np.zeros_like(R1)  # change along rna1 (no change)
V = np.zeros_like(R2)  # change along rna2 (no change)
W = dp                # change along p

# 3D quiver plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(R1, R2, P, U, V, W, length = 5)

ax.set_xlabel('rna1')
ax.set_ylabel('rna2')
ax.set_zlabel('Protein (p)')

# ax.set_xlim(0, 4)
# ax.set_ylim(0, 4)
# ax.set_zlim(10, 30)

# Nullcline surface
rna1_fine = np.linspace(0, 6, 10)
rna2_fine = np.linspace(0, 6, 10)
R1_fine, R2_fine = np.meshgrid(rna1_fine, rna2_fine)
p_null = k_tl * (R1_fine + R2_fine) / k_deg

#Plotting Null Surface
null_surface = ax.plot_surface(R1_fine, R2_fine, p_null, alpha=0.5, cmap='viridis')

fig.colorbar(null_surface, ax=ax, label='Protein Concentration')

plt.show()
