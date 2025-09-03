# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 17:14:57 2025

@author: zacha
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
k_tx = 0.05
k_tl = 0.05
k_deg = 0.0075

Vmax = 0.5
n = 2 
K = 2  

# Protein equation
def f(r1, r2, p):
    
    # uses semi-approx behavior hill equation for dna repressor binding behavoir validated in untitled5 in folder 0xF6 
    inactive_dna = np.where(
    (r1 >= 1) & (r2 >= 1),  # condition
    ( (Vmax*((r1-1)**n))/((K**n)+((r1-1)**n)) ) + 
    ( (Vmax*((r2-1)**n))/((K**n)+((r2-1)**n)) ),
    0) #zero is else term
        
    active_dna = 1 - inactive_dna
    rna = k_tx*(active_dna)
    
    return k_tl*(rna) - k_deg*p

# Create a grid in the (rna1, rna2) plane
r1_values = np.linspace(0, 65, 5)
r2_values = np.linspace(0, 65, 5)
p_values = np.linspace(0, 0.4, 4)
R1, R2, P = np.meshgrid(r1_values, r2_values, p_values)

# Compute the derivative on the grid
dp = f(R1, R2, P)

# For quiver, we need vectors. We'll show change along R1 and R2 as zeros (since only p changes)
U = np.zeros_like(R1)  # change along rna1 (no change)
V = np.zeros_like(R2)  # change along rna2 (no change)
W = dp                # change along p

# 3D quiver plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(R1, R2, P, U, V, W, length = 20)

ax.set_xlabel('r1')
ax.set_ylabel('r2')
ax.set_zlabel('Protein (p)')

# ax.set_xlim(0, 4)
# ax.set_ylim(0, 4)
# ax.set_zlim(10, 30)

# Nullcline surface
r1_fine = np.linspace(1, 65, 100)
r2_fine = np.linspace(1, 65, 100)
R1_fine, R2_fine = np.meshgrid(r1_fine, r2_fine)

inactive_dna = ( (Vmax*((R1_fine-1)**n))/((K**n)+((R1_fine-1)**n)) ) + \
( (Vmax*((R2_fine-1)**n))/((K**n)+((R2_fine-1)**n)) )

active_dna = 1 - inactive_dna
rna = k_tx*(active_dna)

p_null = k_tl * (rna) / k_deg

#Plotting Null Surface
null_surface = ax.plot_surface(R1_fine, R2_fine, p_null, alpha=0.5, cmap='viridis')

fig.colorbar(null_surface, ax=ax, label='Protein Concentration')

plt.show()
