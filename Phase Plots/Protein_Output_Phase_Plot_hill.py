# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 15:18:02 2025

@author: zacha
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
ko = 0.05 #open complex production rate
ko_f = 0.033 # forward RNAP binding rate
ko_r = 1 #reverse RNAP binding rate
nr = 30 #Initial RNAP count
kr_f = 0.5 #Forward repression binding rate
kr_r = 1.0 #reverse repression binding rate
nc = 2 #stochiometry of binding
pAmeR = 2 #promoter
pHlyIIR = 2 #promoter
kd = 0.0075 #protein degradation rate

# Protein equation
def f(AmeR, HlyIIR, BetI):
    
    BetI_1 = pAmeR*ko* ko_f/ko_r *nr / (1+ko_f/ko_r*nr+(kr_f/kr_r*AmeR)**nc)
    BetI_2 = pHlyIIR*ko*ko_f/ko_r*nr/(1+ko_f/ko_r*nr+(kr_f/kr_r*HlyIIR)**nc)
    
    return 10*BetI_1 + 10*BetI_2 - kd*BetI

# Create a grid in the (rna1, rna2) plane
AmeR_values = np.linspace(0, 65, 5)
HlyIIR_values = np.linspace(0, 65, 5)
BetI_values = np.linspace(0, 120, 4)
A, H, B = np.meshgrid(AmeR_values, HlyIIR_values, BetI_values)

# Compute the derivative on the grid
dB = f(A, H, B)

# For quiver
U = np.zeros_like(A)
V = np.zeros_like(H) 
W = dB

# 3D quiver plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(A, H, B, U, V, W, length = 20)

ax.set_xlabel('AmeR')
ax.set_ylabel('HlyIIR')
ax.set_zlabel('BetI')

# ax.set_xlim(0, 4)
# ax.set_ylim(0, 4)
# ax.set_zlim(10, 30)

# Nullcline surface
AmeR_fine = np.linspace(0, 65, 100)
HlyIIR_fine = np.linspace(0, 65, 100)
A_fine, H_fine = np.meshgrid(AmeR_fine, HlyIIR_fine)

BetI_1 = pAmeR*ko* ko_f/ko_r *nr / (1+ko_f/ko_r*nr+(kr_f/kr_r*A_fine)**nc)
BetI_2 = pHlyIIR*ko*ko_f/ko_r*nr/(1+ko_f/ko_r*nr+(kr_f/kr_r*H_fine)**nc)

B_null = (10*BetI_1 + 10*BetI_2) / kd

#Plotting Null Surface
null_surface = ax.plot_surface(A_fine, H_fine, B_null, alpha=0.5, cmap='viridis')

fig.colorbar(null_surface, ax=ax, label='Protein Concentration')

plt.show()
