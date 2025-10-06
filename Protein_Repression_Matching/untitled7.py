# -*- coding: utf-8 -*-
"""
Created on Sat Oct  4 18:24:08 2025

@author: zacha
"""

import matplotlib.pyplot as plt
import pandas as pd
from utils import EQN_fit
from scipy.integrate import solve_ivp
from biocrnpyler import *
import numpy as np

x = [0, 5, 10, 20, 30, 40, 50, 70, 100, 150, 200, 1000]
y_exp = [0.004847121, 0.00742964, 0.012730302, 0.034061053, 0.062861841, 0.100932785, 
         0.144974067, 0.24383541, 0.417799333, 0.742963951, 1, 2.101748011]

kb = 100
ku = 10
coop = 2
ktx_leak = 0.00005
ktx_free = 0.0225
kb_C = 1
ku_C = 1000
coop_I = 2

y_sim = []

for I in x:

    C = 0 
    R = 10
    D = 0.1 #Initial DNA conc
    B = 0 # bound DNA
    M = 0 # mRNA
    A = 0 # protein A
    
    y0 = [I, C, R, D, B, M, A]
    t_span = (0, 1500)
    
    sol = solve_ivp(EQN_fit.crn_odes_PTAC_YFP, t_span, y0, method='LSODA', 
                    args=(kb_C, ku_C, coop_I, kb, ku, coop, ktx_leak, ktx_free))
    y_sim.append(sol.y[6,-1])

ss_res = np.sum((np.array(y_exp) - np.array(y_sim))**2)
ss_tot = np.sum((np.array(y_exp) - np.mean(y_exp))**2)
r_squared = 1 - (ss_res/ss_tot)

plt.figure()
plt.loglog(x, y_exp, marker='o', linestyle='', label = 'Experimental')
plt.loglog(x, y_sim, marker='o', linestyle='', label='Simulation')
plt.title(f'R² = {r_squared:.4f}')
plt.xlabel('Time')
plt.ylabel(f'')
plt.legend()
plt.show()
