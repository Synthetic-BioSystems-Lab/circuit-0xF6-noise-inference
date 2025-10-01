# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 15:56:31 2025

@author: zacha
"""

import matplotlib.pyplot as plt
import pandas as pd
from EQN_fit import EQN_fit
from scipy.integrate import solve_ivp

x_exp = [0, 5, 10, 20, 30, 40, 50, 70, 100, 150, 200, 1000]
y_exp = [0.004847121, 0.00742964, 0.012730302, 0.034061053, 0.062861841, 0.100932785, 
         0.144974067, 0.24383541, 0.417799333, 0.742963951, 1, 2.101748011]

fit = EQN_fit(False, x=x_exp, y=y_exp)

params = fit.data_fit()

kb, ku, coop, ktx_leak, ktx_free, kb_C, ku_C, coop_I, R, r_squared = params

I = 10
C = 0 
D = 0.1 #Initial DNA conc
B = 0 # bound DNA
M = 0 # mRNA
A = 0 # protein A

y0 = [I, C, R, D, B, M, A]
t_span = (0, 1500)

sol = solve_ivp(EQN_fit.crn_odes_PTAC_YFP, t_span, y0, method='LSODA', 
                args=(kb_C, ku_C, coop_I, kb, ku, coop, ktx_leak, ktx_free))

plt.figure()
plt.plot(sol.t, sol.y[6:].T, label='')
plt.xlabel('Time')
plt.ylabel(f'')
plt.legend()
plt.show()
