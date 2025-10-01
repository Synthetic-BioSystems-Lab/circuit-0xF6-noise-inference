# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 14:49:20 2025

@author: zacha
"""

import matplotlib.pyplot as plt
import pandas as pd
from EQN_fit import EQN_fit
from scipy.integrate import solve_ivp

protein_lst = ['AmeR-F1']

exp_values_dict = {'AmeR-F1':[0.2, 3.8, 0.09, 1.4], 'AmtR-A1':[0.06, 3.8, 0.07, 1.6], 
                   'BetI-E1':[0.07, 3.8, 0.41, 2.4], 'HlyIIR-H1':[0.07, 2.5, 0.19, 2.6], 
                   'PhlF-P1':[0.01, 3.9, 0.03, 4.0], 'PhlF-P2':[0.02, 4.1, 0.13, 3.9], 
                   'PhlF-P3':[0.02, 6.8, 0.23, 4.2], 'SrpR-S1':[0.003, 1.3, 0.01, 2.9], 
                   'SrpR-S2':[0.003, 2.1, 0.04, 2.6], 'SrpR-S3':[0.004, 2.1, 0.06, 2.8], 
                   'SrpR-S4':[0.007, 2.1, 0.1, 2.8]}

fits = {}

for i in range(len(protein_lst)):
    
    # Experimantal Values
    y_min = exp_values_dict[protein_lst[i]][0]
    y_max = exp_values_dict[protein_lst[i]][1]
    K = exp_values_dict[protein_lst[i]][2]
    n = exp_values_dict[protein_lst[i]][3]
    
    fits[f'{protein_lst[i]}_fit'] = EQN_fit(True, y_min, y_max, K, n)
         
    kb, ku, coop, ktx_leak, ktx_free, r_squared = fits[f'{protein_lst[i]}_fit'].data_fit_hill(protein_lst[i])
    
    print(f'{protein_lst[i]}: kb = {kb}, ku = {ku}, Cooperativity = {coop}, ktx_leak = {ktx_leak}, '
          f'ktx_free = {ktx_free}')
    
    R = 0 # initial Repressor concentration
    D = 0.1
    B = 0 # bound DNA
    M = 0 # mRNA
    A = 0 # protein A
    
    y0 = [R, D, B, M, A]
    t_span = (0, 1500)
    
    sol = solve_ivp(EQN_fit.crn_odes, t_span, y0, method='LSODA', 
                    args=(kb, ku, coop, ktx_leak, ktx_free))
    
    R = 10 # initial Repressor concentration
    
    y0 = [R, D, B, M, A]
    t_span = (0, 1500)
    
    sol2 = solve_ivp(EQN_fit.crn_odes, t_span, y0, method='LSODA', 
                    args=(kb, ku, coop, ktx_leak, ktx_free))
    
    plt.figure()
    plt.title(f'{protein_lst[i]}')
    plt.plot(sol.t, sol.y[4:].T, label='R = 0')
    plt.plot(sol2.t, sol2.y[4:].T, label='R = 10')
    plt.xlabel('Time')
    plt.ylabel(f'[{protein_lst[i]}]')
    plt.legend()
    plt.show()
    
    df = pd.DataFrame([{
        'Protein':protein_lst[i],
        "kb": kb,
        "ku": ku,
        "Cooperativity": coop,
        "ktx_leak": ktx_leak,
        "ktx_free": ktx_free,
        "R²": r_squared
    }])

    EQN_fit.save_to_excel(protein_lst[i], df, 'Outputs/Protein Parameters.xlsx', i + 2 )
