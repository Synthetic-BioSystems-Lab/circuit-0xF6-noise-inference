# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 15:56:31 2025

@author: zacha
"""

import matplotlib.pyplot as plt
import pandas as pd
from utils import EQN_fit
from scipy.integrate import solve_ivp
from biocrnpyler import *
import numpy as np

x_exp = [0, 5, 10, 20, 30, 40, 50, 70, 100, 150, 200, 1000]
y_exp = [0.004847121, 0.00742964, 0.012730302, 0.034061053, 0.062861841, 0.100932785, 
         0.144974067, 0.24383541, 0.417799333, 0.742963951, 1, 2.101748011]

fit = EQN_fit(False, x=x_exp, y=y_exp)

params = fit.data_fit()

kb, ku, coop, ktx_leak, ktx_free, kb_C, ku_C, coop_I, R, r_squared = params

print(f'kb = {kb}, ku = {ku}, Cooperativity = {coop}, ktx_leak = {ktx_leak}, '
      f'ktx_free = {ktx_free}, kb_C = {kb_C}, ku_C = {ku_C}, Cooperativity I = {coop_I}, '
      f'R = {R}')

I = 10
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

# Parameters and Global Mechanisms
parameters={"cooperativity":coop_I,"kb":100, "ku":10, "ktx":.05, "ktl":0.05, "kdeg":0.001, "kdil":0.0075}
complex_parameters = {'kb':kb_C, 'ku':ku_C}
component_parameters = {
    
    #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):kb, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):ku, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):coop, 
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'simple_transcription', part_id = None, name = "ktx"): ktx_leak,
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Tac_leak', name = "ktx"): ktx_free
    
}

dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)

global_mechanisms = {"dilution":dilution_mechanism}

# Species

IPTG = Species('IPTG',  material_type='protein', attributes=['input']) 
LacI = Species('LacI',  material_type='protein', attributes=['input'])
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

# DNA parts

P_Tac = RegulatedPromoter('P_Tac',  regulators = [LacI], leak=True, 
                          parameters = component_parameters)

rbs = RBS('rbs')

CDS_YFP = CDS('YFP', 'YFP')
CDS_YFP.protein = Species('YFP', material_type='protein', attributes=['degtagged'])

t = Terminator('t')

YFP_construct = DNA_construct([P_Tac, rbs, CDS_YFP, t])

# Mixture and CRN creation

M = SimpleTxTlExtract('simtxtl', parameters = parameters, global_mechanisms=global_mechanisms, 
                      components=[YFP_construct, IPTG_LacI])

CRN = M.compile_crn()

# Simulation and Plotting  
x0 = {YFP_construct.get_species():D, LacI:R, IPTG:I}
timepoints = np.linspace(0, 1500, 500)
Res = CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0)

plt.figure()
plt.plot(Res['time'], Res['protein_YFP_degtagged'], label = 'YFP')
plt.plot(sol.t, sol.y[6:].T, label='Protein A')
plt.xlabel('Time')
plt.ylabel(f'')
plt.legend()
plt.show()
