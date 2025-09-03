# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 18:54:54 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from GCSim import GCSim

''' degredation and production for all componenets except IPTG and aTc'''

parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.05, "kdeg":0.001, "kdil":0.0075}
complex_parameters = {'kb':100, 'ku':10}
component_parameters = {
    #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):10, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):2.0,  
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'kb'):10, 
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'ku'):100, 
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = "ktx"): 0.05, 
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_leak', name = "ktx"): 0.05,
    
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_leak', name = "ktx"): 0.05
    
}


#Species

protease = Species('protease')

IPTG = Species('IPTG',  material_type='protein', attributes=['input']) #Input A
LacI = Species('LacI',  material_type='protein', attributes=['input'])
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

aTc = Species('aTc',  material_type='protein', attributes=['input']) #Input B
TetR = Species('TetR',  material_type='protein', attributes=['input'])
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR], parameters = complex_parameters)

rbs = RBS('UTR1')
CDS_SrpR = CDS('CDS_SrpR', 'SrpR')
t16 = Terminator('t16')

CDS_SrpR.protein = Species('SrpR', material_type='protein', attributes=['degtagged'])

#Promoters

P_Tac = RegulatedPromoter('P_Tac',  regulators = [LacI], leak=True, 
                          parameters = component_parameters)
P_Tet = RegulatedPromoter('P_Tet', regulators = [TetR], leak=True, 
                          parameters = component_parameters)

#DNA_constructs


SrpR_construct = DNA_construct([P_Tac, P_Tet, rbs, CDS_SrpR, t16], name='SrpR_construct')

dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)

global_mechanisms = {"dilution":dilution_mechanism}

# degredation_mechanism = Deg_Tagged_Degredation(protease)

# global_mechanisms = {"degredation":degredation_mechanism}

M = TxTlExtract(name="txtl", parameters = parameters, global_mechanisms=global_mechanisms, 
                components=[SrpR_construct, IPTG_LacI, aTc_TetR])

CRN = M.compile_crn()

with open('temp_CRN_EQNs.txt', 'w') as f:
    f.write(CRN.pretty_print(show_rates = True, show_keys = True))

# num_val = 5
# max_conc = 10
# x0 = {"protein_RNAP":0.5, "protein_Ribo":50., 'protein_RNAase':15, 
#       SrpR_construct.get_species():1, LacI:10, TetR:10}
# timepoints = np.linspace(0, 1000, 1000)

# sim = GCSim(CRN)

# sim.heatmap(x0, timepoints, max_conc, num_val, IPTG, aTc,'protein_SrpR_degtagged', title = 'SrpR Gate Output', 
#             xlabel = 'aTc', ylabel = 'IPTG')

# sim.inputswitch(x0, 5000, 'protein_SrpR_degtagged', IPTG, aTc, 10, LacI, TetR, IPTG_LacI, aTc_TetR)

# for a in [0, 10]:
#     for b in [0, 10]:
#         x0 = {SrpR_construct.get_species():1, LacI:10, TetR:10, IPTG:a, aTc:b,
#               "protein_RNAP":0.5, "protein_Ribo":50., 'protein_RNAase':15}
#         timepoints = np.linspace(0, 21000, 1000)
#         R = sim.basicsim(x0, timepoints, ['protein_SrpR_degtagged'], title = f'IPTG = {a}, aTc = {b}')
        
#         R.to_excel(f'simulation_results_IPTG_{a}_aTc_{b}.xlsx', index=False)