# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 19:26:39 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from GCSim import GCSim

#Parameters #!!! parameters mostly default currently #!!! add in degradation for complexes? try overall circuit? 

''' degredation and production for all componenets except IPTG and aTc'''

complex_parameters = {'kb':100.0, 'ku':0.001}
# component_parameters = {
#     #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
#     ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
#     ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):5.0, 
#     ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):2.0,  
    
#     #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
#     ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'kb'):1, 
#     ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'ku'):100, 
#     ParameterKey(mechanism = 'transcription_mm', part_id = None, name = "ktx"): 1, 
    
#     #Leak Parameters for transcription
#     #These regulate expression of an unbound promoter
#     ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_Tet_leak', name = "kb"): 100,
#     ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_Tet_leak', name = "ku"): 0.001,
#     ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_Tet_leak', name = "ktx"): 1
# }


#Species

protease = Species('protease')

IPTG = Species('IPTG', attributes=['input']) #Input A
LacI = Species('LacI', attributes=['repressor'])
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

aTc = Species('aTc', attributes=['input']) #Input B
TetR = Species('TetR', attributes=['repressor'])
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR], parameters = complex_parameters)

rbs = RBS('UTR1')

CDS_HlyIIR= CDS('HlyIIR', 'HlyIIR')
CDS_HlyIIR.protein = Species('HlyIIR', attributes=['degtagged'])

CDS_SrpR = CDS('SrpR', 'SrpR')
CDS_SrpR.protein = Species('SrpR', attributes=['degtagged'])

CDS_AmeR = CDS('AmeR', 'AmeR')
CDS_AmeR.protein = Species('AmeR', attributes=['degtagged'])

CDS_BetI = CDS('BetI', 'BetI')
CDS_BetI.protein = Species('BetI', attributes=['degtagged'])

CDS_PhlF = CDS('PhlF', 'PhlF')
CDS_PhlF.protein = Species('PhlF', attributes=['degtagged'])

CDS_AmtR = CDS('AmtR', 'AmtR')
CDS_AmtR.protein = Species('AmtR', attributes=['degtagged'])

CDS_YFP = CDS('YFP', 'YFP')
CDS_YFP.protein = Species('YFP', attributes=['degtagged'])

t16 = Terminator('t16')
#Promoters

P_Tac_Tet = CombinatorialPromoter('P_Tac_Tet', regulators = [TetR, LacI], 
                                  tx_capable_list= [[LacI], [TetR]], leak = True)

#DNA_constructs

# parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":1, "ktl":.2, "kdeg":2, "kint":.05, 'kdil':0.0075}
# mechanisms = {"transcription":Transcription_MM(Species("RNAP",material_type="protein", attributes=['machinery'])), 
#               "translation":Translation_MM(Species("Ribo",material_type="protein", attributes=['machinery'])), 
#               "binding":One_Step_Cooperative_Binding()}

SrpR_assembly = DNAassembly("SrpR_assembly", promoter=P_Tac_Tet, rbs="medium", protein='SrpR')

# dilution_mechanism = Dilution(filter_dict = {"input":False, "machinery":False}, default_on = True)

# global_mechanisms = {"dilution":dilution_mechanism}

degredation_mechanism = Deg_Tagged_Degredation(protease)

global_mechanisms = {"degredation":degredation_mechanism}

M = SimpleTxTlDilutionMixture(name="txtl", parameter_file = 'default_parameters.txt',  
                components=[SrpR_assembly, IPTG_LacI, aTc_TetR])

CRN = M.compile_crn()
print(CRN.pretty_print(show_rates = True, show_keys = True))

num_val = 5
max_conc = 10
x0 = {SrpR_assembly.dna:1, LacI:10, TetR:10}
timepoints = np.linspace(0, 5000, 1000)

sim = GCSim(CRN)

sim.heatmap(x0, timepoints, max_conc, num_val, IPTG, aTc,'protein_SrpR', title = 'SrpR Gate Output', 
            xlabel = 'aTc', ylabel = 'IPTG')

# sim.inputswitch(x0, 1000, 'BetI_degtagged', IPTG, aTc, 10, LacI, TetR, IPTG_LacI, aTc_TetR)

for a in [0, 100]:
    for b in [0, 100]:
        x0 = {SrpR_assembly.dna:1, LacI:10, TetR:10, IPTG:a, aTc:b}
        timepoints = np.linspace(0, 5000, 1000)
        R = sim.basicsim(x0, timepoints, 
                         ['protein_SrpR'], 
                         title = f'IPTG = {a}, aTc = {b}')
