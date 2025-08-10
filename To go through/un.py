# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 20:25:08 2025

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

hill_parameters = {"k":1.0, "n":4, "K":20, "kleak":.01}
complex_parameters = {'kb':1.0, 'ku':0.01}
component_parameters = {
    #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):5.0, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):2.0,  
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'kb'):1, 
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'ku'):100, 
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = "ktx"): 1, 
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_leak', name = "ku"): 0.001,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_leak', name = "ktx"): 1,
    
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_leak', name = "ku"): 0.001,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_leak', name = "ktx"): 1
    
}


#Species

protease = Species('protease')

IPTG = Species('IPTG', attributes=['input']) #Input A
LacI = Species('LacI', attributes=['repressor'])
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

aTc = Species('aTc', attributes=['input']) #Input B
TetR = Species('TetR', attributes=['repressor'])
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR], parameters = complex_parameters)

rbs = RBS('UTR1')
CDS_SrpR = CDS('CDS_SrpR', 'SrpR')
t16 = Terminator('t16')

CDS_SrpR.protein = Species('SrpR', attributes=['degtagged'])

#Promoters

P_Tac = RepressiblePromoter('P_Tac', repressor = LacI, leak= False, 
                            parameters = hill_parameters)
P_Tet = RepressiblePromoter('P_Tet', repressor = TetR, leak= False,
                            parameters = hill_parameters)

#DNA_constructs

parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2, "kint":.05, 'kdil':0.0075}

SrpR_assembly1 = DNAassembly("SrpR_assembly", promoter=P_Tac, rbs="medium", protein='SrpR')



# dilution_mechanism = Dilution(filter_dict = {"input":False, "machinery":False}, default_on = True)

# global_mechanisms = {"dilution":dilution_mechanism}

degredation_mechanism = Deg_Tagged_Degredation(protease)

global_mechanisms = {"degredation":degredation_mechanism}

M = SimpleTxTlExtract(name="txtl", parameters = parameters, global_mechanisms=global_mechanisms, 
                components=[SrpR_construct, IPTG_LacI, aTc_TetR])

CRN = M.compile_crn()
print(CRN.pretty_print(show_rates = True, show_keys = True))

num_val = 5
max_conc = 10
x0 = {SrpR_construct.get_species():5, LacI:10, TetR:10, protease:7}
timepoints = np.linspace(0, 1000, 1000)

sim = GCSim(CRN)

sim.heatmap(x0, timepoints, max_conc, num_val, IPTG, aTc,'SrpR_degtagged', title = 'SrpR Gate Output', 
            xlabel = 'aTc', ylabel = 'IPTG')

sim.inputswitch(x0, 5000, 'SrpR_degtagged', IPTG, aTc, 10, LacI, TetR, IPTG_LacI, aTc_TetR)

for a in [0, 100]:
    for b in [0, 100]:
        x0 = {SrpR_construct.get_species():5, LacI:100, TetR:100, IPTG:a, aTc:b, protease:7}
        timepoints = np.linspace(0, 5000, 1000)
        R = sim.basicsim(x0, timepoints, ['SrpR_degtagged'], title = f'IPTG = {a}, aTc = {b}')
