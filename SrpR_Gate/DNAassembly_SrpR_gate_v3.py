# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 19:27:56 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from GCSim import GCSim

#Parameters #!!! parameters mostly default currently #!!! add in degradation???

component_parameters = {
    #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):10.0, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):2.0,  
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'kb'):1, 
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'ku'):100, 
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = "ktx"): 0.05,  
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_leak', name = "ktx"): 0.05, 
}

#Species

IPTG = Species('IPTG') #Input A
LacI = Species('LacI')
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI])

aTc = Species('aTc') #Input B
TetR = Species('TetR')
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR])

protease = Species('protease')

SrpR = Species('SrpR', attributes=['degtagged'])

#Promoters

P_Tac_Tet = RegulatedPromoter('P_Tac_Tet', regulators = [LacI,TetR], leak=True, 
                              parameters = component_parameters)

#DNAassemblies

SrpR_assembly = DNAassembly(name='SrpR_assembly', promoter=P_Tac_Tet, rbs = 'strong', 
                             protein = SrpR)


M = TxTlExtract(name="e coli extract", parameter_file = 'default_parameters.txt', 
                      components=[SrpR_assembly, IPTG_LacI, aTc_TetR])

degredation_mechanism = Deg_Tagged_Degredation(protease)

rna_deg = Degredation_mRNA_MM(Species("RNAase",material_type="protein"))

dilution_mechanism = Dilution()

M.add_mechanism(degredation_mechanism, mech_type='degredation')
M.add_mechanism(rna_deg, mech_type= 'rna_degredation')
# M.add_mechanism(dilution_mechanism, mech_type='dilution')

CRN = M.compile_crn()
print(CRN.pretty_print(show_rates = True, show_keys = True))

num_val = 11
max_conc = 10
x0 = {SrpR_assembly.dna:1, LacI:10, TetR:10, protease:10}
timepoints = np.linspace(0, 10000, 100000)

sim = GCSim(CRN)

sim.heatmap(x0, timepoints, max_conc, num_val, IPTG, aTc,'SrpR_degtagged', title = 'SrpR Gate Output', 
            xlabel = 'aTc', ylabel = 'IPTG')

# sim.inputswitch(x0, 10000, 'SrpR_degtagged', IPTG, aTc, 10, LacI, TetR, IPTG_LacI, aTc_TetR)

# for a in [0, 10]:
#     for b in [0, 10]:
#         x0 = {SrpR_assembly.dna:1, LacI:10, TetR:10, IPTG:a, aTc:b, protease:10}
#         timepoints = np.linspace(0, 10000, 100000)
#         R = sim.basicsim(x0, timepoints, ['SrpR_degtagged'], title = f'IPTG = {a}, aTc = {b}')
#         print(f"{R['rna_SrpR_assembly'][len(timepoints)-1]}")