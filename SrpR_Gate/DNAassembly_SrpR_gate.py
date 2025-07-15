# -*- coding: utf-8 -*-
"""
Created on Sun Jul 13 17:53:46 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import pylab as plt
import pandas as pd

#Parameters #!!! parameters mostly default currently #!!! add in degradation???

# hill_parameters = {"k":1.0, "n":4, "K":20, "kleak":.01}
complex_parameters = {'kb':0.5, 'ku':1.0}
component_parameters = {
    #Promoter LacI Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'binding', part_id = 'P_Tac_Tet_LacI', name = 'kb'):100, 
    ParameterKey(mechanism = 'binding', part_id = "P_Tac_Tet_LacI", name = 'ku'):5.0, 
    ParameterKey(mechanism = 'binding', part_id = "P_Tac_Tet_LacI", name = 'cooperativity'):4.0, 
    
    #LacI Bound Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_LacI', name = 'kb'):1, 
    ParameterKey(mechanism = 'transcription', part_id = "P_Tac_Tet_LacI", name = 'ku'):100, 
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_LacI', name = "ktx"): 1, 
    
    #Promoter TetR Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'binding', part_id = 'P_Tac_Tet_TetR', name = 'kb'):100,
    ParameterKey(mechanism = 'binding', part_id = "P_Tac_Tet_TetR", name = 'ku'):5.0,
    ParameterKey(mechanism = 'binding', part_id = "P_Tac_Tet_TetR", name = 'cooperativity'):4.0,
    
    #TetR Bound Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_TetR', name = 'kb'):1,
    ParameterKey(mechanism = 'transcription', part_id = "P_Tac_Tet_TetR", name = 'ku'):100,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_TetR', name = "ktx"): 1, 
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_leak', name = "kb"): 2,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_leak', name = "ktx"): 1.0, 
}

#Species

IPTG = Species('IPTG') #Input A
LacI = Species('LacI')
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

aTc = Species('aTc') #Input B
TetR = Species('TetR')
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR], parameters = complex_parameters)

SrpR = Species('SrpR') #Gate Protein

#Promoters

P_Tac_Tet = RegulatedPromoter('P_Tac_Tet', regulators = [LacI,TetR], leak=True, 
                              parameters = component_parameters)

#DNAassemblies

SrpR_assembly = DNAassembly(name='SrpR_assembly', promoter=P_Tac_Tet, rbs = 'strong', 
                             protein = SrpR)


M = TxTlExtract(name="e coli extract", parameter_file = 'default_parameters.txt', 
                      components=[SrpR_assembly, IPTG_LacI, aTc_TetR])
CRN = M.compile_crn()

CRN.write_sbml_file('Circuit_0xF6_bioCRNpyler_sbml.xml') #saving CRN as sbml

print(CRN.pretty_print(show_rates = True, show_keys = True))

plt.figure(figsize = (6, 6))
N = 11 #Number of titrations
max_titration = 10
HM = np.zeros((N, N))
for a_ind, A_c in enumerate(np.linspace(0, max_titration, N)):
    for b_ind, B_c in enumerate(np.linspace(0, max_titration, N)):
        x0 = {SrpR_assembly.dna:1, IPTG:A_c, aTc:B_c, LacI:10, TetR:10}
        timepoints = np.linspace(0, 1000, 1000)
        R1 = CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0)
        HM[a_ind, b_ind] = R1['SrpR'][len(timepoints)-1]
plt.title("")
cb = plt.pcolor(HM)
plt.colorbar(cb)
plt.xlabel("aTc")
plt.ylabel("IPTG")

plt.xticks(np.arange(.5, N+.5, 1), [str(i) for i in np.linspace(0, max_titration, N)])
plt.yticks(np.arange(.5, N+.5, 1), [str(i) for i in np.linspace(0, max_titration, N)])
    
plt.show()