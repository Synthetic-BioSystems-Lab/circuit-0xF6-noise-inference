# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 17:52:13 2025

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

#Species

protease = Species('protease')

SrpR = Species('SrpR', attributes=['degtagged'])

IPTG = Species('IPTG', attributes=['input']) #Input A
LacI = Species('LacI', attributes=['repressor'])
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI])

aTc = Species('aTc', attributes=['input']) #Input B
TetR = Species('TetR', attributes=['repressor'])
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR])

#Promoters

P_Tac = RepressiblePromoter('P_Tac', repressor = LacI, leak= False)
P_Tet = RepressiblePromoter('P_Tet', repressor = TetR, leak= False)
#DNA_constructs

# parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":1, "ktl":.2, "kdeg":2, "kint":.05, 'kdil':0.0075}
# mechanisms = {"transcription":Transcription_MM(Species("RNAP",material_type="protein", attributes=['machinery'])), 
#               "translation":Translation_MM(Species("Ribo",material_type="protein", attributes=['machinery'])), 
#               "binding":One_Step_Cooperative_Binding()}

SrpR_assembly1 = DNAassembly("SrpR_assembly1", promoter=P_Tac, rbs="medium", protein=SrpR)

SrpR_assembly2 = DNAassembly("SrpR_assembly2", promoter=P_Tet, rbs="medium", protein=SrpR)

# dilution_mechanism = Dilution(filter_dict = {"input":False, "machinery":False}, default_on = True)

# global_mechanisms = {"dilution":dilution_mechanism}

degredation_mechanism = Deg_Tagged_Degredation(protease)

# global_mechanisms = {"degredation":degredation_mechanism}

M = SimpleTxTlExtract(name="txtl", parameter_file = 'default_parameters.txt', 
                components=[SrpR_assembly1, SrpR_assembly2, IPTG_LacI, aTc_TetR])

M.add_mechanism(degredation_mechanism, mech_type='degredation')

CRN = M.compile_crn()
print(CRN.pretty_print(show_rates = True, show_keys = True))

num_val = 11
max_conc = 10
x0 = {SrpR_assembly1.dna:10, SrpR_assembly2.dna:10, LacI:10, TetR:10, protease:10}
timepoints = np.linspace(0, 6000, 1000)

sim = GCSim(CRN)

sim.heatmap(x0, timepoints, max_conc, num_val, IPTG, aTc,'SrpR_degtagged', title = 'SrpR Gate Output', 
            xlabel = 'aTc', ylabel = 'IPTG')

# sim.inputswitch(x0, 1000, 'BetI_degtagged', IPTG, aTc, 10, LacI, TetR, IPTG_LacI, aTc_TetR)

for a in [0, 10]:
    for b in [0, 10]:
        x0 = {SrpR_assembly1.dna:10, SrpR_assembly2.dna:10, LacI:10, TetR:10, IPTG:a, aTc:b, protease:10}
        timepoints = np.linspace(0, 6000, 1000)
        R = sim.basicsim(x0, timepoints, 
                         ['SrpR_degtagged'], 
                         title = f'IPTG = {a}, aTc = {b}')