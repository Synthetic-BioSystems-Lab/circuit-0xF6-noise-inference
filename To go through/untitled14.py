# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 15:59:21 2025

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

complex_parameters = {'kb':1.0, 'ku':0.01}
component_parameters = {
    #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):5.0, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):2.0, 
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription', part_id = None, name = 'kb'):0.001, 
    ParameterKey(mechanism = 'transcription', part_id = None, name = 'ku'):100, 
    ParameterKey(mechanism = 'transcription', part_id = None, name = "ktx"): 1, 
    
    #P_HlyIIR Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription', part_id = 'P_HlyIIR_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription', part_id = 'P_HlyIIR_leak', name = "ku"): 1,
    ParameterKey(mechanism = 'transcription', part_id = 'P_HlyIIR_leak', name = "ktx"): 1,
    
    #P_AmeR Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription', part_id = 'P_AmeR_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription', part_id = 'P_AmeR_leak', name = "ku"): 1,
    ParameterKey(mechanism = 'transcription', part_id = 'P_AmeR_leak', name = "ktx"): 1,
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_leak', name = "ku"): 0.001,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_leak', name = "ktx"): 1,
    
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tet_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tet_leak', name = "ku"): 0.001,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tet_leak', name = "ktx"): 1
    
}


#Species

protease = Species('protease')

aTc = Species('aTc', attributes=['input']) #Input B
TetR = Species('TetR', attributes=['repressor'])
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR], parameters = complex_parameters)

rbs = RBS('UTR1')

CDS_AmeR = CDS('AmeR', 'AmeR')
CDS_AmeR.protein = Species('AmeR', attributes=['degtagged'])

CDS_BetI = CDS('BetI', 'BetI')
CDS_BetI.protein = Species('BetI', attributes=['degtagged'])

t16 = Terminator('t16')

#Promoters

P_Tet = RegulatedPromoter('P_Tet', regulators = [TetR], leak=True, 
                          parameters = component_parameters)
P_AmeR = RegulatedPromoter('P_AmeR',  regulators = ['AmeR_degtagged'], leak=True, 
                          parameters = component_parameters)

#DNA_constructs

parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2, "kint":.05, 'kdil':0.0075}
mechanisms = {"transcription":Transcription_MM(Species("RNAP",material_type="protein", attributes=['machinery'])), 
              "translation":Translation_MM(Species("Ribo",material_type="protein", attributes=['machinery'])), 
              "binding":One_Step_Cooperative_Binding()}

AmeR_construct = DNA_construct([P_Tet, rbs, CDS_AmeR, t16], mechanisms = mechanisms)


BetI_construct = DNA_construct([P_AmeR, rbs, CDS_BetI, t16], mechanisms = mechanisms)

# dilution_mechanism = Dilution(filter_dict = {"input":False, "machinery":False}, default_on = True)

# global_mechanisms = {"dilution":dilution_mechanism}

degredation_mechanism = Deg_Tagged_Degredation(protease)

global_mechanisms = {"degredation":degredation_mechanism}

M = TxTlExtract(name="txtl", parameters = parameters, global_mechanisms=global_mechanisms, 
                components=[BetI_construct, AmeR_construct, aTc_TetR])

CRN = M.compile_crn()
print(CRN.pretty_print(show_rates = True, show_keys = True))

sim = GCSim(CRN)

for b in [0, 100]:
    x0 = {BetI_construct.get_species():5, 
          AmeR_construct.get_species():5, TetR:10, aTc:b, protease:1,
          "protein_RNAP_machinery":5, "protein_Ribo_machinery":50., 'protein_RNAase':100}
    timepoints = np.linspace(0, 5000, 1000)
    R = sim.basicsim(x0, timepoints, ['BetI_degtagged', 'AmeR_degtagged'], 
                     title = f'aTc = {b}')