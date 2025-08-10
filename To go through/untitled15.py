# -*- coding: utf-8 -*-
"""
Created on Tue Aug  5 16:17:22 2025

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

P_Tac = RegulatedPromoter('P_Tac',  regulators = [LacI], leak=True, 
                          parameters = component_parameters)
P_Tet = RegulatedPromoter('P_Tet', regulators = [TetR], leak=True, 
                          parameters = component_parameters)
P_HlyIIR = RegulatedPromoter('P_HlyIIR',  regulators = ['HlyIIR_degtagged'], leak=True, 
                          parameters = component_parameters)
P_AmeR = RegulatedPromoter('P_AmeR',  regulators = ['AmeR_degtagged'], leak=True, 
                          parameters = component_parameters)

#DNA_constructs

parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2, "kint":.05, 'kdil':0.0075}
mechanisms = {"transcription":Transcription_MM(Species("RNAP",material_type="protein", attributes=['machinery'])), 
              "translation":Translation_MM(Species("Ribo",material_type="protein", attributes=['machinery'])), 
              "binding":One_Step_Cooperative_Binding()}

SrpR_construct = DNA_construct([P_Tac, P_Tet, rbs, CDS_SrpR, t16], mechanisms=mechanisms)

AmeR_construct = DNA_construct([P_Tet, rbs, CDS_AmeR, t16], mechanisms = mechanisms)

HlyIIR_construct = DNA_construct([P_Tac, rbs, CDS_HlyIIR, t16], mechanisms = mechanisms)

BetI_construct = DNA_construct([P_HlyIIR, P_AmeR, rbs, CDS_BetI, t16], mechanisms = mechanisms)

# dilution_mechanism = Dilution(filter_dict = {"input":False, "machinery":False}, default_on = True)

# global_mechanisms = {"dilution":dilution_mechanism}

degredation_mechanism = Deg_Tagged_Degredation(protease)

global_mechanisms = {"degredation":degredation_mechanism}

M = TxTlExtract(name="txtl", parameters = parameters, global_mechanisms=global_mechanisms, 
                components=[SrpR_construct, AmeR_construct, HlyIIR_construct, BetI_construct,
                            IPTG_LacI, aTc_TetR])

CRN = M.compile_crn()
print(CRN.pretty_print(show_rates = True, show_keys = True))

sim = GCSim(CRN)

for a in [0, 100]:
    for b in [0, 100]:
        x0 = {SrpR_construct.get_species():5, AmeR_construct.get_species():5, 
              HlyIIR_construct.get_species():5, BetI_construct.get_species():5,
              LacI:100, TetR:100, IPTG:a, aTc:b, protease:1,
              "protein_RNAP_machinery":5., "protein_Ribo_machinery":50., 'protein_RNAase':100}
        timepoints = np.linspace(0, 5000, 1000)
        R = sim.basicsim(x0, timepoints, 
                         ['SrpR_degtagged', 'AmeR_degtagged', 'HlyIIR_degtagged', 'BetI_degtagged'], 
                         title = f'IPTG = {a}, aTc = {b}')
