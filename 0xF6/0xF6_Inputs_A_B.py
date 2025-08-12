# -*- coding: utf-8 -*-
"""
Created on Sat Aug  9 18:29:44 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import pylab as plt
import pandas as pd
from GCSim import GCSim

#parameters that lead to steady state
# parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2}
# complex_parameters = {'kb':1.0, 'ku':0.01}

parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.05, "kdeg":0.001, "kdil":0.0075}
complex_parameters = {'kb':100, 'ku':10}
component_parameters = {
    #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):10, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):2, 
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'kb'):10, 
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = 'ku'):100, 
    ParameterKey(mechanism = 'transcription_mm', part_id = None, name = "ktx"): 0.05,
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_LacI', name = 'kb'):1, 
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_LacI', name = 'ku'):100, 
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_LacI', name = "ktx"): 0.05,
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_TetR', name = 'kb'):1, 
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_TetR', name = 'ku'):100, 
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_TetR', name = "ktx"): 0.05,
    
    #AraAraC Bound Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_BAD_Ara_2x_AraC_2x', name = 'kb'):100, 
    ParameterKey(mechanism = 'transcription_mm', part_id = "P_BAD_Ara_2x_AraC_2x", name = 'ku'):10, 
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_BAD_Ara_2x_AraC_2x', name = "ktx"): 0.05,
    
    #P_BAD Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_BAD_leak', name = "kb"): 1,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_BAD_leak', name = "ku"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_BAD_leak', name = "ktx"): 0.05,
    
    #P_Tac Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_leak', name = "kb"):100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tac_leak', name = "ktx"): 0.05, 

    #P_Tet Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_Tet_leak', name = "ktx"): 0.05,
    
    #P_SrpR Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_SrpR_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_SrpR_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_SrpR_leak', name = "ktx"): 0.05,

    #P_BetI Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_BetI_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_BetI_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_BetI_leak', name = "ktx"): 0.05,

    #P_HlyIIR Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_HlyIIR_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_HlyIIR_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_HlyIIR_leak', name = "ktx"): 0.05,
    
    #P_AmeR Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_AmeR_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_AmeR_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_AmeR_leak', name = "ktx"): 0.05, 
    
    #P_PhlF Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_PhlF_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_PhlF_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_PhlF_leak', name = "ktx"): 0.05,
    
    #AmtR Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_AmtR_leak', name = "kb"): 100,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_AmtR_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription_mm', part_id = 'P_AmtR_leak', name = "ktx"):0.05 
    
}

#Species

IPTG = Species('IPTG') #Input A
LacI = Species('LacI')
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

aTc = Species('aTc') #Input B
TetR = Species('TetR')
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR], parameters = complex_parameters)

protease = Species('protease')

#DNA parts

rbs = RBS('UTR1')

CDS_HlyIIR= CDS('HlyIIR', 'HlyIIR')
CDS_HlyIIR.protein = Species('HlyIIR', material_type='protein', attributes=['degtagged'])

CDS_SrpR = CDS('SrpR', 'SrpR')
CDS_SrpR.protein = Species('SrpR', material_type='protein', attributes=['degtagged'])

CDS_AmeR = CDS('AmeR', 'AmeR')
CDS_AmeR.protein = Species('AmeR', material_type='protein', attributes=['degtagged'])

CDS_BetI = CDS('BetI', 'BetI')
CDS_BetI.protein = Species('BetI', material_type='protein', attributes=['degtagged'])

CDS_PhlF = CDS('PhlF', 'PhlF')
CDS_PhlF.protein = Species('PhlF', material_type='protein', attributes=['degtagged'])

CDS_AmtR = CDS('AmtR', 'AmtR')
CDS_AmtR.protein = Species('AmtR', material_type='protein', attributes=['degtagged'])

CDS_YFP = CDS('YFP', 'YFP')
CDS_YFP.protein = Species('YFP', material_type='protein', attributes=['degtagged'])

t16 = Terminator('t16')

P_Tac = RegulatedPromoter('P_Tac',  regulators = [LacI], leak=True, 
                          parameters = component_parameters)
P_Tet = RegulatedPromoter('P_Tet', regulators = [TetR], leak=True, 
                          parameters = component_parameters)

P_SrpR = RegulatedPromoter('P_SrpR',  regulators = ['SrpR_degtagged'], leak=True, 
                          parameters = component_parameters)
P_BetI = RegulatedPromoter('P_BetI',  regulators = ['BetI_degtagged'], leak=True, 
                          parameters = component_parameters)
P_HlyIIR = RegulatedPromoter('P_HlyIIR',  regulators = ['HlyIIR_degtagged'], leak=True, 
                          parameters = component_parameters)
P_AmeR = RegulatedPromoter('P_AmeR',  regulators = ['AmeR_degtagged'], leak=True, 
                          parameters = component_parameters)
P_PhlF = RegulatedPromoter('P_PhlF',  regulators = ['PhlF_degtagged'], leak=True, 
                          parameters = component_parameters)

#DNA_constructs

mechanisms = {"transcription":Transcription_MM(Species("RNAP",material_type="protein")), 
              "translation":Translation_MM(Species("Ribo",material_type="protein")), 
              "binding":One_Step_Cooperative_Binding()}

PhlF_construct = DNA_construct([P_SrpR, P_BetI, rbs, CDS_PhlF, t16], mechanisms = mechanisms)

SrpR_construct = DNA_construct([P_Tac, P_Tet, rbs, CDS_SrpR, t16], mechanisms = mechanisms)

BetI_construct = DNA_construct([P_HlyIIR, P_AmeR, rbs, CDS_BetI, t16], mechanisms = mechanisms)

AmeR_construct = DNA_construct([P_Tet, rbs, CDS_AmeR, t16], mechanisms = mechanisms)

HlyIIR_construct = DNA_construct([P_Tac, rbs, CDS_HlyIIR, t16], mechanisms = mechanisms)

YFP_construct = DNA_construct([P_PhlF, rbs, CDS_YFP, t16], mechanisms = mechanisms)

#Mixture and CRN creation

dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)

global_mechanisms = {"dilution":dilution_mechanism}

# degredation_mechanism = Deg_Tagged_Degredation(protease)

# global_mechanisms = {"degredation":degredation_mechanism}

M = TxTlExtract(name="txtl", parameters = parameters, global_mechanisms = global_mechanisms,
                      components=[PhlF_construct, SrpR_construct, BetI_construct, AmeR_construct, 
                                  HlyIIR_construct, YFP_construct, IPTG_LacI, aTc_TetR])
CRN = M.compile_crn()
# CRN.write_sbml_file('Circuit_0xF6_AB_only_sbml.xml') #saving CRN as sbml

print('CRN compiled')

with open('temp_CRN_EQNs.txt', 'w') as f:
    f.write(CRN.pretty_print(show_rates = True, show_keys = True))

sim = GCSim(CRN)

protein_lst = ['protein_SrpR_degtagged', 'protein_BetI_degtagged']

#Plotting
for a in [0, 100]:
    for b in [0, 100]:

        x0 = {PhlF_construct.get_species():5, SrpR_construct.get_species():5, 
              BetI_construct.get_species():5, AmeR_construct.get_species():5, 
              HlyIIR_construct.get_species():5, YFP_construct.get_species():5, 
              IPTG:a, LacI:100, aTc:b, TetR:100, "protein_RNAP":15, 
              "protein_Ribo":150., 'protein_RNAase':45}
        timepoints = np.linspace(0, 4000, 2000)
        R = sim.basicsim(x0, timepoints, protein_lst, title = f'IPTG = {a}, aTc = {b}')
        # print(f"{R['protein_YFP_degtagged'][len(timepoints)-1]}")
        # R.to_excel(f'simulation_results_IPTG_{a}_aTc_{b}.xlsx', index=False)
