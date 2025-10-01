# -*- coding: utf-8 -*-
"""
Created on Sat Sep 20 14:21:49 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import pylab as plt
import pandas as pd
from GCSim import GCSim

#parameters that lead to steady state
# parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.05, "kdeg":0.0075}
# complex_parameters = {'kb':100, 'ku':10}

parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.05, "kdeg":0.001, "kdil":0.0075}
complex_parameters = {'kb':100, 'ku':10}
component_parameters = {
    #Defalt Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):10, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):4,
    ParameterKey(mechanism = 'simple_transcription', part_id = None, name = "ktx"): 1e-15,
    
    # AmeR Parameters
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_AmeR_AmeR_degtagged', name = 'kb'):100.75, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_AmeR_AmeR_degtagged', name = 'ku'):2.50, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_AmeR_AmeR_degtagged', name = 'cooperativity'):1, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_AmeR_leak', name = "ktx"): 0.039,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_AmeR_AmeR_degtagged', name = "ktx"): 0.0023,
    
    # AmtR Parameters
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_AmtR_AmtR_degtagged', name = 'kb'):100.93, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_AmtR_AmtR_degtagged', name = 'ku'):0.75, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_AmtR_AmtR_degtagged', name = 'cooperativity'):1.01, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_AmtR_leak', name = "ktx"):0.037,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_AmtR_AmtR_degtagged', name = "ktx"): 0.00079,
    
    # BetI Parameters
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_BetI_BetI_degtagged', name = 'kb'):100.29, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_BetI_BetI_degtagged', name = 'ku'):6.89, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_BetI_BetI_degtagged', name = 'cooperativity'):2.27, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_BetI_leak', name = "ktx"):0.044,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_BetI_BetI_degtagged', name = "ktx"): 0.00083,
    
    # HlyIIR Parameters
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_HlyIIR_HlyIIR_degtagged', name = 'kb'):102.43, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_HlyIIR_HlyIIR_degtagged', name = 'ku'):0.53, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_HlyIIR_HlyIIR_degtagged', name = 'cooperativity'):2.15, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_HlyIIR_leak', name = "ktx"):0.030,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_HlyIIR_HlyIIR_degtagged', name = "ktx"): 0.00081,
    
    # PhlF Parameters
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_PhlF_PhlF_degtagged', name = 'kb'):91.74, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_PhlF_PhlF_degtagged', name = 'ku'):0.0001, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_PhlF_PhlF_degtagged', name = 'cooperativity'):3.72, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_PhlF_leak', name = "ktx"):0.050,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_PhlF_PhlF_degtagged', name = "ktx"): 0.00024,
    
    # SrpR Parameters
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_SrpR_SrpR_degtagged', name = 'kb'):101.12, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_SrpR_SrpR_degtagged', name = 'ku'):0.25, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = 'P_SrpR_SrpR_degtagged', name = 'cooperativity'):1.02, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_SrpR_leak', name = "ktx"):0.025,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_SrpR_SrpR_degtagged', name = "ktx"): 6.75e-5,
    
    #AraAraC Bound Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name] 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_BAD_protein_AraC_input_2x_protein_Ara_input_2x', name = "ktx"): 0.05,
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Tac_leak', name = "ktx"): 0.05, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Tet_leak', name = "ktx"): 0.05 
    
}

#Species

IPTG = Species('IPTG',  material_type='protein', attributes=['input']) #Input A
LacI = Species('LacI',  material_type='protein', attributes=['input'])
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

aTc = Species('aTc',  material_type='protein', attributes=['input']) #Input B
TetR = Species('TetR',  material_type='protein', attributes=['input'])
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR], parameters = complex_parameters)

Ara = Species('Ara',  material_type='protein', attributes=['input']) #Input C
AraC = Species('AraC',  material_type='protein', attributes=['input'])
AraAraC = ChemicalComplex([Ara, Ara, AraC, AraC], parameters = complex_parameters)

#DNA parts

rbs = RBS('UTR1')
f1 = RBS('F1')

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
P_BAD = RegulatedPromoter('P_BAD', regulators = [AraAraC], leak=True,
                          parameters = component_parameters)

P_SrpR = RegulatedPromoter('P_SrpR',  regulators = ['SrpR_degtagged'], leak=True, 
                          parameters = component_parameters)
P_BetI = RegulatedPromoter('P_BetI',  regulators = ['BetI_degtagged'], leak=True, 
                          parameters = component_parameters)
P_HlyIIR = RegulatedPromoter('P_HlyIIR',  regulators = ['HlyIIR_degtagged'], leak=True, 
                          parameters = component_parameters)
P_AmeR = RegulatedPromoter('P_AmeR',  regulators = ['AmeR_degtagged'], leak=True, 
                          parameters = component_parameters)
P_AmtR = RegulatedPromoter('P_AmtR', regulators = ['AmtR_degtagged'], leak=True, 
                          parameters = component_parameters)
P_PhlF = RegulatedPromoter('P_PhlF',  regulators = ['PhlF_degtagged'], leak=True, 
                          parameters = component_parameters)

#DNA_constructs

PhlF_construct = DNA_construct([P_SrpR, P_BetI, rbs, CDS_PhlF, t16])

SrpR_construct = DNA_construct([P_Tac, P_Tet, rbs, CDS_SrpR, t16])

BetI_construct = DNA_construct([P_HlyIIR, P_AmeR, rbs, CDS_BetI, t16])

AmeR_construct = DNA_construct([P_Tet, rbs, CDS_AmeR, t16])

HlyIIR_construct = DNA_construct([P_Tac, rbs, CDS_HlyIIR, t16])

AmtR_construct = DNA_construct([P_BAD, rbs, CDS_AmtR, t16])

YFP_construct = DNA_construct([P_PhlF, P_AmtR, rbs, CDS_YFP, t16])

#Mixture and CRN creation

dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)

global_mechanisms = {"dilution":dilution_mechanism}

# degredation_mechanism = Deg_Tagged_Degredation(protease)

# global_mechanisms = {"degredation":degredation_mechanism}

M = SimpleTxTlExtract(name="txtl", parameters = parameters, global_mechanisms = global_mechanisms,
                      components=[PhlF_construct, SrpR_construct, BetI_construct, AmeR_construct, 
                                  HlyIIR_construct, AmtR_construct, YFP_construct, IPTG_LacI, 
                                  aTc_TetR, AraAraC])
CRN = M.compile_crn()
# CRN.write_sbml_file('Circuit_0xF6_AB_only_sbml.xml') #saving CRN as sbml

with open('temp_CRN_EQNs.txt', 'w') as f:
    f.write(CRN.pretty_print(show_rates = True, show_keys = True))

print('CRN Compiled')

sim = GCSim(CRN)

protein_lst = ['protein_PhlF_degtagged', 'protein_AmtR_degtagged', 'protein_AmtR_degtagged']

#Plotting
for a in [0,8]:
    for b in [0,8]:
        for c in [0,8]:

            x0 = {PhlF_construct.get_species():0.1, SrpR_construct.get_species():0.1, 
                  BetI_construct.get_species():0.1, AmeR_construct.get_species():0.1, 
                  HlyIIR_construct.get_species():0.1, YFP_construct.get_species():0.1, 
                  AmtR_construct.get_species():0.1, Ara:c, AraC:c,
                  IPTG:a, LacI:5, aTc:b, TetR:5}
            timepoints = np.linspace(0, 10000, 10000)
            R = sim.basicsim(x0, timepoints, protein_lst, title = f'IPTG = {a}, aTc = {b}, Ara = {c}')
            print(R['protein_YFP_degtagged'].iloc[-1])
            # R[['time'] + protein_lst].to_excel(f'{a}, {b}, {c}.xlsx', index=False)
            # R.to_csv(f'simulation_data/IPTG_{a}_aTc_{b}_Ara_{c}.csv', index=False)