# -*- coding: utf-8 -*-
"""
Created on Fri Aug 22 16:45:43 2025

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
    #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):10, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):4, 
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'simple_transcription', part_id = None, name = "ktx"): 1e-15,
    
    #AraAraC Bound Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name] 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_BAD_protein_AraC_input_2x_protein_Ara_input_2x', name = "ktx"): 0.05,
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_BAD_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Tac_leak', name = "ktx"): 0.05, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Tet_leak', name = "ktx"): 0.05,
    
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_SrpR_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_BetI_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_HlyIIR_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_AmeR_leak', name = "ktx"): 0.05, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_PhlF_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Inv1_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Inv2_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_AmtR_leak', name = "ktx"):0.05 
    
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

CDS_Inv1 = CDS('Inv1', 'Inv1')
CDS_Inv1.protein = Species('Inv1', material_type='protein', attributes=['degtagged'])

CDS_Inv2 = CDS('Inv2', 'Inv2')
CDS_Inv2.protein = Species('Inv2', material_type='protein', attributes=['degtagged'])

CDS_YFP = CDS('YFP', 'YFP')
CDS_YFP.protein = Species('YFP', material_type='protein', attributes=['degtagged'])

t16 = Terminator('t16')

P_Tac = RegulatedPromoter('P_Tac',  regulators = [LacI], leak=True, 
                          parameters = component_parameters)
P_Tet = RegulatedPromoter('P_Tet', regulators = [TetR], leak=True, 
                          parameters = component_parameters)
P_BAD = RegulatedPromoter('P_BAD', regulators = [AraAraC], leak=False,
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
P_Inv1 = RegulatedPromoter('P_Inv1', regulators = ['Inv1_degtagged'], leak=True, 
                          parameters = component_parameters)
P_Inv2 = RegulatedPromoter('P_Inv2', regulators = ['Inv2_degtagged'], leak=True, 
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

Inv1_construct = DNA_construct([P_AmtR, rbs, CDS_Inv1, t16])

Inv2_construct = DNA_construct([P_Inv1, rbs, CDS_Inv2, t16])

YFP_construct = DNA_construct([P_PhlF, P_Inv2, rbs, CDS_YFP, t16])

#Mixture and CRN creation

dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)

global_mechanisms = {"dilution":dilution_mechanism}

# degredation_mechanism = Deg_Tagged_Degredation(protease)

# global_mechanisms = {"degredation":degredation_mechanism}

M = SimpleTxTlExtract(name="txtl", parameters = parameters, global_mechanisms = global_mechanisms,
                      components=[PhlF_construct, SrpR_construct, BetI_construct, AmeR_construct, 
                                  HlyIIR_construct, AmtR_construct, YFP_construct, Inv1_construct, 
                                  Inv2_construct, IPTG_LacI, aTc_TetR, AraAraC])
CRN = M.compile_crn()
# CRN.write_sbml_file('Circuit_0xF6_AB_only_sbml.xml') #saving CRN as sbml

with open('temp_CRN_EQNs.txt', 'w') as f:
    f.write(CRN.pretty_print(show_rates = True, show_keys = True))

print('CRN Compiled')

sim = GCSim(CRN)

protein_lst = ['protein_PhlF_degtagged', 'protein_YFP_degtagged', 
               'protein_Inv2_degtagged']

#Plotting
for a in [0,80]:
    for b in [0,80]:
        for c in [0,80]:

            x0 = {PhlF_construct.get_species():1, SrpR_construct.get_species():1, 
                  BetI_construct.get_species():1, AmeR_construct.get_species():1, 
                  HlyIIR_construct.get_species():1, YFP_construct.get_species():1, 
                  AmtR_construct.get_species():1, Inv1_construct.get_species():1, 
                  Inv2_construct.get_species():1, Ara:c, AraC:c,
                  IPTG:a, LacI:50, aTc:b, TetR:50}
            timepoints = np.linspace(0, 6000, 6000)
            R = sim.basicsim(x0, timepoints, protein_lst, title = f'IPTG = {a}, aTc = {b}, Ara = {c}')
            print(R['protein_YFP_degtagged'].iloc[-1])
            # R[['time'] + protein_lst].to_excel(f'{a}, {b}, {c}.xlsx', index=False)
            # R.to_excel(f'simulation_data/IPTG_{a}_aTc_{b}_Ara_{c}.xlsx', index=False)