# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 14:10:00 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import pylab as plt
import pandas as pd
from GCSim import GCSim

'See Excel Sheet Hill EQs Fits for Source for values of V_max, K, and n'

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
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Tac_leak', name = "ktx"): 0.05, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Tet_leak', name = "ktx"): 0.05,
    
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_SrpR_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_BetI_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_HlyIIR_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_AmeR_leak', name = "ktx"): 0.05, 
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_PhlF_leak', name = "ktx"): 0.05,
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_AmtR_leak', name = "ktx"):0.05 
    
}

#Species

IPTG = Species('IPTG',  material_type='protein', attributes=['input']) #Input A
LacI = Species('LacI',  material_type='protein', attributes=['input'])
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

aTc = Species('aTc',  material_type='protein', attributes=['input']) #Input B
TetR = Species('TetR',  material_type='protein', attributes=['input'])
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR], parameters = complex_parameters)

#DNA parts

rbs = RBS('UTR1')

CDS_SrpR = CDS('SrpR', 'SrpR')
CDS_SrpR.protein = Species('SrpR', material_type='protein', attributes=['degtagged'])


t16 = Terminator('t16')

P_Tac = RegulatedPromoter('P_Tac',  regulators = [LacI], leak=True, 
                          parameters = component_parameters)
P_Tet = RegulatedPromoter('P_Tet', regulators = [TetR], leak=True, 
                          parameters = component_parameters)

#DNA_constructs


SrpR_construct = DNA_construct([P_Tac, P_Tet, rbs, CDS_SrpR, t16])


#Mixture and CRN creation

dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)

global_mechanisms = {"dilution":dilution_mechanism}



M = SimpleTxTlExtract(name="txtl", parameters = parameters, global_mechanisms = global_mechanisms,
                      components=[SrpR_construct])
CRN = M.compile_crn()
# CRN.write_sbml_file('Circuit_0xF6_AB_only_sbml.xml') #saving CRN as sbml

with open('temp_CRN_EQNs.txt', 'w') as f:
    f.write(CRN.pretty_print(show_rates = True, show_keys = True))

print('CRN Compiled')

sim = GCSim(CRN)

dna_lst = ['dna_part_P_Tac_forward_part_P_Tet_forward_part_UTR1_forward_part_SrpR_forward_part_t16_forward_',
                'ordered_polymer_complex_part_P_Tac_protein_LacI_input_4x_forward__part_P_Tet_forward_part_UTR1_forward_part_SrpR_forward_part_t16_forward_',
                'ordered_polymer_complex_part_P_Tac_protein_LacI_input_4x_forward__complex_part_P_Tet_protein_TetR_input_4x_forward__part_UTR1_forward_part_SrpR_forward_part_t16_forward_',
                'ordered_polymer_part_P_Tac_forward_complex_part_P_Tet_protein_TetR_input_4x_forward__part_UTR1_forward_part_SrpR_forward_part_t16_forward_']

# dna_lst = [
               
#                'ordered_polymer_complex_part_P_Tac_protein_LacI_input_4x_forward__complex_part_P_Tet_protein_TetR_input_4x_forward__part_UTR1_forward_part_SrpR_forward_part_t16_forward_'
#                ]

IPTG_vals = np.linspace(0, 10, 11)
aTc_vals = np.linspace(0, 10, 11)

heatmaps = {dna: np.zeros((len(IPTG_vals), len(aTc_vals))) for dna in dna_lst}

#Plotting
for i, a in enumerate(IPTG_vals):
    for j, b in enumerate(aTc_vals):

        x0 = {SrpR_construct.get_species():1, 
              LacI:a, TetR:b}
        timepoints = np.linspace(0, 10, 10)
        
        R = CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, safe=True)

        # R[['time'] + protein_lst].to_excel(f'{a}, {b}, {c}.xlsx', index=False)
        # R.to_excel(f'simulation_results_IPTG_{a}_aTc_{b}.xlsx', index=False)
        
        for dna in dna_lst:
                heatmaps[dna][i, j] = R[dna].iloc[-1] 
                
title_lst = ['unbound dna', 'LacI only', 'both repressors', 'TetR only']
# title_lst = ['both repressors']
                
for i,dna in enumerate(dna_lst):
    
    plt.figure()
    plt.title(f'{title_lst[i]}')
    cb = plt.pcolor(heatmaps[dna], cmap='viridis')
    plt.colorbar(cb)
    plt.xlabel('TetR') 
    plt.ylabel('LacI') 
    

    plt.xticks(np.arange(0.5, len(aTc_vals)+0.5, 1), [f"{b:.1f}" for b in aTc_vals])
    plt.yticks(np.arange(0.5, len(IPTG_vals)+0.5, 1), [f"{a:.1f}" for a in IPTG_vals])
        
    plt.show()
    
LacI_vals = np.linspace(0, 30, 31)
TetR_vals = np.linspace(0, 30, 31)
        
heat = np.zeros((len(LacI_vals), len(TetR_vals)))
heat2 = np.zeros((len(LacI_vals), len(TetR_vals)))

Vmax = 0.5
n = 2 
K = 2    

for i, a in enumerate(LacI_vals):
    for j, b in enumerate(TetR_vals):
        
        if a >= 1 and b >= 1:
            repressed_dna_conc = ( (Vmax*((a-1)**n))/((K**n)+((a-1)**n)) ) + \
            ( (Vmax*((b-1)**n))/((K**n)+((b-1)**n)) )
        
        else:
            repressed_dna_conc = 0
            
        unrepressed_dna_conc = 1 - repressed_dna_conc
        
        heat[i,j] = repressed_dna_conc
        heat2[i,j] = unrepressed_dna_conc
        
plt.figure()
plt.title('Hill EQs')
cb = plt.pcolor(heat[:11,:11], cmap='viridis')
plt.colorbar(cb)
plt.xlabel('TetR') 
plt.ylabel('LacI') 


plt.xticks(np.arange(0.5, len(aTc_vals)+0.5, 1), [f"{b:.1f}" for b in aTc_vals])
plt.yticks(np.arange(0.5, len(IPTG_vals)+0.5, 1), [f"{a:.1f}" for a in IPTG_vals])
    
plt.show()   


#3d plot

LacI_mesh, TetR_mesh = np.meshgrid(LacI_vals, TetR_vals) 

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(LacI_mesh, TetR_mesh, heat2, cmap='viridis', alpha=0.9)

ax.view_init(30, 15)
ax.set_xlabel('TetR')
ax.set_ylabel('LacI')
ax.set_zlabel('Repressed DNA Conc')
ax.set_title('3D Hill Equation Surface')

fig.colorbar(surf)

plt.show()
