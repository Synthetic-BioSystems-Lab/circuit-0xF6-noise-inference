# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 16:47:24 2025

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
hill_parameters = {"k":1.0, "n":4, "K":20, "kleak":.01}
component_parameters = {
    #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):10, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):2, 
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'simple_transcription', part_id = None, name = "ktx"): 0.00001,
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'simple_transcription', part_id = 'Promoter_leak', name = "ktx"): 0.05
    
}

#Species

repressor = Species('repressor',  material_type='protein', attributes=['input'])


#DNA parts

rbs = RBS('UTR1')

CDS_A= CDS('A', 'A')
CDS_A.protein = Species('A', material_type='protein', attributes=['degtagged'])


t16 = Terminator('t16')

# Promoter = RepressiblePromoter('Promoter', repressor, parameters= hill_parameters)

Promoter = RegulatedPromoter('Promoter',  regulators = [repressor], leak=True, 
                          parameters = component_parameters)

#DNA_constructs

construct_0 = DNA_construct([Promoter, rbs, CDS_A, t16])


#Mixture and CRN creation

dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)

global_mechanisms = {"dilution":dilution_mechanism}

# degredation_mechanism = Deg_Tagged_Degredation(protease)

# global_mechanisms = {"degredation":degredation_mechanism}

M = SimpleTxTlExtract(name="txtl", parameters = parameters, global_mechanisms = global_mechanisms,
                      components=[construct_0])
CRN = M.compile_crn()
# CRN.write_sbml_file('Circuit_0xF6_AB_only_sbml.xml') #saving CRN as sbml

with open('temp_CRN_EQNs.txt', 'w') as f:
    f.write(CRN.pretty_print(show_rates = True, show_keys = True))

print('CRN Compiled')

sim = GCSim(CRN)

#Plotting

timepoints = np.linspace(0, 6000, 6000)
x_lst = np.linspace(0.01, 1000, 1000)

for i in [0, 10]:
    x0 = {construct_0.get_species():1, 'protein_repressor_input':i}
    sim.basicsim(x0, timepoints, ['protein_repressor_input', 'protein_A_degtagged'])

x0 = {construct_0.get_species():1}

# sim.inout(x0, timepoints, x_lst, 'protein_repressor_input', 'protein_A_degtagged', 
#           loglog = True)

x_lst = np.linspace(0.001, 10, 100)

sim.inout(x0, timepoints, x_lst, 'protein_repressor_input', 'protein_A_degtagged')

