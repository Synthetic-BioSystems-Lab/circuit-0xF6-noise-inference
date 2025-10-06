# -*- coding: utf-8 -*-
"""
Created on Fri Oct  3 16:54:35 2025

@author: zacha
"""

from biocrnpyler import *
import numpy as np
import matplotlib.pyplot as plt
import libsbml

# Parameters and Global Mechanisms
parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.05, "kdeg":0.001, "kdil":0.0075}
complex_parameters = {'kb':100, 'ku':10}
component_parameters = {
    
    #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):100, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):10, 
    ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):4, 
    
    #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'simple_transcription', part_id = None, name = "ktx"): 1e-15,
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Tac_leak', name = "ktx"): 0.05
    
}

dilution_mechanism = Dilution(filter_dict = {"degtagged":True}, default_on = False)

global_mechanisms = {"dilution":dilution_mechanism}

# Species

IPTG = Species('IPTG',  material_type='protein', attributes=['input']) 
LacI = Species('LacI',  material_type='protein', attributes=['input'])
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

# DNA parts

P_Tac = RegulatedPromoter('P_Tac',  regulators = [LacI], leak=True, 
                          parameters = component_parameters)

rbs = RBS('rbs')

CDS_YFP = CDS('YFP', 'YFP')
CDS_YFP.protein = Species('YFP', material_type='protein', attributes=['degtagged'])

t = Terminator('t')

YFP_construct = DNA_construct([P_Tac, rbs, CDS_YFP, t])

# Mixture and CRN creation

M = SimpleTxTlExtract('simtxtl', parameters = parameters, global_mechanisms=global_mechanisms, 
                      components=[YFP_construct, IPTG_LacI])

CRN = M.compile_crn()
CRN.write_sbml_file('temp_sbml_file.xml') #saving CRN as sbml

with open('temp_CRN_EQNs.txt', 'w') as f:
    f.write(CRN.pretty_print(show_rates = True, show_keys = True))
    
timepoints = np.linspace(0, 1500, 500)

x0 = {YFP_construct.get_species():1, LacI:10, IPTG:20}

# Simulation and Plotting    
R = CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0)

plot_lst = ['protein_YFP_degtagged', 'protein_IPTG_input', 'protein_LacI_input']

plt.figure()

for i in range(len(plot_lst)):
    plt.plot(R['time'], R[plot_lst[i]], label = f'{plot_lst[i]}')

plt.legend()
plt.title('')
plt.xlabel('Time')
plt.ylabel('YFP')

plt.show()

print(R['protein_YFP_degtagged'].iloc[-1])

reader = libsbml.SBMLReader()
document = reader.readSBML('temp_sbml_file.xml')

sbml_string = libsbml.writeSBMLToString(document)

with open('filex.txt', 'w') as f:
    f.write(sbml_string)
