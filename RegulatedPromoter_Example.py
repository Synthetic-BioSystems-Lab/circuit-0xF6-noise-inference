# -*- coding: utf-8 -*-
"""
Created on Sun Jul 13 16:49:32 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import pylab as plt
import pandas as pd

#1 Regulated Promoter Needs lots of parameters!
component_parameters = {
    #Promoter repressor2 Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'binding', part_id = 'regulated_promoter_A', name = 'kb'):100, #Promoter - repressor2 Binding
    ParameterKey(mechanism = 'binding', part_id = "regulated_promoter_A", name = 'ku'):5.0, #Unbinding
    ParameterKey(mechanism = 'binding', part_id = "regulated_promoter_A", name = 'cooperativity'):4.0, #Cooperativity
    
    #Activated Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    #These regulate RNAP binding to an activated promoter and transcription
    ParameterKey(mechanism = 'transcription', part_id = 'regulated_promoter_A', name = 'kb'):1, #Promoter - repressor2 Binding
    ParameterKey(mechanism = 'transcription', part_id = "regulated_promoter_A", name = 'ku'):100, #Unbinding
    ParameterKey(mechanism = 'transcription', part_id = 'regulated_promoter_A', name = "ktx"): 1., #Transcription Rate
    
    #Promoter Repressor Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'binding', part_id = 'regulated_promoter_R', name = 'kb'):100,
    ParameterKey(mechanism = 'binding', part_id = "regulated_promoter_R", name = 'ku'):5.0,
    ParameterKey(mechanism = 'binding', part_id = "regulated_promoter_R", name = 'cooperativity'):4.0,
    
    #Repressed Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    #These regulate RNAP binding to a repressed promoter and transcription
    ParameterKey(mechanism = 'transcription', part_id = 'regulated_promoter_R', name = 'kb'):1,
    ParameterKey(mechanism = 'transcription', part_id = "regulated_promoter_R", name = 'ku'):100,
    ParameterKey(mechanism = 'transcription', part_id = 'regulated_promoter_R', name = "ktx"): 1.0, #Transcription Rate
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription', part_id = 'regulated_promoter_leak', name = "kb"): 2,
    ParameterKey(mechanism = 'transcription', part_id = 'regulated_promoter_leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription', part_id = 'regulated_promoter_leak', name = "ktx"): 1.0, #Transcription Rate
}

repressor = Species("R", material_type = "protein")
repressor2 = Species("A", material_type = "protein")
reporter = Species("reporter", material_type = "protein")

#Create a RegulatedPromoter Object named "P_reg" with regulators "repressor2" and "repressor"
#By Loading custom parameters into the promoter, we override the default parameters of the Mixture
P_reg = RegulatedPromoter("regulated_promoter", regulators=[repressor2, repressor], leak=True, 
                          parameters = component_parameters)

#Create a DNA assembly "reporter" with P_reg for its promoter
reg_reporter = DNAassembly(name="reporter", promoter=P_reg, rbs="strong", protein = reporter)

#Use a simple TxTl model with dilution
#M = SimpleTxTlDilutionMixture(name="e coli", parameter_file = "default_parameters.txt", components=[reg_reporter])
M = TxTlExtract(name="e coli extract", parameter_file = "default_parameters.txt", components=[reg_reporter])

CRN = M.compile_crn()
print(CRN.pretty_print(show_rates = True, show_keys = False))

plt.figure(figsize = (6, 6))
N = 11 #Number of titrations
max_titration = 10
HM = np.zeros((N, N))
for a_ind, A_c in enumerate(np.linspace(0, max_titration, N)):
    for b_ind, B_c in enumerate(np.linspace(0, max_titration, N)):
        x0 = {reg_reporter.dna:1, repressor:A_c, repressor2:B_c}
        timepoints = np.linspace(0, 1000, 1000)
        R = CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0)
        HM[a_ind, b_ind] = R["protein_reporter"][len(timepoints)-1]
plt.title("")
cb = plt.pcolor(HM)
plt.colorbar(cb)
plt.xlabel("repressor2")
plt.ylabel("repressor")

plt.xticks(np.arange(.5, N+.5, 1), [str(i) for i in np.linspace(0, max_titration, N)])
plt.yticks(np.arange(.5, N+.5, 1), [str(i) for i in np.linspace(0, max_titration, N)])
    
plt.show()
