# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 13:17:44 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import pylab as plt
import pandas as pd

#Parameters #!!! parameters mostly default currently #!!! add in degradation???

hill_parameters = {"k":1.0, "n":4, "K":20, "kleak":.01}
complex_parameters = {'kb':100, 'ku':1.0}
parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
component_parameters = {
    #Promoter AraAraC Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'binding', part_id = 'P_BAD_AraAraC', name = 'kb'):100, 
    ParameterKey(mechanism = 'binding', part_id = "P_BAD_AraAraC", name = 'ku'):5.0, 
    ParameterKey(mechanism = 'binding', part_id = "P_BAD_AraAraC", name = 'cooperativity'):4.0, 
    
    #AraAraC Bound Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription', part_id = 'P_BAD_AraAraC', name = 'kb'):10, 
    ParameterKey(mechanism = 'transcription', part_id = "P_BAD_AraAraC", name = 'ku'):1, 
    ParameterKey(mechanism = 'transcription', part_id = 'P_BAD_AraAraC', name = "ktx"): 1,
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription', part_id = 'P_AraAraC_Leak', name = "kb"): 2,
    ParameterKey(mechanism = 'transcription', part_id = 'P_AraAraC_Leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription', part_id = 'P_AraAraC_Leak', name = "ktx"): 1.0,
    
    #Promoter AmtR Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'binding', part_id = 'P_AmtR_AmtR', name = 'kb'):100,
    ParameterKey(mechanism = 'binding', part_id = "P_AmtR_AmtR", name = 'ku'):5.0,
    ParameterKey(mechanism = 'binding', part_id = "P_AmtR_AmtR", name = 'cooperativity'):4.0,
    
    #AmtR Bound Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription', part_id = 'P_AmtR_AmtR', name = 'kb'):1,
    ParameterKey(mechanism = 'transcription', part_id = "P_AmtR_AmtR", name = 'ku'):10,
    ParameterKey(mechanism = 'transcription', part_id = 'P_AmtR_AmtR', name = "ktx"):1, 
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription', part_id = 'P_AmtR_Leak', name = "kb"): 2,
    ParameterKey(mechanism = 'transcription', part_id = 'P_AmtR_Leak', name = "ku"): 10,
    ParameterKey(mechanism = 'transcription', part_id = 'P_AmtR_Leak', name = "ktx"):1.0 
}


#Species

Ara = Species('Ara') #Input C
AraC = Species('AraC')
AraAraC = ChemicalComplex([Ara, Ara, AraC, AraC], parameters = complex_parameters)

AmtR = Species('AmtR')
YFP = Species('YFP')

rbs = RBS('UTR1')
CDS_AmtR = CDS('AmtR', 'AmtR')
CDS_YFP = CDS('YFP', 'YFP')
t16 = Terminator('t16')

#Promoters

P_BAD = RegulatedPromoter('P_BAD', regulators = [AraAraC], leak=False,
                          parameters = complex_parameters)

P_AmtR = RegulatedPromoter('P_AmtR', regulators = [AmtR], leak=True, 
                          parameters = complex_parameters)

#DNA_constructs
mechanisms = {"transcription":Transcription_MM(Species("RNAP",material_type="protein")), 
              "translation":Translation_MM(Species("Ribo",material_type="protein")),
              "binding":One_Step_Cooperative_Binding()}

AmtR_construct = DNA_construct([P_BAD, rbs, CDS_AmtR, t16], mechanisms = mechanisms
                               ,parameters=parameters)
YFP_construct = DNA_construct([P_AmtR, rbs, CDS_YFP, t16], mechanisms = mechanisms
                               ,parameters=parameters)

#Mixture and CRN creation

M = TxTlExtract(name='TxTl', parameters=parameters, 
                      components=[AmtR_construct, YFP_construct, AraAraC])

CRN = M.compile_crn()
CRN.write_sbml_file('Ara_Signal_Inverter_sbml.xml') #saving CRN as sbml

print(CRN.pretty_print(show_rates = True, show_keys = True))

#Plotting

for a_c in [0, 25, 50, 100, 200, 1000]:
    x0 = {AmtR_construct.get_species():20, YFP_construct.get_species():20, Ara:a_c, AraC:a_c, 
          "protein_RNAP":10., "protein_Ribo":50., 'protein_RNAase':10}
    timepoints = np.linspace(0, 100, 100)
    R = CRN.simulate_with_bioscrape_via_sbml( timepoints, initial_condition_dict = x0)
    plt.plot(R['time'], R['protein_YFP'], label = 'Initial [Ara], [AraC] = ' + str(a_c) )

plt.ylabel('[YFP]')
plt.xlabel('Time')
plt.legend()
plt.show()

for mech_type, mech in M.mechanisms.items():
    print(f"{mech_type}: {mech}")