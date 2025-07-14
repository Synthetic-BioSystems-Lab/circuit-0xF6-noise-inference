# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 16:57:24 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import pylab as plt
import pandas as pd

#Parameters #!!! parameters mostly default currently #!!! add in degradation???

hill_parameters = {"k":1.0, "n":4, "K":20, "kleak":.01}
complex_parameters = {'kb':0.5, 'ku':1.0}
component_parameters = {
    #Promoter LacI Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'binding', part_id = 'P_Tac_Tet_LacI', name = 'kb'):100, 
    ParameterKey(mechanism = 'binding', part_id = "P_Tac_Tet_LacI", name = 'ku'):5.0, 
    ParameterKey(mechanism = 'binding', part_id = "P_Tac_Tet_LacI", name = 'cooperativity'):4.0, 
    
    #LacI Bound Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_LacI', name = 'kb'):1, 
    ParameterKey(mechanism = 'transcription', part_id = "P_Tac_Tet_LacI", name = 'ku'):100, 
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_LacI', name = "ktx"): 1., 
    
    #Promoter TetR Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'binding', part_id = 'P_Tac_Tet_TetR', name = 'kb'):100,
    ParameterKey(mechanism = 'binding', part_id = "P_Tac_Tet_TetR", name = 'ku'):5.0,
    ParameterKey(mechanism = 'binding', part_id = "P_Tac_Tet_TetR", name = 'cooperativity'):4.0,
    
    #TetR Bound Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_TetR', name = 'kb'):1,
    ParameterKey(mechanism = 'transcription', part_id = "P_Tac_Tet_TetR", name = 'ku'):100,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_TetR', name = "ktx"): 1.0, 
    
    #Leak Parameters for transcription
    #These regulate expression of an unbound promoter
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_leak', name = "kb"): 2.,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_leak', name = "ku"): 100,
    ParameterKey(mechanism = 'transcription', part_id = 'P_Tac_Tet_leak', name = "ktx"): 1.0, 
}

#Species

IPTG = Species('IPTG') #Input A
LacI = Species('LacI')
IPTG_LacI = ChemicalComplex([IPTG, IPTG, LacI, LacI], parameters = complex_parameters)

aTc = Species('aTc') #Input B
TetR = Species('TetR')
aTc_TetR = ChemicalComplex([aTc, aTc, TetR, TetR], parameters = complex_parameters)

Ara = Species('Ara') #Input C
AraC = Species('AraC')
AraAraC = ChemicalComplex([Ara, Ara, AraC, AraC], parameters = complex_parameters)

dummy_species = Species('dummy_species')

HlyIIR = Species('HlyIIR') #Gate Protein
SrpR = Species('SrpR') #Gate Protein
AmeR = Species('AmeR') #Gate Protein
BetI = Species('BetI') #Gate Protein
PhlF = Species('PhlF') #Gate Protein
AmtR = Species('AmtR') #Gate Protein
YFP = Species('YFP') #Reporter Protein

#Promoters

P_Tac = RepressiblePromoter('P_Tac', repressor = LacI, leak= False, 
                            parameters = hill_parameters)
P_Tet = RepressiblePromoter('P_Tet', repressor = TetR, leak= False,
                            parameters = hill_parameters)
P_Tac_Tet = RegulatedPromoter('P_Tac_Tet', regulators = [LacI,TetR], leak=True, 
                              parameters = component_parameters)
P_BAD = ActivatablePromoter('P_BAD', activator = AraAraC, leak= False, 
                            parameters = hill_parameters)

P_SrpR_BetI = CombinatorialPromoter('P_SrpR_BetI', regulators= [SrpR, BetI, dummy_species], 
                              tx_capable_list= 
                              [[dummy_species]], 
                              parameters = hill_parameters , leak= False)
P_HlyIIR_AmeR = CombinatorialPromoter('P_HlyIIR_AmeR', regulators= [HlyIIR, AmeR, dummy_species], 
                              tx_capable_list= 
                              [[dummy_species], [HlyIIR, dummy_species], [AmeR, dummy_species]], 
                              parameters = hill_parameters , leak= False)

P_PhlF_AmtR = CombinatorialPromoter('P_PhlF_AmtR', regulators= [PhlF, AmtR, dummy_species], 
                              tx_capable_list= 
                              [[dummy_species], [dummy_species, PhlF], [dummy_species, AmtR]], 
                              parameters = hill_parameters , leak= False)

#DNAassemblies

PhlF_assembly = DNAassembly(name='PhlF_assembly', promoter=P_SrpR_BetI, rbs = 'strong', 
                            protein = PhlF)
SrpR_assembly = DNAassembly(name='SrpR_assembly', promoter=P_Tac_Tet, rbs = 'strong', 
                             protein = SrpR)
BetI_assembly = DNAassembly(name='BetI_assembly', promoter=P_HlyIIR_AmeR, rbs = 'strong', 
                             protein = BetI)
AmeR_assembly = DNAassembly(name='AmeR_assembly', promoter=P_Tet, rbs = 'strong', 
                             protein = AmeR)
HlyIIR_assembly = DNAassembly(name='HlyIIR_assembly', promoter=P_Tac, rbs = 'strong', 
                             protein = HlyIIR)
AmtR_assembly = DNAassembly(name='AmtR_assembly', promoter=P_BAD, rbs = 'strong', 
                            protein = AmtR)
reporter_assembly = DNAassembly(name='reporter_assembly', promoter=P_PhlF_AmtR, 
                                rbs = 'strong', protein = YFP)

#Mixture and CRN creation

# parameters={"cooperativity":2, "kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2, "kint":.05, 'kdil':0.0075}
M = SimpleTxTlDilutionMixture(name='SimpleTxTl', parameter_file = 'default_parameters.txt', 
                      components=[PhlF_assembly, SrpR_assembly, BetI_assembly, AmeR_assembly, HlyIIR_assembly, 
                                  AmtR_assembly, reporter_assembly, AraAraC, aTc_TetR, IPTG_LacI])

CRN = M.compile_crn()
CRN.write_sbml_file('Circuit_0xF6_bioCRNpyler_sbml.xml') #saving CRN as sbml

print(CRN.pretty_print(show_rates = True, show_keys = True))

#Plotting
for a in [0, 1]:
    for b in [0, 1]:
        for c in [0, 1]:
            x0 = {AmtR_assembly.dna:1, reporter_assembly.dna:1, PhlF_assembly.dna:1, SrpR_assembly.dna:1, 
                  BetI_assembly.dna:1, AmeR_assembly.dna:1, HlyIIR_assembly.dna:1, 
                  IPTG:a, LacI:1, aTc:b, TetR:1, Ara:c, AraC:c, dummy_species:1000}
            timepoints = np.linspace(0, 100, 100)
            R = CRN.simulate_with_bioscrape_via_sbml( timepoints, initial_condition_dict = x0)
            plt.plot(R['time'], R['SrpR'], label = 'SrpR')
            # plt.plot(R['time'], R['PhlF'], label = 'PhlF')
            plt.ylabel(f'Concentration')
            plt.xlabel('Time')
            plt.title(f'A = {a}, B = {b}, C = {c}')
            plt.legend()
            plt.show()

# for mech_type, mech in M.mechanisms.items():
#     print(f"{mech_type}: {mech}")
