# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 19:47:16 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import pylab as plt
import pandas as pd

#Parameters #!!! parameters mostly default currently #!!! add in degradation???

hill_parameters = {"k":1.0, "n":4, "K":20, "kleak":.01}
complex_parameters = {'kb':0.05, 'ku':1.0}

#Species

Ara = Species('Ara') #Input C
AraC = Species('AraC')
AraAraC = ChemicalComplex([Ara, Ara, AraC, AraC], parameters = complex_parameters)

AmtR = Species('AmtR')
YFP = Species('YFP')

#Promoters

P_BAD = ActivatablePromoter('P_BAD', activator = AraAraC, leak = True, 
                            parameters = hill_parameters)

P_AmtR = RepressiblePromoter('P_AmtR', AmtR, leak = True, 
                            parameters = hill_parameters)

#DNAassemblies

main_assembly = DNAassembly(name='main_assembly', promoter=P_BAD, rbs = 'strong', 
                            protein = AmtR)
reporter_assembly = DNAassembly(name='reporter_assembly', promoter=P_AmtR, 
                                rbs = 'strong', protein = YFP)

#Mixture and CRN creation

parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05, 'kdil':0.001}
M = SimpleTxTlExtract(name='SimpleTxTl', parameters=parameters, 
                      components=[main_assembly, reporter_assembly, AraAraC])

CRN = M.compile_crn()
CRN.write_sbml_file('Ara_Signal_Inverter_sbml.xml') #saving CRN as sbml

print(CRN.pretty_print(show_rates = True, show_keys = True))

#Plotting

for a_c in [0, 25, 50, 100, 200, 1000]:
    x0 = {main_assembly.dna:1, reporter_assembly.dna: 1, Ara:a_c, AraC:a_c}
    timepoints = np.linspace(0, 100, 100)
    R = CRN.simulate_with_bioscrape_via_sbml( timepoints, initial_condition_dict = x0)
    plt.plot(R['time'], R['YFP'], label = 'Initial [Ara], [AraC] = ' + str(a_c) )

plt.ylabel(f'[{YFP}]')
plt.xlabel('Time')
plt.legend()
plt.show()

for mech_type, mech in M.mechanisms.items():
    print(f"{mech_type}: {mech}")
