# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 14:59:04 2025

@author: zacha
"""

from biocrnpyler import *
import bioscrape
import numpy as np
import matplotlib.pyplot as plt
import libsbml 

class GCSim:
    
    def __init__(self, CRN):
        self.CRN = CRN
        
    def basicsim(self, x0, timepoints, protein_lst, title = '', xlabel = '', ylabel = ''):
        
        plt.figure()
        
        R = self.CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, safe=True)
        
        for i in range(len(protein_lst)):
            plt.plot(R['time'], R[protein_lst[i]], label = protein_lst[i])
        
        plt.legend()
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
        plt.show()
        
        return R
        
    def heatmap(self, x0, timepoints, max_conc, num_val, input_a, input_b, output_protein, title = '', 
                xlabel = '', ylabel = ''):
        
        plt.figure(figsize = (6, 6))
        array = np.zeros((num_val, num_val))
        
        for a_ind, A_c in enumerate(np.linspace(0, max_conc, num_val)):
            for b_ind, B_c in enumerate(np.linspace(0, max_conc, num_val)):
                x0[input_a] = A_c
                x0[input_b] = B_c
                R = self.CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, safe=True)
                array[a_ind, b_ind] = R[output_protein][len(timepoints)-1]
                
        plt.title(title)
        cb = plt.pcolor(array, cmap='plasma')
        plt.colorbar(cb)
        plt.xlabel(xlabel) #corresponds to b index
        plt.ylabel(ylabel) #corresponds to a index

        plt.xticks(np.arange(.5, num_val+.5, 1), [str(i) for i in np.linspace(0, max_conc, num_val)])
        plt.yticks(np.arange(.5, num_val+.5, 1), [str(i) for i in np.linspace(0, max_conc, num_val)])
            
        plt.show()
    
    def inputswitch(self, x0, period, output_protein, input_a, input_b, input_high, 
                    repressor_a, repressor_b, complex_a, complex_b, title = '', 
                    xlabel = '', ylabel = ''):
        
        plt.figure()
        x_lst = []
        y_lst = []
        i = 0
        
        for a in [0, input_high]:
            for b in [0, input_high]:
                x0[f'{input_a}'] = a
                x0[f'{input_b}'] = b
                # x0[f'{repressor_a}'] = input_high
                # x0[f'{repressor_b}'] = input_high
                timepoints = np.linspace(period*i, period*(i+1), 1000)
                R = self.CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0, safe= True)
                x_lst.extend(R['time'])
                y_lst.extend(R[f'{output_protein}'])
                x0[output_protein] = R[output_protein][len(timepoints)-1]
                # for key in R:
                #     if key not in [input_a.name, input_b.name, 'time', complex_a.name,
                #                    complex_b.name]:
                #         x0[key] = R[key][len(timepoints)-1]
                i += 1
        
        plt.plot(x_lst, y_lst, label = f'{output_protein}')
        plt.axvspan(0, period, color='gray', alpha=0.3, 
                    label=f'{input_a} = 0, {input_b} = 0')
        plt.axvspan(period, period*2, color='gray', alpha=0.01, 
                    label=f'{input_a} = 0, {input_b} = {input_high}')
        plt.axvspan(period*2, period*3, color='gray', alpha=0.3, 
                    label=f'{input_a} = {input_high}, {input_b} = 0')
        plt.axvspan(period*3, period*4, color='gray', alpha=0.01, 
                    label=f'{input_a} = {input_high}, {input_b} = {input_high}')
        plt.legend()
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        
        plt.show()
