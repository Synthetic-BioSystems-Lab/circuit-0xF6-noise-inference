# -*- coding: utf-8 -*-
"""
Created on Thu Sep  4 13:28:43 2025

@author: zacha
"""

import pandas as pd
import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=3, ncols=1, sharey=True)

fig.supxlabel('Time', fontsize=14)
fig.supylabel('Percent Output DNA Bound to Repressors', fontsize=14)

#Lists for Selecting Data and Plotting
state_lst = ['simulation_data/IPTG_0_aTc_0_Ara_0.csv', 'simulation_data/IPTG_0_aTc_80_Ara_0.csv', 
             'simulation_data/IPTG_0_aTc_80_Ara_80.csv']
title_lst = ['000', '010', '011']
label_lst = ['Unbound DNA', 'PhlF Only', 'AmtR Only', 'Both Repressors']
color_lst = ['green', 'purple', 'C0', 'red']
linestyles = ['-', '--', '-.', ':']

for i, state in enumerate(state_lst):
#convert data sets into a DataFrames

    df = pd.read_csv(state)

    #Selecting data to plot
    
    dna_plot_lst = ['time', 'dna_part_P_PhlF_forward_part_P_AmtR_forward_part_UTR1_forward_part_YFP_forward_part_t16_forward_',
                    'ordered_polymer_complex_part_P_PhlF_protein_PhlF_degtagged_4x_forward__part_P_AmtR_forward_part_UTR1_forward_part_YFP_forward_part_t16_forward_',
                    'ordered_polymer_part_P_PhlF_forward_complex_part_P_AmtR_protein_AmtR_degtagged_4x_forward__part_UTR1_forward_part_YFP_forward_part_t16_forward_', 
                    'ordered_polymer_complex_part_P_PhlF_protein_PhlF_degtagged_4x_forward__complex_part_P_AmtR_protein_AmtR_degtagged_4x_forward__part_UTR1_forward_part_YFP_forward_part_t16_forward_']
    
    dna_columns = df[dna_plot_lst]
    
    plot_columns = dna_columns.columns[1:]
    
    for j in range(len(plot_columns)):
    
        ax[i].plot(dna_columns['time'], dna_columns[plot_columns[j]] , 
                   color=color_lst[j], linestyle=linestyles[j], label=label_lst[j])
        
    ax[i].set_xlim(0, 2500)
    ax[i].set_ylim(0, 1.1)
    ax[i].set_title(title_lst[i])
        
        
fig.tight_layout()
#fig.savefig('FL.png', dpi=300, bbox_inches='tight')
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 2.65))
plt.show()


