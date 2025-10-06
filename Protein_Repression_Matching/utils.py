# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 15:37:22 2025

@author: zacha
"""

from biocrnpyler import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
from joblib import Parallel, delayed


class EQN_fit:
    
    def __init__(self, hill, y_min=None, y_max=None, K=None, n=None, x=None, y=None):
        
        self.hill = hill
        
        if hill:
            self.x = np.linspace(0,10,500)
            self.y = EQN_fit.exp_hill(self.x, y_min, y_max, K, n)
        else:
            self.x = x
            self.y = y             
    
    def exp_hill(x, y_min, y_max, K, n):
        '''
        '''
        
        y = y_min + ( (y_max-y_min)*(K**n) ) / ( (x**n) + (K**n) )
        
        return y

    def crn_odes(t, y, kb, ku, coop, ktx_leak, ktx_free, ktl=0.05, kdil=0.0075):
        '''
        '''
        
        R, D, B, M, A = y 
        
        dRdt = coop * ( (ku * B) - (kb * R**coop * D) )
        dDdt = (ku * B) - (kb * R**coop * D)
        dBdt = (kb * R**coop * D) - (ku * B)
        dMdt = (ktx_free * D) + (ktx_leak * B) - (kdil * M)
        dAdt = (ktl * M) - (kdil * A)
        
        return [dRdt, dDdt, dBdt, dMdt, dAdt]
    
    def crn_odes_PTAC_YFP(t, y, kb_C, ku_C, coop_I, kb, ku, coop, ktx_leak, 
                          ktx_free, ktl=0.05, kdil=0.0075):
        '''
        '''
        
        I, C, R, D, B, M, A = y 
        
        dIdt = coop_I * ( (ku_C * C) - (kb_C * I**coop_I * R**coop_I) )
        dCdt = (kb_C * I**coop_I * R**coop_I) - (ku_C * C)
        dRdt = coop * ( (ku * B) - (kb * R**coop * D) ) + coop_I *( (ku_C * C) - (kb_C * I**coop_I * R**coop_I) )
        dDdt = (ku * B) - (kb * R**coop * D)
        dBdt = (kb * R**coop * D) - (ku * B)
        dMdt = (ktx_free * D) + (ktx_leak * B) - (kdil * M)
        dAdt = (ktl * M) - (kdil * A)
        
        return [dIdt, dCdt, dRdt, dDdt, dBdt, dMdt, dAdt]
    
    def crn_ode_bioCRNpyler_sim(timepoints, kb_C, ku_C, coop_I, kb, ku, coop, ktx_leak, 
                          ktx_free, D, R, I, ktl=0.05, kdil=0.0075):
        # Parameters and Global Mechanisms
        parameters={"cooperativity":coop_I,"kb":100, "ku":10, "ktx":.05, "ktl":ktl, "kdeg":0.001, "kdil":kdil}
        complex_parameters = {'kb':kb_C, 'ku':ku_C}
        component_parameters = {
            
            #Defalt Promoter Binding Parameters. Note the part_id = [promoter_name]_[regulator_name]
            ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'kb'):kb, 
            ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'ku'):ku, 
            ParameterKey(mechanism = 'one_step_cooperative_binding', part_id = None, name = 'cooperativity'):coop, 
            
            #Default Promoter Transcription. Note the part_id = [promoter_name]_[regulator_name]
            ParameterKey(mechanism = 'simple_transcription', part_id = None, name = "ktx"): ktx_leak,
            
            #Leak Parameters for transcription
            #These regulate expression of an unbound promoter
            ParameterKey(mechanism = 'simple_transcription', part_id = 'P_Tac_leak', name = "ktx"): ktx_free
            
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
        
        # Simulation and Plotting  
        x0 = {YFP_construct.get_species():D, LacI:R, IPTG:I}
        Res = CRN.simulate_with_bioscrape_via_sbml(timepoints, initial_condition_dict = x0)
        
        return Res['protein_YFP_degtagged'].iloc[-1]

    def crn_sol_parallel_hill(x, kb, ku, coop, ktx_leak, ktx_free, ktl=0.05, kdil=0.0075):
        '''
        '''
        
        def crn_sol_hill(R):
            '''
            '''
                
            D = 0.1 #Initial DNA conc
            B = 0 # bound DNA
            M = 0 # mRNA
            A = 0 # protein A
            
            y0 = [R, D, B, M, A]
            t_span = (0, 1500)
            
            sol = solve_ivp(EQN_fit.crn_odes, t_span, y0, method='LSODA', 
                            args=(kb, ku, coop, ktx_leak, ktx_free))
            
            return sol.y[4,-1]
        
        y_crn_ode = Parallel(n_jobs=-16)(delayed(crn_sol_hill)(R) for R in x)
        
        return y_crn_ode
    
    def crn_sol_parallel(x, kb, ku, coop, ktx_leak, ktx_free, kb_C, ku_C, coop_I,
                         R, ktl=0.05, kdil=0.0075):
        '''
        '''
        
        # def crn_sol(I):
        #     '''
        #     '''
            
        #     C = 0 
        #     D = 0.1 #Initial DNA conc
        #     B = 0 # bound DNA
        #     M = 0 # mRNA
        #     A = 0 # protein A
            
        #     y0 = [I, C, R, D, B, M, A]
        #     t_span = (0, 1500)
            
        #     sol = solve_ivp(EQN_fit.crn_odes_PTAC_YFP, t_span, y0, method='LSODA', 
        #                     args=(kb_C, ku_C, coop_I, kb, ku, coop, ktx_leak, ktx_free))
            
        #     with open(f'debug_{I}.txt', 'w') as f:
        #         f.write(f"coop_I={coop_I}, coop={coop}\n")
        #         f.write(f"[I, C, R, D, B, M, A]={[I, C, R, D, B, M, A]}\n")
        #         f.write(f"sol.y = {sol.y}\n")
            
        #     return sol.y[6,-1]
        
        # y_crn_ode = Parallel(n_jobs=-16)(delayed(crn_sol)(I) for I in x)
        
        y_crn_ode = []
        timepoints = np.linspace(0, 1500, 500)
        D = 0.1 # DNA conc
        
        for I in x:    
            y_crn_ode.append(EQN_fit.crn_ode_bioCRNpyler_sim(timepoints, kb_C, ku_C, coop_I, kb, ku, coop, 
                                                             ktx_leak, ktx_free, D, R, I, ktl=0.05, kdil=0.0075))
        
        return y_crn_ode

    def loss(self, params):
        '''
        least squares approach
        '''
        
        y = np.array(self.y)
        
        if self.hill:
            kb, ku, coop, ktx_leak, ktx_free = params
            y_crn_ode = EQN_fit.crn_sol_parallel_hill(self.x, kb, ku, coop, ktx_leak, ktx_free)
        else:
            kb, ku, coop, ktx_leak, ktx_free, kb_C, ku_C, coop_I, R = params
            y_crn_ode = EQN_fit.crn_sol_parallel(self.x, kb, ku, coop, ktx_leak, ktx_free, 
                                                      kb_C, ku_C, coop_I, R)
        
        y_crn_ode = np.array(y_crn_ode)
         
        return np.sum((y_crn_ode - y)**2)

    def data_fit_hill(self, title='', xlabel='[Repressor]', ylabel='[Output Protein]'):
        '''
        '''

        # Normal plot
        plt.figure()
        
        plt.plot(self.x, self.y, label='Hill Eqn')
        
        
        # Data Fitting kb, ku, coop, ktx_leak, ktx_free
        
        params_init = [100, 10, 2, 0, 0.05]
        
        res = minimize(self.loss, params_init,
                        bounds=[(0.1, None), (0.0001, None), (1, 10), (0, 0.05), (0, 0.1)])
        self.res = res
        
        y_crn_ode = EQN_fit.crn_sol_parallel_hill(self.x, res.x[0], res.x[1], res.x[2], res.x[3], res.x[4])
        
        #R^2 calc
        ss_res = np.sum((np.array(self.y) - np.array(y_crn_ode))**2)
        ss_tot = np.sum((np.array(self.y) - np.mean(self.y))**2)
        r_squared = 1 - (ss_res/ss_tot)
        
        # Mass Action Plot
        plt.plot(self.x, y_crn_ode, label='Mass Action Eqn')
        plt.title(f'{title} R² = {r_squared:.4f}')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend()
        plt.show()
        
        return res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], r_squared
    
    def data_fit(self, title='', xlabel='[Input Protein]', ylabel='[Output Protein]'):
        '''
        

        Parameters
        ----------
        title : TYPE, optional
            DESCRIPTION. The default is ''.
        xlabel : TYPE, optional
            DESCRIPTION. The default is '[Input Protein]'.
        ylabel : TYPE, optional
            DESCRIPTION. The default is '[Output Protein]'.

        Returns
        -------

        '''
        
        # Plot Experimental Data
        plt.figure()
        plt.loglog(self.x, self.y, marker='o', linestyle='', label = 'Experimental Data')
        
        # Data Fitting kb, ku, coop, ktx_leak, ktx_free, kb_C, ku_C, coop_I, R
        
        params_init = [100, 10, 2, 1e-15, 0.05, 100, 10, 2, 10]
        
        res = minimize(self.loss, params_init,
                        bounds=[(0.1, None), (0.0001, None), (1, 10), (1e-15, 0.05), 
                                (1e-15, 0.1), (0.1, None), (0.0001, None), (1, 10), (0,1000)])
        self.res = res
        
        y_crn_ode = EQN_fit.crn_sol_parallel(self.x, res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], 
                                             res.x[5], res.x[6], res.x[7], res.x[8])
        self.y_crn_ode = y_crn_ode
        
        #R^2 calc
        ss_res = np.sum((np.array(self.y) - np.array(y_crn_ode))**2)
        ss_tot = np.sum((np.array(self.y) - np.mean(self.y))**2)
        r_squared = 1 - (ss_res/ss_tot)
        
        # Mass Action Plot
        plt.loglog(self.x, y_crn_ode, marker='o', linestyle='', label='Mass Action Eqn')
        plt.title(f'{title} R² = {r_squared:.4f}')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend()
        plt.show()
        
        return res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], res.x[5], res.x[6], res.x[7], res.x[8], r_squared

    def save_to_excel(protein, df, filename, startrow, startcol=1):
        '''

        Parameters
        ----------
        protein : TYPE
            DESCRIPTION.
        df : TYPE
            DESCRIPTION.
        filename : TYPE
            DESCRIPTION.
        startrow : TYPE
            DESCRIPTION.
        startcol : TYPE, optional
            DESCRIPTION. The default is 1.

        Returns
        -------
        None.

        '''
        
        if protein == 'AmeR-F1':
            header = True
            startrow -= 1
        else:
            header = False
        
        with pd.ExcelWriter(filename, mode='a', if_sheet_exists='overlay', engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Sheet1', startrow=startrow, startcol=startcol, header=header, 
                        index=False)
