# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:25:32 2020

@author: Vuyo-Minenhle.Hlongw
"""

#%% Imports statments
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd


#%% reactions

#rxn1: NO + 2/3NH3 -> 5/6N2 + H20
#rxn3: 2NO -> N2 + O2
#rxn3: O2 + 1//2N2 -> NO2

#%% Constants and initial conditions

species=['NO','NH3','N2','H2O','O2','NO2']
rxn1=np.array([-1,-2/3,5/6,1,0,0])
rxn2=np.array([-2,0,1,0,1,0])
rxn3=np.array([0,0,-1/2,0,-1,1])

stoic_table_df=pd.DataFrame({
        'Species':species,
        'rxn1':rxn1,
        'rxn2':rxn2,
        'rxn3':rxn3,
        })

stoic_table_df.set_index('Species',inplace=True)

#stoic_table_df.info()

k1=0.43 #(dm^3/mol)^1.5/s
k2=2.7 #2.7 (dm^3/mol).s
k3=1

P0=101500 #Pa
T0=298 #K
R=8.314 
V=np.linspace(0,10) #m^3



#%% solving the mass balance

def f(F,V):
    fNO,fNH3,fN2,fH2O,fO2,fNO2=F
    
    f_total=fNO+fNH3+fN2+fH2O+fO2+fNO2
    CT_0=P0/(R*T0)
    
    r_1=k1*((fNH3/f_total)*CT_0)*((fNO/f_total)*CT_0)**(1.5)
    r_2=k2*((fNO/f_total)*CT_0)**(2)
    r_3=k3*((fN2/f_total)*CT_0)*((fO2/f_total)*CT_0)**(2)
    
    r=np.hstack([r_1,r_2,r_3])
    
    dfdv=stoic_table_df.values@r
    
    return dfdv
    
y0=np.array([1000,300,0,0,0,0])

soln = odeint(f, y0, V)

#%% Plotting the solution
fig_1, charts = plt.subplots(1,2,figsize=(20,5))   
charts[0].plot(V,soln)
charts[0].set_xlabel('Volume [m^3]')
charts[0].set_ylabel('Flowrates [mol/s]')
charts[0].legend(species)
charts[1].plot(V,(soln[0,0]-soln[:,0])/soln[0,0])
charts[1].set_xlabel('Volume [m^3]')   
charts[1].set_ylabel('Converion []')

Mole_flows_df=pd.DataFrame({
        'Volume':list(V),
        str('F_'+species[0]):list(soln[:,0]),
        str('F_'+species[1]):list(soln[:,1]),
        str('F_'+species[2]):list(soln[:,2]),
        str('F_'+species[3]):list(soln[:,3]),
        str('F_'+species[4]):list(soln[:,4]),
        str('F_'+species[5]):list(soln[:,5]),
        })    

Mole_flows_df.set_index('Volume',inplace=True)
print(Mole_flows_df)

Mole_flows_df.plot()


    
    
    

    
    
    
    
    
    

