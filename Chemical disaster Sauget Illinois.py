# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 12:35:17 2020

@author: Vuyo-Minenhle.Hlongw
"""

import numpy as np
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt




#%%
#Reaction of interest

# A + 2B -> C + D


#%%


tStop = 121.
tInc = 0.5
t = np.arange(0., tStop, tInc)#time vector running from 0 to 120

k_0=0.00017 #m^3/kmol.min @ 188 deg C
v0=np.array([5.119, 3.26]) # inital volume one for each scenario
FA_0=np.array([9.044, 3.17]) #kmol
FB_0=np.array([33, 43]) #kmol
FW_0=np.array([103.7,103.6]) #kmol
F_0=np.vstack([FA_0,FB_0,FW_0])

Ta_W=25 # deg C

dH_rxn=-5.9e5 #heat of reaction for the production of C
Ea=11273 #cal/mol #activation energy of reaction
Cp_A=40 #cal/mol.K 
Cp_W=18 #cal/mol/K
Cp_B=8.38 #cal/mol.K
Cp=np.vstack([Cp_A,Cp_B,Cp_W])
UA=35.85 #cal/mol.K

def f(y,t):
    T,X=y
    
    k=k_0*(math.exp((Ea/1.987)*(1/461-1/T)))
    r=-k*((F_0[0,0]/v0[0])**2)*(1-X)*(F_0[1,0]/F_0[0,0]-2*X)
    
    dxdt=(r*v0[0])/F_0[0,0]
    dTdt=(-UA*(Ta_W-T)+r*v0[0]*dH_rxn)/(F_0[:,0]*1000@Cp[:,0])

    return [dxdt,dTdt]


y0=[467.992, 0.0423866]

soln = odeint(f, y0, t)
FA=FA_0[0]-FA_0[0]*soln[:,1]
FB=FA_0[0]-2*FA_0[0]*soln[:,1]
FC=FA_0[0]*soln[:,1]
FD=FA_0[0]*soln[:,1]
FW=FW_0[0]-0*soln[:,1]

F=np.column_stack((FA,FB,FC,FD))


#stop


fig_1, charts = plt.subplots(2,2,figsize=(20,10))

charts[0,0].plot(t,soln[:,0])
charts[0,0].set_xlabel('Time [min]')
charts[0,0].set_ylabel('Temperature [F]')
charts[0,1].plot(t,soln[:,1])
charts[0,1].set_xlabel('Time [min]')
charts[0,1].set_ylabel('Conversion []')
charts[1,0].plot(t,F)









