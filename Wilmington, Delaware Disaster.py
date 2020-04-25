# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 23:44:18 2020

@author: Vuyo-Minenhle.Hlongw
"""

#%% Imports
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import newton_krylov as fsolve
exec(open('quickplot.py').read())


#%% Reactions & stoichiometry

# NH4NO3 -> N2O+2H2O 
#A -> B + 2C
stoic_table=np.array([-1,1,2])
Molar_mass=np.array([80,44,18])



#%% System Constants

F_0=310 #lb/h
T0=200 #F
Amonium_feed=0.83
theta=np.array([1,0/Amonium_feed,(1-Amonium_feed)/Amonium_feed])
M=500 #lb
X=0.9999
dHrxn_0=-336 #Btu/lb @959.67 R @500 F

H_H2O_l=168 #Btu/lb
H_H2O_g=1202 #Btu/lb
Cp_NH4NO3=0.38 #Btu/lb
Cp_H2O=0.47 #Btu/lb



rate_constants=pd.DataFrame({
        'Temperature':[969.67, 1019.67], #Temperature in R
        'rate':[0.307, 2.912]      #rate constants in /h
        })
R=1.986 #Btu/lbmol/R


#find the activation energy of the reaction
E=(math.log(rate_constants['rate'][0]/rate_constants['rate'][1])*
(R*rate_constants['Temperature'][0]*rate_constants['Temperature'][1]/
(rate_constants['Temperature'][1]-rate_constants['Temperature'][0]))) #Btu/lbmol



#use standard heats of formation to find CP_N2O
dHF_0=np.array([-339.4,81.55,-241.8])*430.6#Btu/lbmol
dHF_0=(dHF_0/Molar_mass)#Btu/lb

#standard heat of reaction at T=77 F:
dHrxn_T25=dHF_0@stoic_table

Cp_N2O=(dHrxn_0-dHrxn_T25)/(500-77)+(Cp_NH4NO3)-(2*Cp_H2O) #Btu/lb





#%% Reaction system

#Find the temperature
dCp=stoic_table@np.array([Cp_NH4NO3,Cp_N2O,Cp_H2O])
Cp=np.array([Cp_NH4NO3,Cp_N2O,Cp_H2O])
T=(-dHrxn_0*X-theta[2]*(H_H2O_g-H_H2O_l)+theta@Cp*T0+dCp*X*500)/(theta@Cp+dCp*X)


k=rate_constants['rate'][0]*math.exp((E/R)*(1/(T+460)-1/rate_constants['Temperature'][0]))
mass=(F_0*X*Amonium_feed)/(k*M)



tau=1

print(tau)

def f(T):
    k=rate_constants['rate'][0]*math.exp((E/R)*(1/(T+460)-1/rate_constants['Temperature'][0]))
    X_MB=k*tau/(1+k*tau)
    X_EB=(theta@Cp*(T-T0)+theta[2]*(H_H2O_g-H_H2O_l))/-(dHrxn_0+dCp*(T-500))
    G_T=-(dHrxn_0)*X_MB
    R_T=Cp@theta*(T-T0)     
    
    return [X_MB,G_T,R_T]



H=lambda T:f(T)[1]-f(T)[2]

T_i = np.linspace(T0-100,1000)


fig_1 = plt.figure()
chart=plt.subplot(111)
chart.plot(T_i,np.array(list(map(f,T_i)))[:,1],label='Heat generated G(T)')
chart.plot(T_i,np.array(list(map(f,T_i)))[:,2],label='Heat removed R(T)')

chart.annotate('unstable steady-state', xy=(532.77,f(532.77)[2]), xytext=(150,f(532.77)[2]+100),
               arrowprops=dict(facecolor='black', shrink=0.05))
chart.annotate('Upper steady-state',xy=(fsolve(H,1000),f(fsolve(H,1000))[2]), xytext=(fsolve(H,1000)-350,f(fsolve(H,1000))[2]+50),
               arrowprops=dict(facecolor='black',shrink=0.05))
chart.annotate('Lower steady-state',xy=(T0,f(T0)[2]), xytext=(T0-50,f(T0)[2]+100),
               arrowprops=dict(facecolor='black',shrink=0.05))

chart.plot(T,f(T)[1],'o',label='Heat at operating T')
chart.plot(T0,f(T0)[2],'o',label='Heat at T0')

chart.set_xlabel('Temperature [F]')
chart.set_ylabel('Heat [Btu/lb]')
chart.set_ylim(-100,500)
plt.title('Fig.1: Reaction system prior to feed cut off')
chart.legend(loc=4)
plt.show()


tau_i=list(np.linspace(1,10,4))
G_tau=[]
for item in tau_i:
    tau=item
    G_tau.append(np.array(list(map(f,T_i)))[:,1])


fig2=plt.figure()
chart2=plt.subplot(111)

chart2.plot(T_i,np.array(list(map(f,T_i)))[:,2],label='Heat removed R(T)')
chart2.plot(T_i,G_tau[0],label='G(T) tau='+str(tau_i[0]))
chart2.plot(T_i,G_tau[1],label='G(T) tau='+str(tau_i[1]))
chart2.plot(T_i,G_tau[2],label='G(T) tau='+str(tau_i[2]))
chart2.plot(T_i,G_tau[3],label='G(T) tau='+str(tau_i[3]))
chart2.set_xlabel('Temperature [F]')
chart2.set_ylabel('Heat [Btu/lb]')
chart2.annotate('Ignition Temperature',xy=(fsolve(H,T_i[21]),f(fsolve(H,T_i[21]))[1]), 
                xytext=(fsolve(H,T_i[21])+100,f(fsolve(H,T_i[21])+100)[2]), arrowprops=dict(facecolor='black', shrink=0.05))
plt.title('Fig.2: Reaction system as feed flow is shut-off')

chart2.legend(loc=4)
plt.show()

    
    
    
    
    







