"""Calculations for Example 12.7.5 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
V = 10 # gal
Ve = 1.25 # gal
CB0 = 3 # mol/gal
T0 = 20 + 273.15 # K
VFR = 0.5 # gal/min
CAin = 7.5 # mol/gal
CBin = 3 # mol/gal
Tin = 50 + 273.15 # K
Tein = 20 + 273.15 # K
Te0 = 20 + 273.15 # K
mein = 250 # g/min
U = 190 # cal/ft^2/min/K
A = 4 # ft^2
rhoe = 1 # g/cc
Cpe = 1 # cal/g/K
Cp = 1600 # cal/gal/K
tf = 30 # min
k01 = 4.3E8 # gal/mol/min
E1 = 14200 # cal/mol
dH1 = -11000 # cal/mol
k02 = 2.7E8 # gal/mol/min
E2 = 16100 # cal/mol
dH2 = -11600 # cal/mol
k03 = 3.9E8 # gal/mol/min
E3 = 14800 # cal/mol
dH3 = -12100 # cal/mol
# known
R = 1.987 # cal/mol/K
# calculated
nAin = CAin*VFR
nBin = CBin*VFR

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    nW = dep[2]
    nX = dep[3]
    nY = dep[4]
    nZ = dep[5]
    T = dep[6]
    Te = dep[7]

	# calculate the rate and rate of heat transfer
    r1 = k01*np.exp(-E1/R/T)*nA/VFR*nB/VFR
    r2 = k02*np.exp(-E2/R/T)*nA/VFR*nB/VFR
    r3 = k03*np.exp(-E3/R/T)*nA/VFR*nB/VFR
    Q = U*A*(Te-T)
    
	# evaluate the derivatives
    dnAdt = VFR/V*(nAin - nA - V*(r1 + r2 + r3))
    dnBdt = VFR/V*(nBin - nB - V*(r1 + r2 + r3))
    dnWdt = VFR/V*(-nW + r1*V)
    dnXdt = VFR/V*(-nX + r2*V)
    dnYdt = VFR/V*(-nY + r2*V)
    dnZdt = VFR/V*(-nZ + V*(r1 + r2 + r3))
    dTdt = (Q - VFR*Cp*(T-Tin) - V*(r1*dH1 + r2*dH2 + r3*dH3))/V/Cp
    dTedt = (-Q-mein*Cpe*(Te-Tein))/(rhoe*Ve*Cpe)

	# return the derivatives
    return dnAdt, dnBdt, dnWdt, dnXdt, dnYdt, dnZdt, dTdt, dTedt

# reactor model
def profiles():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([0.0, VFR*CB0, 0.0, 0.0, 0.0, 0.0, T0, Te0])

	# define the stopping criterion
    f_var = 0
    f_val = tf
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    # return all profiles
    return t, dep[0,:], dep[1,:], dep[2,:], dep[3,:], dep[4,:], dep[5,:]\
            , dep[6,:], dep[7,:]

# perform the analysis
def perform_the_analysis():

    # solve the reactor design equations
    [t, nA, nB, nW, nX, nY, nZ, T, Tex] = profiles()

    # calculate the other quantities of interest
    CB = nB/VFR
    T = T - 273.15
    Tex = Tex - 273.15

    # display and save the graphs
    plt.figure(1) # reactor temperature vs. time
    plt.plot(t,T)
    plt.xlabel("Time (min)")
    plt.ylabel("Reacting Fluid Temperature (°C)")
    plt.savefig('reb_12_7_5/python/T_vs_t.png')
    plt.show()

    # display and save the graphs
    plt.figure(2) # reactor temperature vs. time
    plt.plot(t,Tex)
    plt.xlabel("Time (min)")
    plt.ylabel("Shell Temperature (°C)")
    plt.savefig('reb_12_7_5/python/Te_vs_t.png')
    plt.show()

    # display and save the graphs
    plt.figure(1) # concentration of B vs. time
    plt.plot(t,CB)
    plt.xlabel("Time (min)")
    plt.ylabel("Concentration of B (mol gal$^-1$)")
    plt.savefig('reb_12_7_5/python/CB_vs_t.png')
    plt.show()

    return

if __name__=="__main__":
    perform_the_analysis()
