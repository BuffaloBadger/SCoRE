"""Calculations for Example 10.5.1 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import pandas as pd
from score_utils import solve_ivodes
import matplotlib.pyplot as plt

# given, known, and calculated constants available to all functions
# given
CA_0 = 15.0 # mol/L
T_0 = 20 + 273.15 # K
CB_in = 15.0 # mol/L
T_in = T_0
V_0 = 5.0 # L
Vdot_in = 0.25/60 # L/s
V_B = 5.0 # L
V_ex = 2.5 # L
Tex_in = T_0
Vdot_ex = 2.5/60 # L/min
A = 2150 # cm^2
U = 73/60/60 # cal/cm^2/s/K
Tex_0 = T_0
rho = 1000.0 # g/L
Cp = 1.0 # cal/g/K
dH_1 = -13700.0 # cal/mol
k0_1 = 8.11E12 # L/mol/s
E_1 = 17700.0 # cal/mol
P = 1.0 # atm
# known
Re = 1.987 # cal/mol/K
Rw = 0.08206 # L-atm/mol/K
M_B = 40.0 # g/mol
# calculated
nDotB_in = Vdot_in*CB_in
mDotB_in = nDotB_in*M_B
nA_0 = CA_0*V_0
t_f = V_B/Vdot_in
mDot_ex = Vdot_ex*rho

# derivatives function
def derivatives(t, dep):
	# extract necessary dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    T = dep[4]
    Tex = dep[5]

    # calculate the reacting fluid volume and the rate coefficient
    V = Vdot_in*t + V_0
    k_1 = k0_1*np.exp(-E_1/Re/T)

    # calculate the rate
    CA = nA/V
    CB = nB/V
    r = k_1*CA*CB

    # calculate the rate of heat exchange
    Qdot = U*A*(Tex - T)

    # evaluate the derivatives
    dnAdt = -V*r
    dnBdt = nDotB_in -V*r
    dnSdt = V*r
    dnWdt = V*r
    dTdt = (Qdot - mDotB_in*Cp*(T - T_in) - V*r*dH_1 + P*Vdot_in*Re/Rw) \
        /(rho*V*Cp)
    dTedt = (-Qdot - mDot_ex*Cp*(Tex - Tex_in))/(rho*V_ex*Cp)

    # return the derivatives
    return [dnAdt, dnBdt, dnSdt, dnWdt, dTdt, dTedt]

# reactor model function
def profiles():
    # set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0, 0.0, 0.0, T_0, Tex_0])

    # define the stopping criterion
    f_var = 0
    f_val = t_f
     
    # solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, False)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nS = dep[2,:]
    nW = dep[3,:]
    T = dep[4,:]
    Tex = dep[5,:]

    # return all profiles
    return t, nA, nB, nS, nW, T, Tex

# function that performs the analysis
def perform_the_analysis():
    # solve the reactor design equations
    t, nA, nB, nS, nW, T, Tex = profiles()

    # calculate the other quantities of interest
    V = Vdot_in*t + V_0
    CA = nA/V
    T_C = T - 273.15
    t_min = t/60
    
    # display and save the graphs
    plt.figure() # concentration profiles
    plt.plot(t_min,CA)
    plt.xlabel("$Time \; (min)$")
    plt.ylabel("$Concentration of A \; (mol \; L^{-1})$")
    plt.savefig('reb_10_5_1/python/concentration_profile.png')
    plt.show()

    plt.figure() # temperature profile
    plt.plot(t_min,T_C)
    plt.xlabel("$Time \; (min)$")
    plt.ylabel("$Temperature \; (°C)$")
    plt.savefig('reb_10_5_1/python/temperature_profile.png')
    plt.show()

    return

if __name__=="__main__":
    perform_the_analysis()