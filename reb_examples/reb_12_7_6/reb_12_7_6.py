"""Calculations for Example 12.7.6 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
V = 500.0 # cm^3
Vdot_in = 1.0 # cm^3 /s
CA_in = 0.015 # mol /cm^3
CB_in = 0.015 # mol /cm^3
T_in_u = 50 + 273.15 # K
T_in_ign = 91 + 273.15 # K
T_in_ext = 3.0 + 273.15 # K
T_0_u = 140 + 273.15 # K
T_0_ign = 96 + 273.15 # K
T_0_ext = 200 + 273.15 # K
fA_0_u = 0.399
fA_0_ign = 0.032
fA_0_ext = 0.881
Cp = 0.35*4.184 # J /g /K
rho = 0.93 # g /cm^3
dH = -20000.0 # J /mol
k0 = 3.24E12 # cm^3 /mol /s
E = 105000.0 # J /mol
# known
R = 8.314 # J /mol /K
# calculated
nA_in = CA_in*Vdot_in
nB_in = CB_in*Vdot_in

# make T_in, T_0, and fA_0 available to all functions
T_in = -1.0
T_0 = -1.0
fA_0 = -1.0

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    nY = dep[2]
    nZ = dep[3]
    T = dep[4]

	# calculate the rate
    k = k0*np.exp(-E/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r = k*CA*CB
    
	# evaluate the derivatives
    dnAdt = Vdot_in/V*(nA_in - nA - V*r)
    dnBdt = Vdot_in/V*(nB_in - nB - V*r)
    dnYdt = Vdot_in/V*( -nY + V*r)
    dnZdt = Vdot_in/V*( -nZ + V*r)
    dTdt = -(Cp*Vdot_in*rho*(T-T_in) + V*r*dH)/(V*Cp*rho)

	# return the derivatives
    return dnAdt, dnBdt, dnYdt, dnZdt, dTdt

# reactor model
def profiles(t_f):
	# set the initial values
    ind_0 = 0.0
    extent = fA_0*nA_in
    dep_0 = np.array([nA_in - extent, nB_in - extent, extent, extent, T_0])

	# define the stopping criterion
    f_var = 0
    f_val = t_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]
    T = dep[4,:]

    # return all profiles
    return t, nA, nB, nY, nZ, T

# perform the analysis
def perform_the_analysis():
    # allow this function to set T_in, T_0, and fA_0
    global T_in, T_0, fA_0

    # analyze the perturbation from the unsteady state
    T_in = T_in_u
    T_0 = T_0_u
    fA_0 = fA_0_u
    t_f = 500.0 # s
    [t_u, nA, nB, nY, nZ, T_u] = profiles(t_f)

    # analyze the perturbation near the ignition point
    T_in = T_in_ign
    T_0 = T_0_ign
    fA_0 = fA_0_ign
    t_f = 10000.0 # s
    [t_ign, nA, nB, nY, nZ, T_ign] = profiles(t_f)

    # analyze the perturbation near the extinction point
    T_in = T_in_ext
    T_0 = T_0_ext
    fA_0 = fA_0_ext
    t_f = 5000.0 # s
    [t_ext, nA, nB, nY, nZ, T_ext] = profiles(t_f)

    # create, display and save the graphs
    plt.figure(1) # perturbation from unsteady state
    plt.plot(t_u/60.0,T_u - 273.15)
    plt.xlabel("Time (min)")
    plt.ylabel("Reacting Fluid Temperature (K")
    plt.savefig('reb_12_7_6/python/unsteady.png')
    plt.show()

    plt.figure(2) # perturbation near the ignition point
    plt.plot(t_ign/60.0,T_ign - 273.15)
    plt.xlabel("Time (min)")
    plt.ylabel("Reacting Fluid Temperature (K")
    plt.savefig('reb_12_7_6/python/ignition.png')
    plt.show() 

    plt.figure(3) # perturbation near the extinction point
    plt.plot(t_ext/60.0,T_ext - 273.15)
    plt.xlabel("Time (min)")
    plt.ylabel("Reacting Fluid Temperature (K")
    plt.savefig('reb_12_7_6/python/extinction.png')
    plt.show()    

    return

if __name__=="__main__":
    perform_the_analysis()
