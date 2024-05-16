"""Calculations for Example 13.7.2 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# constants available to all functions
# given
D = 10 # cm
L = 500 # cm
CA_in = 1E-3 # mol /cm^3
CB_in = 1.2E-3 # mol /cm^3
Vdot_in = 75E3 # cm^3 /min
k0_1 = 8.72E8 # cm^3 /mol /min
E_1 = 7200 # cal /mol
dH_1 = -10700 # cal /mol
Cp = 1.0 # cal /g /K
rho = 1.0 # g /cm^3
f_A = 0.95
# known
R = 1.987 # cal /mol /K
# calculated
nA_in = Vdot_in*CA_in
nB_in = Vdot_in*CB_in
nA_out = nA_in*(1-f_A)

# make missing initial value or IVODE constant available to all functions
T_in = float('nan')

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    T = dep[4]

	# calculate the rate
    k_1 = k0_1*np.exp(-E_1/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r_1 = k_1*CA*CB

	# evaluate the derivatives
    dnAdz = -np.pi*D**2/4*r_1
    dnBdz = -np.pi*D**2/4*r_1
    dnYdz = np.pi*D**2/4*r_1
    dnZdz = np.pi*D**2/4*r_1
    dTdz = -np.pi*D**2/4*r_1*dH_1/(Vdot_in*rho*Cp)

	# return the derivatives
    return dnAdz, dnBdz, dnYdz, dnZdz, dTdz

# reactor model
def profiles():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_in, nB_in, 0.0, 0.0, T_in])

	# define the stopping criterion
    f_var = 0
    f_val = L
     
	# solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
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
    return z, nA, nB, nY, nZ, T

# implicit equation for IVODE initial value as residual
def residual(guess):
    # make the guess available to all functions
    global T_in
    T_in = guess[0]

    # solve the reactor design equations
    z, nA, nB, nY, nZ, T = profiles()

    # extract the calculated final value of nA
    nA_f = nA[-1]

    # evaluate the residual
    residual = nA_out - nA_f

    # return the residual
    return residual

# perform the analysis
def perform_the_analysis():
    # allow this function to set T_in
    global T_in

    # set initial guess for T_in
    initial_guess = 25 + 273.15

    # solve the implict equation for [missing constant]
    soln = sp.optimize.root(residual,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The initial temperature was NOT found: {soln.message}")

    # extract the result
    T_in = soln.x[0]

    # solve the reactor design equations
    z, nA, nB, nY, nZ, T = profiles()

    # calculate the other quantities of interest
    T_f = T[-1]

    # tabulate the results
    data =[['Inlet T',T_in - 273.15,'°C'],['Outlet T',T_f - 273.15,'°C']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    # save the results
    result.to_csv('reb_13_7_2/python/results.csv', index=False)
    # display the results
    print(' ')
    print(result)
    return

if __name__=="__main__":
    perform_the_analysis()
