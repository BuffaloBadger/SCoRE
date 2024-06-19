"""Calculations for Example 16.3.1 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# constants available to all functions
# given
Vdot = 750.0 # L /min
CA_0 = 3.8 # mol /L
T_0 = 25 + 273.15 # K
UA = 1500.0 # kJ /min /K
dH = -79.8 # kJ/mol
Cp = 987.0 * 4.184e-3 # kJ /L /K
fA = 0.8
k_0 = 3.38E6 # /min
E = 50.0 # kJ ;mol
# known
R = 8.31446e-3 # kJ /mol /K
# calculated
nDotA_0 = Vdot*CA_0

# residuals function
def residuals(guess):
    # make the guess available to all functions
    T_1 = guess[0]
    T_3 = guess[1]
    
    #solve the PFR design equations and extract T_2
    V, nDotA, nDotZ, T = profiles(T_1)
    T_2 = T[-1]

    # calculate the heat transfer rate
    if T_3 - T_0 == T_2 - T_1:
        LMTD = T_3 - T_0
    else:
        LMTD = ((T_3 - T_0) - (T_2 - T_1))/np.log((T_3 - T_0)/(T_2 - T_1))
    Qdot = UA*LMTD

    # evaluate the residuals
    epsilon_1 = Vdot*Cp*(T_1-T_0) - Vdot*Cp*(T_2-T_3)
    epsilon_2 = Vdot*Cp*(T_1-T_0) - Qdot

    # return the residual
    return epsilon_1, epsilon_2

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nDotA = dep[0]
    T = dep[2]

	# calculate the rate
    r = k_0*np.exp(-E/R/T)*nDotA/Vdot

	# evaluate the derivatives
    dnDotAdV = -r
    dnDotZdV = r
    dTdV = -r*dH/Vdot/Cp

	# return the derivatives
    return dnDotAdV, dnDotZdV, dTdV

# pfr model
def profiles(T_1):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nDotA_0, 0.0, T_1])

	# define the stopping criterion
    f_var = 1
    f_val = nDotA_0*(1-fA)
     
	# solve the IVODEs
    V, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nDotA = dep[0,:]
    nDotZ = dep[1,:]
    T = dep[2,:]

    # return all profiles
    return V, nDotA, nDotZ, T

# perform the analysis
def perform_the_analysis():
    # case a, PFR only
    # solve the PFR design equations
    V, nDotA, nDotZ, T = profiles(T_0)
    V_pfr_only = V[-1]
    T_pfr_only = T[-1] - 273.15

    # case b, thermally backmixed PFR
    # initial guess for the unknowns
    initial_guess = np.array([T_0 + 30, T_0 + 60])

    # solve the heat exchanger energy balances
    soln = sp.optimize.root(residuals,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The heat exchanger equations were not solved: {soln.message}")

    # extract the result
    T_1 = soln.x[0]

    # solve the PFR design equations
    V, nDotA, nDotZ, T = profiles(T_1)
    V_tb_pfr = V[-1]
    T_tb_pfr = T[-1] - 273.15

    # tabulate the results
    data = [['PFR Volume', V_pfr_only, 'L']
            ,['PFR Outlet Temperature', T_pfr_only, '°C']
            ,['Thermally Backmixed PFR Volume', V_tb_pfr, 'L']
            ,['Thermally Backmixed PFR Outlet Temperature', T_tb_pfr, '°C']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)

    # save the results
    results_df.to_csv('reb_16_3_1/python/results.csv',index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
