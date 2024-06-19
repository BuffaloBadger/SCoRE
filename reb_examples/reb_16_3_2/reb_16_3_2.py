"""Calculations for Example 16.3.2 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# constants available to all functions
# given
Vpfr = 4.0 # m^3
nDotTot_0 = 1.25 # mol/s
yA_0 = 0.5
yB_0 = 0.5
T_0 = 300 # K
P = 2.5 # atm
Cp_i = 25.8 # cal/mol/K
dH = -8600 # cal/mol
k_0 = 8.12e2 # /s
E = 9500.0 # cal/mol
UA = 13.6 # cal /K /s
nDotY_0 = 0
nDotZ_0 = 0
# known
Re = 1.987 # cal/mol/K
Rw = 0.08206 * 1e-3 # m^3*atm/mol/K
# calculated
nDotA_0 = yA_0*nDotTot_0
nDotB_0 = yB_0*nDotTot_0

# residuals function
def residuals(guess):
    # make the guess available to all functions
    T_1 = guess[0]
    T_3 = guess[1]
    
    #solve the PFR design equations
    V, nDotA, nDotB, nDotY, nDotZ, T = profiles(T_1)
    nDotA_2 = nDotA[-1]
    nDotB_2 = nDotB[-1]
    nDotY_2 = nDotY[-1]
    nDotZ_2 = nDotZ[-1]
    T_2 = T[-1]

    # calculate the heat transfer rate
    AMTD = 0.5*((T_3-T_0) + (T_2-T_1))
    Qdot = UA*AMTD

    # evaluate the residuals
    epsilon_1 = (nDotA_0 + nDotB_0 + nDotY_0 + nDotZ_0)*Cp_i*(T_1-T_0) \
            - (nDotA_2+nDotB_2+nDotY_2+nDotZ_2)*Cp_i*(T_2-T_3)
    epsilon_2 = (nDotA_0+nDotB_0+nDotY_0+nDotZ_0)*Cp_i*(T_1-T_0) - Qdot

    # return the residual
    return epsilon_1, epsilon_2

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nDotA = dep[0]
    nDotB = dep[1]
    nDotY = dep[2]
    nDotZ = dep[3]
    T = dep[4]

	# calculate the rate
    nDotTot = nDotA + nDotB + nDotY + nDotZ
    CA = nDotA*P/nDotTot/Rw/T
    r = k_0*np.exp(-E/Re/T)*CA

	# evaluate the derivatives
    dnDotAdV = -r
    dnDotBdV = -r
    dnDotYdV = r
    dnDotZdV = r
    dTdV = -r*dH/nDotTot/Cp_i

	# return the derivatives
    return dnDotAdV, dnDotBdV, dnDotYdV, dnDotZdV, dTdV

# pfr model
def profiles(T_1):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nDotA_0, nDotB_0, nDotY_0, nDotZ_0, T_1])

	# define the stopping criterion
    f_var = 0
    f_val = Vpfr
     
	# solve the IVODEs
    V, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nDotA = dep[0,:]
    nDotB = dep[1,:]
    nDotY = dep[2,:]
    nDotZ = dep[3,:]
    T = dep[4,:]

    # return all profiles
    return V, nDotA, nDotB, nDotY, nDotZ, T

# perform the analysis
def perform_the_analysis():
    # solution for PFR only
    V, nDotA, nDotB, nDotY, nDotZ, T = profiles(T_0)
    fA_PFR_only = 100*(nDotA_0 - nDotA[-1])/nDotA_0

    # initial guess for the unknowns for a low conversion steady state
    initial_guess = np.array([T_0 + 1, T_0 + 2])

    # solve the heat exchanger energy balances
    soln = sp.optimize.root(residuals,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The heat exchanger equations were not solved: {soln.message}")

    # extract the result
    T_1 = soln.x[0]

    # solve the PFR design equations
    V, nDotA, nDotB, nDotY, nDotZ, T = profiles(T_1)
    fA_low_conversion = 100*(nDotA_0 - nDotA[-1])/nDotA_0

    # repeat with a guess for a medium conversion steady state
    initial_guess = np.array([T_0 + 40, T_0 + 80])
    soln = sp.optimize.root(residuals,initial_guess)
    if not(soln.success):
        print(f"The heat exchanger equations were not solved: {soln.message}")
    T_1 = soln.x[0]
    V, nDotA, nDotB, nDotY, nDotZ, T = profiles(T_1)
    fA_medium_conversion = 100*(nDotA_0 - nDotA[-1])/nDotA_0

    # repeat with a guess for a high conversion steady state
    initial_guess = np.array([T_0 + 100, T_0 + 200])
    soln = sp.optimize.root(residuals,initial_guess)
    if not(soln.success):
        print(f"The heat exchanger equations were not solved: {soln.message}")
    T_1 = soln.x[0]
    V, nDotA, nDotB, nDotY, nDotZ, T = profiles(T_1)
    fA_high_conversion = 100*(nDotA_0 - nDotA[-1])/nDotA_0

    # tabulate the results
    data = [['Conversion without Thermal Backmixing', fA_PFR_only, '%']
            ,['Low Steady-State Conversion', fA_low_conversion, '%']
            ,['Medium Steady-State Conversion', fA_medium_conversion, '%']
            ,['High Steady-State Conversion', fA_high_conversion, '%']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)

    # save the results
    results_df.to_csv('reb_16_3_2/python/results.csv',index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
