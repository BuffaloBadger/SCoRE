"""Calculations for Example 15.4.1 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# constants available to all functions
# given
V_cstr = 350.0 # L
V_pfr = 350.0 # L
CA_0 = 2.5 # mol /L
Vdot = 100.0 # L /min
T_0 = 38 + 273.15 # K
dH1 = -21500 # cal /mol
dH2 = -24000 # cal /mol
Cp = 1.0E3 # cal /L /K
k0_1 = 1.2E5 # /min
E_1 = 9100 # cal/mol
k0_2 = 2.17E7 # L /mol /min
E_2 = 13400 # cal /mol
# known
R = 1.987 # cal /mol
# calculated
nDotA_0 = Vdot*CA_0;

# make CSTR inlet stream available to all functions
nDotA_in = float('nan')
nDotD_in = float('nan')
nDotU_in = float('nan')
T_in = float('nan')

# residuals function
def residuals(guess):
    # make the guess available to all functions
    nDotA_out = guess[0]
    nDotD_out = guess[1]
    nDotU_out = guess[2]
    T_out = guess[3]

    # calculate the rates
    r1 = k0_1*np.exp(-E_1/R/T_out)*nDotA_out/Vdot
    r2 = k0_2*np.exp(-E_2/R/T_out)*(nDotA_out/Vdot)**2

    # evaluate the residuals
    epsilon_1 = nDotA_in - nDotA_out + (-r1 -r2)*V_cstr
    epsilon_2 = nDotD_in - nDotD_out + r1*V_cstr
    epsilon_3 = nDotU_in - nDotU_out + r2*V_cstr
    epsilon_4 = -Vdot*Cp*(T_out - T_in) - V_cstr*(r1*dH1 + r2*dH2)

    # return the residual
    return epsilon_1, epsilon_2, epsilon_3, epsilon_4

# cstr model
def unknowns(initial_guess):
     
    # solve the CSTR design equations
    soln = sp.optimize.root(residuals,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The CSTR solution was NOT found: {soln.message}")

    # return the solution
    return soln.x

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nDotA = dep[0]
    T = dep[3]

	# calculate the rates
    r1 = k0_1*np.exp(-E_1/R/T)*nDotA/Vdot
    r2 = k0_2*np.exp(-E_2/R/T)*(nDotA/Vdot)**2

	# evaluate the derivatives
    dnAdV = -r1 -r2
    dnDdV = r1
    dnUdV = r2
    dTdV = -(r1*dH1 + r2*dH2)/Vdot/Cp

	# return the derivatives
    return dnAdV, dnDdV, dnUdV, dTdV

# pfr model
def profiles(dep_0):
	# set the initial values
    ind_0 = 0.0

	# define the stopping criterion
    f_var = 0
    f_val = V_pfr
     
	# solve the IVODEs
    V, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nDotA = dep[0,:]
    nDotD = dep[1,:]
    nDotU = dep[2,:]
    T = dep[3,:]

    # return all profiles
    return V, nDotA, nDotD, nDotU, T

# perform the analysis
def perform_the_analysis():
    # allow this function to set [missing constant]
    global nDotA_in, nDotD_in, nDotU_in, T_in

    # case a
    print('Starting case a')
    print(' ')

    # set the CSTR inlet stream
    nDotA_in = nDotA_0
    nDotD_in = 0.0
    nDotU_in = 0.0
    T_in = T_0

    # set initial guess for the CSTR outlet stream
    initial_guess = nDotA_in/4, nDotA_in/4, nDotA_in/4, T_in + 40

    # solve the CSTR design equations
    cstr_soln = sp.optimize.root(residuals,initial_guess)

    # check that the solution is converged
    if not(cstr_soln.success):
        print(f"The initial temperature was NOT found: {cstr_soln.message}")

    # extract the result
    dep_0 = cstr_soln.x

    # solve the PFR design equations
    V, nDotA, nDotD, nDotU, T = profiles(dep_0)

    # calculate and save the quantities of interest
    fA_case_a = 100.0*(nDotA_0 - nDotA[-1])/nDotA_0
    sel_case_a = nDotD[-1]/nDotU[-1]
    T_case_a = T[-1] - 273.15

    # case b
    print('Starting case b')
    print(' ')

    # solve the PFR design equations
    dep_0 = np.array([nDotA_0, 0.0, 0.0, T_0])
    V, nDotA, nDotD, nDotU, T = profiles(dep_0)

    # set the CSTR inlet stream
    nDotA_in = nDotA[-1]
    nDotD_in = nDotD[-1]
    nDotU_in = nDotU[-1]
    T_in = T[-1]

    # set initial guess for the CSTR outlet stream
    initial_guess = nDotA_in/4, nDotA_in/4, nDotA_in/4, T_in + 40

    # solve the CSTR design equations
    cstr_soln = sp.optimize.root(residuals,initial_guess)

    # extract the results
    nDotA_2 = cstr_soln.x[0]
    nDotD_2 = cstr_soln.x[1]
    nDotU_2 = cstr_soln.x[2]
    T_2 = cstr_soln.x[3]

    # calculate the other quantities of interest
    fA_case_b = 100.0*(nDotA_0 - nDotA_2)/nDotA_0
    sel_case_b = nDotD_2/nDotU_2
    T_case_b = T_2 - 273.15

    # tabulate the results
    data = [['case a conversion', fA_case_a, '%']
            ,['case a selectivity', sel_case_a, 'mol D per mol U']
            ,['T case a', T_case_a, '°C']
            ,['case b conversion', fA_case_b, '%']
            ,['case b selectivity', sel_case_b, 'mol D per mol U']
            ,['T case b', T_case_b, '°C']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('reb_15_4_1/python/results.csv',index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
