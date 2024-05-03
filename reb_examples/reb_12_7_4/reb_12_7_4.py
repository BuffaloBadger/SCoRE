"""Calculations for Example 12.7.4 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd

# constants available to all functions
# given
V = 500 # cm^3
Vdot_in = 1.0 # cm^3 /s
CA_in = 0.015 # mol /cm^3
CB_in = 0.015 # mol /cm^3
T_in = 50 + 273.15 # K
Cp = 0.35*4.184 # J /g /K
rho = 0.93 # g /cm^3
dH = -20000 # J /mol
k0 = 3.24E12 # cm^3 /mol /s
E = 105000 # J /mol
# known
R = 8.314 # J /mol /K
# calculated
nA_in = Vdot_in*CA_in
nB_in = Vdot_in*CB_in

# reactor model
def unknowns(init_guess):
     
	# solve the ATEs
    soln = sp.optimize.root(residuals,init_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The initial temperature was NOT found: {soln.message}")

    # return the solution
    return soln.x

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nA = guess[0]
    nB = guess[1]
    nY = guess[2]
    nZ = guess[3]
    T = guess[4]

    # calculate the rate
    k = k0*np.exp(-E/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r = k*CA*CB

    # evaluate the residuals
    residual_1 = nA_in - nA - V*r
    residual_2 = nB_in - nB - V*r
    residual_3 = -nY + V*r
    residual_4 = -nZ + V*r
    residual_5 = -Vdot_in*rho*Cp*(T - T_in) - V*r*dH

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4
                     , residual_5])

# perform the analysis
def perform_the_analysis():
	# set an initial guess for the low-conversion steady state
    init_guess = np.array([0.99*nA_in, 0.99*nB_in, 0.01*nA_in, 0.01*nA_in
            ,T_in + 1.])

    # solve the reactor design equations
    solution = unknowns(init_guess)

    # save the conversion and temperature
    fA_low = 100*(nA_in - solution[0])/nA_in
    T_low = solution[4] - 273.15

	# set an initial guess for the high-conversion steady state
    init_guess = np.array([0.01*nA_in, 0.01*nB_in, 0.99*nA_in, 0.99*nA_in
            ,T_in + 200.])

    # solve the reactor design equations
    solution = unknowns(init_guess)

    # save the conversion and temperature
    fA_high = 100*(nA_in - solution[0])/nA_in
    T_high = solution[4] - 273.15

	# set an initial guess for the mid-conversion steady state
    fA_guess = (fA_high + fA_low)/100/2
    nA_guess = nA_in*(1-fA_guess)
    nY_guess = nA_in*fA_guess
    init_guess = np.array([nA_guess, nA_guess, nY_guess, nY_guess
            ,(T_low + T_high)/2 + 273.15])

    # solve the reactor design equations
    solution = unknowns(init_guess)

    # save the conversion and temperature
    fA_mid = 100*(nA_in - solution[0])/nA_in
    T_mid = solution[4] - 273.15

    # calculate the other quantities of interest
    # tabulate the results
    data = [[fA_low, T_low],[fA_mid, T_mid],[fA_high, T_high]]
    results_df = pd.DataFrame(data, columns=['Conversion (%)'
                                             ,'Temperature (°C)'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('reb_12_7_4/python/results.csv',index=False)

    # display and save the graphs
    return

if __name__=="__main__":
    perform_the_analysis()
