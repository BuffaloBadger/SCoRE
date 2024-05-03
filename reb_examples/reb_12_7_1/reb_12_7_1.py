"""Calculations for Example 12.7.1 of Reaction Engineering Basics"""

# import libraries
import math
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
k0_1 = 10.2 # gal /mol /min
k0_2 = 17.0 # gal /mol /min
E_1 = 15300 # J /mol
E_2 = 23700 # J /mol
CA_in = 10 # mol /gal
CB_in = 12 # mol /gal
T_in = 350 # K
V = 25 # gal
Vdot_in = 12.5 # gal /min
dH_1_298 = -12000 # J /mol
dH_2_298 = -21300 # J /mol
Cp_A = 85 # J /mol /K
Cp_B = 125 # J /mol /K
Cp_D = 200 # J /mol /K
Cp_U = 170 # J /mol /K
# known
R = 8.314 # J /mol /K
# calculated
nA_in = CA_in*Vdot_in
nB_in = CB_in*Vdot_in

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
    nD = guess[2]
    nU = guess[3]
    T = guess[4]

    # calculate the rates
    k_1 = k0_1*np.exp(-E_1/R/T)
    k_2 = k0_2*np.exp(-E_2/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r_1 = k_1*CA*CB
    r_2 = k_2*CA*CB

    # heats of reaction
    dH_1 = dH_1_298 + (Cp_D - Cp_A - Cp_B)*(T - 298)
    dH_2 = dH_2_298 + (Cp_U - Cp_A - Cp_B)*(T - 298)

    # evaluate the residuals
    residual_1 = nA_in - nA - V*(r_1 + r_2)
    residual_2 = nB_in - nB - V*(r_1 + r_2)
    residual_3 = -nD + V*r_1
    residual_4 = -nU + V*r_2
    residual_5 = (nA_in*Cp_A + nB_in*Cp_B)*(T-T_in) + V*(r_1*dH_1 + r_2*dH_2)

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4 
                     , residual_5])

# perform the analysis
def perform_the_analysis():
	# set the initial guess
    init_guess = np.array([0.5*nA_in, 0.5*nA_in, 0.1*nA_in, 0.1*nA_in
            , T_in + 10.0])
    
    # solve the reactor design equations
    solution = unknowns(init_guess)

    # extract individual unknowns
    nA = solution[0]
    nD = solution[2]
    nU = solution[3]
    T = solution[4]

    # calculate the other quantities of interest
    fA = 100*(nA_in - nA)/nA_in
    S_D_U = nD/nU

    # tabulate the results
    data = [["Conversion",f"{fA}","%"]
            ,["Selectivity",f"{S_D_U}","mol D per mol U"]
            ,["Temperature",f"{T}","K"]]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('reb_12_7_1/python/results.csv',index = False)
    # display and save the graphs
    return

if __name__=="__main__":
    perform_the_analysis()
