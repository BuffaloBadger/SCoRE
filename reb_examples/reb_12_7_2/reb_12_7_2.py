"""Calculations for Example 12.7.2 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd

# constants available to all functions
# given
yA_in = 0.1
yB_in = 0.65
yI_in = 0.25
T_in = 165 + 273.15 # K
P = 5.0 # atm
yA_max = 0.001
tau = 0.5 # min
k0 = 1.37E5 # m^3 /mol /min
E = 11100.0 # cal /mol
dH = -7200.0 # cal /mol
Cp_A = 7.6 # cal /mol /K
Cp_B = 8.2 # cal /mol /K
Cp_I = 4.3 # cal /mol /K
# known
Re = 1.987 # cal /mol /K
Rw = 8.206E-5 # m^3 atm /mol /K
# basis
V = 1.0 # m^3
# calculated
Vdot_in = V/tau
nA_in = yA_in*P*Vdot_in/Rw/T_in
nB_in = yB_in*P*Vdot_in/Rw/T_in
nI_in = yI_in*P*Vdot_in/Rw/T_in

# reactor model
def unknowns(init_guess):
     
	# solve the ATEs
    soln = sp.optimize.root(residuals,init_guess)

    # check that the solution is converged
    if not(soln.success):
        print("")
        print(f"The solver did NOT converge: {soln.message}")

    # return the solution
    return soln.x

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nA = guess[0]
    nB = guess[1]
    nI = guess[2]
    nZ = guess[3]
    T = guess[4]

    # calculate the rate
    k = k0*np.exp(-E/Re/T)
    CA = nA/(nA + nB + nI + nZ)*P/Rw/T
    CB = nB/(nA + nB + nI + nZ)*P/Rw/T
    r = k*CA*CB

    # evaluate the residuals
    residual_1 = nA_in - nA - V*r
    residual_2 = nB_in - nB - V*r
    residual_3 = nI_in - nI
    residual_4 = -nZ + V*r
    residual_5 = -(nA_in*Cp_A + nB_in*Cp_B + nI_in*Cp_I)*(T - T_in) - V*r*dH

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4, residual_5])

# perform the analysis
def perform_the_analysis():
	# set the initial guess
    init_guess = np.array([0.01*nA_in, nB_in - 0.01*nA_in, nI_in, 0.01*nA_in
                           , T_in + 10.0])

    # solve the reactor design equations
    solution = unknowns(init_guess)

    # extract the individual results
    nA = solution[0]
    nB = solution[1]
    nI = solution[2]
    nZ = solution[3]
    T = solution[4] - 273.15

    # calculate the other quantities of interest
    nTot = nA + nB + nI + nZ
    yA = nA/nTot
    yB = nB/nTot
    yI = nI/nTot
    yZ = nZ/nTot

    # tabulate the results
    data = [['yA',f'{100*yA}','%'],['yB',f'{100*yB}','%'],['yI',f'{100*yI}','%']
         ,['yZ',f'{100*yZ}','%'],['T',f'{T}','°C']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('reb_12_7_2/python/results.csv',index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
