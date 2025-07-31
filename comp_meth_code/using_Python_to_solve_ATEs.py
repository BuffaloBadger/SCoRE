""" Calculations for Using Python to Solve ATEs"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd

# global constants
CA_0 = 1.0 # mol /L
V = 50 # L
CB_0 = 1.2 # mol /L
rho = 1.0E3 # g /L
Cp = 1.0 # cal /g /K
T_0 = 303 # K
dH = -10700 # cal /mol
k0 = 8.72E5 # L /mol /min
E = 7200 # cal /mol
R = 1.987 # cal /mol /K

# parameter
Vdot_value = np.array([75,100]) # L /min

# global variables
g_Vdot = float('nan')

# residuals function
def CSTR_residuals(guess):
    # extract the individual guesses
    nDotA_1 = guess[0]
    nDotB_1 = guess[1]
    nDotY_1 = guess[2]
    nDotZ_1 = guess[3]
    T_1 = guess[4]

    # calculate r
    k = k0*np.exp(-E/R/T_1)
    CA = nDotA_1/g_Vdot
    CB = nDotB_1/g_Vdot
    r = k*CA*CB

    # evaluate the residuals
    epsilon_1 = g_Vdot*CA_0 - nDotA_1 - r*V
    epsilon_2 = g_Vdot*CB_0 - nDotB_1 - r*V
    epsilon_3 = - nDotY_1 + r*V
    epsilon_4 = - nDotZ_1 + r*V
    epsilon_5 = rho*g_Vdot*Cp*(T_1 - T_0) + r*V*dH

    # return the residuals as an array
    epsilon = np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4, epsilon_5])
    return epsilon

# CSTR function
def CSTR_variables(Vdot):
    # define guesses for the CSTR variables
    nA1_guess = 0.9*Vdot*CA_0
    nB1_guess = 0.9*Vdot*CB_0
    nY1_guess = 0.0
    nZ1_guess = 0.0
    T1_guess = T_0 + 5.0

    # create an array containing the individual guesses
    guess = np.array([nA1_guess, nB1_guess, nY1_guess, nZ1_guess, T1_guess])

    # make Vdot globally available
    global g_Vdot
    g_Vdot = Vdot

    # solve the ATEs
    soln = sp.optimize.root(CSTR_residuals,guess)

    # check that a solution was found
    if not(soln.success):
        print("")
        print(f"The solver did NOT converge: {soln.message}")

    # extract the CSTR variables to be returned from the solution
    nA = soln.x[0]
    nB = soln.x[1]
    nY = soln.x[2]
    nZ = soln.x[3]
    T = soln.x[4]
    return nA, nB, nY, nZ, T

# deliverables function definition
def deliverables():
    # perform the calculations for each of the parameter values
    for i in range(0, len(Vdot_value)):
        # set Vdot
        Vdot = Vdot_value[i]

        # solve the ATEs
        nA, nB, nY, nZ, T = CSTR_variables(Vdot)

        # tabulate and display the results
        data = [['nA',f'{nA:.1f}','mol/min'],
                ['nB',f'{nB:.1f}','mol/min'],
                ['nY',f'{nY:.1f}','mol/min'],
                ['nZ',f'{nZ:.1f}','mol/min'],
                ['T',f'{T:.1f}','K']]
        results_df = pd.DataFrame(data, columns=['item','value','units'])
        print("")
        print(f"Results for Vdot = {Vdot:.0f} L/min")
        print(results_df)
    return

# execution command
if __name__=="__main__":
    deliverables()