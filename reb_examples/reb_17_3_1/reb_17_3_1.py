"""Calculations for Example 17.3.1 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# constants available to all functions
# given
dH_1 = -14000. # cal /mol
k0_1 = 4.2E15 # cm^3 /mol /min
E_1 = 18000. # cal /mol
Cp = 1.3 # cal /cm^3 /K
CA_0 = 2.0E-3 # mol /cm^3
CZ_0 = 0.0 # mol /cm^3
Vdot_0 = 500. # cm^3 /min
T_0 = 300. # K
R_R = 1.3
D = 5. # cm
L = 50. # cm
# known
R = 1.987 # cal /mol /K

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nDot_A = dep[0]
    nDot_Z = dep[1]
    T = dep[2]

	# calculate other unknown quantities
    Vdot_3 = Vdot_0
    Vdot_4 = R_R*Vdot_3
    Vdot = Vdot_3 + Vdot_4
    k_1 = k0_1*np.exp(-E_1/R/T)
    CA = nDot_A/Vdot
    CZ = nDot_Z/Vdot
    r_1 = k_1*CA*CZ

	# evaluate the derivatives
    dnDotAdz = -np.pi*D**2/4*r_1
    dnDotZdz = np.pi*D**2/4*r_1
    dTdz = -np.pi*D**2/4*r_1*dH_1/Vdot/Cp

	# return the derivatives
    return dnDotAdz, dnDotZdz, dTdz

# reactor model
def profiles(dep_0):
	# set the initial values
    ind_0 = 0.0

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
    nDotA = dep[0,:]
    nDotZ = dep[1,:]
    T = dep[2,:]

    # return all profiles
    return z, nDotA, nDotZ, T

# stream mixer residuals function
def residuals(guess):
    # extract the individual guesses
    nDotA_1 = guess[0]
    nDotZ_1 = guess[1]
    T_1 = guess[2]

    # calculate other unknown constants
    nDotA_0 = Vdot_0*CA_0
    nDotZ_0 = Vdot_0*CZ_0
    Vdot_3 = Vdot_0
    Vdot_4 = R_R*Vdot_3

    # solve the reactor design equations
    dep_0 = np.array([nDotA_1, nDotZ_1, T_1])
    z, nDotA, nDotZ, T = profiles(dep_0)

    # calculate the other unknown quantities
    nDotA_2 = nDotA[-1]
    nDotZ_2 = nDotZ[-1]
    T_4 = T[-1]

    # calculate molar recycle flows
    nDotA_4 = R_R/(1 + R_R)*nDotA_2
    nDotZ_4 = R_R/(1 + R_R)*nDotZ_2

    # evaluate the residual
    eps_1 = nDotA_1 - nDotA_0 - nDotA_4
    eps_2 = nDotZ_1 - nDotZ_0 - nDotZ_4
    eps_3 = Vdot_0*Cp*(T_1 - T_0) + Vdot_4*Cp*(T_1 - T_4)

    # return the residualw
    return eps_1, eps_2, eps_3

# stream mixer model
def unknowns(initial_guess):
    # solve the other equipment mole and energy balances
    soln = sp.optimize.root(residuals,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The initial temperature was NOT found: {soln.message}")
    
    # extract the results
    nDotA_1 = soln.x[0]
    nDotZ_1 = soln.x[1]
    T_1 = soln.x[2]

    # return the results
    return nDotA_1, nDotZ_1, T_1

# perform the analysis
def perform_the_analysis():
    # set initial guess for the unknowns
    nDotA_0 = Vdot_0*CA_0
    initial_guess = np.array([1.2*nDotA_0, nDotA_0, T_0 + 10])

    # solve the other equipment mole and energy balances
    nDotA_1, nDotZ_1, T_1 = unknowns(initial_guess)

    # solve the reactor design equations
    dep_0 = np.array([nDotA_1, nDotZ_1, T_1])
    z, nDotA, nDotZ, T = profiles(dep_0)

    # extract the outlet values
    nDotA_2 = nDotA[-1]
    nDotZ_2 = nDotZ[-1]
    T_3 = T[-1]

    # calculate molar recycle flows
    nDotA_4 = R_R/(1 + R_R)*nDotA_2
    nDotZ_4 = R_R/(1 + R_R)*nDotZ_2

    # calculate the product molar flows
    nDotA_3 = nDotA_2 - nDotA_4
    nDotZ_3 = nDotZ_2 - nDotZ_4

    # calculate the product concentrations
    Vdot_3 = Vdot_0
    CA_3 = nDotA_3/Vdot_3*1000
    CZ_3 = nDotZ_3/Vdot_3*1000
    
    # tabulate the results
    data = [['nDotA 1',nDotA_1,'mol min^-1^']
            ,['nDotZ 1',nDotZ_1,'mol min^-1^']
            ,['T 1',T_1,'K']
            ,['CA 3',CA_3,'M']
            ,['CZ 3',CZ_3,'M']
            ,['T 3',T_3,'K']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)

    # save the results
    results_df.to_csv('reb_17_3_1/python/results.csv', index=False)

    return

if __name__=="__main__":
    perform_the_analysis()
