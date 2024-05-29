"""Calculations for Example K.4.1 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
D = 5. # cm
dH = -14000. # cal /mol
Vdot_feed = 500. # cm^3 /min
Cp = 1.3 # cal /cm^3 /K
R_R = 1.3 #
nDotA_feed = 1.0 # mol /min
nDotZ_feed = 0.0
L = 50. # cm
T_feed = 300. # K
k_0 = 4.2E15 # cm^3 /mol /min
E = 18000. # cal /mol
# known
R = 1.987 # cal /mol /K
# calculated
Vdot_prod = Vdot_feed
Vdot_r = R_R*Vdot_prod
Vdot_in = Vdot_feed + Vdot_r

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nDot_A = dep[0]
    nDot_Z = dep[1]
    T = dep[2]

	# calculate the rate
    r = k_0*np.exp(-E/R/T)*nDot_A*nDot_Z/Vdot_in**2

	# evaluate the derivatives
    dnDotAdz = -np.pi*D**2/4*r
    dnDotZdz = np.pi*D**2/4*r
    dTdz = -np.pi*D**2/4*r*dH/Vdot_in/Cp

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

# residuals function
def residuals(guess):
    # extract the individual guesses
    nDotA_in = guess[0]
    nDotZ_in = guess[1]
    T_in = guess[2]
    T_r = guess[3]
    nDotA_prod = guess[4]
    nDotZ_prod = guess[5]
    T_prod = guess[6]

    # solve the reactor design equations
    dep_0 = np.array([nDotA_in, nDotZ_in, T_in])
    z, nDotA, nDotZ, T = profiles(dep_0)

    # extract the calculated final values
    nDotA_out = nDotA[-1]
    nDotZ_out = nDotZ[-1]
    T_out = T[-1]

    # calculate molar recycle flows
    nDotA_r = R_R*nDotA_prod
    nDotZ_r = R_R*nDotZ_prod

    # evaluate the residual
    eps_1 = nDotA_feed + nDotA_r - nDotA_in
    eps_2 = nDotZ_feed + nDotZ_r - nDotZ_in
    eps_3 = Vdot_feed*Cp*(T_in - T_feed) + Vdot_r*Cp*(T_in - T_r)
    eps_4 = nDotA_out - nDotA_r - nDotA_prod
    eps_5 = nDotZ_out - nDotZ_r - nDotZ_prod
    eps_6 = T_out - T_r
    eps_7 = T_out - T_prod

    # return the residualw
    return eps_1, eps_2, eps_3, eps_4, eps_5, eps_6, eps_7

# other equipment model
def unknowns(initial_guess):
    # solve the other equipment mole and energy balances
    soln = sp.optimize.root(residuals,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The initial temperature was NOT found: {soln.message}")
    
    # extract the results
    nDotA_in = soln.x[0]
    nDotZ_in = soln.x[1]
    T_in = soln.x[2]
    T_r = soln.x[3]
    nDotA_prod = soln.x[4]
    nDotZ_prod = soln.x[5]
    T_prod = soln.x[6]

    # return the results
    return nDotA_in, nDotZ_in, T_in, T_r, nDotA_prod, nDotZ_prod, T_prod

# perform the analysis
def perform_the_analysis():
    # set initial guess for the unknowns
    initial_guess = np.array([0.9*nDotA_feed, 0.1*nDotA_feed, T_feed + 10,
            T_feed + 10, 0.1*nDotA_feed, 0.5*nDotA_feed, T_feed + 10])

    # solve the other equipment mole and energy balances
    nDotA_in, nDotZ_in, T_in, T_r, nDotA_prod, nDotZ_prod, T_prod \
        = unknowns(initial_guess)
    
    # tabulate the results
    data = [['A in',nDotA_in,'mol min^-1^']
            ,['Z in',nDotZ_in,'mol min^-1^']
            ,['T in',T_in,'K']
            ,['T r',T_r,'K']
            ,['A prod',nDotA_prod,'mol min^-1^']
            ,['A prod',nDotZ_prod,'mol min^-1^']
            ,['T prod',T_prod,'K']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)

    # save the results
    results_df.to_csv('reb_K_4_1/python/results.csv', index=False)

    return

if __name__=="__main__":
    perform_the_analysis()
