"""Calculations for Reaction Engineering Basics Example 20.5.3"""

#import libraries
import pandas as pd
import numpy as np
import scipy as sp
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
V = 3.0 # gal

# globally available variables
CA_1 = float('nan')
Vdot = float('nan')
CA_0 = float('nan')
CY_0 = float('nan')
CZ_0 = float('nan')
i_expt_current = 0
kf_current = float('nan')
kr_current = float('nan')

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nA_1 = guess[0]
    nY_1 = guess[1]
    nZ_1 = guess[2]

    # calculate the other unknown quantities
    nA_0 = CA_0[i_expt_current]*Vdot[i_expt_current]
    nY_0 = CY_0[i_expt_current]*Vdot[i_expt_current]
    nZ_0 = CZ_0[i_expt_current]*Vdot[i_expt_current]
    CA = nA_1/Vdot[i_expt_current]
    CY = nY_1/Vdot[i_expt_current]
    CZ = nZ_1/Vdot[i_expt_current]
    r = kf_current*CA**2 - kr_current*CY*CZ

    # evaluate the residuals
    epsilon_1 = nA_0 - nA_1 - V*r
    epsilon_2 = nY_0 - nY_1 + V*r
    epsilon_3 = nZ_0 - nZ_1 + V*r

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3])

# CSTR model function
def unknowns():
    # guess the solution
    nA_guess = CA_1[i_expt_current]*Vdot[i_expt_current]
    xi = CA_0[i_expt_current]*Vdot[i_expt_current] - nA_guess
    nY_guess = CY_0[i_expt_current]*Vdot[i_expt_current] + xi
    nZ_guess = CZ_0[i_expt_current]*Vdot[i_expt_current] + xi
    initial_guess = [nA_guess, nY_guess, nZ_guess]
     
	# solve the ATEs
    soln = sp.optimize.root(residuals,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The ATE solver did not converge: {soln.message}")

    # extract the unknowns
    nA_1 = soln.x[0]
    nY_1 = soln.x[1]
    nZ_1 = soln.x[2]

    return nA_1, nY_1, nZ_1

# predicted responses function
def predicted_responses(adj_inputs, log_kf_guess, log_kr_guess):
    # make the parameters globally available
    global kf_current, kr_current
    kf_current = 10.**log_kf_guess
    kr_current = 10.**log_kr_guess

    # allocate storage for the responses
    nExpts = len(adj_inputs)
    CA_1_model = np.zeros(nExpts)

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # make the experiment index globally available
        global i_expt_current
        i_expt_current = i

        # solve the reactor design equations
        nA_1, nY_1, nZ_1 = unknowns()
        
        # calculate the model-predicted response
        CA_1_model[i] = nA_1/Vdot[i]

    # return the responses
    return CA_1_model

# function that performs the calculations
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_20_5_3/reb_20_5_3_data.csv')
        # columns: Vdot, CA_0, CY_0, CZ_0, CA_1

    # extract the data as arrays
    global Vdot, CA_0, CY_0, CZ_0, CA_1
    Vdot = df['Vdot'].to_numpy()
    CA_0 = df['CA_0'].to_numpy()
    CY_0 = df['CY_0'].to_numpy()
    CZ_0 = df['CZ_0'].to_numpy()
    CA_1 = df['CA_1'].to_numpy()

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([Vdot, CA_0, CY_0, CZ_0]))

    # make a guess for log_10 of kf and kr
    par_guess = [0.0, 0.0]

    # estimate log_10 of kf and kr
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs 
        , CA_1, predicted_responses, False)
    
    # extract the results
    kf = 10.**beta[0]
    kf_ll = 10.**beta_ci[0,0]
    kf_ul = 10.**beta_ci[0,1]
    kr = 10.**beta[1]
    kr_ll = 10.**beta_ci[1,0]
    kr_ul = 10.**beta_ci[1,1]
        
    # calculate the model-predicted response and the residuals
    CA_1_model = predicted_responses(adjusted_inputs, beta[0], beta[1])
    residual = CA_1 - CA_1_model
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(CA_1, CA_1_model, color = 'k', marker='o', ls='')
    plt.plot([min(CA_1),max(CA_1)],[min(CA_1),max(CA_1)], color = 'r', ls = '-')
    plt.xlabel("$C_{A,1}$ (gal lbmol^-1^ min^-1^)")
    plt.ylabel("$C_{A,1, model}$ (gal lbmol^-1^ min^-1^)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(Vdot, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\dot{V}$ (gal/min)")
    plt.ylabel("Residual (lbmol/gal)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_Vdot_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(CA_0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{A,0}$ (lbmol/gal)")
    plt.ylabel("Residual (lbmol/gal)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_CA0_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(CY_0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Y,0}$ (lbmol/gal)")
    plt.ylabel("Residual (lbmol/gal)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_CY0_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(CZ_0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Z,0}$ (lbmol/gal)")
    plt.ylabel("Residual (lbmol/gal)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_CZ0_residual.png')
    plt.show()

    # save the results to a .csv file
    data = [['kf', f'{kf:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['kf_lower_limit', f'{kf_ll:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['kf_upper_limit', f'{kf_ul:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['kr', f'{kr:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['kr_lower_limit', f'{kr_ll:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['kr_upper_limit', f'{kr_ul:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    print(" ")
    print(result)
    result.to_csv("reb_20_5_3/python/reb_20_5_3_results.csv", index=False)

if __name__=="__main__":
    perform_the_calculations()