"""Calculations for Reaction Engineering Basics Example 20.5.1"""

#import libraries
import pandas as pd
import numpy as np
import scipy as sp
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
V = 0.1 # L

# globally available variables
Vdot_current = float('nan')
CA_0_current = float('nan')
CB_0_current = float('nan')
CY_0_current = float('nan')
CZ_0_current = float('nan')
k_current = float('nan')

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nA_1 = guess[0]
    nB_1 = guess[1]
    nY_1 = guess[2]
    nZ_1 = guess[3]

    # calculate the other unknown quantities
    nA_0 = CA_0_current*Vdot_current
    nB_0 = CB_0_current*Vdot_current
    nY_0 = CY_0_current*Vdot_current
    nZ_0 = CZ_0_current*Vdot_current
    CA_1 = nA_1/Vdot_current
    CB_1 = nB_1/Vdot_current
    r = k_current*CA_1*CB_1

    # evaluate the residuals
    epsilon_1 = nA_0 - nA_1 - V*r
    epsilon_2 = nB_0 - nB_1 - V*r
    epsilon_3 = nY_0 - nY_1 + V*r
    epsilon_4 = nZ_0 - nZ_1 + V*r

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4])

# CSTR model function
def unknowns():
    # guess the solution
    initial_guess = np.array([0.25, 0.25, 0.25, 0.25])
     
	# solve the ATEs
    soln = sp.optimize.root(residuals,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The ATE solver did not converge: {soln.message}")

    # extract the unknowns
    nA_1 = soln.x[0]
    nB_1 = soln.x[1]
    nY_1 = soln.x[2]
    nZ_1 = soln.x[3]

    return nA_1, nB_1, nY_1, nZ_1

# predicted responses function
def predicted_responses(adj_inputs, k_log_10):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    CY_1_model = np.zeros(nExpts)

    # make the rate coefficient globally available
    global k_current
    k_current = 10.**k_log_10

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # make the adjusted input variables globally available
        global Vdot_current, CA_0_current, CB_0_current, CY_0_current \
            , CZ_0_current
        Vdot_current = input[0]
        CA_0_current = input[1]
        CB_0_current = input[2]
        CY_0_current = input[3]
        CZ_0_current = input[4]

        # solve the reactor design equations
        nA_1, nB_1, nY_1, nZ_1 = unknowns()
        
        # calculate the model-predicted response
        CY_1_model[i] = nY_1/Vdot_current

    # return the responses
    return CY_1_model

# function that performs the calculations
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_20_5_1/reb_20_5_1_data.csv')
        # columns: Vdot, CA_0, CB_0, CY_0, CZ_0, CY_1

    # extract the data as arrays
    global Vdot, CA_0, CB_0, CY_0, CZ_0, CY_1
    Vdot = (df['Vdot'].to_numpy())*1.0E-3
    CA_0 = df['CA_0'].to_numpy()
    CB_0 = df['CB_0'].to_numpy()
    CY_0 = df['CY_0'].to_numpy()
    CZ_0 = df['CZ_0'].to_numpy()
    CY_1 = df['CY_1'].to_numpy()

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([Vdot, CA_0, CB_0, CY_0, CZ_0]))

    # make a guess for log_10 of k
    par_guess = [0.0]

    # estimate log_10 of k
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs 
        , CY_1, predicted_responses, False)
    
    # extract the results
    k = 10.**beta[0]
    k_ll = 10.**beta_ci[0,0]
    k_ul = 10.**beta_ci[0,1]
        
    # calculate the model-predicted y and the residuals
    CY_1_model = predicted_responses(adjusted_inputs, beta[0])
    residual = CY_1 - CY_1_model
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(CY_1, CY_1_model, color = 'k', marker='o', ls='')
    plt.plot([min(CY_1),max(CY_1)],[min(CY_1),max(CY_1)], color = 'r', ls = '-')
    plt.xlabel("$C_{Y,1}$ (M)")
    plt.ylabel("$C_{Y,1, model}$ (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(Vdot*1.0E3, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\dot{V}$ (cm^3^ min^-1^)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_Vdot_residual.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(CA_0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{A,0}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_CA0_residual.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(CB_0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{B,0}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_CB0_residual.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(CY_0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Y,0}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_CY0_residual.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(CZ_0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Z,0}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_CZ0_residual.png')
    plt.show()

    # save the results to a .csv file
    data = [['k', f'{k:.3g}', 'L mol^-1^ min^-1^'],
        ['k_lower_limit', f'{k_ll:.3g}', 'L mol^-1^ min^-1^'],
        ['k_upper_limit', f'{k_ul:.3g}', 'L mol^-1^ min^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    result.to_csv("reb_20_5_1/python/reb_20_5_1_results.csv", index=False)

if __name__=="__main__":
    perform_the_calculations()