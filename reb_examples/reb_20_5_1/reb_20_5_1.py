"""Calculations for Reaction Engineering Basics Example 20.5.1"""

#import libraries
import pandas as pd
import numpy as np
import scipy as sp
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
V = 0.1 # L
R = 1.987E-3 # kcal/mol/K

# globally available variables
# full set of measured responses
gMeasResp = float('nan')
# adjusted inputs for the current experiment
gT = float('nan')
gVdot = float('nan')
gCA_0 = float('nan')
gCB_0 = float('nan')
gCY_0 = float('nan')
gCZ_0 = float('nan')
# current rate expression parameters
gk0 = float('nan')
gE = float('nan')

# residuals function
def residuals(guess):
    # extract the guesses
    nA_1 = guess[0]
    nB_1 = guess[1]
    nY_1 = guess[2]
    nZ_1 = guess[3]

    # calculate the other unknown quantities
    nA_0 = gCA_0*gVdot
    nB_0 = gCB_0*gVdot
    nY_0 = gCY_0*gVdot
    nZ_0 = gCZ_0*gVdot
    CA_1 = nA_1/gVdot
    CB_1 = nB_1/gVdot
    k = gk0*np.exp(-gE/R/gT)
    r = k*CA_1*CB_1

    # evaluate the residuals
    epsilon_1 = nA_0 - nA_1 - V*r
    epsilon_2 = nB_0 - nB_1 - V*r
    epsilon_3 = nY_0 - nY_1 + V*r
    epsilon_4 = nZ_0 - nZ_1 + V*r

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4])

# CSTR model function
def unknowns(T, Vdot, CA0, CB0, CY0, CZ0, CY1, k0, E):
    # make adjusted inputs and rate expression parameters globally available
    global gT, gVdot, gCA_0, gCB_0, gCY_0, gCZ_0, gk0, gE
    gT = T
    gVdot = Vdot
    gCA_0 = CA0
    gCB_0 = CB0
    gCY_0 = CY0
    gCZ_0 = CZ0
    gk0 = k0
    gE = E

    # guess the solution
    nA_0 = CA0*Vdot
    nB_0 = CB0*Vdot
    nY_0 = CY0*Vdot
    nZ_0 = CZ0*Vdot
    nYguess = CY1*Vdot
    extent = nYguess - nY_0
    nAguess = nA_0 - extent
    nBguess = nB_0 - extent
    nZguess = nZ_0 + extent
    initial_guess = np.array([nAguess, nBguess, nYguess, nZguess])
     
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
def predicted_responses(adj_inputs, beta, E):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    CY_1_model = np.zeros(nExpts)

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # solve the reactor design equations
        nA_1, nB_1, nY_1, nZ_1 = unknowns(input[0],input[1],input[2],input[3]
                                          ,input[4],input[5],gMeasResp[i]
                                          ,10**beta,E)
        
        # calculate the model-predicted response
        CY_1_model[i] = nY_1/input[1]

    # return the responses
    return CY_1_model

# quantities of interest function
def quantities_of_interest(adj_inputs):
    # guess the parameters
    par_guess = [6.0, 9.0]

    # estimate the parameters
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adj_inputs 
        , gMeasResp, predicted_responses, False)
    
    # extract the results
    k0 = 10**beta[0]
    k0_CI = 10**beta_ci[0,:]
    E = beta[1]
    E_CI = beta_ci[1,:]
        
    # calculate the model-predicted response and the residuals
    CY_1_model = predicted_responses(adj_inputs, beta[0], beta[1])
    epsilon_expt = gMeasResp - CY_1_model

    return k0, k0_CI, E, E_CI, r_squared, CY_1_model, epsilon_expt

# master function
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_20_5_1/reb_20_5_1_data.csv')
        # columns: T, Vdot, CA_0, CB_0, CY_0, CZ_0, CY_1

    # extract the data as arrays
    T = df['T'].to_numpy()
    Vdot = (df['Vdot'].to_numpy())*1.0E-3
    CA_0 = df['CA_0'].to_numpy()
    CB_0 = df['CB_0'].to_numpy()
    CY_0 = df['CY_0'].to_numpy()
    CZ_0 = df['CZ_0'].to_numpy()
    CY_1 = df['CY_1'].to_numpy()
    global gMeasResp
    gMeasResp = df['CY_1'].to_numpy()

    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([T, Vdot, CA_0, CB_0, CY_0, CZ_0]))

    # calculate the quantities of interest
    k0, k0_CI, E, E_CI, r_squared, CY_1_model, epsilon_expt \
        = quantities_of_interest(adj_inputs)
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(gMeasResp, CY_1_model, color = 'k', marker='o', ls='', label='Data')
    plt.plot([min(gMeasResp),max(gMeasResp)],[min(gMeasResp),max(gMeasResp)]
             , color = 'r', ls = '-', label = 'Parity Line')
    plt.xlabel("$C_{Y,1}$ (M)")
    plt.ylabel("$C_{Y,1, model}$ (M)")
    plt.legend
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(Vdot*1.0E3, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\dot{V}$ (cm^3^ min^-1^)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_Vdot_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(T, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("T (K)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_T_residual.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(CA_0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{A,0}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_CA0_residual.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(CB_0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{B,0}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_CB0_residual.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(CY_0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Y,0}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_CY0_residual.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(CZ_0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Z,0}$ (M)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_1/python/reb_20_5_1_CZ0_residual.png')
    plt.show()

    # save the results to a .csv file
    data = [['k0', f'{k0:.3g}', 'L mol^-1^ min^-1^'],
        ['k_lower_limit', f'{k0_CI[0]:.3g}', 'L mol^-1^ min^-1^'],
        ['k_upper_limit', f'{k0_CI[1]:.3g}', 'L mol^-1^ min^-1^'],
        ['E', f'{E:.3g}', 'kcal mol^-1^'],
        ['E_lower_limit', f'{E_CI[0]:.3g}', 'kcal mol^-1^'],
        ['E_upper_limit', f'{E_CI[1]:.3g}', 'kcal mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    print(" ")
    print(result)
    result.to_csv("reb_20_5_1/python/reb_20_5_1_results.csv", index=False)

if __name__=="__main__":
    perform_the_calculations()