"""Calculations for Reaction Engineering Basics Example 20.5.1"""

#import libraries
import pandas as pd
import numpy as np
import scipy as sp
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
P = 3.0 # atm
V = 1.0 # L (basis)
R = 0.08206 # L*atm/mol/K
Ren = 8.314E-3 # kJ/mol/K

# globally available variables
# full set of measured responses
gMeasResp = float('nan')
# adjusted inputs for the current experiment
gT = float('nan')
gtau = float('nan')
gyA_0 = float('nan')
# current rate expression parameters
gk0 = float('nan')
gE = float('nan')
galphaA = float('nan')
galphaB = float('nan')

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nA_1 = guess[0]
    nB_1 = guess[1]
    nZ_1 = guess[2]

    # calculate the other unknown quantities
    Vdot_0 = V/gtau
    nA_0 = gyA_0*P*Vdot_0/R/gT
    nB_0 = (1-gyA_0)*P*Vdot_0/R/gT
    PA = nA_1*P/(nA_1+nB_1+nZ_1)
    PB = nB_1*P/(nA_1+nB_1+nZ_1)
    k = gk0*np.exp(-gE/Ren/gT)
    r = k*PA**galphaA*PB**galphaB

    # evaluate the residuals
    epsilon_1 = nA_0 - nA_1 - V*r
    epsilon_2 = nB_0 - nB_1 - V*r
    epsilon_3 = -nZ_1 + V*r

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3])

# CSTR model function
def unknowns(T, tau, yA_0, CZ_1, k0, E, alphaA, alphaB):
    # make adjusted inputs and rate expression parameters globally available
    global gT, gtau, gyA_0, gk0, gE, galphaA, galphaB
    gT =T
    gtau = tau
    gyA_0 = yA_0
    gk0 = k0
    gE = E
    galphaA = alphaA
    galphaB = alphaB

    # guess the solution
    Vdot_0 = V/tau
    nA_0 = yA_0*P*Vdot_0/R/T
    nB_0 = (1-yA_0)*P*Vdot_0/R/T
    nZ_guess = CZ_1*Vdot_0
    nA_guess = nA_0 - nZ_guess
    nB_guess = nB_0 - nZ_guess
    initial_guess = [nA_guess, nB_guess, nZ_guess]
     
	# solve the ATEs
    soln = sp.optimize.root(residuals,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The ATE solver did not converge: {soln.message}")

    # extract the unknowns
    nA_1 = soln.x[0]
    nB_1 = soln.x[1]
    nZ_1 = soln.x[2]

    return nA_1, nB_1, nZ_1

# predicted responses function
def predicted_responses(adj_inputs, beta, E, alphaA, alphaB):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    CZ_1_model = np.zeros(nExpts)

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # solve the reactor design equations
        nA_1, nB_1, nZ_1 = unknowns(input[0], input[1], input[2]
                                    , gMeasResp[i], 10**beta, E, alphaA, alphaB)
        
        # calculate the model-predicted response
        CZ_1_model[i] = nZ_1*P/(nA_1 + nB_1 + nZ_1)/R/input[0]

    # return the responses
    return CZ_1_model

# quantities of interest function
def quantities_of_interest(adj_inputs):
    # guess the parameters
    #par_guess = [0.0, 75.0, 1.0, 1.0] # original guess
    par_guess = [np.log10(7260), 114, 1.4, 0.6]

    # estimate the parameters
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adj_inputs 
        , gMeasResp, predicted_responses, False)

    # extract the results
    k0 = 10**beta[0]
    k0_CI = 10**beta_ci[0,:]
    E = beta[1]
    E_CI = beta_ci[1,:]
    alphaA = beta[2]
    alphaA_CI = beta_ci[2,:]
    alphaB = beta[3]
    alphaB_CI = beta_ci[3,:]

    # calculate the model-predicted responses and the experiment residuals
    CZ_1_model = predicted_responses(adj_inputs, beta[0], beta[1], beta[2]
                                     , beta[3])
    epsilon_expt = gMeasResp - CZ_1_model
    
    return k0, k0_CI, E, E_CI, alphaA, alphaA_CI, alphaB, alphaB_CI\
        , r_squared, CZ_1_model, epsilon_expt

# master function
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_20_5_2/reb_20_5_2_data.csv')
        # columns: T, tau, yA_0, CZ_1

    # extract the data as arrays
    T = (df['T (°C)']).to_numpy() + 273.15
    tau = df['tau (s)'].to_numpy()
    yA_0 = df['yA_0'].to_numpy()
    global gMeasResp
    gMeasResp = (df['CZ_1 (mmol/L)'].to_numpy())/1000

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([T, tau, yA_0]))

    # calculate the quantities of interest
    k0, k0_CI, E, E_CI, alphaA, alphaA_CI, alphaB, alphaB_CI, r_squared\
        , CZ_1_model, epsilon_expt = quantities_of_interest(adjusted_inputs)
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(gMeasResp, CZ_1_model, color = 'k', marker='o', ls=''
             , label = 'Data')
    plt.plot([min(gMeasResp),max(gMeasResp)],[min(gMeasResp),max(gMeasResp)]
             , color = 'r', ls = '-', label = 'Parity Line')
    plt.xlabel("$C_{Z,1}$ (M)")
    plt.ylabel("$C_{Z,1, model}$ (M)")
    plt.legend
    plt.tight_layout()
    plt.savefig('reb_20_5_2/python/reb_20_5_2_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(T-273.15, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("T (°C)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_2/python/reb_20_5_2_T_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(tau, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\\tau$ (s)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_2/python/reb_20_5_2_tau_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(yA_0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{A,0}$")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_2/python/reb_20_5_2_yA0_residual.png')
    plt.show()

    # save the results to a .csv file
    data = [['k0', f'{k0:.3g}', 'mol L^-1^ atm^-2^ s^-1^'],
        ['k0_lower_limit', f'{k0_CI[0]:.3g}', 'mol L^-1^ atm^-2^ s^-1^'],
        ['k_upper_limit', f'{k0_CI[1]:.3g}', 'mol L^-1^ atm^-2^ s^-1^'],
        ['E', f'{E:.3g}', 'kJ mol^-1^'],
        ['E_lower_limit', f'{E_CI[0]:.3g}', 'kJ mol^-1^'],
        ['E_upper_limit', f'{E_CI[1]:.3g}', 'kJ mol^-1^'],
        ['alphaA', f'{alphaA:.3g}', ''],
        ['alphaA_lower_limit', f'{alphaA_CI[0]:.3g}', ''],
        ['alphaA_upper_limit', f'{alphaA_CI[1]:.3g}', ''],
        ['alphaB', f'{alphaB:.3g}', ''],
        ['alphaB_lower_limit', f'{alphaB_CI[0]:.3g}', ''],
        ['alphaB_upper_limit', f'{alphaB_CI[1]:.3g}', ''],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    result.to_csv("reb_20_5_2/python/reb_20_5_2_results.csv", index=False)

    print(' ')
    print(result)

if __name__=="__main__":
    perform_the_calculations()