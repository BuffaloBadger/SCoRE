"""Calculations for Reaction Engineering Basics Example 20.5.1"""

#import libraries
import pandas as pd
import numpy as np
import scipy as sp
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
P = 3.0 # atm
T = 450 + 273.15 # K
V = 1.0 # L (basis)
R = 0.08206 # L*atm/mol/K

# globally available variables
CZ_1 = float('nan')
tau = float('nan')
yA_0 = float('nan')
i_expt_current = 0
k_current = float('nan')
alphaA_current = float('nan')
alphaB_current = float('nan')

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nA_1 = guess[0]
    nB_1 = guess[1]
    nZ_1 = guess[2]

    # calculate the other unknown quantities
    Vdot_0 = V/tau[i_expt_current]
    nA_0 = yA_0[i_expt_current]*P*Vdot_0/R/T
    nB_0 = (1-yA_0[i_expt_current])*P*Vdot_0/R/T
    PA = nA_1*P/(nA_1+nB_1+nZ_1)
    PB = nB_1*P/(nA_1+nB_1+nZ_1)
    r = k_current*PA**alphaA_current*PB**alphaB_current

    # evaluate the residuals
    epsilon_1 = nA_0 - nA_1 - V*r
    epsilon_2 = nB_0 - nB_1 - V*r
    epsilon_3 = -nZ_1 + V*r

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3])

# CSTR model function
def unknowns():
    # guess the solution
    Vdot_0 = V/tau[i_expt_current]
    nA_0 = yA_0[i_expt_current]*P*Vdot_0/R/T
    nB_0 = (1-yA_0[i_expt_current])*P*Vdot_0/R/T
    nZ_guess = CZ_1[i_expt_current]*Vdot_0
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
def predicted_responses(adj_inputs, log_k_guess, alphaA_guess
                        , alphaB_guess):
    # make the parameters globally available
    global k_current, alphaA_current, alphaB_current
    k_current = 10.**log_k_guess
    alphaA_current = alphaA_guess
    alphaB_current = alphaB_guess

    # allocate storage for the responses
    nExpts = len(adj_inputs)
    CZ_1_model = np.zeros(nExpts)

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # make the experiment index globally available
        global i_expt_current
        i_expt_current = i

        # solve the reactor design equations
        nA_1, nB_1, nZ_1 = unknowns()
        
        # calculate the model-predicted response
        CZ_1_model[i] = nZ_1*P/(nA_1 + nB_1 + nZ_1)/R/T

    # return the responses
    return CZ_1_model

# function that performs the calculations
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_20_5_2/reb_20_5_2_data.csv')
        # columns: tau, yA_0, CZ_1

    # extract the data as arrays
    global CZ_1, tau, yA_0
    tau = df['tau'].to_numpy()
    yA_0 = df['yA_0'].to_numpy()
    CZ_1 = (df['CZ_1'].to_numpy())/1000

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([tau, yA_0]))

    # make a guess for log_10 of k
    par_guess = [-5.0, 1.3, 0.5]

    # estimate log_10 of k
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs 
        , CZ_1, predicted_responses, False)
    
    # extract the results
    k = 10.**beta[0]
    k_ll = 10.**beta_ci[0,0]
    k_ul = 10.**beta_ci[0,1]
    alphaA = beta[1]
    alphaA_ll = beta_ci[1,0]
    alphaA_ul = beta_ci[1,1]
    alphaB = beta[2]
    alphaB_ll = beta_ci[2,0]
    alphaB_ul = beta_ci[2,1]
        
    # calculate the model-predicted response and the residuals
    CZ_1_model = predicted_responses(adjusted_inputs, beta[0], alphaA, alphaB)
    residual = CZ_1 - CZ_1_model
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(CZ_1, CZ_1_model, color = 'k', marker='o', ls='')
    plt.plot([min(CZ_1),max(CZ_1)],[min(CZ_1),max(CZ_1)], color = 'r', ls = '-')
    plt.xlabel("$C_{Z,1}$ (M)")
    plt.ylabel("$C_{Z,1, model}$ (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_2/python/reb_20_5_2_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(tau, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\tau$ (s)")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_2/python/reb_20_5_2_tau_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(yA_0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$y_{A,0}$")
    plt.ylabel("Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_20_5_2/python/reb_20_5_2_yA0_residual.png')
    plt.show()

    # save the results to a .csv file
    data = [['k', f'{k:.3g}', 'L mol^-1^ s^-1^'],
        ['k_lower_limit', f'{k_ll:.3g}', 'L mol^-1^ s^-1^'],
        ['k_upper_limit', f'{k_ul:.3g}', 'L mol^-1^ s^-1^'],
        ['alphaA', f'{alphaA:.3g}', ''],
        ['alphaA_lower_limit', f'{alphaA_ll:.3g}', ''],
        ['alphaA_upper_limit', f'{alphaA_ul:.3g}', ''],
        ['alphaB', f'{alphaB:.3g}', ''],
        ['alphaB_lower_limit', f'{alphaB_ll:.3g}', ''],
        ['alphaB_upper_limit', f'{alphaB_ul:.3g}', ''],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    result.to_csv("reb_20_5_2/python/reb_20_5_2_results.csv", index=False)

if __name__=="__main__":
    perform_the_calculations()