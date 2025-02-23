"""Calculations for Reaction Engineering Basics Example 20.5.3"""

#import libraries
import pandas as pd
import numpy as np
import scipy as sp
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
V = 3.0 # gal
Ren = 1.987 # BTU/lbmol/°R

# globally available variables
# full set of measured responses
gMeasResp = float('nan')
# adjusted inputs for the current experiment
gT = float('nan')
gVdot = float('nan')
gCA_0 = float('nan')
gCY_0 = float('nan')
gCZ_0 = float('nan')
# current rate expression parameters
gk0f = float('nan')
gEf = float('nan')
gk0r = float('nan')
gEr = float('nan')

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nA_1 = guess[0]
    nY_1 = guess[1]
    nZ_1 = guess[2]

    # calculate the other unknown quantities
    nA_0 = gCA_0*gVdot
    nY_0 = gCY_0*gVdot
    nZ_0 = gCZ_0*gVdot
    CA = nA_1/gVdot
    CY = nY_1/gVdot
    CZ = nZ_1/gVdot
    kf = gk0f*np.exp(-gEf/Ren/gT)
    kr = gk0r*np.exp(-gEr/Ren/gT)
    r = kf*CA**2 - kr*CY*CZ

    # evaluate the residuals
    epsilon_1 = nA_0 - nA_1 - V*r
    epsilon_2 = nY_0 - nY_1 + V*r
    epsilon_3 = nZ_0 - nZ_1 + V*r

    # return the residuals
    return np.array([epsilon_1, epsilon_2, epsilon_3])

# CSTR model function
def unknowns(T, Vdot, CA_0, CY_0, CZ_0, CA_1, k0f, Ef, k0r, Er):
    # make the adjusted inputs and rate expression parameters globally available
    global gT, gVdot, gCA_0, gCY_0, gCZ_0, gk0f, gEf, gk0r, gEr
    gT = T
    gVdot = Vdot
    gCA_0 = CA_0
    gCY_0 = CY_0
    gCZ_0 = CZ_0
    gk0f = k0f
    gEf = Ef
    gk0r = k0r
    gEr = Er

    # guess the solution
    nA_guess = CA_1*Vdot
    xi = CA_0*Vdot - nA_guess
    nY_guess = CY_0*Vdot + xi
    nZ_guess = CZ_0*Vdot + xi
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
def predicted_responses(adj_inputs, betaf, Ef, betar, Er):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    CA_1_model = np.zeros(nExpts)

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # solve the reactor design equations
        nA_1, nY_1, nZ_1 = unknowns(input[0], input[1], input[2], input[3]
                                    , input[4], gMeasResp[i], 10**betaf
                                    , Ef, 10**betar, Er)
        
        # calculate the model-predicted response
        CA_1_model[i] = nA_1/input[1]

    # return the responses
    return CA_1_model

# quantities of interest function
def quantities_of_interest(adj_inputs):
    # make a guess for log_10 of kf and kr
    #par_guess = [6.0, 15000.0, 6.0, 15000.0]
    par_guess = [np.log10(5.34E5), 11800.0, np.log10(6.44E7), 18500.0]

    # estimate log_10 of kf and kr
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adj_inputs 
        , gMeasResp, predicted_responses, False)
    
    # extract the results
    k0f = 10.**beta[0]
    k0f_CI = 10.**beta_ci[0,:]
    Ef = beta[1]
    Ef_CI = beta_ci[1,:]
    k0r = 10.**beta[2]
    k0r_CI = 10.**beta_ci[2,:]
    Er = beta[3]
    Er_CI = beta_ci[3,:]
        
    # calculate the model-predicted response and the residuals
    CA_1_model = predicted_responses(adj_inputs, beta[0], beta[1], beta[2]
                                     , beta[3])
    epsilon_expt = gMeasResp - CA_1_model

    return k0f, k0f_CI, Ef, Ef_CI, k0r, k0r_CI, Er, Er_CI, r_squared\
        , CA_1_model, epsilon_expt

# master function
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_20_5_3/reb_20_5_3_data.csv')
        # columns: T, Vdot, CA_0, CY_0, CZ_0, CA_1

    # extract the data as arrays
    T = df['T (°F)'].to_numpy() + 459.7
    Vdot = df['Vdot (gal/min)'].to_numpy()
    CA_0 = df['CA_0 (lbmol/gal)'].to_numpy()
    CY_0 = df['CY_0 (lbmol/gal)'].to_numpy()
    CZ_0 = df['CZ_0 (lbmol/gal)'].to_numpy()
    global gMeasResp
    gMeasResp = df['CA_1 (lbmol/gal)'].to_numpy()

    # combine the adjusted inputs as a matrix
    adj_inputs = np.transpose(np.array([T, Vdot, CA_0, CY_0, CZ_0]))

    # calculate the quantities of interest
    k0f, k0f_CI, Ef, Ef_CI, k0r, k0r_CI, Er, Er_CI, r_squared, CA_1_model\
        , epsilon_expt = quantities_of_interest(adj_inputs)
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(gMeasResp, CA_1_model, color = 'k', marker='o', ls=''
             , label = 'Data')
    plt.plot([min(gMeasResp),max(gMeasResp)],[min(gMeasResp),max(gMeasResp)]
             , color = 'r', ls = '-', label = 'Parity Line')
    plt.xlabel("$C_{A,1}$ (gal lbmol$^{-1}$ min$^{-1}$)")
    plt.ylabel("$C_{A,1,model}$ (gal lbmol$^{-1}$ min$^{-1}$)")
    plt.legend
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(T - 459.7, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("T (°F)")
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_T_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(Vdot, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\dot{V}$ (gal min$^{-1}$)")
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_Vdot_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(CA_0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{A,0}$ (lbmol gal$^{-1}$)")
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_CA0_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(CY_0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Y,0}$ (lbmol gal$^{-1}$)")
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_CY0_residual.png')
    plt.show()

    plt.figure() 
    plt.plot(CZ_0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{Z,0}$ (lbmol gal$^{-1}$)")
    plt.ylabel("Residual (lbmol gal$^{-1}$)")
    plt.tight_layout()
    plt.savefig('reb_20_5_3/python/reb_20_5_3_CZ0_residual.png')
    plt.show()

    # save the results to a .csv file
    data = [['k0f', f'{k0f:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0f_lower_limit', f'{k0f_CI[0]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0f_upper_limit', f'{k0f_CI[1]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['Ef', f'{Ef:.3g}', 'BTU/lbmol'],
        ['Ef_lower_limit', f'{Ef_CI[0]:.3g}', 'BTU/lbmol'],
        ['Ef_upper_limit', f'{Ef_CI[1]:.3g}', 'BTU/lbmol'],
        ['k0r', f'{k0r:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0r_lower_limit', f'{k0r_CI[0]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['k0r_upper_limit', f'{k0r_CI[1]:.3g}', 'gal lbmol^-1^ min^-1^'],
        ['Er', f'{Er:.3g}', 'BTU/lbmol'],
        ['Er_lower_limit', f'{Er_CI[0]:.3g}', 'BTU/lbmol'],
        ['Er_upper_limit', f'{Er_CI[1]:.3g}', 'BTU/lbmol'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    print(" ")
    print(result)
    result.to_csv("reb_20_5_3/python/reb_20_5_3_results.csv", index=False)

if __name__=="__main__":
    perform_the_calculations()