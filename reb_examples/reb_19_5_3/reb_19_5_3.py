"""Calculations for Reaction Engineering Basics Example 19.5.3"""

# import libraries
import pandas as pd
import numpy as np
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
V = 100.0 # cm^3
P0 = 6.0 # atm
T = 275 + 273.15 # K
R = 82.06 # cm^3 atm/mol/K

# make current value of k available to all functions
global k_current
k_current = float('nan')

# derivatives function
def derivatives(t,dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]

    # calculate the rate
    PA = nA*R*T/V
    PB = nB*R*T/V
    if PA <= 0:
        r = 0.0
    elif PB <=0:
        r = 0.0
    else:
        r = k_current*PA*np.sqrt(PB)

    # evaluate the derivatives
    dnAdt = -r*V
    dnBdt = -r*V
    dnZdt = r*V

    # collect and return the derivatives
    ddt = np.array([dnAdt, dnBdt, dnZdt])
    return ddt

# BSTR model function
def profiles(PA0, tf):
    # initial values
    ind0 = 0
    nA0 = PA0*V/R/T
    PB0 = P0 - PA0
    nB0 = PB0*V/R/T
    dep0 = [nA0, nB0, 0.0]
    stop_var = 0

    # solve the BSTR reactor design equations
    t, dep, success, message = solve_ivodes(ind0, dep0, stop_var, tf, derivatives
                                            ,True)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")
    
    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nZ = dep[2,:]

    # return the profiles
    return t, nA, nB, nZ

# predicted responses function
def predicted_responses(adj_inputs, k_log_10):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # make the rate coefficient available to all functions
    global k_current
    k_current = 10.**k_log_10

    # loop through all of the experiments
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs and make T available to all functions
        PA0 = input[0]
        tf = input[1]

        # solve the BSTR design equations
        t, nA, nB, nZ = profiles(PA0,tf)

        # calculate the response
        nAout = nA[-1]
        nBout = nB[-1]
        nZout = nZ[-1]
        resp[i] = (nAout + nBout + nZout)*R*T/V
    
    # return the responses
    return resp

# function that performs the calculations
def perform_the_calculations():
    # Read the experimental data into a dataframe
    df = pd.read_csv('reb_19_5_3/reb_19_5_3_data.csv')
            # columns: PA0, tf, P

    # extract the data as arrays
    PA0 = df['PA0'].to_numpy()
    tf = df['tf'].to_numpy()
    Pf = df['Pf'].to_numpy()

    # combine the adjusted inputs in a matrix
    adjusted_inputs = np.transpose(np.array([PA0, tf]))

    # make a guess for log_10 of k
    par_guess = [-6.0]

    # estimate log_10 of k
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs, Pf, 
            predicted_responses, False)

    # extract the results
    k = 10.**beta[0]
    k_ll = 10.**beta_ci[0,0]
    k_ul = 10.**beta_ci[0,1]

    # Display the results
    print(' ')
    print(f'k: {k:.3g} mol/cc/min/atm^1.5, 95% CI [{k_ll:.3g}, {k_ul:.3g}]')
    print(f'R-squared: {r_squared:.3g}')
    print(' ')

    # Save the results
    data = [['k', f'{k:.3g}', 'mol cm^-3^ min^-1^ atm^-1.5^'],
        ['k_lower_limit', f'{k_ll:.3g}', 'mol cm^-3^ min^-1^ atm^-1.5^'],
        ['k_upper_limit', f'{k_ul:.3g}', 'mol cm^-3^ min^-1^ atm^-1.5^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    result.to_csv('reb_19_5_3/python/reb_19_5_3_results.csv', index=False)

    # calculate the model-predicted responses and the residuals
    P_model = predicted_responses(adjusted_inputs, beta[0])
    residuals = Pf - P_model

    # create, display and save the parity plot
    plt.figure(1) 
    plt.plot(Pf, P_model, color = 'k', marker='o', ls='')
    plt.plot([np.min(Pf), np.max(Pf)], [np.min(Pf), np.max(Pf)], color = 'r')
    plt.xlabel("Measured Pressure (atm)")
    plt.ylabel("Predicted Pressure (atm)")
    plt.savefig('reb_19_5_3/python/reb_19_5_3_parity.png')
    plt.show()

    # create, display and save the residuals plots
    # residuals vs. PA0
    plt.figure(2) 
    plt.plot(PA0, residuals, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("Initial Partial Pressure of A (atm)")
    plt.ylabel("Residual (atm)")
    plt.savefig('reb_19_5_3/python/reb_19_5_3_PA0_residuals.png')
    plt.show()

    # residuals vs. tRxn
    plt.figure(2) 
    plt.plot(tf, residuals, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("Reaction time (min)")
    plt.ylabel("Residual (atm)")
    plt.savefig('reb_19_5_3/python/reb_19_5_3_tf_residuals.png')
    plt.show()

if __name__=="__main__":
    perform_the_calculations()
