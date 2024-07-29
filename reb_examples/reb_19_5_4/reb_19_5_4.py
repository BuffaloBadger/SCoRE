"""Calculations for Reaction Engineering Basics Example 19.5.4"""

# import libraries
import pandas as pd
import numpy as np
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
V = 50.0E-3 # L

# make current values of Vmax and Km available to all functions
global Vmax_current, Km_current
Vmax_current = float('nan')
Km_current = float('nan')

# derivatives function
def derivatives(t,dep):
    # extract the dependent variables
    nS = dep[0]

    # calculate the rate
    CS = nS/V
    r = Vmax_current*CS/(Km_current + CS)

    # evaluate the derivatives
    dnSdt = -r*V
    dnPdt = r*V

    # collect and return the derivatives
    ddt = np.array([dnSdt, dnPdt])
    return ddt

# BSTR model function
def profiles(CS0, tf):
    # initial values
    ind0 = 0
    nS0 = CS0*V
    dep0 = [nS0, 0.0]
    stop_var = 0

    # solve the BSTR reactor design equations
    t, dep, success, message = solve_ivodes(ind0, dep0, stop_var, tf, derivatives
                                            ,True)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")
    
    # extract the dependent variable profiles
    nS = dep[0,:]
    nP = dep[1,:]

    # return the profiles
    return t, nS, nP

# predicted responses function
def predicted_responses(adj_inputs, log_Vmax_guess, log_Km_guess):
    # allocate storage for the responses
    shape = np.shape(adj_inputs)
    n_data = shape[0]
    resp = np.ones(n_data)*float('nan')

    # make current values of Vmax and Km available to all functions
    global Vmax_current, Km_current
    Vmax_current = 10**log_Vmax_guess
    Km_current = 10**log_Km_guess

    # loop through all of the experiments
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs and make T available to all functions
        CS0 = input[0]
        tf = input[1]

        # solve the BSTR design equations
        t, nS, nP = profiles(CS0,tf)

        # calculate the response
        nPf = nP[-1]
        resp[i] = nPf/V
    
    # return the responses
    return resp

# function that performs the calculations
def perform_the_calculations():
    # Read the experimental data into a dataframe
    df = pd.read_csv('reb_19_5_4/reb_19_5_4_data.csv')
            # columns: CS0, tf, CPf

    # extract the data as arrays
    CS0 = df['CS0'].to_numpy()
    tf = df['tf'].to_numpy()
    CPf = df['CPf'].to_numpy()

    # combine the adjusted inputs in a matrix
    adjusted_inputs = np.transpose(np.array([CS0, tf]))

    # make a guess for log_10 of k
    par_guess = [0.0, 0.0]

    # estimate the kinetics parameters
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs, CPf, 
            predicted_responses, False)

    # extract the results
    Vmax = 10**beta[0]
    Vmax_lower = 10**beta_ci[0,0]
    Vmax_upper = 10**beta_ci[0,1]
    Km = 10**beta[1]
    Km_lower = 10**beta_ci[1,0]
    Km_upper = 10**beta_ci[1,1]

    # report the results
    print(' ')
    print(f'Vmax: {Vmax:.3g} mmol/L/min, 95% CI ['\
        + f'{Vmax_lower:.3g}, {Vmax_upper:.3g}]')
    print(f'Km: {Km:.3g} mmol/L, 95% CI [{Km_lower:.3g}, {Km_upper:.3g}]')
    print(f'R-squared: {r_squared:.3g}')
    print(' ')

    # save the results to a .csv file
    data = [['Vmax', f'{Vmax:.3g}', 'mmol L^-1^ min^-1^'],
        ['Vmax_lower_limit', f'{Vmax_lower:.3g}', 'mmol L^-1^ min^-1^'],
        ['Vmax_upper_limit', f'{Vmax_upper:.3g}', 'mmol L^-1^ min^-1^'],
        ['Km', f'{Km:.3g}', 'mmol L^-1^'],
        ['Km_lower_limit', f'{Km_lower:.3g}', 'mmol L^-1^'],
        ['Km_upper_limit', f'{Km_upper:.3g}', 'mmol L^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    result.to_csv("reb_19_5_4/python/reb_19_5_4_results.csv", index=False)

    # calculate the model-predicted responses
    CP_model = predicted_responses(adjusted_inputs,beta[0],beta[1])

    # calculate the residuals
    residuals = CPf - CP_model

    # create a parity plot
    plt.figure(1) 
    plt.plot(CPf, CP_model, color = 'k', marker='o', ls='')
    plt.plot([np.min(CPf), np.max(CPf)], [np.min(CPf), np.max(CPf)], color = 'r')
    plt.xlabel("experimental response (mmol L$^{-1}$)")
    plt.ylabel("model-predicted response (mmol L$^{-1}$)")

    # save and show the parity plot
    plt.savefig('reb_19_5_4/python/reb_19_5_4_parity.png')
    plt.show()

    # create a residuals plot for the reaction time
    plt.figure(2) 
    plt.plot(tf, residuals, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("Reaction time (min)")
    plt.ylabel("Residual (mmol L$^{-1}$)")

    # save and show the residuals plot for the reaction time
    plt.savefig('reb_19_5_4/python/reb_19_5_4_tf_residuals.png')
    plt.show()

    # create a residuals plot for the initial substrate concentration
    plt.figure(3) 
    plt.plot(CS0, residuals, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("Initial Substrate Concentration (mmol L$^{-1}$)")
    plt.ylabel("Residual (mmol L$^{-1}$)")

    # save and show the residuals plot for the initial substrate concentrationfilename = filepath + 'reb_19_4_residuals_vs_CS.png'
    plt.savefig('reb_19_5_4/python/reb_19_5_4_CS0_residuals.png')
    plt.show()

if __name__=="__main__":
    perform_the_calculations()
