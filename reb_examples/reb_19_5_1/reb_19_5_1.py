"""Calculations for Reaction Engineering Basics Example 19.5.1"""

#import libraries
import pandas as pd
import numpy as np
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
from score_utils import Arrhenius_parameters
import matplotlib.pyplot as plt

# given and known constants
V = 1.0 # L
R = 8.314E-3 # kJ/mol/K

# globally available variables
k0_current = float('nan')
E_current = float('nan')
T_current = float('nan')

# derivatives function
def derivatives(t,dep):
    # get the dependent variables
    nA = dep[0]

    # calculate the rate coefficient
    k = k0_current*np.exp(-E_current/R/T_current)

    # calculate the concentration
    CA = nA/V

    # calculate the rate
    r = k*CA

    # calculate the time derivatives of the dependent variables
    dnAdt = -r*V
    dnZdt = r*V

    # return an array containing the derivatives
    ddt = np.array([dnAdt, dnZdt])
    return ddt

# BSTR model function
def profiles(k0,E,T,CA0,tf):
    # make the rate expression parameters and adjusted inputs available
    global k0_current, E_current, T_current
    k0_current = k0
    E_current = E
    T_current = T

    # set initial values and stopping criterion
    t0 = 0
    nA0 = CA0*V
    n0 = np.array([nA0, 0.0])
    stop_var = 0

    # solve the BSTR reactor design equations
    t, dep, success, message = solve_ivodes(t0, n0, stop_var, tf, derivatives
                                            ,True)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")
    
    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]

    # return the profiles
    return t, nA, nZ

# predicted responses function
def predicted_responses(adj_inputs, k0_log_10, E):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # make the rate coefficient available to other functions
    k0 = 10.**k0_log_10

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs
        T = input[0]
        CA0 = input[1]
        tf = input[2]        

        # solve the reactor design equations
        t, nA, nZ = profiles(k0, E, T, CA0, tf)
        
        # calculate the model-predicted response
        nAf = nA[-1]
        CAf = nAf/V
        resp[i] = CAf

    # return the responses
    return resp

# quantities of interest function
def quantities_of_interest(adjusted_inputs, CAf):
    # make a guess for log_10 of k0 and E
    par_guess = [2.0, 20.0]

    # estimate the parameters
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs 
            , CAf,  predicted_responses, False)
    
    # extract the results
    k0 = 10.**beta[0]
    k0_CI = 10.**beta_ci[0,:]
    E = beta[1]
    E_CI = beta_ci[1,:]

    # calculate the model-predicted y and the residuals
    CAf_model = predicted_responses(adjusted_inputs, beta[0], beta[1])
    epsilon_expt = CAf - CAf_model

    return k0, k0_CI, E, E_CI, CAf_model, r_squared, epsilon_expt


# master function
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_19_5_1/reb_19_5_1_data.csv')
        # columns: Experiment, T, CA0, tf, CAf
    
    # extract the data as arrays
    T = df['T'].to_numpy()
    CA0 = df['CA0'].to_numpy()
    tf = df['tf'].to_numpy()
    CAf = df['CAf'].to_numpy()
    T_K = T + 273.15

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([T_K, CA0, tf]))

    # calculate the quantities of interest
    k0, k0_CI, E, E_CI, CAf_model, r_squared, epsilon_expt = \
        quantities_of_interest(adjusted_inputs, CAf)
    
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(CAf, CAf_model, color = 'k', marker='o', ls='', label = 'data')
    plt.plot([min(CAf),max(CAf)],[min(CAf),max(CAf)], color = 'r', ls = '-'
             , label = 'parity line')
    plt.xlabel("$C_{A, expt}$ (M)")
    plt.ylabel("$C_{A, model}$ (M)")
    plt.legend()
    plt.tight_layout()
    plt.savefig('reb_19_5_1/python/reb_19_5_1_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(CA0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{A,0}$ (M)")
    plt.ylabel("Experiment Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_19_5_1/python/reb_19_5_1_residuals_CA0.png')
    plt.show()

    plt.figure() 
    plt.plot(tf, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$t_f$ (min)")
    plt.ylabel("Experiment Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_19_5_1/python/reb_19_5_1_residuals_tf.png')
    plt.show()

    plt.figure() 
    plt.plot(T, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$T$ (°C)")
    plt.ylabel("Experiment Residual (M)")
    plt.tight_layout()
    plt.savefig('reb_19_5_1/python/reb_19_5_1_residuals_T.png')
    plt.show()
    
    # show and save the fitting results
    data = [['k0', f'{k0:.3g}', 'min^-1^'],
        ['k0_lower_limit', f'{k0_CI[0]:.3g}', 'min^-1^'],
        ['k0_upper_limit', f'{k0_CI[1]:.3g}', 'min^-1^'],
        ['E', f'{E:.3g}', 'kJ mol^-1^'],
        ['E_lower_limit', f'{E_CI[0]:.3g}', 'kJ mol^-1^'],
        ['E_upper_limit', f'{E_CI[1]:.3g}', 'kJ mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    print(result)
    result.to_csv("reb_19_5_1/python/reb_19_5_1_results.csv", 
                index=False)

if __name__=="__main__":
    perform_the_calculations()