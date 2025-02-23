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
R = 82.06 # cm^3 atm/mol/K
Ren = 1.987E-3 # kcal/mol

# globally available variables
global k0_current, E_current, T_current
k0_current = float('nan')
E_current = float('nan')
T_current = float('nan')

# derivatives function
def derivatives(t,dep):
    # extract the dependent variables
    nA = dep[0]
    nB = dep[1]

    # calculate the rate
    k = k0_current*np.exp(-E_current/Ren/T_current)
    PA = nA*R*T_current/V
    PB = nB*R*T_current/V
    if PA <= 0:
        r = 0.0
    elif PB <=0:
        r = 0.0
    else:
        r = k*PA*np.sqrt(PB)

    # evaluate the derivatives
    dnAdt = -r*V
    dnBdt = -r*V
    dnZdt = r*V

    # collect and return the derivatives
    ddt = np.array([dnAdt, dnBdt, dnZdt])
    return ddt

# BSTR model function
def profiles(T, PA0, tf):
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
def predicted_responses(adj_inputs, k0, E):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # make the rate coefficient available to all functions
    global k0_current, E_current, T_current
    k0_current = k0
    E_current = E

    # loop through all of the experiments
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs and make T available to all functions
        T_current = input[0]
        PA0 = input[1]
        tf = input[2]

        # solve the BSTR design equations
        t, nA, nB, nZ = profiles(T_current, PA0, tf)

        # calculate the response
        nAout = nA[-1]
        nBout = nB[-1]
        nZout = nZ[-1]
        resp[i] = (nAout + nBout + nZout)*R*T_current/V
    
    # return the responses
    return resp

# quantities of interest function
def quantities_of_interest(adjusted_inputs, Pf):
    # make a guess for k0 and E
    par_guess = [1.0, 14.0]

    # estimate k0 and E
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs, Pf, 
            predicted_responses, False)

    # extract the results
    k0 = beta[0]
    k0_CI = beta_ci[0,:]
    E = beta[1]
    E_CI = beta_ci[1,:]

    # calculate the model-predicted responses and the residuals
    Pf_model = predicted_responses(adjusted_inputs, k0, E)
    epsilon_expt = Pf - Pf_model

    return k0, k0_CI, E, E_CI, r_squared, Pf_model, epsilon_expt

# master function
def perform_the_calculations():
    # Read the experimental data into a dataframe
    df = pd.read_csv('reb_19_5_3/reb_19_5_3_data.csv')
            # columns: PA0, tf, P

    # extract the data as arrays
    T = df['T'].to_numpy() + 273.15
    PA0 = df['PA0'].to_numpy()
    tf = df['tf'].to_numpy()
    Pf = df['Pf'].to_numpy()

    # combine the adjusted inputs in a matrix
    adjusted_inputs = np.transpose(np.array([T, PA0, tf]))

    # calculate the quantities of interest
    k0, k0_CI, E, E_CI, r_squared, Pf_model, epsilon_expt \
        = quantities_of_interest(adjusted_inputs,Pf)

    # Display the results
    print(' ')
    print(f'k0: {k0:.3g} mol/cc/min/atm^1.5, 95% CI [{k0_CI[0]:.3g}, {k0_CI[1]:.3g}]')
    print(f'E: {E:.3g} kcal/mol, 95% CI [{E_CI[0]:.3g}, {E_CI[1]:.3g}]')
    print(f'R-squared: {r_squared:.3g}')
    print(' ')

    # Save the results
    data = [['k0', f'{k0:.3g}', 'mol cm^-3^ min^-1^ atm^-1.5^'],
        ['k0_lower_limit', f'{k0_CI[0]:.3g}', 'mol cm^-3^ min^-1^ atm^-1.5^'],
        ['k0_upper_limit', f'{k0_CI[1]:.3g}', 'mol cm^-3^ min^-1^ atm^-1.5^'],
        ['E', f'{E:.3g}', 'kcal mol^-1^'], 
        ['E_lower_limit', f'{E_CI[0]:.3g}', 'kcal mol^-1^'], 
        ['E_upper_limit', f'{E_CI[1]:.3g}', 'kcal mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    result.to_csv('reb_19_5_3/python/reb_19_5_3_results.csv', index=False)

    # create, display and save the parity plot
    plt.figure(1) 
    plt.plot(Pf, Pf_model, color = 'k', marker='o', ls='', label="Data")
    plt.plot([np.min(Pf), np.max(Pf)], [np.min(Pf), np.max(Pf)], color = 'r'
             , label = "Parity Line")
    plt.xlabel("Measured Pressure (atm)")
    plt.ylabel("Predicted Pressure (atm)")
    plt.legend
    plt.savefig('reb_19_5_3/python/reb_19_5_3_parity.png')
    plt.show()

    # create, display and save the residuals plots
    # residuals vs. T
    plt.figure(2) 
    plt.plot(T-273.15, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Residual (atm)")
    plt.savefig('reb_19_5_3/python/reb_19_5_3_T_residuals.png')
    plt.show()

    # residuals vs. PA0
    plt.figure(3) 
    plt.plot(PA0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("Initial Partial Pressure of A (atm)")
    plt.ylabel("Residual (atm)")
    plt.savefig('reb_19_5_3/python/reb_19_5_3_PA0_residuals.png')
    plt.show()

    # residuals vs. tRxn
    plt.figure(4) 
    plt.plot(tf, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("Reaction time (min)")
    plt.ylabel("Residual (atm)")
    plt.savefig('reb_19_5_3/python/reb_19_5_3_tf_residuals.png')
    plt.show()

if __name__=="__main__":
    perform_the_calculations()
