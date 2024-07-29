"""Calculations for Reaction Engineering Basics Example 19.5.2"""

#import libraries
import pandas as pd
import numpy as np
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# given and known constants
V = 500.0 # cc
Re = 1.987E-3 # kcal/mol/K
Rpv = 82.06 # cc atm/mol/K

# make current values of k0, E and T available to all functions
k0_current = float('nan')
E_current = float('nan')
T_current = float('nan')

# derivatives function
def derivatives(t,dep):
    # get the dependent variables that are needed
    nA = dep[0]
    nB = dep[1]

    # calculate the partial pressures
    PA = nA*Rpv*T_current/V
    PB = nB*Rpv*T_current/V

    # calculate the rate
    k = k0_current*np.exp(-E_current/Re/T_current)
    r = k*PA*PB

    # calculate the time derivatives of the dependent variables
    dnAdt = -r*V
    dnBdt = -r*V
    dnYdt = r*V
    dnZdt = r*V

    # return an array containing the derivatives
    ddt = np.array([dnAdt, dnBdt, dnYdt, dnZdt])
    return ddt

# BSTR model function
def profiles(PA0, PB0, tf):
    # set initial values and stopping criterion
    t0 = 0
    nA0 = PA0*V/Rpv/T_current
    nB0 = PB0*V/Rpv/T_current
    dep0 = np.array([nA0, nB0, 0.0, 0.0])
    stop_var = 0

    # solve the BSTR reactor design equations
    t, dep, success, message = solve_ivodes(t0, dep0, stop_var, tf, derivatives
                                            ,True)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")
    
    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]

    # return the profiles
    return t, nA, nB, nY, nZ

# predicted responses function
def predicted_responses(adj_inputs, k0, E):
    # allocate storage for the responses
    nExpts = len(adj_inputs)
    resp = np.zeros(nExpts)

    # make the rate coefficient available to other functions
    global k0_current, E_current
    k0_current = k0
    E_current = E

    # loop through the experiments in the data set
    for i, input in enumerate(adj_inputs):
        # get the adjusted inputs
        global T_current
        T_current = input[0]
        PA0 = input[1]
        PB0 = input[2]
        tf = input[3]

        # solve the reactor design equations
        t, nA, nB, nY, nZ = profiles(PA0, PB0, tf)
        
        # calculate the model-predicted response
        nAf = nA[-1]
        nA0 = PA0*V/Rpv/T_current
        fA = (nA0 - nAf)/nA0
        resp[i] = fA

    # return the responses
    return resp

# function that performs the calculations
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_19_5_2/reb_19_5_2_data.csv')

    # extract the data as arrays
    T = df['T'].to_numpy()
    PA0 = df['PA0'].to_numpy()
    PB0 = df['PB0'].to_numpy()
    tf = df['tf'].to_numpy()
    fA = df['fA'].to_numpy()
    T_K = T + 273.15

    # combine the adjusted inputs as a matrix
    adjusted_inputs = np.transpose(np.array([T_K, PA0, PB0, tf]))

    # guess the parameters
    par_guess = [1.0, 20.0]

    # estimate the parameters
    beta, beta_ci, r_squared = fit_to_SR_data(par_guess, adjusted_inputs 
            , fA,  predicted_responses, False)
    
    # extract the results
    k0 = beta[0]
    k0_ll = beta_ci[0,0]
    k0_ul = beta_ci[0,1]
    E = beta[1]
    E_ll = beta_ci[1,0]
    E_ul = beta_ci[1,1]
        
    # calculate the model-predicted y and the residuals
    fA_model = predicted_responses(adjusted_inputs, k0, E)
    residual = fA - fA_model
        
    # create, show, and save a parity plot
    plt.figure() 
    plt.plot(fA, fA_model, color = 'k', marker='o', ls='')
    plt.plot([min(fA),max(fA)],[min(fA),max(fA)], \
        color = 'r', ls = '-')
    plt.xlabel("$f_{A, expt}$")
    plt.ylabel("$f_{A, model}$")
    plt.tight_layout()
    plt.savefig('reb_19_5_2/python/reb_19_5_2_parity.png')
    plt.show()

    # create, show and save residuals plots
    plt.figure() 
    plt.plot(PA0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$P_{A,0}$ (atm)")
    plt.ylabel("Residual")
    plt.tight_layout()
    plt.savefig('reb_19_5_2/python/reb_19_5_2_residuals_PA0.png')
    plt.show()

    plt.figure() 
    plt.plot(PB0, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$P_{B,0}$ (atm)")
    plt.ylabel("Residual")
    plt.tight_layout()
    plt.savefig('reb_19_5_2/python/reb_19_5_2_residuals_PB0.png')
    plt.show()

    plt.figure() 
    plt.plot(tf, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$t_f$ (min)")
    plt.ylabel("Residual")
    plt.tight_layout()
    plt.savefig('reb_19_5_2/python/reb_19_5_2_residuals_tf.png')
    plt.show()

    plt.figure() 
    plt.plot(T, residual, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("T (°C)")
    plt.ylabel("Residual")
    plt.tight_layout()
    plt.savefig('reb_19_5_2/python/reb_19_5_2_residuals_T.png')
    plt.show()

    # report the results
    print(' ')
    print(f'k0: {k0:.3g} mol/cc/min/atm^2, 95% CI [{k0_ll:.3g}, {k0_ul:.3g}]')
    print(f'E: {E:.3g} kcal/mol, 95% CI [{E_ll:.3g}, {E_ul:.3g}]')
    print(f'R-squared: {r_squared:.3g}')
    print(' ')

    # save the results to a .csv file
    data = [['k0', f'{k0:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
        ['k0_lower_limit', f'{k0_ll:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
        ['k0_upper_limit', f'{k0_ul:.3g}', 'mol cm^-3^ min^-1^ atm^-2^'],
        ['E', f'{E:.3g}', 'kcal mol^-1^'],
        ['E_lower_limit', f'{E_ll:.3g}', 'kcal mol^-1^'],
        ['E_upper_limit', f'{E_ul:.3g}', 'kcal mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    result.to_csv("reb_19_5_2/python/reb_19_5_2_results.csv", index=False)

if __name__=="__main__":
    perform_the_calculations()