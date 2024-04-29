""" Calculations for Reaction Engineering Basics Example I.6"""

# import libraries
import scipy as sp
import numpy as np
import pandas as pd

# given and known constants
n_A_in = 500. # mol/h
V_fluid = 500. # L
n_Z_in = 0. # mol/h
Cp_vol = 1170. # cal/L/K
V_flow = 250. # L/h
T_in_K = 423. # K
dH_rxn = 18200. # cal/mol
k_0 = 1.14E9 # L/mol/h
E = 16200. # cal/mol
Re = 1.987 # cal/mol/K

# residuals function
def residuals(guess):
    # Extract the guess values
    n_A_guess = guess[0]
    n_Z_guess = guess[1]
    temp_guess_K = guess[2]

    # Calculate the concentration of A
    C_A = n_A_guess/V_flow

    # Calculate the rate
    r = k_0*np.exp(-E/Re/temp_guess_K)*C_A**2
    
    # Evaluate the residuals
    residual_1 = n_A_in - n_A_guess - V_fluid*r
    residual_2 = n_Z_in - n_Z_guess + V_fluid*r
    residual_3 = Cp_vol*V_flow*(temp_guess_K - T_in_K) + V_fluid*r*dH_rxn

    # Return the residuals
    return [residual_1, residual_2, residual_3]

# reactor model
def unknowns():
    # Initial guess
    init_guess = [n_A_in/2.0, n_A_in/2.0, T_in_K - 10.0]

    # Solve the ATEs
    soln = sp.optimize.root(residuals,init_guess)

    # Check that the solution is converged
    if not(soln.success):
        print(f"A solution was NOT obtained: {soln.message}")

    # return the solution
    return soln


# perform the analysis
def perform_the_analysis():
    # solve the reactor design equations
    soln = unknowns()
    
    # Extract the solution
    n_A_out = soln.x[0]
    n_Z_out = soln.x[1]
    temp_out_K = soln.x[2]

    # calculate the other quantities of interest
    temp_out_C = temp_out_K - 273.15

    # tabulate the results
    data = [['Flow Rate of A', f'{n_A_out}', ' mol/h'],
    ['Flow Rate of Z', f'{n_Z_out}', ' mol/h'],
    ['Temperature',f'{temp_out_K - 273.15}', ' °C']]
    result = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(result)

    # save the results
    result.to_csv('reb_I_6/python/results.csv', index=False)
    return

if __name__=="__main__":
    perform_the_analysis()






