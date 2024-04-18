"""Calculations for Example J.6.3 of Reaction Engineering Basics"""

# import libraries
import math
import numpy as np
import scipy as sp
import pandas as pd
from score_utils import solve_ivodes

# given and known constants
T_f = 325. # K
D = 5. # cm
dH = -14000. # cal /mol
Cp = 1.3 # cal /cm^3 /K
C_A_in = 0.0025 # mol /cm^3
C_Z_in = 0.0 # mol /min
T_in = 300. # K
L = 50.0 # cm
k_0 = 4.2E15 # cm^3 /mol /min
E = 18000. # cal /mol
Re = 1.987 # cal /mol /K

# make Vdot available to all functions
Vdot = float('nan')

# derivatives function
def derivatives(z,dep):
    # extract the dependent variables for this integration step
    nDot_A = dep[0]
    nDot_Z = dep[1]
    T = dep[2]

    # calculate the rate
    r = k_0*math.exp(-E/Re/T)*nDot_A**2/Vdot**2
    
    # evaluate the derivatives
    dnDotAdz = -math.pi*D**2/4*r
    dnDotZdz = math.pi*D**2/4*r
    dTdz = -math.pi*D**2/4*r*dH/Vdot/Cp

    # return the derivatives
    return [dnDotAdz, dnDotZdz, dTdz]

# residual function
def residual(guess):
    # allow this function to set Vdot
    global Vdot
    Vdot = guess[0]
  
    # solve the reactor design equations
    z, n_A, n_Z, T = profiles()

    # extract the calculated final T
    T_f_calc = T[-1]

    # evaluate the residual
    resid = T_f - T_f_calc

    # return the residual
    return resid

# reactor model function
def profiles():
    # set the initial values
    nDot_A_in = Vdot*C_A_in
    nDot_Z_in = Vdot*C_Z_in
    ind_0 = 0.0
    dep_0 = np.array([nDot_A_in, nDot_Z_in, T_in])

    # define the stopping criterion
    f_var = 0
    f_val = L
 
    # solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0
                            , f_var, f_val, derivatives)

    # check that a solution was found
    if not(success):
        print(f"The profiles were NOT found: {message}")

    # extract dependent variable profiles
    nDot_A = dep[0,:]
    nDot_Z = dep[1,:]
    T = dep[2,:]

    # return all profiles
    return z, nDot_A, nDot_Z, T

# function that performs the analysis
def perform_the_analysis():
    # calculate Vdot and make it available to all functions
    global Vdot

    # initial guess for Vdot
    initial_guess = 100.0
    
    # solve the implicit equation for Vdot
    soln = sp.optimize.root(residual,initial_guess)

    # check that the solution converged
    if not(soln.success):
        print(f"Vdot was NOT found: {soln.message}")

    # set Vdot
    Vdot = soln.x[0]

    # solve the reactor design equations
    z, nDot_A, nDot_Z, T = profiles()

    # tabulate the results
    Vdot_results_df = pd.DataFrame(columns=['item','value','units'])
    Vdot_results_df.loc[0] = ['Vdot' , Vdot, 'cm^3^ min^-1^']
    profile_results_df = pd.DataFrame({'z':z , 'nDot_A':nDot_A
                                       , 'nDot_Z':nDot_Z, 'T':T})
        
    # display the results
    print(' ')
    print(f'Inlet Volumetric Flow Rate: {T_in:.3g} cm^3/min')
    print(' ')
    print('Molar Flow and Temperature Profiles')
    print(profile_results_df.to_string(index=False))

    # save the results
    Vdot_results_file_spec = \
            './reb_J_6_3/python/Vdot_results.csv'
    profile_results_file_spec = \
            './reb_J_6_3/python/profile_results.csv'
    Vdot_results_df.to_csv(Vdot_results_file_spec, index=False)
    profile_results_df.to_csv(profile_results_file_spec, index=False)

if __name__=="__main__":
    perform_the_analysis()