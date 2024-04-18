"""Calculations for Reaction Engineering Basics Example 9.6.2"""

# import libraries
import math
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# given, known and calculated constants available to all functions
# given
V = 10.0E3 # cm^3
Vex = 1.4E3 # cm^3
Uex = 138. # cal /ft^2 /min /K
Aex = 1200./929. # ft^2
Tex_0 = 40. + 273.15 # K
Tex_in = 40. + 273.15 # K
mDot_ex = 100. # g /min
rho = 1.0 # g /cm^3
rho_ex = 1.0 # g /cm^3
Cp = 1.0 # cal /g /K
Cp_ex = 1.0 # cal /g /K
CA_0 = 5.0E-3 # mol /cm^3
CB_0 = 7.0E-3 # mol /cm^3
dH_1 = -16.7E3 # cal /mol
dH_2 = -14.3E3 # cal /mol
k0_1 = 9.74E12 # cm^3 /mol /min
E_1 = 20.1E3 # cal /mol
k0_2 = 2.38E13 # /min
E_2 = 25.3E3 # cal /mol
t_f = 30. # min
fA_f = 0.45
# known
Re = 1.987 # cal /mol /K
# calculated
nA_0 = CA_0*V
nB_0 = CB_0*V
nA_f = nA_0*(1 - fA_f)

# current value of T0 available to all functions
T_0 = float('nan')

# derivatives function
def derivatives(ind, dep):
	# extract necessary dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    T = dep[5]
    Tex = dep[6]

	# calculate the rate
    CA = nA/V
    CB = nB/V
    k_1 = k0_1*math.exp(-E_1/Re/T)
    r_1 = k_1*CA*CB
    k_2 = k0_2*math.exp(-E_2/Re/T)
    r_2 = k_2*CA

    # calculate the rate of heat exchange
    Qdot = Uex*Aex*(Tex - T)

	# evaluate the derivatives
    dnAdt = -(r_1 + r_2)*V
    dnBdt = -r_1*V
    dnXdt = r_1*V
    dnYdt = r_1*V
    dnZdt = r_2*V
    dTdt = (Qdot-(r_1*dH_1 + r_2*dH_2)*V)/rho/V/Cp
    dTexdt = (-Qdot - mDot_ex*Cp_ex*(Tex - Tex_in))/rho_ex/Vex/Cp_ex

	# return the derivatives
    return [dnAdt, dnBdt, dnXdt, dnYdt, dnZdt, dTdt, dTexdt]

# reactor model function
def profiles():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, nB_0, 0.0, 0.0, 0.0, T_0, Tex_0])

	# define the stopping criterion
    f_var = 0
    f_val = t_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nX = dep[2,:]
    nY = dep[3,:]
    nZ = dep[4,:]
    T = dep[5,:]
    Tex = dep[6,:]

    # return all profiles
    return t, nA, nB, nX, nY, nZ, T, Tex

# residual function
def residual(guess):
    # make the guess available to all functions
    global T_0
    T_0 = guess[0]

    # solve the reactor design equations
    _, nA, _, _, _, _, _, _ = profiles()

    # extract the calculated molar amount
    nA_f_calc = nA[-1]

    # evaluate the residual
    resid = nA_f_calc - nA_f

    # return the residual
    return resid

# function that performs the analysis
def perform_the_analysis():
    # allow this function to set T0
    global T_0

    # set initial guess for T0
    initial_guess = Tex_in + 20

    # solve the implict equation for T0
    soln = sp.optimize.root(residual,initial_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The initial temperature was NOT found: {soln.message}")

    # extract the result
    T_0 = soln.x[0]

    # solve the reactor design equations
    _, nA, _, nX, _, nZ, T, Tex = profiles()

    # calculate the other quantities of interest
    T_0 = T_0 - 273.15
    T_f = T[-1] - 273.15
    Tex_out = Tex[-1] - 273.15
    sel_X_Z = nX[-1]/nZ[-1]

    # tabulate the results
    data = [['T0', f'{T_0}', '°C'],
    ['Tf', f'{T_f}', '°C'],
    ['Te_out',f'{Tex_out}', '°C'],
    ['sel_X_Z',f'{sel_X_Z}','mol X per mol Z']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(f'Initial Temperature: {T_0:.3g} °C')   
    print(f'Final Temperature: {T_f:.3g} °C')
    print(f'Outlet Coolant Temperature: {Tex_out:.3g} °C')
    print(f'Selectivity: {sel_X_Z:.3g} mol X per mol Z')

    # save the results
    results_df.to_csv('reb_9_6_2/python/results.csv'
                      , index=False)

    return

if __name__=="__main__":
    perform_the_analysis()
