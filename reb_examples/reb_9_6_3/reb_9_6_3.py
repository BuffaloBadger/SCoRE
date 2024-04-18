"""Calculations for Reaction Engineering Basics Example 9.6.3"""

# import libraries
import math
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# given, known, and calculated constants available to all functions
# given
PA_0 = 1. # atm
PB_0 = 2. # atm
T_0 = 25 + 273.14 # K
V = 2. # L
Aex = 600. # cm^2
Uex = 0.6 # cal /cm^2 /min /K
Tex = 30 + 273.15 # K
k_0_1 =  3.34E9 # mol /cm^3 /min /atm^2
k_0_2 = 4.99E9 # mol /cm^3 /min /atm^2
E_1 = 20.5E3 # cal /mol
E_2 = 21.8E3 # cal /mol
dH_1 = -6300. # cal /mol
dH_2 = -6900. # cal /mol
Cp_A = 7.4 # cal /mol /K
Cp_B = 8.6 # cal /mol /K
Cp_D = 10.7 # cal /mol /K
Cp_Z = 5.2 # cal /mol /K
Cp_U = 10.3 # cal /mol /K
# known
Re = 1.987 # cal /mol /K
Rw = 82.057 # cm^3 atm /mol /K
# calculated
nA_0 = PA_0*V/Rw/T_0
nB_0 = PB_0*V/Rw/T_0
P_0 = PA_0 + PB_0

# derivatives function
def derivatives(ind, dep):
	# extract necessary dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    nD = dep[2]
    nZ = dep[3]
    nU = dep[4]
    T = dep[5]

	# calculate the rate
    PA = nA*Rw*T/V
    PB = nB*Rw*T/V
    PD = nD*Rw*T/V
    k_1 = k_0_1*math.exp(-E_1/Re/T)
    k_2 = k_0_2*math.exp(-E_2/Re/T)
    r_1 = k_1*PA*PB
    r_2 = k_2*PD*PB

    # calculate the rate of heat exchange
    Qdot = Uex*Aex*(Tex - T)

	# Create mass matrix, setting all elements to zero
    mass_matrix = np.zeros((7,7))

    # Add 1 on the diagonal for the first 5 rows
    mass_matrix[0,0] = 1.0
    mass_matrix[1,1] = 1.0
    mass_matrix[2,2] = 1.0
    mass_matrix[3,3] = 1.0
    mass_matrix[4,4] = 1.0

    # Add the elements for the energy balance
    mass_matrix[5,5] = nA*Cp_A + nB*Cp_B + nD*Cp_D + nZ*Cp_Z + nU*Cp_U
    mass_matrix[5,6] = -V*Re/Rw

    # Add the elements for the ideal gas law equation
    mass_matrix[6,0] = Rw*T
    mass_matrix[6,1] = Rw*T
    mass_matrix[6,2] = Rw*T
    mass_matrix[6,3] = Rw*T
    mass_matrix[6,4] = Rw*T
    mass_matrix[6,5] = Rw*(nA + nB + nD + nZ +nU)
    mass_matrix[6,6] = -V

    # Create right side vector
    rhs1 = -r_1*V
    rhs2 = (-r_1 -r_2)*V
    rhs3 = (r_1-r_2)*V
    rhs4 = (r_1+r_1)*V
    rhs5 = r_2*V
    rhs6 = Qdot -(r_1*dH_1 + r_2*dH_2)*V
    rhs7 = 0.0
    rhs = np.array([rhs1, rhs2, rhs3, rhs4, rhs5, rhs6, rhs7])

    # Evaluate the derivatives
    derivs = sp.linalg.solve(mass_matrix, rhs)

    # Return the derivatives
    return derivs

# reactor model function
def profiles(t_f):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, nB_0, 0.0, 0.0, 0.0, T_0, P_0])

	# define the stopping criterion
    f_var = 0
    f_val = t_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: message")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nD = dep[2,:]
    nZ = dep[3,:]
    nU = dep[4,:]
    T = dep[5,:]
    P = dep[6,:]

    # return all profiles
    return t, nA, nB, nD, nZ, nU, T, P

# function that performs the analysis
def perform_the_analysis():
    # set a rango of reaction times
    times = np.linspace(1,60.0,1000)

    # create storage for the yields
    yields = np.zeros(1000)

    # calculate the yield at each reaction time
    for i in range(0, 1000):
        t, nA, nB, nD, nZ, nU, T, P = profiles(times[i])
        yields[i] = nD[-1]/nA_0

    # find the reaction time where the yield is maximum
    i_max = np.argmax(yields)
    t_max = times[i_max]

    # solve the reactor design equations at the optimum reaction time
    _, nA, _, nD, _, _, _, _ = profiles(t_max)

    # calculate the other quantities of interest
    yield_max = nD[-1]/nA_0
    fA_pct = 100*(nA_0 - nA[-1])/nA_0

    # tabulate the results
    data = [['t_max', f'{t_max}', 'min'],
    ['Y_D/A', f'{yield_max}', 'mol D per initial mol A'],
    ['Conversion',f'{fA_pct}', '%']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(f'Optimum reaction time: {t_max:.3g} min')   
    print(f'Yield: {yield_max:.3g} mol D per initial mol A')
    print(f'Conversion of A: {fA_pct:.3g} %')
    print(' ')

    # save the results
    results_df.to_csv('reb_9_6_3/python/results.csv'
                      , index=False)
    
    # plot the yield vs. the reaction time
    plt.figure()
    plt.plot(times, yields)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Yield (mol D per initial mol A)")

    # show and save the graph
    plt.savefig('reb_9_6_3/python/yield_vs_t.png')
    plt.show()
    return

if __name__=="__main__":
    perform_the_analysis()
