"""Calculations for Example 12.7.3 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
V = 0.5 # m^3
nA_in = 70 # mol /s
nB_in = 1500 # mol /s
Vdot_in = 40E-3 # m^3 /s
k0 = 1.2E9 # m^3 /mol /s
E = 25800*4.184 # J /mol
K0 = 4.2E-18 # m^3 /mol
dH = -22400*4.184 # J /mol
Cp_A = 412 # J /mol /K
Cp_B = 75.5 # J /mol /K
# known
R = 8.314 # J /mol /K

# make T_in available to all functions
T_in = -100.0

# reactor model
def unknowns(init_guess):
     
	# solve the ATEs
    soln = sp.optimize.root(residuals,init_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The initial temperature was NOT found: {soln.message}")

    # return the solution
    return soln.x

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nA = guess[0]
    nB = guess[1]
    nZ = guess[2]
    T = guess[3]

    # calculate the rate
    k = k0*np.exp(-E/R/T)
    K = K0*np.exp(-dH/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    CZ = nZ/Vdot_in
    r = k*CA*CB*(1 - CZ/(K*CA*CB))

    # evaluate the residuals
    residual_1 = nA_in - nA - V*r
    residual_2 = nB_in - nB - V*r
    residual_3 = -nZ + V*r
    residual_4 = -(nA_in*Cp_A + nB_in*Cp_B)*(T-T_in) - V*r*dH

    # return the residuals
    return np.array([residual_1, residual_2, residual_3, residual_4])

# perform the analysis
def perform_the_analysis():
    # allow this function to set T_in
    global T_in

    # set a range for T_in
    T_in_range = np.linspace(75.,125.,100) + 273.15

	# set the initial guess
    init_guess = np.array([0.5*nA_in, nB_in, 0.5*nA_in, T_in_range[0] + 10.0])

    # allocate storage for conversions
    fA_range = np.zeros(100)

    # calculate fA for each T_in
    for iT in range(0,100):

        # set T_in
        T_in = T_in_range[iT]

        # solve the reactor design equations
        solution = unknowns(init_guess)

        # calculate the conversion
        fA_range[iT] = 100*(nA_in - solution[0])/nA_in
        init_guess = solution

    # plot the results
    plt.figure(1) 
    plt.plot(T_in_range-273.15, fA_range)
    plt.xlabel("Inlet Temperature (°C)")
    plt.ylabel("Conversion of A (%)")

    # save and show the figure
    plt.savefig('reb_12_7_3/python/fA_vs_Tin.png')
    plt.show()

    # find the optimum T_in
    i_max = np.argmax(fA_range)

    # solve the reactor design equations using the optimum T_in
    T_in = T_in_range[i_max]
    fA_guess = fA_range[i_max]
    init_guess = np.array([(1-fA_guess)*nA_in, nB_in, fA_guess*nA_in, T_in])
    solution = unknowns(init_guess)

    # calculate the quantities of interest
    fA = 100*(nA_in - solution[0])/nA_in
    T = solution[3] - 273.15

    # tabulate the results
    data = [['opt T_in',f'{T_in - 273.15}','°C'],['fA',f'{fA}','%']
         ,['T out',f'{T}','°C']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('reb_12_7_3/python/results.csv',index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
