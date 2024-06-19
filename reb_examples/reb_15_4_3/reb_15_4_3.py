"""Calculations for Example 15.4.3 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
T_0 = 30 + 273.15 # K
CA_0 = 1.0 # mol /l
CB_0 = 1.2 # mol /l
nDotY_0 = 0.0
nDotZ_0 = 0.0
Vdot = 75 # l /min
k0 = 8.72E5 # l /mol /min
E = 7200 # cal /mol
dH = -10700 # cal /mol
Cp = 1.0 # cal /g /K
rho = 1.0E3 # g /l
fA = 0.9
# known
R = 1.987 # cal /mol /K
# calculated
nDotA_0 = CA_0*Vdot
nDotB_0 = CB_0*Vdot
nDotA_2 = nDotA_0*(1-fA)

# make V_R2 available to all functions
V_R2 = -1

# reactor model
def unknowns(init_guess):
     
	# solve the ATEs
    soln = sp.optimize.root(residuals,init_guess)

    # check that the solution is converged
    if not(soln.success):
        print(f"The ATE solver did not converge: {soln.message}")

    # return the solution
    return soln.x

# residuals function
def residuals(guess):
    # extract the indiviaual guesses
    nDotA_1 = guess[0]
    nDotB_1 = guess[1]
    nDotY_1 = guess[2]
    nDotZ_1 = guess[3]
    T_1 = guess[4]
    V_R1 = guess[5]
    nDotB_2 = guess[6]
    nDotY_2 = guess[7]
    nDotZ_2 = guess[8]
    T_2 = guess[9]

    # calculate the rates
    r_R1 = k0*np.exp(-E/R/T_1)*nDotA_1*nDotB_1/Vdot**2
    r_R2 = k0*np.exp(-E/R/T_2)*nDotA_2*nDotB_2/Vdot**2

    # evaluate the residuals
    eps1 = nDotA_0 - nDotA_1 - V_R1*r_R1
    eps2 = nDotB_0 - nDotB_1 - V_R1*r_R1
    eps3 = nDotY_0 - nDotY_1 + V_R1*r_R1
    eps4 = nDotZ_0 - nDotZ_1 + V_R1*r_R1
    eps5 = rho*Vdot*Cp*(T_1 - T_0) + V_R1*r_R1*dH
    eps6 = nDotA_1 - nDotA_2 - V_R2*r_R2
    eps7 = nDotB_1 - nDotB_2 - V_R2*r_R2
    eps8 = nDotY_1 - nDotY_2 + V_R2*r_R2
    eps9 = nDotZ_1 - nDotZ_2 + V_R2*r_R2
    eps10 = rho*Vdot*Cp*(T_2 - T_1) + V_R2*r_R2*dH

    # return the residuals
    return np.array([eps1, eps2, eps3, eps4, eps5, eps6, eps7, eps8, eps9
            , eps10])

# perform the analysis
def perform_the_analysis():
    # allow this function to set V_R2
    global V_R2

    # choose a range of values for V_R2
    V_R2_range  = np.linspace(50.0,70.0,100)

    # set the initial guess
    init_guess = np.array([nDotA_0*0.11,
        nDotB_0 - nDotA_0*0.89,
        nDotA_0*0.89,
        nDotA_0*0.89,
        T_0 + 10.0,
        100.,
        nDotB_0 - nDotA_0*0.9,
        nDotA_0*0.9,
        nDotA_0*0.9,
        T_0 + 10.0])

    # allocate storage for the corresponding values of V_R1
    V_R1 = np.zeros(100)

    # calculate the corresponding values of V_R1
    for i in range(0,100):
        # set the volume of reactor R2
        V_R2 = V_R2_range[i]

        # solve the reactor design equations
        solution = unknowns(init_guess)

        # save the result
        V_R1[i] = solution[5]

        # use the result as the next initial guess
        init_guess = solution
        
    # plot the results
    plt.figure(1) 
    plt.plot(V_R2_range,V_R1 + V_R2_range)
    plt.xlabel("Reactor 2 Volume (l)")
    plt.ylabel("Total Reactor Volume (l)")

    # save and show the figure
    plt.savefig('reb_15_4_3/python/volume_plot.png')
    plt.show()

    # find the index of the minimum volume
    i_min = np.argmin(V_R1 + V_R2_range)

    # calculate the other quantities of interest
    opt_V_R1 = V_R1[i_min]
    opt_V_R2 = V_R2_range[i_min]
    min_total_V =  opt_V_R1 + opt_V_R2

    # tabulate the results
    data = [['Minimum Total Volume', min_total_V, 'L']
         ,['Optimum Reactor 1 Volume', opt_V_R1, 'L']
         ,['Optimum Reactor 2 Volume', opt_V_R2, 'L']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('reb_15_4_3/python/results.csv',index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
