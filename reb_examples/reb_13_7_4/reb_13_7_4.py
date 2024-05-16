"""Calculations for Example 13.7.4 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
dH = -9120 # cal /mol
K0 = 0.132
Cp_CO = 29.3/4.184 # cal /mol /K
Cp_H2O = 34.3/4.184 # cal /mol /K
Cp_CO2 = 41.3/4.184 # cal /mol /K
Cp_H2 = 29.2/4.184 # cal /mol /K
Cp_I = 40.5/4.184 # cal /mol /K
k0 = 0.0354 # mol /cm^3 /min /atm^2
E = 9740 # cal /mol
P = 26 # atm
Tin = 320 + 273.15 # K
nCO_in = 1 # mol /h
nCO2_in = 0.359 # mol /h
nH2_in = 4.44 # mol /h
nI_in = 0.18 # mol /h
fCO = 0.55
# known
Re = 1.987 # cal /mol /K

# make the inlet H2O flow available to all functions
nH2O_in = float('nan')

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nCO = dep[0]
    nH2O = dep[1]
    nCO2 = dep[2]
    nH2 = dep[3]
    nI = dep[4]
    T = dep[5]

	# calculate the rate
    k = k0*np.exp(-E/Re/T)
    K = K0*np.exp(-dH/Re/T)
    n = nCO + nH2O + nCO2 + nH2 + nI
    P_CO = nCO/n*P
    P_H2O = nH2O/n*P
    P_CO2 = nCO2/n*P
    P_H2 = nH2/n*P
    r = k*(P_CO*P_H2O - P_CO2*P_H2/K)

	# evaluate the derivatives
    dnCOdz = -r
    dnH2Odz = -r
    dnCO2dz = r
    dnH2dz = r
    dnIdz = 0.0
    dTdz = -r*dH/(nCO*Cp_CO + nH2O*Cp_H2O + nCO2*Cp_CO2 + nH2*Cp_H2 + nI*Cp_I)

	# return the derivatives
    return dnCOdz, dnH2Odz, dnCO2dz, dnH2dz, dnIdz, dTdz

# reactor model
def profiles():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nCO_in, nH2O_in, nCO2_in, nH2_in, nI_in, Tin])

	# define the stopping criterion
    f_var = 1
    f_val = nCO_in*(1-fCO)
     
	# solve the IVODEs
    V, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nCO = dep[0,:]
    nH2O = dep[1,:]
    nCO2 = dep[2,:]
    nH2 = dep[3,:]
    nI = dep[4,:]
    T = dep[5,:]

    # return all profiles
    return V, nCO, nH2O, nCO2, nH2, nI, T

# perform the analysis
def perform_the_analysis():
    # allow this function to set the inlet H2O flow
    global nH2O_in

    # set a range of inlet H2O flow rates
    nH2O_range = np.linspace(3*nCO_in, 6*nCO_in, 100)

    # allocate storage for the corresponding volumes and outlet temperatures
    V_range = np.zeros(100)
    T_range = np.zeros(100)

    # solve the reactor design equations for each inlet H2O flow
    for i in range(0,100):
        # set the inlet H2O flow rate
        nH2O_in = nH2O_range[i]

        # solve the reactor design equations
        V, nCO, nH2O, nCO2, nH2, nI, T = profiles()

        # save the volume and temperature
        V_range[i] = V[-1]
        T_range[i] = T[-1] - 273.15

    # find the index of the minimum volume
    i_min = np.argmin(V_range)

    # get the optimum H2O flow, minimum volume and temperature
    nH2O_opt = nH2O_range[i_min]
    Vmin = V_range[i_min]
    Tout = T_range[i_min]

    # tabulate the results
    data =[["Optimum H2O Feed Rate",nH2O_opt,"mol h^-1^"]
           ,["Minimum Volume",Vmin,"cm^3^"]
           ,["Outlet Temperature",Tout,"°C"]]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)

    # save the results
    results_df.to_csv('reb_13_7_4/python/results.csv', index=False)

    # display and save the graphs
    plt.figure(1) 
    plt.plot(nH2O_range,V_range)
    plt.xlabel("H$_2$O Feed Rate (mol min$^{-1}$)")
    plt.ylabel("Reactor Volume (cm$^3$)")
    plt.savefig('reb_13_7_4/python/V_vs_H2O_in.png')
    plt.show()

    plt.figure(2) 
    plt.plot(nH2O_range,T_range)
    plt.xlabel("H$_2$O Feed Rate (mol min$^{-1}$)")
    plt.ylabel("Outlet Temperature (°C)")
    plt.savefig('reb_13_7_4/python/T_vs_H2O_in.png')
    plt.show()
    return

if __name__=="__main__":
    perform_the_analysis()
