"""Calculations for Example 15.4.2 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# constants available to all functions
# given
dH = -9120 # cal /mol
K0 = 0.132
CpCO = 29.3/4.184 # cal /mol /K
CpH2O = 34.3/4.184 # cal /mol /K
CpCO2 = 41.3/4.184 # cal /mol /K
CpH2 = 29.1/4.184 # cal /mol /K
CpI = 40.5/4.184 # cal /mol /K
nDotCO_in = 1.0 # mol /h
nDotCO2_in = 0.359 # mol /h
nDotH2_in = 4.44 # mol /h
nDotI_in = 0.18 # mol/h
nDotH2O_in = 9.32 # mol /h
P = 26.0 # atm
T_in = 445 + 273.15 # K
V_1 = 685 # cm^3
k0_1 = 3.54E-2 # mol /cm^3 /min /atm^2
E_1 = 9740 # cal /mol
Tex_in = 20 + 273.15 # K
mDotEx = 1100 # g /h
Tex_out = 50 + 273.15 # K
CpEx = 1.0 # cal /g /K
V_2 = 3950 # cm^3
k0_2 = 1.77E-3 # mol /cm^3 /min /atm^2
E_2 = 3690 # cal /mol
# known
R = 1.987 # cal /mol /K

# make the current reactor available to all functions
reactor_number = float('nan')

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nDotCO = dep[0]
    nDotH2O = dep[1]
    nDotCO2 = dep[2]
    nDotH2 = dep[3]
    nDotI = dep[4]
    T = dep[5]

	# calculate the rate
    nDotTotal = nDotCO + nDotH2O + nDotCO2 + nDotH2 + nDotI
    PCO = P*nDotCO/nDotTotal
    PH2O = P*nDotH2O/nDotTotal
    PCO2 = P*nDotCO2/nDotTotal
    PH2 = P*nDotH2/nDotTotal
    K = K0*np.exp(-dH/R/T)
    if reactor_number == 1:
        r = k0_1*np.exp(-E_1/R/T)*(PCO*PH2O - PCO2*PH2/K)
    else:
        r = k0_2*np.exp(-E_2/R/T)*(PCO*PH2O - PCO2*PH2/K)
    
	# evaluate the derivatives
    dCOdV = -r
    dH2OdV = -r
    dCO2dV = r
    dH2dV = r
    dIdV = 0
    dTdV = -r*dH/(nDotCO*CpCO + nDotH2O*CpH2O + nDotCO2*CpCO2 \
        + nDotH2*CpH2 + nDotI*CpI)
    
	# return the derivatives
    return dCOdV, dH2OdV, dCO2dV, dH2dV, dIdV, dTdV

# reactor model
def profiles(dep_0, f_val):
	# set the initial values
    ind_0 = 0.0

	# define the stopping criterion
    f_var = 0
     
	# solve the IVODEs
    V, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nDotCO = dep[0,:]
    nDotH2O = dep[1,:]
    nDotCO2 = dep[2,:]
    nDotH2 = dep[3,:]
    nDotI = dep[4,:]
    T = dep[5,:]

    # return all profiles
    return V, nDotCO, nDotH2O, nDotCO2, nDotH2, nDotI, T

# perform the analysis
def perform_the_analysis():
    # allow this function to set [missing constant]
    global reactor_number

    # reactor 1
    reactor_number = 1
    dep0 = np.array([nDotCO_in, nDotH2O_in, nDotCO2_in, nDotH2_in, nDotI_in,
        T_in])
    f_val = V_1

    # solve the reactor design equations
    V, nDotCO, nDotH2O, nDotCO2, nDotH2, nDotI, T = profiles(dep0, f_val)
    T_1 = T[-1] - 273.15
    fCO_1 = 100*(nDotCO_in - nDotCO[-1])/nDotCO_in

    # heat exchanger
    T_in_2 = T[-1] - mDotEx*CpEx*(Tex_out - Tex_in) \
        /(nDotCO[-1]*CpCO + nDotH2O[-1]*CpH2O \
        + nDotCO2[-1]*CpCO2 + nDotH2[-1]*CpH2 + nDotI[-1]*CpI)

    # reactor 2
    reactor_number = 2
    dep0 = np.array([nDotCO[-1], nDotH2O[-1], nDotCO2[-1], nDotH2[-1],
        nDotI[-1], T_in_2])
    f_val = V_2

    # solve the reactor design equations
    V, nDotCO, nDotH2O, nDotCO2, nDotH2, nDotI, T = profiles(dep0, f_val)

    # calculate the other quantities of interest
    T_in_2 = T_in_2 - 273.15
    T_2 = T[-1] - 273.15
    fCO_2 = 100*(nDotCO_in - nDotCO[-1])/nDotCO_in

    # tabulate the results
    data = [['Reactor 1 Outlet T',T_1,'°C']
            ,['Reactor 1 Conversion',fCO_1,'%']
            ,['Reactor 2 Inlet T',T_in_2,'°C']
            ,['Reactor 2 Outlet T',T_2,'°C']
            ,['Reactor 2 Conversion',fCO_2,'%']]
    results_df = pd.DataFrame(data,columns=['item','value','units'])

    # display the results
    print(results_df)

    # save the results
    results_df.to_csv('reb_15_4_2/python/results.csv',index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
