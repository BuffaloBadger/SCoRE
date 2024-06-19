"""Calculations for Example X.Y.Z of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# constants available to all functions
# given
V_pfr1 = 60.0 # L
V_pfr2 = 40.0 # L
CA_0 = 1.0 # mol /L
T_0 = 60 + 273.15 # K
Vdot_0 = 0.55 # L /min
k0 = 2.63E7 # L /mol /min
E = 62000 # J /mol
dH = -35000 # J /mol
Cp = 800 # J /L /K
# known
R = 8.314 # J /mol /K
# calculated
nDotA_0 = Vdot_0*CA_0

# make Vdot available to all functions
Vdot = float('nan')

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nDot_A = dep[0]
    T = dep[3]
    
	# calculate the rate
    r = k0*np.exp(-E/R/T)*nDot_A**2/Vdot**2

	# evaluate the derivatives
    dnAdV = -2*r
    dnYdV = r
    dnZdV = r
    dTdV = -r*dH/Vdot/Cp

	# return the derivatives
    return dnAdV, dnYdV, dnZdV, dTdV

# reactor model
def profiles(f_val):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([Vdot*CA_0, 0.0, 0.0, T_0])

	# define the stopping criterion
    f_var = 0
     
	# solve the IVODEs
    V, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nDotA = dep[0,:]
    nDotY = dep[1,:]
    nDotZ = dep[2,:]
    T = dep[3,:]

    # return all profiles
    return V, nDotA, nDotY, nDotZ, T

# perform the analysis
def perform_the_analysis():
    # allow this function to set Vdot
    global Vdot

    # case a, equal flow rates
    Vdot = Vdot_0/2.0

    # solve the reactor design equations for reactor R1
    V, nDotA, nDotY, nDotZ, T = profiles(V_pfr1)
    nDotA_3 = nDotA[-1]

    # solve the reactor design equations for reactor R2
    V, nDotA, nDotY, nDotZ, T = profiles(V_pfr2)
    nDotA_4 = nDotA[-1]

    # calculate the conversion
    nDotA_5 = nDotA_3 + nDotA_4
    fA_a = 100.0*(nDotA_0 - nDotA_5)/nDotA_0

    # case b, equal space times
    Vdot_1 = Vdot_0*V_pfr1/(V_pfr1 + V_pfr2)
    Vdot_2 = Vdot_0 - Vdot_1 
    
    # solve the reactor design equations for reactor R1
    Vdot = Vdot_1;
    V, nDotA, nDotY, nDotZ, T = profiles(V_pfr1)
    nDotA_3 = nDotA[-1]

    # solve the reactor design equations for reactor R2
    Vdot = Vdot_2
    V, nDotA, nDotY, nDotZ, T = profiles(V_pfr2)
    nDotA_4 = nDotA[-1]

    # calculate the conversion
    nDotA_5 = nDotA_3 + nDotA_4
    fA_b = 100.0*(nDotA_0 - nDotA_5)/nDotA_0

    # tabulate the results
    data = [["Conversion for Equal Flow Rates",fA_a,'%']
            ,['Conversion for Equal Space Times',fA_b,'%']]
    results_df = pd.DataFrame(data,columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)

    # save the results
    results_df.to_csv('reb_15_4_4/python/results.csv',index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
