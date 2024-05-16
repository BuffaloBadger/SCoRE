"""Calculations for Example 13.7.6 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# constants available to all functions
# given
T = 400 + 273.15 # K
k = 0.15 # /s
eps = 0.4
Vdot_in = 1.0 # cm^3 /s
mDot = 0.054 # g /s
P_in = 30.0 # atm
Rpart = 0.2 # cm
Deff = 0.0029 # cm^2 /s
mu = 0.022E-2 # g /cm /s
D = 1.5 # cm
fA = 0.85
sph = 1.0
# known
R = 82.06 # cm^3 atm /mol /K
Pconv = 9.872E-7 # atm cm^2 /dyne
# calculated
nA_in = Vdot_in*P_in/R/T
nA_out = nA_in*(1-fA)
G = 4*mDot/np.pi/D**2
phi = Rpart*np.sqrt(k/Deff)
eta = 3/phi*(1/np.tanh(phi) - 1/phi)

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nZ = dep[1]
    P = dep[2]

	# calculate the rate
    CA = nA/(nA + nZ)*P/R/T
    r = k*CA

    # calculate the density
    rho = mDot*P/(nA + nZ)/R/T

	# evaluate the derivatives
    dnAdz = -np.pi*D**2/4*(1-eps)*eta*r
    dnZdz = np.pi*D**2/4*(1-eps)*eta*r
    dPdz = -Pconv*(1-eps)/eps**3*G**2/rho/sph/2/Rpart \
        *(150*(1-eps)*mu/sph/2/Rpart/G + 1.75)
    
	# return the derivatives
    return dnAdz, dnZdz, dPdz

# reactor model
def profiles():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_in, 0.0, P_in])

	# define the stopping criterion
    f_var = 1
    f_val = nA_out
     
	# solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    P = dep[2,:]

    # return all profiles
    return z, nA, nZ, P

# perform the analysis
def perform_the_analysis():
    # solve the reactor design equations
    z, nA, nZ, P = profiles()

    # calculate the other quantities of interest
    L = z[-1]
    P_out = P[-1]

    # tabulate the results
    data = [['effectiveness factor', eta, '']
            ,['reactor length', L, 'cm']
            ,['outlet pressure', P_out, 'atm']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])
    
    # display the results
    print(' ')
    print(results_df)

    # save the results
    results_df.to_csv('reb_13_7_6/python/results.csv', index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
