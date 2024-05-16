"""Calculations for Example 13.7.3 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import pandas as pd

# constants available to all functions
# given
yA_in = 0.76
yI_in = 0.25
Vdot_in = 100 # cm^3 /s
mDot_in = 0.44 # g /s
P_in = 3.0 # atm
T_in = 400 + 273.15 # K
D = 2.5 # cm
L = 800. # cm
T_ex = 375 + 273.15 # K
U = 187E-1/3600 # J /s /cm^2 /K
Dp = 0.25 # cm
phi = 0.7
eps = 0.6
k0f = 9E17 # mol /cm^3 /s /atm
Ef = 285E3 # J /mol
k0r = 4.09E-4 # mol /cm^3 /s /atm^4
Er = 85E3 # J /mol
dH = 200E3 # J /mol
CpA = 11.7*4.184 # J /mol /K
CpY = 8.3*4.184 # J /mol /K
CpZ = 4.2*4.184 # J /mol /K
CpI = 5.8*4.184 # J /mol /K
mu = 0.027E-2 # g /cm /s
# known
Re = 8.314 # J /mol /K
Rw = 82.06 # cm^3 atm /mol /K
P_conv = 9.872E-7 # atm cm^2 /dyne
# calculated
G = mDot_in/(np.pi*D**2/4)
nA_in = yA_in*Vdot_in*P_in/Rw/T_in
nI_in = yI_in*Vdot_in*P_in/Rw/T_in

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nY = dep[1]
    nZ = dep[2]
    nI = dep[3]
    T = dep[4]
    P = dep[5]

	# calculate the rate
    kf = k0f*np.exp(-Ef/Re/T)
    kr = k0r*np.exp(-Er/Re/T)
    ntot = nA + nY + nZ + nI
    PA = nA/ntot*P
    PY = nY/ntot*P
    PZ = nZ/ntot*P
    r = kf*PA - kr*PY*PZ**3

    # calculate the density
    Vdot = ntot*Rw*T/P
    rho = mDot_in/Vdot

	# evaluate the derivatives
    dnAdz = -np.pi*D**2/4*r
    dnYdz = np.pi*D**2/4*r
    dnZdz = 3*np.pi*D**2/4*r
    dnIdz = 0.0
    dTdz = (np.pi*D*U*(T_ex - T) - np.pi*D**2/4*r*dH)/(nA*CpA + nY*CpY 
                                                       + nZ*CpZ + nI*CpI)
    dPdz = -P_conv*(1-eps)/eps**3*G**2/rho/phi/Dp*(150*(1-eps)*mu/phi/Dp/G 
                                                   + 1.75)

	# return the derivatives
    return dnAdz, dnYdz, dnZdz, dnIdz, dTdz, dPdz

# reactor model
def profiles():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_in,0.0,0.0,nI_in,T_in,P_in])

	# define the stopping criterion
    f_var = 0
    f_val = L
     
	# solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nY = dep[1,:]
    nZ = dep[2,:]
    nI = dep[3,:]
    T = dep[4,:]
    P = dep[5,:]

    # return all profiles
    return z,nA,nY,nZ,nI,T,P

# perform the analysis
def perform_the_analysis():
    # solve the reactor design equations
    z, nA, nY, nZ, nI, T, P = profiles()

    # calculate the other quantities of interest
    T_out = T[-1] - 273.15
    P_out = P[-1]
    fA = 100*(nA_in - nA[-1])/nA_in

    # tabulate the results
    data =[['Outlet T',T_out,'°C'], ['Outlet P',P_out,'atm']
           ,['Conversion',fA,'%']]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(' ')
    print(results_df)

    # save the results
    results_df.to_csv('reb_13_7_3/python/results.csv',index=False)
    return

if __name__=="__main__":
    perform_the_analysis()
