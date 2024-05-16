"""Calculations for Example 13.7.5 of Reaction Engineering Basics"""

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

# constants available to all functions
# given
D = 10 # cm
L = 500 # cm
Vdot_in = 75E3 # cm^3 /min
T_before = 25 + 273.15 # K
CA_in = 1E-3 # mol /cm^3
CB_in = 1.2E-3 # mol /cm^3
T_in = 30 + 273.15 # K
k0_1 = 8.72E8 # cm^3 /mol /min
E_1 = 7200 # cal /mol
dH_1 = -10700 # cal /mol
Cp = 1.0 # cal /g /K
rho = 1.0 # g /cm^3
t = np.array([0.131, 0.262, 0.393, 0.524])
# known
R = 1.987 # cal /mol /K
# calculated
nA_in = Vdot_in*CA_in
nB_in = Vdot_in*CB_in

# derivatives function
def derivatives(ind, dep):
	# extract the dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    T = dep[4]

	# calculate the rate
    k_1 = k0_1*np.exp(-E_1/R/T)
    CA = nA/Vdot_in
    CB = nB/Vdot_in
    r_1 = k_1*CA*CB

	# evaluate the derivatives
    dnAdz = -np.pi*D**2/4*r_1
    dnBdz = -np.pi*D**2/4*r_1
    dnYdz = np.pi*D**2/4*r_1
    dnZdz = np.pi*D**2/4*r_1
    dTdz = -np.pi*D**2/4*r_1*dH_1/(Vdot_in*rho*Cp)

	# return the derivatives
    return dnAdz, dnBdz, dnYdz, dnZdz, dTdz

# reactor model
def profiles(z_front):
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_in, nB_in, 0.0, 0.0, T_in])

	# define the stopping criterion
    f_var = 0
    f_val = z_front
     
	# solve the IVODEs
    z, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nB = dep[1,:]
    nY = dep[2,:]
    nZ = dep[3,:]
    T = dep[4,:]
    # return all profiles
    return z, nA, nB, nY, nZ, T

# perform the analysis
def perform_the_analysis():
    # analyze each of the four specified elapsed times
    for i in range(0,4):
        # calculate the position of the front
        V_front = Vdot_in*t[i]
        z_front = V_front/(np.pi*D**2/4)

        # get the profiles up to the front
        z, nA, nB, nY, nZ, T = profiles(z_front)

        if (z[-1] < L):
            # add the front after the front
            z = np.concatenate((z,np.array([z_front, L])))
            nA = np.concatenate((nA,np.array([0.0, 0.0])))
            T = np.concatenate((T,np.array([T_before, T_before])))

        # plot, display and save the fronts
        fig, ax1 = plt.subplots()

        color = 'tab:blue'
        ax1.set_xlabel('Axial Distance from Inlet (cm)')
        ax1.set_ylabel('Flow of A (mol min$^{-1}$)', color=color)
        ax1.plot(z, nA, color=color)
        ax1.set_ylim(-5.0,77.5)

        ax2 = ax1.twinx()

        color = 'tab:red'
        ax2.set_ylabel('Outlet Temperature °C', color=color)
        ax2.plot(z, T-273.15, color=color)
        ax2.set_ylim(23.0, 41.0)

        fig.tight_layout()
        plt.savefig(f'reb_13_7_5/python/t_{1000*t[i]:.0f}_min_profiles.png')
        plt.show()
    return

if __name__=="__main__":
    perform_the_analysis()
