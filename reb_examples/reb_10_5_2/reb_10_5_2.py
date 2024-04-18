"""Calculations for Example 10.5.2 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import pandas as pd
from score_utils import solve_ivodes

# given, known, and calculated constants available to all functions
# given
k0_1 = 1.83E12 # L/mol/h
E_1 = 18000.0 # cal/mol
k0_2 = 5.08E13 # L/mol/h
E_2 = 20500.0 # cal/mol
dH_1 = -9000.0 # cal/mol
dH_2 = -7800.0 # cal/mol
Cp = 863.0# cal/L/K
CA_0 = 2.0 # mol/L
CB_in = 0.5 # mol/L
T_0 = 40 + 273.15 # K
T_in = T_0
P = 1.0 # atm
V_0 = 2000.0 # L
V_B = 8000.0 # L
t_f = 8.0 # h
# known
Re = 1.987 # cal/mol/K
Rw = 0.08206 # L-atm/mol/K
# calculated
nA_0 = CA_0*V_0

# make t_add available to all functions
t_add = float('nan')

# derivatives function
def derivatives(t, dep):
    # extract necessary dependent variables for this integration step
    nA = dep[0]
    nB = dep[1]
    T = dep[4]

    # calculate the reacting fluid volume
    if (t<t_add):
        Vdot_in = V_B/t_add
        V = Vdot_in*t + V_0
    else:
        Vdot_in = 0.0
        V = V_0 + V_B

    # calculate the inlet molar flow of B
    nDotB_in = Vdot_in*CB_in
    
    # calculate the rate coefficients
    k_1 = k0_1*np.exp(-E_1/Re/T)
    k_2 = k0_2*np.exp(-E_2/Re/T)

    # calculate the rate
    CA = nA/V
    CB = nB/V
    r_1 = k_1*CA*CB
    r_2 = k_2*CB**2

    # evaluate the derivatives
    dnAdt = -V*r_1
    dnBdt = nDotB_in - V*(r_1 + 2*r_2)
    dnDdt = V*r_1
    dnUdt = V*r_2
    dTdt = (-Vdot_in*Cp*(T - T_in) - V*r_1*dH_1 - V*r_2*dH_2 \
            + P*Vdot_in*Re/Rw)/(V*Cp)

    # return the derivatives
    return [dnAdt, dnBdt, dnDdt, dnUdt, dTdt]

# reactor model function
def profiles():
    # set the initial values for stage 1
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0, 0.0, 0.0, T_0])

    # define the stopping criterion
    f_var = 0
    f_val = t_add
     
    # solve the IVODEs
    t1, dep1, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, False)

    # check that a solution was found
    if not(success):
        print(f"A stage 1 IVODE solution was NOT obtained: {message}")
    
    # set the initial values for stage 2
    ind_0 = t1[-1]
    dep_0 = np.array([dep1[0,-1], dep1[1,-1], dep1[2,-1], dep1[3,-1]
                      , dep1[4,-1]])

    # define the stopping criterion
    f_val = t_f
     
    # solve the IVODEs for stage 2
    t2, dep2, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, False)

    # check that a solution was found
    if not(success):
        print(f"A stage 2 IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    t = np.concatenate((t1, t2))
    nA = np.concatenate((dep1[0,:], dep2[0,:]))
    nB = np.concatenate((dep1[1,:], dep2[1,:]))
    nD = np.concatenate((dep1[2,:], dep2[2,:]))
    nU = np.concatenate((dep1[3,:], dep2[3,:]))
    T = np.concatenate((dep1[4,:], dep2[4,:]))

    # return all profiles
    return t, nA, nB, nD, nU, T

# function that performs the analysis
def perform_the_analysis():
    # allow this function to set t_add
    global t_add

    # set the add times and allocate storage for the quantities of interest
    add_times = np.array([1.0, 3.0, 5.0, 7.0])
    f_B = np.zeros(4)
    S_DoverU = np.zeros(4)
    Y_DfromB = np.zeros(4)

    for iAdd in range(0,4):
        # set t_add
        t_add = add_times[iAdd]

        # solve the reactor design equations
        t, nA, nB, nD, nU, T = profiles()

        # calculate the quantities of interest
        f_B[iAdd] = (V_B*CB_in - nB[-1])/(V_B*CB_in)
        S_DoverU[iAdd] = nD[-1]/nU[-1]
        Y_DfromB[iAdd] = nD[-1]/(V_B*CB_in)


    # convert to percents
    f_B = 100.0*f_B
    Y_DfromB = 100.0*Y_DfromB
    
    # tabulate, display and save the results
    results_df = pd.DataFrame({'t_add':add_times , 'f_B':f_B
                               ,'S_DoverU':S_DoverU, 'Y_DfromB':Y_DfromB})
    print(results_df)
    results_df.to_csv("reb_10_5_2/python/results.csv", index=False)

    return

if __name__=="__main__":
    perform_the_analysis()