"""Calculations for Reaction Engineering Basics Example 9.6.4"""

# import libraries
import math
import numpy as np
from score_utils import solve_ivodes
import pandas as pd
import matplotlib.pyplot as plt

# given, known, and calculated constants available in all functions
# given
k_0_1 = 2.59E9 # /min
E_1 = 16500. # cal /mol
dH_1 = -22200. # cal /mol
CA_0 = 2. # mol /L
T_0 = 23  + 273.15 # K
Cp = 440. # cal /L /K
V = 4.0 # L
V_shell = 0.5 # L
A_shell = 0.6 # ft^2
U_shell = 1.13E4/60 # cal /ft^2 /min /K
Tex_in = 20 + 273.15 # K
rho_ex = 1.0 # g /cm^3
Cp_ex = 1.0 # cal /g /K
U_coil = 3.8E4/60 # cal /ft^2 /min /K
A_coil = 0.23 # ft^2
T_coil = 120 + 273.15 # K
Tex_0 = 23 + 273.15 # K
T_1_f = 50 + 273.15 # K
T_f = 25 + 273.15 # K
t_turn = 25 # min
# known
Re = 1.987 # cal /mol /K
# calculated
nA_0 = CA_0*V

# make the current coolant flow rate available to all functions
mDot_ex = float('nan')

# derivatives function for the first stage of operation
def first_stage_derivatives(ind, dep):
	# extract necessary dependent variables for this integration step
    nA = dep[0]
    T = dep[2]
    Tex = dep[3]

	# calculate the rate
    CA = nA/V
    k_1 = k_0_1*math.exp(-E_1/Re/T)
    r_1 = k_1*CA
    
    # calculate the rates of heat exchange
    Qdot_shell = U_shell*A_shell*(Tex - T)
    Qdot_coil = U_coil*A_coil*(T_coil - T)

	# evaluate the derivatives
    dnAdt = -V*r_1
    dnZdt = V*r_1
    dTdt = (Qdot_shell + Qdot_coil - V*r_1*dH_1)/V/Cp
    dTexdt = -Qdot_shell/rho_ex/V_shell/Cp_ex
    
	# return the derivatives
    return [dnAdt, dnZdt, dTdt, dTexdt]

# derivatives function for the second stage of operation
def second_stage_derivatives(ind, dep):
	# extract necessary dependent variables for this integration step
    nA = dep[0]
    T = dep[2]
    Tex = dep[3]

	# calculate the rate
    CA = nA/V
    k_1 = k_0_1*math.exp(-E_1/Re/T)
    r_1 = k_1*CA

    # calculate the rates of heat exchange
    Qdot_shell = U_shell*A_shell*(Tex - T)

	# evaluate the derivatives
    dnAdt = -V*r_1
    dnZdt = V*r_1
    dTdt = (Qdot_shell - V*r_1*dH_1)/V/Cp
    dTexdt = -(Qdot_shell + mDot_ex*Cp_ex*(Tex-Tex_in))/rho_ex/V_shell/Cp_ex
    
	# return the derivatives
    return [dnAdt, dnZdt, dTdt, dTexdt]

# reactor model for the first stage of operation
def first_stage_profiles():
	# set the initial values
    ind_0 = 0.0
    dep_0 = np.array([nA_0, 0.0, T_0, Tex_0])

	# define the stopping criterion
    f_var = 3
    f_val = T_1_f
    
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , first_stage_derivatives)
    
    # check that a solution was found
    if not(success):
        print(f"First stage IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]
    Tex = dep[3,:]

    # return all profiles
    return t, nA, nZ, T, Tex

# reactor model for the second stage of operation
def second_stage_profiles(t_0, nA_0, nZ_0, T_0, Tex_0):
	# set the initial values
    ind_0 = t_0
    dep_0 = np.array([nA_0, nZ_0, T_0, Tex_0])

	# define the stopping criterion
    f_var = 3
    f_val = T_f
     
	# solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , second_stage_derivatives)

    # check that a solution was found
    if not(success):
        print(f"Second stage IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    nA = dep[0,:]
    nZ = dep[1,:]
    T = dep[2,:]
    Tex = dep[3,:]

    # return all profiles
    return t, nA, nZ, T, Tex

# function that performs the analysis
def perform_the_analysis():
    # allow this function to set mDot_ex
    global mDot_ex

    # choose a range of coolant flow rates
    coolant_flows = np.linspace(100.0, 250.0, 100)

    # calculate the net rate for each coolant flow rate
    net_rates = np.zeros(100)
    for i in range(0,100):
        # make the coolant flow rate available
        mDot_ex = coolant_flows[i]

        # solve the reactor design equations
        t1, nA1, nZ1, T1, Tex1 = first_stage_profiles()
        t, _, nZ, _, _ = second_stage_profiles(t1[-1], nA1[-1], nZ1[-1]
                                                , T1[-1], Tex1[-1])
        
        # calculate the net rate
        net_rates[i] = nZ[-1]/(t[-1] + t_turn)

    # find the coolant flow where the net rate is maximized
    i_max = np.argmax(net_rates)
    mDot_max = coolant_flows[i_max]

    # solve the reactor design equations using the optimum coolant flow
    mDot_ex = mDot_max
    t_1, nA_1, nZ_1, T_1, Tex_1 = first_stage_profiles()
    t_2, nA_2, nZ_2, T_2, _ = second_stage_profiles(t_1[-1], nA_1[-1], nZ_1[-1]
                                                , T_1[-1], Tex_1[-1])
    
    # combine the profiles
    t = np.concatenate((t_1, t_2))
    nA = np.concatenate((nA_1, nA_2))
    nZ = np.concatenate((nZ_1, nZ_2))
    T = np.concatenate((T_1, T_2))
    
    # calculate the conversion vs time at the optimum coolant flow rate
    pct_conversion = 100*(nA_0 - nA)/nA_0
    
    # tabulate the results
    max_net_rate = nZ[-1]/(t[-1] + t_turn)
    data = [["Optimum Coolant Flow", f"{mDot_max}", "g min^-1^"],
            ["Maximum Net Rate", f"{max_net_rate}", "mol min^-1^"]]
    results_df = pd.DataFrame(data, columns=['item','value','units'])

    # display the results
    print(" ")
    print(f"Optimum Coolant Flow: {mDot_max} g /min")
    print(f"Maximum Net Rate: {max_net_rate} mol /min")

    # save the results
    results_df.to_csv('reb_9_6_4/python/results.csv'
                      , index=False)

    # display and save the graphs
    plt.figure() # net rate vs coolant flow
    plt.plot(coolant_flows, net_rates)
    plt.xlabel("Coolant Flow (g/min)")
    plt.ylabel("Net Rate (mol/min)")
    plt.savefig(
        'reb_9_6_4/python/net_rate_vs_coolant_flow.png')
    plt.show()

    plt.figure() # conversion profile
    plt.plot(t,pct_conversion)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Conversion (%)")
    plt.savefig('reb_9_6_4/python/conversion_profile.png')
    plt.show()

    plt.figure() # temperature profile
    plt.plot(t,T - 273.15)
    plt.xlabel("Reaction Time (min)")
    plt.ylabel("Temperature (°C)")
    plt.savefig('reb_9_6_4/python/temperature_profile.png')
    plt.show()

    return

if __name__=="__main__":
    perform_the_analysis()
