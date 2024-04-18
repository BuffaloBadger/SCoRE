"""Calculations for Example 10.5.3 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import pandas as pd
from score_utils import solve_ivodes
import matplotlib.pyplot as plt

# constants available to all functions
# given
k0_1 = 1.192E15 # /min
E_1 = 97600 # J/mol
dH_1 =-58615 # J/mol
V_ex = 300 # cm^3
Tex_in = 60 + 273.15 # K
Vdot_ex = 250 # cm^3/min
UA = 260*4.184 # J/min/K
P = 1.0 # atm
VW_0 = 67 # cm^3
VZ_0 = 283 # cm^3
Vacid = 0.3 # cm^3
T_0 = 60 + 273.15 # K
Tex_0 = 60 + 273.15 # K
Cp = 2.68 # J/cm^3/K
V_A = 350 # cm^3
T_in = 21 + 273.15 # K
T_f = 65 + 273.15 # K
T_max = 95 + 273.15 #K
rho_A = 1.082 # g/cm^3
rho_W = 1.0 # g/cm^3
rho_Z = 1.049 # g/cm^3
M_A = 102 # g/mol
M_W = 18 # g/mol
M_Z = 60 # g/mol
Cp_A = 168.2 # J/mol/K
rho_ex = 1.0 # g/cm^3
Cp_ex = 1.0*4.184 # cal/g/K
# known
Re = 1.987*4.184 # J/mol/K
Rw = 82.06 # cm^3-atm/mol/K
# calculated
V_0 = VW_0 + VZ_0 + Vacid
m_ex = Vdot_ex*rho_ex
nW_0 = VW_0*rho_W/M_W
nZ_0 = VZ_0*rho_Z*M_Z

# make Vdot_in available to all functions
Vdot_in = float('nan')

# derivatives function
def derivatives(t, dep):
    # extract necessary dependent variables for this integration step
    nA = dep[0]
    T = dep[3]
    Tex = dep[4]

    # calculate the batch processing time
    t_sb = V_A/Vdot_in

    # calculate the molar feed rate and reacting fluid volume
    if (t<t_sb):
        nA_in = Vdot_in*rho_A/M_A
        V = Vdot_in*t + V_0
        W_exp = P*Vdot_in*Re/Rw
    else:
        nA_in = 0.0
        V = V_0 + V_A
        W_exp = 0
    
    # calculate the rate coefficient
    k_1 = k0_1*np.exp(-E_1/Re/T)

    # calculate the rate
    CA = nA/V
    r_1 = k_1*CA

    # calculate the rate of heat exchange
    Qdot = UA*(Tex-T)

    # evaluate the derivatives
    dnAdt = nA_in - V*r_1
    dnWdt = - V*r_1
    dnZdt = V*r_1
    dTdt = (Qdot - nA_in*Cp_A*(T - T_in) - V*r_1*dH_1 + W_exp)\
        /(Cp*V)
    dTexdt = (-Qdot - m_ex*Cp_ex*(Tex-Tex_in))/(rho_ex*V_ex*Cp_ex)

    # return the derivatives
    return [dnAdt, dnWdt, dnZdt, dTdt, dTexdt]

# reactor model function
def profiles():
    # set the initial values for stage 1
    ind_0 = 0.0
    dep_0 = np.array([0.0, nW_0, nZ_0, T_0, Tex_0])

    # define the stopping criterion
    t_sb = V_A/Vdot_in
    f_var = 0
    f_val = t_sb
     
    # solve the IVODEs
    t, dep, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                        , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"A stage 1 IVODE solution was NOT obtained: {message}")
    
    nA = dep[0,:]
    nW = dep[1,:]
    nZ = dep[2,:]
    T = dep[3,:]
    Tex = dep[4,:]

    # set the initial values for stage 2
    ind_0 = t[-1]
    dep_0 = np.array([dep[0,-1], dep[1,-1], dep[2,-1], dep[3,-1]\
                        , dep[4,-1]])

    # define the stopping criterion
    f_var = 4
    f_val = T_f
    
    # solve the IVODEs for stage 2
    t2, dep2, success, message = solve_ivodes(ind_0, dep_0, f_var, f_val
                                    , derivatives, True)

    # check that a solution was found
    if not(success):
        print(f"A stage 2 IVODE solution was NOT obtained: {message}")

    # extract the dependent variable profiles
    t = np.concatenate((t, t2))
    nA = np.concatenate((nA, dep2[0,:]))
    nW = np.concatenate((nW, dep2[1,:]))
    nZ = np.concatenate((nZ, dep2[2,:]))
    T = np.concatenate((T, dep2[3,:]))
    Tex = np.concatenate((Tex, dep2[4,:]))

    # return all profiles
    return t, nA, nW, nZ, T, Tex

# function that performs the analysis
def perform_the_analysis():
    # allow this function to set Vdot_in
    global Vdot_in

    # set feed rates and allocate storage
    feed_rates = np.linspace(37,41,100)
    t_proc = 1E20*np.ones(100)

    for iFeed in range(0,100):
        # set Vdot_in
        Vdot_in = feed_rates[iFeed]

        # solve the reactor design equations
        t, nA, nW, nZ, T, Tex = profiles()

        # save the processing time if the temperature constraint is satisfied
        if np.max(T) <= T_max:
            t_proc[iFeed] = t[-1]

    # find the optimum
    t_opt = np.min(t_proc)
    i_opt = np.argmin(t_proc)
    Vdot_opt = feed_rates[i_opt]
    t_sb = V_A/Vdot_opt
    
    # solve the reactor design equations using the optimum feed rate
    Vdot_in = Vdot_opt
    t, nA, nW, nZ, T, Tex = profiles()
    f_A = 100*(V_A*rho_A/M_A - nA[-1])/(V_A*rho_A/M_A)

    # tabulate results
    data = [["Minimum Processing Time",f"{t_opt}","min"]
            ,["Semi-Batch Processing Time",f"{t_sb}","min"]
            ,["Optimum Feed Rate",f"{Vdot_opt}","cm^3^ min^-1^"]
            ,["Maximum Temperature",f"{np.max(T)-273.15}","°C"]
            ,["Maximum Cooling Water Temperature",f"{np.max(Tex)-273.15}","°C"]
            ,["Acetic Anhydride Conversion",f"{f_A}","%"]]
    results_df = pd.DataFrame(data, columns=['item','value','units'])
    print(results_df)
    results_df.to_csv('reb_10_5_3/python/results.csv',index=False)

    # create, display and save graphs
    plt.figure(1) 
    plt.plot(t, T-273.15)
    plt.xlabel("Time (min)")
    plt.ylabel("Temperature (°C)")
    plt.savefig('reb_10_5_3/python/temperature_profile.png')
    plt.show()

    # calculate the concentration profile
    CA = np.zeros(t.size)
    for i in range(0,t.size):
        if t[i] < t_sb:
            V = V_0 + Vdot_in*t[i]
        else:
            V = V_0 + V_A
        CA[i] = nA[i]/V

    plt.figure(2) 
    plt.plot(t, CA)
    plt.xlabel("Time (min)")
    plt.ylabel("Acetic Anhydride Concentration (mol cm$^{-3}$)")
    plt.savefig('reb_10_5_3/python/concentration_profile.png')
    plt.show()
    return

if __name__=="__main__":
    perform_the_analysis()