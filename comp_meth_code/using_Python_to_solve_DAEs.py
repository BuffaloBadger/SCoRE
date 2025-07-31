""" Calculations for Using Python to Solve DAEs """

# import libraries
import numpy as np
from score_utils import solve_ivodes
import scipy as sp
import pandas as pd

# global constants
nA0 = 2400 # mol
T0 = 300 # K
nZ0 = 0 # mol
k0 = 1.44E10 # /min
E = 15.3 # kcal /mol
dH = -7 # kcal /mol
Cp = 1 # kcal /L /K
tf = 9 # min
nAf = 480 # mol
R = 1.987E-3 # kcal /mol /K

# global variable
g_V = float('nan')

# derivatives function
def IVODE_derivatives(ind, dep):
    # extract variables needed to evaluate the derivatives
    nA = dep[0]
    nZ = dep[1]
    T = dep[2]

    # calculate r
    k = k0*np.exp(-E/R/T)
    CA = nA/g_V
    r = k*CA

    # evaluate the derivatives
    dnAdt = -r*g_V
    dnZdt = r*g_V
    dTdt = -r*dH/Cp

    # return an array containing the derivatives
    ddt = np.array([dnAdt, dnZdt, dTdt])
    return ddt

# residuals function
def ATE_residuals(V):
    # set the initial values and stopping criterion for the IVODEs
    ind_0 = 0
    dep_0 = np.array([nA0, nZ0, T0])
    stop_var = 0
    stop_val = tf

    # make V available to the derivatives function
    global g_V
    g_V = V[0]

    # solve the IVODEs
    odes_are_stiff = False
    ind, dep, success, message = solve_ivodes(ind_0, dep_0, stop_var, stop_val
        , IVODE_derivatives, odes_are_stiff)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print("WARNING: The ODE solution may not be accurate!")
        print(f"         {message}")

    # extract the final value of nA
    nA = dep[0, :]
    nAfinal = nA[-1]

    # evaluate and return the residual
    epsilon = nAfinal - nAf
    return epsilon

# BSTR model function
def BSTR_variables():
    # guess the fluid volume
    V_guess = 2400 # L

    # calculate the fluid volume
    soln = sp.optimize.root(ATE_residuals,V_guess)

    # check that the solution is converged
    if not(soln.success):
        print("")
        print(f"The solver did NOT converge: {soln.message}")

    # extract the fluid volume
    V = soln.x[0]

    # set the initial values and stopping criterion for the IVODEs
    ind_0 = 0
    dep_0 = np.array([nA0, nZ0, T0])
    stop_var = 0
    stop_val = tf

    # make V available to the derivatives function
    global g_V
    g_V = V

    # solve the IVODEs
    odes_are_stiff = False
    ind, dep, success, message = solve_ivodes(ind_0, dep_0, stop_var, stop_val
        , IVODE_derivatives, odes_are_stiff)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print("WARNING: The ODE solution may not be accurate!")
        print(f"         {message}")

    # extract the profiles
    t = ind
    nA = dep[0, :]
    nZ = dep[1, :]
    T = dep[2, :]

    # return the BSTR variables
    return V, t, nA, nZ, T

# deliverable function
def deliverables():
    # calculate the BSTR variables
    V, t, nA, nZ, T = BSTR_variables()

    # extract the final values
    t_final = t[-1]
    nA_final = nA[-1]
    T_final = T[-1]

    # tabulate the results
    data = [["Fluid volume (L)", f"{V:.2f}","L"],
            ["Final time (min)", f"{t_final:.2f}","min"],
            ["Final nA (mol)", f"{nA_final:.2f}","mol"],
            ["Final T (K)", f"{T_final:.2f}","K"]]
    results_df = pd.DataFrame(data, columns=["Item", "Value", "Units"])

    # print the results
    print(" ")
    print(results_df)

    # save the results to a CSV file
    results_df.to_csv("results.csv", index=False)

# execution command
if __name__ == "__main__":
    deliverables()