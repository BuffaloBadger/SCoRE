""" Calculations for Solving IVODEs using Python """

# import libraries
import numpy as np
from score_utils import solve_ivodes
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi=300)

# global constants
V = 1.0 # L
tf = 30.0 # min
CA0 = 0.5 # mol/L
CZ0 = 0.0 # mol/L
T = 338 # K
k0 = 3.6E8 # L/mol/min
E = 67.5 # kJ/mol
R = 8.314E-3 # kJ/mol/K

# derivatives function
def BSTR_derivatives(t, dep):
    # extract the dependent variables
    nA = dep[0]
    nZ = dep[1]

    # calculate the rate
    k = k0*np.exp(-E/R/T)
    CA = nA/V
    r = k*CA

    # calculate the derivatives of the dependent variables
    dnAdt = -r*V
    dnZdt = r*V

    # return an array containing the derivatives
    ddt = np.array([dnAdt, dnZdt])
    return ddt

# BSTR function
def BSTR_variables():
    # set the initial values and stopping criterion
    ind_0 = 0
    nA0 = CA0*V
    nZ0 = CZ0*V
    dep_0 = np.array([nA0, nZ0])
    stop_var = 0
    stop_val = tf

    # solve the BSTR reactor design equations
    odes_are_stiff = False
    ind, dep, success, message = solve_ivodes(ind_0, dep_0, stop_var, stop_val
            , BSTR_derivatives, odes_are_stiff)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print("WARNING: The ODE solution may not be accurate!")
        print(f"         {message}")
    
    # extract the corresponding sets of values for each of the variables
    t = ind
    nA = dep[0,:]
    nZ = dep[1,:]

    # return the solution
    return t, nA, nZ

# deliverables function
def deliverables():
    # get the BSTR variables
    t, nA, nZ = BSTR_variables()

    # plot, display, and save the results
    plt.figure()
    plt.plot(t,nA,label='A')
    plt.plot(t,nZ,label='Z')
    plt.xlabel('Time (min)')
    plt.ylabel('Amount (moles)')
    plt.legend()
    plt.savefig('profiles.png')
    plt.show()

if __name__=="__main__":
    deliverables()