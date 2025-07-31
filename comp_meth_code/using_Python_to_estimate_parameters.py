""" Calculations for Solving IVODEs using Python """

# import libraries
import numpy as np
import pandas as pd
from score_utils import solve_ivodes
from score_utils import fit_to_SR_data
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi=300)

# global constants
L = 30.0 # cm
D = 2.0 # cm
R = 8.3145E-3 # kJ mol^-1^ K^-1^

# read the data from the .csv file
df = pd.read_csv('estimating_parameters_data.csv')

# extract the data as arrays
T = df['T'].to_numpy() + 273.15  # convert to °C to K
VFR = df['VFR'].to_numpy()
CA0 = df['CA0'].to_numpy()/1000.0  # convert to mol/cc
CB0 = df['CB0'].to_numpy()/1000.0  # convert to mol/cc
exptResp = df['fA'].to_numpy()

# combine the adjusted inputs as a matrix
adjExptInputs = np.transpose(np.array([T, VFR, CA0, CB0]))

# determine the number of experiments
nExpts = len(exptResp)

# globally available variable
global g_k0, E, iExpt
g_k0 = float('nan')
g_E = float('nan')
g_iExpt = -1

# derivatives function
def derivatives(z,dep):
    # extract the dependent variables
    nDot_A = dep[0]
    nDot_B = dep[1]
    nDot_Z = dep[2]

    # calculate the rate
    k=g_k0 * np.exp(-g_E/(R*T[g_iExpt]))
    CA = nDot_A/VFR[g_iExpt]
    CB = nDot_B/VFR[g_iExpt]
    r = k * CA * CB

    # evaluate the derivatives
    dnAdz = -np.pi * (D**2/4)*r
    dnBdz = -np.pi * (D**2/4)*r
    dnZdz = 2*np.pi * (D**2/4)*r

    # collect and return the derivatives
    ddz = np.array([dnAdz, dnBdz, dnZdz])
    return ddz

# PFR function
def PFR_variables(iExpt):
    # define initial values
    ind0 = 0
    nDotA0 = CA0[iExpt] * VFR[iExpt]  # mol/min
    nDotB0 = CB0[iExpt] * VFR[iExpt]  # mol/min
    nDotZ0 = 0.0  # mol/min
    dep0 =[nDotA0, nDotB0, nDotZ0]

    # define stopping criterion
    stop_var = 0.0
    stop_val = L

    # make iExpt available to the derivatives function
    global g_iExpt
    g_iExpt = iExpt

    # solve the PFR design equations
    z, dep, success, message = solve_ivodes(ind0, dep0, stop_var, stop_val
                                            , derivatives, False)

    # print a warning if there was a problem solving the design equations
    if not(success):
        print(f"An IVODE solution was NOT obtained: {message}")
    
    # extract the dependent variable profiles
    nDotA = dep[0,:]
    nDotB = dep[1,:]
    nDotZ = dep[2,:]

    # return the profiles
    return z, nDotA, nDotB, nDotZ

# predicted responses function
def predicted_responses(adjExptInputs, k0, E):
    # extract the parameters and make them globally available
    global g_k0, g_E
    g_k0 = k0
    g_E = E

    # allocate storage for the predicted responses
    modelResp = np.zeros(nExpts)

    # loop through the data points
    for i, input in enumerate(adjExptInputs):
        # Solve the PFR design equations
        z, nDotA, nDotB, nDotZ = PFR_variables(i)

        # Calculate the response
        nDotA0 = CA0[i] * VFR[i]
        nDotA_final = nDotA[-1]
        modelResp[i] = (nDotA0 - nDotA_final)/nDotA0
    
    # return the responses
    return modelResp

# parameter estimation function
def parameter_estimates():
    # define guesses for the parameters
    k0guess = 1.0E9
    Eguess = 30.0
    guess = [k0guess, Eguess]

    beta, beta_ci, r_squared = fit_to_SR_data(guess, adjExptInputs 
        , exptResp, predicted_responses, False)
    
    # extract the results
    k0 = beta[0]
    k0_CI = beta_ci[0,:]
    E = beta[1]
    E_CI = beta_ci[1,:]

    return k0, k0_CI, E, E_CI, r_squared

# deliverables function
def deliverables():
    # estimate the parameters
    k0, k0_CI, E, E_CI, r_squared = parameter_estimates()

    # save the results to a .csv file
    data = [['k0', f'{k0:.4g}', 'cm^3^ mol^-1^ min^-1^'],
        ['k0_lower_limit', f'{k0_CI[0]:.4g}', 'cm^3^ mol^-1^ min^-1^'],
        ['k0_upper_limit', f'{k0_CI[1]:.4g}', 'cm^3^ mol^-1^ min^-1^'],
        ['E', f'{E:.4g}', 'kJ mol^-1^'],
        ['E_lower_limit', f'{E_CI[0]:.4g}', 'kJ mol^-1^'],
        ['E_upper_limit', f'{E_CI[1]:.4g}', 'kJ mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    print(" ")
    print(result)
    result.to_csv("results.csv", index=False)

    # create, show, and save a parity plot
    global g_k0, g_E
    g_k0 = k0
    g_E = E
    modelResp = predicted_responses(adjExptInputs, g_k0, g_E)
    plt.figure() 
    plt.plot(exptResp, modelResp, color = 'k', marker='o', ls=''
             , label = 'Data')
    plt.plot([min(exptResp),max(exptResp)],[min(exptResp),max(exptResp)]
             , color = 'r', ls = '-', label = 'Parity Line')
    plt.xlabel("$f_A$")
    plt.ylabel("$f_{A,model}$")
    plt.legend
    plt.tight_layout()
    plt.savefig('parity.png')
    plt.show()

    # create, show and save residuals plots
    epsilon_expt = (modelResp - exptResp)
    plt.figure() 
    plt.plot(T, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$T$ (K)")
    plt.ylabel("Residual")
    plt.tight_layout()
    plt.savefig('residual_T.png')
    plt.show()

    plt.figure() 
    plt.plot(VFR, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$\dot{V}$ (cm$^3$ min$^{-1}$)")
    plt.ylabel("Residual")
    plt.tight_layout()
    plt.savefig('residual_Vdot.png')
    plt.show()

    plt.figure() 
    plt.plot(CA0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{A,0}$ (mol cm$^{-3}$)")
    plt.ylabel("Residual")
    plt.tight_layout()
    plt.savefig('residual_CA0.png')
    plt.show()

    plt.figure() 
    plt.plot(CB0, epsilon_expt, color = 'k', marker='o', ls='')
    plt.axhline(y=0, color = 'r')
    plt.xlabel("$C_{B,0}$ (mol cm$^{-3}$)")
    plt.ylabel("Residual")
    plt.tight_layout()
    plt.savefig('residual_CB0.png')
    plt.show()

if __name__=="__main__":
    deliverables()
