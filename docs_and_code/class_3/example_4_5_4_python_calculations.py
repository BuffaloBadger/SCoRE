""" Calculations for Example 4.5.4 of Reaction Engineering Basics"""

# import libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from score_utils import Arrhenius_parameters

# set the dpi for figures
plt.rc("savefig", dpi=300)

# global constants
R = 1.987E-3 # kcal/mol/K

# Read the data from the .csv file
df = pd.read_csv('example_4_5_4_data.csv')

# extract the data as arrays
T = df['T (°C)'].to_numpy() + 273.15  # convert to °C to K
k = df['k (L/mol/min)'].to_numpy()  # rate constant in L/mol/K

# calculate the Arrhenius expression parameters
k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k, T, R)

# print the results
print(" ")
print(f"k0 = {k0:.2e}, 95% CI [{k0_ci[0]:.2e}, {k0_ci[1]:.2e}] L/mol/min")
print(f"E = {E:.2f}, 95% CI [{E_ci[0]:.2f}, {E_ci[1]:.2f}] kcal/mol")
print(f"R^2 = {r_squared:.4f}")

# save the results to a .csv file
data = [['k0', f'{k0:.4g}', 'cm^3^ mol^-1^ min^-1^'],
        ['k0_lower_limit', f'{k0_ci[0]:.4g}', 'cm^3^ mol^-1^ min^-1^'],
        ['k0_upper_limit', f'{k0_ci[1]:.4g}', 'cm^3^ mol^-1^ min^-1^'],
        ['E', f'{E:.4g}', 'kJ mol^-1^'],
        ['E_lower_limit', f'{E_ci[0]:.4g}', 'kJ mol^-1^'],
        ['E_upper_limit', f'{E_ci[1]:.4g}', 'kJ mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
result = pd.DataFrame(data, columns=['item','value','units'])
result.to_csv("example_4_5_4_python_results.csv", index=False)

# create, show and save an Arrhenius plot
x=1/R/T
y_model = np.log(k0*np.exp(-E/R/T))
y = np.log(k)
plt.figure()
plt.plot(x, y, color = 'k', marker='o', ls=''
            , label = 'Experimental Data')
plt.plot(x,y_model, color = 'r', ls = '-', label = 'Arrhenius Expression')
plt.xlabel("$T^{-1}$ K$^{-1}$")
plt.ylabel("$\log{k}$")
plt.legend
plt.tight_layout()
plt.savefig('example_4_5_4_python_model_plot.png')
plt.show()