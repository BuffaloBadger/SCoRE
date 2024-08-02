"""Calculations for Reaction Engineering Basics Example 19.5.1 using differential analysis"""

# import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from score_utils import Arrhenius_parameters

# given and known constants
V = 1.0 # L
R = 8.314E-3 # kJ/mol/K

# function that performs the calculations
def perform_the_calculations():
    # read the data from the .csv file
    df = pd.read_csv('reb_19_5_1/reb_19_5_1_data.csv')
        # columns: Experiment, T, CA0, tf, CAf
    
    # get the temperatures of the data blocks
    block_temperatures = df['T'].unique()

    # create a dataframe for the fitting results from the same-temperature blocks
    parameter_estimates = pd.DataFrame(columns=['T', 'k', 'k_ll', 'k_ul', 'R_sq'])

    # process the data blocks
    for T in block_temperatures:
        # create the block
        block = df[df['T'] == T]

        # extract the data as arrays
        T_expts = block['T'].to_numpy()
        CA0 = block['CA0'].to_numpy()
        tf = block['tf'].to_numpy()
        CAf = block['CAf'].to_numpy()

        # calculate x
        x = block['CAf'].to_numpy()*V

        # allocate storage for y
        y = np.zeros(len(block.index))

        # calculate y
        for i in range(0,len(x)):
            if tf[i] == 5.0:
                x_minus = CA0[i]*V
                t_minus = 0
            else:
                x_minus = x[i-1]
                t_minus = tf[i-1]
            y[i] = (x_minus - x[i])/(tf[i] - t_minus)
        
        # fit y = mx to the data
        model = sm.OLS(y,x)
        res = model.fit()
        m = res.params[0]
        m_ci = res.conf_int(alpha=0.05)

        # add the fitting results to the dataframe
        parameter_estimates.loc[len(parameter_estimates)] = [T, m, m_ci[0,0]
                    , m_ci[0,1], res.rsquared]
    
        # calculate the model-predicted y
        y_model = m*x
        
        # create, show, and save a model plot
        plt.figure() 
        plt.plot(x, y_model, color = 'r', ls='-' , label = 'model')
        plt.plot(x, y, color = 'k', marker='o', ls='', label = 'experiment')
        plt.xlabel("x")
        plt.ylabel("y")
        T_as_text = format(T,'.0f')
        plt.legend(title='T = ' + T_as_text + '°C')
        plt.tight_layout()
        f_name = 'reb_19_5_1_model_' + T_as_text + '.png'
        plt.savefig('reb_19_5_1/python/' + f_name)
        plt.show()

    # show and save the fitting results
    print('\nParameter Estimates:\n')
    print(parameter_estimates)
    parameter_estimates.to_csv(
        'reb_19_5_1/python/reb_19_5_1_diff_params.csv', index=False)
    
    # fit the Arrhenius expression to the k vs T data
    T_block = parameter_estimates['T'].to_numpy() + 273.15
    k = parameter_estimates['k'].to_numpy()
    # fit the Arrhenius expression to the data
    k0, k0_ci, E, E_ci, r_squared = Arrhenius_parameters(k,T_block,R)

    # report the results
    print(' ')
    print(f'k0: {k0:.3g} /min, 95% CI [{k0_ci[0]:.3g}, {k0_ci[1]:.3g}]')
    print(f'E: {E:.3g} kJ/mol, 95% CI [{E_ci[0]:.3g}, {E_ci[1]:.3g}]')
    print(f'R-squared: {r_squared:.3g}')
    print(' ')

    # save the results to a .csv file
    data = [['k0', f'{k0:.3g}', 'min^-1^'],
        ['k0_lower_limit', f'{k0_ci[0]:.3g}', 'min^-1^'],
        ['k0_upper_limit', f'{k0_ci[1]:.3g}', 'min^-1^'],
        ['E', f'{E:.3g}', 'kJ mol^-1^'],
        ['E_lower_limit', f'{E_ci[0]:.3g}', 'kJ mol^-1^'],
        ['E_upper_limit', f'{E_ci[1]:.3g}', 'kJ mol^-1^'],
        ['R_squared', f'{r_squared:.3g}', '']]
    result = pd.DataFrame(data, columns=['item','value','units'])
    result.to_csv("reb_19_5_1/python/reb_19_5_1_Arrhenius_diff.csv", 
                index=False)

    # create, show, and save an Arrhenius plot
    y_pred = k0*np.exp(-E/R/T_block)
    plt.figure()
    plt.semilogy(1/T_block,k,color='k',marker='o', ls='none')
    plt.semilogy(1/T_block,y_pred,color='r')
    plt.xlabel('T$^{-1}$ (K$^{-1}$)')
    plt.ylabel('k (min$^{-1}$)')
    plt.xticks(rotation=25)
    plt.tight_layout()
    plt.savefig('reb_19_5_1/python/reb_19_5_1_Arrhenius_diff.png')
    plt.show()

if __name__=="__main__":
    perform_the_calculations()