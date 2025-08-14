function example_4_5_4_matlab_calculations
% Calculations for Example 4.5.4 of Reaction Engineering Basics
    % Given and known
    R = 1.987E-3; % kcal/mol/K

    % Read the data into a table
    tbl = readtable('example_4_5_4_data.csv','VariableNamingRule',...
        'preserve');

    % Extract the columns from the table
    T = table2array(tbl(:,1)); % °C
    k = table2array(tbl(:,2)); % L/mol/min

    % Convert to absolute temperatures
    T = T + 273.15;

    % Calculate the Arrhenius expression parameters
    [k0, k0_ci, E, E_ci, r_squared] = Arrhenius_parameters(k,T,R);
    
    % Report the results
    disp('')
    disp(['k0 = ',num2str(k0,3),' L mol^-^1 min^-^1, 95% CI ['...
        ,num2str(k0_ci(1),3),', ',num2str(k0_ci(2),3),']'])
    disp(['E = ',num2str(E,3),' kcal mol^-^1, 95% CI ['...
        ,num2str(E_ci(1),3),', ',num2str(E_ci(2),3),']'])
    disp(['R = ',num2str(r_squared,3)])

    % Save the results to a .csv file
    item = ["k0"; "k0_lower_limit"; "k0_upper_limit"; "E";...
        "E_lower_limit"; "E_upper_limit"; "R-squared"];
    value = [k0; k0_ci(1); k0_ci(2); E; E_ci(1); E_ci(2); r_squared];
    units = ["s^-1^"; "s^-1^"; ...
        "s^-1^"; "kcal mol^-1^"; "kcal mol^-1^"; ...
        "kcal mol^-1^"; " "];
    results_table = table(item,value,units);
    writetable(results_table,'example_4_5_4_matlab_results.csv');

    % Generate a model plot
    y_model = log(k0*exp(-E/R./T));
    y = log(k);
    figure;
    plot(1/R./T,y_model,'r',1/R./T,y,'ok','LineWidth',2)
    set(gca, 'FontSize', 14);
    xlabel('1/(RT) (/K)','FontSize', 14)
    ylabel('ln(k)','FontSize', 14)
    legend({'Arrhenius Expression','Experimental Data'},'Location'...
        ,'northeast','FontSize',14)

    % Save the model plot to a .png file    
    saveas(gcf,"example_4_5_4_matlab_model_plot.png")
end