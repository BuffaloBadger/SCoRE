function reb_19_5_1()
%reb_16_5_1 Reaction Engineering Basics Example 19.5.1

    % given and known constants
    V = 1.0; % L
    R = 8.314E-3; % kJ /mol /K

    % globally available variables
    k0_current = nan;
    E_current = nan;
    T_current = nan;

    % derivatives function
    function derivs = derivatives(~,dep)
        % get nA
        nA = dep(1);

        % calculate the rate coefficient
        k = k0_current*exp(-E_current/R/T_current);

        % calculate the concentration
        CA = nA/V;

        % calculate the rate
        r = k*CA;

        % evaluate the derivatives
        dnAdt = -r*V;
        dnZdt = r*V;

        % return the derivatives
        derivs = [dnAdt; dnZdt];
    end

    % BSTR model function
    function [t, nA, nZ] = profiles(k0,E,T,CA0,tf)
        % make rate expression parameters and adjusted inputs available
        k0_current = k0;
        E_current = E;
        T_current = T;

        % initial values and stopping criterion
        ind_0 = 0.0;
        nA0 = CA0*V;
        dep_0 = [nA0; 0.0];

        % stopping criterion
        f_var = 0;
        f_val = tf;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract the dependent variable profiles
        nA = dep(:,1);
        nZ = dep(:,2);
    end

    % predicted responses function
    function CAf_model = predicted_responses(guesses, expt_inputs)
        % set the current value of k
        k0 = 10^guesses(1);
        E = guesses(2);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CAf_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % solve the BSTR design equations
            T = expt_inputs(i,1);
            CA0 = expt_inputs(i,2);
            tf = expt_inputs(i,3);
            [~, nA, ~] = profiles(k0,E,T,CA0,tf);

            % calculate the response
            CAf_model(i) = nA(end)/V;
        end
    end

    % quantities of interest function
    function [k0, k0_CI, E, E_CI, r_squared, CAf_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs, CAf)
        % guess the log 10 of k0 and E
        guess = [2.0; 20.0];

        % estimate the parameters
        useRelErr = false;
        [beta, betaCI, r_squared] = fit_to_SR_data(guess...
            , adj_inputs, CAf, @predicted_responses, useRelErr);

        % extract the results
        k0 = 10^beta(1);
        k0_CI = 10.^betaCI(1,:);
        E = beta(2);
        E_CI = betaCI(2,:);

        % calculate the predicted responses and experimental residuals
        CAf_model = predicted_responses(beta, adj_inputs);
        epsilon_expt = CAf - CAf_model;
        
    end

    % master function
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../reb_19_5_1_data.csv'...
            , 'VariableNamingRule', 'preserve');
        adj_inputs = table2array(data_table(:,2:4));
        adj_inputs(:,1) = adj_inputs(:,1) + 273.15;
        CAf = table2array(data_table(:,5));

        % calculate the quantities of interest
        [k0, k0_CI, E, E_CI, r_squared, CAf_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs, CAf);

        % create, show, and save a parity plot
        figure
        hold on
        plot([min(CAf),max(CAf)], [min(CAf),max(CAf)],'r'...
            ,'LineWidth',2)
        plot(CAf, CAf_model,'ok','MarkerSize',10,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental C_A (M)','FontSize', 14)
        ylabel('Predicted C_A (M)','FontSize', 14)
        saveas(gcf,'reb_19_5_1_parity.png')

        % create show, and save residuals plots
        figure
        tf = adj_inputs(:,3);
        plot(tf, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Elapsed Time (min)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_19_5_1_residuals_tf.png')

        figure
        CA0 = adj_inputs(:,2);
        plot(CA0,epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('CA0 (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_19_5_1_residuals_CA0.png')

        figure
        T = adj_inputs(:,1) - 273.15;
        plot(T,epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('T (°C)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_19_5_1_residuals_T.png')

        % tabulate, show, and save the results
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit"; "E"; 
            "E_lower_limit"; "E_upper_limit"; "R_squared"];
        value = [k0; k0_CI(1); k0_CI(2); E; E_CI(1); E_CI(2)...
            ; r_squared];
        units = ["min^-1^"; "min^-1^"; "min^-1^"; "kJ mol^-1^";
            "kJ mol^-1^"; "kJ mol^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_19_5_1_results.csv')
    end

    % perform the calculations
    perform_the_calculations()
end