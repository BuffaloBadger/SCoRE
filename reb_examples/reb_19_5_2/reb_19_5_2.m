function reb_19_5_2()
%reb_16_5_1 Reaction Engineering Basics Example 19.5.2

    % given and known 
    V = 500.0; % cc
    Re = 1.987E-3; % kcal /mol /K
    Rpv = 82.06; % cc atm /mol /K

    % make current k0, E, and T available to all functions
    k0_current = nan;
    E_current = nan;
    T_current = nan;

    % derivatives function
    function derivs = derivatives(~,dep)
        % get dependent variables that are needed
        nA = dep(1);
        nB = dep(2);

        % calculate the partial pressures
        PA = nA*Rpv*T_current/V;
        PB = nB*Rpv*T_current/V;

        % calculate the rate
        k = k0_current*exp(-E_current/Re/T_current);
        r = k*PA*PB;

        % evaluate the derivatives
        dnAdt = -r*V;
        dnBdt = -r*V;
        dnYdt = r*V;
        dnZdt = r*V;

        % return the derivatives
        derivs = [dnAdt; dnBdt; dnYdt; dnZdt];
    end

    % BSTR model function
    function [t, nA, nB, nY, nZ] = profiles(PA0, PB0, tf)
        % initial values and stopping criterion
        ind_0 = 0.0;
        nA0 = PA0*V/Rpv/T_current;
        nB0 = PB0*V/Rpv/T_current;
        dep_0 = [nA0; nB0; 0.0; 0.0];

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
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
    end

    % predicted responses function
    function fA_model = predicted_responses(params, expt_inputs)
        % set the current value of k
        k0_current = params(1);
        E_current = params(2);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        fA_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % get the inputs
            T_current = expt_inputs(i,1);
            PA0 = expt_inputs(i,2);
            PB0 = expt_inputs(i,3);
            tf = expt_inputs(i,4);

            % solve the BSTR design equations
            [~, nA, ~, ~, ~] = profiles(PA0, PB0, tf);

            % calculate the response
            nA0 = PA0*V/Rpv/T_current;
            nAf = nA(end);
            fA_model(i) = (nA0 - nAf)/nA0;
        end
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_file = '../reb_19_5_2_data.csv';
        data_table = readtable(data_file, 'VariableNamingRule'...
            , 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        adj_inputs = data(:,1:4);
        adj_inputs(:,1) = adj_inputs(:,1) + 273.15;
        fA = data(:,5);

        % guess the parameters
        par_guess = [0.0; 20.0];

        % estimate the parameters
        useRelErr = false;
        [beta, betaCI, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, fA, @predicted_responses, useRelErr);

        % extract the results
        k0 = beta(1);
        k0_ll = betaCI(1,1);
        k0_ul = betaCI(1,2);
        E = beta(2);
        E_ll = betaCI(2,1);
        E_ul = betaCI(2,2);
        
        % calculate the model-predicted response and the residuals
        fA_model = predicted_responses(beta, adj_inputs);
        residual = fA - fA_model;

        % create, show, and save a parity plot
        figure
        hold on
        plot([min(fA),max(fA)], [min(fA),max(fA)],'r'...
            ,'LineWidth',2)
        plot(fA, fA_model,'ok','MarkerSize',10,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Conversion','FontSize', 14)
        ylabel('Predicted Conversion','FontSize', 14)
        saveas(gcf,'reb_19_5_2_parity.png')

        % create show, and save residuals plots
        figure
        tf = adj_inputs(:,4);
        plot(tf, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Elapsed Time (min)','FontSize', 14)
        ylabel('Residual','FontSize', 14)
        saveas(gcf,'reb_19_5_2_residuals_tf.png')

        figure
        PA0 = adj_inputs(:,2);
        plot(PA0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('PA0 (atm)','FontSize', 14)
        ylabel('Residual','FontSize', 14)
        saveas(gcf,'reb_19_5_2_residuals_PA0.png')

        figure
        PB0 = adj_inputs(:,3);
        plot(PB0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('PB0 (atm)','FontSize', 14)
        ylabel('Residual','FontSize', 14)
        saveas(gcf,'reb_19_5_2_residuals_PB0.png')

        figure
        T = adj_inputs(:,1);
        plot(T, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('T (C)','FontSize', 14)
        ylabel('Residual','FontSize', 14)
        saveas(gcf,'reb_19_5_2_residuals_T.png')

        % tabulate, show, and save the results
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit"; "E"; 
            "E_lower_limit"; "E_upper_limit"; "R_squared"];
        value = [k0; k0_ll; k0_ul; E; E_ll; E_ul;
             r_squared];
        units = ["mol cm^-3^ min^-1^ atm^-2^"; 
            "mol cm^-3^ min^-1^ atm^-2^"; "mol cm^-3^ min^-1^ atm^-2^^";
            "kcal mol^-1^"; "kcal mol^-1^"; "kcal mol^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_19_5_2_results.csv')
    end

    % perform the calculations
    perform_the_calculations()
end