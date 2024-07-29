function reb_19_5_3()
%reb_16_5_3 Reaction Engineering Basics Example 19.5.3

    % given and known 
    V = 100.0; %cm^3
    P0 = 6.0; % atm
    T = 275 + 273.15; % K
    R = 82.06; % cm^3 atm/mol/K

    % make current k available to all functions
    k_current = nan;

    % derivatives function
    function derivs = derivatives(~,dep)
        % get dependent variables that are needed
        nA = dep(1);
        nB = dep(2);

        % calculate the partial pressures
        PA = nA*R*T/V;
        PB = nB*R*T/V;

        % calculate the rate
        if PA < 0
            r = 0.0;
        elseif PB<0
            r = 0.0;
        else
            r =k_current*PA*sqrt(PB);
        end

        % evaluate the derivatives
        dnAdt = -r*V;
        dnBdt = -r*V;
        dnZdt = r*V;

        % return the derivatives
        derivs = [dnAdt; dnBdt; dnZdt];
    end

    % BSTR model function
    function [t, nA, nB, nZ] = profiles(PA0, tf)
        % initial values and stopping criterion
        ind_0 = 0.0;
        nA0 = PA0*V/R/T;
        PB0 = P0 - PA0;
        nB0 = PB0*V/R/T;
        dep_0 = [nA0; nB0; 0.0];

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
        nZ = dep(:,3);
    end

    % predicted responses function
    function Pf_model = predicted_responses(params, expt_inputs)
        % set the current value of k
        k_current = 10^params(1);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        Pf_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % get the inputs
            PA0 = expt_inputs(i,1);
            tf = expt_inputs(i,2);

            % solve the BSTR design equations
            [~, nA, nB, nZ] = profiles(PA0, tf);

            % calculate the response
            Pf_model(i) = (nA(end) + nB(end) + nZ(end))*R*T/V;
        end
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_file = '../reb_19_5_3_data.csv';
        data_table = readtable(data_file, 'VariableNamingRule'...
            , 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        adj_inputs = data(:,1:2);
        Pf = data(:,3);

        % guess the parameters
        par_guess = -6.0;

        % estimate the parameters
        useRelErr = false;
        [beta, betaCI, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, Pf, @predicted_responses, useRelErr);

        % extract the results
        k = beta(1);
        k_ll = betaCI(1,1);
        k_ul = betaCI(1,2);
        
        % tabulate, show, and save the results
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit"; "R_squared"];
        value = [k; k_ll; k_ul; r_squared];
        units = ["mol cm^-3^ min^-1^ atm^-1.5^"; 
            "mol cm^-3^ min^-1^ atm^-1.5^";
            "mol cm^-3^ min^-1^ atm^-1.5^^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_19_5_3_results.csv')
        
        % calculate the model-predicted response and the residuals
        Pf_model = predicted_responses(beta, adj_inputs);
        residual = Pf - Pf_model;

        % create, show, and save a parity plot
        figure
        hold on
        plot([min(Pf),max(Pf)], [min(Pf),max(Pf)],'r'...
            ,'LineWidth',2)
        plot(Pf, Pf_model,'ok','MarkerSize',10,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Final Pressure (atm)','FontSize', 14)
        ylabel('Predicted Final Pressure (atm)','FontSize', 14)
        saveas(gcf,'reb_19_5_3_parity.png')

        % create show, and save residuals plots
        figure
        PA0 = adj_inputs(:,1);
        plot(PA0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Initial Pressure of A (atm)','FontSize', 14)
        ylabel('Residual (atm)','FontSize', 14)
        saveas(gcf,'reb_19_5_3_PA0_residuals.png')

        figure
        tf = adj_inputs(:,2);
        plot(tf, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('final time (min)','FontSize', 14)
        ylabel('Residual (atm)','FontSize', 14)
        saveas(gcf,'reb_19_5_3_tf_residuals.png')

    end

    % perform the calculations
    perform_the_calculations()
end