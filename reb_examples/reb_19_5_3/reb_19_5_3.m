function reb_19_5_3()
%reb_16_5_3 Reaction Engineering Basics Example 19.5.3

    % given and known constants
    V = 100.0; %cm^3
    P0 = 6.0; % 
    Rpv = 82.06; % cm^3 atm/mol/K
    Ren = 1.987E-3; % kcal/mol

    % globally available variables
    k0_current = nan;
    E_current = nan;
    T_current = nan;

    % derivatives function
    function derivs = derivatives(~,dep)
        % get dependent variables that are needed
        nA = dep(1);
        nB = dep(2);

        % calculate the rate coefficient
        k = k0_current*exp(-E_current/Ren/T_current);

        % calculate the partial pressures
        PA = nA*Rpv*T_current/V;
        PB = nB*Rpv*T_current/V;

        % calculate the rate
        if PA < 0
            r = 0.0;
        elseif PB<0
            r = 0.0;
        else
            r =k*PA*sqrt(PB);
        end

        % evaluate the derivatives
        dnAdt = -r*V;
        dnBdt = -r*V;
        dnZdt = r*V;

        % return the derivatives
        derivs = [dnAdt; dnBdt; dnZdt];
    end

    % BSTR model function
    function [t, nA, nB, nZ] = profiles(T, PA0, tf)
        % initial values and stopping criterion
        ind_0 = 0.0;
        nA0 = PA0*V/Rpv/T;
        PB0 = P0 - PA0;
        nB0 = PB0*V/Rpv/T;
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
        % make the rate expression parameters globally available
        k0_current = params(1);
        E_current = params(2);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        Pf_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % get the inputs
            T = expt_inputs(i,1);
            PA0 = expt_inputs(i,2);
            tf = expt_inputs(i,3);

            % make the temperature globally available
            T_current = T;

            % solve the BSTR design equations
            [~, nA, nB, nZ] = profiles(T, PA0, tf);

            % calculate the response
            Pf_model(i) = (nA(end) + nB(end) + nZ(end))*Rpv*T/V;
        end
    end

    % quantities of interest function
    function [k0, k0_CI, E, E_CI, r_squared, Pf_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs, Pf)
        % guess the parameters
        par_guess = [1.0; 15.0];

        % estimate the parameters
        useRelErr = false;
        [beta, betaCI, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, Pf, @predicted_responses, useRelErr);

        % extract the results
        k0 = beta(1);
        k0_CI = betaCI(1,:);
        E = beta(2);
        E_CI = betaCI(2,:);

        % calculate the model-predicted response and the residuals
        Pf_model = predicted_responses(beta, adj_inputs);
        epsilon_expt = Pf - Pf_model;
    end

    % master function
    function perform_the_calculations()
        % read the experimental data
        data_file = '../reb_19_5_3_data.csv';
        data_table = readtable(data_file, 'VariableNamingRule'...
            , 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        adj_inputs = data(:,1:3);
        adj_inputs(:,1) = adj_inputs(:,1) + 273.15;
        Pf = data(:,4);

        % calculate the quantities of interest
        [k0, k0_CI, E, E_CI, r_squared, Pf_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs, Pf);
        
        % tabulate, show, and save the results
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit";...
            "E"; "E_lower_limit"; "E_upper_limit"; "R_squared"];
        value = [k0; k0_CI(1); k0_CI(2); E; E_CI(1); E_CI(2)...
            ; r_squared];
        units = ["mol cm^-3^ min^-1^ atm^-1.5^"; 
            "mol cm^-3^ min^-1^ atm^-1.5^";
            "mol cm^-3^ min^-1^ atm^-1.5^^"; "kcal/mol"; "kcal/mol";...
            "kcal/mol"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_19_5_3_results.csv')

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
        legend({'Parity Line','Data'},'Location','northwest',...
            'FontSize',14)
        saveas(gcf,'reb_19_5_3_parity.png')

        % create show, and save residuals plots
        figure
        T = adj_inputs(:,1) - 273.15;
        plot(T, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Temperature (°C)','FontSize', 14)
        ylabel('Residual (atm)','FontSize', 14)
        saveas(gcf,'reb_19_5_3_T_residuals.png')

        figure
        PA0 = adj_inputs(:,2);
        plot(PA0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Initial Pressure of A (atm)','FontSize', 14)
        ylabel('Residual (atm)','FontSize', 14)
        saveas(gcf,'reb_19_5_3_PA0_residuals.png')

        figure
        tf = adj_inputs(:,3);
        plot(tf, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('final time (min)','FontSize', 14)
        ylabel('Residual (atm)','FontSize', 14)
        saveas(gcf,'reb_19_5_3_tf_residuals.png')

    end

    % perform the calculations
    perform_the_calculations()
end