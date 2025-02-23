function reb_19_5_4()
%reb_16_5_4 Reaction Engineering Basics Example 19.5.4

    % given and known constants
    V = 50.0E-3; % L

    % globally available variables
    Vmax_current = nan;
    Km_current = nan;

    % derivatives function
    function derivs = derivatives(~,dep)
        % get dependent variables that are needed
        nS = dep(1);

        % calculate the rate
        CS = nS/V;
        r = Vmax_current*CS/(Km_current + CS);

        % evaluate the derivatives
        dnSdt = -r*V;
        dnPdt = r*V;

        % return the derivatives
        derivs = [dnSdt; dnPdt];
    end

    % BSTR model function
    function [t, nS, nP] = profiles(CS0, tf)
        % initial values and stopping criterion
        ind_0 = 0.0;
        nS0 = CS0*V;
        dep_0 = [nS0; 0.0];

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
        nS = dep(:,1);
        nP = dep(:,2);
    end

    % predicted responses function
    function CPf_model = predicted_responses(params, expt_inputs)
        % set the current value of k
        Vmax_current = 10^params(1);
        Km_current = 10^params(2);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CPf_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % get the inputs
            CS0 = expt_inputs(i,1);
            tf = expt_inputs(i,2);

            % solve the BSTR design equations
            [~, ~, nP] = profiles(CS0, tf);

            % calculate the response
            CPf_model(i) = nP(end)/V;
        end
    end

    % quantities of interest function
    function [Vmax, Vmax_CI, Km, Km_CI, r_squared, CPf_model...
            , epsilon_expt] = quantities_of_interest(adj_inputs, CPf)
        % guess the base 10 log of the parametersparameters
        par_guess = [0.0; 0.0];

        % estimate the parameters
        useRelErr = false;
        [beta, betaCI, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, CPf, @predicted_responses, useRelErr);

        % extract the results
        Vmax = 10^beta(1);
        Vmax_CI = 10.^betaCI(1,:);
        Km = 10^beta(2);
        Km_CI = 10.^betaCI(2,:);
        
        % calculate the model-predicted response and the residuals
        CPf_model = predicted_responses(beta, adj_inputs);
        epsilon_expt = CPf - CPf_model;
    end

    % master function
    function perform_the_calculations()
        % read the experimental data
        data_file = '../reb_19_5_4_data.csv';
        data_table = readtable(data_file, 'VariableNamingRule'...
            , 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        adj_inputs = data(:,1:2);
        CPf = data(:,3);

        [Vmax, Vmax_CI, Km, Km_CI, r_squared, CPf_model...
            , epsilon_expt] = quantities_of_interest(adj_inputs, CPf);
        
        % tabulate, show, and save the results
        item = ["Vmax"; "Vmax_lower_limit"; "Vmax_upper_limit";
            "Km"; "Km_lower_limit"; "Km_upper_limit"; "R_squared"];
        value = [Vmax; Vmax_CI(1); Vmax_CI(2); Km; Km_CI(1); Km_CI(2)...
            ; r_squared];
        units = ["mmol L^-1^ min^-1^"; "mmol L^-1^ min^-1^"; 
            "mmol L^-1^ min^-1^"; "mmol L^-1^"; "mmol L^-1^"; 
            "mmol L^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_19_5_4_results.csv')

        % create, show, and save a parity plot
        figure
        hold on
        plot([min(CPf),max(CPf)], [min(CPf),max(CPf)],'r'...
            ,'LineWidth',2)
        plot(CPf, CPf_model,'ok','MarkerSize',10,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental P Concentration (mmol/L)','FontSize', 14)
        ylabel('Predicted P Concentration (mmol/L)','FontSize', 14)
        legend({'Parity Line','Data'},'Location','northwest'...
            ,'FontSize',14)
        saveas(gcf,'reb_19_5_4_parity.png')

        % create show, and save residuals plots
        figure
        CS0 = adj_inputs(:,1);
        plot(CS0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Initial S Concentration (mmol/L','FontSize', 14)
        ylabel('Residual (mmol/L)','FontSize', 14)
        saveas(gcf,'reb_19_5_4_CS0_residuals.png')

        figure
        tf = adj_inputs(:,2);
        plot(tf, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Final time (min)','FontSize', 14)
        ylabel('Residual (mmol/L)','FontSize', 14)
        saveas(gcf,'reb_19_5_4_tf_residuals.png')

    end

    % perform the calculations
    perform_the_calculations()
end