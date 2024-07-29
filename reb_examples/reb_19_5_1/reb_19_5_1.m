function reb_19_5_1()
%reb_16_5_1 Reaction Engineering Basics Example 19.5.1

    % given and known 
    V = 1.0; % L
    R = 8.314E-3; % kJ /mol /K

    % make k available to all functions
    k_current = nan;

    % derivatives function
    function derivs = derivatives(~,dep)
        % get nA
        nA = dep(1);

        % calculate the concentration
        CA = nA/V;

        % calculate the rate
        r = k_current*CA;

        % evaluate the derivatives
        dnAdt = -r*V;
        dnZdt = r*V;

        % return the derivatives
        derivs = [dnAdt; dnZdt];
    end

    % BSTR model function
    function [t, nA, nZ] = profiles(CA0,tf)
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
    function CAf_model = predicted_responses(log_k_guess, expt_inputs)
        % set the current value of k
        k_current = 10^log_k_guess;

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CAf_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % solve the BSTR design equations
            CA0 = expt_inputs(i,2);
            tf = expt_inputs(i,3);
            [~, nA, ~] = profiles(CA0, tf);

            % calculate the response
            CAf_model(i) = nA(end)/V;
        end
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_file = '../reb_19_5_1_data.csv';
        data_table = readtable(data_file, 'VariableNamingRule'...
            , 'preserve');
        data = table2array(data_table(:,2:end));

        % get the temperatures of the data blocks
        block_temperatures = unique(data(:,1));
        n_blocks = size(block_temperatures,1);

        % create vectors to store the fitting results
        k_fit = nan(n_blocks,1);
        kll_fit = nan(n_blocks,1);
        kul_fit = nan(n_blocks,1);
        rsq_fit = nan(n_blocks,1);

        % process the data blocks
        for iBlock = 1:n_blocks
            % create the block
            idx = table2array(data_table(:,2) == block_temperatures(iBlock));
            data_block = data(idx,:);

            % extract the adjusted inputs and response
            adj_inputs = data_block(:,1:3);
            CAf = data_block(:,4);

            % guess log 10 of k
            par_guess = -2.0;

            % estimate log 10 of k
            useRelErr = false;
            [beta, betaCI, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, CAf, @predicted_responses, useRelErr);

            % extract the results
            k_fit(iBlock) = 10^beta;
            kll_fit(iBlock) = 10^betaCI(1);
            kul_fit(iBlock) = 10^betaCI(2);
            rsq_fit(iBlock) = r_squared;

            % calculate the model-predicted response and the residuals
            CAf_model = predicted_responses(beta, adj_inputs);
            residual = CAf - CAf_model;
            T_as_text = int2str(data_block(1,1));

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
            fname = strcat('reb_19_5_1_resp_fcn_parity_',T_as_text...
                ,'.png');
            saveas(gcf,fname)

            % create show, and save residuals plots
            figure
            tf = adj_inputs(:,3);
            plot(tf, residual,'ok','MarkerSize',10,'LineWidth',2)
            yline(0.0,'r','LineWidth',2)
            set(gca, 'FontSize', 14);
            xlabel('Elapsed Time (min)','FontSize', 14)
            ylabel('Residual (M)','FontSize', 14)
            fname = strcat('reb_19_5_1_resp_fcn_residual_tf_'...
                ,T_as_text,'.png');
            saveas(gcf,fname)

            figure
            CA0 = adj_inputs(:,2);
            plot(CA0, residual,'ok','MarkerSize',10,'LineWidth',2)
            yline(0.0,'r','LineWidth',2)
            set(gca, 'FontSize', 14);
            xlabel('CA0 (M)','FontSize', 14)
            ylabel('Residual (M)','FontSize', 14)
            fname = strcat('reb_19_5_1_resp_fcn_residual_CA0_'...
                ,T_as_text,'.png');
            saveas(gcf,fname)
            
        end

        % tabulate, show, and save the results
        results_file ="reb_19_5_1_resp_fcn_params.csv";
        results_table = table(block_temperatures, k_fit, kll_fit...
            , kul_fit, rsq_fit);
        results_table.Properties.VariableNames = ["T", "k", "k_ll"...
            , "k_ul", "R_sq"];
        disp(results_table)
        writetable(results_table,results_file);

        % fit the Arrhenius expression to the k vs. T data
        T = block_temperatures + 273.15;
        [k0, k0_ci, E, E_ci, r_squared] = Arrhenius_parameters(k_fit...
            ,T,R);

        % tabulate, show, and save the results
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit"; "E"; 
            "E_lower_limit"; "E_upper_limit"; "R_squared"];
        value = [k0; k0_ci(1); k0_ci(2); E; E_ci(1); E_ci(2);
             r_squared];
        units = ["min^-1^"; "min^-1^"; "min^-1^"; "kJ mol^-1^";
            "kJ mol^-1^"; "kJ mol^-1^"; ""];
        Arrhenius_table = table(item, value, units);
        disp(Arrhenius_table)
        writetable(Arrhenius_table,'reb_19_5_1_Arrhenius_resp_fcn.csv')

        % create, show and save an Arrhenius plot
        k_pred = k0*exp(-E/R./T);
        x = 1./T;
        figure
        semilogy(x, k_fit,'ok', x, k_pred,'r', 'MarkerSize' ,10 ...
            ,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('1/T (/K)','FontSize', 14)
        ylabel('k (/min)','FontSize', 14)
        saveas(gcf,'reb_19_5_1_Arrhenius_resp_fcn.png')
    end

    % perform the calculations
    perform_the_calculations()
end