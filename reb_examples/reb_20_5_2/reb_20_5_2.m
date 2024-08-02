function reb_20_5_2()
%reb_20_5_2 Reaction Engineering Basics Example 20.5.2

    % given and known constants
    P = 3.0; % atm
    T = 450 + 273.15; % K
    V = 1.0; % L (basis)
    R = 0.08206; % L*atm/mol/K

    % globally available variables
    CZ_1 = nan;
    tau = nan;
    yA_0 = nan;
    i_expt_current = -1;
    k_current = nan;
    alphaA_current = nan;
    alphaB_current = nan;

    % residuals function
    function resids = residuals(guess)
        % get dependent variables that are needed
        nA_1 = guess(1);
        nB_1 = guess(2);
        nZ_1 = guess(3);

        % calculate the other unknown quantities
        Vdot_0 = V/tau(i_expt_current);
        nA_0 = yA_0(i_expt_current)*P*Vdot_0/R/T;
        nB_0 = (1-yA_0(i_expt_current))*P*Vdot_0/R/T;
        PA = nA_1*P/(nA_1+nB_1+nZ_1);
        PB = nB_1*P/(nA_1+nB_1+nZ_1);
        r = k_current*PA^alphaA_current*PB^alphaB_current;
    
        % evaluate the residuals
        epsilon_1 = nA_0 - nA_1 - V*r;
        epsilon_2 = nB_0 - nB_1 - V*r;
        epsilon_3 = -nZ_1 + V*r;

        % return the derivatives
        resids = [epsilon_1; epsilon_2; epsilon_3];
    end

    % CSTR model function
    function [nA_1, nB_1, nZ_1] = unknowns()
        % guess the solution
        Vdot_0 = V/tau(i_expt_current);
        nA_0 = yA_0(i_expt_current)*P*Vdot_0/R/T;
        nB_0 = (1-yA_0(i_expt_current))*P*Vdot_0/R/T;
        nZ_guess = CZ_1(i_expt_current)*Vdot_0;
        nA_guess = nA_0 - nZ_guess;
        nB_guess = nB_0 - nZ_guess;
        initial_guess = [nA_guess, nB_guess, nZ_guess];
         
	    % solve the ATEs
        [soln, flag, message] = solve_ates(@residuals, initial_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    
        % extract the unknowns
        nA_1 = soln(1);
        nB_1 = soln(2);
        nZ_1 = soln(3);
    end

    % predicted responses function
    function CZ_1_model = predicted_responses(params, expt_inputs)
        % set the current value of the parameters
        k_current = 10^params(1);
        alphaA_current = params(2);
        alphaB_current = params(3);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CZ_1_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % get the inputs and make them globally available
            i_expt_current = i;

            % solve the BSTR design equations
            [nA_1, nB_1, nZ_1] = unknowns();

            % calculate the response
            CZ_1_model(i) = nZ_1*P/(nA_1 + nB_1 + nZ_1)/R/T;
        end
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../reb_20_5_2_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        tau = data(:,1);
        yA_0 = data(:,2);
        CZ_1 = data(:,3)/1000.0;
        adj_inputs = data(:,1:2);

        % guess the parameters
        par_guess = [-5.0; 1.3; 0.5];

        % estimate the parameters
        useRelErr = false;
        [beta, beta_ci, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, CZ_1, @predicted_responses, useRelErr);

        % extract the results
        k = 10^beta(1);
        k_ll = 10^beta_ci(1,1);
        k_ul = 10^beta_ci(1,2);
        alphaA = beta(2);
        alphaA_ll = beta_ci(2,1);
        alphaA_ul = beta_ci(2,2);
        alphaB = beta(3);
        alphaB_ll = beta_ci(3,1);
        alphaB_ul = beta_ci(3,2);
        
        % tabulate, show, and save the results
        item = ["k"; "k_lower_limit"; "k_upper_limit"; 
            "alphaA"; "alphaA_lower_limit"; "alphaA_upper_limit";
            "alphaB"; "alphaB_lower_limit"; "alphaB_upper_limit";
            "R_squared"];
        value = [k; k_ll; k_ul; alphaA; alphaA_ll; alphaA_ul; alphaB;
            alphaB_ll; alphaB_ul; r_squared];
        units = ["L mol^-1^ s^-1^"; "L mol^-1^ s^-1^";
            "L mol^-1^ s^-1^"; "L mol^-1^ s^-1^"; "L mol^-1^ s^-1^";
            "L mol^-1^ s^-1^"; "L mol^-1^ s^-1^"; "L mol^-1^ s^-1^";
            "L mol^-1^ s^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_20_5_2_results.csv')
        
        % calculate the model-predicted response and the residuals
        CZ_1_model = predicted_responses(beta, adj_inputs);
        residual = CZ_1 - CZ_1_model;

        % create, show, and save a parity plot
        figure
        hold on
        plot([min(CZ_1),max(CZ_1)], [min(CZ_1),max(CZ_1)],'r'...
            ,'LineWidth',2)
        plot(CZ_1, CZ_1_model,'ok','MarkerSize',10,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Outlet Z Concentration (M)','FontSize', 14)
        ylabel('Predicted Outlet Z Concentration (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_2_parity.png')

        % create show, and save residuals plots
        figure
        plot(tau, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Space Time (s)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_2_tau_residual.png')

        figure
        plot(yA_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Mole Fraction of A','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_2_yA0_residual.png')
    end

    % perform the calculations
    perform_the_calculations()
end