function reb_20_5_2()
%reb_20_5_2 Reaction Engineering Basics Example 20.5.2

    % given and known constants
    P = 3.0; % atm
    V = 1.0; % L (basis)
    R = 0.08206; % L*atm/mol/K
    Ren = 8.314E-3; % kJ/mol/K

    % globally available variables
    % full set of measured responses
    gMeasResp = nan;
    % adjusted inputs for the current experiment
    gT = nan;
    gtau = nan;
    gyA_0 = nan;
    % current rate expression parameters
    gk0 = nan;
    gE = nan;
    galphaA = nan;
    galphaB = nan;

    % residuals function
    function resids = residuals(guess)
        % get dependent variables that are needed
        nA_1 = guess(1);
        nB_1 = guess(2);
        nZ_1 = guess(3);

        % calculate the other unknown quantities
        Vdot_0 = V/gtau;
        nA_0 = gyA_0*P*Vdot_0/R/gT;
        nB_0 = (1-gyA_0)*P*Vdot_0/R/gT;
        PA = nA_1*P/(nA_1+nB_1+nZ_1);
        PB = nB_1*P/(nA_1+nB_1+nZ_1);
        k = gk0*exp(-gE/Ren/gT);
        r = k*PA^galphaA*PB^galphaB;
    
        % evaluate the residuals
        epsilon_1 = nA_0 - nA_1 - V*r;
        epsilon_2 = nB_0 - nB_1 - V*r;
        epsilon_3 = -nZ_1 + V*r;

        % return the derivatives
        resids = [epsilon_1; epsilon_2; epsilon_3];
    end

    % CSTR model function
    function [nA_1, nB_1, nZ_1] = unknowns(T, tau, yA_0, CZ_1, k0...
            , E, alphaA, alphaB)
        % make the adjusted inputs globally available
        gT = T;
        gtau = tau;
        gyA_0 = yA_0;
        gk0 = k0;
        gE = E;
        galphaA = alphaA;
        galphaB = alphaB;

        % guess the solution
        Vdot_0 = V/tau;
        nA_0 = yA_0*P*Vdot_0/R/T;
        nB_0 = (1-yA_0)*P*Vdot_0/R/T;
        nZ_guess = CZ_1*Vdot_0;
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
        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CZ_1_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % solve the BSTR design equations
            [nA_1, nB_1, nZ_1] = unknowns(expt_inputs(i,1)...
                , expt_inputs(i,2), expt_inputs(i,3), gMeasResp(i)...
                , 10^params(1), params(2), params(3), params(4));

            % calculate the response
            CZ_1_model(i) = nZ_1*P/(nA_1 + nB_1 + nZ_1)/R/expt_inputs(i,1);
        end
    end

    % quantities of interest function
    function [k0, k0_CI, E, E_CI, alphaA, alphaA_CI, alphaB...
            , alphaB_CI, r_squared, CZ_1_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs)
        % guess the parameters
        %par_guess = [0.0; 75.0; 1.0; 1.0];
        par_guess = [log10(7208.6); 113.65; 1.3911; 0.59762];

        % estimate the parameters
        useRelErr = false;
        [beta, beta_ci, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, gMeasResp, @predicted_responses...
                , useRelErr);

        % extract the results
        k0 = 10.^beta(1);
        k0_CI = 10.^beta_ci(1,:);
        E = beta(2);
        E_CI = beta_ci(2,:);
        alphaA = beta(3);
        alphaA_CI = beta_ci(3,:);
        alphaB = beta(4);
        alphaB_CI = beta_ci(4,:);

        % calculate the model-predicted response and the residuals
        CZ_1_model = predicted_responses(beta, adj_inputs);
        epsilon_expt = gMeasResp - CZ_1_model;
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../reb_20_5_2_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        gMeasResp = data(:,4)/1000.0;
        adj_inputs = data(:,1:3);
        adj_inputs(:,1) = adj_inputs(:,1) + 273.15;

        % calculate the quantities of interest
        [k0, k0_CI, E, E_CI, alphaA, alphaA_CI, alphaB...
            , alphaB_CI, r_squared, CZ_1_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs);
        
        % tabulate, show, and save the results
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit";
            "E"; "E_lower_limit"; "E_upper_limit";
            "alphaA"; "alphaA_lower_limit"; "alphaA_upper_limit";
            "alphaB"; "alphaB_lower_limit"; "alphaB_upper_limit";
            "R_squared"];
        value = [k0; k0_CI(1); k0_CI(2); E; E_CI(1); E_CI(2); alphaA...
            ; alphaA_CI(1); alphaA_CI(2); alphaB; alphaB_CI(1)...
            ; alphaB_CI(2); r_squared];
        units = ["L mol^-1^ s^-1^"; "L mol^-1^ s^-1^"...
            ; "L mol^-1^ s^-1^"; "kJ/mol"; "kJ/mol"; "kJ/mol"...
            ; " "; " "; " "; " "; " "; " "; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_20_5_2_results.csv')

        % create, show, and save a parity plot
        figure
        hold on
        plot(gMeasResp, CZ_1_model,'ok','MarkerSize',10,'LineWidth',2)
        plot([min(gMeasResp),max(gMeasResp)], [min(gMeasResp) ...
            ,max(gMeasResp)],'r','LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Outlet Z Concentration (M)','FontSize', 14)
        ylabel('Predicted Outlet Z Concentration (M)','FontSize', 14)
        legend({'Data','Parity Line'},'Location', 'northwest'...
            ,'FontSize',14)
        saveas(gcf,'reb_20_5_2_parity.png')

        % create show, and save residuals plots
        figure
        T = adj_inputs(:,1) - 273.15;
        plot(T, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('T (°C)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_2_T_residual.png')

        figure
        tau = adj_inputs(:,2);
        plot(tau, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Space Time (s)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_2_tau_residual.png')

        figure
        yA_0 = adj_inputs(:,3);
        plot(yA_0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Mole Fraction of A','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_2_yA0_residual.png')
    end

    % perform the calculations
    perform_the_calculations()
end