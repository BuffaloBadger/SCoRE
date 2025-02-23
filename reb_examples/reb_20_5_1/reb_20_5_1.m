function reb_20_5_1()
%reb_20_5_1 Reaction Engineering Basics Example 20.5.1

    % given and known constants
    V = 0.1; % L
    R = 1.987E-3; % kcal/mol/K

    % globally available variables
    % full set of measured responses
    gMeasResp = nan;
    % adjusted inputs for the current experiment
    gT = nan;
    gVdot = nan;
    gCA_0 = nan;
    gCB_0 = nan;
    gCY_0 = nan;
    gCZ_0 = nan;
    % current rate expression parameters
    gk0 = nan;
    gE = nan;

    % residuals function
    function epsilon = residuals(guess)
        % extract the guesses
        nA_1 = guess(1);
        nB_1 = guess(2);
        nY_1 = guess(3);
        nZ_1 = guess(4);

        % calculate the other unknown quantities
        nA_0 = gCA_0*gVdot;
        nB_0 = gCB_0*gVdot;
        nY_0 = gCY_0*gVdot;
        nZ_0 = gCZ_0*gVdot;
        CA_1 = nA_1/gVdot;
        CB_1 = nB_1/gVdot;
        k = gk0*exp(-gE/R/gT);
        r = k*CA_1*CB_1;
    
        % evaluate the residuals
        epsilon_1 = nA_0 - nA_1 - V*r;
        epsilon_2 = nB_0 - nB_1 - V*r;
        epsilon_3 = nY_0 - nY_1 + V*r;
        epsilon_4 = nZ_0 - nZ_1 + V*r;

        % return the derivatives
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end

    % CSTR model function
    function [nA_1, nB_1, nY_1, nZ_1] = unknowns(T, Vdot, CA_0, CB_0...
            , CY_0, CZ_0, CY_1, k0, E)
        % make the adjusted inputs and parameters globally available
        gT = T;
        gVdot = Vdot;
        gCA_0 = CA_0;
        gCB_0 = CB_0;
        gCY_0 = CY_0;
        gCZ_0 = CZ_0;
        gk0 = k0;
        gE = E;

        % guess the solution
        nA_0 = CA_0*Vdot;
        nB_0 = CB_0*Vdot;
        nY_0 = CY_0*Vdot;
        nZ_0 = CZ_0*Vdot;
        nY_guess = CY_1*Vdot;
        extent = nY_guess - nY_0;
        nA_guess = nA_0 - extent;
        nB_guess = nB_0 - extent;
        nZ_guess = nZ_0 + extent;
        initial_guess = [nA_guess; nB_guess; nY_guess; nZ_guess];
         
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
        nY_1 = soln(3);
        nZ_1 = soln(4);
    end

    % predicted responses function
    function CY_1_model = predicted_responses(params, expt_inputs)
        % extract the rate expression parameters
        k0 = 10^params(1);
        E = params(2);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CY_1_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % extract the adjusted inputs for this experiment
            T = expt_inputs(i,1);
            Vdot = expt_inputs(i,2);
            CA_0 = expt_inputs(i,3);
            CB_0 = expt_inputs(i,4);
            CY_0 = expt_inputs(i,5);
            CZ_0 = expt_inputs(i,6);

            % extract the measured response for this experiment
            CY_1 = gMeasResp(i);

            % solve the BSTR design equations
            [~, ~, nY_1,~] = unknowns(T, Vdot, CA_0, CB_0...
            , CY_0, CZ_0, CY_1, k0, E);

            % calculate the response
            CY_1_model(i) = nY_1/Vdot;
        end
    end

    % quantities of interest function
    function [k0, k0_CI, E, E_CI, r_squared, CY_1_model...
            , epsilon_expt] = quantities_of_interest(adj_inputs)
        % guess the parameters
        par_guess = [4.0; 10.0];

        % estimate the parameters
        useRelErr = false;
        [beta, betaCI, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, gMeasResp, @predicted_responses...
                , useRelErr);

        % extract the results
        k0 = 10.^beta(1);
        k0_CI = 10.^betaCI(1,:);
        E = beta(2);
        E_CI = betaCI(2,:);

        % calculate the model-predicted response and the residuals
        CY_1_model = predicted_responses(beta, adj_inputs);
        epsilon_expt = gMeasResp - CY_1_model;
    end

    % master function
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../reb_20_5_1_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and measured responses
        adj_inputs = data(:,1:6);
        adj_inputs(:,2) = adj_inputs(:,2)*1.0E-3;
        gMeasResp = data(:,7);

        % calculate the quantities of interest
        [k0, k0_CI, E, E_CI, r_squared, CY_1_model, epsilon_expt] ...
            = quantities_of_interest(adj_inputs);
        
        % tabulate, show, and save the results
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit"; "E"...
            ; "E_lower_limit"; "E_upper_limit"; "R_squared"];
        value = [k0; k0_CI(1); k0_CI(2); E; E_CI(1); E_CI(2)...
            ; r_squared];
        units = ["L mol^-1^ min^-1^"; "L mol^-1^ min^-1^";
            "L mol^-1^ min^-1^"; "kcal mol^-1^"; "kcal mol^-1^"...
            ; "kcal mol^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_20_5_1_results.csv')

        % create, show, and save a parity plot
        figure
        hold on
        plot(gMeasResp, CY_1_model,'ok','MarkerSize',10,'LineWidth',2)
        plot([min(gMeasResp),max(gMeasResp)], [min(gMeasResp)...
            ,max(gMeasResp)],'r'...
            ,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Outlet Y Concentration (M)','FontSize', 14)
        ylabel('Predicted Outlet Y Concentration (M)','FontSize', 14)
        legend({'Data','Parity Line'},'Location','northwest'...
            ,'FontSize',14)
        saveas(gcf,'reb_20_5_1_parity.png')

        % create show, and save residuals plots
        figure
        plot(adj_inputs(:,1), epsilon_expt,'ok','MarkerSize',10 ...
            ,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Temperature (K)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_T_residual.png')

        figure
        plot(adj_inputs(:,2)*1.0E3, epsilon_expt,'ok','MarkerSize',10 ...
            ,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Volumetric Flow Rate (cm^3 min^-^1)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_Vdot_residual.png')

        figure
        plot(adj_inputs(:,3), epsilon_expt,'ok','MarkerSize',10 ...
            ,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of A (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CA0_residual.png')

        figure
        plot(adj_inputs(:,4), epsilon_expt,'ok','MarkerSize',10 ...
            ,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of B (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CB0_residual.png')

        figure
        plot(adj_inputs(:,5), epsilon_expt,'ok','MarkerSize',10 ...
            ,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of Y (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CY0_residual.png')

        figure
        plot(adj_inputs(:,6), epsilon_expt,'ok','MarkerSize',10 ...
            ,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of Z (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CZ0_residual.png')
    end

    % perform the calculations
    perform_the_calculations()
end