function reb_20_5_1()
%reb_20_5_1 Reaction Engineering Basics Example 20.5.1

    % given and known constants
    V = 0.1; % L

    % globally available variables
    Vdot_current = nan;
    CA_0_current = nan;
    CB_0_current = nan;
    CY_0_current = nan;
    CZ_0_current = nan;
    k_current = nan;

    % residuals function
    function resids = residuals(guess)
        % get dependent variables that are needed
        nA_1 = guess(1);
        nB_1 = guess(2);
        nY_1 = guess(3);
        nZ_1 = guess(4);

        % calculate the other unknown quantities
        nA_0 = CA_0_current*Vdot_current;
        nB_0 = CB_0_current*Vdot_current;
        nY_0 = CY_0_current*Vdot_current;
        nZ_0 = CZ_0_current*Vdot_current;
        CA_1 = nA_1/Vdot_current;
        CB_1 = nB_1/Vdot_current;
        r = k_current*CA_1*CB_1;
    
        % evaluate the residuals
        epsilon_1 = nA_0 - nA_1 - V*r;
        epsilon_2 = nB_0 - nB_1 - V*r;
        epsilon_3 = nY_0 - nY_1 + V*r;
        epsilon_4 = nZ_0 - nZ_1 + V*r;

        % return the derivatives
        resids = [epsilon_1; epsilon_2; epsilon_3; epsilon_4];
    end

    % CSTR model function
    function [nA_1, nB_1, nY_1, nZ_1] = unknowns()
        % guess the solution
        initial_guess = [0.25; 0.25; 0.25; 0.25];
         
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
        % set the current value of k
        k_current = 10^params(1);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CY_1_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % get the inputs and make them globally available
            Vdot_current = expt_inputs(i,1);
            CA_0_current = expt_inputs(i,2);
            CB_0_current = expt_inputs(i,3);
            CY_0_current = expt_inputs(i,4);
            CZ_0_current = expt_inputs(i,5);

            % solve the BSTR design equations
            [~, ~, nY_1,~] = unknowns();

            % calculate the response
            CY_1_model(i) = nY_1/Vdot_current;
        end
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../reb_20_5_1_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        adj_inputs = data(:,1:5);
        adj_inputs(:,1) = adj_inputs(:,1)*1.0E-3;
        CY_1 = data(:,6);

        % guess the parameters
        par_guess = 0.0;

        % estimate the parameters
        useRelErr = false;
        [beta, betaCI, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, CY_1, @predicted_responses, useRelErr);

        % extract the results
        k = 10^beta(1);
        k_ll = 10^betaCI(1,1);
        k_ul = 10^betaCI(1,2);
        
        % tabulate, show, and save the results
        item = ["k"; "k_lower_limit"; "k_upper_limit"; "R_squared"];
        value = [k; k_ll; k_ul; r_squared];
        units = ["L mol^-1^ min^-1^"; "L mol^-1^ min^-1^";
            "L mol^-1^ min^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_20_5_1_results.csv')
        
        % calculate the model-predicted response and the residuals
        CY_1_model = predicted_responses(beta, adj_inputs);
        residual = CY_1 - CY_1_model;

        % create, show, and save a parity plot
        figure
        hold on
        plot([min(CY_1),max(CY_1)], [min(CY_1),max(CY_1)],'r'...
            ,'LineWidth',2)
        plot(CY_1, CY_1_model,'ok','MarkerSize',10,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Outlet Y Concentration (M)','FontSize', 14)
        ylabel('Predicted Outlet Y Concentration (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_parity.png')

        % create show, and save residuals plots
        figure
        Vdot = adj_inputs(:,1);
        plot(Vdot, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Volumetric Flow Rate (cm^3 min^-^1)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_Vdot_residual.png')

        % create show, and save residuals plots
        figure
        CA_0 = adj_inputs(:,2);
        plot(CA_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of A (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CA0_residual.png')

        % create show, and save residuals plots
        figure
        CB_0 = adj_inputs(:,3);
        plot(CB_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of B (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CB0_residual.png')

        % create show, and save residuals plots
        figure
        CY_0 = adj_inputs(:,4);
        plot(CY_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of Y (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CY0_residual.png')

        % create show, and save residuals plots
        figure
        CZ_0 = adj_inputs(:,5);
        plot(CZ_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet Concentration of Z (M)','FontSize', 14)
        ylabel('Residual (M)','FontSize', 14)
        saveas(gcf,'reb_20_5_1_CZ0_residual.png')

    end

    % perform the calculations
    perform_the_calculations()
end