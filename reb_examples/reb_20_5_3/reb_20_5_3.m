function reb_20_5_3()
%reb_20_5_3 Reaction Engineering Basics Example 20.5.3

    % given and known constants
    V = 3.0; % gal

    % globally available variables
    Vdot = nan;
    CA_0 = nan;
    CY_0 = nan;
    CZ_0 = nan;
    CA_1 = nan;
    i_expt_current = -1;
    kf_current = nan;
    kr_current = nan;

    % residuals function
    function resids = residuals(guess)
        % get dependent variables that are needed
        nA_1 = guess(1);
        nY_1 = guess(2);
        nZ_1 = guess(3);

        % calculate the other unknown quantities
        nA_0 = CA_0(i_expt_current)*Vdot(i_expt_current);
        nY_0 = CY_0(i_expt_current)*Vdot(i_expt_current);
        nZ_0 = CZ_0(i_expt_current)*Vdot(i_expt_current);
        CA = nA_1/Vdot(i_expt_current);
        CY = nY_1/Vdot(i_expt_current);
        CZ = nZ_1/Vdot(i_expt_current);
        r = kf_current*CA^2 - kr_current*CY*CZ;
    
        % evaluate the residuals
        epsilon_1 = nA_0 - nA_1 - V*r;
        epsilon_2 = nY_0 - nY_1 + V*r;
        epsilon_3 = nZ_0 - nZ_1 + V*r;

        % return the derivatives
        resids = [epsilon_1; epsilon_2; epsilon_3];
    end

    % CSTR model function
    function [nA_1, nY_1, nZ_1] = unknowns()
        % guess the solution
        nA_0 = CA_0(i_expt_current)*Vdot(i_expt_current);
        nY_0 = CY_0(i_expt_current)*Vdot(i_expt_current);
        nZ_0 = CZ_0(i_expt_current)*Vdot(i_expt_current);
        nA_guess = CA_1(i_expt_current)*Vdot(i_expt_current);
        xi = nA_0 - nA_guess;
        nY_guess = nY_0 + xi;
        nZ_guess = nZ_0 + xi;
        initial_guess = [nA_guess, nY_guess, nZ_guess];
         
	    % solve the ATEs
        [soln, flag, message] = solve_ates(@residuals, initial_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    
        % extract the unknowns
        nA_1 = soln(1);
        nY_1 = soln(2);
        nZ_1 = soln(3);
    end

    % predicted responses function
    function CA_1_model = predicted_responses(params, expt_inputs)
        % set the current value of the parameters
        kf_current = 10^params(1);
        kr_current = 10^params(2);

        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CA_1_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % make the experiment index globally available
            i_expt_current = i;

            % solve the BSTR design equations
            [nA_1, ~, ~] = unknowns();

            % calculate the response
            CA_1_model(i) = nA_1/Vdot(i);
        end
    end

    % function that performs the calculations
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../reb_20_5_3_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        Vdot = data(:,1);
        CA_0 = data(:,2);
        CY_0 = data(:,3);
        CZ_0 = data(:,4);
        CA_1 = data(:,5);
        adj_inputs = data(:,1:4);

        % guess the parameters
        par_guess = [0.0; 0.0];

        % estimate the parameters
        useRelErr = false;
        [beta, beta_ci, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, CA_1, @predicted_responses, useRelErr);

        % extract the results
        kf = 10^beta(1);
        kf_ll = 10^beta_ci(1,1);
        kf_ul = 10^beta_ci(1,2);
        kr = 10^beta(2);
        kr_ll = 10^beta_ci(2,1);
        kr_ul = 10^beta_ci(2,2);
        
        % tabulate, show, and save the results
        item = ["kf"; "kf_lower_limit"; "kf_upper_limit"; 
            "kr"; "kr_lower_limit"; "kr_upper_limit"; "R_squared"];
        value = [kf; kf_ll; kf_ul; kr; kr_ll; kr_ul; r_squared];
        units = ["gal lbmol^-1^ min^-1^"; "gal lbmol^-1^ min^-1^";
            "gal lbmol^-1^ min^-1^"; "gal lbmol^-1^ min^-1^"; 
            "gal lbmol^-1^ min^-1^"; "gal lbmol^-1^ min^-1^"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_20_5_3_results.csv')
        
        % calculate the model-predicted response and the residuals
        CA_1_model = predicted_responses(beta, adj_inputs);
        residual = CA_1 - CA_1_model;

        % create, show, and save a parity plot
        figure
        hold on
        plot([min(CA_1),max(CA_1)], [min(CA_1),max(CA_1)],'r'...
            ,'LineWidth',2)
        plot(CA_1, CA_1_model,'ok','MarkerSize',10,'LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Outlet C_A (lbmol/gal)','FontSize', 14)
        ylabel('Predicted Outlet C_A (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_parity.png')

        % create show, and save residuals plots
        figure
        plot(Vdot, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('VFR (gal/min)','FontSize', 14)
        ylabel('Residual (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_Vdot_residual.png')

        figure
        plot(CA_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet C_A (lbmol/gal)','FontSize', 14)
        ylabel('Residual (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_CA0_residual.png')

        figure
        plot(CY_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet C_Y (lbmol/gal)','FontSize', 14)
        ylabel('Residual (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_CY0_residual.png')

        figure
        plot(CZ_0, residual,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet C_Z (lbmol/gal)','FontSize', 14)
        ylabel('Residual (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_CZ0_residual.png')
    end

    % perform the calculations
    perform_the_calculations()
end