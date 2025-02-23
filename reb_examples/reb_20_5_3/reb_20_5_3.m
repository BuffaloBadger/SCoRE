function reb_20_5_3()
%reb_20_5_3 Reaction Engineering Basics Example 20.5.3

    % given and known constants
    V = 3.0; % gal
    Ren = 1.987; % BTU/lbmol

    % globally available variables
    % full det of measured responses
    gMeasResp = nan;
    % adjusted inputs for the current experiment
    gT = nan;
    gVdot = nan;
    gCA_0 = nan;
    gCY_0 = nan;
    gCZ_0 = nan;
    % current rate expression parameters
    gk0f = nan;
    gEf = nan;
    gk0r = nan;
    gEr = nan;

    % residuals function
    function resids = residuals(guess)
        % get dependent variables that are needed
        nA_1 = guess(1);
        nY_1 = guess(2);
        nZ_1 = guess(3);

        % calculate the other unknown quantities
        nA_0 = gCA_0*gVdot;
        nY_0 = gCY_0*gVdot;
        nZ_0 = gCZ_0*gVdot;
        CA = nA_1/gVdot;
        CY = nY_1/gVdot;
        CZ = nZ_1/gVdot;
        kf = gk0f*exp(-gEf/Ren/gT);
        kr = gk0r*exp(-gEr/Ren/gT);
        r = kf*CA^2 - kr*CY*CZ;
    
        % evaluate the residuals
        epsilon_1 = nA_0 - nA_1 - V*r;
        epsilon_2 = nY_0 - nY_1 + V*r;
        epsilon_3 = nZ_0 - nZ_1 + V*r;

        % return the derivatives
        resids = [epsilon_1; epsilon_2; epsilon_3];
    end

    % CSTR model function
    function [nA_1, nY_1, nZ_1] = unknowns(T, Vdot, CA_0, CY_0, CZ_0...
            , CA_1, k0f, Ef, k0r, Er)
        % make the adjusted inputs and parameters globally available
        gT = T;
        gVdot = Vdot;
        gCA_0 = CA_0;
        gCY_0 = CY_0;
        gCZ_0 = CZ_0;
        gk0f = k0f;
        gEf = Ef;
        gk0r = k0r;
        gEr = Er;

        % guess the solution
        nA_guess = CA_1*Vdot;
        xi = CA_0*Vdot - nA_guess;
        nY_guess = CY_0*Vdot + xi;
        nZ_guess = CZ_0*Vdot + xi;
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
        % get the number of experiments
        n_expt = size(expt_inputs,1);

        % allocate storage for the predicted responses
        CA_1_model = nan(n_expt,1);

        % loop through the experiments
        for i = 1:n_expt
            % solve the BSTR design equations
            [nA_1, ~, ~] = unknowns(expt_inputs(i,1)...
                , expt_inputs(i,2), expt_inputs(i,3)...
                , expt_inputs(i,4), expt_inputs(i,5), gMeasResp(i)...
                , 10^params(1), params(2), 10^params(3), params(4));

            % calculate the response
            CA_1_model(i) = nA_1/expt_inputs(i,2);
        end
    end

    % quantities of interest function
    function [k0f, k0f_CI, Ef, Ef_CI, k0r, k0r_CI, Er, Er_CI...
            , r_squared, CA_1_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs)
        % guess the parameters
        par_guess = [6.0; 15000.0; 6.0; 15000.0];
        %par_guess = [log10(5.34E5); 11800; log10(6.44E7); 18500];

        % estimate the parameters
        useRelErr = false;
        [beta, beta_ci, r_squared] = fit_to_SR_data(par_guess...
                , adj_inputs, gMeasResp, @predicted_responses...
                , useRelErr);

        % extract the results
        k0f = 10.^beta(1);
        k0f_CI = 10.^beta_ci(1,:);
        Ef = beta(2);
        Ef_CI = beta_ci(2,:);
        k0r = 10.^beta(3);
        k0r_CI = 10.^beta_ci(3,:);
        Er = beta(4);
        Er_CI = beta_ci(4,:);

        % calculate the model-predicted response and the residuals
        CA_1_model = predicted_responses(beta, adj_inputs);
        epsilon_expt = gMeasResp - CA_1_model;
    end

    % master function
    function perform_the_calculations()
        % read the experimental data
        data_table = readtable('../reb_20_5_3_data.csv'...
            , 'VariableNamingRule', 'preserve');
        data = table2array(data_table(:,:));

        % extract the adjusted inputs and response
        T = data(:,1) + 459.7;
        Vdot = data(:,2);
        CA_0 = data(:,3);
        CY_0 = data(:,4);
        CZ_0 = data(:,5);
        gMeasResp = data(:,6);
        adj_inputs = data(:,1:5);
        adj_inputs(:,1) = adj_inputs(:,1) + 459.7;

        % calculate the quantities of interest
        [k0f, k0f_CI, Ef, Ef_CI, k0r, k0r_CI, Er, Er_CI...
            , r_squared, CA_1_model, epsilon_expt]...
            = quantities_of_interest(adj_inputs);

        % tabulate, show, and save the results
        item = ["k0f"; "k0f_lower_limit"; "k0f_upper_limit";
            "Ef"; "Ef_lower_limit"; "Ef_upper_limit";
            "k0r"; "k0r_lower_limit"; "k0r_upper_limit"; 
            "Er"; "Er_lower_limit"; "Er_upper_limit"; "R_squared"];
        value = [k0f; k0f_CI(1); k0f_CI(2); Ef; Ef_CI(1); Ef_CI(2);
            k0r; k0r_CI(1); k0r_CI(2); Er; Er_CI(1); Er_CI(2);
            r_squared];
        units = ["gal lbmol^-1^ min^-1^"; "gal lbmol^-1^ min^-1^";
            "gal lbmol^-1^ min^-1^"; "BTU/lbmol"; "BTU/lbmol"; 
            "BTU/lbmol"; "gal lbmol^-1^ min^-1^"; 
            "gal lbmol^-1^ min^-1^"; "gal lbmol^-1^ min^-1^";
            "BTU/lbmol"; "BTU/lbmol"; "BTU/lbmol"; ""];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'reb_20_5_3_results.csv')
        
        % create, show, and save a parity plot
        figure
        hold on
        plot(gMeasResp, CA_1_model,'ok','MarkerSize',10,'LineWidth',2)
        plot([min(gMeasResp),max(gMeasResp)], [min(gMeasResp) ...
            ,max(gMeasResp)],'r','LineWidth',2)
        hold off
        set(gca, 'FontSize', 14);
        xlabel('Experimental Outlet C_A (lbmol/gal)','FontSize', 14)
        ylabel('Predicted Outlet C_A (lbmol/gal)','FontSize', 14)
        legend({'Data','Parity Line'},'Location', 'northwest' ...
            ,'FontSize',14)
        saveas(gcf,'reb_20_5_3_parity.png')

        % create show, and save residuals plots
        figure
        plot(T - 459.7, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('T (°F)','FontSize', 14)
        ylabel('Residual (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_T_residual.png')

        figure
        plot(Vdot, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('VFR (gal/min)','FontSize', 14)
        ylabel('Residual (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_Vdot_residual.png')

        figure
        plot(CA_0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet C_A (lbmol/gal)','FontSize', 14)
        ylabel('Residual (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_CA0_residual.png')

        figure
        plot(CY_0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet C_Y (lbmol/gal)','FontSize', 14)
        ylabel('Residual (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_CY0_residual.png')

        figure
        plot(CZ_0, epsilon_expt,'ok','MarkerSize',10,'LineWidth',2)
        yline(0.0,'r','LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Inlet C_Z (lbmol/gal)','FontSize', 14)
        ylabel('Residual (lbmol/gal)','FontSize', 14)
        saveas(gcf,'reb_20_5_3_CZ0_residual.png')
    end

    % perform the calculations
    perform_the_calculations()
end