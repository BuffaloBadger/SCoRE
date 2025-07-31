function using_Matlab_to_estimate_parameters
% Calculations for Using Matlab to Estimate Parameters
    % global constants
    L = 30.0; % cm
    D = 2.0; % cm
    R = 8.3145E-3; % kJ/mol/K

    % read the data and separate into adjusted inputs and experimental
    % responses
    dataTable = readtable('estimating_parameters_data.csv'...
        , 'VariableNamingRule','preserve');
    adjExptInputs = table2array(dataTable(:,1:4));
    adjExptInputs(:,1) = adjExptInputs(:,1) + 273.15; % convert °C to K
    adjExptInputs(:,3) = adjExptInputs(:,3)/1000.0; % convert to mol/cc
    adjExptInputs(:,4) = adjExptInputs(:,4)/1000.0; % convert to mol/cc
    exptResp = table2array(dataTable(:,5));

    % extract individual adjusted inputs
    T = adjExptInputs(:,1);
    VFR = adjExptInputs(:,2);
    CA0 = adjExptInputs(:,3);
    CB0 = adjExptInputs(:,4);

    % determine the number of experiments
    nExpt = length(exptResp);

    % global variables
    g_k0 = nan;
    g_E = nan;
    g_iExpt = -1;

    % derivatives function
    function ddz = derivatives(~,dep)
        % get dependent variables that are needed
        nDotA = dep(1);
        nDotB = dep(2);
        nDotZ = dep(3);

        % calculate the rate
        k = g_k0*exp(-g_E/R/T(g_iExpt));
        CA = nDotA/VFR(g_iExpt);
        CB = nDotB/VFR(g_iExpt);
        r = k*CA*CB;

        % evaluate the derivatives
        dnAdz = -pi*D^2/4*r;
        dnBdz = -pi*D^2/4*r;
        dnZdz = 2*pi*D^2/4*r;

        % return the derivatives as a column vector
        ddz = [dnAdz; dnBdz; dnZdz];
    end

    % PFR function
    function [z, nDotA, nDotB, nDotZ] = PFR_variables(iExpt)
        % define initial values
        ind_0 = 0.0;
        nDotA0 = CA0(iExpt)*VFR(iExpt);
        nDotB0 = CB0(iExpt)*VFR(iExpt);
        nDotZ0 = 0.0;
        dep_0 = [nDotA0; nDotB0; nDotZ0];
        
        % define stopping criterion
        stop_var = 0.0;
        stop_val = L;

        % make iExpt available to derivatives function
        g_iExpt = iExpt;

        % solve the design equations
        odes_are_stiff = false;
        [ind, dep, flag, message] = solve_ivodes(ind_0, dep_0...
            , stop_var, stop_val, @derivatives, odes_are_stiff);
    
        % print a warning if there was a problem solving the design 
        % equations
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
            disp(['         ',message])
        end
        
        % extract the corresponding sets of values for each of the
        % variables
        z = ind;
        nDotA = dep(:,1);
        nDotB = dep(:,2);
        nDotZ = dep(:,3);
    end

    % predicted responses function
    function modelResp = predicted_responses(params,~)
        % extract the parameters and make them globally available
        g_k0 = params(1);
        g_E = params(2);

        % Allocate storage for the predicted response
        modelResp = nan(nExpt,1);

        % loop through all of the experiments
        for iExpt = 1:nExpt
            % solve the PFR design equations
            [~, nDotA, ~, ~] = PFR_variables(iExpt);

            % extract the final value of nDotA
            nDotA_1 = nDotA(end);

            % calculate the predicted response
            nDotA_0 = CA0(iExpt)*VFR(iExpt);
            modelResp(iExpt) = (nDotA_0 - nDotA_1)/nDotA_0;
        end
    end

    % parameter estimation function
    function [k0, k0_CI, E, E_CI, rSquared] = parameter_estimates()
        % define guesses for the parameters
        k0guess = 1.0E9;
        Eguess = 30.0;
        guess = [k0guess; Eguess];

        % fit the model to the data
        useRelErr = false;
        [beta, betaCI, rSquared] = fit_to_SR_data(guess, adjExptInputs,...
            exptResp, @predicted_responses, useRelErr);

        % extract the parameters
        k0 = beta(1);
        k0_CI = betaCI(1,:);
        E = beta(2);
        E_CI = betaCI(2,:);
    end

    % deliverables function
    function deliverables()
        % estimate the parameters
        [k0, k0_CI, E, E_CI, rSquared] = parameter_estimates();

        % extract the confidence intervals
        k0_ll = k0_CI(1);
        k0_ul = k0_CI(2);
        E_ll = E_CI(1);
        E_ul = E_CI(2);

        % report the results
        disp(' ')
        disp(['k0 = ',num2str(k0,4),' cc/mol/min [',num2str(k0_ll,4)...
            ,', ', num2str(k0_ul,4),']'])
        disp(['E = ',num2str(E,4),' kJ/mol [',num2str(E_ll,4),', ',...
            num2str(E_ul,4),']'])
        disp(['R = ',num2str(rSquared,3)])
    
        % Save the results to a .csv file
        item = ["k0"; "k0_lower_limit"; "k0_upper_limit"; "E";...
            "E_lower_limit"; "E_upper_limit"; "R-squared"];
        value = [k0; k0_ll; k0_ul; E; E_ll; E_ul; rSquared];
        units = ["cm^3^ mol^-1^ min^-1^"; "cm^3^ mol^-1^ min^-1^"; ...
            "cm^3^ mol^-1^ min^-1^"; "kJ mol^-1^"; "kJ mol^-1^"; ...
            "kJ mol^-1^"; " "];
        results_table = table(item,value,units);
        writetable(results_table,'results.csv');

        % generate, show, and save a parity plot
        estParams = [k0, E];
        modelResp = predicted_responses(estParams,adjExptInputs);
        parityRange = [min(exptResp),max(exptResp)];
        figure;
        plot(exptResp,modelResp,'ko',parityRange,parityRange,'r'...
            ,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Experimental Conversion','FontSize', 14)
        ylabel('Predicted Conversion','FontSize', 14)
        legend({'Data','Parity Line'},'Location','northwest'...
            ,'FontSize',14)
        saveas(gcf,"parity.pdf")

        % generate, show, and save residuals plots
        epsilonExpt = modelResp - exptResp;
        figure;
        plot(T,epsilonExpt,'ko','LineWidth',2)
        yline(0, 'r', 'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('T (K)','FontSize', 14)
        ylabel('Residual','FontSize', 14)
        saveas(gcf,"residual_T.pdf")

        figure;
        plot(VFR,epsilonExpt,'ko','LineWidth',2)
        yline(0, 'r', 'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Volumetric Flow Rate (cm^3 min^-^1)','FontSize', 14)
        ylabel('Residual','FontSize', 14)
        saveas(gcf,"residual_VFR.pdf")

        figure;
        plot(CA0,epsilonExpt,'ko','LineWidth',2)
        yline(0, 'r', 'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('C_A_0 (mol cm^-^3)','FontSize', 14)
        ylabel('Residual','FontSize', 14)
        saveas(gcf,"residual_CA0.pdf")

        figure;
        plot(CB0,epsilonExpt,'ko','LineWidth',2)
        yline(0, 'r', 'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('C_B_0 (mol cm^-^3)','FontSize', 14)
        ylabel('Residual','FontSize', 14)
        saveas(gcf,"residual_CB0.pdf")
    end

    % execution command
    deliverables();
end