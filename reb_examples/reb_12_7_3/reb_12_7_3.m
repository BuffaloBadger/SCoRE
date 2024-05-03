function reb_12_7_3
%REB_12_7_3 Reaction Engineering Basics Example 12.7.3

    % constants available to all functions
    % given
    V = 0.5; % /m^3
    nA_in =  70; % mol /s
    nB_in = 1500; % mol /s
    Vdot_in = 40E-3; % m^3 /s
    k0 = 1.2E9; % m^3 /mol /s
    E = 25800*4.184; % J /mol
    K0 = 4.2E-18; % m^3 /mol
    dH = -22400*4.184; % J /mol
    Cp_A = 412; % J /mol /K
    Cp_B = 75.5; % J /mol /K
    % known
    R = 8.314; % J /mol /K

    % make T_in available in all functions
    T_in = nan;

    % reactor model
    function soln = unknowns(init_guess)
        % solve the ATEs
        [soln, flag, message] = solve_ates(@residuals, init_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(num2str(T_in - 273.15))
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    end

    % residuals function
    function resids = residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nB = guess(2);
        nZ = guess(3);
        T = guess(4);

        % rate
        k = k0*exp(-E/R/T);
        K = K0*exp(-dH/R/T);
        CA = nA/Vdot_in;
        CB = nB/Vdot_in;
        CZ = nZ/Vdot_in;
        r_1 = k*CA*CB*(1 - CZ/(K*CA*CB));

        % evaluate and return the residuals
        residual_1 = nA_in - nA - V*r_1;
        residual_2 = nB_in - nB - V*r_1;
        residual_3 = -nZ + V*r_1;
        residual_4 = -(nA_in*Cp_A + nB_in*Cp_B)*(T-T_in) - V*r_1*dH;
        resids = [residual_1;residual_2;residual_3;residual_4];
    end

    % function that performs the analysis
    function perform_the_analysis()

        % set a range of inlet temperatures
        T_in_range = linspace(75.,125.)' + 273.15;

        % set the initial guess
        init_guess = [0.5*nA_in; nB_in; 0.5*nA_in; T_in_range(1) + 5.];
    
        % allocate storage for fA vs. Tin
        fA_range = nan(100,1);
    
        % calculate fA for each Tin
        for iT = 1:100
            % make tau avialable
            T_in = T_in_range(iT);
    
            % solve the reactor design equations
            soln = unknowns(init_guess);
    
            % save nZ and use the result as the new guess
            fA_range(iT) = 100*(nA_in - soln(1))/nA_in;
            init_guess = soln;
        end
    
        % find tau that maximizes nZ
        [~, iMax] = max(fA_range);
    
        % plot fA vs Tin
        figure; 
        plot(T_in_range - 273.15,fA_range,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('T_in (K)','FontSize', 14)
        ylabel('Conversion of A (%)','FontSize', 14)
        saveas(gcf,"fA_vs_Tin.png")
    
        % solve the design equations using the optimum tau
        T_in = T_in_range(iMax);
        fA_guess = fA_range(iMax)/100.;
        init_guess = [(1-fA_guess)*nA_in; nB_in; fA_guess*nA_in;
            T_in_range(1) + 10.];
        soln = unknowns(init_guess);
    
        % extract the individual values
        nA = soln(1);
        T = soln(4) - 273.15;
        fA = 100*(nA_in - nA)/nA_in;
    
        % tabulate the results
        item = ["Opt Tin";"fA";"T out"];
        value = [T_in - 273.15; fA; T];
        units = ["°C";"%";"°C"];
        results_table = table(item,value,units);
    
        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,'results.csv');
    end

    % perform the analysis
    perform_the_analysis();
end