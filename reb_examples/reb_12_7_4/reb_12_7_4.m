function reb_12_7_4
%REB_12_7_4 Reaction Engineering Basics Example 12.7.4

    % constants available to all functions
    % given
    V = 500; % cm^3
    Vdot_in = 1.0; % cm^3 /s
    CA_in = 0.015; % mol /cm^3
    CB_in = 0.015; % mol /cm^3
    T_in = 50 + 273.15; % K
    Cp = 0.35*4.184; % J /g /K
    rho = 0.93; % g /cm^3
    dH = -20000; % J /mol
    k0 = 3.24E12; % cm^3 /mol /s
    E = 105000; % J/mol
    % known
    R = 8.314; % J /mol /K
    % calculated
    nA_in = Vdot_in*CA_in;
    nB_in = Vdot_in*CB_in;

    % reactor model
    function soln = unknowns(init_guess)
        % solve the ATEs
        [soln, flag, message] = solve_ates(@residuals, init_guess);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp(['WARNING: The ATE solver did not converge: ',message])
        end
    end

    % residuals function
    function resids = residuals(guess)
        % extract the individual guesses
        nA = guess(1);
        nB = guess(2);
        nY = guess(3);
        nZ = guess(4);
        T = guess(5);

        % rate
        k = k0*exp(-E/R/T);
        CA = nA/Vdot_in;
        CB = nB/Vdot_in;
        r_1 = k*CA*CB;

        % evaluate and return the residuals
        residual_1 = nA_in - nA - V*r_1;
        residual_2 = nB_in - nB - V*r_1;
        residual_3 = -nY + V*r_1;
        residual_4 = -nZ + V*r_1;
        residual_5 = -Vdot_in*rho*Cp*(T-T_in) - V*r_1*dH;
        resids = [residual_1;residual_2;residual_3;residual_4;
            residual_5];
    end

    % function that performs the analysis
    function perform_the_analysis()

        % set an initial guess for the low-conversion steady state
        init_guess = [0.99*nA_in; 0.99*nB_in; 0.01*nA_in; 0.01*nA_in;
            T_in + 1.];
    
        % solve the reactor design equations
        soln = unknowns(init_guess);

        % save the conversion and temperature
        fA_low = 100*(nA_in - soln(1))/nA_in;
        T_low = soln(5) - 273.15;

        % set an initial guess for the high-conversion steady state
        init_guess = [0.01*nA_in; 0.01*nB_in; 0.99*nA_in; 0.99*nA_in;
            T_in + 220.];
    
        % solve the reactor design equations
        soln = unknowns(init_guess);

        % save the conversion and temperature
        fA_high = 100*(nA_in - soln(1))/nA_in;
        T_high = soln(5) - 273.15;

        % set an initial guess for the mid-conversion steady state
        fA_guess = (fA_high + fA_low)/100/2;
        nA_guess = nA_in*(1-fA_guess);
        nY_guess = nA_in*fA_guess;
        init_guess = [nA_guess; nA_guess; nY_guess; nY_guess;
            (T_low + T_high)/2 + 273.15 - 25];
    
        % solve the reactor design equations
        soln = unknowns(init_guess);

        % save the conversion and temperature
        fA_mid = 100*(nA_in - soln(1))/nA_in;
        T_mid = soln(5) - 273.15;
    
        % tabulate the results
        conversion = [fA_low; fA_mid; fA_high];
        temperature = [T_low; T_mid; T_high];
        results_table = table(conversion, temperature);
    
        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,'results.csv');
    end

    % perform the analysis
    perform_the_analysis();
end