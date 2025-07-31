function using_Matlab_to_solve_ATEs
% Calculations for Using Matlab to Solve ATEs
    % global constants
    CA_0 = 1.0; % mol /L
    V = 50; % L
    CB_0 = 1.2; % mol /L
    rho = 1.0E3; % g /L
    Cp = 1.0; % cal /g /K
    T_0 = 303; % K
    dH = -10700; % cal /mol
    k0 = 8.72E5; % L /mol /min
    E = 7200; % cal /mol
    R = 1.987; % cal /mol /K

    % parameter
    Vdot_value = [75, 100]; % L /min

    % global variables
    g_Vdot = nan;

    % residuals function
    function epsilon = CSTR_residuals(guess)
        % extract the individual guesses
        nDotA_1 = guess(1);
        nDotB_1 = guess(2);
        nDotY_1 = guess(3);
        nDotZ_1 = guess(4);
        T_1 = guess(5);

        % calculate r
        k = k0*exp(-E/R/T_1);
        CA = nDotA_1/g_Vdot;
        CB = nDotB_1/g_Vdot;
        r = k*CA*CB;

        % evaluate the residuals
        epsilon_1 = g_Vdot*CA_0 - nDotA_1 - r*V;
        epsilon_2 = g_Vdot*CB_0 - nDotB_1 - r*V;
        epsilon_3 = - nDotY_1 + r*V;
        epsilon_4 = - nDotZ_1 + r*V;
        epsilon_5 = rho*g_Vdot*Cp*(T_1 - T_0) + r*V*dH;

        % return the residuals as a vector
        epsilon = [epsilon_1; epsilon_2; epsilon_3; epsilon_4; epsilon_5];
    end

    % CSTR function
    function [nA1, nB1, nY1, nZ1, T1] = CSTR_variables(Vdot)
        % define guesses for the CSTR variables
        nA1_guess = 0.9*Vdot*CA_0;
        nB1_guess = 0.9*Vdot*CB_0;
        nY1_guess = 0.0;
        nZ1_guess = 0.0;
        T1_guess = T_0 + 5.0;

        % create a vector containing the individual guesses
        guess = [nA1_guess; nB1_guess; nY1_guess; nZ1_guess; T1_guess];

        % make vDot globally available
        g_Vdot = Vdot;
        
        % solve the ATEs
        options = optimoptions('fsolve','Display','off');
        [soln, ~, flag, details] = fsolve(@CSTR_residuals, guess, options);
    
        % check that a solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ATE solver did not converge')
            disp(['         ' , details.message])
        end

        % extract the CSTR variables to be returned from the solution
        nA1 = soln(1);
        nB1 = soln(2);
        nY1 = soln(3);
        nZ1 = soln(4);
        T1 = soln(5);
    end

    % deliverables function
    function deliverables()
        % perform the calculations for each of the parameter values
        for i = 1:length(Vdot_value)
            % set Vdot
            Vdot = Vdot_value(i);

            % solve the ATEs
            [nA, nB, nY, nZ, T] = CSTR_variables(Vdot);

            % display the results
            disp(' ')
            disp(['Results for Vdot = ',num2str(Vdot,3), 'L/min:'])
            item = ["nA"; "nB"; "nY"; "nZ"; "T"];
            value = [nA; nB; nY; nZ; T];
            units = ["mol/min"; "mol/min"; "mol/min"; "mol/min"; "K"];
            results_table = table(item, value, units);
            disp(results_table);
        end
    end
        
    % execution command
    deliverables();
end