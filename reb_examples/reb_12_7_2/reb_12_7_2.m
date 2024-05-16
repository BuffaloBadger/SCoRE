function reb_12_7_2
%REB_12_7_2 Reaction Engineering Basics Example 12.7.2
    % constants available to all functions
    % given
    yA_in = 0.1;
    yB_in = 0.65;
    yI_in = 0.25;
    T_in = 165 + 273.15; % K
    P = 5; % atm
    yA = 0.001;
    k0_1 = 1.37E5; % m^3 /mol /min
    E_1 = 11100; % cal /mol
    dH_1 = -7200; % cal /mol
    Cp_A = 7.6; % cal /mol K
    Cp_B = 8.2; % cal /mol K
    Cp_I = 4.3; % cal /mol K
    % known
    Re = 1.987; % cal /mol K 
    Rw = 8.206E-5; % m^3 atm /mol /K
    % basis
    Vdot_in = 1.0; % m^3 /min
    % calculated
    nA_in = yA_in*P*Vdot_in/Rw/T_in;
    nB_in = yB_in*P*Vdot_in/Rw/T_in;
    nI_in = yI_in*P*Vdot_in/Rw/T_in;

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
        V = guess(1);
        nB = guess(2);
        nI = guess(3);
        nZ = guess(4);
        T = guess(5);
        
        % calculate nA
        nA = yA*(nB + nI + nZ)/(1 - yA);

        % rates
        k_1 = k0_1*exp(-E_1/Re/T);
        CA = nA/(nA + nB + nI + nZ)*P/Rw/T;
        CB = nB/(nA + nB + nI + nZ)*P/Rw/T;
        r_1 = k_1*CA*CB;

        % evaluate and return the residuals
        residual_1 = nA_in - nA - V*r_1;
        residual_2 = nB_in - nB - V*r_1;
        residual_3 = nI_in - nI;
        residual_4 = -nZ + V*r_1;
        residual_5 = -(nA_in*Cp_A + nB_in*Cp_B + nI_in*Cp_I)...
            *(T - T_in) - V*r_1*dH_1;
        resids = [residual_1;residual_2;residual_3;residual_4;
            residual_5];
    end

    % function that performs the analysis
    function perform_the_analysis()
    
        % set the initial guess
        init_guess = [1.0, nB_in - 0.01*nA_in, nI_in, 0.01*nA_in...
            , T_in + 50.0];
    
        % solve the reactor design equations
        soln = unknowns(init_guess);
    
        % extract the individual values
        V = soln(1);
        T = soln(5) - 273.15;
    
        % calculate the other quantities of interest
        tau = V/Vdot_in;
    
        % tabulate the results
        item = ["tau";"T"];
        value = [tau; T];
        units = ["min";"°C"];
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