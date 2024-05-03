function reb_12_7_1
%REB_12_7_1 Reaction Engineering Basics Example 12.7.1
    % constants available to all functions
    % given
    k0_1 = 10.2; %  gal /mol /min
    k0_2 = 17.0; % gal /mol /min
    E_1 = 15300; % J /mol
    E_2 = 23700; % kJ /mol
    CA_in = 10; % mol /gal
    CB_in = 12; % mol /gal
    T_in = 350; % K
    V = 25; % gal
    Vdot_in = 12.5; % gal /min
    dH_1_298 = -12000;% J /mol
    dH_2_298 = -21300; % J /mol
    Cp_A = 85; % J /mol /K
    Cp_B = 125; % J /mol /K
    Cp_D = 200; % J /mol /K
    Cp_U = 170; % J /mol /K
    % known
    R = 8.314; % J /mol /K
    % calculated
    nA_in = CA_in*Vdot_in;
    nB_in = CB_in*Vdot_in;

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
        nD = guess(3);
        nU = guess(4);
        T = guess(5);

        % rates
        k_1 = k0_1*exp(-E_1/R/T);
        k_2 = k0_2*exp(-E_2/R/T);
        CA = nA/Vdot_in;
        CB = nB/Vdot_in;
        r_1 = k_1*CA*CB;
        r_2 = k_2*CA*CB;

        % heats of reaction
        dH1 = dH_1_298 + (Cp_D - Cp_A - Cp_B)*(T - 298.);
        dH2 = dH_2_298 + (Cp_U - Cp_A - Cp_B)*(T - 298.);

        % evaluate and return the residuals
        residual_1 = nA_in - nA + V*(-r_1 -r_2);
        residual_2 = nB_in - nB + V*(-r_1 -r_2);
        residual_3 = -nD + V*r_1;
        residual_4 = -nU + V*r_2;
        residual_5 = (nA_in*Cp_A + nB_in*Cp_B)*(T - T_in) ...
            + V*(r_1*dH1 + r_2*dH2);
        resids = [residual_1;residual_2;residual_3;residual_4;
            residual_5];
    end

    % function that performs the analysis
    function perform_the_analysis()
        % set the initial guess
        init_guess = [0.5*nA_in; 0.5*nA_in; 0.1*nA_in; 0.1*nA_in;
        T_in + 10.0];
    
        % solve the reactor design equations
        soln = unknowns(init_guess);
    
        % extract the individual values
        nA = soln(1);
        nD = soln(3);
        nU = soln(4);
        T = soln(5);
    
        % calculate the other quantities of interest
        fA = 100*(nA_in - nA)/nA_in;
        S_D_U = nD/nU;
    
        % tabulate the results
        item = ["Conversion";"Selectivity";"Temperature"];
        value = [fA; S_D_U; T];
        units = ["%";"mol D per mol U";"K"];
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