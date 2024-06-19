function reb_16_3_2
    %REB_16_3_2 Calculations for Example 16.3.2 of Reaction Engineering Basics

    % constants available to all functions
	% given
    Vpfr = 4.0;
    nDotTot_0 = 1.25; % mol/s
    yA_0 = 0.5;
    yB_0 = 0.5;
    T_0 = 300; % K
    P = 2.5; % atm
    Cp_i = 25.8; % cal/mol/K
    dH = -8600; % cal/mol
    k_0 = 8.12E2; % /s
    E = 9500.0;
    UA = 13.6; % cal /K /s
    nDotY_0 = 0;
    nDotZ_0 = 0;
    % known
    Re = 1.987; % cal/mol/K
    Rw = 0.08206 * 1e-3; % m^3*atm/mol/K
	% calculated
    nDotA_0 = yA_0*nDotTot_0;
    nDotB_0 = yB_0*nDotTot_0;

    % residuals function for the ATEs
    function resid = residuals(guess)
        T_1 = guess(1);
        T_3 = guess(2);

        % solve the PFR design equations
        [~, nDotA, nDotB, nDotY, nDotZ, T] = profiles(T_1);
        nDotA_2 = nDotA(end);
        nDotB_2 = nDotB(end);
        nDotY_2 = nDotY(end);
        nDotZ_2 = nDotZ(end);
        T_2 = T(end);

        % calculate the heat transfer rate
        AMTD = 0.5*((T_3-T_0) + (T_2-T_1));
        Qdot = UA*AMTD;

        % evaluate the residuals
        eps1 = (nDotA_0 + nDotB_0 + nDotY_0 + nDotZ_0)*Cp_i*(T_1-T_0) ...
            - (nDotA_2+nDotB_2+nDotY_2+nDotZ_2)*Cp_i*(T_2-T_3);
        eps2 = (nDotA_0+nDotB_0+nDotY_0+nDotZ_0)*Cp_i*(T_1-T_0) - Qdot;
        resid =[eps1; eps2];
    end

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract ind and dep vars for the current integration step
        nDot_A = dep(1);
        nDot_B = dep(2);
        nDot_Y = dep(3);
        nDot_Z = dep(4);
        T = dep(5);

        % calculate rate
        nDotTot = nDot_A + nDot_B + nDot_Y + nDot_Z;
        CA = nDot_A*P/nDotTot/Rw/T;
        r = k_0*exp(-E/Re/T)*CA;

        % evaluate the derivatives
        dnDotAdV = -r;
        dnDotBdV = -r;
        dnDotYdV = r;
        dnDotZdV = r;
        dTdV = -r*dH/nDotTot/Cp_i;

        % return the derivatives
        derivs = [dnDotAdV; dnDotBdV; dnDotYdV; dnDotZdV; dTdV];
    end

    % reactor model
    function [V, nDotA, nDotB, nDotY, nDotZ, T] = profiles(T_1)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nDotA_0, nDotB_0, nDotY_0, nDotZ_0, T_1];

        % define the stopping criterion
        f_var = 0;
        f_val = Vpfr;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [V, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nDotA = dep(:,1);
        nDotB = dep(:,2);
        nDotY = dep(:,3);
        nDotZ = dep(:,4);
        T = dep(:,5);
    end

    % function that performs the analysis
	function perform_the_analysis()
        % solve the PFR reactor design equations
        [~, nDotA, ~, ~, ~, ~] = profiles(T_0);
        fA_PFR_only = 100*(nDotA_0 - nDotA(end))/nDotA_0;

        % initial guess for the unknowns for a low conversion steady state
        initial_guess = [T_0 + 1; T_0 + 2];

        % solve the heat exchanger energy balances
        [soln, flag, message] = solve_ates(@residuals, initial_guess);

        % check that the solution converged
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end

        % extract the results
        T_1 = soln(1);

        % solve the PFR reactor design equations
        [~, nDotA, ~, ~, ~, ~] = profiles(T_1);
        fA_low_conversion = 100*(nDotA_0 - nDotA(end))/nDotA_0;

        % repeat with a guess for a medium conversion steady state
        initial_guess = [T_0 + 40; T_0 + 80];
        [soln, flag, message] = solve_ates(@residuals, initial_guess);
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end
        T_1 = soln(1);
        [~, nDotA, ~, ~, ~, ~] = profiles(T_1);
        fA_medium_conversion = 100*(nDotA_0 - nDotA(end))/nDotA_0;

        % repeat with a guess for a high conversion steady state
        initial_guess = [T_0 + 100; T_0 + 200];
        [soln, flag, message] = solve_ates(@residuals, initial_guess);
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end
        T_1 = soln(1);
        [~, nDotA, ~, ~, ~, ~] = profiles(T_1);
        fA_high_conversion = 100*(nDotA_0 - nDotA(end))/nDotA_0;

        % tabulate the results
        item = ["Conversion without Thermal Backmixing"
            "Low Steady-State Conversion"; "Medium Steady-State Conversion"
            "High Steady-State Conversion"];
        value = [fA_PFR_only;fA_low_conversion; fA_medium_conversion; fA_high_conversion];
        units = ["%"; "%"; "%"; "%"];
        results_table = table(item,value,units);

        % display the results
        disp(results_table)

        % save the results
        writetable(results_table,"results.csv");
    end

    % perform the analysis
    perform_the_analysis()

end