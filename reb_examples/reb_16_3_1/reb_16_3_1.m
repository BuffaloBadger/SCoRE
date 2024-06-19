function reb_16_3_1
    %REB_16_3_1 Calculations for Example 16.1.3 of Reaction Engineering Basics

    % constants available to all functions
	% given
    Vdot = 750.0; % L /min
    CA_0 = 3.8; % mol /L
    T_0 = 25 + 273.15; % K
    UA = 1500.0; % kJ /min /K
    dH = -79.8; % kJ/mol
    Cp = 987.0 * 4.184e-3; % kJ /L /K
    fA = 0.8;
    k_0 = 3.38E6; % /min
    E = 50.0; % kJ ;mol
    % known
    R = 8.31446e-3; % kJ /mol /K
	% calculated
    nDotA_0 = Vdot*CA_0;

    % residuals function for the ATEs
    function resid = residuals(guess)
        T_1 = guess(1);
        T_3 = guess(2);

        % solve the PFR design equations and extract T_2
        [~, ~, ~, T] = profiles(T_1);
        T_2 = T(end);

        % calculate the heat transfer rate
        if T_3 - T_0 == T_2 - T_1
            LMTD = T_3 - T_0;
        else
            LMTD = ((T_3 - T_0) - (T_2 - T_1))/log((T_3 - T_0)/(T_2 - T_1));
        end
        Qdot = UA*LMTD;

        % evaluate the residuals
        resid =[Vdot*Cp*(T_1-T_0) - Vdot*Cp*(T_2-T_3)
            Vdot*Cp*(T_1-T_0) - Qdot
            ];
    end

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract ind and dep vars for the current integration step
        nDot_A = dep(1);
        T = dep(3);

        % calculate rate
        r = k_0*exp(-E/R/T)*nDot_A/Vdot;

        % evaluate the derivatives
        dnDotAdV = -r;
        dnDotZdV = r;
        dTdV = -r*dH/Vdot/Cp;

        % return the derivatives
        derivs = [dnDotAdV; dnDotZdV; dTdV];
    end

    % reactor model
    function [V, nDotA, nDotZ, T] = profiles(T_1)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nDotA_0, 0.0, T_1];

        % define the stopping criterion
        f_var = 1;
        f_val = nDotA_0*(1-fA);
        
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
        nDotZ = dep(:,2);
        T = dep(:,3);
    end

    % function that performs the analysis
	function perform_the_analysis()

        % case a, PFR only
        % solve the PFR reactor design equations
        [V, ~, ~, T] = profiles(T_0);
        V_pfr_only = V(end);
        T_pfr_only = T(end) - 273.15;

        % case b, thermally backmixed PFR
        % initial guess for the unknowns
        initial_guess = [T_0 + 25; T_0 + 50];

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
        [V, ~, ~, T] = profiles(T_1);
        V_tb_pfr = V(end);
        T_tb_pfr = T(end) - 273.15;

        % tabulate the results
        item = ["PFR Volume"; "PFR Outlet Temperature"
            "Thermally Backmixed PFR Volume"
            "Thermally Backmixed PFR Outlet Temperature"];
        value = [V_pfr_only; T_pfr_only; V_tb_pfr; T_tb_pfr];
        units = ["L"; "°C"; "L"; "°C"];
        results_table = table(item,value,units);

        % display the results
        disp(results_table)

        % save the results
        writetable(results_table,"results.csv");
    end

    % perform the analysis
    perform_the_analysis()

end