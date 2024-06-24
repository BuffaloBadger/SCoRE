function REB_17_3_1
    %REB_17_3_1 Calculations for Example 17.3.1 of Reaction Engineering Basics

    % constants available to all functions
	% given
    dH_1 = -14000.; % cal /mol
    k0_1 = 4.2E15; % cm^3 /mol /min
    E_1 = 18000.; % cal /mol
    Cp = 1.3; % cal /cm^3 /K
    CA_0 = 2.0E-3; % mol /cm^3
    CZ_0 = 0.0; % mol /cm^3
    Vdot_0 = 500.; % cm^3 /min
    T_0 = 300.; % K
    R_R = 1.3; %
    D = 5.; % cm
    L = 50.; % cm
    % known
    R = 1.987; % cal /mol /K

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract ind and dep vars for the current integration step
        nDot_A = dep(1);
        nDot_Z = dep(2);
        T = dep(3);

        % calculate other unknown quantities
        Vdot_3 = Vdot_0;
        Vdot_4 = R_R*Vdot_3;
        Vdot = Vdot_3 + Vdot_4;
        k_1 = k0_1*exp(-E_1/R/T);
        C_A = nDot_A/Vdot;
        C_Z = nDot_Z/Vdot;
        r_1 = k_1*C_A*C_Z;

        % evaluate the derivatives
        dnDotAdz = -pi()*D^2/4*r_1;
        dnDotZdz = pi()*D^2/4*r_1;
        dTdz = -pi()*D^2/4*r_1*dH_1/Vdot/Cp;

        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];
    end

    % reactor model
    function [z, nDotA, nDotZ, T] = profiles(dep_0)
        % set the initial values
        ind_0 = 0.0;

        % define the stopping criterion
        f_var = 0;
        f_val = L;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [z, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
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

    % stream mixer residuals function
    function resid = residuals(guess)
        % extract the individual guesses
        nDotA_1 = guess(1);
        nDotZ_1 = guess(2);
        T_1 = guess(3);

        % calculate other unknown constants
        nDotA_0 = Vdot_0*CA_0;
        nDotZ_0 = Vdot_0*CZ_0;
        Vdot_3 = Vdot_0;
        Vdot_4 = R_R*Vdot_3;

        % solve the reactor design equations
        initial_values = [nDotA_1; nDotZ_1; T_1];
        [~, nDotA, nDotZ, T] = profiles(initial_values);

        % extract final values
        nDotA_2 = nDotA(end);
        nDotZ_2 = nDotZ(end);
        T_4 = T(end);

        % calculate molar recycle flows
        nDotA_4 = R_R/(1 + R_R)*nDotA_2;
        nDotZ_4 = R_R/(1 + R_R)*nDotZ_2;

        % evaluate the residuals
        eps_1 = nDotA_1 - nDotA_0 - nDotA_4;
        eps_2 = nDotZ_1 - nDotZ_0 - nDotZ_4;
        eps_3 = Vdot_0*Cp*(T_1 - T_0) + Vdot_4*Cp*(T_1 - T_4);

        % return the residuals
        resid = [eps_1; eps_2; eps_3];
    end

    % stream mixer model
    function [nDotA_1, nDotZ_1, T_1] = unknowns(initial_guess)

        % solve the stream mixer mole and energy balances
        [soln, flag, message] = solve_ates(@residuals, initial_guess);

        % check that the solution converged
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end

        % extract and return the results
        nDotA_1 = soln(1);
        nDotZ_1 = soln(2);
        T_1 = soln(3);
    end

    % function that performs the analysis
	function perform_the_analysis()

        % initial guess for the unknowns
        nDotA_0 = Vdot_0*CA_0;
        initial_guess = [1.2*nDotA_0; nDotA_0; T_0 + 10];

        % solve the stream mixer balances
        [nDotA_1, nDotZ_1, T_1] = unknowns(initial_guess);

        % solve the PFR design equations
        dep_0 = [nDotA_1; nDotZ_1; T_1];
        [~, nDotA, nDotZ, T] = profiles(dep_0);
        
        % extract the outlet values
        nDotA_2 = nDotA(end);
        nDotZ_2 = nDotZ(end);
        T_3 = T(end);

        % calculate molar recycle flows
        nDotA_4 = R_R/(1 + R_R)*nDotA_2;
        nDotZ_4 = R_R/(1 + R_R)*nDotZ_2;

        % calculate the product molar flows
        nDotA_3 = nDotA_2 - nDotA_4;
        nDotZ_3 = nDotZ_2 - nDotZ_4;

        % calculate the product concentrations
        Vdot_3 = Vdot_0;
        CA_3 = nDotA_3/Vdot_3*1000;
        CZ_3 = nDotZ_3/Vdot_3*1000;

        % tabulate the results
        item = ["nDotA 1"; "nDotZ 1"; "T 1"; "CA 3"; "CZ 3"; "T 3"];
        value = [nDotA_1; nDotZ_1; T_1; CA_3; CZ_3; T_3];
        units = ["mol/min"; "mol/min"; "K"; "M"; "M"; "K"];
        results_table = table(item,value,units);

        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,"results.csv");
    end

    % perform the analysis
    perform_the_analysis()
end