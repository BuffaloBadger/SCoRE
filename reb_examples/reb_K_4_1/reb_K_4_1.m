function reb_K_4_1
    %REB_K_4_1 Calculations for Example K.4.1 of Reaction Engineering Basics

    % constants available to all functions
	% given
    D = 5.; % cm
    dH = -14000.; % cal /mol
    Vdot_feed = 500.; % cm^3 /min
    Cp = 1.3; % cal /cm^3 /K
    R_R = 1.3; %
    nDotA_feed = 1.0; % mol /min
    nDotZ_feed = 0.0;
    L = 50.; % cm
    T_feed = 300.; % K
    k_0 = 4.2E15; % cm^3 /mol /min
    E = 18000.; % cal /mol
    % known
    R = 1.987; % cal /mol /K
	% calculated
    Vdot_prod = Vdot_feed;
    Vdot_r = R_R*Vdot_prod;
    Vdot_in = Vdot_feed + Vdot_r;

    % other equipment model
    function [nDotA_in, nDotZ_in, T_in, T_r, nDotA_prod, nDotZ_prod...
            , T_prod] = unknowns(initial_guess)

        % solve the other equipment mole and energy balances
        [soln, flag, message] = solve_ates(@residuals, initial_guess);

        % check that the solution converged
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end

        % extract and return the results
        nDotA_in = soln(1);
        nDotZ_in = soln(2);
        T_in = soln(3);
        T_r = soln(4);
        nDotA_prod = soln(5);
        nDotZ_prod = soln(6);
        T_prod = soln(7);
    end

    % residuals function for the ATEs
    function resid = residuals(guess)
        % extract the individual guesses
        nDotA_in = guess(1);
        nDotZ_in = guess(2);
        T_in = guess(3);
        T_r = guess(4);
        nDotA_prod = guess(5);
        nDotZ_prod = guess(6);
        T_prod = guess(7);

        % solve the reactor design equations
        initial_values = [nDotA_in; nDotZ_in; T_in];
        [~, nDotA, nDotZ, T] = profiles(initial_values);

        % extract final values
        nDotA_out = nDotA(end);
        nDotZ_out = nDotZ(end);
        T_out = T(end);

        % calculate molar recycle flows
        nDotA_r = R_R*nDotA_prod;
        nDotZ_r = R_R*nDotZ_prod;

        % evaluate the residuals
        eps_1 = nDotA_feed + nDotA_r - nDotA_in;
        eps_2 = nDotZ_feed + nDotZ_r - nDotZ_in;
        eps_3 = Vdot_feed*Cp*(T_in - T_feed) + Vdot_r*Cp*(T_in - T_r);
        eps_4 = nDotA_out - nDotA_r - nDotA_prod;
        eps_5 = nDotZ_out - nDotZ_r - nDotZ_prod;
        eps_6 = T_out - T_r;
        eps_7 = T_out - T_prod;

        % return the residuals
        resid = [eps_1; eps_2; eps_3; eps_4; eps_5; eps_6; eps_7];
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

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract ind and dep vars for the current integration step
        nDot_A = dep(1);
        nDot_Z = dep(2);
        T = dep(3);

        % calculate rate
        r = k_0*exp(-E/R/T)*nDot_A*nDot_Z/Vdot_in^2;

        % evaluate the derivatives
        dnDotAdz = -pi()*D^2/4*r;
        dnDotZdz = pi()*D^2/4*r;
        dTdz = -pi()*D^2/4*r*dH/Vdot_in/Cp;

        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];
    end

    % function that performs the analysis
	function perform_the_analysis()

        % initial guess for the unknowns
        initial_guess = [0.9*nDotA_feed; 0.1*nDotA_feed; T_feed + 5
            T_feed + 10; 0.1*nDotA_feed; 0.1*nDotA_feed; T_feed + 10];

        % calculate the unknowns
        [nDotA_in, nDotZ_in, T_in, T_r, nDotA_prod, nDotZ_prod...
            , T_prod] = unknowns(initial_guess);

        % tabulate the results
        item = ["A in"; "Z in"; "T in"; "T r"; "A prod"; "Z prod"
            "T prod"];
        value = [nDotA_in; nDotZ_in; T_in; T_r; nDotA_prod; nDotZ_prod
            T_prod];
        units = ["mol min^-1^"; "mol min^-1^"; "K"; "K"; "mol min^-1^"
            "mol min^-1^"; "K"];
        results_table = table(item,value,units);

        % display the results
        disp(results_table)

        % save the results
        writetable(results_table,"results.csv");
    end

    % perform the analysis
    perform_the_analysis()

end