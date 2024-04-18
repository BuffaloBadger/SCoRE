function reb_J_6_2
%REB_J_6_2 Calculations for Example J.6.2 of Reaction Engineering Basics

    % given and known constants available to all functions
    T_out = 400.; % K
    D = 1. ; % in
    k_0 = 7.49E9*61.02; % in^3 /mol /min
    E = 15300.; % cal /mol
    P = 4.; % atm
    dH = -14500.; % cal /mol
    Cp_A = 10.9; % cal /mol /K
    Cp_Z = 21.8; % cal /mol /K
    L = 100.; % in
    nDot_A_in = 1.5; % mol /min
    nDot_Z_in = 0.0;
    Re = 1.987; % cal /mol /K
    Rw = 0.08206*61.02; % in^3 atm /mol /K

    % make T_in available in all functions
    T_in = nan;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract ind and dep vars for this integration step
        nDot_A = dep(1);
        nDot_Z = dep(2);
        T = dep(3);

        % calculate rate
        r = k_0*exp(-E/Re/T)*(nDot_A*P/Rw/T/(nDot_A + nDot_Z))^2;

        % evaluate the derivatives
        dnDotAdz = -2*pi()*D^2/4*r;
        dnDotZdz = pi()*D^2/4*r;
        dTdz = -pi()*D^2/4*r*dH/(nDot_A*Cp_A + nDot_Z*Cp_Z);

        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];

    end

    % residual function
    function resid = residual(guess)
        % make the guess available to all functions
        T_in = guess;

        % solve the reactor design equations
        [~, ~, ~, T] = profiles();

        % extract the calculated final temperature
        T_f = T(end);

        % evaluate and return the residual
        resid = T_f - T_out;
    end

    % reactor model function
    function [z, nDot_A, nDot_Z, T] = profiles()
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nDot_A_in; nDot_Z_in; T_in];

        % stopping criterion
        f_var = 0;
        f_val = L;

        % solve the IVODEs
        odes_are_stiff = false;
        [z, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % Check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nDot_A = dep(:,1);
        nDot_Z = dep(:,2);
        T = dep(:,3);
    end

    % function that perform the analysis
    function perform_the_analysis()

        % calculate the inlet temperature
        % initial guess
        initial_guess = T_out - 100.0;
    
        % solve the implicit equation for Tin
        [T_in, flag, message] = solve_ates(@residual, initial_guess);
    
        % check that the solution converged
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end
    
        % solve the reactor design equations
        [z, nDot_A, nDot_Z, T] = profiles();
    
        % tabulate the results
        item = "T_in";
        value = T_in;
        units = "K";
        Tin_results_table = table(item,value,units);
        profile_results_table = table(z, nDot_A, nDot_Z, T);
    
        % display the results
        disp(' ')
        disp(['Inlet Temperature: ', num2str(T_in,3), ' K'])
        disp(' ')
        disp('Molar flow and Temperature Profiles')
        disp(profile_results_table)
    
        % Save the results
        Tin_results_file = "Tin_results.csv";
        writetable(Tin_results_table,Tin_results_file);
        profile_results_file ...
            = "profile_results.csv";
        writetable(profile_results_table,profile_results_file);
    end

    % perform the analysis
    perform_the_analysis()
end