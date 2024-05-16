function reb_13_7_2
%REB_13_7_2 Calculations for Example 13.7.2 of Reaction Engineering Basics
    % constants available to all functions
	% given
    D = 10.0; % cm
    L = 500.0; % cm
    CA_in = 1.0E-3; % mol /cm^3
    CB_in = 1.2E-3; % mol /cm^3
    Vdot_in = 75E3; % mv^3 /min
    k0_1 = 8.72E8; % cm^3 /mol /min
    E_1 = 7200; % cal /mol
    dH_1 = -10700; % cal /mol
    Cp = 1.0; % cal /g /K
    rho = 1.0; % g \cm^3
    f_A = 0.95;
	% known
    R = 1.987; % cal /mol /K
	% calculated
    nA_in = Vdot_in*CA_in;
    nB_in = Vdot_in*CB_in;
    nA_out = nA_in*(1-f_A);

    % make [missing constant] available to all functions
    T_in = nan;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        T = dep(5);

        % calculate [quantities in the derivatives]
        k_1 = k0_1*exp(-E_1/R/T);
        CA = nA/Vdot_in;
        CB = nB/Vdot_in;
        r_1 = k_1*CA*CB;

        % evaluate the derivatives
        dnAdz = - pi()*D^2/4*r_1;
        dnBdz = - pi()*D^2/4*r_1;
        dnYdz = pi()*D^2/4*r_1;
        dnZdz = pi()*D^2/4*r_1;
        dTdz = - pi()*D^2/4*r_1*dH_1/(Vdot_in*rho*Cp);

        % return the derivatives
        derivs = [dnAdz; dnBdz; dnYdz; dnZdz; dTdz];
    end

    % reactor model
    function [z, nA, nB, nY, nZ, T] = profiles()
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_in; nB_in; 0.0; 0.0; T_in];

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
        nA = dep(:,1);
        nB = dep(:,2);
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
    end

    % implicit equation for IVODE initial value as residual
    function resid = residual(guess)
        % make the guess available to all functions
        T_in = guess;

        % solve the reactor design equations
        [z, nA, nB, nY, nZ, T] = profiles();

        % extract the calculated final [value]
        nA_f = nA(end);

        % evaluate and return the residual
        resid = nA_f - nA_out;
    end

    % function that performs the analysis
	function perform_the_analysis()

        % initial guess for [missing constant]
        initial_guess = 25 + 273.15;

        % calculate [missing constant]
        [T_in, flag, message] = solve_ates(@residual, initial_guess);

        % check that the solution converged
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end

        % solve the reactor design equations
        [z, nA, nB, nY, nZ, T] = profiles();
        
        % tabulate the results
        item = ["Inlet T";"Outlet T"];
        value = [T_in - 273.15; T(end) - 273.15];
        units = ["°C";"°C"];
        results_table = table(item,value,units);

        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,'results.csv');

        % display and save the graphs
    end

    % perform the analysis
    perform_the_analysis()
end