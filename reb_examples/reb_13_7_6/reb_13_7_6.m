function reb_13_7_6
%REB_13_7_6 Calculations for Example 13.7.6 of Reaction Engineering Basics
    % constants available to all functions
	% given
    T = 400 + 273.15; % K
    k = 0.15; % /s
    eps = 0.4;
    Vdot_in = 1.0; % cm^3 /s
    mDot = 0.054; % g /s
    P_in = 30.0; % atm
    Rpart = 0.2; % cm
    Deff = 0.0029; % cm^2 /s
    mu = 0.022E-2; % g /s /cm
    D = 1.5; % cm
    fA = 0.85;
    sph = 1.0;
	% known
    R = 82.06; % cm^3 atm /mol /K
    Pconv = 9.872E-7; % atm cm^2 /dyne
	% calculated
    nA_in = Vdot_in*P_in/R/T;
    nA_out = nA_in*(1-fA);
    G = 4*mDot/pi()/D^2;
    phi = Rpart*sqrt(k/Deff);
    eta = 3/phi*(1/tanh(phi) - 1/phi);

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nZ = dep(2);
        P = dep(3);

        % calculate the rate
        CA = nA/(nA + nZ)*P/R/T;
        r = k*CA;

        % calculate the fluid density
        rho = mDot*P/(nA + nZ)/R/T;

        % evaluate the derivatives
        dnAdz = - pi()*D^2/4*(1-eps)*eta*r;
        dnZdz = pi()*D^2/4*(1-eps)*eta*r;
        dPdz = -Pconv*(1-eps)/eps^3*G^2/rho/sph/2/Rpart...
            *(150*(1-eps)*mu/sph/2/Rpart/G + 1.75);

        % return the derivatives
        derivs = [dnAdz; dnZdz; dPdz];
    end

    % reactor model
    function [z, nA, nZ, P] = profiles()
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_in; 0.0; P_in];

        % define the stopping criterion
        f_var = 1;
        f_val = nA_out;
        
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
        nZ = dep(:,2);
        P = dep(:,3);
    end

    % function that performs the analysis
	function perform_the_analysis()
        % solve the reactor design equations
        [z, ~, ~, P] = profiles();

        % calculate the other quantities of interest
        L = z(end);
        P_out = P(end);

        % tabulate the results
        item = ["effectiveness factor";"length";"outlet pressure"];
        value = [eta;L; P_out];
        units = ["";"cm";"atm"];
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