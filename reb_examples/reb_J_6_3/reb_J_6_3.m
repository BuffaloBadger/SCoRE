function reb_J_6_3
%REB_J_6_3 Calculations for Example J.6.3 of Reaction Engineering Basics

    % given and known constants available in all functions
    T_f = 325.; % K
    D = 5.; % cm
    dH = -14000.; % cal /mol
    Cp = 1.3; % cal /cm^3 /K
    C_A_in = 0.0025; % mol /cm^3
    C_Z_in = 0.0; % mol /min
    T_in = 300.; % K
    L = 50.0; % cm
    k_0 = 4.2E15; % cm^3 /mol /min
    E = 18000.; % cal /mol
    Re = 1.987; % cal /mol /K
    
    % make Vdot available to all functions
    Vdot = nan;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nDot_A = dep(1);
        % nDot_Z = dep(2); not needed to evaluate derivatives
        T = dep(3);
    
        % calculate rate
        r = k_0*exp(-E/Re/T)*(nDot_A^2)/Vdot^2;
    
        % evaluate the derivatives
        dnDotAdz = -pi()*D^2/4*r;
        dnDotZdz = pi()*D^2/4*r;
        dTdz = -pi()*D^2/4*r*dH/Vdot/Cp;
    
        % return the derivatives
        derivs = [dnDotAdz; dnDotZdz; dTdz];
    end

    % residual function
    function resid = residual(guess)
        % set Vdot
        Vdot = guess;

        % solve the reactor design equations
        [~, ~, ~, T] = profiles();

        % extract the calculated final T
        T_f_calc = T(end);

        % evaluate and return the residual
        resid = T_f - T_f_calc;
    end

    % reactor model function
    function [z, nDot_A, nDot_Z, T] = profiles()
        % set the intial values
        nDot_A_in = Vdot*C_A_in;
        nDot_Z_in = Vdot*C_Z_in;
        ind_0 = 0.0;
        dep_0 = [nDot_A_in; nDot_Z_in; T_in];

        % define the stopping criterion
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

    % function that performs the analysis
    function perform_the_analysis()
        % initial guess for Vdot
        initial_guess = 1150.;
    
        % solve the implicit equation for Vdot
        [Vdot, flag, message] = solve_ates(@residual, initial_guess);
    
        % check that the solution converged
        if flag <= 0
            disp(' ')
            disp(['The ATE solver did not converge: ',message])
        end
    
        % solve the reactor design equations
        [z, nDot_A, nDot_Z, T] = profiles();
    
        % tabulate the results
        item = "Vdot";
        value = Vdot;
        units = "cm^3^ min^-1^";
        Tin_results_table = table(item,value,units);
        profile_results_table = table(z, nDot_A, nDot_Z, T);
    
        % display the results
        disp(' ')
        disp(['Volumetric Flow Rate: ', num2str(Vdot,3), ' cm^3/min'])
        disp(' ')
        disp('Molar flow and Temperature Profiles')
        disp(profile_results_table)
    
        % Save the results
        Vdot_results_file = "Vdot_results.csv";
        writetable(Tin_results_table,Vdot_results_file);
        profile_results_file ...
            = "profile_results.csv";
        writetable(profile_results_table,profile_results_file);
    end

    % perform the analysis
    perform_the_analysis();
end