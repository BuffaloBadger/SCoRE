function using_Matlab_to_solve_DAEs()
% Calculations for Using Matlab to Solve DAEs
    % global constants
    nA0 = 2400; % mol
    T0 = 300; % K
    nZ0 = 0; % mol
    k0 = 1.44E10; % /min
    E = 15.3; % kcal /mol
    dH = -7; % kcal /mol
    Cp = 1; % kcal /L /K
    tf = 9; % min
    nAf = 480; % mol
    R = 1.987E-3; % kcal /mol /K

    % global variable
    g_V = nan;

    % derivatives function
    function ddt = IVODE_derivatives(ind, dep)
        % extract the dependent variables
        nA = dep(1);
        nZ = dep(2);
        T = dep(3);

        % calculate r
        k = k0*exp(-E/R/T);
        CA = nA/g_V;
        r = k*CA;

        % evaluate the derivatives
        dnAdt = -r*g_V;
        dnZdt = r*g_V;
        dTdt = -r*dH/Cp;

        % return an array containing the derivatives
        ddt = [dnAdt; dnZdt; dTdt];
    end

    % residuals function
    function epsilon = ATE_residuals(V)
        % set the initial values and stopping criterion for the IVODEs
        ind_0 = 0;
        dep_0 = [nA0; nZ0; T0];
        stop_var = 0;
        stop_val = tf;

        % make V available to the derivatives function
        g_V = V;

        % solve the IVODEs
        odes_are_stiff = false;
        [ind, dep, flag, message] = solve_ivodes(ind_0, dep_0...
            , stop_var, stop_val, @IVODE_derivatives, odes_are_stiff);
    
        % print a warning if there was a problem solving the design 
        % equations
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
            disp(['         ',message])
        end

        % extract the final value of nA
        nA = dep(:,1);
        nAfinal = nA(end);

        % evaluate the residual
        epsilon = nAfinal - nAf;
    end

    % BSTR model function
    function [V, t, nA, nZ, T] = BSTR_variables()
        % guess the fluid volume
        V_guess = 2400; % L

        % calculate the fluid volume
        options = optimoptions('fsolve','Display','off');
        [soln, ~, flag, details] = fsolve(@ATE_residuals, V_guess, options);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ATE solver did not converge')
            disp(['         ' , details.message])
        end

        % extract the fluid volume from the solution
        V = soln;

        % set the initial values and stopping criterion for the IVODEs
        ind_0 = 0;
        dep_0 = [nA0; nZ0; T0];
        stop_var = 0;
        stop_val = tf;

        % make V available to the derivatives function
        g_V = V;

        % solve the IVODEs
        odes_are_stiff = false;
        [ind, dep, flag, message] = solve_ivodes(ind_0, dep_0...
            , stop_var, stop_val, @IVODE_derivatives, odes_are_stiff);
    
        % print a warning if there was a problem solving the design 
        % equations
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
            disp(['         ',message])
        end

        % extract the profiles
        t = ind;
        nA = dep(:,1);
        nZ = dep(:,2);
        T = dep(:,3);
    end

    % deliverables function
    function deliverables()
        % solve the BSTR model equations
        [V, t, nA, ~, T] = BSTR_variables();

        % extract the final values
        t_final = t(end);
        nA_final = nA(end);
        T_final = T(end);

        % tabulate the results
        item = ["V";"t final";"nA final";"T final"];
        value = [V; t_final; nA_final; T_final];
        units = ["L";"min";"moles";"K"];
        resultsTable = table(item,value,units);

        % display the results
        disp(' ')
        disp(resultsTable)

        % save the results
        writetable(resultsTable,'results.csv')
    end

    % execution command
    deliverables();
end