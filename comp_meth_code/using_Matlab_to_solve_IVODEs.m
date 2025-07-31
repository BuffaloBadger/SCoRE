function using_Matlab_to_solve_IVODEs
% Calculations for Using Matlab to Solve IVODEs
    % global constants
    V = 1.0; % L
    tf = 30.0; % min
    CA0 = 0.5; % mol/L
    CZ0 = 0.0; % mol/L
    T = 338; % K
    k0 = 3.6E8; % L/mol/min
    E = 67.5; % kJ/mol
    R = 8.314E-3; % kJ/mol/K

    % derivatives function
    function ddt = BSTR_derivatives(t, dep)
        % extract the dependent variables
        nA = dep(1);
        nZ = dep(2);

        % calculate the rate
        k = k0*exp(-E/R/T);
        CA = nA/V;
        r = k*CA;
    
        % evaluate the derivatives of the dependent variables
        dnAdt = -r*V;
        dnZdt = r*V;
    
        % return an array containing the derivatives
        ddt = [dnAdt; dnZdt];
    end

    % BSTR function
    function [t, nA, nZ] = BSTR_variables()
        % set the initial values and stopping criterion
        ind_0 = 0;
        nA0 = CA0*V;
        nZ0 = CZ0*V;
        dep_0 = [nA0; nZ0];
        stop_var = 0;
        stop_val = tf;
    
        % solve the BSTR reactor design equations
        odes_are_stiff = false;
        [ind, dep, flag, message] = solve_ivodes(ind_0, dep_0...
            , stop_var, stop_val, @BSTR_derivatives, odes_are_stiff);
    
        % print a warning if there was a problem solving the design 
        % equations
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
            disp(['         ',message])
        end
        
        % extract the corresponding sets of values for each of the
        % variables
        t = ind;
        nA = dep(:,1);
        nZ = dep(:,2);
    end

    % deliverables function
    function deliverables()
        % get the BSTR variables
        [t, nA, nZ] = BSTR_variables();
    
        % plot, display, and save the results
        figure;
        plot(t,nA,t,nZ,'LineWidth',2)
        set(gca, 'FontSize', 14);
        legend({'A','Z'},'Location','northeast','FontSize',14)
        xlabel('Time (min)','FontSize', 14)
        ylabel('Amount (moless)','FontSize', 14)
        saveas(gcf,"profiles.pdf")
    end
    
    % execution statement
    deliverables();
end