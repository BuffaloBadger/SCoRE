function reb_15_4_4
%REB_15_4_4 Calculations for Example 15.4.4 of Reaction Engineering Basics
    % constants available to all functions
	% given
    V_pfr1 = 60.0; % L
    V_pfr2 = 40.0; % L
    CA_0 = 1.0; % mol /L
    T_0 = 60 + 273.15; % K
    Vdot_0 = 0.55; % L /min
    k0 = 2.63E7; % L /mol /min
    E = 62000; % J /mol
    dH = -35000; % J /mol
    Cp = 800; % J /L /K
	% known
    R = 8.314; % J /mol /K
	% calculated
    nDotA_0 = Vdot_0*CA_0;

    % make Vdot available to all functions
    Vdot = nan;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nDot_A = dep(1);
        T = dep(4);

        % calculate the rate
        r = k0*exp(-E/R/T)*nDot_A^2/Vdot^2;

        % evaluate the derivatives
        dnAdV = -2*r;
        dnYdV = r;
        dnZdV = r;
        dTdV = -r*dH/Vdot/Cp;

        % return the derivatives
        derivs = [dnAdV; dnYdV; dnZdV; dTdV];
    end

    % reactor model
    function [V, nDotA, nDotY, nDotZ, T] = profiles(f_val)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [Vdot*CA_0; 0.0; 0.0; T_0];

        % define the stopping criterion
        f_var = 0;
        
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
        nDotY = dep(:,2);
        nDotZ = dep(:,3);
        T = dep(:,4);
    end

    % function that performs the analysis
	function perform_the_analysis()
        % case 1, equal flow rates
        % set Vdot
        Vdot = Vdot_0/2.0;

        % solve the reactor design equations for reactor R1
        [~, nDotA, ~, ~, ~] = profiles(V_pfr1);
        nDotA_3 = nDotA(end);

        % solve the reactor design equations for reactor R2
        [~, nDotA, ~, ~, ~] = profiles(V_pfr2);
        nDotA_4 = nDotA(end);

        % calculate the conversion
        nDotA_5 = nDotA_3 + nDotA_4;
        fA_a = 100.0*(nDotA_0 - nDotA_5)/nDotA_0;

        % case 2, equal space velocities
        Vdot_1 = Vdot_0*V_pfr1/(V_pfr1 + V_pfr2);
        Vdot_2 = Vdot_0 - Vdot_1;

        % solve the reactor design equations for reactor R1
        Vdot = Vdot_1;
        [~, nDotA, ~, ~, ~] = profiles(V_pfr1);
        nDotA_3 = nDotA(end);

        % solve the reactor design equations for reactor R2
        Vdot = Vdot_2;
        [~, nDotA, ~, ~, ~] = profiles(V_pfr2);
        nDotA_4 = nDotA(end);

        % calculate the conversion
        nDotA_5 = nDotA_3 + nDotA_4;
        fA_b = 100.0*(nDotA_0 - nDotA_5)/nDotA_0;

        % tabulate the results
        item = ["Conversion for Equal Flow Rates"
            "Conversion for Equal Space Times"];
        value = [fA_a; fA_b];
        units = ["%";"%"];
        results_table = table(item,value,units);

        % display the results
        disp(' ')
        disp(results_table)

        % save the results
        writetable(results_table,'results.csv');
    end

    % perform the analysis
    perform_the_analysis()
end