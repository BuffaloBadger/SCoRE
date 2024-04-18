function reb_10_5_2
%REB_10_5_2 Calculations for Reaction Engineering Basics Example 10.5.2

    % given, known and calculated constants available to all functions
    % given
    k0_1 = 1.83E12; % L/mol/h
    E_1 = 18000; % cal/mol
    k0_2 = 5.08E13; % L/mol/h
    E_2 = 20500; % cal/mol
    dH_1 = -9000; % cal/mol
    dH_2 = -7800; % cal/mol
    Cp = 863; % cal/L/K
    CA_0 = 2.0; % mol/L
    CB_in = 0.5; % mol/L
    T_0 = 40 + 273.15; % K
    T_in = T_0;
    P = 1.0; % atm
    V_0 = 2000; % L
    V_B = 8000; % L
    t_f = 8.0; % h
    % known
    Re = 1.987; % cal/mol/K
    Rw = 0.08206; % L*atm/mol/K
    % calculated
    nA_0 = CA_0*V_0;

    % make t_add available to all functions
    t_add = nan(1,1);

    % derivatives function
    function derivs = derivatives(t, dep)
        % extract necessary dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        T = dep(5);
        
        % calculate the reacting fluid volume
        if t<t_add
            Vdot_in = V_B/t_add;
            V = V_0 + Vdot_in*t;
        else
            Vdot_in = 0.0;
            V = V_0 + V_B;
        end

        % calculate the inlet molar flow rate of B
        nDotB_in = Vdot_in*CB_in;

        % calculate rate coefficients
        k_1 = k0_1*exp(-E_1/Re/T);
        k_2 = k0_2*exp(-E_2/Re/T);

        % calculate the rate
        CA = nA/V;
        CB = nB/V;
        r_1 = k_1*CA*CB;
        r_2 = k_2*CB^2;

        % evaluate the derivatives
        dAdt = -V*r_1;
        dBdt = nDotB_in - V*(r_1 + 2*r_2);
        dDdt = V*r_1;
        dUdt = V*r_2;
        dTdt = (-Vdot_in*Cp*(T-T_in) - r_1*V*dH_1 -r_2*V*dH_2 ...
            + P*Vdot_in*Re/Rw)/V/Cp;

        % return the derivatives
        derivs = [dAdt; dBdt; dDdt; dUdt; dTdt];
    end

    % reactor model function
    function [t, nA, nB, nD, nU, T] = profiles()
        % set the initial values for the first stage
        ind_0 = 0.0;
        dep_0 = [nA_0; 0.0; 0.0; 0.0; T_0];

        % define the stopping criterion for the first stage
        f_var = 0;
        f_val = t_add;
        
        % solve the IVODEs for the first stage
        odes_are_stiff = false;
        [t1, dep1, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: Stage 1 solution may not be accurate!')
        end

        % set the initial values for the second stage
        ind_0 = t1(end);
        dep_0 = dep1(end,1:5)';

        % define the stopping criterion for the second stage
        f_val = t_f;

        % solve the design equations for the second stage
        [t2, dep2, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);

        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: Stage 2 solution may not be accurate!')
        end

        % combine the stages and return the profiles
        t = [t1; t2];
        nA = [dep1(:,1); dep2(:,1)];
        nB = [dep1(:,2); dep2(:,2)];
        nD = [dep1(:,3); dep2(:,3)];
        nU = [dep1(:,4); dep2(:,4)];
        T = [dep1(:,5); dep2(:,5)];
    end

    % function that performs the analysis
    function perform_the_analysis()
        % set the add times
        add_times = [1.0; 3.0; 5.0; 7.0];

        % allocate storage for the results
        f_B = nan(4,1);
        S_DoverU = nan(4,1);
        Y_DfromB = nan(4,1);

        % loop through the add times
        for iAdd = 1:4
            % set t_add
            t_add = add_times(iAdd);

            % solve the reactor design equations
            [~, ~, nB, nD, nU, ~] = profiles();

            % calculate the quantities of interest
            f_B(iAdd) = (V_B*CB_in - nB(end))/V_B/CB_in;
            S_DoverU(iAdd) = nD(end)/nU(end);
            Y_DfromB(iAdd) = nD(end)/V_B/CB_in;
        end

        % convert to percents
        f_B = 100.0*f_B;
        Y_DfromB = 100.0*Y_DfromB;

        % tabulate, display and save the results
        results_table = table(add_times, f_B, S_DoverU, Y_DfromB);
        disp(results_table)
        writetable(results_table,'results.csv');
    end

    % perform the analysis
    perform_the_analysis();
end