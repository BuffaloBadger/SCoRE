function reb_J_6_1
%REB_J_6_1 Calculations for Example J.6.1 of Reaction Engineering Basics

    % given and known constants available to all functions
    n_A_0 = 190.; % mol
    n_B_0 = 190.; % mol
    n_Y_0 = 0.; % mol
    n_Z_0 = 0.; % mol
    T_0 = 450.; % K
    P_0 = 7.0; % atm
    t_f = 2.0; % h
    V = 2.0; % m^3
    dH_1 = -3470; % cal/mol
    Cp_A = 7.5; % cal/mol/K
    Cp_B = 8.5; % cal/mol/K
    Cp_Y = 12.1; % cal/mol/K
    Cp_Z = 5.7; % cal/mol/K
    k_0_1 = 83; % m^3 /mol /h
    E_1 = 10200; % cal/mol
    Re = 1.987; % cal/mol/K
    Rw = 8.206E-5; % m^3 atm/mol/K

    % reactor design equations as derivative expressions
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        n_A = dep(1);
        n_B = dep(2);
        n_Y = dep(3);
        n_Z = dep(4);
        T = dep(5);
        % P = dep(6); not needed to evaluate the derivatives

        % create empty mass matrix
        mass_matrix = zeros(6,6);

        % add 1 on the diagonal for the first 4 rows
        mass_matrix(1,1) = 1.0;
        mass_matrix(2,2) = 1.0;
        mass_matrix(3,3) = 1.0;
        mass_matrix(4,4) = 1.0;

        % Add the elements for the energy balance
        mass_matrix(5,5) = n_A*Cp_A + n_B*Cp_B + n_Y*Cp_Y + n_Z*Cp_Z;
        mass_matrix(5,6) = -V*Re/Rw;

        % Add the elements for the ideal gas law equation
        mass_matrix(6,1) = Rw*T;
        mass_matrix(6,2) = Rw*T;
        mass_matrix(6,3) = Rw*T;
        mass_matrix(6,4) = Rw*T;
        mass_matrix(6,5) = Rw*(n_A + n_B + n_Y + n_Z);
        mass_matrix(6,6) = -V;

        % calculate the rate
        C_A = n_A/V;
        C_B = n_B/V;
        r = k_0_1*exp(-E_1/Re/T)*C_A*C_B;

        % Create right side vector
        rhs = [-V*r; -V*r; V*r; V*r; -V*r*dH_1; 0];

        % Evaluate the derivatives
        derivs = mass_matrix\rhs;
    end

    % reactor model
    function [t, n_A, n_B, n_Y, n_Z, T, P] = profiles()
        % set the initial and final values
        ind_0 = 0.0;
        dep_0 = [n_A_0; n_B_0; n_Y_0; n_Z_0; T_0; P_0];
        f_var = 0;
        f_val = t_f;

        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);

        % Check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract the dependent variable profiles
        n_A = dep(:,1);
        n_B = dep(:,2);
        n_Y = dep(:,3);
        n_Z = dep(:,4);
        T = dep(:,5);
        P = dep(:,6);
    end

    % function that performs the analysis
    function perform_the_analysis()
        % solve the reactor design equations
        [t, n_A, n_B, n_Y, n_Z, T, P] = profiles();
     
        % Tabulate the results
        results_table = table(t,n_A,n_B,n_Y,n_Z,T,P);
    
        % Display the results
        disp(' ')
        disp(results_table)
    
        % Save the results
        results_file = "results.csv";
        writetable(results_table,results_file);
    end
    
    % perform the analysis
    perform_the_analysis();
end