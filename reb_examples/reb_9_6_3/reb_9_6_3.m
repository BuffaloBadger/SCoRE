function reb_9_6_3
%REB_9_6_3 Calculations for Example 9.6.3 of Reaction Engineering Basics

    % given, known and calculated constants available to all functions
    % given
    PA_0 = 1.; % atm
    PB_0 = 2.; % atm
    T_0 = 25 + 273.15; % K
    V = 2.; % L
    Aex = 600.; % cm^2
    Uex = 0.6; % cal /cm^2 /min /K
    Tex = 30 + 273.15; % K
    k_0_1 =  3.34E9; % mol /cm^3 /min /atm^2
    k_0_2 = 4.99E9; % mol /cm^3 /min /atm^2
    E_1 = 20.5E3; % cal /mol
    E_2 = 21.8E3; % cal /mol
    dH_1 = -6300.; % cal /mol
    dH_2 = -6900.; % cal /mol
    Cp_A = 7.4; % cal /mol /K
    Cp_B = 8.6; % cal /mol /K
    Cp_D = 10.7; % cal /mol /K
    Cp_Z = 5.2; % cal /mol /K
    Cp_U = 10.3; % cal /mol /K
    % known
    Re = 1.987; % cal /mol /K
    Rw = 82.057; % cm^3 atm /mol /K
    % calculated
    nA_0 = PA_0*V/Rw/T_0;
    nB_0 = PB_0*V/Rw/T_0;
    P_0 = PA_0 + PB_0;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract necessary dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        nD = dep(3);
        nZ = dep(4);
        nU = dep(5);
        T = dep(6);

        % calculate the rate
        PA = nA*Rw*T/V;
        PB = nB*Rw*T/V;
        PD = nD*Rw*T/V;
        k_1 = k_0_1*exp(-E_1/Re/T);
        k_2 = k_0_2*exp(-E_2/Re/T);
        r_1 = k_1*PA*PB;
        r_2 = k_2*PD*PB;

        % calculate the rate of heat exchange
        Qdot = Uex*Aex*(Tex-T);

        % create empty mass matrix
        mass_matrix = zeros(7,7);

        % add 1 on the diagonal for the first 4 rows
        mass_matrix(1,1) = 1.0;
        mass_matrix(2,2) = 1.0;
        mass_matrix(3,3) = 1.0;
        mass_matrix(4,4) = 1.0;
        mass_matrix(5,5) = 1.0;

        % Add the elements for the energy balance
        mass_matrix(6,6) = nA*Cp_A + nB*Cp_B + nD*Cp_D + nZ*Cp_Z + nU*Cp_U;
        mass_matrix(6,7) = -V*Re/Rw;

        % Add the elements for the ideal gas law equation
        mass_matrix(7,1) = Rw*T;
        mass_matrix(7,2) = Rw*T;
        mass_matrix(7,3) = Rw*T;
        mass_matrix(7,4) = Rw*T;
        mass_matrix(7,5) = Rw*T;
        mass_matrix(7,6) = Rw*(nA + nB + nD + nZ + nU);
        mass_matrix(7,7) = -V;

        % Create right side vector
        rhs1 = -r_1*V;
        rhs2 = (-r_1 -r_2)*V;
        rhs3 = (r_1-r_2)*V;
        rhs4 = (r_1+r_1)*V;
        rhs5 = r_2*V;
        rhs6 = Qdot -(r_1*dH_1 + r_2*dH_2)*V;
        rhs7 = 0.0;
        rhs = [rhs1; rhs2; rhs3; rhs4; rhs5; rhs6; rhs7];

        % Evaluate the derivatives
        derivs = mass_matrix\rhs;
    end

    % reactor model function
    function [t, nA, nB, nD, nZ, nU, T, P] = profiles(t_f)
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_0; nB_0; 0.0; 0.0; 0.0; T_0; P_0];

        % define the stopping criterion
        f_var = 0;
        f_val = t_f;
        
        % solve the IVODEs
        odes_are_stiff = true;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nA = dep(:,1);
        nB = dep(:,2);
        nD = dep(:,3);
        nZ = dep(:,4);
        nU = dep(:,5);
        T = dep(:,6);
        P = dep(:,7);
    end

    % function that performs the analysis
    function perform_the_analysis()
    
        % set a range of reaction times
        times = linspace(1, 60.0, 1000);
    
        % create storage for the yields
        yields = nan(1000,1);
    
        % calculate the yield at each reaction time
        for i=1:1000
            [~, ~, ~, nD, ~, ~, ~, ~] = profiles(times(i));
            yields(i) = nD(end)/nA_0;
        end
    
        % find the reaction time where the yield is maximum
        [~, i_max] = max(yields);
        t_max = times(i_max);
    
        % solve the reactor design equations at the optimum reaction time
        [~, nA, ~, nD, ~, ~, ~, ~] = profiles(t_max);
    
        % calculate the other quantities of interest
        yield_max = nD(end)/nA_0;
        fA_pct = 100*(nA_0 - nA(end))/nA_0;
    
        % tabulate the results
        item = ["t_max";"Y_D_A";"Conversion"];
        value = [t_max; yield_max; fA_pct];
        units = ["min"; "mol D per initial mol A"; "%"];
        results_table = table(item,value,units);
    
        % display the results
        disp(' ')
        disp(['Optimum reaction time: ', num2str(t_max,3), ' min'])
        disp(['Yield: ', num2str(yield_max,3), ' mol D per initial mol A'])
        disp(['Conversion of A: ', num2str(fA_pct,3), ' %'])
        disp(' ')
    
        % save the results
        writetable(results_table,'results.csv');
    
        % plot the yield vs. the reaction time
        figure; 
        plot(times, yields,'LineWidth',2)
        set(gca, 'FontSize', 14);
        xlabel('Reaction Time (min)', 'FontSize', 14)
        ylabel('Yield (mol D per initial mol A)', 'FontSize', 14)
    
        % save the graph
        saveas(gcf,"yield_vs_t.png")
    end

    % perform the analysis
    perform_the_analysis();
end