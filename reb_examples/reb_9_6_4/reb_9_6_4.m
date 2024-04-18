function reb_9_6_4
%REB_9_6_4 Calculations for Example 9.6.4 of Reaction Engineering Basics

    % given, known, and calculated constants available to all functions
    % given
    k_0_1 = 2.59E9; % /min
    E_1 = 16500; % cal/mol
    dH_1 = -22200; % cal/mol
    CA_0 = 2.0; % mol /L
    T_0 = 23 + 273.15; % K
    Cp = 440.0; % cal/L/K
    V = 4.0; % L
    V_shell = 0.5; % L
    A_shell = 0.6; % ft^2
    U_shell = 1.13E4/60.0; % cal/ft^2/min/K
    Tex_in = 20 + 273.15; % K
    rho_ex = 1.0; % g/cm^3
    Cp_ex = 1.0; % cal/g/K
    U_coil = 3.8E4/60.0; % cal/ft^2/min/K
    A_coil = 0.23; % ft^2
    T_coil = 120 + 273.15; % K
    Tex_0 = 23 + 273.15; % K
    T_1_f = 50 + 273.15; % K
    T_f = 25 + 273.15; % K
    t_turn = 25.0; % min
    % known
    Re = 1.987; % cal/mol/K
    % calculated
    nA_0 = CA_0*V;

    % make the current coolant flow rate available to all functions
    mDot_ex = nan;

    % derivatives function for the first stage of operation
    function derivs = first_stage_derivatives(~, dep)
        % extract necessary dependent variables for this integration step
        nA = dep(1);
        T = dep(3);
        Tex = dep(4);

        % calculate the rate
        CA = nA/V;
        k_1 = k_0_1*exp(-E_1/Re/T);
        r_1 = k_1*CA;

        % calculate the heat exhange rates
        Qdot_shell = U_shell*A_shell*(Tex-T);
        Qdot_coil = U_coil*A_coil*(T_coil-T);

        % evaluate the derivatives
        dnAdt = -V*r_1;
        dnZdt = V*r_1;
        dTdt = (Qdot_shell + Qdot_coil - V*r_1*dH_1)/V/Cp;
        dTexdt = -Qdot_shell/rho_ex/V_shell/Cp_ex;

        % return the derivatives
        derivs = [dnAdt; dnZdt; dTdt; dTexdt];
    end

    % derivatives function for the second stage of operation
    function derivs = second_stage_derivatives(~, dep)
        % extract necessary dependent variables for this integration step
        nA = dep(1);
        T = dep(3);
        Tex = dep(4);

        % calculate the rate
        CA = nA/V;
        k_1 = k_0_1*exp(-E_1/Re/T);
        r_1 = k_1*CA;

        % calculate the heat exhange rates
        Qdot_shell = U_shell*A_shell*(Tex-T);

        % evaluate the derivatives
        dnAdt = -V*r_1;
        dnZdt = V*r_1;
        dTdt = (Qdot_shell - V*r_1*dH_1)/V/Cp;
        dTexdt = -(Qdot_shell + mDot_ex*Cp_ex*(Tex-Tex_in))/rho_ex/V_shell/Cp_ex;

        % return the derivatives
        derivs = [dnAdt; dnZdt; dTdt; dTexdt];
    end

    % reactor model for the first stage of operation
    function [t, nA, nZ, T, Tex] = first_stage_profiles()
        % set the initial values
        ind_0 = 0.0;
        dep_0 = [nA_0; 0.0; T_0; Tex_0];

        % define the stopping criterion
        f_var = 3;
        f_val = T_1_f;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @first_stage_derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nA = dep(:,1);
        nZ = dep(:,2);
        T = dep(:,3);
        Tex = dep(:,4);
    end

    % reactor model for the second stage of operation
    function [t, nA, nZ, T, Tex] = second_stage_profiles(t_0, nA_0...
            , nZ_0, T_0, Tex_0)
        % set the initial values
        ind_0 = t_0;
        dep_0 = [nA_0; nZ_0; T_0; Tex_0];

        % define the stopping criterion
        f_var = 3;
        f_val = T_f;
        
        % solve the IVODEs
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @second_stage_derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: The ODE solution may not be accurate!')
        end

        % extract and return the dependent variable profiles
        nA = dep(:,1);
        nZ = dep(:,2);
        T = dep(:,3);
        Tex = dep(:,4);
    end

    % function that performs the analysis
        function perform_the_analysis()
    
        % choose a range of coolant flow rates
        coolant_flows = linspace(100.0, 250.0, 100);
    
        % calculate the net rate for each coolant flow rate
        net_rates = nan(100,1);
        for i = 1:100

            % set the coolant flow rate to be used in the calculations
            mDot_ex = coolant_flows(i);

            % solve the reactor design equations for the first stage
            [t_1, nA_1, nZ_1, T_1, Te_1] = first_stage_profiles();

            % solve the reactor design equations for the second stage
            [t, ~, nZ, ~, ~] = second_stage_profiles(t_1(end)...
                , nA_1(end), nZ_1(end), T_1(end), Te_1(end));

            % calculate the net rate
            net_rates(i) = nZ(end)/(t(end) + t_turn);
        end

        % find the coolant flow where the net rate is maximized
        [r_net_max,i_max] = max(net_rates);
        mDot_max = coolant_flows(i_max);

        % solve the reactor design equations using the optimum coolant flow
        mDot_ex = mDot_max;
        [t_1, nA_1, nZ_1, T_1, Te_1] = first_stage_profiles();
        [t_2, nA_2, ~, T_2, ~] = second_stage_profiles(t_1(end)...
                , nA_1(end), nZ_1(end), T_1(end), Te_1(end));

        % combine the profiles
        t = [t_1 ; t_2];
        nA = [nA_1; nA_2];
        T = [T_1; T_2];

        % calculate the conversion vs time at the optimum cooland flow rate
        pct_conversion = 100*(nA_0 - nA)/nA_0;
    
        % tabulate the results
        item = ["Optimum Coolant Flow Rate";"Maximum Net Rate"];
        value = [mDot_max; r_net_max];
        units = ["g min^-1^"; "mol L^-1^ min^-1^"];
        results_table = table(item, value, units);

        % display the results
        disp(' ')
        disp(['Optimum Coolant Flow Rate: ',num2str(mDot_max,3)...
            ,' g/min'])
        disp(['Maximum Net Rate: ', num2str(r_net_max,3), ' mol/L/min'])

        % save the results
        writetable(results_table,'results.csv');
    
        % display and save the graphs
        figure; % net rate vs coolant flow
        plot(coolant_flows, net_rates, 'LineWidth', 2)
        set(gca, 'FontSize', 14);
        xlabel('Coolant Flow Rate (g/min)','FontSize', 14)
        ylabel('Net Rate (mol/min)','FontSize', 14)
        saveas(gcf,"net_rate_vs_coolant_flow.png")

        figure; % conversion profile
        plot(t, pct_conversion, 'LineWidth', 2)
        set(gca, 'FontSize', 14);
        xlabel('Reaction Time (min)','FontSize', 14)
        ylabel('Conversion (%)','FontSize', 14)
        saveas(gcf,"conversion_profile.png")

        figure; % temperature profile
        plot(t, T - 273.15, 'LineWidth', 2)
        set(gca, 'FontSize', 14);
        xlabel('Reaction Time (min)','FontSize', 14)
        ylabel('Temperature (°C)','FontSize', 14)
        saveas(gcf,"temperature_profile.png")
   
    end

    % perform the analysis
    perform_the_analysis();
end