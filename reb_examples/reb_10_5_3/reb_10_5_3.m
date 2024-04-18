function reb_10_5_3
%REB_10_5_3 Calculations for Reaction Engineering Basics Example 10.5.3

    % constants available to all functions
    % given
    k0_1 = 1.192E15; % /min
    E_1 = 97600.0; % J/mol
    dH_1 = -58615.0; % J/mol
    V_ex = 300.0; % cm^3
    Tex_in = 60 + 273.15; % K
    Vdot_ex = 250.0; % cm^3/min
    UA = 260.0*4.184; % J/min/K
    P = 1.0; % atm
    VW_0 = 67.0; % cm^3
    VZ_0 = 283.0; % cm^3
    Vacid = 0.3; % cm^3
    T_0 = 60 + 273.15; % K
    Tex_0 = 60 + 273.15; %K
    Cp = 2.68; % J/cm^3/K
    V_A = 350.0; % cm^3
    T_in = 21 + 273.15; % K
    T_f = 65 + 273.15; % K
    T_max = 95 + 273.15; % K
    rho_A = 1.082; % g/cm^3
    rho_W = 1.0; % g/cm^3
    rho_Z = 1.049; % g/cm^3
    M_A = 102; % g/mol
    M_W = 18; % g/mol
    M_Z = 60; % g/mol
    Cp_A = 168.2; % J/mol/K
    rho_ex = 1.0; % g/cm^3
    Cp_ex = 1.0*4.184; % J/g/K
    % known
    Re = 1.987*4.184; % J/mol/K
    Rw = 82.06; % cm^3*atm/mol/K
    % calculated
    V_0 = VW_0 + VZ_0 + Vacid;
    m_ex = Vdot_ex*rho_ex;
    nW_0 = VW_0*rho_W/M_W;
    nZ_0 = VZ_0*rho_Z/M_Z;
    

    % make Vdot_in available to all functions
    Vdot_in = nan(1,1);

    % derivatives function
    function derivs = derivatives(t, dep)
        % extract necessary dependent variables for this integration step
        nA = dep(1);
        T = dep(4);
        T_ex = dep(5);

        % calculate the time when all of the A has been added
        t_sb = V_A/Vdot_in;
        
        % calculate the reacting fluid volume and set the molar feed rate
        if t<t_sb
            V = V_0 + Vdot_in*t;
            nDotA_in = Vdot_in*rho_A/M_A;
            W_exp = P*Vdot_in*Re/Rw;
        else
            V = V_0 + V_A;
            nDotA_in = 0.0;
            W_exp = 0.0;
        end

        % calculate rate coefficient
        k_1 = k0_1*exp(-E_1/Re/T);

        % calculate the rate
        CA = nA/V;
        r_1 = k_1*CA;

        % calculate the rate of heat exchange
        Qdot = UA*(T_ex - T);

        % evaluate the derivatives
        dAdt = nDotA_in - V*r_1;
        dWdt = -V*r_1;
        dZdt = V*r_1;
        dTdt = (Qdot - nDotA_in*Cp_A*(T-T_in) - V*r_1*dH_1 + W_exp)...
            /(V*Cp);
        dTexdt = (-Qdot - m_ex*Cp_ex*(T_ex - Tex_in))...
            /(rho_ex*V_ex*Cp_ex);

        % return the derivatives
        derivs = [dAdt; dWdt; dZdt; dTdt; dTexdt];
    end

    % reactor model function
    function [t, nA, nW, nZ, T, T_ex] = profiles()
        % set the initial values for the first stage
        ind_0 = 0.0;
        dep_0 = [0.0; nW_0; nZ_0; T_0; Tex_0];

        % define the stopping criterion for the first stage
        t_sb = V_A/Vdot_in;
        f_var = 0;
        f_val = t_sb;
        
        % solve the IVODEs for the first stage
        odes_are_stiff = false;
        [t, dep, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);
    
        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: Stage 1 solution may not be accurate!')
        end

        % extract the individual profiles
        nA = dep(:,1);
        nW = dep(:,2);
        nZ = dep(:,3);
        T = dep(:,4);
        T_ex = dep(:,5);

        % set the initial values for the second stage
        ind_0 = t(end);
        dep_0 = dep(end,1:5)';

        % define the stopping criterion for the second stage
        f_var = 4;
        f_val = T_f;

        % solve the design equations for the second stage
        [t2, dep2, flag] = solve_ivodes(ind_0, dep_0, f_var, f_val...
            , @derivatives, odes_are_stiff);

        % check that the solution was found
        if flag <= 0
            disp(' ')
            disp('WARNING: Stage 2 solution may not be accurate!')
        end

        % combine the stages and return the profiles
        t = [t; t2];
        nA = [nA; dep2(:,1)];
        nW = [nW; dep2(:,2)];
        nZ = [nZ; dep2(:,3)];
        T = [T; dep2(:,4)];
        T_ex = [T_ex; dep2(:,5)];
    end

    % function that performs the analysis
    function perform_the_analysis()
        % set a range of feed rates
        %feed_rates = linspace(30,300,100);
        feed_rates = linspace(37,41,100);

        % allocate storage for the results
        t_proc = nan(100,1);

        % loop through the feed rates
        for iFeedRate = 1:100
            % set Vdot_in
            Vdot_in = feed_rates(iFeedRate);

            % solve the reactor design equations
            [t, ~, ~, ~, T, ~] = profiles();

            % save the processing time if the temperature constraint is
            % satisfied
            if max(T) < T_max
                t_proc(iFeedRate) = t(end);
            end
        end

        % find the optimum
        [t_opt,i_opt] = min(t_proc);
        Vdot_opt = feed_rates(i_opt);
        t_sb = V_A/Vdot_opt;

        % solve the reactor design equations using the optimum feed rate
        Vdot_in = Vdot_opt;
        [t, nA, ~, ~, T, Tex] = profiles();
        f_A = 100*(V_A*rho_A/M_A - nA(end))/(V_A*rho_A/M_A);

        % tabulate the results
        item = ["Minimum Processing Time"...
            ; "Semi-Batch Processing Time"...
            ; "Optimum Feed Rate"...
            ; "Maximum Temperature"...
            ; "Maximum Cooling Water Temperature"...
            ; "Acetic Anhydride Conversion"];
        value = [t_opt...
            ; t_sb...
            ; Vdot_in...
            ; max(T - 273.15)...
            ; (max(Tex - 273.15))...
            ; f_A];
        units = ["min"...
            ; "min"...
            ; "cm^3^ min^-1^"...
            ; "°C"...
            ; "°C"...
            ; "%"];
        results_table = table(item, value, units);
        disp(results_table)
        writetable(results_table,'results.csv')

        % calculate the acetic acid concentration profile
        nData = length(t);
        CA = nan(nData,1);
        for i = 1:nData
            if t(i) < t_sb
                V = V_0 + Vdot_in*t(i);
            else
                V = V_0 + V_A;
            end
        CA(i) = nA(i)/V;
        end

        % create, show and save the graphs
        figure
        plot(t,T-273.15,'LineWidth',2)
        xlabel("Time (min)","FontSize",14)
        ylabel("Temperture (°C)","FontSize",14)
        set(gca,'FontSize',14)
        saveas(gcf,'temperature_profile.png')

        figure
        plot(t,CA,'LineWidth',2)
        xlabel("Time (min)","FontSize",14)
        ylabel("Acetic Anhydride Concentration (mol cm^-^3)","FontSize",14)
        set(gca,'FontSize',14)
        saveas(gcf,'concentration_profile.png')
    end

    % perform the analysis
    perform_the_analysis();
end