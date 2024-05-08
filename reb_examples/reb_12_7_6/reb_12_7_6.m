function reb_12_7_6
%REB_12_7_6 Calculations for Example 12.7.6 of Reaction Engineering Basics
    
    % constants available to all functions
    % given
    V = 500; % cm^3
    Vdot_in = 1.0; % cm^3 /s
    CA_in = 0.015; % mol /cm^3
    CB_in = 0.015; % mol /cm^3
    T_in_u = 50 + 273.15; % K
    T_in_ign = 91 + 273.15; % K
    T_in_ext = 3 + 273.15; % K
    T_0_u = 140 + 273.15; % K
    T_0_ign = 96 + 273.15; % K
    T_0_ext = 200 + 273.15; % K
    fA_0_u = 0.399;
    fA_0_ign = 0.032;
    fA_0_ext = 0.881;
    Cp = 0.35*4.184; % J /g /K
    rho = 0.93; % g /cm^3
    dH = -20000; % J /mol
    k0 = 3.24E12; % cm^3 /mol /s
    E = 105000; % J/mol
    % known
    R = 8.314; % J /mol /K
    % calculated
    nA_in = Vdot_in*CA_in;
    nB_in = Vdot_in*CB_in;

    % make T_in, T_0, and fA_0 available to all functions
    T_in = nan;
    T_0 = nan;
    fA_0 = nan;

    % derivatives function
    function derivs = derivatives(~, dep)
        % extract the dependent variables for this integration step
        nA = dep(1);
        nB = dep(2);
        nY = dep(3);
        nZ = dep(4);
        T = dep(5);

        % calculate the rate
        k = k0*exp(-E/R/T);
        CA = nA/Vdot_in;
        CB = nB/Vdot_in;
        r_1 = k*CA*CB;

        % evaluate the derivatives
        dnAdt = Vdot_in/V*(nA_in - nA - V*r_1);
        dnBdt = Vdot_in/V*(nB_in - nB - V*r_1);
        dnYdt = Vdot_in/V*(-nY + V*r_1);
        dnZdt = Vdot_in/V*(-nZ + V*r_1);
        dTdt = -(Cp*Vdot_in*rho*(T - T_in)+V*r_1*dH)/(V*Cp*rho);

        % return the derivatives
        derivs = [dnAdt; dnBdt; dnYdt; dnZdt; dTdt];
    end

    % reactor model
    function [t, nA, nB, nY, nZ, T] = profiles(t_f)
        % set the initial values
        ind_0 = 0.0;
        extent = fA_0*nA_in;
        dep_0 = [nA_in - extent; nB_in - extent; extent; extent;
            T_0];

        % define the stopping criterion
        f_var = 0;
        f_val = t_f;
        
        % solve the IVODEs
        odes_are_stiff = false;
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
        nY = dep(:,3);
        nZ = dep(:,4);
        T = dep(:,5);
    end

    % function that performs the analysis
	function perform_the_analysis()

        % analysis for perturbation from the unsteady state

        T_0 = T_0_u;
        T_in = T_in_u;
        fA_0 = fA_0_u;

        % solve the reactor design equations
        t_f = 500.0; % s
        [t_u, ~, ~, ~, ~, T_u] = profiles(t_f);

        % analysis for perturbation from the steady state near
        % the ignition point

        T_0 = T_0_ign;
        T_in = T_in_ign;
        fA_0 = fA_0_ign;

        % solve the reactor design equations
        t_f = 10000.0; % s
        [t_ign, ~, ~, ~, ~, T_ign] = profiles(t_f);

        % analysis for perturbation from the steady state near 
        % the extinction point

        T_0 = T_0_ext;
        T_in = T_in_ext;
        fA_0 = fA_0_ext;

        % solve the reactor design equations
        t_f = 4000; % s
        [t_ext, ~, ~, ~, ~, T_ext] = profiles(t_f);

        % create, display and save the graphs
        figure; % unsteady state scenario
        plot(t_u/60.0, T_u - 273.15,'LineWidth',2)
        set(gca, 'FontSize', 14);
        title('Perturbation from an Unsteady State','FontSize',14)
        xlabel('Time (min)','FontSize', 14)
        ylabel('Reacting Fluid Temperature (°C)','FontSize', 14)
        saveas(gcf,"unsteady.png")

        figure; % near ignition scenario
        plot(t_ign/60.0, T_ign - 273.15,'LineWidth',2)
        set(gca, 'FontSize', 14);
        title('Ignition','FontSize',14)
        xlabel('Time (min)','FontSize', 14)
        ylabel('Reacting Fluid Temperature (°C)','FontSize', 14)
        saveas(gcf,"ignition.png")

        figure; % near extinction scenario
        plot(t_ext/60.0, T_ext - 273.15,'LineWidth',2)
        set(gca, 'FontSize', 14);
        title('Extinction','FontSize',14)
        xlabel('Time (min)','FontSize', 14)
        ylabel('Reacting Fluid Temperature (°C)','FontSize', 14)
        saveas(gcf,"extinction.png")
    end

    % perform the analysis
    perform_the_analysis()
end